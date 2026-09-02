"""FRUTON -> AMBER atom-name normalization.

FRUTON's protonation stage (gmx pdb2gmx) writes hydrogens with GROMACS-style
naming, which sometimes differs from AMBER ff14SB conventions.  Bare tleap
then errors out (missing atom / unknown residue) on the mismatch.

Six categories of transformation are needed (per rename_h_atoms.py that we
historically ran as a post-processing script on MMBSA_200):

1. Methylene hydrogens: HX1/HX2 -> HX2/HX3 (AMBER numbers from 2, not 1)
   Applied to: HB (many residues except GLY/ALA/ILE/VAL/THR),
               HG (ASN/ASP/CYS/GLN/GLU/LYS/MET/PRO/ARG),
               HD (LYS/PRO/ARG),
               HE (LYS)

2. ILE special: CD -> CD1, HD1/HD2/HD3 -> HD11/HD12/HD13, HG11/HG12 -> HG12/HG13.

3. C-terminal oxygens: OC1/OC2 deleted (tleap re-adds O/OXT from topology).
   NOTE: FRUTON's ``prepared_structure.py`` already handles OC1/OC2 for
   internal fragments; this module handles the C-terminal case that
   escapes that path.

4. CYM (deprotonated Cys) HG: deleted (thiolate anion has no thiol H).

5. HIS residue name normalization: rename to HID / HIE / HIP based on
   which imidazole H atoms are present.
      HD1 only            -> HID
      HE2 only            -> HIE
      HD1 AND HE2         -> HIP

6. N-terminal PRO: H1/H2 -> H2/H3 (secondary amine numbering).
"""
from __future__ import annotations

from pathlib import Path

# Residues with beta-methylene (HB1/HB2 -> HB2/HB3). Excludes:
#   GLY (no Cbeta), ALA (methyl HB1/HB2/HB3 unchanged),
#   ILE/THR/VAL (single HB).
_HB_METHYLENE = {
    "ARG", "ASN", "ASP", "CYS", "CYM", "CYX", "GLN", "GLU", "GLH", "HIS",
    "HID", "HIE", "HIP", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "TRP",
    "TYR",
    # N-terminal variants that gmx/pdb2gmx sometimes emits
    "NARG", "NASN", "NASP", "NCYS", "NCYM", "NCYX", "NGLN", "NGLU", "NGLH",
    "NHIS", "NHID", "NHIE", "NHIP", "NLEU", "NLYS", "NMET", "NPHE",
    "NPRO", "NSER", "NTRP", "NTYR",
}

# Residues with gamma-methylene (HG1/HG2 -> HG2/HG3). ILE-Cbeta HG is handled
# separately.
_HG_METHYLENE = {
    "ASN", "ASP", "CYS", "GLN", "GLU", "GLH", "LYS", "MET", "PRO", "ARG",
    "NASN", "NASP", "NCYS", "NGLN", "NGLU", "NGLH", "NLYS", "NMET", "NPRO",
    "NARG",
}

# Residues with delta-methylene HD (chain methylene, NOT ring HD which stays
# named HD1/HD2 in HIS/PHE/TYR/TRP).
_HD_METHYLENE = {
    "LYS", "PRO", "ARG",
    "NLYS", "NPRO", "NARG",
}

# Residues with epsilon-methylene HE (LYS side-chain).
_HE_METHYLENE = {"LYS", "NLYS"}

_ILE_VARIANTS = {"ILE", "NILE"}
_NPRO_VARIANTS = {"PRO"}  # only N-terminal Pro has H1/H2 on N

_DELETE_RULES = (
    # (atom_name, resname or None for any)
    ("OC1", None),  # C-terminal oxygens dropped; tleap re-adds O/OXT
    ("OC2", None),
    ("HG", "CYM"),  # thiolate anion has no thiol H
)


def _build_rename_map(resname: str) -> dict[str, str]:
    r = resname.upper().strip()
    out: dict[str, str] = {}

    if r in _HB_METHYLENE:
        out["HB1"] = "HB2"
        out["HB2"] = "HB3"
    if r in _HG_METHYLENE:
        out["HG1"] = "HG2"
        out["HG2"] = "HG3"
    if r in _HD_METHYLENE:
        out["HD1"] = "HD2"
        out["HD2"] = "HD3"
    if r in _HE_METHYLENE:
        out["HE1"] = "HE2"
        out["HE2"] = "HE3"

    if r in _ILE_VARIANTS:
        out["CD"] = "CD1"
        out["HD1"] = "HD11"
        out["HD2"] = "HD12"
        out["HD3"] = "HD13"
        # gmx sometimes emits HG11/HG12 for the CG methylene; AMBER wants
        # HG12/HG13.
        out["HG11"] = "HG12"
        out["HG12"] = "HG13"

    if r in _NPRO_VARIANTS:
        # Applied only to the N-terminal Pro (mid-chain Pro has no H1/H2 on N).
        # Applying globally is safe because mid-chain Pro won't have H1 to
        # begin with.
        out["H1"] = "H2"
        out["H2"] = "H3"

    return out


def _should_delete(atom_name: str, resname: str) -> bool:
    resname_u = resname.upper().strip()
    for a, r in _DELETE_RULES:
        if atom_name != a:
            continue
        if r is None or r == resname_u:
            return True
    return False


def _classify_his_by_h(atoms: set[str]) -> str:
    """HID/HIE/HIP from imidazole H presence."""
    has_hd1 = "HD1" in atoms
    has_he2 = "HE2" in atoms
    if has_hd1 and has_he2:
        return "HIP"
    if has_hd1:
        return "HID"
    if has_he2:
        return "HIE"
    # No labile H at all — leave as HIS (tleap will assign default).
    return "HIS"


def rename_atoms_in_pdb(pdb_path: str | Path, *, in_place: bool = True) -> Path:
    """Apply FRUTON -> AMBER atom-name normalization to a PDB file.

    All 6 transformation categories above are applied.  Deletions do not
    renumber remaining atom serials (tleap re-numbers on load).

    Returns the output path (same as input if in_place=True).
    """
    p = Path(pdb_path)
    if not p.is_file():
        raise FileNotFoundError(p)

    lines_in = p.read_text().splitlines()

    # First pass: for HIS residues, collect the set of atom names to classify.
    his_key_atoms: dict[tuple[str, int, str], set[str]] = {}  # (chain, resi, icode) -> set of atom names
    for line in lines_in:
        if not (line.startswith("ATOM") or line.startswith("HETATM")):
            continue
        resn = line[17:20].strip().upper()
        if resn not in ("HIS", "HID", "HIE", "HIP"):
            continue
        try:
            resi = int(line[22:26])
        except ValueError:
            continue
        chain = line[21]
        icode = line[26].strip()
        atom = line[12:16].strip()
        his_key_atoms.setdefault((chain, resi, icode), set()).add(atom)

    his_new_name: dict[tuple[str, int, str], str] = {
        k: _classify_his_by_h(atoms) for k, atoms in his_key_atoms.items()
    }

    # Second pass: rename + delete.
    lines_out: list[str] = []
    for line in lines_in:
        if not (line.startswith("ATOM") or line.startswith("HETATM")):
            lines_out.append(line)
            continue

        resn = line[17:20].strip().upper()
        atom = line[12:16].strip()

        # 4. Deletion rules.
        if _should_delete(atom, resn):
            continue

        # 1, 2, 6. Rename map.
        rename = _build_rename_map(resn)
        if atom in rename:
            new_atom = rename[atom]
            # PDB atom name field is 4-char, columns 13-16 (1-indexed).
            if len(new_atom) < 4:
                atom_field = f" {new_atom:<3s}"
            else:
                atom_field = new_atom[:4]
            line = line[:12] + atom_field + line[16:]

        # 5. HIS -> HID/HIE/HIP.
        if resn in ("HIS", "HID", "HIE", "HIP"):
            try:
                resi = int(line[22:26])
            except ValueError:
                lines_out.append(line)
                continue
            chain = line[21]
            icode = line[26].strip()
            new_resn = his_new_name.get((chain, resi, icode))
            if new_resn and new_resn != resn:
                line = line[:17] + f"{new_resn:>3s}" + line[20:]

        lines_out.append(line)

    if in_place:
        out_path = p
    else:
        out_path = p.with_name(p.stem + "_renamed" + p.suffix)
    out_path.write_text("\n".join(lines_out) + "\n")
    return out_path
