"""Convert protein termini to AMBER force-field naming conventions."""
from __future__ import annotations

import re
from dataclasses import dataclass
from pathlib import Path

from Bio.PDB import PDBIO, PDBParser
from Bio.PDB.Atom import Atom
from Bio.PDB.Residue import Residue
from Bio.Vector import Vector

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

_CAP_RESNAMES = {"ACE", "NME"}
_SKIP_RESNAMES = {"HOH", "WAT", "SOL"}
_STANDARD_AA = frozenset({
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS",
    "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP",
    "TYR", "VAL",
    # AMBER protonation-state variants
    "HID", "HIE", "HIP", "ASH", "GLH", "CYM", "CYX", "LYN",
    # Common non-standard residues that have backbone
    "MSE", "SEC", "SEP", "TPO", "PTR",
})


# ---------------------------------------------------------------------------
# Public helpers
# ---------------------------------------------------------------------------


def sanitize_variant_label(label: str) -> str:
    """Replace spaces/slashes/special chars with '_', strip, raise ValueError if empty."""
    cleaned = re.sub(r"[^A-Za-z0-9_\-]", "_", label.strip())
    cleaned = cleaned.strip("_")
    if not cleaned:
        raise ValueError(f"Invalid empty variant label derived from {label!r}")
    return cleaned


def build_default_terminus_output_path(pdb_directory: Path, pdb_id: str) -> Path:
    """Return the default AMBER terminus output path."""
    return Path(pdb_directory) / "components" / f"{pdb_id}_protein_amber_termini.pdb"


def build_variant_terminus_output_path(
    pdb_directory: Path, pdb_id: str, variant_label: str
) -> Path:
    """Return the variant-specific AMBER terminus output path."""
    safe = sanitize_variant_label(variant_label)
    return Path(pdb_directory) / "components" / f"{pdb_id}_protein_amber_termini_{safe}.pdb"


# ---------------------------------------------------------------------------
# Result dataclass
# ---------------------------------------------------------------------------


@dataclass
class TerminusConversionResult:
    chain_id: str
    first_resseq: int
    last_resseq: int
    first_old_resname: str
    first_new_resname: str
    last_old_resname: str
    last_new_resname: str
    single_residue_chain: bool


# ---------------------------------------------------------------------------
# Private helpers
# ---------------------------------------------------------------------------


def _is_polymer_residue(residue: Residue) -> bool:
    """Return True if the residue is a standard polymer residue (not HOH/cap/HETATM)."""
    hetflag = residue.get_id()[0]
    resname = residue.get_resname().strip().upper()
    if hetflag not in {" ", ""}:
        return False
    if resname in _SKIP_RESNAMES:
        return False
    if resname in _CAP_RESNAMES:
        return False
    return True


def _ensure_n_terminus_hydrogens(residue: Residue) -> None:
    """Ensure H1, H2, H3 atoms exist on the N-terminus residue."""
    existing = {atom.get_name().strip() for atom in residue.get_atoms()}
    # Get N atom position as reference.
    n_atom = None
    for atom in residue.get_atoms():
        if atom.get_name().strip() == "N":
            n_atom = atom
            break

    if n_atom is None:
        ref_coord = [0.0, 0.0, 0.0]
    else:
        ref_coord = list(n_atom.get_vector())

    serial = 9000
    for h_name in ("H1", "H2", "H3"):
        if h_name not in existing:
            offset = {"H1": [0.3, 0.3, 0.0], "H2": [-0.3, 0.3, 0.0], "H3": [0.0, 0.3, 0.3]}[h_name]
            coord = [ref_coord[i] + offset[i] for i in range(3)]
            new_atom = Atom(
                name=h_name,
                coord=coord,
                bfactor=20.0,
                occupancy=1.0,
                altloc=" ",
                fullname=f" {h_name}",
                serial_number=serial,
                element="H",
            )
            residue.add(new_atom)
            serial += 1


def _ensure_c_terminus_oxt(residue: Residue) -> None:
    """Ensure OXT atom exists on the C-terminus residue."""
    existing = {atom.get_name().strip() for atom in residue.get_atoms()}
    if "OXT" in existing:
        return

    # Get C atom position as reference.
    c_atom = None
    for atom in residue.get_atoms():
        if atom.get_name().strip() == "C":
            c_atom = atom
            break

    if c_atom is None:
        coord = [0.0, 1.25, 0.0]
    else:
        ref = list(c_atom.get_vector())
        coord = [ref[0], ref[1] + 1.25, ref[2]]

    new_atom = Atom(
        name="OXT",
        coord=coord,
        bfactor=20.0,
        occupancy=1.0,
        altloc=" ",
        fullname=" OXT",
        serial_number=9999,
        element="O",
    )
    residue.add(new_atom)


# ---------------------------------------------------------------------------
# Core conversion
# ---------------------------------------------------------------------------


def convert_protein_termini_to_amber(
    input_pdb_path: Path,
    output_pdb_path: Path,
) -> list[TerminusConversionResult]:
    """Convert protein chain termini to AMBER naming conventions.

    The N-prefix (NALA etc.) and C-prefix (CSER etc.) are recorded in the
    returned TerminusConversionResult objects, but the PDB file uses the
    original three-letter residue names so BioPython can parse it normally.

    Parameters
    ----------
    input_pdb_path
        Source PDB file (must exist).
    output_pdb_path
        Destination PDB file.

    Returns
    -------
    list[TerminusConversionResult]
        One entry per processed chain.

    Raises
    ------
    FileNotFoundError
        If *input_pdb_path* does not exist.
    """
    input_pdb_path = Path(input_pdb_path)
    output_pdb_path = Path(output_pdb_path)

    if not input_pdb_path.exists():
        raise FileNotFoundError(f"Input PDB not found: {input_pdb_path}")

    parser = PDBParser(QUIET=True)
    structure = parser.get_structure(input_pdb_path.stem, str(input_pdb_path))

    results: list[TerminusConversionResult] = []

    for model in structure.get_models():
        for chain in model.get_chains():
            chain_id = chain.get_id()

            # Find polymer residues (skip ACE/NME caps and water).
            polymer_residues = [r for r in chain.get_residues() if _is_polymer_residue(r)]

            if not polymer_residues:
                continue

            first_res = polymer_residues[0]
            last_res = polymer_residues[-1]

            first_old = first_res.get_resname().strip()
            last_old = last_res.get_resname().strip()

            single = len(polymer_residues) == 1

            if single:
                first_new = first_old
                last_new = last_old
                # Ensure terminal atoms even for single residue chains.
                _ensure_n_terminus_hydrogens(first_res)
                _ensure_c_terminus_oxt(last_res)
            else:
                first_new = f"N{first_old}"
                last_new = f"C{last_old}"
                # Ensure terminal atoms.
                _ensure_n_terminus_hydrogens(first_res)
                _ensure_c_terminus_oxt(last_res)
                # Note: we do NOT rename the residue in the PDB structure;
                # NALA/CSER naming is only reported in the result object.

            results.append(
                TerminusConversionResult(
                    chain_id=chain_id,
                    first_resseq=first_res.get_id()[1],
                    last_resseq=last_res.get_id()[1],
                    first_old_resname=first_old,
                    first_new_resname=first_new,
                    last_old_resname=last_old,
                    last_new_resname=last_new,
                    single_residue_chain=single,
                )
            )

    output_pdb_path.parent.mkdir(parents=True, exist_ok=True)
    io = PDBIO()
    io.set_structure(structure)
    io.save(str(output_pdb_path))

    return results


def convert_variant_protein_termini_to_amber(
    pdb_id: str,
    pdb_directory: Path,
    variant_label: str,
    input_pdb_path: Path,
) -> dict:
    """Convert termini for a variant structure and return summary dict.

    Returns
    -------
    dict
        Keys: amber_termini_success, amber_termini_variant,
        amber_termini_input_path, amber_termini_output_path,
        amber_termini_chain_count.
    """
    pdb_directory = Path(pdb_directory)
    input_pdb_path = Path(input_pdb_path)
    output_pdb_path = build_variant_terminus_output_path(pdb_directory, pdb_id, variant_label)

    try:
        chain_results = convert_protein_termini_to_amber(
            input_pdb_path=input_pdb_path,
            output_pdb_path=output_pdb_path,
        )
        success = True
        chain_count = len(chain_results)
    except Exception:
        success = False
        chain_count = 0

    return {
        "amber_termini_success": success,
        "amber_termini_variant": variant_label,
        "amber_termini_input_path": str(input_pdb_path.resolve()),
        "amber_termini_output_path": str(output_pdb_path.resolve()),
        "amber_termini_chain_count": chain_count,
    }
