from __future__ import annotations

"""Metal-center parametrization preparation module.

Purpose
-------
Prepare per-site MCPB scaffold folders for every transition-metal atom found
in a FRUTON protein directory.  For each site the module:

1. Detects transition-metal atoms using the same deterministic logic as
   :mod:`stack_protein_preparation.metalls_check`.
2. Creates a numbered site folder under ``metall_params/`` with four
   sub-directories: ``00_input/``, ``01_hydrogen_cleanup/``, ``02_contacts/``,
   and ``03_mcpb/``.
3. Runs :func:`~stack_protein_preparation.metal_hydrogen_cleanup.remove_metal_coordinating_hydrogens`
   to strip spurious hydrogens from the environment before parametrization.
4. Writes a deterministic ``deterministic_metal_contacts.tsv`` (Python only,
   no Chimera dependency).
5. Writes an MCPB.py input scaffold with ``TODO`` placeholders for the user to
   fill in before running Gaussian.
6. Optionally writes and runs a ChimeraX inspection script (never used for
   decisions; for visual inspection only).
7. Writes a top-level ``manifest.json`` and ``all_sites_summary.tsv``.

Public API
----------
:func:`run_metal_parametrization_for_protein_dir`
"""

import json
import os
import shutil
import subprocess
from pathlib import Path
from typing import Any

from stack_protein_preparation.metal_hydrogen_cleanup import (
    STANDARD_AMINO_ACID_NAMES as _STD_AA_NAMES,
    remove_metal_coordinating_hydrogens,
)
from stack_protein_preparation.metalls_check import (
    DONOR_ELEMENTS,
    METAL_RESIDUE_IDENTITY,
    TRANSITION_ION_GEOMETRY_PREFERENCES,
    TRANSITION_METAL_ELEMENTS,
    classify_coordination_geometry,
    distance_between_atoms,
    find_metal_contacts,
    is_true_metal_atom,
    read_pdb_atom_records,
)
from stack_protein_preparation.pdb_components import WATER_NAMES as _WATER_NAMES

# Spin multiplicity is unambiguous only for closed-shell d10 ions.
# All others depend on oxidation state and ligand field — left as TODO in the scaffold.
_UNAMBIGUOUS_SPIN: dict[str, int] = {
    "ZN": 1,  # d10 Zn2+: singlet (closed shell, no unpaired electrons)
    "CD": 1,  # d10 Cd2+: singlet
    "HG": 1,  # d10 Hg2+: singlet
}

# ---------------------------------------------------------------------------
# ChimeraX executable candidates (inspection only)
# ---------------------------------------------------------------------------

DEFAULT_CHIMERAX_EXECUTABLE_CANDIDATES = [
    "chimerax",
    "ChimeraX",
    "/usr/bin/chimerax",
    "/usr/local/bin/chimerax",
    "/opt/UCSF/ChimeraX/bin/ChimeraX",
    "/Applications/ChimeraX.app/Contents/MacOS/ChimeraX",
]

# Legacy UCSF Chimera fallbacks (no longer used for decisions, kept for discovery)
DEFAULT_CHIMERA_EXECUTABLE_CANDIDATES = [
    "chimera",
    "/usr/bin/chimera",
    "/home/grheco/.local/bin/chimera",
    "/opt/UCSF/Chimera64-1.17/bin/chimera",
    "/opt/UCSF/Chimera64-1.16/bin/chimera",
]


# ---------------------------------------------------------------------------
# Input selection
# ---------------------------------------------------------------------------

def _choose_best_protein_input(components_dir: Path, pdb_id: str) -> Path | None:
    """Return the best available protonated protein structure for this site.

    Prefers the fully assembled prepared variant (prepared/<variant>/<PDBID>.pdb)
    over per-fragment component files because it contains all chain segments with
    unique chain IDs and continuous serial numbers — required for correct MCPB
    parametrization of gap-containing proteins.  Falls back to per-fragment
    component files when no prepared variant is present.
    """
    protein_dir = components_dir.parent
    prepared_variants = ["single", "best_complete", "gaps", "large_gap_complete", "complete"]
    for variant in prepared_variants:
        prepared_path = protein_dir / "prepared" / variant / f"{pdb_id}.pdb"
        if prepared_path.is_file() and prepared_path.stat().st_size > 0:
            return prepared_path

    candidates = [
        components_dir / f"{pdb_id}_proteinH_single.pdb",
        components_dir / f"{pdb_id}_proteinH_best_complete.pdb",
        components_dir / f"{pdb_id}_proteinH_gaps.pdb",
        components_dir / f"{pdb_id}_protein_monomer.pdb",
        components_dir / f"{pdb_id}_protein.pdb",
        components_dir / f"{pdb_id}_protein_final.pdb",
        components_dir / f"{pdb_id}_protein_as_Amber.pdb",
        components_dir / f"{pdb_id}_proteinH.pdb",
    ]
    for path in candidates:
        if path.is_file():
            return path
    return None


def _get_optional_component_paths(
    components_dir: Path, pdb_id: str
) -> dict[str, Path | None]:
    """Return optional environment component paths if the files are present."""
    water_path = components_dir / f"{pdb_id}_water.pdb"
    ligand_path = components_dir / f"{pdb_id}_ligand.pdb"
    metal_path = components_dir / f"{pdb_id}_metals.pdb"
    return {
        "water": water_path if water_path.is_file() else None,
        "ligand": ligand_path if ligand_path.is_file() else None,
        "metal": metal_path if metal_path.is_file() else None,
    }


def _has_metal_file(components_dir: Path, pdb_id: str) -> bool:
    metal_path = components_dir / f"{pdb_id}_metals.pdb"
    return metal_path.is_file() and metal_path.stat().st_size > 0


# ---------------------------------------------------------------------------
# PDB file helpers
# ---------------------------------------------------------------------------

def _read_pdb_payload_lines(pdb_path: Path) -> list[str]:
    """Read all lines except END / ENDMDL terminators."""
    payload: list[str] = []
    with pdb_path.open("r", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.strip() in {"END", "ENDMDL"}:
                continue
            payload.append(line)
    return payload


_METAL_RESNAMES: frozenset[str] = frozenset(METAL_RESIDUE_IDENTITY.keys())

# Residue names that are never "non-standard organic ligands" in the MCPB
# context — they are either standard amino acids, water molecules, or metal
# ions, all of which are handled by MCPB natively.
_NON_NAA_RESNAMES: frozenset[str] = (
    frozenset(_STD_AA_NAMES) | frozenset(_WATER_NAMES) | _METAL_RESNAMES
)

# ---------------------------------------------------------------------------
# GROMACS → AMBER atom name mapping
# ---------------------------------------------------------------------------
# GROMACS (amber99sb-ildn) numbers methylene hydrogens from 1 (HB1/HB2),
# but AMBER (ff99SB) numbers them from 2 (HB2/HB3).  The mismatch causes
# MCPB.py to raise KeyError during the residue libdict lookup.
# ILE additionally uses CD/HD[1-3] in GROMACS vs CD1/HD1[1-3] in AMBER.
_GROMACS_TO_AMBER_ATOM: dict[str, dict[str, str]] = {
    "ARG": {"HB1": "HB2", "HB2": "HB3", "HG1": "HG2", "HG2": "HG3",
            "HD1": "HD2", "HD2": "HD3"},
    "ASN": {"HB1": "HB2", "HB2": "HB3"},
    "ASP": {"HB1": "HB2", "HB2": "HB3"},
    "CYS": {"HB1": "HB2", "HB2": "HB3"},
    "GLN": {"HB1": "HB2", "HB2": "HB3", "HG1": "HG2", "HG2": "HG3"},
    "GLU": {"HB1": "HB2", "HB2": "HB3", "HG1": "HG2", "HG2": "HG3"},
    "GLY": {"HA1": "HA2", "HA2": "HA3"},
    "HID": {"HB1": "HB2", "HB2": "HB3"},
    "HIE": {"HB1": "HB2", "HB2": "HB3"},
    "HIP": {"HB1": "HB2", "HB2": "HB3"},
    "HIS": {"HB1": "HB2", "HB2": "HB3"},
    "ILE": {"HG11": "HG12", "HG12": "HG13",
            "HD1": "HD11", "HD2": "HD12", "HD3": "HD13", "CD": "CD1"},
    "LEU": {"HB1": "HB2", "HB2": "HB3"},
    "LYS": {"HB1": "HB2", "HB2": "HB3", "HG1": "HG2", "HG2": "HG3",
            "HD1": "HD2", "HD2": "HD3", "HE1": "HE2", "HE2": "HE3"},
    "MET": {"HB1": "HB2", "HB2": "HB3", "HG1": "HG2", "HG2": "HG3"},
    "PHE": {"HB1": "HB2", "HB2": "HB3"},
    "PRO": {"HB1": "HB2", "HB2": "HB3", "HG1": "HG2", "HG2": "HG3",
            "HD1": "HD2", "HD2": "HD3"},
    "SER": {"HB1": "HB2", "HB2": "HB3"},
    "TRP": {"HB1": "HB2", "HB2": "HB3"},
    "TYR": {"HB1": "HB2", "HB2": "HB3"},
}


def _find_naa_resnames_in_contacts(contacts: list[Any]) -> list[str]:
    """Return unique non-standard organic ligand residue names from *contacts*.

    A contact is a non-standard ligand if its *donor_residue_name* is not a
    standard amino acid, not water, and not a metal ion.  The returned list
    preserves the order of first appearance and has no duplicates.
    """
    seen: list[str] = []
    for c in contacts:
        rn = c.donor_residue_name.strip().upper()
        if rn not in _NON_NAA_RESNAMES and rn not in seen:
            seen.append(rn)
    return seen


def _extract_residue_pdb(
    mcpb_pdb: Path,
    resname: str,
    output_pdb: Path,
) -> bool:
    """Extract all ATOM/HETATM lines matching *resname* from *mcpb_pdb* into
    *output_pdb*.  Returns True if at least one line was written."""
    lines = [
        ln for ln in mcpb_pdb.read_text(encoding="utf-8", errors="replace").splitlines(keepends=True)
        if ln.startswith(("ATOM  ", "HETATM")) and ln[17:20].strip().upper() == resname.upper()
    ]
    if not lines:
        return False
    output_pdb.write_text("".join(lines) + "END\n", encoding="utf-8")
    return True


def _generate_naa_mol2(
    resname: str,
    pdb_path: Path,
    mol2_path: Path,
    net_charge: int = 0,
) -> bool:
    """Try to generate a GAFF2 mol2 for *resname* using antechamber.

    Returns True if antechamber produced the mol2 file, False otherwise
    (antechamber not found, failed, or no pdb to work from).
    """
    if not pdb_path.is_file() or pdb_path.stat().st_size == 0:
        return False
    antechamber = shutil.which("antechamber")
    if antechamber is None:
        return False
    try:
        result = subprocess.run(
            [
                antechamber,
                "-i", str(pdb_path),
                "-fi", "pdb",
                "-o", str(mol2_path),
                "-fo", "mol2",
                "-c", "bcc",
                "-nc", str(net_charge),
                "-at", "gaff2",
                "-s", "2",
                "-pf", "y",
            ],
            capture_output=True,
            text=True,
            timeout=300,
            cwd=str(mol2_path.parent),
        )
        return mol2_path.is_file() and mol2_path.stat().st_size > 0
    except Exception:
        return False


def _is_metal_hetatm_line(line: str) -> bool:
    """Return True if *line* is a HETATM record for a bare metal ion."""
    if not line.startswith("HETATM"):
        return False
    resname = line[17:20].strip().upper()
    return resname in _METAL_RESNAMES


def _write_combined_tmp_param_pdb(
    output_path: Path,
    protein_path: Path,
    water_path: Path | None,
    ligand_path: Path | None,
) -> dict[str, Any]:
    """Concatenate protein + optional water + optional ligand into one PDB.

    Metal HETATM records are stripped so the environment is a clean donor
    context for hydrogen cleanup and MCPB input.  This is important when
    *protein_path* is the fully assembled prepared variant (which already
    embeds the metal) rather than a component-only file.
    """
    source_paths: list[Path] = [protein_path]
    if water_path is not None:
        source_paths.append(water_path)
    if ligand_path is not None:
        source_paths.append(ligand_path)

    line_count = 0
    with output_path.open("w", encoding="utf-8") as out:
        for src in source_paths:
            lines = [
                ln for ln in _read_pdb_payload_lines(src)
                if not _is_metal_hetatm_line(ln)
            ]
            for line in lines:
                out.write(line)
                line_count += 1
            if lines and not lines[-1].startswith("TER"):
                out.write("TER\n")
                line_count += 1
        out.write("END\n")
        line_count += 1

    return {
        "source_paths": [str(p) for p in source_paths],
        "line_count_written": line_count,
    }


def _copy_metal_only_pdb(metal_input_path: Path, metal_only_output_path: Path) -> None:
    """Copy the metal component to a local file for this site."""
    text = metal_input_path.read_text(encoding="utf-8", errors="replace")
    metal_only_output_path.write_text(text, encoding="utf-8")


def _write_single_metal_pdb(
    metal_atom_raw_line: str,
    output_path: Path,
) -> None:
    """Write a minimal PDB containing only one metal atom."""
    with output_path.open("w", encoding="utf-8") as fh:
        fh.write(metal_atom_raw_line.rstrip("\n") + "\n")
        fh.write("END\n")


def _combined_mcpb_pdb(
    environment_pdb: Path,
    metal_only_pdb: Path,
    output_path: Path,
) -> None:
    """Concatenate the cleaned environment and the metal-only PDB for MCPB."""
    env_lines = _read_pdb_payload_lines(environment_pdb)
    metal_lines = _read_pdb_payload_lines(metal_only_pdb)
    with output_path.open("w", encoding="utf-8") as fh:
        for line in env_lines:
            fh.write(line)
        if env_lines and not env_lines[-1].startswith("TER"):
            fh.write("TER\n")
        for line in metal_lines:
            fh.write(line)
        fh.write("END\n")


def _renumber_pdb_atoms(pdb_path: Path) -> None:
    """Rewrite *pdb_path* in-place with sequential atom serial numbers starting
    from 1.

    Only ATOM and HETATM lines are renumbered; all other lines (REMARK, TER,
    CONECT, END, …) are written through unchanged.  This is needed when
    *pdb_path* is assembled from multiple source PDB files whose serial
    counters all start at 1, producing duplicate serials that confuse
    MCPB.py -s 1.
    """
    lines = pdb_path.read_text(encoding="utf-8", errors="replace").splitlines(keepends=True)
    serial = 1
    out_lines: list[str] = []
    for line in lines:
        if line.startswith(("ATOM  ", "HETATM")):
            # PDB columns 7-11 (0-based 6:11) are the serial number field
            out_lines.append(f"{line[:6]}{serial:5d}{line[11:]}")
            serial += 1
        else:
            out_lines.append(line)
    pdb_path.write_text("".join(out_lines), encoding="utf-8")


def _rename_gromacs_to_amber_atoms(pdb_path: Path) -> None:
    """Rename GROMACS (amber99sb-ildn) atom names to AMBER (ff99SB) in-place.

    GROMACS numbers methylene H's from 1 (HB1/HB2); MCPB.py's ff99SB library
    expects them from 2 (HB2/HB3).  Also fixes ILE CD→CD1 and HD[1-3]→HD1[1-3].
    Only ATOM/HETATM lines for residues listed in _GROMACS_TO_AMBER_ATOM are
    touched; all other lines pass through unchanged.
    """
    lines = pdb_path.read_text(encoding="utf-8", errors="replace").splitlines(keepends=True)
    out_lines: list[str] = []
    for line in lines:
        if line.startswith(("ATOM  ", "HETATM")):
            resname = line[17:20].strip()
            rename_map = _GROMACS_TO_AMBER_ATOM.get(resname)
            if rename_map:
                atom_name = line[12:16].strip()
                new_name = rename_map.get(atom_name)
                if new_name is not None:
                    # PDB atom-name field is 4 chars (cols 12-15, 0-indexed).
                    # 1-3 char names carry a leading space; 4-char names fill the field.
                    if len(new_name) < 4:
                        new_field = f" {new_name:<3s}"
                    else:
                        new_field = f"{new_name:<4s}"
                    line = line[:12] + new_field + line[16:]
        out_lines.append(line)
    pdb_path.write_text("".join(out_lines), encoding="utf-8")


def _rename_his_by_protonation(pdb_path: Path) -> None:
    """Rename HIS residues to HID/HIE/HIP in-place based on which H atoms are present.

    GROMACS outputs all histidines as ``HIS`` regardless of protonation state.
    MCPB.py's AMBER ff99SB chargedict expects ``HID`` (δ-protonated), ``HIE``
    (ε-protonated), or ``HIP`` (doubly protonated).

    Detection: for each HIS residue, check for HD1 (δ N-H) and HE2 (ε N-H).
    """
    lines = pdb_path.read_text(encoding="utf-8", errors="replace").splitlines(keepends=True)

    # Collect per-residue H atom sets: key = (chain, resseq, icode)
    his_atoms: dict[tuple[str, str, str], set[str]] = {}
    for line in lines:
        if line.startswith(("ATOM  ", "HETATM")) and line[17:20].strip() == "HIS":
            key = (line[21].strip(), line[22:26].strip(), line[26].strip())
            atom = line[12:16].strip()
            his_atoms.setdefault(key, set()).add(atom)

    # Determine new residue name for each HIS key
    his_rename: dict[tuple[str, str, str], str] = {}
    for key, atoms in his_atoms.items():
        has_hd1 = "HD1" in atoms
        has_he2 = "HE2" in atoms
        if has_hd1 and has_he2:
            his_rename[key] = "HIP"
        elif has_hd1:
            his_rename[key] = "HID"
        else:
            his_rename[key] = "HIE"  # default: ε-protonated (most common)

    if not his_rename:
        return

    out_lines: list[str] = []
    for line in lines:
        if line.startswith(("ATOM  ", "HETATM")) and line[17:20].strip() == "HIS":
            key = (line[21].strip(), line[22:26].strip(), line[26].strip())
            new_rn = his_rename.get(key, "HIE")
            line = line[:17] + f"{new_rn:<3s}" + line[20:]
        out_lines.append(line)
    pdb_path.write_text("".join(out_lines), encoding="utf-8")


def _generate_water_mol2(mol2_path: Path, water_pdb: Path) -> bool:
    """Generate a mol2 for a coordinating water using antechamber (GAFF2/BCC).

    Extracts the first HOH record from *water_pdb* and runs antechamber.
    Returns True on success.
    """
    if not water_pdb.is_file() or water_pdb.stat().st_size == 0:
        return False
    antechamber = shutil.which("antechamber")
    if antechamber is None:
        return False
    try:
        result = subprocess.run(
            [
                antechamber,
                "-i", str(water_pdb),
                "-fi", "pdb",
                "-o", str(mol2_path),
                "-fo", "mol2",
                "-c", "bcc",
                "-nc", "0",
                "-at", "gaff2",
                "-s", "2",
                "-pf", "y",
            ],
            capture_output=True,
            text=True,
            timeout=120,
            cwd=str(mol2_path.parent),
        )
        return mol2_path.is_file() and mol2_path.stat().st_size > 0
    except Exception:
        return False


# ---------------------------------------------------------------------------
# Chimera executable resolution
# ---------------------------------------------------------------------------

def _resolve_chimera_executable(chimera_executable: str | None = None) -> str | None:
    """Resolve a ChimeraX (or legacy Chimera) executable.

    Priority:
    1. Explicit argument
    2. CHIMERAX_EXECUTABLE / CHIMERA_EXECUTABLE environment variables
    3. Candidate lists (ChimeraX preferred)
    """
    def _find(name: str) -> str | None:
        resolved = shutil.which(name)
        if resolved:
            return resolved
        path = Path(name)
        if path.is_file():
            return str(path)
        return None

    if chimera_executable:
        found = _find(chimera_executable)
        if found:
            return found

    for env_var in ("CHIMERAX_EXECUTABLE", "CHIMERA_EXECUTABLE"):
        env_val = os.environ.get(env_var)
        if env_val:
            found = _find(env_val)
            if found:
                return found

    for candidate in DEFAULT_CHIMERAX_EXECUTABLE_CANDIDATES:
        found = _find(candidate)
        if found:
            return found

    for candidate in DEFAULT_CHIMERA_EXECUTABLE_CANDIDATES:
        found = _find(candidate)
        if found:
            return found

    return None


# ---------------------------------------------------------------------------
# ChimeraX inspection script (02_contacts)
# ---------------------------------------------------------------------------

def _write_chimerax_cxc_script(
    script_path: Path,
    analysis_pdb_path: Path,
    contact_cutoff_angstrom: float,
    pdb_id: str,
    site_folder_name: str,
) -> None:
    """Write a ChimeraX ``.cxc`` inspection script for one metal site."""
    atom_spec_metals = "@" + ",".join(sorted(TRANSITION_METAL_ELEMENTS))
    script_lines = [
        f"open {analysis_pdb_path}",
        "delete solvent false",
        (
            f"contacts {atom_spec_metals} restrict @N*,O*,S*,SE* "
            f"distance {contact_cutoff_angstrom:.2f} reveal true "
            f"name metal_contacts color gold"
        ),
        f"# site: {site_folder_name}  pdb_id: {pdb_id}",
        "exit",
        "",
    ]
    script_path.write_text("\n".join(script_lines), encoding="utf-8")


def _run_chimerax_script(
    script_path: Path,
    log_path: Path,
    chimera_executable: str,
) -> str:
    """Run the ChimeraX script and return a status string."""
    try:
        completed = subprocess.run(
            [chimera_executable, "--nogui", str(script_path)],
            check=False,
            capture_output=True,
            text=True,
            timeout=180,
        )
    except Exception as exc:
        log_path.write_text(f"ChimeraX execution failed: {exc!r}\n", encoding="utf-8")
        return "failed"

    log_path.write_text(completed.stdout + "\n" + completed.stderr, encoding="utf-8")
    return "success" if completed.returncode == 0 else f"returncode_{completed.returncode}"


# ---------------------------------------------------------------------------
# Deterministic contacts TSV (02_contacts)
# ---------------------------------------------------------------------------

def _write_deterministic_contacts_tsv(
    output_path: Path,
    metal_atom: Any,
    all_records: list[Any],
    contact_cutoff_angstrom: float,
) -> None:
    """Write metal donor contacts computed purely in Python."""
    contacts = find_metal_contacts(
        metal_atom=metal_atom,
        atom_records=all_records,
        contact_cutoff_angstrom=contact_cutoff_angstrom,
    )
    contacts_sorted = sorted(contacts, key=lambda c: c.distance_angstrom)

    header = "\t".join([
        "metal_element", "metal_chain", "metal_resseq",
        "donor_element", "donor_residue_name", "donor_chain",
        "donor_resseq", "donor_atom_name", "distance_angstrom",
    ])
    rows = []
    for c in contacts_sorted:
        rows.append("\t".join([
            c.metal_element,
            metal_atom.chain_id,
            str(metal_atom.residue_number),
            c.donor_element,
            c.donor_residue_name,
            c.donor_chain_id,
            str(c.donor_residue_number),
            c.donor_atom_name.strip(),
            f"{c.distance_angstrom:.3f}",
        ]))

    output_path.write_text(header + "\n" + "\n".join(rows) + "\n", encoding="utf-8")


# ---------------------------------------------------------------------------
# MCPB scaffold writers (03_mcpb)
# ---------------------------------------------------------------------------

def _write_mcpb_in(
    output_path: Path,
    pdb_id: str,
    element: str,
    chain: str,
    resseq: int,
    max_contact_dist: float,
    ion_serial: int | None,
    ion_mol2file: str,
    formal_charge: int | None,
    spin_guess: int | None,
    naa_mol2files: list[str] | None = None,
) -> None:
    """Write a first-pass MCPB.py ``.in`` scaffold.

    Fields that FRUTON can resolve automatically (serial number, mol2 filename,
    metal formal charge, spin for d10 ions) are filled in directly.  Fields
    that require manual chemistry input (spin for open-shell metals, ligand
    partial charges) are written as ``TODO`` with an explanatory comment.
    """
    group_name = f"{pdb_id}_{element}_{chain}_{resseq}"
    cut_off = max_contact_dist + 0.5

    # ion_ids — resolved from actual serial in merged PDB
    if ion_serial is not None:
        ion_ids_line = f"ion_ids {ion_serial}"
        ion_ids_comment = f"# Serial of {element} in {pdb_id}_mcpb.pdb (auto-resolved by FRUTON)."
    else:
        ion_ids_line = "ion_ids TODO_metal_serial_in_mcpb_pdb"
        ion_ids_comment = f"# TODO: Serial number of {element} in {pdb_id}_mcpb.pdb."

    # ion_mol2files — generated by FRUTON when charge is known
    if not ion_mol2file.startswith("TODO"):
        mol2_comment = "# mol2 file for the metal ion (auto-generated by FRUTON from known formal charge)."
    else:
        mol2_comment = "# TODO: Provide mol2 file for the metal ion and any non-standard ligands."

    # charge_m_sm / charge_m_lm
    if formal_charge is not None and spin_guess is not None:
        charge_sm_line = f"charge_m_sm {formal_charge} {spin_guess}"
        charge_lm_line = f"charge_m_lm {formal_charge} {spin_guess}"
        charge_comment = (
            f"# Metal formal charge {formal_charge}+ and spin {spin_guess} "
            f"(unambiguous for {element}, auto-filled).\n"
            "# NOTE: Add net charge of any ionised ligand residues in the QM region."
        )
    elif formal_charge is not None:
        charge_sm_line = f"charge_m_sm {formal_charge} TODO_spin"
        charge_lm_line = f"charge_m_lm {formal_charge} TODO_spin"
        charge_comment = (
            f"# Metal formal charge {formal_charge}+ auto-filled.\n"
            f"# TODO: Set spin multiplicity — depends on oxidation state and\n"
            f"#       ligand-field splitting of {element}. "
            "1=singlet, 2=doublet, 3=triplet, …\n"
            "# NOTE: Add net charge of any ionised ligand residues in the QM region."
        )
    else:
        charge_sm_line = "charge_m_sm TODO_charge TODO_spin"
        charge_lm_line = "charge_m_lm TODO_charge TODO_spin"
        charge_comment = (
            "# TODO: Set total charge and spin multiplicity for the QM region.\n"
            "# Spin: 1=singlet, 2=doublet, 3=triplet, …"
        )

    # naa_mol2files — non-standard organic ligands in the coordination shell
    naa_block = ""
    if naa_mol2files:
        naa_list = " ".join(naa_mol2files)
        naa_block = (
            "\n# mol2 file(s) for non-standard organic ligand(s) in the QM region.\n"
            f"naa_mol2files {naa_list}\n"
        )

    content = f"""\
# MCPB.py input file — first-pass scaffold generated by FRUTON
# Review carefully before running. Remaining TODOs are marked explicitly.

original_pdb {pdb_id}_mcpb.pdb
group_name {group_name}

# cut_off derived from max donor-metal distance + 0.5 Å buffer. Verify.
cut_off {cut_off:.1f}

{ion_ids_comment}
{ion_ids_line}

{mol2_comment}
ion_mol2files {ion_mol2file}
{naa_block}
# 1 = optimise only H positions in the large model (recommended).
large_opt 1

force_field ff99SB

# Verify Gaussian version on your system (g09 or g16).
software_version g16

water_model tip3p

{charge_comment}
{charge_sm_line}
{charge_lm_line}
"""
    output_path.write_text(content, encoding="utf-8")


def _write_commands_sh(
    output_path: Path,
    pdb_id: str,
    group_name: str,
    mol2_filename: str,
    naa_mol2files: list[str] | None = None,
) -> None:
    """Write the step-1 script: creates step01_gen_inputs/, runs MCPB.py -s 1,
    then stages .com files into step02_gaussian/."""
    # Build copy line for naa mol2 files (if any)
    naa_copy_lines = ""
    naa_note = ""
    if naa_mol2files:
        naa_files_str = " ".join(naa_mol2files)
        naa_copy_lines = f"cp {naa_files_str} step01_gen_inputs/\n"
        naa_note = (
            "\n# NOTE: Non-standard ligand mol2 file(s) must exist before running:\n"
            + "".join(f"#   {f}\n" for f in naa_mol2files)
            + "# Generate with antechamber if needed (see comments in the .in file).\n"
        )

    content = f"""\
#!/usr/bin/env bash
# =============================================================================
# MCPB.py — Step 1  |  {pdb_id}  |  site {group_name}
# =============================================================================
# Creates step01_gen_inputs/, copies input files, runs MCPB.py -s 1.
# Stages the three Gaussian .com files into step02_gaussian/.
#
# After this script:
#   1. Review .com files in step02_gaussian/.
#   2. Submit Gaussian jobs:  bash submit_gaussian.sh
#   3. Once all Gaussian jobs complete:  bash run_after_gaussian.sh
# =============================================================================
set -euo pipefail
{naa_note}
# -- Create step01 directory and copy all required inputs --------------------
mkdir -p step01_gen_inputs
cp {pdb_id}.in {pdb_id}_mcpb.pdb {mol2_filename} step01_gen_inputs/
{naa_copy_lines}
# -- Run MCPB.py step 1 (generates Gaussian input files and model PDBs) ------
cd step01_gen_inputs
MCPB.py -i {pdb_id}.in -s 1
cd ..

# -- Stage Gaussian input files into step02_gaussian/ ------------------------
mkdir -p step02_gaussian
cp step01_gen_inputs/{group_name}_small_opt.com step02_gaussian/
cp step01_gen_inputs/{group_name}_small_fc.com  step02_gaussian/
cp step01_gen_inputs/{group_name}_large_mk.com  step02_gaussian/

echo ""
echo "Step 1 done.  Gaussian inputs staged in step02_gaussian/:"
echo "  {group_name}_small_opt.com   geometry optimisation (small model)"
echo "  {group_name}_small_fc.com    force constants / Hessian (small model)"
echo "  {group_name}_large_mk.com    RESP charges (large model)"
echo ""
echo "Review each .com file, then run:  bash submit_gaussian.sh"
"""
    output_path.write_text(content, encoding="utf-8")
    output_path.chmod(0o755)


def _write_submit_gaussian_sh(output_path: Path, group_name: str) -> None:
    """Write a CESGA-compatible MCPB Gaussian submission script.

    MCPB requires three Gaussian jobs in a specific order:
      - small_opt  (geometry optimisation → writes .chk)
      - large_mk   (MK ESP, independent of small_opt)
      - small_fc   (frequency/force-constants → reads small_opt.chk)

    small_opt and large_mk are submitted immediately; small_fc is submitted
    with ``--dependency=afterok:<small_opt_job_id>`` so SLURM holds it until
    small_opt succeeds.  This avoids the "submit all, let one fail, resubmit"
    anti-pattern.
    """
    content = f"""\
#!/usr/bin/env bash
# =============================================================================
# Gaussian HPC submission — {group_name}
# =============================================================================
# Submits the three MCPB Gaussian jobs with the correct SLURM dependency:
#   small_opt  ──┐ (parallel)
#   large_mk     │
#   small_fc  ←──┘ afterok:small_opt  (needs small_opt.chk)
#
# Run ONLY after bash commands.sh has completed successfully.
# After all jobs complete:  bash run_after_gaussian.sh
# =============================================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${{BASH_SOURCE[0]}}")" && pwd)"
# Locate the central CESGA submission script.
# Priority: $FRUTON_SCRIPTS env var (set this when working from $LUSTRE
# where the git repo is not present), then git-based discovery, then
# directory walk-up.
if [[ -n "${{FRUTON_SCRIPTS:-}}" && -f "${{FRUTON_SCRIPTS}}/submit_gaussian.sh" ]]; then
    CENTRAL="${{FRUTON_SCRIPTS}}/submit_gaussian.sh"
elif REPO_ROOT="$(git -C "$SCRIPT_DIR" rev-parse --show-toplevel 2>/dev/null)"; then
    CENTRAL="${{REPO_ROOT}}/scripts/CESGA_SLURM/submit_gaussian.sh"
else
    CENTRAL=""
    DIR="$SCRIPT_DIR"
    for _ in 1 2 3 4 5 6 7 8; do
        if [[ -f "${{DIR}}/scripts/CESGA_SLURM/submit_gaussian.sh" ]]; then
            CENTRAL="${{DIR}}/scripts/CESGA_SLURM/submit_gaussian.sh"
            break
        fi
        DIR="$(dirname "$DIR")"
    done
fi
[[ -n "${{CENTRAL:-}}" && -f "$CENTRAL" ]] || {{ echo "ERROR: central submit_gaussian.sh not found. Set FRUTON_SCRIPTS=/path/to/scripts/CESGA_SLURM" >&2; exit 2; }}

GAUSS_DIR="$SCRIPT_DIR/step02_gaussian"
SLURM_OPTS=(--partition short --time 06:00:00 --mem 4G --cpus 32 --gpus 1)

# 1. small_opt — geometry optimisation (produces .chk needed by small_fc)
OPT_OUT=$(bash "$CENTRAL" "$GAUSS_DIR/{group_name}_small_opt.com" "${{SLURM_OPTS[@]}}" "$@")
echo "$OPT_OUT"
OPT_JOB=$(echo "$OPT_OUT" | grep -oP '(?<=Submitted batch job )\\d+' || true)
if [[ -z "$OPT_JOB" ]]; then
    echo "WARNING: could not parse small_opt job ID — small_fc will be submitted without dependency" >&2
fi

# 2. large_mk — ESP calculation, independent of small_opt
bash "$CENTRAL" "$GAUSS_DIR/{group_name}_large_mk.com" "${{SLURM_OPTS[@]}}" "$@"

# 3. small_fc — force constants, depends on small_opt.chk
if [[ -n "$OPT_JOB" ]]; then
    bash "$CENTRAL" "$GAUSS_DIR/{group_name}_small_fc.com" "${{SLURM_OPTS[@]}}" \
        --dependency "afterok:$OPT_JOB" "$@"
else
    bash "$CENTRAL" "$GAUSS_DIR/{group_name}_small_fc.com" "${{SLURM_OPTS[@]}}" "$@"
fi

echo ""
echo "Submitted: small_opt (job $OPT_JOB), large_mk, small_fc (afterok:$OPT_JOB)"
echo "Once all complete:  bash run_after_gaussian.sh"
"""
    output_path.write_text(content, encoding="utf-8")
    output_path.chmod(0o755)


def _write_run_after_gaussian_sh(output_path: Path, pdb_id: str, group_name: str) -> None:
    """Write steps 2-4 script: creates step03_amber_params/, merges step01 +
    step02 outputs, runs MCPB.py -s 2,3,4 and tleap."""
    content = f"""\
#!/usr/bin/env bash
# =============================================================================
# MCPB.py — Steps 2-4  |  {pdb_id}  |  site {group_name}
# =============================================================================
# Creates step03_amber_params/, assembles all required inputs, then runs
# MCPB.py -s 2, -s 3, -s 4, and tleap.
#
# Run ONLY after all three Gaussian jobs in step02_gaussian/ have completed.
#
# Required files in step02_gaussian/ after Gaussian:
#   {group_name}_small_opt.fchk   (from: formchk _small_opt.chk)
#   {group_name}_small_fc.log
#   {group_name}_large_mk.log
#
# Final outputs in step03_amber_params/:
#   {group_name}.frcmod            AMBER bond/angle force constants
#   {group_name}_addmed.mol2       AMBER partial charges (RESP)
#   {group_name}_tleap.in / {group_name}_tleap.out
# =============================================================================
set -euo pipefail

# -- Create step03 directory and assemble all required inputs ----------------
mkdir -p step03_amber_params

# All step01 outputs (model PDBs, .in, mol2, MCPB intermediate files)
cp step01_gen_inputs/* step03_amber_params/

# Gaussian outputs from step02
cp step02_gaussian/{group_name}_small_opt.fchk step03_amber_params/
cp step02_gaussian/{group_name}_small_fc.log   step03_amber_params/
cp step02_gaussian/{group_name}_large_mk.log   step03_amber_params/

# -- Run MCPB.py steps 2-4 and tleap ----------------------------------------
cd step03_amber_params

# Step 2: Derive bond/angle force constants (Seminario method).
MCPB.py -i {pdb_id}.in -s 2

# Step 3: Fit RESP charges.
MCPB.py -i {pdb_id}.in -s 3

# Step 4: Generate LEaP topology input.
MCPB.py -i {pdb_id}.in -s 4
tleap -s -f {group_name}_tleap.in > {group_name}_tleap.out

echo ""
echo "Parametrization complete.  Key outputs in step03_amber_params/:"
echo "  {group_name}.frcmod"
echo "  {group_name}_addmed.mol2"
echo "  {group_name}_tleap.in / {group_name}_tleap.out"
"""
    output_path.write_text(content, encoding="utf-8")
    output_path.chmod(0o755)


def _write_metal_mol2(
    output_path: Path,
    element: str,
    formal_charge: int,
    x: float,
    y: float,
    z: float,
) -> None:
    """Write a minimal TRIPOS mol2 for a bare transition-metal ion.

    The atom type uses the AMBER/GAFF convention (element symbol in title case,
    e.g. ``Zn``, ``Fe``).  The formal charge is taken from
    :data:`METAL_RESIDUE_IDENTITY` and is written as a float in the charge
    column so that MCPB.py can read it directly.
    """
    atom_type = element.capitalize()
    charge_int = formal_charge
    # subst_name must match the PDB residue name (e.g. "ZN") so that MCPB.py's
    # chargedict lookup (keyed by residue name) finds the charge from this mol2.
    subst_name = element
    content = (
        "@<TRIPOS>MOLECULE\n"
        f"{element}\n"
        " 1 0 0 0 0\n"
        "SMALL\n"
        "NO_CHARGES\n"
        "\n"
        "\n"
        "@<TRIPOS>ATOM\n"
        f"      1 {element:<8s}{x:10.4f}{y:10.4f}{z:10.4f} "
        f"{atom_type:<8s}  1  {subst_name:<8s}{float(formal_charge):10.4f}\n"
        "@<TRIPOS>BOND\n"
    )
    output_path.write_text(content, encoding="utf-8")


def _find_metal_serial_in_pdb(pdb_path: Path, element: str) -> int | None:
    """Return the serial number of the first metal atom matching *element*."""
    try:
        for rec in read_pdb_atom_records(pdb_path, source_role="serial_lookup"):
            if is_true_metal_atom(rec) and rec.element == element:
                return rec.serial
    except Exception:
        pass
    return None


def _write_mcpb_readme(output_path: Path, pdb_id: str, group_name: str) -> None:
    """Write a README.md explaining the full MCPB multi-step workflow."""
    content = f"""\
# MCPB Parametrization — {pdb_id}  |  site {group_name}

## Overview

MCPB.py (Metal Center Parameter Builder) derives bonded-model AMBER parameters
for transition-metal centers using quantum-chemical calculations.

Reference: doi:10.1021/acs.jcim.5b00674

## Directory layout

```
03_mcpb/                         ← run all scripts from here
├── {pdb_id}.in
├── {pdb_id}_mcpb.pdb
├── *.mol2
├── commands.sh           Step 1
├── submit_gaussian.sh    Gaussian HPC submission
├── run_after_gaussian.sh Steps 2-4
│
├── step01_gen_inputs/    created by commands.sh
│   ├── {group_name}_small_opt.com
│   ├── {group_name}_small_fc.com
│   ├── {group_name}_large_mk.com
│   └── {group_name}_small.pdb / _large.pdb  (model PDBs)
│
├── step02_gaussian/      .com files staged here; Gaussian logs written here
│   ├── {group_name}_small_opt.log / .fchk
│   ├── {group_name}_small_fc.log
│   └── {group_name}_large_mk.log
│
└── step03_amber_params/  created by run_after_gaussian.sh; final outputs here
    ├── {group_name}.frcmod
    ├── {group_name}_addmed.mol2
    └── {group_name}_tleap.in / _tleap.out
```

## Workflow (3 scripts, 4 MCPB.py invocations)

```
bash commands.sh
  → mkdir step01_gen_inputs/, run MCPB.py -s 1
  → mkdir step02_gaussian/, copy .com files

bash submit_gaussian.sh
  → sbatch all three Gaussian jobs from step02_gaussian/

bash run_after_gaussian.sh
  → mkdir step03_amber_params/, merge step01 + step02 outputs
  → MCPB.py -s 2, -s 3, -s 4, tleap
```

## Before running Step 1

Review ``{pdb_id}.in``.  Fields auto-filled by FRUTON are marked in comments.
Fields that require manual chemistry input carry ``TODO`` placeholders:

- ``software_version`` — Gaussian version on your cluster (g09 or g16)
- ``charge_m_sm / charge_m_lm`` — total charge + spin multiplicity of the QM
  region; spin is unambiguous only for closed-shell d10 ions (Zn²⁺, Cu⁺)

## Gaussian

Gaussian is a licensed external program; FRUTON does not run it automatically.
``submit_gaussian.sh`` provides a SLURM template — adjust ``--cpus``,
``--gpus``, ``--mem``, and ``--time`` to match your cluster's policies.

All three Gaussian jobs are independent and can be submitted simultaneously.
Steps 2-4 require all three to finish before starting.
"""
    output_path.write_text(content, encoding="utf-8")


# ---------------------------------------------------------------------------
# Per-site orchestration
# ---------------------------------------------------------------------------

def _process_one_site(
    *,
    site_index: int,
    metal_atom: Any,
    all_records: list[Any],
    protein_dir: Path,
    pdb_id: str,
    protein_input: Path,
    water_input: Path | None,
    ligand_input: Path | None,
    metal_input: Path,
    metall_params_dir: Path,
    chimera_executable: str | None,
    contact_cutoff_angstrom: float,
) -> dict[str, Any]:
    """Build and populate one site folder; return a per-site result dict."""

    element = metal_atom.element
    chain = metal_atom.chain_id or "X"
    resseq = metal_atom.residue_number

    site_folder_name = f"site_{site_index:03d}_{element}_{chain}_{resseq}"
    site_dir = metall_params_dir / site_folder_name

    input_dir = site_dir / "00_input"
    cleanup_dir = site_dir / "01_hydrogen_cleanup"
    contacts_dir = site_dir / "02_contacts"
    mcpb_dir = site_dir / "03_mcpb"
    logs_dir = site_dir / "logs"

    for d in (input_dir, cleanup_dir, contacts_dir, mcpb_dir, logs_dir):
        d.mkdir(parents=True, exist_ok=True)

    site_result: dict[str, Any] = {
        "site_index": site_index,
        "site_folder": site_folder_name,
        "element": element,
        "chain_id": chain,
        "residue_number": resseq,
        "status": "failed",
        "hydrogen_cleanup_status": "failed",
        "removed_h_count": 0,
        "kept_h_count": 0,
        "manual_review_count": 0,
        "coordination_number": 0,
        "found_geometry": "",
        "geometry_ok": False,
        "mcpb_scaffold_generated": False,
        "mcpb_folder": str(mcpb_dir.relative_to(metall_params_dir)),
        "message": "",
    }

    # ------------------------------------------------------------------
    # 00_input — environment PDB (no metal) + metal_only.pdb
    # ------------------------------------------------------------------
    env_before_cleanup = input_dir / "environment_before_cleanup.pdb"
    metal_only_pdb = input_dir / "metal_only.pdb"
    prepared_variant_pdb = input_dir / "prepared_variant.pdb"

    _write_combined_tmp_param_pdb(
        output_path=env_before_cleanup,
        protein_path=protein_input,
        water_path=water_input,
        ligand_path=ligand_input,
    )

    # metal_only.pdb: write the raw line for this specific metal atom
    if metal_atom.raw_line:
        _write_single_metal_pdb(metal_atom.raw_line, metal_only_pdb)
    else:
        _copy_metal_only_pdb(metal_input, metal_only_pdb)

    # prepared_variant = env + metal together
    _combined_mcpb_pdb(env_before_cleanup, metal_only_pdb, prepared_variant_pdb)

    # ------------------------------------------------------------------
    # 01_hydrogen_cleanup
    # Pass prepared_variant_pdb (env + metal) so the metal is visible and
    # remove_metal_coordinating_hydrogens can actually find transition-metal
    # contacts.  The cleaned output is the full structure with bad H removed.
    # ------------------------------------------------------------------
    cleanup_result = remove_metal_coordinating_hydrogens(
        input_pdb_path=prepared_variant_pdb,
        output_dir=cleanup_dir,
        pdb_id=pdb_id,
        variant_label=site_folder_name,
        contact_cutoff_angstrom=contact_cutoff_angstrom,
    )

    site_result["hydrogen_cleanup_status"] = cleanup_result.status
    site_result["removed_h_count"] = cleanup_result.removed_count
    site_result["kept_h_count"] = cleanup_result.kept_count
    site_result["manual_review_count"] = cleanup_result.manual_review_count

    # cleaned_pdb_path is prepared_variant (env + metal) with bad H removed
    env_after_cleanup = cleanup_result.cleaned_pdb_path

    # ------------------------------------------------------------------
    # 02_contacts — deterministic TSV + optional ChimeraX script
    # ------------------------------------------------------------------
    contacts_tsv = contacts_dir / "deterministic_metal_contacts.tsv"
    analysis_records = read_pdb_atom_records(env_before_cleanup, source_role="contacts")
    # Include the metal record for distance calculations
    metal_records_extra = read_pdb_atom_records(metal_only_pdb, source_role="contacts")
    all_analysis_records = analysis_records + metal_records_extra

    # Find the metal atom in analysis records (matched by position/element)
    metal_in_analysis = next(
        (
            a for a in metal_records_extra
            if is_true_metal_atom(a) and a.element == element
        ),
        metal_atom,
    )

    _write_deterministic_contacts_tsv(
        output_path=contacts_tsv,
        metal_atom=metal_in_analysis,
        all_records=all_analysis_records,
        contact_cutoff_angstrom=contact_cutoff_angstrom,
    )

    # Coordination number + geometry classification
    contacts_list = find_metal_contacts(
        metal_atom=metal_in_analysis,
        atom_records=all_analysis_records,
        contact_cutoff_angstrom=contact_cutoff_angstrom,
    )
    site_result["coordination_number"] = len(contacts_list)
    max_contact_dist = max((c.distance_angstrom for c in contacts_list), default=2.5)

    found_geometry = classify_coordination_geometry(
        metal_atom=metal_in_analysis,
        contacts=tuple(contacts_list),
    )
    expected_geometries = TRANSITION_ION_GEOMETRY_PREFERENCES.get(element, [])
    geometry_ok = bool(expected_geometries) and found_geometry in expected_geometries
    site_result["found_geometry"] = found_geometry
    site_result["geometry_ok"] = geometry_ok

    # Optional ChimeraX script
    chimerax_script = contacts_dir / "chimera_contacts.cxc"
    _write_chimerax_cxc_script(
        script_path=chimerax_script,
        analysis_pdb_path=prepared_variant_pdb.resolve(),
        contact_cutoff_angstrom=contact_cutoff_angstrom,
        pdb_id=pdb_id,
        site_folder_name=site_folder_name,
    )

    resolved_chimera = _resolve_chimera_executable(chimera_executable)
    if resolved_chimera:
        chimera_log = logs_dir / f"chimerax_{site_folder_name}.log"
        _run_chimerax_script(chimerax_script, chimera_log, resolved_chimera)

    # ------------------------------------------------------------------
    # 03_mcpb — only generated when geometry matches standard expectations
    # ------------------------------------------------------------------
    mcpb_pdb = mcpb_dir / f"{pdb_id}_mcpb.pdb"

    if not geometry_ok:
        # Geometry is non-standard — MCPB parametrization is not appropriate
        # without manual review.  Write a note instead of the scaffold.
        (mcpb_dir / "GEOMETRY_MISMATCH.txt").write_text(
            f"Site {site_folder_name}\n"
            f"Found geometry   : {found_geometry}\n"
            f"Expected for {element:2s}  : {', '.join(expected_geometries) or 'unknown'}\n"
            "\n"
            "MCPB scaffold was NOT generated because the coordination geometry\n"
            "does not match the standard geometry for this metal ion.  Possible\n"
            "causes: incomplete coordination sphere, crystal contacts, or wrong\n"
            "protonation state.  Inspect the site manually before proceeding.\n",
            encoding="utf-8",
        )
    else:
        # Geometry matches — generate the full MCPB scaffold.
        # mcpb_pdb = env_after_cleanup already contains env + metal (bad H removed)
        shutil.copy(env_after_cleanup, mcpb_pdb)

        # Renumber atom serials sequentially from 1 to avoid duplicate-serial
        # errors in MCPB.py -s 1 that arise when multiple source PDB files are
        # concatenated (each starting at serial 1).
        _renumber_pdb_atoms(mcpb_pdb)

        # Rename GROMACS-style atom names (HB1/HB2 methylene H's) to AMBER
        # convention (HB2/HB3) so MCPB.py's ff99SB libdict lookups succeed.
        _rename_gromacs_to_amber_atoms(mcpb_pdb)

        # Rename HIS → HID/HIE/HIP based on H atoms present (GROMACS writes
        # all histidines as HIS; MCPB.py's ff99SB chargedict uses HID/HIE/HIP).
        _rename_his_by_protonation(mcpb_pdb)

        ion_serial = _find_metal_serial_in_pdb(mcpb_pdb, element)

        metal_identity = METAL_RESIDUE_IDENTITY.get(element)
        formal_charge = metal_identity[1] if metal_identity else None
        spin_guess = _UNAMBIGUOUS_SPIN.get(element)

        mol2_filename: str
        if formal_charge is not None:
            mol2_filename = f"{element}.mol2"
            _write_metal_mol2(
                output_path=mcpb_dir / mol2_filename,
                element=element,
                formal_charge=formal_charge,
                x=metal_atom.x,
                y=metal_atom.y,
                z=metal_atom.z,
            )
        else:
            mol2_filename = "TODO.mol2"

        # Detect non-standard organic ligands in the coordination shell and
        # attempt to generate GAFF2 mol2 files via antechamber (BCC charges).
        naa_resnames = _find_naa_resnames_in_contacts(contacts_list)
        naa_mol2files: list[str] = []
        for naa_rn in naa_resnames:
            naa_pdb = mcpb_dir / f"{naa_rn}.pdb"
            naa_mol2 = mcpb_dir / f"{naa_rn}.mol2"
            extracted = _extract_residue_pdb(mcpb_pdb, naa_rn, naa_pdb)
            if extracted:
                _generate_naa_mol2(naa_rn, naa_pdb, naa_mol2)
            naa_mol2files.append(f"{naa_rn}.mol2")

        # Coordinating water is handled by MCPB.py via the water_model keyword
        # (tip3p in the .in scaffold).  It must NOT appear in naa_mol2files —
        # that would override tip3p with incompatible GAFF2 atom types and break
        # MCPB.py steps 2-4.

        group_name = f"{pdb_id}_{element}_{chain}_{resseq}"
        mcpb_in = mcpb_dir / f"{pdb_id}.in"
        _write_mcpb_in(
            output_path=mcpb_in,
            pdb_id=pdb_id,
            element=element,
            chain=chain,
            resseq=resseq,
            max_contact_dist=max_contact_dist,
            ion_serial=ion_serial,
            ion_mol2file=mol2_filename,
            formal_charge=formal_charge,
            spin_guess=spin_guess,
            naa_mol2files=naa_mol2files if naa_mol2files else None,
        )
        _write_commands_sh(
            mcpb_dir / "commands.sh",
            pdb_id,
            group_name,
            mol2_filename,
            naa_mol2files=naa_mol2files if naa_mol2files else None,
        )
        _write_submit_gaussian_sh(mcpb_dir / "submit_gaussian.sh", group_name)
        _write_run_after_gaussian_sh(mcpb_dir / "run_after_gaussian.sh", pdb_id, group_name)
        _write_mcpb_readme(mcpb_dir / "README.md", pdb_id, group_name)
        site_result["mcpb_scaffold_generated"] = True
        site_result["formal_charge"] = formal_charge
        site_result["spin_multiplicity"] = spin_guess

    # ------------------------------------------------------------------
    # logs — pipeline log for this site
    # ------------------------------------------------------------------
    log_path = logs_dir / f"metall_params_{pdb_id}.log"
    log_lines = [
        f"site              : {site_folder_name}",
        f"pdb_id            : {pdb_id}",
        f"element           : {element}",
        f"chain             : {chain}",
        f"residue_number    : {resseq}",
        f"coordination      : {site_result['coordination_number']}",
        f"found_geometry    : {found_geometry}",
        f"geometry_ok       : {geometry_ok}",
        f"mcpb_scaffold     : {site_result['mcpb_scaffold_generated']}",
        f"h_removed         : {cleanup_result.removed_count}",
        f"h_kept            : {cleanup_result.kept_count}",
        f"h_manual_review   : {cleanup_result.manual_review_count}",
        f"cleanup_status    : {cleanup_result.status}",
        f"cleanup_message   : {cleanup_result.message}",
    ]
    log_path.write_text("\n".join(log_lines) + "\n", encoding="utf-8")

    site_result["status"] = "success"
    site_result["message"] = (
        f"Site {site_folder_name} prepared; H removed={cleanup_result.removed_count}."
    )
    return site_result


# ---------------------------------------------------------------------------
# Public entry point
# ---------------------------------------------------------------------------

def run_metal_parametrization_for_protein_dir(
    protein_dir: str | Path,
    chimera_executable: str | None = None,
) -> dict[str, Any]:
    """Run the metal parametrization preparation step for one protein directory.

    For every transition-metal atom found in ``<protein_dir>/components/``
    this function creates a numbered site folder under
    ``<protein_dir>/metall_params/`` containing:

    - ``00_input/``         environment PDB and metal-only PDB
    - ``01_hydrogen_cleanup/`` environment with bad H removed
    - ``02_contacts/``      deterministic contacts TSV and optional ChimeraX script
    - ``03_mcpb/``          MCPB.py scaffold (input file, commands, README)
    - ``logs/``             per-site log

    A top-level ``manifest.json`` and ``all_sites_summary.tsv`` are written to
    ``<protein_dir>/metall_params/``.

    Parameters
    ----------
    protein_dir:
        Path to the FRUTON protein directory, e.g. ``data/proteins/2AFX``.
    chimera_executable:
        Optional explicit ChimeraX/Chimera executable for optional visual
        inspection scripts.  If ``None``, the function auto-discovers ChimeraX
        from well-known locations.  Chimera output is never used for decisions.

    Returns
    -------
    dict
        Structured result for pipeline-state integration with keys:
        ``status``, ``pdb_id``, ``metall_params_dir``,
        ``transition_metal_site_count``, ``manifest_path``,
        ``site_results``, ``message``.
    """
    protein_dir = Path(protein_dir).resolve()
    pdb_id = protein_dir.name
    components_dir = protein_dir / "components"
    metall_params_dir = protein_dir / "metall_params"

    result: dict[str, Any] = {
        "status": "failed",
        "pdb_id": pdb_id,
        "metall_params_dir": str(metall_params_dir),
        "transition_metal_site_count": 0,
        "manifest_path": None,
        "site_results": [],
        "message": "",
    }

    if not components_dir.is_dir():
        result["message"] = f"Missing components directory: {components_dir}"
        return result

    if not _has_metal_file(components_dir, pdb_id):
        result["status"] = "skipped"
        result["message"] = "No metal file found. metall_params step skipped."
        return result

    protein_input = _choose_best_protein_input(components_dir, pdb_id)
    optional = _get_optional_component_paths(components_dir, pdb_id)
    water_input = optional["water"]
    ligand_input = optional["ligand"]
    metal_input = optional["metal"]

    if protein_input is None:
        result["message"] = (
            "No suitable protein input found. Expected one of: "
            f"{pdb_id}_proteinH_single.pdb, {pdb_id}_proteinH_best_complete.pdb, "
            f"{pdb_id}_proteinH_gaps.pdb, {pdb_id}_protein_monomer.pdb, "
            f"{pdb_id}_protein.pdb, {pdb_id}_protein_final.pdb, "
            f"{pdb_id}_protein_as_Amber.pdb, {pdb_id}_proteinH.pdb"
        )
        return result

    if metal_input is None:
        result["message"] = "Metal file missing despite positive metal detection."
        return result

    # Read metal records to enumerate individual transition-metal atoms
    try:
        metal_records = read_pdb_atom_records(metal_input, source_role="metal_detection")
    except Exception as exc:
        result["message"] = f"Could not read metal PDB: {exc!r}"
        return result

    transition_metals = [
        a for a in metal_records
        if is_true_metal_atom(a) and a.element in TRANSITION_METAL_ELEMENTS
    ]

    if not transition_metals:
        result["status"] = "skipped"
        result["message"] = (
            "Metal file present but contains no transition-metal atoms. "
            "metall_params step skipped."
        )
        return result

    metall_params_dir.mkdir(parents=True, exist_ok=True)

    # Also read environment for contact analysis
    env_records_for_analysis: list[Any] = []
    try:
        env_records_for_analysis = read_pdb_atom_records(
            protein_input, source_role="env_analysis"
        )
        if water_input:
            env_records_for_analysis += read_pdb_atom_records(
                water_input, source_role="env_analysis"
            )
        if ligand_input:
            env_records_for_analysis += read_pdb_atom_records(
                ligand_input, source_role="env_analysis"
            )
    except Exception:
        pass  # Non-fatal; contacts step will re-read from file

    site_results: list[dict[str, Any]] = []
    contact_cutoff = 3.5  # fixed default; could be parameterised later

    for idx, metal_atom in enumerate(transition_metals, start=1):
        try:
            site_res = _process_one_site(
                site_index=idx,
                metal_atom=metal_atom,
                all_records=env_records_for_analysis + metal_records,
                protein_dir=protein_dir,
                pdb_id=pdb_id,
                protein_input=protein_input,
                water_input=water_input,
                ligand_input=ligand_input,
                metal_input=metal_input,
                metall_params_dir=metall_params_dir,
                chimera_executable=chimera_executable,
                contact_cutoff_angstrom=contact_cutoff,
            )
        except Exception as exc:
            element = getattr(metal_atom, "element", "?")
            chain = getattr(metal_atom, "chain_id", "?") or "?"
            resseq = getattr(metal_atom, "residue_number", 0)
            site_res = {
                "site_index": idx,
                "site_folder": f"site_{idx:03d}_{element}_{chain}_{resseq}",
                "element": element,
                "chain_id": chain,
                "residue_number": resseq,
                "status": "failed",
                "hydrogen_cleanup_status": "failed",
                "removed_h_count": 0,
                "kept_h_count": 0,
                "manual_review_count": 0,
                "coordination_number": 0,
                "mcpb_folder": f"site_{idx:03d}_{element}_{chain}_{resseq}/03_mcpb",
                "message": f"Exception: {exc!r}",
            }
        site_results.append(site_res)

    result["transition_metal_site_count"] = len(transition_metals)
    result["site_results"] = site_results

    # ------------------------------------------------------------------
    # manifest.json
    # ------------------------------------------------------------------
    manifest_sites = [
        {
            "site_index": s["site_index"],
            "site_folder": s["site_folder"],
            "element": s["element"],
            "chain_id": s["chain_id"],
            "residue_number": s["residue_number"],
            "coordination_number": s["coordination_number"],
            "found_geometry": s["found_geometry"],
            "geometry_ok": s["geometry_ok"],
            "hydrogen_cleanup_status": s["hydrogen_cleanup_status"],
            "removed_h_count": s["removed_h_count"],
            "mcpb_scaffold_generated": s["mcpb_scaffold_generated"],
            "mcpb_folder": s["mcpb_folder"],
        }
        for s in site_results
    ]

    all_statuses = {s["status"] for s in site_results}
    if all_statuses == {"success"}:
        overall_status = "success"
    elif "success" in all_statuses:
        overall_status = "partial"
    elif all_statuses == {"skipped"}:
        overall_status = "skipped"
    else:
        overall_status = "failed"

    manifest = {
        "pdb_id": pdb_id,
        "status": overall_status,
        "transition_metal_site_count": len(transition_metals),
        "sites": manifest_sites,
    }
    manifest_path = metall_params_dir / "manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2), encoding="utf-8")

    # ------------------------------------------------------------------
    # all_sites_summary.tsv
    # ------------------------------------------------------------------
    summary_tsv = metall_params_dir / "all_sites_summary.tsv"
    summary_header = "\t".join([
        "site_index", "element", "chain_id", "residue_number",
        "coordination_number", "found_geometry", "geometry_ok",
        "h_removed", "h_kept", "h_manual_review",
        "mcpb_scaffold_generated", "mcpb_folder", "status",
    ])
    summary_rows = []
    for s in site_results:
        summary_rows.append("\t".join([
            str(s["site_index"]),
            s["element"],
            s["chain_id"],
            str(s["residue_number"]),
            str(s["coordination_number"]),
            s.get("found_geometry", ""),
            str(s.get("geometry_ok", False)),
            str(s["removed_h_count"]),
            str(s["kept_h_count"]),
            str(s["manual_review_count"]),
            str(s.get("mcpb_scaffold_generated", False)),
            s["mcpb_folder"],
            s["status"],
        ]))
    summary_tsv.write_text(
        summary_header + "\n" + "\n".join(summary_rows) + "\n",
        encoding="utf-8",
    )

    result["status"] = overall_status
    result["manifest_path"] = str(manifest_path)
    result["message"] = (
        f"{overall_status.capitalize()}: processed {len(transition_metals)} "
        f"transition-metal site(s) for {pdb_id}."
    )
    return result


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description=(
            "Prepare per-site MCPB scaffold folders for metal-containing "
            "protein directories."
        )
    )
    parser.add_argument(
        "protein_dir",
        type=str,
        help="Path to one protein directory, e.g. data/proteins/2AFX",
    )
    parser.add_argument(
        "--chimera",
        dest="chimera_executable",
        type=str,
        default=None,
        help="Optional explicit ChimeraX executable path or name.",
    )

    args = parser.parse_args()
    cli_result = run_metal_parametrization_for_protein_dir(
        protein_dir=args.protein_dir,
        chimera_executable=args.chimera_executable,
    )
    import json as _json
    print(_json.dumps(cli_result, indent=2))
