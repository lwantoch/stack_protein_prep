"""Internal helper functions for metall_params.py.

All private/helper functions and constants are kept here so metall_params.py
can be a thin orchestrator exposing only the public API.
"""

from __future__ import annotations

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
    filter_bidentate_carboxylate_contacts,
    find_metal_contacts,
    is_true_metal_atom,
    read_pdb_atom_records,
)
from stack_protein_preparation.pdb_components import WATER_NAMES as _WATER_NAMES
from stack_protein_preparation._metall_params_mcpb import (  # noqa: F401 (re-exported)
    _UNAMBIGUOUS_SPIN,
    _get_spin_candidates,
    _write_mcpb_in,
    _write_commands_sh,
    _write_submit_gaussian_sh,
    _write_run_after_gaussian_sh,
    _write_metal_mol2,
    _find_metal_serial_in_pdb,
    _write_mcpb_readme,
)

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


# ---------------------------------------------------------------------------
# Input selection
# ---------------------------------------------------------------------------

def _choose_best_protein_input(components_dir: Path, pdb_id: str) -> Path | None:
    """Return the best available protonated protein structure for metal contact analysis.

    Crystal-coordinate variants (single, gaps) are strongly preferred over
    AlphaFold/MODELLER complete models.  Complete models have predicted loop
    conformations that do not match the crystal zinc environment; using them
    causes the contact detection to find wrong backbone atoms near the metal.
    """
    protein_dir = components_dir.parent
    # Crystal variants first, then complete models as fallback only
    prepared_variants = ["single", "gaps", "large_gap_complete", "best_complete", "complete"]
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
    cofactor_path = components_dir / f"{pdb_id}_cofactor.pdb"
    return {
        "water": water_path if water_path.is_file() else None,
        "ligand": ligand_path if ligand_path.is_file() else None,
        "metal": metal_path if metal_path.is_file() else None,
        "cofactor": cofactor_path if cofactor_path.is_file() else None,
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


def _find_naa_resnames_in_contacts(contacts: list[Any]) -> list[str]:
    """Return unique non-standard organic ligand residue names from *contacts*."""
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


def _add_hydrogens_to_residue_pdb(input_pdb: Path, output_pdb: Path) -> bool:
    """Add hydrogens to a non-standard residue PDB using OpenBabel.

    Crystal PDB ligand records typically contain only heavy atoms.  Passing
    such an H-less PDB to ``antechamber`` produces a GAFF2 mol2 with the same
    missing-H valence structure, which then forces MCPB.py to compute an
    incorrect total electron count (often odd, forcing a wrong-state doublet
    optimisation that does not converge in Berny).  The general remedy is to
    add hydrogens from standard valences before antechamber sees the file.

    OpenBabel's ``obabel -h`` infers H count from heavy-atom connectivity and
    standard valences and writes the H-added structure to ``output_pdb``.
    Returns True when the output file is non-empty.
    """
    if not input_pdb.is_file() or input_pdb.stat().st_size == 0:
        return False
    obabel = shutil.which("obabel") or shutil.which("babel")
    if obabel is None:
        # No OpenBabel — caller can still proceed with antechamber on the
        # H-less PDB and accept the consequences (will likely produce a wrong
        # electron count for non-rigid organic ligands).
        return False
    try:
        subprocess.run(
            [obabel, str(input_pdb), "-O", str(output_pdb), "-h"],
            capture_output=True,
            text=True,
            timeout=120,
        )
    except Exception:
        return False
    return output_pdb.is_file() and output_pdb.stat().st_size > 0


def _generate_naa_mol2(
    resname: str,
    pdb_path: Path,
    mol2_path: Path,
    net_charge: int = 0,
) -> bool:
    """Try to generate a GAFF2 mol2 for *resname* using antechamber.

    Pre-step: ensure the input PDB has hydrogens.  Crystal ligand PDB records
    almost never include hydrogens, so we route the file through OpenBabel
    (``obabel -h``) first when available.  The H-added intermediate is written
    next to ``mol2_path`` and reused as antechamber's input.  Without this
    pre-step, antechamber happily produces a mol2 with the missing-H valence
    structure, and MCPB.py downstream computes an incorrect total electron
    count (typically odd, forcing a wrong-multiplicity doublet calculation
    that does not converge under Berny optimisation in Gaussian).
    """
    if not pdb_path.is_file() or pdb_path.stat().st_size == 0:
        return False
    antechamber = shutil.which("antechamber")
    if antechamber is None:
        return False

    antechamber_input = pdb_path
    h_added_pdb = mol2_path.parent / f"{resname}_H.pdb"
    if _add_hydrogens_to_residue_pdb(pdb_path, h_added_pdb):
        antechamber_input = h_added_pdb

    try:
        subprocess.run(
            [
                antechamber,
                "-i", str(antechamber_input),
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
    cofactor_path: Path | None = None,
) -> dict[str, Any]:
    """Concatenate protein + optional water + optional ligand + optional
    cofactor into one PDB.

    The cofactor file is important for metal-contact detection when the metal
    is embedded in a cofactor residue (heme HEM/HEC/HEB, iron-sulfur clusters
    FES/SF4/F3S, chlorophyll, etc.).  The porphyrin nitrogens or bridging
    sulfides ARE the metal's primary donors, so leaving the cofactor out of
    the environment file causes ``find_metal_contacts`` to see only the axial
    protein donor (e.g. Cys-SG for cytochromes) plus water, misclassifying an
    octahedral heme Fe as ``bent`` two-coordinate.
    """
    source_paths: list[Path] = [protein_path]
    if water_path is not None:
        source_paths.append(water_path)
    if ligand_path is not None:
        source_paths.append(ligand_path)
    if cofactor_path is not None:
        source_paths.append(cofactor_path)

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
    """Rewrite *pdb_path* in-place with sequential atom serial numbers starting from 1."""
    lines = pdb_path.read_text(encoding="utf-8", errors="replace").splitlines(keepends=True)
    serial = 1
    out_lines: list[str] = []
    for line in lines:
        if line.startswith(("ATOM  ", "HETATM")):
            out_lines.append(f"{line[:6]}{serial:5d}{line[11:]}")
            serial += 1
        else:
            out_lines.append(line)
    pdb_path.write_text("".join(out_lines), encoding="utf-8")


def _renumber_residues_globally(pdb_path: Path) -> None:
    """Renumber residues with a single sequential counter across all chains.

    pdb4amber renumbers residues within each chain independently (chain A:
    1..N, chain B: 1..M, ...).  MCPB.py ignores chain IDs and resolves
    residues by number alone, so two residues at position 1 in different
    chains cause KeyErrors such as 'NPRO-C1'.  This pass assigns a globally
    unique number to every (chain, resseq, icode) group in file order.
    """
    lines = pdb_path.read_text(encoding="utf-8", errors="replace").splitlines(keepends=True)
    out: list[str] = []
    counter = 0
    resmap: dict[tuple[str, str, str], int] = {}
    for ln in lines:
        if not ln.startswith(("ATOM  ", "HETATM")):
            out.append(ln)
            continue
        chain = ln[21] if len(ln) > 21 else " "
        orig_res = ln[22:26] if len(ln) > 26 else "   1"
        icode = ln[26] if len(ln) > 26 else " "
        key = (chain, orig_res, icode)
        if key not in resmap:
            counter += 1
            resmap[key] = counter
        new_res = resmap[key]
        ln = ln[:22] + f"{new_res:4d}" + ln[26:]
        out.append(ln)
    pdb_path.write_text("".join(out), encoding="utf-8")


def _normalize_mcpb_pdb(pdb_path: Path) -> None:
    """Normalize a mcpb PDB for AMBER/MCPB.py compatibility using pdb4amber.

    pdb4amber handles all GROMACS→AMBER atom naming differences and correctly
    identifies HIS protonation states (HID/HIE/HIP) from which H atoms are
    present.  A follow-up pass (_renumber_residues_globally) assigns a single
    sequential counter across all chains so MCPB.py — which ignores chain IDs
    — sees unique residue numbers everywhere.

    Preprocessing: strip MODEL/ENDMDL records that pdb2gmx sometimes emits.
    If pdb4amber is not found in PATH, falls back to the legacy per-residue
    rename table + HIS detection.
    """
    # Strip MODEL/ENDMDL records — pdb4amber (via parmed) rejects multi-model files
    text = pdb_path.read_text(encoding="utf-8", errors="replace")
    cleaned = "".join(
        ln for ln in text.splitlines(keepends=True)
        if not ln.startswith(("MODEL ", "ENDMDL"))
    )
    pdb_path.write_text(cleaned, encoding="utf-8")

    pdb4amber_exe = shutil.which("pdb4amber")
    if pdb4amber_exe is None:
        # Legacy fallback — keeps pipeline working without AmberTools in PATH
        _rename_gromacs_to_amber_atoms(pdb_path)
        _rename_his_by_protonation(pdb_path)
        _renumber_residues_globally(pdb_path)
        return

    tmp = pdb_path.with_suffix(".pdb4amber_tmp.pdb")
    try:
        result = subprocess.run(
            [pdb4amber_exe, "-i", str(pdb_path), "-o", str(tmp),
             "--no-conect", "--noter"],
            capture_output=True, text=True, timeout=120,
            cwd=str(pdb_path.parent),
        )
        if result.returncode == 0 and tmp.is_file() and tmp.stat().st_size > 0:
            shutil.move(str(tmp), str(pdb_path))
        else:
            # pdb4amber failed — fall back
            _rename_gromacs_to_amber_atoms(pdb_path)
            _rename_his_by_protonation(pdb_path)
    except Exception:
        _rename_gromacs_to_amber_atoms(pdb_path)
        _rename_his_by_protonation(pdb_path)
    finally:
        if tmp.is_file():
            tmp.unlink(missing_ok=True)
    # pdb4amber handles residue names (HIS→HID/HIE/HIP) but does not map
    # GROMACS-style atom names (e.g. HB1/HB2 on HID/HIE) to AMBER names
    # (HB2/HB3).  Run the systematic rename table for all residues after
    # pdb4amber so these are caught regardless of which residue type they
    # appear in.
    _rename_gromacs_to_amber_atoms(pdb_path)
    # Resolve any HIS left by pdb4amber when H atoms are absent (heavy-atom
    # crystal structures) → default HIE.
    _rename_his_by_protonation(pdb_path)
    # Always renumber globally after pdb4amber — it renumbers within each
    # chain independently, which leaves cross-chain collisions intact.
    _renumber_residues_globally(pdb_path)


def _rename_gromacs_to_amber_atoms(pdb_path: Path) -> None:
    """Rename GROMACS (amber99sb-ildn) atom names to AMBER (ff99SB) in-place."""
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
                    if len(new_name) < 4:
                        new_field = f" {new_name:<3s}"
                    else:
                        new_field = f"{new_name:<4s}"
                    line = line[:12] + new_field + line[16:]
        out_lines.append(line)
    pdb_path.write_text("".join(out_lines), encoding="utf-8")


def _rename_his_by_protonation(pdb_path: Path) -> None:
    """Rename HIS residues to HID/HIE/HIP in-place based on which H atoms are present."""
    lines = pdb_path.read_text(encoding="utf-8", errors="replace").splitlines(keepends=True)

    his_atoms: dict[tuple[str, str, str], set[str]] = {}
    for line in lines:
        if line.startswith(("ATOM  ", "HETATM")) and line[17:20].strip() == "HIS":
            key = (line[21].strip(), line[22:26].strip(), line[26].strip())
            atom = line[12:16].strip()
            his_atoms.setdefault(key, set()).add(atom)

    his_rename: dict[tuple[str, str, str], str] = {}
    for key, atoms in his_atoms.items():
        has_hd1 = "HD1" in atoms
        has_he2 = "HE2" in atoms
        if has_hd1 and has_he2:
            his_rename[key] = "HIP"
        elif has_hd1:
            his_rename[key] = "HID"
        else:
            his_rename[key] = "HIE"

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
    """Generate a mol2 for a coordinating water using antechamber (GAFF2/BCC)."""
    if not water_pdb.is_file() or water_pdb.stat().st_size == 0:
        return False
    antechamber = shutil.which("antechamber")
    if antechamber is None:
        return False
    try:
        subprocess.run(
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


# Standard TIP3P water mol2 used as fallback when antechamber is unavailable.
# Charges match the TIP3P model; atom types are GAFF2-compatible.
_TIP3P_HOH_MOL2 = """\
@<TRIPOS>MOLECULE
HOH
 3 2 1 0 0
SMALL
USER_CHARGES

@<TRIPOS>ATOM
      1 OW          0.0000    0.0000    0.0000 OW        1 HOH    -0.834000
      2 HW1         0.9572    0.0000    0.0000 HW        1 HOH     0.417000
      3 HW2        -0.2399    0.9267    0.0000 HW        1 HOH     0.417000
@<TRIPOS>BOND
     1     1     2    1
     2     1     3    1
@<TRIPOS>SUBSTRUCTURE
     1 HOH         1 TEMP              0 ****  ****    0 ROOT
"""


def _write_standard_hoh_mol2(mol2_path: Path) -> None:
    """Write a standard TIP3P HOH mol2 used when a water molecule coordinates a metal."""
    mol2_path.write_text(_TIP3P_HOH_MOL2, encoding="utf-8")


def _has_water_contacts(contacts: list[Any]) -> bool:
    """Return True if any contact donor is a water molecule."""
    return any(c.donor_residue_name.strip().upper() in _WATER_NAMES for c in contacts)


# ---------------------------------------------------------------------------
# Chimera executable resolution
# ---------------------------------------------------------------------------

def _resolve_chimera_executable(chimera_executable: str | None = None) -> str | None:
    """Resolve a ChimeraX (or legacy Chimera) executable."""
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
    cofactor_input: Path | None,
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
        cofactor_path=cofactor_input,
    )

    if metal_atom.raw_line:
        _write_single_metal_pdb(metal_atom.raw_line, metal_only_pdb)
    else:
        _copy_metal_only_pdb(metal_input, metal_only_pdb)

    _combined_mcpb_pdb(env_before_cleanup, metal_only_pdb, prepared_variant_pdb)

    # ------------------------------------------------------------------
    # 01_hydrogen_cleanup
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

    env_after_cleanup = cleanup_result.cleaned_pdb_path

    # ------------------------------------------------------------------
    # 02_contacts — deterministic TSV + optional ChimeraX script
    # ------------------------------------------------------------------
    contacts_tsv = contacts_dir / "deterministic_metal_contacts.tsv"
    analysis_records = read_pdb_atom_records(env_before_cleanup, source_role="contacts")
    metal_records_extra = read_pdb_atom_records(metal_only_pdb, source_role="contacts")
    all_analysis_records = analysis_records + metal_records_extra

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

    contacts_list = filter_bidentate_carboxylate_contacts(
        find_metal_contacts(
            metal_atom=metal_in_analysis,
            atom_records=all_analysis_records,
            contact_cutoff_angstrom=contact_cutoff_angstrom,
        )
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

    # Reference-oracle cross-check (2026-08-23): the PDB Chemical Component
    # Dictionary usually knows the intended ion identity + oxidation state
    # (e.g. FE2 vs FE3, CU1 vs CU2, MN vs MN3).  Look the raw metal resname
    # up in the 89-row TS reference table and validate the observed coord
    # sphere against known chemistry.  This catches assignment errors like
    # "resname says ZN but coord sphere is 6-octahedral all-O" (looks more
    # like Mg or Ca; the crystallographer may have picked the wrong ion).
    try:
        from stack_protein_preparation._metal_reference import (
            lookup_metal, validate_coordination,
        )
        _donor_elements = tuple(
            (c.donor_element or "").strip().upper() for c in contacts_list
            if getattr(c, "donor_element", None)
        )
        _ref_entry = lookup_metal(residue_name)
        _oracle_ok, _oracle_reasons = validate_coordination(
            pdb_resname=residue_name,
            coord_number=len(contacts_list),
            geometry=found_geometry,
            donor_elements=_donor_elements,
        )
        site_result["reference_ion"] = _ref_entry.pdb_resname if _ref_entry else None
        site_result["reference_consistent"] = _oracle_ok
        site_result["reference_advisory"] = _oracle_reasons
    except Exception:  # noqa: BLE001 -- oracle is advisory, never blocks
        site_result["reference_ion"] = None
        site_result["reference_consistent"] = True
        site_result["reference_advisory"] = []

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
        # Water-proposal helper: when the coord number is BELOW the expected
        # coord number for this metal AND there is a nearby crystal HOH (2.0-
        # 3.0 A from the metal), suggest promoting the water to a ligand slot
        # OR placing a fresh water at the geometric void.  Written as a
        # user-review file next to GEOMETRY_MISMATCH.txt so the user can
        # decide + rerun with active_site_overrides or by copying the
        # proposed water into the coord PDB.
        _EXPECTED_COORD = {
            "ZN": 4, "CU": 4, "FE": 6, "MN": 6, "CO": 6, "NI": 6,
            "CR": 6, "V": 6, "MO": 6, "W": 6, "CD": 4, "HG": 2,
        }
        _propose_water_msg = _propose_water_for_incomplete_coordination(
            metal_atom=metal_in_analysis,
            contacts=contacts_list,
            expected_coord=_EXPECTED_COORD.get(element),
            expected_geometries=expected_geometries,
            found_geometry=found_geometry,
            analysis_pdb_path=prepared_variant_pdb,
        )
        _ref_block = ""
        if site_result.get("reference_ion"):
            _ref_advisories = site_result.get("reference_advisory") or []
            _ref_block = (
                f"\nReference oracle for '{site_result['reference_ion']}':\n"
                + "\n".join(f"  - {r}" for r in _ref_advisories)
                + "\n"
            )
        (mcpb_dir / "GEOMETRY_MISMATCH.txt").write_text(
            f"Site {site_folder_name}\n"
            f"Found geometry   : {found_geometry}\n"
            f"Expected for {element:2s}  : {', '.join(expected_geometries) or 'unknown'}\n"
            f"Coord number     : found={len(contacts_list)}, expected~{_EXPECTED_COORD.get(element, '?')}\n"
            "\n"
            "MCPB scaffold was NOT generated because the coordination geometry\n"
            "does not match the standard geometry for this metal ion.  Possible\n"
            "causes: incomplete coordination sphere, crystal contacts, or wrong\n"
            "protonation state.  Inspect the site manually before proceeding.\n"
            + _ref_block
            + (f"\n{_propose_water_msg}\n" if _propose_water_msg else ""),
            encoding="utf-8",
        )
    else:
        shutil.copy(env_after_cleanup, mcpb_pdb)

        _renumber_pdb_atoms(mcpb_pdb)
        _normalize_mcpb_pdb(mcpb_pdb)

        ion_serial = _find_metal_serial_in_pdb(mcpb_pdb, element)

        metal_identity = METAL_RESIDUE_IDENTITY.get(element)
        formal_charge = metal_identity[1] if metal_identity else None
        # Spin selection: prefer the (element, formal_charge) lookup which covers
        # Cu(II), Ti(III), Mn(II), Fe(II/III), etc.  The first candidate is the
        # biological default (HS for first-row open-shell, LS for strong-field
        # d6); remaining candidates are recorded in the .in comment for the user.
        # Falls back to the legacy element-only view for back-compat when no
        # candidates are tabulated.
        spin_candidates = _get_spin_candidates(element, formal_charge)
        if spin_candidates:
            spin_guess = spin_candidates[0]
            spin_alternatives = spin_candidates[1:]
        else:
            spin_guess = _UNAMBIGUOUS_SPIN.get(element)
            spin_alternatives = ()

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

        naa_resnames = _find_naa_resnames_in_contacts(contacts_list)
        naa_mol2files: list[str] = []
        for naa_rn in naa_resnames:
            naa_pdb = mcpb_dir / f"{naa_rn}.pdb"
            naa_mol2 = mcpb_dir / f"{naa_rn}.mol2"
            extracted = _extract_residue_pdb(mcpb_pdb, naa_rn, naa_pdb)
            if extracted:
                # Attempt mol2 generation now (succeeds if antechamber is in PATH).
                # If it fails, commands.sh will regenerate it in the AmberTools env.
                _generate_naa_mol2(naa_rn, naa_pdb, naa_mol2)
                naa_mol2files.append(f"{naa_rn}.mol2")

        # Water molecules coordinating the metal require HOH.mol2 for MCPB.py.
        # Provide a standard TIP3P mol2 (antechamber-independent).
        if _has_water_contacts(contacts_list):
            hoh_mol2 = mcpb_dir / "HOH.mol2"
            _write_standard_hoh_mol2(hoh_mol2)
            naa_mol2files.append("HOH.mol2")

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
            spin_alternatives=spin_alternatives,
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
        site_result["spin_alternatives"] = list(spin_alternatives)

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


def _propose_water_for_incomplete_coordination(
    metal_atom,
    contacts: list,
    expected_coord: int | None,
    expected_geometries: list[str],
    found_geometry: str,
    analysis_pdb_path: Path,
) -> str:
    """Return a user-readable proposal block for GEOMETRY_MISMATCH.txt when
    the metal coordination sphere is under-populated versus the expected
    coord number.

    Two proposal modes:
      1. **Promote existing HOH**: scan the analysis PDB for HOH oxygens
         within 2.0-3.0 A of the metal.  These are ideal candidates to
         re-label as coordinated waters (WAT / HOH bound) that MCPB.py
         includes in the QM cluster.
      2. **Place new water at geometric void**: compute a unit vector
         opposite the average ligand direction (approximate void direction)
         and suggest placing a fresh water at 2.1 A along that vector.
         Written as ``PROPOSED_WATER_XYZ`` lines the user can eyeball in
         PyMOL / ChimeraX before accepting.

    Empty string when: no expected_coord known, or coord is already full,
    or there is no clean void direction.  Never mutates the analysis PDB
    -- the user reviews and re-runs FRUTON with the accepted water added.
    """
    if expected_coord is None:
        return ""
    n_missing = expected_coord - len(contacts)
    if n_missing <= 0:
        return ""

    lines: list[str] = [
        "=" * 62,
        "PROPOSED WATER COMPLETION (auto-suggested; user must review)",
        "=" * 62,
        f"The site has {len(contacts)} ligand(s) but {expected_coord} are",
        f"expected for {found_geometry or 'this geometry'}.  Missing {n_missing} ligand(s).",
        "",
    ]

    metal_coord = (
        float(metal_atom.x_orthogonal_coord),
        float(metal_atom.y_orthogonal_coord),
        float(metal_atom.z_orthogonal_coord),
    )

    # Mode 1: scan for nearby HOH oxygens (2.0-3.0 A -- typical M-OW bond)
    import math
    nearby_hoh: list[tuple[str, int, float, tuple[float, float, float]]] = []
    try:
        for raw in analysis_pdb_path.read_text(
            encoding="utf-8", errors="replace"
        ).splitlines():
            if not (raw.startswith("HETATM") or raw.startswith("ATOM")):
                continue
            if len(raw) < 54:
                continue
            resname = raw[17:20].strip().upper()
            if resname not in ("HOH", "WAT", "H2O"):
                continue
            atom_name = raw[12:16].strip()
            if atom_name not in ("O", "OW"):
                continue
            try:
                x = float(raw[30:38]); y = float(raw[38:46]); z = float(raw[46:54])
                resnum = int(raw[22:26])
            except ValueError:
                continue
            d = math.sqrt(
                (x - metal_coord[0]) ** 2
                + (y - metal_coord[1]) ** 2
                + (z - metal_coord[2]) ** 2
            )
            if 2.0 <= d <= 3.0:
                nearby_hoh.append((raw[21], resnum, d, (x, y, z)))
    except Exception:  # noqa: BLE001
        pass

    if nearby_hoh:
        nearby_hoh.sort(key=lambda t: t[2])
        lines.append(f"Mode 1 -- Promote existing HOH ({len(nearby_hoh)} candidate(s)):")
        for chain, resnum, dist, xyz in nearby_hoh[:5]:
            lines.append(
                f"  HOH {chain}:{resnum} at ({xyz[0]:.2f}, {xyz[1]:.2f}, {xyz[2]:.2f}) "
                f"is {dist:.2f} A from metal -> good coordination-water candidate"
            )
        lines.append("")
        lines.append("Action: verify the water is real (not a symmetry mate),")
        lines.append("then either (a) accept AS IS -- MCPB.py includes HETATMs in")
        lines.append("the QM cluster automatically -- or (b) if MCPB.py did not pick")
        lines.append("it up, add the resname to the FRUTON active-site override list.")
        lines.append("")

    # Mode 2: geometric void direction (opposite the average ligand direction)
    if contacts:
        vx = vy = vz = 0.0
        for c in contacts:
            try:
                # contacts are PDBAtomRecord-like objects
                lx = float(c.atom.x_orthogonal_coord)
                ly = float(c.atom.y_orthogonal_coord)
                lz = float(c.atom.z_orthogonal_coord)
            except AttributeError:
                continue
            dx = lx - metal_coord[0]
            dy = ly - metal_coord[1]
            dz = lz - metal_coord[2]
            norm = math.sqrt(dx * dx + dy * dy + dz * dz) or 1.0
            vx += dx / norm; vy += dy / norm; vz += dz / norm
        void_norm = math.sqrt(vx * vx + vy * vy + vz * vz)
        if void_norm > 0.1:
            ux = -vx / void_norm; uy = -vy / void_norm; uz = -vz / void_norm
            wat_xyz = (
                metal_coord[0] + 2.1 * ux,
                metal_coord[1] + 2.1 * uy,
                metal_coord[2] + 2.1 * uz,
            )
            lines.append("Mode 2 -- Place new water at geometric void:")
            lines.append(
                f"  PROPOSED_WATER_XYZ  {wat_xyz[0]:8.3f} {wat_xyz[1]:8.3f} {wat_xyz[2]:8.3f}"
            )
            lines.append(
                f"  (direction opposite average ligand vector, 2.1 A from metal)"
            )
            lines.append("")
            lines.append("Action: open the analysis PDB in PyMOL / ChimeraX, verify")
            lines.append("the void is a physically reasonable water pocket (no clashes")
            lines.append("with sidechains), then add HETATM record with resname HOH")
            lines.append("at the proposed XYZ and rerun FRUTON.")
        else:
            lines.append("Mode 2 -- no clean void direction (ligands nearly symmetric).")
    return "\n".join(lines)
