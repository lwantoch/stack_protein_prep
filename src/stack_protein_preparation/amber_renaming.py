"""Rename BioPython structure residues to AMBER force-field names based on hydrogen atoms."""
from __future__ import annotations

import json
import math
import re
from pathlib import Path
from typing import Iterator

from Bio.PDB import PDBIO, PDBParser
from Bio.PDB.Residue import Residue


# ---------------------------------------------------------------------------
# Private helpers
# ---------------------------------------------------------------------------


def _iter_atoms_altloc_ok(residue: Residue) -> Iterator:
    """Yield atoms whose altloc is blank, 'A', or ''."""
    for atom in residue.get_atoms():
        altloc = atom.get_altloc()
        if altloc in {" ", "A", ""}:
            yield atom


def _atom_names(residue: Residue) -> set[str]:
    """Return atom name strings for the residue, filtering by altloc."""
    return {atom.get_name().strip() for atom in _iter_atoms_altloc_ok(residue)}


def _distance(coord_a: tuple[float, float, float], coord_b: tuple[float, float, float]) -> float:
    """Euclidean distance between two 3D coordinates."""
    return math.sqrt(
        (coord_a[0] - coord_b[0]) ** 2
        + (coord_a[1] - coord_b[1]) ** 2
        + (coord_a[2] - coord_b[2]) ** 2
    )


# ---------------------------------------------------------------------------
# Public path helpers
# ---------------------------------------------------------------------------


def sanitize_variant_label(label: str) -> str:
    """Replace spaces/slashes/special chars with '_', strip, raise ValueError if empty."""
    cleaned = re.sub(r"[^A-Za-z0-9_\-]", "_", label.strip())
    cleaned = cleaned.strip("_")
    if not cleaned:
        raise ValueError(f"Invalid empty variant label derived from {label!r}")
    return cleaned


def build_standard_amber_output_path(pdb_id: str, protein_dir: Path) -> Path:
    """Return the standard AMBER output PDB path."""
    return Path(protein_dir) / "components" / f"{pdb_id}_protein_as_Amber.pdb"


def build_variant_amber_output_path(pdb_id: str, protein_dir: Path, variant_label: str) -> Path:
    """Return the variant-specific AMBER output PDB path."""
    safe = sanitize_variant_label(variant_label)
    return Path(protein_dir) / "components" / f"{pdb_id}_protein_as_Amber_{safe}.pdb"


# ---------------------------------------------------------------------------
# Renaming logic
# ---------------------------------------------------------------------------

_POLYMER_RECORD_TYPES = {"ATOM"}


def rename_structure_by_hydrogens(
    structure,
    *,
    strict_his: bool = False,
    disulf_min: float = 1.8,
    disulf_max: float = 2.2,
) -> dict[str, int]:
    """Rename residues in-place based on hydrogen atoms present.

    Returns
    -------
    dict
        Counts keyed by transition strings such as "HIS->HID", "CYS->CYX", etc.
    """
    counts: dict[str, int] = {
        "HIS->HID": 0,
        "HIS->HIE": 0,
        "HIS->HIP": 0,
        "HIS_no_protons": 0,
        "ASP->ASH": 0,
        "GLU->GLH": 0,
        "CYS->CYM": 0,
        "CYS kept": 0,
        "CYS->CYX": 0,
        "CYM->CYX": 0,
        "CYM->CYS": 0,
    }

    # First pass: rename HIS, ASP, GLU, CYS based on protonation state.
    for residue in structure.get_residues():
        resname = residue.get_resname().strip()
        # Skip water and non-ATOM records (HETATM)
        hetflag = residue.get_id()[0]
        if resname == "HOH" or hetflag not in {" ", ""}:
            continue

        names = _atom_names(residue)

        if resname == "HIS":
            has_hd1 = "HD1" in names
            has_he2 = "HE2" in names
            if has_hd1 and has_he2:
                residue.resname = "HIP"
                counts["HIS->HIP"] += 1
            elif has_hd1:
                residue.resname = "HID"
                counts["HIS->HID"] += 1
            elif has_he2:
                residue.resname = "HIE"
                counts["HIS->HIE"] += 1
            else:
                if strict_his:
                    counts["HIS_no_protons"] += 1
                # Leave as HIS

        elif resname == "ASP":
            if "HD2" in names:
                residue.resname = "ASH"
                counts["ASP->ASH"] += 1

        elif resname == "GLU":
            if "HE2" in names:
                residue.resname = "GLH"
                counts["GLU->GLH"] += 1

        elif resname == "CYS":
            if "HG" not in names:
                residue.resname = "CYM"
                counts["CYS->CYM"] += 1
            else:
                counts["CYS kept"] += 1

    # Second pass: detect disulfide bonds between CYS/CYM residues.
    # Collect all CYS and CYM residues with their SG coordinates.
    cys_cym_residues: list[tuple[Residue, list]] = []
    for residue in structure.get_residues():
        resname = residue.get_resname().strip()
        hetflag = residue.get_id()[0]
        if resname in {"CYS", "CYM"} and hetflag in {" ", ""}:
            sg_atoms = [atom for atom in residue.get_atoms() if atom.get_name().strip() == "SG"]
            if sg_atoms:
                cys_cym_residues.append((residue, sg_atoms[0].get_vector()))

    # Check pairwise distances.
    disulfide_residues: set[int] = set()
    for i, (res_i, vec_i) in enumerate(cys_cym_residues):
        for j, (res_j, vec_j) in enumerate(cys_cym_residues):
            if j <= i:
                continue
            coord_i = (vec_i[0], vec_i[1], vec_i[2])
            coord_j = (vec_j[0], vec_j[1], vec_j[2])
            dist = _distance(coord_i, coord_j)
            if disulf_min <= dist <= disulf_max:
                disulfide_residues.add(id(res_i))
                disulfide_residues.add(id(res_j))

    for residue, _vec in cys_cym_residues:
        if id(residue) in disulfide_residues:
            old_name = residue.get_resname().strip()
            residue.resname = "CYX"
            if old_name == "CYS":
                counts["CYS->CYX"] += 1
            elif old_name == "CYM":
                counts["CYM->CYX"] += 1

    return counts


# ---------------------------------------------------------------------------
# PDB I/O
# ---------------------------------------------------------------------------


def write_amber_renamed_pdb(
    input_pdb_path: Path,
    output_pdb_path: Path,
    log_path: Path | None = None,
) -> dict[str, int]:
    """Parse, rename, and write a PDB file with AMBER residue names.

    Also writes a stats JSON file alongside the output PDB.

    Parameters
    ----------
    input_pdb_path
        Source PDB file (must exist).
    output_pdb_path
        Destination PDB file.
    log_path
        Optional log file path. If None, no log is written.

    Returns
    -------
    dict
        Renaming statistics counts.

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

    stats = rename_structure_by_hydrogens(structure)

    output_pdb_path.parent.mkdir(parents=True, exist_ok=True)
    io = PDBIO()
    io.set_structure(structure)
    io.save(str(output_pdb_path))

    # Write stats JSON alongside output.
    stats_path = output_pdb_path.with_name(f"{output_pdb_path.stem}_stats.json")
    stats_path.write_text(json.dumps(stats, indent=2), encoding="utf-8")

    # Write log if requested.
    if log_path is not None:
        log_path = Path(log_path)
        log_path.parent.mkdir(parents=True, exist_ok=True)
        log_lines = [
            f"input={input_pdb_path.resolve()}",
            f"output={output_pdb_path.resolve()}",
        ]
        for key, value in stats.items():
            log_lines.append(f"{key}={value}")
        log_path.write_text("\n".join(log_lines) + "\n", encoding="utf-8")

    return stats


# ---------------------------------------------------------------------------
# High-level API
# ---------------------------------------------------------------------------

_STATS_KEY_MAP = {
    "HIS->HID": "his_to_hid",
    "HIS->HIE": "his_to_hie",
    "HIS->HIP": "his_to_hip",
    "ASP->ASH": "asp_to_ash",
    "GLU->GLH": "glu_to_glh",
    "CYS->CYM": "cys_to_cym",
    "CYS->CYX": "cys_to_cyx",
    "CYM->CYX": "cym_to_cyx",
    "CYM->CYS": "cym_to_cys",
}


def amber_rename_selected_structure(
    input_pdb_path: Path,
    output_pdb_path: Path,
    *,
    strict_his: bool = False,
    disulf_min: float = 1.8,
    disulf_max: float = 2.2,
    variant_label: str,
    log_path: Path | None = None,
) -> dict:
    """Rename residues in the selected structure and return summary dict.

    Returns
    -------
    dict
        Keys: amber_input_path, amber_output_path, amber_log_path,
        amber_renaming_success, amber_variant, amber_stats_json_path,
        asp_to_ash, glu_to_glh, his_to_hid, his_to_hie, his_to_hip,
        cys_to_cym, cym_to_cys, cys_to_cyx, cym_to_cyx.
    """
    input_pdb_path = Path(input_pdb_path)
    output_pdb_path = Path(output_pdb_path)

    try:
        if not input_pdb_path.exists():
            raise FileNotFoundError(f"Input PDB not found: {input_pdb_path}")

        parser = PDBParser(QUIET=True)
        structure = parser.get_structure(input_pdb_path.stem, str(input_pdb_path))
        stats = rename_structure_by_hydrogens(
            structure, strict_his=strict_his, disulf_min=disulf_min, disulf_max=disulf_max
        )

        output_pdb_path.parent.mkdir(parents=True, exist_ok=True)
        io = PDBIO()
        io.set_structure(structure)
        io.save(str(output_pdb_path))

        stats_path = output_pdb_path.with_name(f"{output_pdb_path.stem}_stats.json")
        stats_path.write_text(json.dumps(stats, indent=2), encoding="utf-8")

        if log_path is not None:
            log_path = Path(log_path)
            log_path.parent.mkdir(parents=True, exist_ok=True)
            log_lines = [
                f"input={input_pdb_path.resolve()}",
                f"output={output_pdb_path.resolve()}",
            ]
            for key, value in stats.items():
                log_lines.append(f"{key}={value}")
            log_path.write_text("\n".join(log_lines) + "\n", encoding="utf-8")

        success = True
    except Exception:
        stats = {k: 0 for k in _STATS_KEY_MAP}
        stats_path = output_pdb_path.with_name(f"{output_pdb_path.stem}_stats.json")
        success = False

    return {
        "amber_input_path": str(input_pdb_path.resolve()),
        "amber_output_path": str(output_pdb_path.resolve()),
        "amber_log_path": str(Path(log_path).resolve()) if log_path is not None else "",
        "amber_renaming_success": success,
        "amber_variant": variant_label,
        "amber_stats_json_path": str(
            output_pdb_path.with_name(f"{output_pdb_path.stem}_stats.json").resolve()
        ),
        "asp_to_ash": stats.get("ASP->ASH", 0),
        "glu_to_glh": stats.get("GLU->GLH", 0),
        "his_to_hid": stats.get("HIS->HID", 0),
        "his_to_hie": stats.get("HIS->HIE", 0),
        "his_to_hip": stats.get("HIS->HIP", 0),
        "cys_to_cym": stats.get("CYS->CYM", 0),
        "cym_to_cys": stats.get("CYM->CYS", 0),
        "cys_to_cyx": stats.get("CYS->CYX", 0),
        "cym_to_cyx": stats.get("CYM->CYX", 0),
    }


def amber_rename_variant_structure(
    pdb_id: str,
    protein_dir: Path,
    variant_label: str,
    input_pdb_path: Path,
) -> dict:
    """Rename residues for a variant structure using standard path conventions.

    Output is written to:
        protein_dir/components/{pdb_id}_protein_as_Amber_{sanitized_label}.pdb

    Log is written to:
        protein_dir/components/amber_renamed.log
    """
    protein_dir = Path(protein_dir)
    output_pdb_path = build_variant_amber_output_path(pdb_id, protein_dir, variant_label)
    log_path = protein_dir / "components" / "amber_renamed.log"

    return amber_rename_selected_structure(
        input_pdb_path=input_pdb_path,
        output_pdb_path=output_pdb_path,
        variant_label=variant_label,
        log_path=log_path,
    )


def amber_rename_protein_structure(
    pdb_id: str,
    protein_dir: Path,
) -> dict:
    """Rename residues for the standard protein structure using convention paths.

    Input:  protein_dir/components/{pdb_id}_proteinH.pdb
    Output: protein_dir/components/{pdb_id}_protein_as_Amber.pdb
    Log:    protein_dir/components/amber_renamed.log
    """
    protein_dir = Path(protein_dir)
    input_pdb_path = protein_dir / "components" / f"{pdb_id}_proteinH.pdb"
    output_pdb_path = build_standard_amber_output_path(pdb_id, protein_dir)
    log_path = protein_dir / "components" / "amber_renamed.log"

    return amber_rename_selected_structure(
        input_pdb_path=input_pdb_path,
        output_pdb_path=output_pdb_path,
        variant_label="standard",
        log_path=log_path,
    )
