# /home/grheco/repositorios/stack_protein_prep/src/stack_protein_preparation/amber_renaming.py

from __future__ import annotations

import datetime as _dt
import json
import math
import traceback
from pathlib import Path
from typing import Dict, Iterable, Set, Tuple, Union

from Bio.PDB import PDBIO, PDBParser
from Bio.PDB.Atom import Atom
from Bio.PDB.Residue import Residue

from stack_protein_preparation.pipeline_logging import (
    append_exception_traceback,
    append_key_value_block,
    append_log_text,
    environment_snapshot_text,
    print_double_header_box,
    print_light_table_box,
    print_mixed_box,
)

ALTLOC_OK: Set[str] = {" ", "A"}  # keep primary or 'A' conformer
MODULE_NAME = "amber_renaming"


def _iter_atoms_altloc_ok(residue: Residue) -> Iterable[Atom]:
    """
    Yield atoms from one residue, keeping only accepted altloc states.
    """
    for atom in residue.get_atoms():
        if atom.get_altloc() in ALTLOC_OK:
            yield atom


def _atom_names(residue: Residue) -> Set[str]:
    """
    Return a set of accepted atom names for one residue.
    """
    return {atom.get_name().strip() for atom in _iter_atoms_altloc_ok(residue)}


def _distance(a, b) -> float:
    """
    Return Euclidean distance between two 3D coordinates.
    """
    dx = float(a[0] - b[0])
    dy = float(a[1] - b[1])
    dz = float(a[2] - b[2])
    return math.sqrt(dx * dx + dy * dy + dz * dz)


def sanitize_variant_label(variant_label: str) -> str:
    """
    Convert a variant label into a filesystem-safe token.
    """
    cleaned = "".join(
        character if character.isalnum() or character in {"_", "-"} else "_"
        for character in variant_label.strip()
    )
    cleaned = cleaned.strip("_")

    if not cleaned:
        raise ValueError(f"Invalid empty variant label derived from {variant_label!r}")

    return cleaned


def build_standard_amber_output_path(
    pdb_id: str,
    protein_dir: str | Path,
) -> Path:
    """
    Standard single-path AMBER renaming output.

    Example
    -------
    data/proteins/1ABC/components/1ABC_protein_as_Amber.pdb
    """
    protein_dir = Path(protein_dir)
    components_dir = protein_dir / "components"
    return components_dir / f"{pdb_id}_protein_as_Amber.pdb"


def build_variant_amber_output_path(
    pdb_id: str,
    protein_dir: str | Path,
    variant_label: str,
) -> Path:
    """
    Variant-specific AMBER renaming output.

    Example
    -------
    data/proteins/1ABC/components/1ABC_protein_as_Amber_gaps.pdb
    data/proteins/1ABC/components/1ABC_protein_as_Amber_best_complete.pdb
    """
    protein_dir = Path(protein_dir)
    components_dir = protein_dir / "components"
    safe_variant_label = sanitize_variant_label(variant_label)
    return components_dir / f"{pdb_id}_protein_as_Amber_{safe_variant_label}.pdb"


def build_standard_amber_log_path(
    protein_dir: str | Path,
) -> Path:
    """
    Shared standard log path for single-path workflow.
    """
    protein_dir = Path(protein_dir)
    components_dir = protein_dir / "components"
    return components_dir / "amber_renamed.log"


def build_variant_amber_log_path(
    protein_dir: str | Path,
    variant_label: str,
) -> Path:
    """
    Variant-specific log path to avoid mixed logs between variants.
    """
    protein_dir = Path(protein_dir)
    components_dir = protein_dir / "components"
    safe_variant_label = sanitize_variant_label(variant_label)
    return components_dir / f"amber_renamed_{safe_variant_label}.log"


def build_stats_json_path(output_pdb_path: str | Path) -> Path:
    """
    Build the JSON stats sidecar path for one AMBER-renamed output file.
    """
    output_pdb_path = Path(output_pdb_path)
    return output_pdb_path.with_name(f"{output_pdb_path.stem}_stats.json")


def _collect_summary_counts(stats: Dict[str, int]) -> dict[str, int]:
    return {
        "his_to_hid": int(stats.get("HIS->HID", 0)),
        "his_to_hie": int(stats.get("HIS->HIE", 0)),
        "his_to_hip": int(stats.get("HIS->HIP", 0)),
        "cys_to_cym": int(stats.get("CYS->CYM", 0)),
        "cym_to_cys": int(stats.get("CYM->CYS", 0)),
        "cys_to_cyx": int(stats.get("CYS->CYX", 0)),
        "cym_to_cyx": int(stats.get("CYM->CYX", 0)),
        "asp_to_ash": int(stats.get("ASP->ASH", 0)),
        "glu_to_glh": int(stats.get("GLU->GLH", 0)),
    }


def _print_start_box(
    *,
    input_pdb_path: Path,
    output_pdb_path: Path,
    log_path: Path,
    stats_path: Path,
    variant_label: str | None,
    strict_his: bool,
    disulf_min: float,
    disulf_max: float,
) -> None:
    line_list = [
        f"Input       : {input_pdb_path}",
        f"Output      : {output_pdb_path}",
        f"Stats JSON  : {stats_path}",
        f"Log file    : {log_path}",
        f"Variant     : {variant_label or 'standard'}",
        f"strict_his  : {strict_his}",
        f"disulf_min  : {disulf_min}",
        f"disulf_max  : {disulf_max}",
    ]
    print_mixed_box(
        title=f"{MODULE_NAME} · START",
        line_list=line_list,
    )


def _print_result_table(
    *,
    stats: Dict[str, int],
) -> None:
    ordered_row_list = [
        ("HIS -> HID", str(int(stats.get("HIS->HID", 0)))),
        ("HIS -> HIE", str(int(stats.get("HIS->HIE", 0)))),
        ("HIS -> HIP", str(int(stats.get("HIS->HIP", 0)))),
        ("ASP -> ASH", str(int(stats.get("ASP->ASH", 0)))),
        ("GLU -> GLH", str(int(stats.get("GLU->GLH", 0)))),
        ("CYS -> CYM", str(int(stats.get("CYS->CYM", 0)))),
        ("CYM -> CYS", str(int(stats.get("CYM->CYS", 0)))),
        ("CYS -> CYX", str(int(stats.get("CYS->CYX", 0)))),
        ("CYM -> CYX", str(int(stats.get("CYM->CYX", 0)))),
        ("CYS kept", str(int(stats.get("CYS kept", 0)))),
        ("CYM kept", str(int(stats.get("CYM kept", 0)))),
        ("CYX kept", str(int(stats.get("CYX kept", 0)))),
        ("HIS no protons", str(int(stats.get("HIS_no_protons", 0)))),
    ]

    print_light_table_box(
        title=f"{MODULE_NAME} · residue renaming summary",
        row_list=ordered_row_list,
        left_header="Transformation",
        right_header="Count",
    )


def _print_end_box_success(
    *,
    output_pdb_path: Path,
    stats_json_path: Path,
    n_total_events: int,
) -> None:
    print_double_header_box(
        title=f"{MODULE_NAME} · END",
        line_list=[
            "Status      : success",
            f"Output PDB  : {output_pdb_path}",
            f"Stats JSON  : {stats_json_path}",
            f"Events      : {n_total_events}",
        ],
    )


def _print_end_box_failure(
    *,
    input_pdb_path: Path,
    output_pdb_path: Path,
    error_message: str,
) -> None:
    print_double_header_box(
        title=f"{MODULE_NAME} · END",
        line_list=[
            "Status      : failed",
            f"Input       : {input_pdb_path}",
            f"Output      : {output_pdb_path}",
            f"Error       : {error_message}",
        ],
    )


def rename_structure_by_hydrogens(
    structure,
    *,
    disulf_min: float = 1.8,
    disulf_max: float = 2.2,
    strict_his: bool = False,
) -> Dict[str, int]:
    """
    Rename residues in-place for AMBER-style naming.

    Rules
    -----
    - HIS -> HID / HIE / HIP based on hydrogens
    - ASP -> ASH if protonated
    - GLU -> GLH if protonated
    - CYS -> CYM if deprotonated
    - CYS/CYM in disulfide -> CYX

    Notes
    -----
    - no renumbering is performed
    - waters are ignored
    - hetero residues (ligands, ions, cofactors) are ignored
    """
    stats: Dict[str, int] = {}

    sulfur_residue_list: list[Tuple[Residue, Tuple[float, float, float]]] = []

    for model in structure:
        for chain in model:
            for residue in chain:
                resname = residue.get_resname().strip()

                if resname not in {"CYS", "CYM"}:
                    continue

                atoms = {
                    atom.get_name().strip(): atom
                    for atom in _iter_atoms_altloc_ok(residue)
                }
                sg_atom = atoms.get("SG")

                if sg_atom is None:
                    continue

                sulfur_residue_list.append((residue, tuple(sg_atom.get_coord())))

    disulfide_residue_id_set: Set[int] = set()

    for i in range(len(sulfur_residue_list)):
        residue_i, coord_i = sulfur_residue_list[i]

        for j in range(i + 1, len(sulfur_residue_list)):
            residue_j, coord_j = sulfur_residue_list[j]
            sulfur_distance = _distance(coord_i, coord_j)

            if disulf_min <= sulfur_distance <= disulf_max:
                disulfide_residue_id_set.add(id(residue_i))
                disulfide_residue_id_set.add(id(residue_j))

    for model in structure:
        for chain in model:
            for residue in chain:
                resname = residue.get_resname().strip()

                if resname in {"HOH", "WAT"}:
                    continue

                if str(residue.id[0]).startswith("H_"):
                    continue

                atom_name_set = _atom_names(residue)

                if resname in {"CYS", "CYM"}:
                    if id(residue) in disulfide_residue_id_set:
                        if resname != "CYX":
                            residue.resname = "CYX"
                            stats[f"{resname}->CYX"] = (
                                stats.get(f"{resname}->CYX", 0) + 1
                            )
                        else:
                            stats["CYX kept"] = stats.get("CYX kept", 0) + 1
                        continue

                    has_hg = ("HG" in atom_name_set) or any(
                        atom_name.endswith("HG") for atom_name in atom_name_set
                    )

                    if has_hg:
                        if resname != "CYS":
                            residue.resname = "CYS"
                            stats[f"{resname}->CYS"] = (
                                stats.get(f"{resname}->CYS", 0) + 1
                            )
                        else:
                            stats["CYS kept"] = stats.get("CYS kept", 0) + 1
                    else:
                        if resname != "CYM":
                            residue.resname = "CYM"
                            stats[f"{resname}->CYM"] = (
                                stats.get(f"{resname}->CYM", 0) + 1
                            )
                        else:
                            stats["CYM kept"] = stats.get("CYM kept", 0) + 1
                    continue

                if resname == "HIS":
                    has_hd1 = ("HD1" in atom_name_set) or any(
                        atom_name.endswith("HD1") for atom_name in atom_name_set
                    )
                    has_he2 = ("HE2" in atom_name_set) or any(
                        atom_name.endswith("HE2") for atom_name in atom_name_set
                    )

                    if has_hd1 and has_he2:
                        residue.resname = "HIP"
                        stats["HIS->HIP"] = stats.get("HIS->HIP", 0) + 1
                    elif has_hd1:
                        residue.resname = "HID"
                        stats["HIS->HID"] = stats.get("HIS->HID", 0) + 1
                    elif has_he2:
                        residue.resname = "HIE"
                        stats["HIS->HIE"] = stats.get("HIS->HIE", 0) + 1
                    else:
                        if strict_his:
                            stats["HIS_no_protons"] = stats.get("HIS_no_protons", 0) + 1
                    continue

                if resname == "ASP":
                    has_hd2 = ("HD2" in atom_name_set) or any(
                        atom_name.endswith("HD2") for atom_name in atom_name_set
                    )
                    if has_hd2:
                        residue.resname = "ASH"
                        stats["ASP->ASH"] = stats.get("ASP->ASH", 0) + 1
                    continue

                if resname == "GLU":
                    has_he2 = ("HE2" in atom_name_set) or any(
                        atom_name.endswith("HE2") for atom_name in atom_name_set
                    )
                    if has_he2:
                        residue.resname = "GLH"
                        stats["GLU->GLH"] = stats.get("GLU->GLH", 0) + 1
                    continue

    return stats


def _append_renamed_log(
    log_path: Path,
    *,
    input_pdb_path: Path,
    output_pdb_path: Path,
    stats: Dict[str, int],
    disulf_min: float,
    disulf_max: float,
    strict_his: bool,
    variant_label: str | None,
) -> None:
    """
    Append one plain-text AMBER renaming log line.
    """
    timestamp = _dt.datetime.now().isoformat(timespec="seconds")
    parts = [f"{key}={stats[key]}" for key in sorted(stats)] if stats else ["(none)"]

    line = (
        f"[{timestamp}] "
        f"variant={variant_label or 'standard'} "
        f"input={input_pdb_path} "
        f"output={output_pdb_path} "
        f"disulf_min={disulf_min} "
        f"disulf_max={disulf_max} "
        f"strict_his={strict_his} "
        f"stats={' '.join(parts)}\n"
    )

    log_path.parent.mkdir(parents=True, exist_ok=True)
    with log_path.open("a", encoding="utf-8") as handle:
        handle.write(line)

    append_key_value_block(
        log_path=log_path,
        title="RUN SUMMARY",
        key_value_pairs=[
            ("timestamp", timestamp),
            ("module", MODULE_NAME),
            ("variant", variant_label or "standard"),
            ("input_pdb_path", str(input_pdb_path)),
            ("output_pdb_path", str(output_pdb_path)),
            ("disulf_min", str(disulf_min)),
            ("disulf_max", str(disulf_max)),
            ("strict_his", str(strict_his)),
        ],
    )

    append_key_value_block(
        log_path=log_path,
        title="RENAME COUNTS",
        key_value_pairs=[(key, str(stats[key])) for key in sorted(stats)]
        or [("none", "0")],
    )


def _append_failure_log(
    log_path: Path,
    *,
    input_pdb_path: Path,
    output_pdb_path: Path,
    disulf_min: float,
    disulf_max: float,
    strict_his: bool,
    variant_label: str | None,
    error_message: str,
) -> None:
    timestamp = _dt.datetime.now().isoformat(timespec="seconds")
    line = (
        f"[{timestamp}] "
        f"variant={variant_label or 'standard'} "
        f"input={input_pdb_path} "
        f"output={output_pdb_path} "
        f"disulf_min={disulf_min} "
        f"disulf_max={disulf_max} "
        f"strict_his={strict_his} "
        f"status=FAILED "
        f"error={error_message}\n"
    )
    log_path.parent.mkdir(parents=True, exist_ok=True)
    with log_path.open("a", encoding="utf-8") as handle:
        handle.write(line)

    append_key_value_block(
        log_path=log_path,
        title="FAILED RUN",
        key_value_pairs=[
            ("timestamp", timestamp),
            ("module", MODULE_NAME),
            ("variant", variant_label or "standard"),
            ("input_pdb_path", str(input_pdb_path)),
            ("output_pdb_path", str(output_pdb_path)),
            ("disulf_min", str(disulf_min)),
            ("disulf_max", str(disulf_max)),
            ("strict_his", str(strict_his)),
            ("error_message", error_message),
        ],
    )


def write_amber_renamed_pdb(
    input_pdb_path: Union[str, Path],
    output_pdb_path: Union[str, Path],
    *,
    log_path: Union[str, Path, None] = None,
    disulf_min: float = 1.8,
    disulf_max: float = 2.2,
    strict_his: bool = False,
    variant_label: str | None = None,
) -> Dict[str, int]:
    """
    Read a protonated PDB and write a new AMBER-renamed PDB.

    Typical usage
    -------------
    input : <PDBID>_proteinH.pdb
    output: <PDBID>_protein_as_Amber.pdb

    Also writes
    -----------
    <output_stem>_stats.json
    """
    input_pdb_path = Path(input_pdb_path).resolve()
    output_pdb_path = Path(output_pdb_path).resolve()

    if not input_pdb_path.exists():
        raise FileNotFoundError(input_pdb_path)

    if log_path is None:
        log_path = output_pdb_path.parent / "amber_renamed.log"
    else:
        log_path = Path(log_path).resolve()

    stats_path = build_stats_json_path(output_pdb_path)

    _print_start_box(
        input_pdb_path=input_pdb_path,
        output_pdb_path=output_pdb_path,
        log_path=log_path,
        stats_path=stats_path,
        variant_label=variant_label,
        strict_his=strict_his,
        disulf_min=disulf_min,
        disulf_max=disulf_max,
    )

    append_key_value_block(
        log_path=log_path,
        title="ENVIRONMENT",
        key_value_pairs=[
            (line.split("=", 1)[0], line.split("=", 1)[1])
            for line in environment_snapshot_text()
        ],
    )

    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("prot", str(input_pdb_path))

    stats = rename_structure_by_hydrogens(
        structure,
        disulf_min=disulf_min,
        disulf_max=disulf_max,
        strict_his=strict_his,
    )

    output_pdb_path.parent.mkdir(parents=True, exist_ok=True)

    io = PDBIO()
    io.set_structure(structure)
    io.save(str(output_pdb_path))

    _append_renamed_log(
        log_path=log_path,
        input_pdb_path=input_pdb_path,
        output_pdb_path=output_pdb_path,
        stats=stats,
        disulf_min=disulf_min,
        disulf_max=disulf_max,
        strict_his=strict_his,
        variant_label=variant_label,
    )

    stats_path.write_text(
        json.dumps(stats, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )

    _print_result_table(stats=stats)
    _print_end_box_success(
        output_pdb_path=output_pdb_path,
        stats_json_path=stats_path,
        n_total_events=sum(int(value) for value in stats.values()),
    )

    return stats


def amber_rename_selected_structure(
    *,
    input_pdb_path: str | Path,
    output_pdb_path: str | Path,
    strict_his: bool = False,
    disulf_min: float = 1.8,
    disulf_max: float = 2.2,
    variant_label: str | None = None,
    log_path: str | Path | None = None,
) -> dict[str, str | bool | int | float]:
    """
    Rename one explicitly selected protonated structure.

    This is the variant-friendly core function. It does not decide which input
    should be used; it only renames the path it receives.

    On failure, it returns structured diagnostics instead of only raising.
    """
    input_pdb_path = Path(input_pdb_path).resolve()
    output_pdb_path = Path(output_pdb_path).resolve()

    resolved_log_path = (
        Path(log_path).resolve()
        if log_path is not None
        else (output_pdb_path.parent / "amber_renamed.log").resolve()
    )
    stats_json_path = build_stats_json_path(output_pdb_path)

    try:
        stats = write_amber_renamed_pdb(
            input_pdb_path=input_pdb_path,
            output_pdb_path=output_pdb_path,
            log_path=resolved_log_path,
            disulf_min=disulf_min,
            disulf_max=disulf_max,
            strict_his=strict_his,
            variant_label=variant_label,
        )

        amber_renaming_success = (
            output_pdb_path.is_file() and output_pdb_path.stat().st_size > 0
        )

        output_dict: dict[str, str | bool | int | float] = {
            "amber_renaming_success": amber_renaming_success,
            "amber_input_path": str(input_pdb_path),
            "amber_output_path": str(output_pdb_path),
            "amber_log_path": str(resolved_log_path),
            "amber_strict_his": strict_his,
            "amber_disulf_min": disulf_min,
            "amber_disulf_max": disulf_max,
            "amber_stats_json_path": str(stats_json_path),
            "amber_error_message": "",
            "amber_exception_type": "",
            "amber_stats_raw": json.dumps(stats, sort_keys=True),
        }

        output_dict.update(_collect_summary_counts(stats))

        if variant_label is not None:
            output_dict["amber_variant"] = variant_label

        return output_dict

    except Exception as exc:
        error_message = (
            f"AMBER renaming failed for input {input_pdb_path}: "
            f"{type(exc).__name__}: {exc}"
        )

        try:
            _append_failure_log(
                log_path=resolved_log_path,
                input_pdb_path=input_pdb_path,
                output_pdb_path=output_pdb_path,
                disulf_min=disulf_min,
                disulf_max=disulf_max,
                strict_his=strict_his,
                variant_label=variant_label,
                error_message=error_message,
            )
            append_exception_traceback(resolved_log_path, exc)
            append_log_text(resolved_log_path, traceback.format_exc().rstrip())
        except Exception:
            pass

        _print_end_box_failure(
            input_pdb_path=input_pdb_path,
            output_pdb_path=output_pdb_path,
            error_message=error_message,
        )

        output_dict = {
            "amber_renaming_success": False,
            "amber_input_path": str(input_pdb_path),
            "amber_output_path": str(output_pdb_path),
            "amber_log_path": str(resolved_log_path),
            "amber_strict_his": strict_his,
            "amber_disulf_min": disulf_min,
            "amber_disulf_max": disulf_max,
            "amber_stats_json_path": str(stats_json_path),
            "amber_error_message": error_message,
            "amber_exception_type": type(exc).__name__,
            "amber_traceback": traceback.format_exc(),
            "amber_stats_raw": "{}",
            "his_to_hid": 0,
            "his_to_hie": 0,
            "his_to_hip": 0,
            "cys_to_cym": 0,
            "cym_to_cys": 0,
            "cys_to_cyx": 0,
            "cym_to_cyx": 0,
            "asp_to_ash": 0,
            "glu_to_glh": 0,
        }

        if variant_label is not None:
            output_dict["amber_variant"] = variant_label

        return output_dict


def amber_rename_variant_structure(
    pdb_id: str,
    protein_dir: str | Path,
    *,
    variant_label: str,
    input_pdb_path: str | Path,
    strict_his: bool = False,
    disulf_min: float = 1.8,
    disulf_max: float = 2.2,
    use_variant_specific_log: bool = True,
) -> dict[str, str | bool | int | float]:
    """
    Rename one named variant and write a variant-specific output file.

    Example output path
    -------------------
    components/<PDBID>_protein_as_Amber_gaps.pdb
    components/<PDBID>_protein_as_Amber_best_complete.pdb
    """
    protein_dir = Path(protein_dir)

    output_pdb_path = build_variant_amber_output_path(
        pdb_id=pdb_id,
        protein_dir=protein_dir,
        variant_label=variant_label,
    )

    if use_variant_specific_log:
        log_path = build_variant_amber_log_path(
            protein_dir=protein_dir,
            variant_label=variant_label,
        )
    else:
        log_path = build_standard_amber_log_path(protein_dir=protein_dir)

    return amber_rename_selected_structure(
        input_pdb_path=input_pdb_path,
        output_pdb_path=output_pdb_path,
        strict_his=strict_his,
        disulf_min=disulf_min,
        disulf_max=disulf_max,
        variant_label=variant_label,
        log_path=log_path,
    )


def amber_rename_protein_structure(
    pdb_id: str,
    protein_dir: str | Path,
    *,
    strict_his: bool = False,
    disulf_min: float = 1.8,
    disulf_max: float = 2.2,
) -> dict[str, str | bool | int | float]:
    """
    Backward-compatible single-path AMBER renaming workflow.

    Expected input path:
        components/<PDBID>_proteinH.pdb

    Standardized output path:
        components/<PDBID>_protein_as_Amber.pdb
    """
    protein_dir = Path(protein_dir)
    components_dir = protein_dir / "components"

    input_pdb_path = components_dir / f"{pdb_id}_proteinH.pdb"
    output_pdb_path = build_standard_amber_output_path(
        pdb_id=pdb_id,
        protein_dir=protein_dir,
    )
    log_path = build_standard_amber_log_path(protein_dir=protein_dir)

    return amber_rename_selected_structure(
        input_pdb_path=input_pdb_path,
        output_pdb_path=output_pdb_path,
        strict_his=strict_his,
        disulf_min=disulf_min,
        disulf_max=disulf_max,
        log_path=log_path,
    )
