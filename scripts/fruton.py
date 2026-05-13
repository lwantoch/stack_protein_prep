#!/usr/bin/env python3

from __future__ import annotations

import re
import sys
import traceback
from datetime import datetime
from pathlib import Path
from typing import Any

SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_ROOT_DIR = SCRIPT_DIR.parent
SRC_DIR = PROJECT_ROOT_DIR / "src"

if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

from stack_protein_preparation.amber_renaming import amber_rename_variant_structure
from stack_protein_preparation.fasta_files import create_fasta_files_for_pdb_directory
from stack_protein_preparation.filler import (
    choose_filler_chain_id_from_range_and_alignments,
    run_filler_for_chain,
)
from stack_protein_preparation.fragment_split import split_gap_variant_structure
from stack_protein_preparation.gaps import summarize_gaps
from stack_protein_preparation.insertion_codes import (
    STATUS_FAILED as INSERTION_STATUS_FAILED,
)
from stack_protein_preparation.insertion_codes import (
    STATUS_NONE as INSERTION_STATUS_NONE,
)
from stack_protein_preparation.insertion_codes import (
    STATUS_SUCCESS as INSERTION_STATUS_SUCCESS,
)
from stack_protein_preparation.insertion_codes import (
    find_input_pdb_for_protein,
    process_pdb_for_delinsertion,
)
from stack_protein_preparation.pdb_components import split_pdb_components
from stack_protein_preparation.pdb_sync import (
    read_pdb_records_from_csv,
    sync_pdb_csv_and_directories,
)
from stack_protein_preparation.pipeline_state import (
    ALIGNMENT_DIRECTORY_COLUMN_NAME,
    AMBER_INPUT_PATH_COLUMN_NAME,
    AMBER_OUTPUT_PATH_COLUMN_NAME,
    AMBER_RENAMING_STATUS_COLUMN_NAME,
    AMBER_TERMINI_INPUT_PATH_COLUMN_NAME,
    AMBER_TERMINI_OUTPUT_PATH_COLUMN_NAME,
    AMBER_TERMINI_STATUS_COLUMN_NAME,
    AVAILABLE_MODELS_COLUMN_NAME,
    COMPONENTS_DIRECTORY_COLUMN_NAME,
    FASTA_DIRECTORY_COLUMN_NAME,
    FASTA_FILES_DONE_COLUMN_NAME,
    FILLER_DIRECTORY_COLUMN_NAME,
    FILLER_MODEL_PATH_COLUMN_NAME,
    FILLER_MODEL_SOURCE_COLUMN_NAME,
    FILLER_STATUS_COLUMN_NAME,
    GAP_SIZES_COLUMN_NAME,
    HAS_GAPS_COLUMN_NAME,
    HAS_LIGANDS_COLUMN_NAME,
    HAS_METALS_COLUMN_NAME,
    HAS_NONSTANDARD_RESIDUES_COLUMN_NAME,
    INSERTION_CODES_DONE_COLUMN_NAME,
    INTERNAL_CAPPING_INPUT_PATH_COLUMN_NAME,
    INTERNAL_CAPPING_OUTPUT_PATH_COLUMN_NAME,
    INTERNAL_CAPPING_STATUS_COLUMN_NAME,
    N_GAPS_COLUMN_NAME,
    PDB_DIRECTORY_COLUMN_NAME,
    PDB_ID_COLUMN_NAME,
    PDB_SYNC_DONE_COLUMN_NAME,
    PREPARED_ALPHAFOLD_OUTPUT_PATH_COLUMN_NAME,
    PREPARED_DIRECTORY_COLUMN_NAME,
    PREPARED_GAPS_OUTPUT_PATH_COLUMN_NAME,
    PREPARED_MODELLER_OUTPUT_PATH_COLUMN_NAME,
    PREPARED_STRUCTURE_OUTPUT_PATH_COLUMN_NAME,
    PREPARED_STRUCTURE_PROTEIN_INPUT_PATH_COLUMN_NAME,
    PREPARED_STRUCTURE_STATUS_COLUMN_NAME,
    PREPARED_STRUCTURE_VARIANT_COLUMN_NAME,
    PROTONATION_INPUT_PATH_COLUMN_NAME,
    PROTONATION_INPUT_SOURCE_COLUMN_NAME,
    PROTONATION_OUTPUT_PATH_COLUMN_NAME,
    PROTONATION_STATUS_COLUMN_NAME,
    RANGE_COLUMN_NAME,
    RETAIN_ALPHAFOLD_VARIANT_COLUMN_NAME,
    RETAIN_GAPS_VARIANT_COLUMN_NAME,
    RETAIN_MODELLER_VARIANT_COLUMN_NAME,
    SEQUENCE_ALIGNMENT_DONE_COLUMN_NAME,
    STATUS_FAILED,
    STATUS_REQUIRED,
    STATUS_SKIPPED,
    STATUS_SUCCESS,
    STATUS_WARNING,
    UNIPROT_ID_COLUMN_NAME,
    VARIANT_POLICY_COLUMN_NAME,
    create_protein_record,
)
from stack_protein_preparation.pipeline_table import (
    get_record_by_pdb_id,
    load_pipeline_table,
    save_pipeline_table,
)
from stack_protein_preparation.pipeline_xlsx import write_pipeline_to_xlsx
from stack_protein_preparation.prepared_structure import (
    build_prepared_structure_for_variant,
)
from stack_protein_preparation.protonation import (
    DEFAULT_GROMACS_FORCE_FIELD,
    DEFAULT_GROMACS_WATER_MODEL,
    protonate_variant_structure,
)
from stack_protein_preparation.sequence_alignment import (
    run_alignments_for_pdb_directory,
)
from stack_protein_preparation.terminus import convert_protein_termini_to_amber

FRUTON_LOGO = r"""
███████╗██████╗ ██╗   ██╗████████╗ ██████╗ ███╗   ██╗
██╔════╝██╔══██╗██║   ██║╚══██╔══╝██╔═══██╗████╗  ██║
█████╗  ██████╔╝██║   ██║   ██║   ██║   ██║██╔██╗ ██║
██╔══╝  ██╔══██╗██║   ██║   ██║   ██║   ██║██║╚██╗██║
██║     ██║  ██║╚██████╔╝   ██║   ╚██████╔╝██║ ╚████║
╚═╝     ╚═╝  ╚═╝ ╚═════╝    ╚═╝    ╚═════╝ ╚═╝  ╚═══╝
""".strip("\n")

FRUTON_EXPANSION = (
    "Framework for Reconstruction, UniProt alignment, and "
    "Topology-Oriented protein Normalization"
)

FRUTON_LOG_DIR = PROJECT_ROOT_DIR / "data" / "proteins" / "logs" / "fruton"
FRUTON_LOG_PATH = FRUTON_LOG_DIR / "fruton.log"

FRUTON_PROTONATION_PH = 7.4
FRUTON_PROTONATION_METHOD = "gmx pdb2gmx -ignh"


def _print_logo() -> None:
    print(FRUTON_LOGO)
    print(FRUTON_EXPANSION)
    print()


def _build_box_table(
    title: str,
    rows: list[tuple[str, str]],
    left_header: str = "Metric",
    right_header: str = "Value",
) -> str:
    left_width = max(len(left_header), *(len(left) for left, _ in rows))
    right_width = max(len(right_header), *(len(right) for _, right in rows))

    inner_width = left_width + right_width + 7
    title_text = f" {title} "
    if len(title_text) > inner_width:
        title_text = title_text[: inner_width - 1]

    fill_total = inner_width - len(title_text)
    fill_left = fill_total // 2
    fill_right = fill_total - fill_left

    top = f"┌{'─' * fill_left}{title_text}{'─' * fill_right}┐"
    header_sep = f"├{'─' * (left_width + 2)}┬{'─' * (right_width + 2)}┤"
    row_sep = f"├{'─' * (left_width + 2)}┼{'─' * (right_width + 2)}┤"
    bottom = f"└{'─' * (left_width + 2)}┴{'─' * (right_width + 2)}┘"

    lines = [
        top,
        f"│ {left_header.ljust(left_width)} │ {right_header.rjust(right_width)} │",
        header_sep,
    ]

    for index, (left, right) in enumerate(rows):
        lines.append(f"│ {left.ljust(left_width)} │ {right.rjust(right_width)} │")
        if index != len(rows) - 1:
            lines.append(row_sep)

    lines.append(bottom)
    return "\n".join(lines)


def _print_summary(summary: dict[str, int]) -> None:
    print()
    _print_logo()

    rows = [
        ("PDB records read", str(summary["total_records"])),
        ("Proteins with max gap > 5", str(summary["max_gap_gt_5"])),
        ("Proteins with max gap >= 8", str(summary["max_gap_ge_8"])),
        (
            "Single best-available prepared variants written",
            str(summary["single_best_available_written"]),
        ),
        ("Gaps prepared variants written", str(summary["gaps_variant_written"])),
        (
            "MODELLER-complete prepared variants written",
            str(summary["modeller_complete_written"]),
        ),
        (
            "AlphaFold-complete prepared variants written",
            str(summary["alphafold_complete_written"]),
        ),
        (
            "Proteins with at least one prepared output",
            str(summary["proteins_with_prepared_output"]),
        ),
        ("Proteins with no prepared output", str(summary["proteins_failed"])),
    ]

    print(_build_box_table(title="FRUTON Summary", rows=rows))
    print()


def _timestamp() -> str:
    return datetime.now().strftime("%Y-%m-%d %H:%M:%S")


def _append_fruton_log(title: str, lines: list[str]) -> None:
    FRUTON_LOG_DIR.mkdir(parents=True, exist_ok=True)
    with FRUTON_LOG_PATH.open("a", encoding="utf-8") as handle:
        handle.write("═" * 100 + "\n")
        handle.write(f"[{_timestamp()}] {title}\n")
        handle.write("─" * 100 + "\n")
        for line in lines:
            handle.write(f"{line}\n")
        handle.write("\n")


def _log_fruton_exception(title: str, exc: Exception) -> None:
    _append_fruton_log(
        title,
        [
            f"exception_type: {type(exc).__name__}",
            f"exception_repr: {exc!r}",
            "traceback:",
            traceback.format_exc().rstrip(),
        ],
    )


def _screen_step(step_number: int, step_name: str) -> None:
    print(f"[FRUTON] Step {step_number}: {step_name}")


def _screen_item(message: str) -> None:
    print(f"[FRUTON]   {message}")


def _safe_log_name_part(text: str) -> str:
    normalized = re.sub(r"[^A-Za-z0-9_.-]+", "_", str(text).strip())
    return normalized.strip("_") or "unknown"


def _path_list_to_text(path_list: list[Path]) -> str:
    return ";".join(str(path) for path in path_list)


def _extract_fragment_paths_from_split_result(split_result: Any) -> list[Path]:
    if isinstance(split_result, dict):
        raw_fragment_paths = split_result.get("fragment_paths", [])
    else:
        raw_fragment_paths = getattr(split_result, "fragment_paths", [])

    fragment_path_list: list[Path] = []
    for raw_path in raw_fragment_paths:
        fragment_path = Path(raw_path)
        if fragment_path.exists():
            fragment_path_list.append(fragment_path)

    return fragment_path_list


def build_pipeline_records_from_input_csv(
    pdb_record_list: list[dict[str, str]],
    protein_data_dir: Path,
) -> list[dict[str, str]]:
    pipeline_record_list: list[dict[str, str]] = []

    for pdb_record in pdb_record_list:
        pdb_id = str(pdb_record.get(PDB_ID_COLUMN_NAME, "")).strip().upper()
        residue_range = str(pdb_record.get(RANGE_COLUMN_NAME, "")).strip()

        if not pdb_id:
            continue

        pdb_directory = protein_data_dir / pdb_id
        fasta_directory = pdb_directory / "fasta"
        alignment_directory = fasta_directory / "alignments"
        components_directory = pdb_directory / "components"
        prepared_directory = pdb_directory / "prepared"

        pipeline_record = create_protein_record(
            pdb_id=pdb_id,
            residue_range=residue_range,
        )
        pipeline_record[PDB_DIRECTORY_COLUMN_NAME] = str(pdb_directory)
        pipeline_record[FASTA_DIRECTORY_COLUMN_NAME] = str(fasta_directory)
        pipeline_record[ALIGNMENT_DIRECTORY_COLUMN_NAME] = str(alignment_directory)
        pipeline_record[COMPONENTS_DIRECTORY_COLUMN_NAME] = str(components_directory)
        pipeline_record[PREPARED_DIRECTORY_COLUMN_NAME] = str(prepared_directory)
        pipeline_record[PDB_SYNC_DONE_COLUMN_NAME] = STATUS_SUCCESS

        pipeline_record_list.append(pipeline_record)

    return pipeline_record_list


def merge_existing_and_new_pipeline_records(
    existing_pipeline_record_list: list[dict[str, str]],
    new_pipeline_record_list: list[dict[str, str]],
) -> list[dict[str, str]]:
    merged_record_list: list[dict[str, str]] = []

    for new_record in new_pipeline_record_list:
        pdb_id = new_record[PDB_ID_COLUMN_NAME]
        existing_record = get_record_by_pdb_id(existing_pipeline_record_list, pdb_id)

        if existing_record is None:
            merged_record_list.append(new_record)
            continue

        merged_record = dict(existing_record)
        merged_record[PDB_ID_COLUMN_NAME] = new_record[PDB_ID_COLUMN_NAME]
        merged_record[RANGE_COLUMN_NAME] = new_record[RANGE_COLUMN_NAME]
        merged_record[PDB_DIRECTORY_COLUMN_NAME] = new_record[PDB_DIRECTORY_COLUMN_NAME]
        merged_record[FASTA_DIRECTORY_COLUMN_NAME] = new_record[
            FASTA_DIRECTORY_COLUMN_NAME
        ]
        merged_record[ALIGNMENT_DIRECTORY_COLUMN_NAME] = new_record[
            ALIGNMENT_DIRECTORY_COLUMN_NAME
        ]
        merged_record[COMPONENTS_DIRECTORY_COLUMN_NAME] = new_record[
            COMPONENTS_DIRECTORY_COLUMN_NAME
        ]
        merged_record[PREPARED_DIRECTORY_COLUMN_NAME] = new_record[
            PREPARED_DIRECTORY_COLUMN_NAME
        ]
        merged_record[PDB_SYNC_DONE_COLUMN_NAME] = STATUS_SUCCESS

        merged_record_list.append(merged_record)

    return merged_record_list


def _find_uniprot_id_for_protein(pdb_dir: Path) -> str:
    fasta_dir = pdb_dir / "fasta"

    if not fasta_dir.exists():
        return ""

    for fasta_path in sorted(fasta_dir.glob("UniProt_*.fasta")):
        match = re.match(r"^UniProt_([A-Z0-9]+)\.fasta$", fasta_path.name)
        if match:
            return match.group(1)

    return ""


def _find_template_pdb_for_filler(
    pdb_id: str,
    pdb_dir: Path,
    components_dir: Path,
) -> Path | None:
    candidate_path_list = [
        components_dir / f"{pdb_id}_protein.pdb",
        pdb_dir / f"{pdb_id}_delins.pdb",
        pdb_dir / f"{pdb_id}.pdb",
    ]

    for candidate_path in candidate_path_list:
        if candidate_path.exists():
            return candidate_path

    return None


def _detect_model_source_from_result(
    *,
    modeller_model_path: Path | None,
    alphafold_model_path: Path | None,
) -> str:
    if alphafold_model_path is not None:
        return "alphafold"
    if modeller_model_path is not None:
        return "modeller"
    return ""


def _parse_gap_sizes(gap_sizes_text: str) -> list[int]:
    text = str(gap_sizes_text).strip().lower()
    if not text or text == "none":
        return []

    result: list[int] = []
    for token in text.split("|"):
        token = token.strip()
        if not token:
            continue
        try:
            result.append(int(token))
        except ValueError:
            continue
    return result


def _get_max_gap_size(pipeline_record: dict[str, str]) -> int:
    gap_sizes = _parse_gap_sizes(pipeline_record.get(GAP_SIZES_COLUMN_NAME, ""))
    return max(gap_sizes) if gap_sizes else 0


def _build_available_models_value(
    *,
    has_initial: bool,
    has_gaps: bool,
    has_modeller: bool,
    has_alphafold: bool,
) -> str:
    model_list: list[str] = []

    if has_initial:
        model_list.append("initial")
    if has_gaps:
        model_list.append("gaps")
    if has_modeller:
        model_list.append("modeller")
    if has_alphafold:
        model_list.append("alphafold")

    return " | ".join(model_list)


def _update_available_models_field(
    pipeline_record: dict[str, str],
) -> None:
    prepared_variant_text = str(
        pipeline_record.get(PREPARED_STRUCTURE_VARIANT_COLUMN_NAME, "")
    ).strip()
    prepared_variant_set = {
        token.strip() for token in prepared_variant_text.split(",") if token.strip()
    }

    has_initial = "single" in prepared_variant_set
    has_gaps = bool(
        str(pipeline_record.get(PREPARED_GAPS_OUTPUT_PATH_COLUMN_NAME, "")).strip()
    ) or ("gaps" in prepared_variant_set)

    has_modeller = bool(
        str(pipeline_record.get(PREPARED_MODELLER_OUTPUT_PATH_COLUMN_NAME, "")).strip()
    )
    has_alphafold = bool(
        str(pipeline_record.get(PREPARED_ALPHAFOLD_OUTPUT_PATH_COLUMN_NAME, "")).strip()
    )

    pipeline_record[AVAILABLE_MODELS_COLUMN_NAME] = _build_available_models_value(
        has_initial=has_initial,
        has_gaps=has_gaps,
        has_modeller=has_modeller,
        has_alphafold=has_alphafold,
    )


def _clear_component_fields(pipeline_record: dict[str, str]) -> None:
    pipeline_record[HAS_METALS_COLUMN_NAME] = ""
    pipeline_record[HAS_LIGANDS_COLUMN_NAME] = ""
    pipeline_record[HAS_NONSTANDARD_RESIDUES_COLUMN_NAME] = ""


def _clear_gap_fields(pipeline_record: dict[str, str]) -> None:
    pipeline_record[N_GAPS_COLUMN_NAME] = ""
    pipeline_record[GAP_SIZES_COLUMN_NAME] = ""
    pipeline_record[HAS_GAPS_COLUMN_NAME] = ""
    pipeline_record[VARIANT_POLICY_COLUMN_NAME] = ""
    pipeline_record[RETAIN_GAPS_VARIANT_COLUMN_NAME] = ""
    pipeline_record[RETAIN_MODELLER_VARIANT_COLUMN_NAME] = ""
    pipeline_record[RETAIN_ALPHAFOLD_VARIANT_COLUMN_NAME] = ""


def _clear_filler_fields(pipeline_record: dict[str, str]) -> None:
    pipeline_record[FILLER_DIRECTORY_COLUMN_NAME] = ""
    pipeline_record[FILLER_MODEL_PATH_COLUMN_NAME] = ""
    pipeline_record[FILLER_MODEL_SOURCE_COLUMN_NAME] = ""
    pipeline_record[FILLER_STATUS_COLUMN_NAME] = ""


def _clear_protonation_fields(pipeline_record: dict[str, str]) -> None:
    pipeline_record[PROTONATION_STATUS_COLUMN_NAME] = ""
    pipeline_record[PROTONATION_INPUT_SOURCE_COLUMN_NAME] = ""
    pipeline_record[PROTONATION_INPUT_PATH_COLUMN_NAME] = ""
    pipeline_record[PROTONATION_OUTPUT_PATH_COLUMN_NAME] = ""


def _clear_amber_renaming_fields(pipeline_record: dict[str, str]) -> None:
    pipeline_record[AMBER_RENAMING_STATUS_COLUMN_NAME] = ""
    pipeline_record[AMBER_INPUT_PATH_COLUMN_NAME] = ""
    pipeline_record[AMBER_OUTPUT_PATH_COLUMN_NAME] = ""


def _clear_internal_capping_fields(pipeline_record: dict[str, str]) -> None:
    pipeline_record[INTERNAL_CAPPING_STATUS_COLUMN_NAME] = ""
    pipeline_record[INTERNAL_CAPPING_INPUT_PATH_COLUMN_NAME] = ""
    pipeline_record[INTERNAL_CAPPING_OUTPUT_PATH_COLUMN_NAME] = ""


def _clear_amber_termini_fields(pipeline_record: dict[str, str]) -> None:
    pipeline_record[AMBER_TERMINI_STATUS_COLUMN_NAME] = ""
    pipeline_record[AMBER_TERMINI_INPUT_PATH_COLUMN_NAME] = ""
    pipeline_record[AMBER_TERMINI_OUTPUT_PATH_COLUMN_NAME] = ""


def _clear_prepared_structure_fields(pipeline_record: dict[str, str]) -> None:
    pipeline_record[PREPARED_STRUCTURE_STATUS_COLUMN_NAME] = ""
    pipeline_record[PREPARED_STRUCTURE_VARIANT_COLUMN_NAME] = ""
    pipeline_record[PREPARED_STRUCTURE_PROTEIN_INPUT_PATH_COLUMN_NAME] = ""
    pipeline_record[PREPARED_STRUCTURE_OUTPUT_PATH_COLUMN_NAME] = ""
    pipeline_record[PREPARED_GAPS_OUTPUT_PATH_COLUMN_NAME] = ""
    pipeline_record[PREPARED_MODELLER_OUTPUT_PATH_COLUMN_NAME] = ""
    pipeline_record[PREPARED_ALPHAFOLD_OUTPUT_PATH_COLUMN_NAME] = ""
    pipeline_record[AVAILABLE_MODELS_COLUMN_NAME] = ""


def _clear_downstream_after_gap_stage(pipeline_record: dict[str, str]) -> None:
    _clear_filler_fields(pipeline_record)
    _clear_protonation_fields(pipeline_record)
    _clear_amber_renaming_fields(pipeline_record)
    _clear_internal_capping_fields(pipeline_record)
    _clear_amber_termini_fields(pipeline_record)
    _clear_prepared_structure_fields(pipeline_record)


def _set_component_flags(
    pipeline_record: dict[str, str],
    component_summary: dict[str, object],
) -> None:
    pipeline_record[HAS_METALS_COLUMN_NAME] = (
        "yes" if bool(component_summary.get("has_metals", False)) else "no"
    )
    pipeline_record[HAS_LIGANDS_COLUMN_NAME] = (
        "yes" if bool(component_summary.get("has_ligands", False)) else "no"
    )
    pipeline_record[HAS_NONSTANDARD_RESIDUES_COLUMN_NAME] = (
        "yes"
        if bool(component_summary.get("has_nonstandard_residues", False))
        else "no"
    )


def _set_variant_policy_fields(
    pipeline_record: dict[str, str],
    n_gaps: int,
    max_gap_size: int,
    complete_model_source: str,
) -> None:
    if n_gaps == 0:
        pipeline_record[VARIANT_POLICY_COLUMN_NAME] = "single_best_available"
        pipeline_record[RETAIN_GAPS_VARIANT_COLUMN_NAME] = "no"
        pipeline_record[RETAIN_MODELLER_VARIANT_COLUMN_NAME] = "no"
        pipeline_record[RETAIN_ALPHAFOLD_VARIANT_COLUMN_NAME] = "no"
        return

    has_complete = bool(complete_model_source)

    if has_complete:
        if max_gap_size > 8:
            pipeline_record[VARIANT_POLICY_COLUMN_NAME] = "gaps_plus_large_gap_complete"
        else:
            pipeline_record[VARIANT_POLICY_COLUMN_NAME] = "gaps_plus_best_complete"
    else:
        pipeline_record[VARIANT_POLICY_COLUMN_NAME] = "gaps_only_no_complete_model"

    pipeline_record[RETAIN_GAPS_VARIANT_COLUMN_NAME] = "yes"
    pipeline_record[RETAIN_MODELLER_VARIANT_COLUMN_NAME] = (
        "yes" if complete_model_source == "modeller" else "no"
    )
    pipeline_record[RETAIN_ALPHAFOLD_VARIANT_COLUMN_NAME] = (
        "yes" if complete_model_source == "alphafold" else "no"
    )


def _choose_variant_structure_input_path(
    pdb_id: str,
    pdb_dir: Path,
    model_source: str,
    model_path: Path | None,
) -> tuple[Path, str]:
    if model_source == "default":
        default_path = pdb_dir / "components" / f"{pdb_id}_protein.pdb"
        if not default_path.exists():
            raise FileNotFoundError(f"Default protein input not found: {default_path}")
        return default_path, "protein"

    if model_path is None or not model_path.exists():
        raise FileNotFoundError(
            f"Complete-model input missing for {pdb_id}: {model_source=} {model_path=}"
        )

    if model_source == "modeller":
        return model_path, "modeller"

    if model_source == "alphafold":
        return model_path, "alphafold"

    raise ValueError(f"Unsupported model source: {model_source!r}")


def _run_fruton_protonation_for_variant(
    *,
    pdb_id: str,
    pdb_dir: Path,
    variant_label: str,
    input_pdb: Path,
    input_source: str,
) -> dict[str, str | bool | float | int]:
    """Run FRUTON-default GROMACS protonation for one variant.

    This helper keeps the runner-level protonation route explicit and avoids
    scattering GROMACS defaults across the linear and fragment chemistry paths.
    The pH value is retained as compatibility metadata because the previous
    protonation API recorded it, but the GROMACS backend does not apply pH
    directly. The actual backend route is ``gmx pdb2gmx -ignh`` using the
    AMBER-family defaults defined by the protonation module.
    """

    return protonate_variant_structure(
        pdb_id=pdb_id,
        protein_dir=pdb_dir,
        variant_label=variant_label,
        input_pdb=input_pdb,
        input_source=input_source,  # type: ignore[arg-type]
        ph=FRUTON_PROTONATION_PH,
        ff=DEFAULT_GROMACS_FORCE_FIELD,
        water_model=DEFAULT_GROMACS_WATER_MODEL,
    )


def _build_fragment_variant_label(
    base_variant_label: str,
    fragment_index: int,
) -> str:
    return f"{base_variant_label}_fragment_{fragment_index:02d}"


def _build_fragment_amber_termini_output_path(amber_output_path: Path) -> Path:
    return amber_output_path.with_name(f"{amber_output_path.stem}_amber_termini.pdb")


def _run_linear_shared_chemistry_for_variant(
    pdb_id: str,
    pdb_dir: Path,
    variant_label: str,
    model_source: str,
    model_path: Path | None,
) -> dict[str, str]:
    _append_fruton_log(
        "shared_chemistry_linear:start",
        [
            f"pdb_id             : {pdb_id}",
            f"variant_label      : {variant_label}",
            f"model_source       : {model_source}",
            f"model_path         : {model_path}",
            "internal_capping    : skipped",
            f"protonation_method : {FRUTON_PROTONATION_METHOD}",
            f"protonation_ph     : {FRUTON_PROTONATION_PH} compatibility metadata",
            f"gromacs_ff         : {DEFAULT_GROMACS_FORCE_FIELD}",
            f"gromacs_water      : {DEFAULT_GROMACS_WATER_MODEL}",
            "route              : input -> gmx pdb2gmx protonation -> amber rename -> amber termini",
        ],
    )

    selected_input_path, input_source_label = _choose_variant_structure_input_path(
        pdb_id=pdb_id,
        pdb_dir=pdb_dir,
        model_source=model_source,
        model_path=model_path,
    )

    protonation_result = _run_fruton_protonation_for_variant(
        pdb_id=pdb_id,
        pdb_dir=pdb_dir,
        variant_label=variant_label,
        input_pdb=selected_input_path,
        input_source=input_source_label,
    )

    if not protonation_result.get("protonation_success", False):
        _append_fruton_log(
            "shared_chemistry_linear:protonation_failed",
            [
                f"pdb_id                    : {pdb_id}",
                f"variant_label             : {variant_label}",
                f"input_path                : {selected_input_path}",
                f"input_source              : {input_source_label}",
                f"protonation_method        : {protonation_result.get('protonation_method', '')}",
                f"protonation_output_path   : {protonation_result.get('protonation_output_path', '')}",
                f"protonation_topology_path : {protonation_result.get('protonation_topology_path', '')}",
                f"protonation_error_message : {protonation_result.get('protonation_error_message', '')}",
            ],
        )
        return {
            "status": STATUS_FAILED,
            "variant_label": variant_label,
            "structure_input_path": str(selected_input_path),
            "structure_input_source": input_source_label,
            "internal_capping_input_path": "",
            "internal_capping_output_path": "",
            "protonation_input_source": str(
                protonation_result.get("protonation_input_source", "")
            ),
            "protonation_input_path": str(
                protonation_result.get("protonation_input_path", "")
            ),
            "protonation_output_path": str(
                protonation_result.get("protonation_output_path", "")
            ),
            "amber_input_path": "",
            "amber_output_path": "",
            "amber_termini_input_path": "",
            "amber_termini_output_path": "",
            "final_protein_input_paths": "",
        }

    protonation_output_path = Path(str(protonation_result["protonation_output_path"]))

    amber_result = amber_rename_variant_structure(
        pdb_id=pdb_id,
        protein_dir=pdb_dir,
        variant_label=variant_label,
        input_pdb_path=protonation_output_path,
        strict_his=False,
        disulf_min=1.8,
        disulf_max=2.2,
    )

    if not amber_result.get("amber_renaming_success", False):
        return {
            "status": STATUS_FAILED,
            "variant_label": variant_label,
            "structure_input_path": str(selected_input_path),
            "structure_input_source": input_source_label,
            "internal_capping_input_path": "",
            "internal_capping_output_path": "",
            "protonation_input_source": str(
                protonation_result.get("protonation_input_source", "")
            ),
            "protonation_input_path": str(
                protonation_result.get("protonation_input_path", "")
            ),
            "protonation_output_path": str(protonation_output_path),
            "amber_input_path": str(amber_result.get("amber_input_path", "")),
            "amber_output_path": str(amber_result.get("amber_output_path", "")),
            "amber_termini_input_path": "",
            "amber_termini_output_path": "",
            "final_protein_input_paths": "",
        }

    amber_output_path = Path(str(amber_result["amber_output_path"]))
    amber_termini_output_path = (
        pdb_dir / "components" / f"{pdb_id}_protein_amber_termini_{variant_label}.pdb"
    )

    convert_protein_termini_to_amber(
        input_pdb_path=amber_output_path,
        output_pdb_path=amber_termini_output_path,
    )

    if (
        not amber_termini_output_path.exists()
        or amber_termini_output_path.stat().st_size == 0
    ):
        return {
            "status": STATUS_FAILED,
            "variant_label": variant_label,
            "structure_input_path": str(selected_input_path),
            "structure_input_source": input_source_label,
            "internal_capping_input_path": "",
            "internal_capping_output_path": "",
            "protonation_input_source": str(
                protonation_result.get("protonation_input_source", "")
            ),
            "protonation_input_path": str(
                protonation_result.get("protonation_input_path", "")
            ),
            "protonation_output_path": str(protonation_output_path),
            "amber_input_path": str(amber_result["amber_input_path"]),
            "amber_output_path": str(amber_output_path),
            "amber_termini_input_path": str(amber_output_path),
            "amber_termini_output_path": str(amber_termini_output_path),
            "final_protein_input_paths": "",
        }

    return {
        "status": STATUS_SUCCESS,
        "variant_label": variant_label,
        "structure_input_path": str(selected_input_path),
        "structure_input_source": input_source_label,
        "internal_capping_input_path": "",
        "internal_capping_output_path": "",
        "protonation_input_source": str(protonation_result["protonation_input_source"]),
        "protonation_input_path": str(protonation_result["protonation_input_path"]),
        "protonation_output_path": str(protonation_output_path),
        "amber_input_path": str(amber_result["amber_input_path"]),
        "amber_output_path": str(amber_output_path),
        "amber_termini_input_path": str(amber_output_path),
        "amber_termini_output_path": str(amber_termini_output_path),
        "final_protein_input_paths": str(amber_termini_output_path),
    }


def _run_gaps_fragment_route(
    pdb_id: str,
    pdb_dir: Path,
    variant_label: str,
) -> dict[str, str]:
    if variant_label != "gaps":
        raise ValueError(
            f"_run_gaps_fragment_route only supports variant_label='gaps', got {variant_label!r}"
        )

    raw_protein_input_path = pdb_dir / "components" / f"{pdb_id}_protein.pdb"
    if not raw_protein_input_path.exists():
        raise FileNotFoundError(
            f"Missing raw protein input for gaps route: {raw_protein_input_path}"
        )

    _append_fruton_log(
        "shared_chemistry_gaps:start",
        [
            f"pdb_id                 : {pdb_id}",
            f"variant_label          : {variant_label}",
            f"raw_protein_input_path : {raw_protein_input_path}",
            "internal_capping: skipped",
            f"protonation_method     : {FRUTON_PROTONATION_METHOD}",
            f"gromacs_ff             : {DEFAULT_GROMACS_FORCE_FIELD}",
            f"gromacs_water          : {DEFAULT_GROMACS_WATER_MODEL}",
            "route: split raw protein -> gmx pdb2gmx fragments -> amber rename fragments -> amber termini fragments",
        ],
    )

    split_result = split_gap_variant_structure(
        pdb_id=pdb_id,
        protein_dir=pdb_dir,
        variant_label=variant_label,
        input_pdb_path=raw_protein_input_path,
    )

    fragment_input_path_list = _extract_fragment_paths_from_split_result(split_result)

    if not fragment_input_path_list:
        return {
            "status": STATUS_FAILED,
            "variant_label": variant_label,
            "structure_input_path": str(raw_protein_input_path),
            "structure_input_source": "protein",
            "internal_capping_input_path": "",
            "internal_capping_output_path": "",
            "protonation_input_source": "",
            "protonation_input_path": "",
            "protonation_output_path": "",
            "amber_input_path": "",
            "amber_output_path": "",
            "amber_termini_input_path": "",
            "amber_termini_output_path": "",
            "final_protein_input_paths": "",
        }

    protonation_input_path_list: list[Path] = []
    protonation_output_path_list: list[Path] = []
    amber_input_path_list: list[Path] = []
    amber_output_path_list: list[Path] = []
    amber_termini_input_path_list: list[Path] = []
    final_fragment_output_path_list: list[Path] = []

    for fragment_index, fragment_input_path in enumerate(
        fragment_input_path_list, start=1
    ):
        fragment_variant_label = _build_fragment_variant_label(
            base_variant_label=variant_label,
            fragment_index=fragment_index,
        )

        protonation_result = _run_fruton_protonation_for_variant(
            pdb_id=pdb_id,
            pdb_dir=pdb_dir,
            variant_label=fragment_variant_label,
            input_pdb=fragment_input_path,
            input_source="protein",
        )

        if not protonation_result.get("protonation_success", False):
            return {
                "status": STATUS_FAILED,
                "variant_label": variant_label,
                "structure_input_path": _path_list_to_text(fragment_input_path_list),
                "structure_input_source": "protein_fragments",
                "internal_capping_input_path": "",
                "internal_capping_output_path": "",
                "protonation_input_source": "protein_fragments",
                "protonation_input_path": _path_list_to_text(fragment_input_path_list),
                "protonation_output_path": _path_list_to_text(
                    protonation_output_path_list
                ),
                "amber_input_path": _path_list_to_text(amber_input_path_list),
                "amber_output_path": _path_list_to_text(amber_output_path_list),
                "amber_termini_input_path": _path_list_to_text(
                    amber_termini_input_path_list
                ),
                "amber_termini_output_path": _path_list_to_text(
                    final_fragment_output_path_list
                ),
                "final_protein_input_paths": "",
            }

        protonated_fragment_path = Path(
            str(protonation_result["protonation_output_path"])
        )

        amber_result = amber_rename_variant_structure(
            pdb_id=pdb_id,
            protein_dir=pdb_dir,
            variant_label=fragment_variant_label,
            input_pdb_path=protonated_fragment_path,
            strict_his=False,
            disulf_min=1.8,
            disulf_max=2.2,
        )

        if not amber_result.get("amber_renaming_success", False):
            return {
                "status": STATUS_FAILED,
                "variant_label": variant_label,
                "structure_input_path": _path_list_to_text(fragment_input_path_list),
                "structure_input_source": "protein_fragments",
                "internal_capping_input_path": "",
                "internal_capping_output_path": "",
                "protonation_input_source": "protein_fragments",
                "protonation_input_path": _path_list_to_text(fragment_input_path_list),
                "protonation_output_path": _path_list_to_text(
                    protonation_output_path_list + [protonated_fragment_path]
                ),
                "amber_input_path": _path_list_to_text(
                    amber_input_path_list + [protonated_fragment_path]
                ),
                "amber_output_path": _path_list_to_text(amber_output_path_list),
                "amber_termini_input_path": _path_list_to_text(
                    amber_termini_input_path_list
                ),
                "amber_termini_output_path": _path_list_to_text(
                    final_fragment_output_path_list
                ),
                "final_protein_input_paths": "",
            }

        amber_fragment_path = Path(str(amber_result["amber_output_path"]))
        amber_termini_fragment_path = _build_fragment_amber_termini_output_path(
            amber_output_path=amber_fragment_path
        )

        convert_protein_termini_to_amber(
            input_pdb_path=amber_fragment_path,
            output_pdb_path=amber_termini_fragment_path,
        )

        if (
            not amber_termini_fragment_path.exists()
            or amber_termini_fragment_path.stat().st_size == 0
        ):
            return {
                "status": STATUS_FAILED,
                "variant_label": variant_label,
                "structure_input_path": _path_list_to_text(fragment_input_path_list),
                "structure_input_source": "protein_fragments",
                "internal_capping_input_path": "",
                "internal_capping_output_path": "",
                "protonation_input_source": "protein_fragments",
                "protonation_input_path": _path_list_to_text(fragment_input_path_list),
                "protonation_output_path": _path_list_to_text(
                    protonation_output_path_list + [protonated_fragment_path]
                ),
                "amber_input_path": _path_list_to_text(
                    amber_input_path_list + [protonated_fragment_path]
                ),
                "amber_output_path": _path_list_to_text(
                    amber_output_path_list + [amber_fragment_path]
                ),
                "amber_termini_input_path": _path_list_to_text(
                    amber_termini_input_path_list + [amber_fragment_path]
                ),
                "amber_termini_output_path": _path_list_to_text(
                    final_fragment_output_path_list
                ),
                "final_protein_input_paths": "",
            }

        protonation_input_path_list.append(fragment_input_path)
        protonation_output_path_list.append(protonated_fragment_path)
        amber_input_path_list.append(protonated_fragment_path)
        amber_output_path_list.append(amber_fragment_path)
        amber_termini_input_path_list.append(amber_fragment_path)
        final_fragment_output_path_list.append(amber_termini_fragment_path)

    return {
        "status": STATUS_SUCCESS,
        "variant_label": variant_label,
        "structure_input_path": _path_list_to_text(fragment_input_path_list),
        "structure_input_source": "protein_fragments",
        "internal_capping_input_path": "",
        "internal_capping_output_path": "",
        "protonation_input_source": "protein_fragments",
        "protonation_input_path": _path_list_to_text(protonation_input_path_list),
        "protonation_output_path": _path_list_to_text(protonation_output_path_list),
        "amber_input_path": _path_list_to_text(amber_input_path_list),
        "amber_output_path": _path_list_to_text(amber_output_path_list),
        "amber_termini_input_path": _path_list_to_text(amber_termini_input_path_list),
        "amber_termini_output_path": _path_list_to_text(
            final_fragment_output_path_list
        ),
        "final_protein_input_paths": _path_list_to_text(
            final_fragment_output_path_list
        ),
    }


def _run_shared_chemistry_for_variant(
    pdb_id: str,
    pdb_dir: Path,
    variant_label: str,
    model_source: str,
    model_path: Path | None,
) -> dict[str, str]:
    if variant_label == "gaps":
        return _run_gaps_fragment_route(
            pdb_id=pdb_id,
            pdb_dir=pdb_dir,
            variant_label=variant_label,
        )

    return _run_linear_shared_chemistry_for_variant(
        pdb_id=pdb_id,
        pdb_dir=pdb_dir,
        variant_label=variant_label,
        model_source=model_source,
        model_path=model_path,
    )


def _build_prepared_structure_for_variant(
    pdb_id: str,
    pdb_dir: Path,
    variant_label: str,
    final_protein_input_path: Path | None = None,
    final_protein_input_paths: list[Path] | None = None,
) -> Path:
    water_input_path = pdb_dir / "components" / f"{pdb_id}_water.pdb"
    ligand_input_path = pdb_dir / "components" / f"{pdb_id}_ligand.pdb"
    metals_input_path = pdb_dir / "components" / f"{pdb_id}_metals.pdb"

    kwargs: dict[str, object] = {
        "pdb_directory": pdb_dir,
        "pdb_id": pdb_id,
        "structure_variant": variant_label,
        "water_input_path": water_input_path if water_input_path.exists() else None,
        "ligand_input_path": ligand_input_path if ligand_input_path.exists() else None,
        "metals_input_path": metals_input_path if metals_input_path.exists() else None,
    }

    if final_protein_input_paths:
        kwargs["protein_input_paths"] = final_protein_input_paths
    else:
        if final_protein_input_path is None:
            raise ValueError(
                "Either final_protein_input_path or final_protein_input_paths must be provided."
            )
        kwargs["protein_input_path"] = final_protein_input_path

    summary = build_prepared_structure_for_variant(**kwargs)
    return summary.output_pdb_path


def _get_retained_variant_plan(
    n_gaps: int,
    max_gap_size: int,
    filler_model_source: str,
    filler_model_path: Path | None,
) -> list[dict[str, str | Path | None]]:
    variant_plan_list: list[dict[str, str | Path | None]] = []

    if n_gaps == 0:
        variant_plan_list.append(
            {
                "variant_label": "single",
                "model_source": "default",
                "model_path": None,
            }
        )
        return variant_plan_list

    variant_plan_list.append(
        {
            "variant_label": "gaps",
            "model_source": "default",
            "model_path": None,
        }
    )

    if filler_model_source and filler_model_path is not None:
        if max_gap_size > 8:
            variant_plan_list.append(
                {
                    "variant_label": "large_gap_complete",
                    "model_source": filler_model_source,
                    "model_path": filler_model_path,
                }
            )
        else:
            variant_plan_list.append(
                {
                    "variant_label": "best_complete",
                    "model_source": filler_model_source,
                    "model_path": filler_model_path,
                }
            )

    return variant_plan_list


def run_pipeline() -> None:
    protein_data_dir = PROJECT_ROOT_DIR / "data" / "proteins"
    pdb_ids_csv_path = protein_data_dir / "pdb_ids.csv"
    pipeline_json_path = protein_data_dir / "pipeline.json"
    pipeline_xlsx_path = protein_data_dir / "pipeline.xlsx"

    FRUTON_LOG_DIR.mkdir(parents=True, exist_ok=True)

    summary = {
        "total_records": 0,
        "max_gap_gt_5": 0,
        "max_gap_ge_8": 0,
        "single_best_available_written": 0,
        "gaps_variant_written": 0,
        "modeller_complete_written": 0,
        "alphafold_complete_written": 0,
        "proteins_with_prepared_output": 0,
        "proteins_failed": 0,
    }

    _print_logo()
    print("[FRUTON] Starting")
    _append_fruton_log(
        "run_pipeline:start",
        [
            f"project_root_dir      : {PROJECT_ROOT_DIR}",
            f"protein_data_dir      : {protein_data_dir}",
            f"pdb_ids_csv_path      : {pdb_ids_csv_path}",
            f"pipeline_json         : {pipeline_json_path}",
            f"pipeline_xlsx         : {pipeline_xlsx_path}",
            "internal_capping_in_fruton: disabled",
            f"protonation_method    : {FRUTON_PROTONATION_METHOD}",
            f"protonation_ph_compat : {FRUTON_PROTONATION_PH}",
            f"gromacs_force_field   : {DEFAULT_GROMACS_FORCE_FIELD}",
            f"gromacs_water_model   : {DEFAULT_GROMACS_WATER_MODEL}",
        ],
    )

    _screen_step(1, "pdb_sync")
    try:
        sync_pdb_csv_and_directories(protein_data_dir)
    except Exception as exc:
        _log_fruton_exception("step_1:pdb_sync", exc)
        raise

    _screen_step(2, "read_input_csv")
    pdb_record_list = read_pdb_records_from_csv(pdb_ids_csv_path)
    summary["total_records"] = len(pdb_record_list)
    _screen_item(f"loaded {len(pdb_record_list)} input records")
    _append_fruton_log(
        "step_2:read_input_csv",
        [f"record_count: {len(pdb_record_list)}"],
    )

    _screen_step(3, "load_existing_pipeline_state")
    existing_pipeline_record_list = load_pipeline_table(pipeline_json_path)
    _screen_item(
        f"loaded {len(existing_pipeline_record_list)} existing pipeline records"
    )
    _append_fruton_log(
        "step_3:load_existing_pipeline_state",
        [f"existing_record_count: {len(existing_pipeline_record_list)}"],
    )

    _screen_step(4, "build_and_merge_pipeline_state")
    new_pipeline_record_list = build_pipeline_records_from_input_csv(
        pdb_record_list=pdb_record_list,
        protein_data_dir=protein_data_dir,
    )
    pipeline_record_list = merge_existing_and_new_pipeline_records(
        existing_pipeline_record_list=existing_pipeline_record_list,
        new_pipeline_record_list=new_pipeline_record_list,
    )
    _screen_item(f"active pipeline records: {len(pipeline_record_list)}")
    _append_fruton_log(
        "step_4:build_and_merge_pipeline_state",
        [
            f"new_record_count    : {len(new_pipeline_record_list)}",
            f"merged_record_count : {len(pipeline_record_list)}",
        ],
    )

    _screen_step(5, "fasta_files")
    for pipeline_record in pipeline_record_list:
        pdb_id = pipeline_record[PDB_ID_COLUMN_NAME]
        pdb_dir = Path(pipeline_record[PDB_DIRECTORY_COLUMN_NAME])
        residue_range = str(pipeline_record.get(RANGE_COLUMN_NAME, "")).strip()

        _screen_item(f"fasta_files -> {pdb_id}")

        try:
            create_fasta_files_for_pdb_directory(
                pdb_directory=pdb_dir,
                residue_range=residue_range,
            )

            pipeline_record[FASTA_FILES_DONE_COLUMN_NAME] = STATUS_SUCCESS
            pipeline_record[UNIPROT_ID_COLUMN_NAME] = _find_uniprot_id_for_protein(
                pdb_dir
            )
        except Exception as error:
            print(f"[FRUTON] fasta_files failed for {pdb_id}: {error!r}")
            _log_fruton_exception(f"step_5:fasta_files:{pdb_id}", error)
            pipeline_record[FASTA_FILES_DONE_COLUMN_NAME] = STATUS_FAILED
            pipeline_record[UNIPROT_ID_COLUMN_NAME] = ""
            pipeline_record[SEQUENCE_ALIGNMENT_DONE_COLUMN_NAME] = STATUS_REQUIRED
            pipeline_record[INSERTION_CODES_DONE_COLUMN_NAME] = STATUS_REQUIRED
            _clear_component_fields(pipeline_record)
            _clear_gap_fields(pipeline_record)
            _clear_downstream_after_gap_stage(pipeline_record)

    _screen_step(6, "sequence_alignment")
    for pipeline_record in pipeline_record_list:
        pdb_id = pipeline_record[PDB_ID_COLUMN_NAME]
        pdb_dir = Path(pipeline_record[PDB_DIRECTORY_COLUMN_NAME])

        if pipeline_record.get(FASTA_FILES_DONE_COLUMN_NAME, "") != STATUS_SUCCESS:
            _screen_item(f"sequence_alignment skipped for {pdb_id}")
            pipeline_record[SEQUENCE_ALIGNMENT_DONE_COLUMN_NAME] = STATUS_SKIPPED
            pipeline_record[INSERTION_CODES_DONE_COLUMN_NAME] = STATUS_REQUIRED
            _clear_component_fields(pipeline_record)
            _clear_gap_fields(pipeline_record)
            _clear_downstream_after_gap_stage(pipeline_record)
            continue

        _screen_item(f"sequence_alignment -> {pdb_id}")

        try:
            run_alignments_for_pdb_directory(pdb_dir)

            pipeline_record[SEQUENCE_ALIGNMENT_DONE_COLUMN_NAME] = STATUS_SUCCESS

            if not str(pipeline_record.get(UNIPROT_ID_COLUMN_NAME, "")).strip():
                pipeline_record[UNIPROT_ID_COLUMN_NAME] = _find_uniprot_id_for_protein(
                    pdb_dir
                )
        except Exception as error:
            print(f"[FRUTON] sequence_alignment failed for {pdb_id}: {error!r}")
            _log_fruton_exception(f"step_6:sequence_alignment:{pdb_id}", error)
            pipeline_record[SEQUENCE_ALIGNMENT_DONE_COLUMN_NAME] = STATUS_FAILED
            pipeline_record[INSERTION_CODES_DONE_COLUMN_NAME] = STATUS_REQUIRED
            _clear_component_fields(pipeline_record)
            _clear_gap_fields(pipeline_record)
            _clear_downstream_after_gap_stage(pipeline_record)

    _screen_step(7, "insertion_codes")
    for pipeline_record in pipeline_record_list:
        pdb_id = pipeline_record[PDB_ID_COLUMN_NAME]
        pdb_dir = Path(pipeline_record[PDB_DIRECTORY_COLUMN_NAME])

        if (
            pipeline_record.get(SEQUENCE_ALIGNMENT_DONE_COLUMN_NAME, "")
            != STATUS_SUCCESS
        ):
            _screen_item(f"insertion_codes skipped for {pdb_id}")
            pipeline_record[INSERTION_CODES_DONE_COLUMN_NAME] = STATUS_SKIPPED
            _clear_component_fields(pipeline_record)
            _clear_gap_fields(pipeline_record)
            _clear_downstream_after_gap_stage(pipeline_record)
            continue

        _screen_item(f"insertion_codes -> {pdb_id}")

        input_pdb_path = find_input_pdb_for_protein(pdb_dir)
        if input_pdb_path is None:
            print(f"[FRUTON] insertion_codes failed for {pdb_id}: input PDB not found")
            _append_fruton_log(
                f"step_7:insertion_codes:{pdb_id}",
                ["input PDB not found"],
            )
            pipeline_record[INSERTION_CODES_DONE_COLUMN_NAME] = STATUS_FAILED
            _clear_component_fields(pipeline_record)
            _clear_gap_fields(pipeline_record)
            _clear_downstream_after_gap_stage(pipeline_record)
            continue

        output_pdb_path = pdb_dir / f"{pdb_id}_delins.pdb"

        try:
            insertion_result = process_pdb_for_delinsertion(
                input_pdb_path=input_pdb_path,
                output_pdb_path=output_pdb_path,
            )
        except Exception as error:
            print(f"[FRUTON] insertion_codes failed for {pdb_id}: {error!r}")
            _log_fruton_exception(f"step_7:insertion_codes:{pdb_id}", error)
            pipeline_record[INSERTION_CODES_DONE_COLUMN_NAME] = STATUS_FAILED
            _clear_component_fields(pipeline_record)
            _clear_gap_fields(pipeline_record)
            _clear_downstream_after_gap_stage(pipeline_record)
            continue

        insertion_status = str(insertion_result.get("status", "")).strip().lower()

        if insertion_status in {INSERTION_STATUS_NONE, INSERTION_STATUS_SUCCESS}:
            pipeline_record[INSERTION_CODES_DONE_COLUMN_NAME] = STATUS_SUCCESS
        elif insertion_status == INSERTION_STATUS_FAILED:
            pipeline_record[INSERTION_CODES_DONE_COLUMN_NAME] = STATUS_FAILED
            _clear_component_fields(pipeline_record)
            _clear_gap_fields(pipeline_record)
            _clear_downstream_after_gap_stage(pipeline_record)
        else:
            pipeline_record[INSERTION_CODES_DONE_COLUMN_NAME] = STATUS_WARNING
            _clear_component_fields(pipeline_record)
            _clear_gap_fields(pipeline_record)
            _clear_downstream_after_gap_stage(pipeline_record)

    _screen_step(8, "component_split")
    for pipeline_record in pipeline_record_list:
        pdb_id = pipeline_record[PDB_ID_COLUMN_NAME]
        pdb_dir = Path(pipeline_record[PDB_DIRECTORY_COLUMN_NAME])
        components_dir = Path(pipeline_record[COMPONENTS_DIRECTORY_COLUMN_NAME])

        if pipeline_record.get(INSERTION_CODES_DONE_COLUMN_NAME, "") != STATUS_SUCCESS:
            _screen_item(f"component_split skipped for {pdb_id}")
            _clear_component_fields(pipeline_record)
            _clear_gap_fields(pipeline_record)
            _clear_downstream_after_gap_stage(pipeline_record)
            continue

        component_input_pdb_path = pdb_dir / f"{pdb_id}_delins.pdb"

        if not component_input_pdb_path.exists():
            _screen_item(f"component_split skipped for {pdb_id}: missing input")
            _clear_component_fields(pipeline_record)
            _clear_gap_fields(pipeline_record)
            _clear_downstream_after_gap_stage(pipeline_record)
            continue

        _screen_item(f"component_split -> {pdb_id}")

        try:
            component_summary = split_pdb_components(
                pdb_path=component_input_pdb_path,
                output_dir=components_dir,
                protein_stem=pdb_id,
            )
        except Exception as error:
            print(f"[FRUTON] component_split failed for {pdb_id}: {error!r}")
            _log_fruton_exception(f"step_8:component_split:{pdb_id}", error)
            _clear_component_fields(pipeline_record)
            _clear_gap_fields(pipeline_record)
            _clear_downstream_after_gap_stage(pipeline_record)
            continue

        pipeline_record[COMPONENTS_DIRECTORY_COLUMN_NAME] = str(
            component_summary["output_dir"]
        )
        _set_component_flags(pipeline_record, component_summary)

    _screen_step(9, "gap_detection")
    for pipeline_record in pipeline_record_list:
        pdb_id = pipeline_record[PDB_ID_COLUMN_NAME]
        components_dir = Path(pipeline_record[COMPONENTS_DIRECTORY_COLUMN_NAME])
        residue_range = str(pipeline_record.get(RANGE_COLUMN_NAME, "")).strip()

        if pipeline_record.get(INSERTION_CODES_DONE_COLUMN_NAME, "") != STATUS_SUCCESS:
            _screen_item(f"gap_detection skipped for {pdb_id}")
            _clear_gap_fields(pipeline_record)
            _clear_downstream_after_gap_stage(pipeline_record)
            continue

        gap_input_pdb_path = components_dir / f"{pdb_id}_protein.pdb"

        if not gap_input_pdb_path.exists():
            _screen_item(f"gap_detection skipped for {pdb_id}: missing protein file")
            _clear_gap_fields(pipeline_record)
            _clear_downstream_after_gap_stage(pipeline_record)
            continue

        _screen_item(f"gap_detection -> {pdb_id}")

        try:
            try:
                gap_summary = summarize_gaps(
                    gap_input_pdb_path,
                    residue_range=residue_range,
                )
            except TypeError:
                gap_summary = summarize_gaps(gap_input_pdb_path)
        except Exception as error:
            print(f"[FRUTON] gap_detection failed for {pdb_id}: {error!r}")
            _log_fruton_exception(f"step_9:gap_detection:{pdb_id}", error)
            _clear_gap_fields(pipeline_record)
            _clear_downstream_after_gap_stage(pipeline_record)
            continue

        n_gaps = int(gap_summary.get("n_gaps", 0))
        gap_sizes = gap_summary.get("gap_sizes", [])

        pipeline_record[N_GAPS_COLUMN_NAME] = str(n_gaps)
        pipeline_record[GAP_SIZES_COLUMN_NAME] = (
            "none" if n_gaps == 0 else "|".join(str(g) for g in gap_sizes)
        )
        pipeline_record[HAS_GAPS_COLUMN_NAME] = "yes" if n_gaps > 0 else "no"

        max_gap_size = max(gap_sizes) if gap_sizes else 0
        if max_gap_size > 5:
            summary["max_gap_gt_5"] += 1
        if max_gap_size >= 8:
            summary["max_gap_ge_8"] += 1

    _screen_step(10, "filler")
    for pipeline_record in pipeline_record_list:
        pdb_id = pipeline_record[PDB_ID_COLUMN_NAME]
        pdb_dir = Path(pipeline_record[PDB_DIRECTORY_COLUMN_NAME])
        alignment_dir = Path(pipeline_record[ALIGNMENT_DIRECTORY_COLUMN_NAME])
        components_dir = Path(pipeline_record[COMPONENTS_DIRECTORY_COLUMN_NAME])

        _clear_filler_fields(pipeline_record)

        n_gaps_text = str(pipeline_record.get(N_GAPS_COLUMN_NAME, "")).strip()
        if not n_gaps_text:
            _screen_item(f"filler skipped for {pdb_id}: missing gap info")
            pipeline_record[FILLER_STATUS_COLUMN_NAME] = STATUS_SKIPPED
            continue

        n_gaps = int(n_gaps_text)

        if n_gaps == 0:
            _screen_item(f"filler skipped for {pdb_id}: no gaps")
            pipeline_record[FILLER_STATUS_COLUMN_NAME] = STATUS_SKIPPED
            continue

        if pipeline_record.get(INSERTION_CODES_DONE_COLUMN_NAME, "") != STATUS_SUCCESS:
            _screen_item(f"filler skipped for {pdb_id}: upstream failure")
            pipeline_record[FILLER_STATUS_COLUMN_NAME] = STATUS_SKIPPED
            continue

        template_pdb_path = _find_template_pdb_for_filler(
            pdb_id=pdb_id,
            pdb_dir=pdb_dir,
            components_dir=components_dir,
        )

        if template_pdb_path is None or not alignment_dir.exists():
            _screen_item(f"filler skipped for {pdb_id}: missing inputs")
            pipeline_record[FILLER_STATUS_COLUMN_NAME] = STATUS_WARNING
            continue

        residue_range = str(pipeline_record.get(RANGE_COLUMN_NAME, "")).strip()

        try:
            chain_id = choose_filler_chain_id_from_range_and_alignments(
                template_pdb_path=template_pdb_path,
                alignment_directory=alignment_dir,
                residue_range=residue_range,
            )

            _screen_item(f"filler -> {pdb_id} [chain {chain_id}]")

            template_id = template_pdb_path.stem
            target_id = f"{pdb_id}_{chain_id}"
            filler_output_dir = alignment_dir / "filler" / chain_id
            uniprot_id = str(pipeline_record.get(UNIPROT_ID_COLUMN_NAME, "")).strip()

            filler_result = run_filler_for_chain(
                alignment_directory=alignment_dir,
                template_pdb_path=template_pdb_path,
                output_dir=filler_output_dir,
                template_id=template_id,
                target_id=target_id,
                chain_id=chain_id,
                final_model_name="final_filled_model.pdb",
                starting_model=1,
                ending_model=1,
                uniprot_id=uniprot_id,
                residue_range=residue_range,
            )
        except Exception as error:
            print(f"[FRUTON] filler failed for {pdb_id}: {error!r}")
            _log_fruton_exception(f"step_10:filler:{pdb_id}", error)
            pipeline_record[FILLER_STATUS_COLUMN_NAME] = STATUS_FAILED
            continue

        pipeline_record[FILLER_DIRECTORY_COLUMN_NAME] = str(filler_result.output_dir)

        representative_model_path = filler_result.final_model_path
        pipeline_record[FILLER_MODEL_PATH_COLUMN_NAME] = (
            str(representative_model_path)
            if representative_model_path is not None
            else ""
        )
        pipeline_record[FILLER_MODEL_SOURCE_COLUMN_NAME] = (
            _detect_model_source_from_result(
                modeller_model_path=filler_result.modeller_model_path,
                alphafold_model_path=filler_result.alphafold_model_path,
            )
        )

        if (
            filler_result.final_model_path is not None
            and filler_result.final_model_path.exists()
        ):
            pipeline_record[FILLER_STATUS_COLUMN_NAME] = STATUS_SUCCESS
        elif filler_result.skipped:
            pipeline_record[FILLER_STATUS_COLUMN_NAME] = STATUS_WARNING
        else:
            pipeline_record[FILLER_STATUS_COLUMN_NAME] = STATUS_WARNING

    _screen_step(11, "shared chemistry and prepared variants")
    for pipeline_record in pipeline_record_list:
        pdb_id = pipeline_record[PDB_ID_COLUMN_NAME]
        pdb_dir = Path(pipeline_record[PDB_DIRECTORY_COLUMN_NAME])

        _clear_protonation_fields(pipeline_record)
        _clear_amber_renaming_fields(pipeline_record)
        _clear_internal_capping_fields(pipeline_record)
        _clear_amber_termini_fields(pipeline_record)
        _clear_prepared_structure_fields(pipeline_record)

        n_gaps_text = str(pipeline_record.get(N_GAPS_COLUMN_NAME, "")).strip()
        if not n_gaps_text:
            _screen_item(f"prepared variants skipped for {pdb_id}: no gap data")
            pipeline_record[PREPARED_STRUCTURE_STATUS_COLUMN_NAME] = STATUS_FAILED
            pipeline_record[AVAILABLE_MODELS_COLUMN_NAME] = ""
            summary["proteins_failed"] += 1
            continue

        n_gaps = int(n_gaps_text)
        max_gap_size = _get_max_gap_size(pipeline_record)

        filler_model_source = str(
            pipeline_record.get(FILLER_MODEL_SOURCE_COLUMN_NAME, "")
        ).strip()
        filler_model_path_text = str(
            pipeline_record.get(FILLER_MODEL_PATH_COLUMN_NAME, "")
        ).strip()

        filler_model_path = None
        if filler_model_path_text:
            candidate = Path(filler_model_path_text)
            if candidate.exists():
                filler_model_path = candidate

        _set_variant_policy_fields(
            pipeline_record=pipeline_record,
            n_gaps=n_gaps,
            max_gap_size=max_gap_size,
            complete_model_source=(
                filler_model_source if filler_model_path is not None else ""
            ),
        )

        variant_plan_list = _get_retained_variant_plan(
            n_gaps=n_gaps,
            max_gap_size=max_gap_size,
            filler_model_source=(
                filler_model_source if filler_model_path is not None else ""
            ),
            filler_model_path=filler_model_path,
        )

        successful_variant_label_list: list[str] = []
        representative_final_protein_path = ""
        representative_prepared_output_path = ""

        for variant_plan in variant_plan_list:
            variant_label = str(variant_plan["variant_label"])
            model_source = str(variant_plan["model_source"])
            model_path = variant_plan["model_path"]
            model_path = model_path if isinstance(model_path, Path) else None

            _screen_item(f"shared_chemistry -> {pdb_id} [{variant_label}]")

            try:
                chemistry_result = _run_shared_chemistry_for_variant(
                    pdb_id=pdb_id,
                    pdb_dir=pdb_dir,
                    variant_label=variant_label,
                    model_source=model_source,
                    model_path=model_path,
                )
            except Exception as error:
                print(
                    f"[FRUTON] shared_chemistry failed for {pdb_id} [{variant_label}]: {error!r}"
                )
                _log_fruton_exception(
                    f"step_11:shared_chemistry:{pdb_id}:{variant_label}",
                    error,
                )
                continue

            if chemistry_result["status"] != STATUS_SUCCESS:
                continue

            pipeline_record[INTERNAL_CAPPING_STATUS_COLUMN_NAME] = STATUS_SKIPPED
            pipeline_record[INTERNAL_CAPPING_INPUT_PATH_COLUMN_NAME] = ""
            pipeline_record[INTERNAL_CAPPING_OUTPUT_PATH_COLUMN_NAME] = ""

            pipeline_record[PROTONATION_STATUS_COLUMN_NAME] = STATUS_SUCCESS
            pipeline_record[PROTONATION_INPUT_SOURCE_COLUMN_NAME] = chemistry_result[
                "protonation_input_source"
            ]
            pipeline_record[PROTONATION_INPUT_PATH_COLUMN_NAME] = chemistry_result[
                "protonation_input_path"
            ]
            pipeline_record[PROTONATION_OUTPUT_PATH_COLUMN_NAME] = chemistry_result[
                "protonation_output_path"
            ]

            pipeline_record[AMBER_RENAMING_STATUS_COLUMN_NAME] = STATUS_SUCCESS
            pipeline_record[AMBER_INPUT_PATH_COLUMN_NAME] = chemistry_result[
                "amber_input_path"
            ]
            pipeline_record[AMBER_OUTPUT_PATH_COLUMN_NAME] = chemistry_result[
                "amber_output_path"
            ]

            pipeline_record[AMBER_TERMINI_STATUS_COLUMN_NAME] = STATUS_SUCCESS
            pipeline_record[AMBER_TERMINI_INPUT_PATH_COLUMN_NAME] = chemistry_result[
                "amber_termini_input_path"
            ]
            pipeline_record[AMBER_TERMINI_OUTPUT_PATH_COLUMN_NAME] = chemistry_result[
                "amber_termini_output_path"
            ]

            final_protein_input_paths_text = chemistry_result.get(
                "final_protein_input_paths",
                "",
            ).strip()
            final_input_path_list = [
                Path(token)
                for token in final_protein_input_paths_text.split(";")
                if token.strip()
            ]

            if not final_input_path_list:
                continue

            try:
                if len(final_input_path_list) == 1:
                    prepared_output_path = _build_prepared_structure_for_variant(
                        pdb_id=pdb_id,
                        pdb_dir=pdb_dir,
                        variant_label=variant_label,
                        final_protein_input_path=final_input_path_list[0],
                        final_protein_input_paths=None,
                    )
                else:
                    prepared_output_path = _build_prepared_structure_for_variant(
                        pdb_id=pdb_id,
                        pdb_dir=pdb_dir,
                        variant_label=variant_label,
                        final_protein_input_path=None,
                        final_protein_input_paths=final_input_path_list,
                    )
            except Exception as error:
                print(
                    f"[FRUTON] prepared_structure failed for {pdb_id} [{variant_label}]: {error!r}"
                )
                _log_fruton_exception(
                    f"step_11:prepared_structure:{pdb_id}:{variant_label}",
                    error,
                )
                continue

            successful_variant_label_list.append(variant_label)

            if not representative_final_protein_path:
                representative_final_protein_path = final_protein_input_paths_text
            if not representative_prepared_output_path:
                representative_prepared_output_path = str(prepared_output_path)

            if variant_label == "single":
                summary["single_best_available_written"] += 1
            elif variant_label == "gaps":
                summary["gaps_variant_written"] += 1
                pipeline_record[PREPARED_GAPS_OUTPUT_PATH_COLUMN_NAME] = str(
                    prepared_output_path
                )
            elif variant_label in {"best_complete", "large_gap_complete"}:
                if model_source == "modeller":
                    summary["modeller_complete_written"] += 1
                    pipeline_record[PREPARED_MODELLER_OUTPUT_PATH_COLUMN_NAME] = str(
                        prepared_output_path
                    )
                elif model_source == "alphafold":
                    summary["alphafold_complete_written"] += 1
                    pipeline_record[PREPARED_ALPHAFOLD_OUTPUT_PATH_COLUMN_NAME] = str(
                        prepared_output_path
                    )

        if successful_variant_label_list:
            summary["proteins_with_prepared_output"] += 1
            pipeline_record[PREPARED_STRUCTURE_STATUS_COLUMN_NAME] = STATUS_SUCCESS
            pipeline_record[PREPARED_STRUCTURE_VARIANT_COLUMN_NAME] = ",".join(
                successful_variant_label_list
            )
            pipeline_record[PREPARED_STRUCTURE_PROTEIN_INPUT_PATH_COLUMN_NAME] = (
                representative_final_protein_path
            )
            pipeline_record[PREPARED_STRUCTURE_OUTPUT_PATH_COLUMN_NAME] = (
                representative_prepared_output_path
            )
            _update_available_models_field(pipeline_record)
        else:
            summary["proteins_failed"] += 1
            pipeline_record[INTERNAL_CAPPING_STATUS_COLUMN_NAME] = STATUS_SKIPPED
            if not pipeline_record[PROTONATION_STATUS_COLUMN_NAME]:
                pipeline_record[PROTONATION_STATUS_COLUMN_NAME] = STATUS_FAILED
            if not pipeline_record[AMBER_RENAMING_STATUS_COLUMN_NAME]:
                pipeline_record[AMBER_RENAMING_STATUS_COLUMN_NAME] = STATUS_FAILED
            if not pipeline_record[AMBER_TERMINI_STATUS_COLUMN_NAME]:
                pipeline_record[AMBER_TERMINI_STATUS_COLUMN_NAME] = STATUS_FAILED
            pipeline_record[PREPARED_STRUCTURE_STATUS_COLUMN_NAME] = STATUS_FAILED
            pipeline_record[AVAILABLE_MODELS_COLUMN_NAME] = ""

    _screen_step(12, "save_pipeline_json")
    save_pipeline_table(pipeline_record_list, pipeline_json_path)

    _screen_step(13, "write_pipeline_xlsx")
    write_pipeline_to_xlsx(pipeline_record_list, pipeline_xlsx_path)

    print(f"[FRUTON] Pipeline JSON written to: {pipeline_json_path}")
    print(f"[FRUTON] Pipeline XLSX written to: {pipeline_xlsx_path}")

    _append_fruton_log(
        "run_pipeline:finished",
        [
            f"pipeline_json_path : {pipeline_json_path}",
            f"pipeline_xlsx_path : {pipeline_xlsx_path}",
            f"summary            : {summary}",
        ],
    )

    _print_summary(summary)
    print("[FRUTON] Finished")


def main() -> None:
    run_pipeline()


if __name__ == "__main__":
    main()