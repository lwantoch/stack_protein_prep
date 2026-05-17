#!/usr/bin/env python3

from __future__ import annotations

import contextlib
import io
import os
import shutil
import sys
import time
import traceback
from datetime import datetime
from pathlib import Path
from typing import Any

SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_ROOT_DIR = SCRIPT_DIR.parent
SRC_DIR = PROJECT_ROOT_DIR / "src"

if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

from stack_protein_preparation.dependencies import ensure_fruton_dependencies

FRUTON_DEPENDENCY_REPORT = ensure_fruton_dependencies(
    install_missing=True,
    include_optional_python=False,
    strict=False,
)

from stack_protein_preparation.fasta_files import create_fasta_files_for_pdb_directory
from stack_protein_preparation.filler import (
    choose_filler_chain_id_from_range_and_alignments,
    run_filler_for_chain,
)
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
    AVAILABLE_MODELS_COLUMN_NAME,
    COMPONENTS_DIRECTORY_COLUMN_NAME,
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
    PREPARED_ALPHAFOLD_OUTPUT_PATH_COLUMN_NAME,
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
    SEQUENCE_ALIGNMENT_DONE_COLUMN_NAME,
    METALL_PARAMS_MANIFEST_PATH_COLUMN_NAME,
    METALL_PARAMS_SITE_COUNT_COLUMN_NAME,
    METALL_PARAMS_STATUS_COLUMN_NAME,
    NONSTD_RESIDUE_PARAMS_MANIFEST_PATH_COLUMN_NAME,
    NONSTD_RESIDUE_PARAMS_N_RESIDUES_COLUMN_NAME,
    NONSTD_RESIDUE_PARAMS_STATUS_COLUMN_NAME,
    STATUS_FAILED,
    STATUS_REQUIRED,
    STATUS_SKIPPED,
    STATUS_SUCCESS,
    STATUS_WARNING,
    UNIPROT_ID_COLUMN_NAME,
)
from stack_protein_preparation.pipeline_table import (
    load_pipeline_table,
    save_pipeline_table,
)
from stack_protein_preparation.pipeline_xlsx import write_pipeline_to_xlsx
from stack_protein_preparation.protonation import (
    DEFAULT_GROMACS_FORCE_FIELD,
    DEFAULT_GROMACS_WATER_MODEL,
)
from stack_protein_preparation.sequence_alignment import (
    run_alignments_for_pdb_directory,
)
from stack_protein_preparation import pipeline_runner
from stack_protein_preparation.pipeline_records import (
    build_pipeline_records_from_input_csv,
    merge_existing_and_new_pipeline_records,
    _clear_component_fields,
    _clear_gap_fields,
    _clear_filler_fields,
    _clear_protonation_fields,
    _clear_internal_capping_fields,
    _clear_prepared_structure_fields,
    _clear_downstream_after_gap_stage,
    _set_component_flags,
    _clear_prepared_variant_output_paths,
    _recompute_prepared_summary_from_state,
)
from stack_protein_preparation.pipeline_variants import (
    _detect_model_source_from_result,
    _get_max_gap_size,
    _update_available_models_field,
    _set_variant_policy_fields,
    _get_retained_variant_plan,
    _choose_representative_successful_variant,
    _blocks_complete_model_generation,
    _select_default_protein_input_path,
)
from stack_protein_preparation.pipeline_runner import (
    _clear_component_output_files,
    _write_representative_monomer_for_protein,
    _copy_representative_protein_component_to_monomer_path,
    _find_uniprot_id_for_protein,
    _find_template_pdb_for_filler,
    _run_sanitize_for_representative_protein,
    _run_sanitize_for_prepared_variant,
    _run_parameter_audit_for_variant,
    _append_parameter_audit_variant_summary,
    _run_metals_check_for_input_inventory,
    _run_metals_check_for_prepared_variant,
    _store_metal_check_result_in_record,
    _variant_metal_check_accepts_end_model,
    _variant_audit_accepts_end_model,
    _apply_accepted_variant_audit_decision,
    _run_protonation_route_for_variant,
    _build_prepared_structure_for_variant,
    _run_metall_params_for_protein,
    _run_nonstd_residue_params_for_protein,
)

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


def _print_dependency_report() -> None:
    required_missing = len(FRUTON_DEPENDENCY_REPORT.missing_required)
    optional_missing = len(FRUTON_DEPENDENCY_REPORT.missing_optional)
    print("┏━ FRUTON dependency check")
    print(f"┃ required missing : {required_missing}")
    print(f"┃ optional missing : {optional_missing}")
    if required_missing:
        for check in FRUTON_DEPENDENCY_REPORT.missing_required:
            print(f"┃ missing required : {check.name} — {check.hint}")
    if optional_missing and _verbose_screen_enabled():
        for check in FRUTON_DEPENDENCY_REPORT.missing_optional:
            print(f"┃ optional missing : {check.name}")
    print()


def _build_summary_box(
    title: str,
    rows: list[tuple[str, str]],
    left_header: str = "Metric",
    right_header: str = "Value",
) -> str:
    left_width = max(len(left_header), *(len(left) for left, _ in rows))
    right_width = max(len(right_header), *(len(right) for _, right in rows))

    logo_width = max(len(line) for line in FRUTON_LOGO.splitlines())
    minimum_table_width = left_width + right_width + 5
    content_width = max(logo_width, minimum_table_width, len(title) + 4)

    left_width = content_width - right_width - 5

    top = f"╔{'═' * content_width}╗"
    bottom = f"╚{'═' * content_width}╝"

    title_text = f" {title} "
    title_fill = content_width - len(title_text)
    title_left = title_fill // 2
    title_right = title_fill - title_left
    title_sep = f"╠{'═' * title_left}{title_text}{'═' * title_right}╣"

    header_sep = f"╠{'═' * (left_width + 2)}╪{'═' * (right_width + 2)}╣"
    light_row_sep = f"╟{'─' * (left_width + 2)}┼{'─' * (right_width + 2)}╢"

    lines = [top]

    for logo_line in FRUTON_LOGO.splitlines():
        lines.append(f"║{logo_line.center(content_width)}║")

    lines.extend(
        [
            title_sep,
            f"║ {left_header.ljust(left_width)} │ {right_header.rjust(right_width)} ║",
            header_sep,
        ]
    )

    for index, (left, right) in enumerate(rows):
        lines.append(f"║ {left.ljust(left_width)} │ {right.rjust(right_width)} ║")
        if index != len(rows) - 1:
            lines.append(light_row_sep)

    lines.append(bottom)
    return "\n".join(lines)


def _print_summary(summary: dict[str, int]) -> None:
    print()

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

    print(_build_summary_box(title="FRUTON Summary", rows=rows))
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


pipeline_runner.init_runner(
    _append_fruton_log,
    force_field=DEFAULT_GROMACS_FORCE_FIELD,
    water_model=DEFAULT_GROMACS_WATER_MODEL,
    ph=FRUTON_PROTONATION_PH,
)


class TerminalProgress:
    """Render compact runner progress without adding a UI dependency.

    The progress line is intentionally simple: it uses a carriage return to
    redraw the current terminal line instead of depending on Rich, curses, or
    tqdm. This keeps the runner usable in normal terminals and copied logs. The
    bar is based on estimated file-changing work units and is updated whenever a
    mutating pipeline call finishes. External module stdout/stderr is captured
    into the FRUTON log so the progress line is not destroyed by verbose module
    output.
    """

    def __init__(self, *, total_steps: int) -> None:
        self.total_steps = total_steps
        self.current_step = 0
        self.current_step_name = "startup"
        self.current_label = "initialising"
        self.current_pdb_id = "-"
        self.total_work = 1
        self.completed_work = 0
        self.started_at = time.monotonic()
        self.enabled = sys.stdout.isatty() and not _verbose_screen_enabled()

    def configure(self, *, total_work: int) -> None:
        self.total_work = max(int(total_work), 1)
        self.completed_work = min(self.completed_work, self.total_work)
        self.render()

    def start_step(self, *, step_number: int, step_name: str) -> None:
        self.current_step = step_number
        self.current_step_name = step_name
        self.current_label = step_name
        self.current_pdb_id = "-"
        self.clear()
        print(f"\n┏━ FRUTON step {step_number:02d}/{self.total_steps:02d} · {step_name}")
        self.render()

    def advance(self, *, label: str, pdb_id: str = "-", status: str = "done") -> None:
        self.completed_work = min(self.completed_work + 1, self.total_work)
        self.current_label = f"{label} · {status}"
        self.current_pdb_id = pdb_id or "-"
        self.render()

    def set_status(self, *, label: str, pdb_id: str = "-", status: str = "running") -> None:
        self.current_label = f"{label} · {status}"
        self.current_pdb_id = pdb_id or "-"
        self.render()

    def finish(self) -> None:
        self.completed_work = self.total_work
        self.current_label = "finished"
        self.current_pdb_id = "-"
        self.render()
        if self.enabled:
            print()

    def clear(self) -> None:
        if not self.enabled:
            return
        width = shutil.get_terminal_size((120, 20)).columns
        print("\r" + " " * (width - 1) + "\r", end="", flush=True)

    def render(self) -> None:
        if not self.enabled:
            return

        width = shutil.get_terminal_size((120, 20)).columns
        fraction = min(max(self.completed_work / max(self.total_work, 1), 0.0), 1.0)
        percent = 100.0 * fraction
        bar_width = 28
        filled = int(round(bar_width * fraction))
        bar = "█" * filled + "░" * (bar_width - filled)
        elapsed = int(time.monotonic() - self.started_at)

        line = (
            f"FRUTON [{bar}] {percent:5.1f}% | "
            f"work {self.completed_work}/{self.total_work} | "
            f"step {self.current_step:02d}/{self.total_steps:02d} | "
            f"{self.current_pdb_id} | {self.current_label} | {elapsed}s"
        )
        if len(line) >= width:
            line = line[: width - 4] + "..."
        print("\r" + line.ljust(width - 1), end="", flush=True)


class _CapturedOutput:
    """Context manager that redirects noisy module output to memory."""

    def __init__(self) -> None:
        self.stdout = io.StringIO()
        self.stderr = io.StringIO()

    def __enter__(self) -> "_CapturedOutput":
        self._stdout_cm = contextlib.redirect_stdout(self.stdout)
        self._stderr_cm = contextlib.redirect_stderr(self.stderr)
        self._stdout_cm.__enter__()
        self._stderr_cm.__enter__()
        return self

    def __exit__(self, exc_type: Any, exc: Any, tb: Any) -> None:
        self._stderr_cm.__exit__(exc_type, exc, tb)
        self._stdout_cm.__exit__(exc_type, exc, tb)

    def captured_lines(self) -> list[str]:
        lines: list[str] = []
        stdout_text = self.stdout.getvalue().strip()
        stderr_text = self.stderr.getvalue().strip()
        if stdout_text:
            lines.append("[captured stdout]")
            lines.extend(stdout_text.splitlines())
        if stderr_text:
            lines.append("[captured stderr]")
            lines.extend(stderr_text.splitlines())
        return lines


def _verbose_screen_enabled() -> bool:
    return os.environ.get("FRUTON_VERBOSE", "").strip().lower() in {"1", "true", "yes", "on"}


_PROGRESS = TerminalProgress(total_steps=17)


def _estimate_initial_work_units(record_count: int) -> int:
    """Estimate mutating work units for the current run.

    The estimate is deliberately practical rather than chemically exact. The
    runner knows the number of records after Step 4, but it does not yet know
    how many retained variants each protein will have after gap detection and
    filler. The estimate therefore counts the normal file-changing stages per
    record and reserves two late work units per record for protonation and final
    prepared-structure writing. The progress bar is finalized to 100% at the end
    of the run even if some records are skipped.
    """

    fixed_units = 4 + 2
    per_record_units = 10  # +1 for nonstd_residue_params (step 15)
    return fixed_units + record_count * per_record_units


def _screen_step(step_number: int, step_name: str) -> None:
    _PROGRESS.start_step(step_number=step_number, step_name=step_name)


def _screen_item(message: str) -> None:
    _append_fruton_log("screen", [message])
    if _verbose_screen_enabled():
        _PROGRESS.clear()
        print(f"┃ {message}")
    _PROGRESS.render()


def _screen_notice(message: str) -> None:
    _PROGRESS.clear()
    print(f"┃ {message}")
    _PROGRESS.render()


def _screen_error(message: str) -> None:
    _PROGRESS.clear()
    print(f"┃ ERROR · {message}")
    _PROGRESS.render()


def _run_mutating_call(
    *,
    log_title: str,
    work_label: str,
    screen_pdb_id: str,
    func: Any,
    **kwargs: Any,
) -> Any:
    """Run one file-changing function while keeping terminal output compact."""

    _PROGRESS.set_status(label=work_label, pdb_id=screen_pdb_id, status="running")
    with _CapturedOutput() as captured:
        result = func(**kwargs)

    captured_lines = captured.captured_lines()
    if captured_lines:
        _append_fruton_log(log_title, captured_lines)

    _PROGRESS.advance(label=work_label, pdb_id=screen_pdb_id, status="done")
    return result



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
    _print_dependency_report()
    _screen_notice("starting pipeline")
    _append_fruton_log(
        "dependency_check",
        FRUTON_DEPENDENCY_REPORT.to_log_lines(),
    )
    _append_fruton_log(
        "run_pipeline:start",
        [
            f"project_root_dir      : {PROJECT_ROOT_DIR}",
            f"protein_data_dir      : {protein_data_dir}",
            f"pdb_ids_csv_path      : {pdb_ids_csv_path}",
            f"pipeline_json         : {pipeline_json_path}",
            f"pipeline_xlsx         : {pipeline_xlsx_path}",
            "protein_unit_policy   : representative unit required",
            "ligand_policy         : split components from representative unit, not full multimer",
            "variant_policy       : GROMACS/prepared variants are built before final audit decision",
            "sanitize             : runs before protonation inputs and again on prepared variants",
            "parameter_audit       : runs on sanitized prepared variants",
            "internal_capping      : skipped for gaps; temporary fragment route",
            f"protonation_method    : {FRUTON_PROTONATION_METHOD}",
            f"protonation_ph_compat : {FRUTON_PROTONATION_PH}",
            f"gromacs_force_field   : {DEFAULT_GROMACS_FORCE_FIELD}",
            f"gromacs_water_model   : {DEFAULT_GROMACS_WATER_MODEL}",
        ],
    )

    _screen_step(1, "pdb_sync")
    try:
        _run_mutating_call(log_title="step_1:pdb_sync:captured_output", work_label="pdb_sync", screen_pdb_id="-", func=sync_pdb_csv_and_directories, protein_data_dir=protein_data_dir)
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
    _screen_notice(f"active pipeline records: {len(pipeline_record_list)}")
    _PROGRESS.configure(total_work=_estimate_initial_work_units(len(pipeline_record_list)))
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
            _run_mutating_call(
                log_title=f"step_5:fasta_files:{pdb_id}:captured_output",
                work_label="fasta_files",
                screen_pdb_id=pdb_id,
                func=create_fasta_files_for_pdb_directory,
                pdb_directory=pdb_dir,
                residue_range=residue_range,
            )

            pipeline_record[FASTA_FILES_DONE_COLUMN_NAME] = STATUS_SUCCESS
            pipeline_record[UNIPROT_ID_COLUMN_NAME] = _find_uniprot_id_for_protein(
                pdb_dir
            )
        except Exception as error:
            _screen_error(f"fasta_files failed for {pdb_id}: {error!r}")
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
            _run_mutating_call(
                log_title=f"step_6:sequence_alignment:{pdb_id}:captured_output",
                work_label="sequence_alignment",
                screen_pdb_id=pdb_id,
                func=run_alignments_for_pdb_directory,
                pdb_directory=pdb_dir,
            )

            pipeline_record[SEQUENCE_ALIGNMENT_DONE_COLUMN_NAME] = STATUS_SUCCESS

            if not str(pipeline_record.get(UNIPROT_ID_COLUMN_NAME, "")).strip():
                pipeline_record[UNIPROT_ID_COLUMN_NAME] = _find_uniprot_id_for_protein(
                    pdb_dir
                )
        except Exception as error:
            _screen_error(f"sequence_alignment failed for {pdb_id}: {error!r}")
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
            _screen_error(f"insertion_codes failed for {pdb_id}: input PDB not found")
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
            insertion_result = _run_mutating_call(
                log_title=f"step_7:insertion_codes:{pdb_id}:captured_output",
                work_label="insertion_codes",
                screen_pdb_id=pdb_id,
                func=process_pdb_for_delinsertion,
                input_pdb_path=input_pdb_path,
                output_pdb_path=output_pdb_path,
            )
        except Exception as error:
            _screen_error(f"insertion_codes failed for {pdb_id}: {error!r}")
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

    _screen_step(8, "representative_unit and component_split")
    for pipeline_record in pipeline_record_list:
        pdb_id = pipeline_record[PDB_ID_COLUMN_NAME]
        pdb_dir = Path(pipeline_record[PDB_DIRECTORY_COLUMN_NAME])
        components_dir = Path(pipeline_record[COMPONENTS_DIRECTORY_COLUMN_NAME])

        if pipeline_record.get(INSERTION_CODES_DONE_COLUMN_NAME, "") != STATUS_SUCCESS:
            _screen_item(f"representative_unit skipped for {pdb_id}")
            _clear_component_fields(pipeline_record)
            _clear_gap_fields(pipeline_record)
            _clear_downstream_after_gap_stage(pipeline_record)
            continue

        component_input_pdb_path = pdb_dir / f"{pdb_id}_delins.pdb"

        if not component_input_pdb_path.exists():
            _screen_item(f"representative_unit skipped for {pdb_id}: missing input")
            _clear_component_fields(pipeline_record)
            _clear_gap_fields(pipeline_record)
            _clear_downstream_after_gap_stage(pipeline_record)
            continue

        _screen_item(f"representative_unit -> {pdb_id}")

        try:
            representative_unit_path = _run_mutating_call(
                log_title=f"step_8:representative_unit:{pdb_id}:captured_output",
                work_label="representative_unit",
                screen_pdb_id=pdb_id,
                func=_write_representative_monomer_for_protein,
                pdb_id=pdb_id,
                pdb_dir=pdb_dir,
                input_pdb_path=component_input_pdb_path,
            )
            _screen_item(
                f"representative_unit -> {pdb_id}: {representative_unit_path.name}"
            )
        except Exception as error:
            _screen_error(f"representative_unit failed for {pdb_id}: {error!r}")
            _log_fruton_exception(f"step_8:representative_unit:{pdb_id}", error)
            _clear_component_fields(pipeline_record)
            _clear_gap_fields(pipeline_record)
            _clear_downstream_after_gap_stage(pipeline_record)
            continue

        try:
            _clear_component_output_files(pdb_id=pdb_id, components_dir=components_dir)

            component_summary = _run_mutating_call(
                log_title=f"step_8:component_split:{pdb_id}:captured_output",
                work_label="component_split",
                screen_pdb_id=pdb_id,
                func=split_pdb_components,
                pdb_path=representative_unit_path,
                output_dir=components_dir,
                protein_stem=pdb_id,
            )

            monomer_path = _run_mutating_call(
                log_title=f"step_8:monomer_copy:{pdb_id}:captured_output",
                work_label="monomer_copy",
                screen_pdb_id=pdb_id,
                func=_copy_representative_protein_component_to_monomer_path,
                pdb_id=pdb_id,
                pdb_dir=pdb_dir,
            )
        except Exception as error:
            _screen_error(f"component_split failed for {pdb_id}: {error!r}")
            _log_fruton_exception(f"step_8:component_split:{pdb_id}", error)
            _clear_component_fields(pipeline_record)
            _clear_gap_fields(pipeline_record)
            _clear_downstream_after_gap_stage(pipeline_record)
            continue

        pipeline_record[COMPONENTS_DIRECTORY_COLUMN_NAME] = str(
            component_summary["output_dir"]
        )
        _set_component_flags(pipeline_record, component_summary)


        _append_fruton_log(
            f"step_8:component_split:{pdb_id}",
            [
                f"representative_unit_path : {representative_unit_path}",
                f"monomer_protein_path     : {monomer_path}",
                f"components_dir           : {component_summary['output_dir']}",
                f"has_metals               : {pipeline_record[HAS_METALS_COLUMN_NAME]}",
                f"has_ligands              : {pipeline_record[HAS_LIGANDS_COLUMN_NAME]}",
                f"has_nonstandard_residues : {pipeline_record[HAS_NONSTANDARD_RESIDUES_COLUMN_NAME]}",
            ],
        )

        _screen_item(f"component_split -> {pdb_id}: from representative unit")
        _screen_item(f"monomer_selection -> {pdb_id}: {monomer_path.name}")

        try:
            metal_inventory_result = _run_mutating_call(
                log_title=f"step_8:metalls_check_input_inventory:{pdb_id}:captured_output",
                work_label="metalls_check:input_inventory",
                screen_pdb_id=pdb_id,
                func=_run_metals_check_for_input_inventory,
                pdb_id=pdb_id,
                pdb_dir=pdb_dir,
                monomer_path=Path(monomer_path),
            )
            _store_metal_check_result_in_record(
                pipeline_record=pipeline_record,
                metal_check_result=metal_inventory_result,
            )
        except Exception as error:
            _screen_error(f"input metal inventory failed for {pdb_id}: {error!r}")
            _log_fruton_exception(
                f"step_8:metalls_check_input_inventory:{pdb_id}",
                error,
            )

    _screen_step(9, "gap_detection")
    for pipeline_record in pipeline_record_list:
        pdb_id = pipeline_record[PDB_ID_COLUMN_NAME]
        pdb_dir = Path(pipeline_record[PDB_DIRECTORY_COLUMN_NAME])
        residue_range = str(pipeline_record.get(RANGE_COLUMN_NAME, "")).strip()

        if pipeline_record.get(INSERTION_CODES_DONE_COLUMN_NAME, "") != STATUS_SUCCESS:
            _screen_item(f"gap_detection skipped for {pdb_id}")
            _clear_gap_fields(pipeline_record)
            _clear_downstream_after_gap_stage(pipeline_record)
            continue

        try:
            gap_input_pdb_path = _select_default_protein_input_path(
                pdb_id=pdb_id,
                pdb_dir=pdb_dir,
            )
        except FileNotFoundError as error:
            _screen_item(f"gap_detection skipped for {pdb_id}: missing monomer")
            _log_fruton_exception(f"step_9:gap_detection:{pdb_id}", error)
            _clear_gap_fields(pipeline_record)
            _clear_downstream_after_gap_stage(pipeline_record)
            continue

        _screen_item(f"gap_detection -> {pdb_id} [representative unit protein]")

        try:
            if residue_range:
                gap_summary = _run_mutating_call(
                    log_title=f"step_9:gap_detection:{pdb_id}:captured_output",
                    work_label="gap_detection",
                    screen_pdb_id=pdb_id,
                    func=summarize_gaps,
                    pdb_path=gap_input_pdb_path,
                    residue_range=residue_range,
                )
            else:
                gap_summary = _run_mutating_call(
                    log_title=f"step_9:gap_detection:{pdb_id}:captured_output",
                    work_label="gap_detection",
                    screen_pdb_id=pdb_id,
                    func=summarize_gaps,
                    pdb_path=gap_input_pdb_path,
                )
        except Exception as error:
            _screen_error(f"gap_detection failed for {pdb_id}: {error!r}")
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

        blocks_complete_model, block_reason = _blocks_complete_model_generation(
            pipeline_record
        )

        if blocks_complete_model:
            _screen_item(f"filler skipped for {pdb_id}: {block_reason}")
            pipeline_record[FILLER_STATUS_COLUMN_NAME] = STATUS_SKIPPED
            pipeline_record[FILLER_MODEL_PATH_COLUMN_NAME] = ""
            pipeline_record[FILLER_MODEL_SOURCE_COLUMN_NAME] = ""
            pipeline_record[FILLER_DIRECTORY_COLUMN_NAME] = ""
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
            _screen_item(f"filler skipped for {pdb_id}: missing monomer/alignment input")
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

            filler_result = _run_mutating_call(
                log_title=f"step_10:filler:{pdb_id}:captured_output",
                work_label="filler",
                screen_pdb_id=pdb_id,
                func=run_filler_for_chain,
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
            _screen_error(f"filler failed for {pdb_id}: {error!r}")
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
    prepared_variant_results_by_pdb: dict[str, list[dict[str, str]]] = {}

    _screen_step(11, "GROMACS protonation + prepared variants")
    for pipeline_record in pipeline_record_list:
        pdb_id = pipeline_record[PDB_ID_COLUMN_NAME]
        pdb_dir = Path(pipeline_record[PDB_DIRECTORY_COLUMN_NAME])

        _clear_protonation_fields(pipeline_record)
        _clear_internal_capping_fields(pipeline_record)
        _clear_prepared_structure_fields(pipeline_record)

        has_metals = str(
            pipeline_record.get(HAS_METALS_COLUMN_NAME, "")
        ).strip().lower()

        if has_metals == "yes":
            _append_fruton_log(
                "step_11:metal_structure_allowed",
                [
                    f"pdb_id      : {pdb_id}",
                    "policy      : metal-containing structures are prepared, then checked by metalls_check",
                ],
            )

        n_gaps_text = str(pipeline_record.get(N_GAPS_COLUMN_NAME, "")).strip()
        if not n_gaps_text:
            _screen_item(f"prepared variants skipped for {pdb_id}: no gap data")
            pipeline_record[PREPARED_STRUCTURE_STATUS_COLUMN_NAME] = STATUS_FAILED
            pipeline_record[AVAILABLE_MODELS_COLUMN_NAME] = ""
            summary["proteins_failed"] += 1
            prepared_variant_results_by_pdb[pdb_id] = []
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
        successful_variant_result_list: list[dict[str, str]] = []

        for variant_plan in variant_plan_list:
            variant_label = str(variant_plan["variant_label"])
            model_source = str(variant_plan["model_source"])
            model_path = variant_plan["model_path"]
            model_path = model_path if isinstance(model_path, Path) else None

            _screen_item(f"protonation/prepared -> {pdb_id} [{variant_label}]")

            try:
                route_result = _run_mutating_call(
                    log_title=f"step_11:protonation_route:{pdb_id}:{variant_label}:captured_output",
                    work_label=f"protonation:{variant_label}",
                    screen_pdb_id=pdb_id,
                    func=_run_protonation_route_for_variant,
                    pdb_id=pdb_id,
                    pdb_dir=pdb_dir,
                    variant_label=variant_label,
                    model_source=model_source,
                    model_path=model_path,
                )
            except Exception as error:
                _screen_error(f"protonation route failed for {pdb_id} [{variant_label}]: {error!r}")
                _log_fruton_exception(
                    f"step_11:protonation_route:{pdb_id}:{variant_label}",
                    error,
                )
                continue

            if route_result["status"] != STATUS_SUCCESS:
                pipeline_record[INTERNAL_CAPPING_STATUS_COLUMN_NAME] = route_result.get(
                    "internal_capping_status",
                    pipeline_record.get(INTERNAL_CAPPING_STATUS_COLUMN_NAME, ""),
                )
                pipeline_record[INTERNAL_CAPPING_INPUT_PATH_COLUMN_NAME] = route_result.get(
                    "internal_capping_input_path",
                    pipeline_record.get(INTERNAL_CAPPING_INPUT_PATH_COLUMN_NAME, ""),
                )
                pipeline_record[INTERNAL_CAPPING_OUTPUT_PATH_COLUMN_NAME] = route_result.get(
                    "internal_capping_output_path",
                    pipeline_record.get(INTERNAL_CAPPING_OUTPUT_PATH_COLUMN_NAME, ""),
                )
                pipeline_record[PROTONATION_STATUS_COLUMN_NAME] = STATUS_FAILED
                pipeline_record[PROTONATION_INPUT_SOURCE_COLUMN_NAME] = route_result.get(
                    "protonation_input_source",
                    "",
                )
                pipeline_record[PROTONATION_INPUT_PATH_COLUMN_NAME] = route_result.get(
                    "protonation_input_path",
                    "",
                )
                pipeline_record[PROTONATION_OUTPUT_PATH_COLUMN_NAME] = route_result.get(
                    "protonation_output_path",
                    "",
                )
                continue

            pipeline_record[INTERNAL_CAPPING_STATUS_COLUMN_NAME] = route_result.get(
                "internal_capping_status",
                STATUS_SKIPPED,
            )
            pipeline_record[INTERNAL_CAPPING_INPUT_PATH_COLUMN_NAME] = route_result.get(
                "internal_capping_input_path",
                "",
            )
            pipeline_record[INTERNAL_CAPPING_OUTPUT_PATH_COLUMN_NAME] = route_result.get(
                "internal_capping_output_path",
                "",
            )

            pipeline_record[PROTONATION_STATUS_COLUMN_NAME] = STATUS_SUCCESS
            pipeline_record[PROTONATION_INPUT_SOURCE_COLUMN_NAME] = route_result[
                "protonation_input_source"
            ]
            pipeline_record[PROTONATION_INPUT_PATH_COLUMN_NAME] = route_result[
                "protonation_input_path"
            ]
            pipeline_record[PROTONATION_OUTPUT_PATH_COLUMN_NAME] = route_result[
                "protonation_output_path"
            ]

            final_protein_input_paths_text = route_result.get(
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
                    prepared_output_path = _run_mutating_call(
                        log_title=f"step_11:prepared_structure:{pdb_id}:{variant_label}:captured_output",
                        work_label=f"prepared:{variant_label}",
                        screen_pdb_id=pdb_id,
                        func=_build_prepared_structure_for_variant,
                        pdb_id=pdb_id,
                        pdb_dir=pdb_dir,
                        variant_label=variant_label,
                        final_protein_input_path=final_input_path_list[0],
                        final_protein_input_paths=None,
                    )
                else:
                    prepared_output_path = _run_mutating_call(
                        log_title=f"step_11:prepared_structure:{pdb_id}:{variant_label}:captured_output",
                        work_label=f"prepared:{variant_label}",
                        screen_pdb_id=pdb_id,
                        func=_build_prepared_structure_for_variant,
                        pdb_id=pdb_id,
                        pdb_dir=pdb_dir,
                        variant_label=variant_label,
                        final_protein_input_path=None,
                        final_protein_input_paths=final_input_path_list,
                    )
            except Exception as error:
                _screen_error(f"prepared_structure failed for {pdb_id} [{variant_label}]: {error!r}")
                _log_fruton_exception(
                    f"step_11:prepared_structure:{pdb_id}:{variant_label}",
                    error,
                )
                continue

            try:
                metal_check_result = _run_mutating_call(
                    log_title=f"step_11:metalls_check:{pdb_id}:{variant_label}:captured_output",
                    work_label=f"metalls_check:{variant_label}",
                    screen_pdb_id=pdb_id,
                    func=_run_metals_check_for_prepared_variant,
                    pdb_id=pdb_id,
                    pdb_dir=pdb_dir,
                    variant_label=variant_label,
                    prepared_output_path=Path(prepared_output_path),
                )
                _store_metal_check_result_in_record(
                    pipeline_record=pipeline_record,
                    metal_check_result=metal_check_result,
                )
            except Exception as error:
                _screen_error(f"metalls_check failed for {pdb_id} [{variant_label}]: {error!r}")
                _log_fruton_exception(
                    f"step_11:metalls_check:{pdb_id}:{variant_label}",
                    error,
                )
                metal_check_result = None

            successful_variant_label_list.append(variant_label)
            successful_variant_result_list.append(
                {
                    "variant_label": variant_label,
                    "model_source": model_source,
                    "final_protein_input_paths": final_protein_input_paths_text,
                    "prepared_output_path": str(prepared_output_path),
                    "metals.status": str(getattr(metal_check_result, "status", "")) if metal_check_result is not None else "failed",
                    "metals.ion_type": str(getattr(metal_check_result, "ion_type_text", "")) if metal_check_result is not None else "",
                    "metals.class": str(getattr(metal_check_result, "metal_class", "")) if metal_check_result is not None else "",
                    "metals.geometry_found": str(getattr(metal_check_result, "found_geometry_text", "")) if metal_check_result is not None else "",
                    "metals.geometry_probable": str(getattr(metal_check_result, "probable_geometry_text", "")) if metal_check_result is not None else "",
                    "metals.geometry_match": str(getattr(metal_check_result, "geometry_match", "")) if metal_check_result is not None else "",
                    "metals.parameter_reference": str(getattr(metal_check_result, "parameter_reference_text", "")) if metal_check_result is not None else "",
                    "metals.model_ready": "yes" if (metal_check_result is None or bool(getattr(metal_check_result, "model_metal_ready", False))) else "no",
                    "metals.check_log_path": str(getattr(metal_check_result, "log_path", "")) if metal_check_result is not None else "",
                    "metals.check_manifest_path": str(getattr(metal_check_result, "manifest_path", "")) if metal_check_result is not None else "",
                }
            )

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

        prepared_variant_results_by_pdb[pdb_id] = successful_variant_result_list

        if successful_variant_label_list:
            representative_variant_result = _choose_representative_successful_variant(
                successful_variant_result_list
            )

            summary["proteins_with_prepared_output"] += 1
            pipeline_record[PREPARED_STRUCTURE_STATUS_COLUMN_NAME] = STATUS_SUCCESS
            pipeline_record[PREPARED_STRUCTURE_VARIANT_COLUMN_NAME] = ",".join(
                successful_variant_label_list
            )
            pipeline_record[PREPARED_STRUCTURE_PROTEIN_INPUT_PATH_COLUMN_NAME] = (
                representative_variant_result["final_protein_input_paths"]
            )
            pipeline_record[PREPARED_STRUCTURE_OUTPUT_PATH_COLUMN_NAME] = (
                representative_variant_result["prepared_output_path"]
            )
            _update_available_models_field(pipeline_record)
        else:
            prepared_variant_results_by_pdb[pdb_id] = []
            summary["proteins_failed"] += 1
            pipeline_record[INTERNAL_CAPPING_STATUS_COLUMN_NAME] = STATUS_SKIPPED
            if not pipeline_record[PROTONATION_STATUS_COLUMN_NAME]:
                pipeline_record[PROTONATION_STATUS_COLUMN_NAME] = STATUS_FAILED
            pipeline_record[PREPARED_STRUCTURE_STATUS_COLUMN_NAME] = STATUS_FAILED
            pipeline_record[AVAILABLE_MODELS_COLUMN_NAME] = ""

    _screen_step(12, "sanitize variants")
    sanitized_variant_results_by_pdb: dict[str, list[dict[str, str]]] = {}

    for pipeline_record in pipeline_record_list:
        pdb_id = pipeline_record[PDB_ID_COLUMN_NAME]
        pdb_dir = Path(pipeline_record[PDB_DIRECTORY_COLUMN_NAME])
        prepared_variant_result_list = prepared_variant_results_by_pdb.get(pdb_id, [])

        if not prepared_variant_result_list:
            sanitized_variant_results_by_pdb[pdb_id] = []
            _screen_item(f"sanitize variants skipped for {pdb_id}: no prepared variants")
            continue

        sanitized_variant_result_list: list[dict[str, str]] = []

        for variant_result in prepared_variant_result_list:
            variant_label = variant_result["variant_label"]
            prepared_output_path = Path(variant_result["prepared_output_path"])

            if not prepared_output_path.exists() or prepared_output_path.stat().st_size == 0:
                _screen_item(
                    f"sanitize variants skipped for {pdb_id} [{variant_label}]: missing prepared output"
                )
                _append_fruton_log(
                    f"step_12:sanitize_variant_skipped:{pdb_id}:{variant_label}",
                    [
                        f"pdb_id               : {pdb_id}",
                        f"variant_label        : {variant_label}",
                        f"prepared_output_path : {prepared_output_path}",
                        "reason               : prepared output missing or empty",
                    ],
                )
                continue

            try:
                sanitize_result = _run_mutating_call(
                    log_title=f"step_12:sanitize_variant:{pdb_id}:{variant_label}:captured_output",
                    work_label=f"sanitize_variant:{variant_label}",
                    screen_pdb_id=pdb_id,
                    func=_run_sanitize_for_prepared_variant,
                    pdb_id=pdb_id,
                    pdb_dir=pdb_dir,
                    variant_label=variant_label,
                    prepared_output_path=prepared_output_path,
                )
            except Exception as error:
                _screen_error(
                    f"sanitize variant failed for {pdb_id} [{variant_label}]: {error!r}"
                )
                _log_fruton_exception(
                    f"step_12:sanitize_variant:{pdb_id}:{variant_label}",
                    error,
                )
                continue

            sanitized_output_path = Path(
                str(getattr(sanitize_result, "output_pdb_path", ""))
            )
            if not sanitized_output_path.exists() or sanitized_output_path.stat().st_size == 0:
                _screen_item(
                    f"sanitize variants failed for {pdb_id} [{variant_label}]: no sanitized output"
                )
                continue

            sanitized_variant_result = dict(variant_result)
            sanitized_variant_result["sanitized_audit_input_path"] = str(
                sanitized_output_path
            )
            sanitized_variant_result["sanitize_status"] = str(
                getattr(sanitize_result, "status", "")
            )
            sanitized_variant_result["sanitize_log_path"] = str(
                getattr(sanitize_result, "log_path", "")
            )
            sanitized_variant_result_list.append(sanitized_variant_result)

            _screen_item(
                f"sanitize variants -> {pdb_id} [{variant_label}]: {sanitized_output_path.name}"
            )

        sanitized_variant_results_by_pdb[pdb_id] = sanitized_variant_result_list

    _screen_step(13, "parameter_audit of variants")

    for pipeline_record in pipeline_record_list:
        pdb_id = pipeline_record[PDB_ID_COLUMN_NAME]
        pdb_dir = Path(pipeline_record[PDB_DIRECTORY_COLUMN_NAME])
        sanitized_variant_result_list = sanitized_variant_results_by_pdb.get(pdb_id, [])

        if not sanitized_variant_result_list:
            _apply_accepted_variant_audit_decision(
                pipeline_record=pipeline_record,
                accepted_variant_result_list=[],
            )
            _screen_item(f"parameter_audit variants skipped for {pdb_id}: no sanitized variants")
            continue

        accepted_variant_result_list: list[dict[str, str]] = []

        for variant_result in sanitized_variant_result_list:
            variant_label = variant_result["variant_label"]
            sanitized_input_path = Path(variant_result["sanitized_audit_input_path"])

            try:
                parameter_audit_result = _run_mutating_call(
                    log_title=f"step_13:parameter_audit_variant:{pdb_id}:{variant_label}:captured_output",
                    work_label=f"parameter_audit:{variant_label}",
                    screen_pdb_id=pdb_id,
                    func=_run_parameter_audit_for_variant,
                    pdb_id=pdb_id,
                    pdb_dir=pdb_dir,
                    variant_label=variant_label,
                    input_pdb_path=sanitized_input_path,
                )
            except Exception as error:
                _screen_error(
                    f"parameter_audit variant failed for {pdb_id} [{variant_label}]: {error!r}"
                )
                _log_fruton_exception(
                    f"step_13:parameter_audit_variant:{pdb_id}:{variant_label}",
                    error,
                )
                continue

            accepted = _variant_audit_accepts_end_model(parameter_audit_result)
            if accepted and not _variant_metal_check_accepts_end_model(variant_result):
                accepted = False
                _screen_notice(
                    f"parameter_audit variants -> {pdb_id} [{variant_label}]: metal parameters required"
                )
            _append_parameter_audit_variant_summary(
                pdb_id=pdb_id,
                variant_label=variant_label,
                audit_result=parameter_audit_result,
                accepted=accepted,
            )

            if accepted:
                accepted_variant_result_list.append(variant_result)
                _screen_item(f"parameter_audit variants -> {pdb_id} [{variant_label}]: accepted")
            elif getattr(parameter_audit_result, "requires_qm_parameters", False):
                _screen_notice(
                    f"parameter_audit variants -> {pdb_id} [{variant_label}]: unsupported residue parameters required"
                )
            elif getattr(parameter_audit_result, "requires_repair", False):
                _screen_notice(
                    f"parameter_audit variants -> {pdb_id} [{variant_label}]: input repair still required"
                )
            else:
                _screen_notice(
                    f"parameter_audit variants -> {pdb_id} [{variant_label}]: rejected"
                )

        _apply_accepted_variant_audit_decision(
            pipeline_record=pipeline_record,
            accepted_variant_result_list=accepted_variant_result_list,
        )

    _recompute_prepared_summary_from_state(
        summary=summary,
        pipeline_record_list=pipeline_record_list,
    )

    _screen_step(14, "metall_params")
    for pipeline_record in pipeline_record_list:
        pdb_id = pipeline_record[PDB_ID_COLUMN_NAME]
        pdb_dir = Path(pipeline_record[PDB_DIRECTORY_COLUMN_NAME])

        has_metals = str(pipeline_record.get(HAS_METALS_COLUMN_NAME, "")).strip().lower()
        if has_metals != "yes":
            _screen_item(f"metall_params skipped for {pdb_id}: no metals")
            continue

        try:
            metall_params_result = _run_mutating_call(
                log_title=f"step_14:metall_params:{pdb_id}:captured_output",
                work_label=f"metall_params:{pdb_id}",
                screen_pdb_id=pdb_id,
                func=_run_metall_params_for_protein,
                pdb_id=pdb_id,
                pdb_dir=pdb_dir,
            )
        except Exception as error:
            _screen_error(f"metall_params failed for {pdb_id}: {error!r}")
            _log_fruton_exception(f"step_14:metall_params:{pdb_id}", error)
            pipeline_record[METALL_PARAMS_STATUS_COLUMN_NAME] = STATUS_FAILED
            continue

        pipeline_record[METALL_PARAMS_STATUS_COLUMN_NAME] = str(
            metall_params_result.get("status", "")
        )
        pipeline_record[METALL_PARAMS_SITE_COUNT_COLUMN_NAME] = str(
            metall_params_result.get("transition_metal_site_count", 0)
        )
        pipeline_record[METALL_PARAMS_MANIFEST_PATH_COLUMN_NAME] = str(
            metall_params_result.get("manifest_path", "") or ""
        )
        site_count = metall_params_result.get("transition_metal_site_count", 0)
        status = metall_params_result.get("status", "")
        _screen_item(f"metall_params -> {pdb_id}: {status} ({site_count} transition-metal site(s))")

    _screen_step(15, "nonstd_residue_params")
    for pipeline_record in pipeline_record_list:
        pdb_id = pipeline_record[PDB_ID_COLUMN_NAME]
        pdb_dir = Path(pipeline_record[PDB_DIRECTORY_COLUMN_NAME])

        has_nonstd = str(pipeline_record.get(HAS_NONSTANDARD_RESIDUES_COLUMN_NAME, "")).strip().lower()
        if has_nonstd != "yes":
            _screen_item(f"nonstd_residue_params skipped for {pdb_id}: no non-standard residues")
            continue

        try:
            nonstd_result = _run_mutating_call(
                log_title=f"step_15:nonstd_residue_params:{pdb_id}:captured_output",
                work_label=f"nonstd_residue_params:{pdb_id}",
                screen_pdb_id=pdb_id,
                func=_run_nonstd_residue_params_for_protein,
                pdb_id=pdb_id,
                pdb_dir=pdb_dir,
            )
        except Exception as error:
            _screen_error(f"nonstd_residue_params failed for {pdb_id}: {error!r}")
            _log_fruton_exception(f"step_15:nonstd_residue_params:{pdb_id}", error)
            pipeline_record[NONSTD_RESIDUE_PARAMS_STATUS_COLUMN_NAME] = STATUS_FAILED
            continue

        pipeline_record[NONSTD_RESIDUE_PARAMS_STATUS_COLUMN_NAME] = str(
            nonstd_result.get("status", "")
        )
        pipeline_record[NONSTD_RESIDUE_PARAMS_N_RESIDUES_COLUMN_NAME] = str(
            nonstd_result.get("n_residues", 0)
        )
        pipeline_record[NONSTD_RESIDUE_PARAMS_MANIFEST_PATH_COLUMN_NAME] = str(
            nonstd_result.get("manifest_path", "") or ""
        )
        n_residues = nonstd_result.get("n_residues", 0)
        status = nonstd_result.get("status", "")
        _screen_item(f"nonstd_residue_params -> {pdb_id}: {status} ({n_residues} non-standard residue(s))")

    _screen_step(16, "save_pipeline_json")
    _run_mutating_call(log_title="step_16:save_pipeline_json:captured_output", work_label="save_pipeline_json", screen_pdb_id="-", func=save_pipeline_table, protein_record_list=pipeline_record_list, json_path=pipeline_json_path)
    _screen_step(17, "write_pipeline_xlsx")
    _run_mutating_call(log_title="step_17:write_pipeline_xlsx:captured_output", work_label="write_pipeline_xlsx", screen_pdb_id="-", func=write_pipeline_to_xlsx, protein_record_list=pipeline_record_list, output_path=pipeline_xlsx_path)

    _screen_notice(f"pipeline JSON written: {pipeline_json_path}")
    _screen_notice(f"pipeline XLSX written: {pipeline_xlsx_path}")

    _append_fruton_log(
        "run_pipeline:finished",
        [
            f"pipeline_json_path : {pipeline_json_path}",
            f"pipeline_xlsx_path : {pipeline_xlsx_path}",
            f"summary            : {summary}",
        ],
    )

    _print_summary(summary)
    _PROGRESS.finish()
    print("FRUTON finished")


def main() -> None:
    run_pipeline()


if __name__ == "__main__":
    main()