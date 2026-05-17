#!/usr/bin/env python3

from __future__ import annotations

import contextlib
import io
import os
import re
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
from stack_protein_preparation.cap import cap_internal_gaps_for_variant
from stack_protein_preparation.fragment_split import split_gap_variant_structure
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
from stack_protein_preparation.monomer import write_single_representative_monomer_unit
from stack_protein_preparation.metalls_check import (
    run_metals_check_for_model,
    run_metals_inventory_for_structure,
)
from stack_protein_preparation.pdb_components import split_pdb_components
from stack_protein_preparation.pdb_sync import (
    read_pdb_records_from_csv,
    sync_pdb_csv_and_directories,
)
from stack_protein_preparation.parameter_audit import (
    audit_protein_residue_parameters,
)
from stack_protein_preparation.pipeline_state import (
    ALIGNMENT_DIRECTORY_COLUMN_NAME,
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
    METALS_CHECK_LOG_PATH_COLUMN_NAME,
    METALS_CHECK_MANIFEST_PATH_COLUMN_NAME,
    METALS_CHECK_STATUS_COLUMN_NAME,
    METALS_CLASS_COLUMN_NAME,
    METALS_GEOMETRY_FOUND_COLUMN_NAME,
    METALS_GEOMETRY_MATCH_COLUMN_NAME,
    METALS_GEOMETRY_PROBABLE_COLUMN_NAME,
    METALS_ION_TYPE_COLUMN_NAME,
    METALS_MODEL_READY_COLUMN_NAME,
    METALS_PARAMETER_REFERENCE_COLUMN_NAME,
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
from stack_protein_preparation.sanitize import sanitize_pdb_for_gromacs
from stack_protein_preparation.standard_residue_repair import (
    repair_standard_residue_heavy_atoms_for_gromacs,
)
from stack_protein_preparation.sequence_alignment import (
    run_alignments_for_pdb_directory,
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


_PROGRESS = TerminalProgress(total_steps=15)


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
    per_record_units = 9
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


def _build_raw_protein_path(pdb_id: str, pdb_dir: Path) -> Path:
    return pdb_dir / "components" / f"{pdb_id}_protein.pdb"


def _build_monomer_protein_path(pdb_id: str, pdb_dir: Path) -> Path:
    return pdb_dir / "components" / f"{pdb_id}_protein_monomer.pdb"


def _build_representative_unit_path(pdb_id: str, pdb_dir: Path) -> Path:
    return pdb_dir / "components" / f"{pdb_id}_representative_unit.pdb"


def _clear_component_output_files(pdb_id: str, components_dir: Path) -> None:
    """Remove stale component files before splitting the representative unit.

    The component splitter writes category files such as protein, ligand, water,
    and metal outputs, but an old output can survive when a new representative
    unit no longer contains that category. That is dangerous here because a
    ligand file from the full multimer can be reintroduced later by prepared
    structure assembly. Removing known component outputs before splitting keeps
    all component files synchronized with the representative unit.
    """

    for suffix in (
        "protein",
        "protein_monomer",
        "ligand",
        "water",
        "metals",
        "nonstandard_residues",
    ):
        (components_dir / f"{pdb_id}_{suffix}.pdb").unlink(missing_ok=True)


def _extract_monomer_output_path(result: Any, fallback_path: Path) -> Path:
    if isinstance(result, str | Path):
        return Path(result)

    output_pdb = getattr(result, "output_pdb", None)
    if output_pdb is not None:
        return Path(output_pdb)

    return fallback_path


def _write_representative_monomer_for_protein(
    pdb_id: str,
    pdb_dir: Path,
    input_pdb_path: Path,
) -> Path:
    """Write the representative unit from the full processed PDB input.

    The monomer selector must run on ``<PDBID>_delins.pdb`` rather than on
    ``components/<PDBID>_protein.pdb``. The protein component file is already
    stripped of ligands, metals, cofactors, and other HETATM records, so a
    representative monomer written from that file can never retain ligand records
    per selected chain. Running on the full processed input lets ``monomer.py``
    keep non-water HETATM records assigned to the selected representative chains.

    This helper only writes ``components/<PDBID>_representative_unit.pdb``. Step 8
    then splits that representative unit into protein, ligand, water, and metal
    component files, so downstream prepared-structure assembly cannot reintroduce
    ligand records from discarded duplicate chains.
    """

    components_dir = Path(pdb_dir) / "components"
    representative_unit_path = _build_representative_unit_path(
        pdb_id=pdb_id,
        pdb_dir=pdb_dir,
    )

    if not input_pdb_path.exists():
        raise FileNotFoundError(
            f"Representative-unit input file not found: {input_pdb_path}"
        )

    components_dir.mkdir(parents=True, exist_ok=True)

    result = write_single_representative_monomer_unit(
        input_pdb_path,
        representative_unit_path,
        keep_non_protein_hetero=True,
        keep_waters=False,
        overwrite=True,
    )

    result_path = _extract_monomer_output_path(
        result=result,
        fallback_path=representative_unit_path,
    )

    if not result_path.exists() or result_path.stat().st_size == 0:
        raise RuntimeError(
            f"Representative unit output was not written or is empty: {result_path}"
        )

    representative_chain_ids = getattr(result, "representative_chain_ids", ())
    represented_chain_groups = getattr(result, "represented_chain_groups", ())

    _append_fruton_log(
        "monomer_selection:success",
        [
            f"pdb_id                     : {pdb_id}",
            f"representative_unit_input  : {input_pdb_path}",
            f"representative_unit_output : {result_path}",
            f"representative_chain_ids   : {'/'.join(representative_chain_ids)}",
            f"represented_chain_groups   : {represented_chain_groups}",
            "ligand_policy              : keep non-water HETATM only for selected representative chains",
        ],
    )

    return result_path


def _copy_representative_protein_component_to_monomer_path(
    pdb_id: str,
    pdb_dir: Path,
) -> Path:
    """Copy the split representative protein component to the monomer path.

    Downstream stages already use ``components/<PDBID>_protein_monomer.pdb`` as
    the strict default protein input. After Step 8 is changed to split the full
    representative unit, ``components/<PDBID>_protein.pdb`` already contains only
    the representative protein chains. This helper mirrors that file to the
    historical monomer path so the rest of the runner does not need a broader
    path-interface change.
    """

    split_protein_path = _build_raw_protein_path(pdb_id=pdb_id, pdb_dir=pdb_dir)
    monomer_protein_path = _build_monomer_protein_path(pdb_id=pdb_id, pdb_dir=pdb_dir)

    if not split_protein_path.exists() or split_protein_path.stat().st_size == 0:
        raise RuntimeError(
            f"Representative protein component was not written: {split_protein_path}"
        )

    monomer_protein_path.write_text(
        split_protein_path.read_text(encoding="utf-8"),
        encoding="utf-8",
    )

    if not monomer_protein_path.exists() or monomer_protein_path.stat().st_size == 0:
        raise RuntimeError(
            f"Representative monomer protein output was not written: {monomer_protein_path}"
        )

    return monomer_protein_path



def _run_parameter_audit_for_protein(
    *,
    pdb_id: str,
    pdb_dir: Path,
    input_pdb_path: Path,
) -> Any:
    """Run the residue parameter audit for one representative protein.

    The audit writes a per-protein JSON manifest and a per-protein text log
    under the normal FRUTON output tree. It does not change the input PDB and it
    does not generate Gaussian files. Its role is only to decide whether the
    current GROMACS force field accepts the concrete representative protein,
    whether the input needs repair, or whether unsupported residues should be
    sent to the separate residue-parameter workflow.
    """

    audit_manifest_path = pdb_dir / "parameter_audit.json"
    audit_log_path = pdb_dir.parent / "logs" / pdb_id / "parameter_audit.log"
    audit_log_path.parent.mkdir(parents=True, exist_ok=True)

    return audit_protein_residue_parameters(
        input_pdb_path=input_pdb_path,
        pdb_id=pdb_id,
        force_field=DEFAULT_GROMACS_FORCE_FIELD,
        water_model=DEFAULT_GROMACS_WATER_MODEL,
        manifest_path=audit_manifest_path,
        log_path=audit_log_path,
    )


def _append_parameter_audit_summary(
    *,
    pdb_id: str,
    audit_result: Any,
) -> None:
    """Write the compact parameter-audit decision into the FRUTON log."""

    _append_fruton_log(
        f"step_12:parameter_audit:{pdb_id}",
        [
            f"status                       : {getattr(audit_result, 'status', '')}",
            f"input_pdb_path               : {getattr(audit_result, 'input_pdb_path', '')}",
            f"manifest_path                : {getattr(audit_result, 'manifest_path', '')}",
            f"log_path                     : {getattr(audit_result, 'log_path', '')}",
            f"gromacs_probe_success        : {getattr(audit_result, 'gromacs_probe_success', '')}",
            f"can_use_current_force_field  : {getattr(audit_result, 'can_use_current_force_field', '')}",
            f"requires_repair              : {getattr(audit_result, 'requires_repair', '')}",
            f"requires_external_parameters : {getattr(audit_result, 'requires_external_parameters', '')}",
            f"requires_qm_parameters       : {getattr(audit_result, 'requires_qm_parameters', '')}",
            f"requires_metal_parameters    : {getattr(audit_result, 'requires_metal_parameters', '')}",
            f"supported_nonstandard        : {', '.join(getattr(audit_result, 'supported_nonstandard_residue_names', ())) or '(none)'}",
            f"unsupported_residues         : {', '.join(getattr(audit_result, 'unsupported_residue_names', ())) or '(none)'}",
        ],
    )

def _select_default_protein_input_path(
    pdb_id: str,
    pdb_dir: Path,
) -> Path:
    """Return the original representative monomer protein input.

    Gap detection and fragment splitting must see the original residue numbering
    from ``components/<PDBID>_protein_monomer.pdb``. Sanitized copies are used
    later for parameter audit and GROMACS-facing protonation inputs, but they do
    not replace the canonical monomer path used for structural decisions.
    """

    monomer_protein_path = _build_monomer_protein_path(pdb_id=pdb_id, pdb_dir=pdb_dir)

    if monomer_protein_path.exists() and monomer_protein_path.stat().st_size > 0:
        return monomer_protein_path

    raw_protein_path = _build_raw_protein_path(pdb_id=pdb_id, pdb_dir=pdb_dir)

    raise FileNotFoundError(
        f"Representative monomer protein input not found for {pdb_id}: "
        f"{monomer_protein_path}. Raw protein path is not used as a fallback: "
        f"{raw_protein_path}."
    )

def _blocks_complete_model_generation(
    pipeline_record: dict[str, str],
) -> tuple[bool, str]:
    """Return whether this protein should skip complete-model generation.

    Complete-model generation means MODELLER or AlphaFold reconstruction. The
    runner only blocks this step for metal-containing structures because the
    current reconstruction route cannot preserve metal coordination chemistry
    reliably. Non-standard residues are audited separately by
    ``parameter_audit.py`` and must not suppress MODELLER/AlphaFold logs by
    themselves. Unsupported residue chemistry should be handled later through
    the residue-parameter workflow rather than by hiding the reconstruction
    attempt.
    """

    has_metals = str(
        pipeline_record.get(HAS_METALS_COLUMN_NAME, "")
    ).strip().lower()

    if has_metals == "yes":
        return True, "metal-containing structure"

    return False, ""


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
    _ = components_dir

    try:
        return _select_default_protein_input_path(pdb_id=pdb_id, pdb_dir=pdb_dir)
    except FileNotFoundError:
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


def _clear_internal_capping_fields(pipeline_record: dict[str, str]) -> None:
    pipeline_record[INTERNAL_CAPPING_STATUS_COLUMN_NAME] = ""
    pipeline_record[INTERNAL_CAPPING_INPUT_PATH_COLUMN_NAME] = ""
    pipeline_record[INTERNAL_CAPPING_OUTPUT_PATH_COLUMN_NAME] = ""


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
    _clear_internal_capping_fields(pipeline_record)
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
        return _select_default_protein_input_path(pdb_id=pdb_id, pdb_dir=pdb_dir), "protein"

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
    """Run the FRUTON-default GROMACS protonation route for one variant.

    The runner treats ``gmx pdb2gmx -ignh`` as the protein protonation and
    GROMACS naming step. The output is passed directly into prepared-structure
    assembly because it should remain in the naming convention produced by the
    selected GROMACS force field. The historical pH value is still passed as
    compatibility metadata because the protonation module records it, but it is
    not applied directly by GROMACS.
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



def _build_sanitized_protonation_path(
    *,
    pdb_id: str,
    pdb_dir: Path,
    variant_label: str,
) -> Path:
    """Return a stable temporary sanitized input path for protonation.

    Sanitized structures are implementation artefacts of the GROMACS route.
    They are kept under ``components/.tmp_sanitized_protonation`` so they do not
    look like final FRUTON products, but they remain inspectable when a route
    fails. The variant label is normalized before it is used as a filename so
    fragment labels and future model-source labels cannot create nested paths.
    """

    safe_variant_label = re.sub(r"[^A-Za-z0-9_.-]+", "_", variant_label).strip("_")
    if not safe_variant_label:
        safe_variant_label = "variant"

    sanitize_dir = pdb_dir / "components" / ".tmp_sanitized_protonation"
    sanitize_dir.mkdir(parents=True, exist_ok=True)
    return sanitize_dir / f"{pdb_id}_{safe_variant_label}_sanitized.pdb"


def _build_sanitize_log_path(
    *,
    pdb_id: str,
    pdb_dir: Path,
    variant_label: str,
) -> Path:
    """Return the per-protein sanitizer log path for one protonation input."""

    safe_variant_label = re.sub(r"[^A-Za-z0-9_.-]+", "_", variant_label).strip("_")
    if not safe_variant_label:
        safe_variant_label = "variant"

    return pdb_dir.parent / "logs" / pdb_id / f"sanitize_{safe_variant_label}.log"


def _summarize_sanitize_issues(sanitize_result: Any) -> str:
    """Return a compact one-line summary of sanitizer issues."""

    issue_list = list(getattr(sanitize_result, "issues", ()) or ())
    if not issue_list:
        return "<none>"

    return " | ".join(
        f"{getattr(issue, 'severity', '')}:{getattr(issue, 'code', '')}"
        for issue in issue_list
    )



def _build_standard_residue_repair_path(
    *,
    pdb_id: str,
    pdb_dir: Path,
    variant_label: str,
) -> Path:
    """Return the temporary MODELLER repair path for one protonation input."""

    safe_variant_label = "".join(
        character if character.isalnum() or character in {"_", "-"} else "_"
        for character in variant_label.strip()
    ).strip("_")
    repair_dir = pdb_dir / "components" / ".tmp_standard_residue_repair"
    repair_dir.mkdir(parents=True, exist_ok=True)
    return repair_dir / f"{pdb_id}_{safe_variant_label}_standard_repaired.pdb"


def _build_standard_residue_repair_log_path(
    *,
    pdb_id: str,
    pdb_dir: Path,
    variant_label: str,
) -> Path:
    """Return the repair log path for one protonation input."""

    safe_variant_label = "".join(
        character if character.isalnum() or character in {"_", "-"} else "_"
        for character in variant_label.strip()
    ).strip("_")
    return pdb_dir.parent / "logs" / pdb_id / f"standard_residue_repair_{safe_variant_label}.log"


def _run_standard_residue_repair_for_protonation_input(
    *,
    pdb_id: str,
    pdb_dir: Path,
    variant_label: str,
    input_pdb: Path,
) -> Any:
    """Run minimal MODELLER missing-heavy-atom repair before sanitizer."""

    repair_result = repair_standard_residue_heavy_atoms_for_gromacs(
        input_pdb_path=input_pdb,
        output_pdb_path=_build_standard_residue_repair_path(
            pdb_id=pdb_id,
            pdb_dir=pdb_dir,
            variant_label=variant_label,
        ),
        log_path=_build_standard_residue_repair_log_path(
            pdb_id=pdb_id,
            pdb_dir=pdb_dir,
            variant_label=variant_label,
        ),
        force_field=DEFAULT_GROMACS_FORCE_FIELD,
        allow_nonstandard_residues=False,
    )
    _append_fruton_log(
        "standard_residue_repair:result",
        [
            f"pdb_id                            : {pdb_id}",
            f"variant_label                     : {variant_label}",
            f"input_pdb                         : {input_pdb}",
            f"used_output_path                  : {getattr(repair_result, 'used_output_path', '')}",
            f"repair_log_path                   : {getattr(repair_result, 'log_path', '')}",
            f"status                            : {getattr(repair_result, 'status', '')}",
            f"repair_attempted                  : {getattr(repair_result, 'repair_attempted', '')}",
            f"repair_success                    : {getattr(repair_result, 'repair_success', '')}",
            f"missing_heavy_atom_count_before   : {getattr(repair_result, 'missing_heavy_atom_count_before', '')}",
            f"missing_heavy_atom_count_after    : {getattr(repair_result, 'missing_heavy_atom_count_after', '')}",
            f"message                           : {getattr(repair_result, 'message', '')}",
        ],
    )
    return repair_result

def _run_sanitize_for_protonation_input(
    *,
    pdb_id: str,
    pdb_dir: Path,
    variant_label: str,
    input_pdb: Path,
) -> Any:
    """Sanitize one structure before passing it to ``pdb2gmx``.

    FRUTON treats sanitizing as a deterministic cleanup gate immediately before
    protonation. The gate selects alternate locations, normalizes occupancies,
    removes non-polymer heterogens for the protein-only route, and writes a
    clean PDB copy. It does not build atoms, replace PTMs, add loops, add
    hydrogens, or generate non-standard residue parameters; unresolved chemistry
    is left visible for parameter audit and the final GROMACS check.
    """

    repair_result = _run_standard_residue_repair_for_protonation_input(
        pdb_id=pdb_id,
        pdb_dir=pdb_dir,
        variant_label=variant_label,
        input_pdb=input_pdb,
    )
    repair_input_path = Path(str(getattr(repair_result, "used_output_path", input_pdb)))

    output_pdb_path = _build_sanitized_protonation_path(
        pdb_id=pdb_id,
        pdb_dir=pdb_dir,
        variant_label=variant_label,
    )
    log_path = _build_sanitize_log_path(
        pdb_id=pdb_id,
        pdb_dir=pdb_dir,
        variant_label=variant_label,
    )

    sanitize_result = sanitize_pdb_for_gromacs(
        input_pdb_path=repair_input_path,
        output_pdb_path=output_pdb_path,
        pdb_id=pdb_id,
        variant_label=variant_label,
        force_field=DEFAULT_GROMACS_FORCE_FIELD,
        log_path=log_path,
        fail_on_missing_heavy_atoms=False,
        fail_on_nonstandard_residues=False,
    )

    _append_fruton_log(
        "protonation_sanitize:result",
        [
            f"pdb_id                            : {pdb_id}",
            f"variant_label                     : {variant_label}",
            f"input_pdb_path                    : {input_pdb}",
            f"sanitized_output_path             : {getattr(sanitize_result, 'output_pdb_path', '')}",
            f"sanitize_log_path                 : {getattr(sanitize_result, 'log_path', '')}",
            f"status                            : {getattr(sanitize_result, 'status', '')}",
            f"can_run_gromacs                   : {getattr(sanitize_result, 'can_run_gromacs', '')}",
            f"atom_count_in                     : {getattr(sanitize_result, 'atom_count_in', '')}",
            f"atom_count_out                    : {getattr(sanitize_result, 'atom_count_out', '')}",
            f"selected_altloc_count             : {getattr(sanitize_result, 'selected_altloc_count', '')}",
            f"normalized_occupancy_count        : {getattr(sanitize_result, 'normalized_occupancy_count', '')}",
            f"missing_heavy_atom_count          : {getattr(sanitize_result, 'missing_heavy_atom_count', '')}",
            f"nonstandard_residue_names         : {', '.join(getattr(sanitize_result, 'nonstandard_residue_names', ())) or '(none)'}",
            f"issues                            : {_summarize_sanitize_issues(sanitize_result)}",
        ],
    )

    return sanitize_result


def _build_sanitized_monomer_protein_path(
    pdb_id: str,
    pdb_dir: Path,
) -> Path:
    """Return the stable sanitized representative protein path.

    This file is created after representative-unit selection and component
    splitting. It is the earliest safe place to sanitize the protein-only input
    because doing it before representative-unit/component splitting would either
    drop ligands/metals too early or sanitize the wrong multimeric structure.
    Downstream default protein selection prefers this file when it exists.
    """

    return pdb_dir / "components" / f"{pdb_id}_protein_monomer_sanitized.pdb"


def _build_representative_sanitize_log_path(
    pdb_id: str,
    pdb_dir: Path,
) -> Path:
    """Return the sanitizer log path for the representative protein input."""

    return pdb_dir.parent / "logs" / pdb_id / "sanitize_representative_protein.log"


def _run_sanitize_for_representative_protein(
    *,
    pdb_id: str,
    pdb_dir: Path,
    input_pdb: Path,
) -> Any:
    """Sanitize the representative monomer protein before gap/filler/protonation.

    This visible sanitize gate runs after gap detection and filler decisions but
    before parameter audit. It writes a cleaned protein-only copy for audit and
    later GROMACS-facing steps without replacing the original monomer file used
    for gap detection and fragment splitting.
    """

    output_pdb_path = _build_sanitized_monomer_protein_path(
        pdb_id=pdb_id,
        pdb_dir=pdb_dir,
    )
    log_path = _build_representative_sanitize_log_path(
        pdb_id=pdb_id,
        pdb_dir=pdb_dir,
    )

    sanitize_result = sanitize_pdb_for_gromacs(
        input_pdb_path=input_pdb,
        output_pdb_path=output_pdb_path,
        pdb_id=pdb_id,
        variant_label="representative_protein",
        force_field=DEFAULT_GROMACS_FORCE_FIELD,
        log_path=log_path,
        fail_on_missing_heavy_atoms=False,
        fail_on_nonstandard_residues=False,
    )

    _append_fruton_log(
        "step_11:sanitize:result",
        [
            f"pdb_id                            : {pdb_id}",
            f"input_pdb_path                    : {input_pdb}",
            f"sanitized_output_path             : {getattr(sanitize_result, 'output_pdb_path', '')}",
            f"sanitize_log_path                 : {getattr(sanitize_result, 'log_path', '')}",
            f"status                            : {getattr(sanitize_result, 'status', '')}",
            f"can_run_gromacs                   : {getattr(sanitize_result, 'can_run_gromacs', '')}",
            f"atom_count_in                     : {getattr(sanitize_result, 'atom_count_in', '')}",
            f"atom_count_out                    : {getattr(sanitize_result, 'atom_count_out', '')}",
            f"selected_altloc_count             : {getattr(sanitize_result, 'selected_altloc_count', '')}",
            f"normalized_occupancy_count        : {getattr(sanitize_result, 'normalized_occupancy_count', '')}",
            f"missing_heavy_atom_count          : {getattr(sanitize_result, 'missing_heavy_atom_count', '')}",
            f"nonstandard_residue_names         : {', '.join(getattr(sanitize_result, 'nonstandard_residue_names', ())) or '(none)'}",
            f"issues                            : {_summarize_sanitize_issues(sanitize_result)}",
        ],
    )

    return sanitize_result


def _get_sanitize_result_output_path(
    sanitize_result: Any,
    fallback_path: Path,
) -> Path:
    """Return a usable sanitizer output path or the original fallback path."""

    raw_output_path = str(getattr(sanitize_result, "output_pdb_path", "")).strip()
    if raw_output_path:
        output_path = Path(raw_output_path)
        if output_path.exists() and output_path.stat().st_size > 0:
            return output_path

    return fallback_path



def _build_sanitized_prepared_variant_path(
    *,
    pdb_id: str,
    pdb_dir: Path,
    variant_label: str,
) -> Path:
    """Return the sanitized protein-only audit path for one prepared variant.

    Step 12 sanitizes the produced prepared-structure model rather than the
    starting representative protein. The output is kept in an implementation
    directory because it is an audit artefact: it removes ligands, waters,
    metals, and other non-polymer heterogens so ``parameter_audit.py`` can test
    the actual produced protein model in a protein-only GROMACS probe.
    """

    safe_variant_label = re.sub(r"[^A-Za-z0-9_.-]+", "_", variant_label).strip("_")
    if not safe_variant_label:
        safe_variant_label = "variant"

    audit_dir = pdb_dir / "components" / ".tmp_sanitized_variant_audit"
    audit_dir.mkdir(parents=True, exist_ok=True)
    return audit_dir / f"{pdb_id}_{safe_variant_label}_prepared_sanitized.pdb"


def _build_prepared_variant_sanitize_log_path(
    *,
    pdb_id: str,
    pdb_dir: Path,
    variant_label: str,
) -> Path:
    """Return the sanitizer log path for a produced prepared variant."""

    safe_variant_label = re.sub(r"[^A-Za-z0-9_.-]+", "_", variant_label).strip("_")
    if not safe_variant_label:
        safe_variant_label = "variant"

    return pdb_dir.parent / "logs" / pdb_id / f"sanitize_prepared_{safe_variant_label}.log"


def _run_sanitize_for_prepared_variant(
    *,
    pdb_id: str,
    pdb_dir: Path,
    variant_label: str,
    prepared_output_path: Path,
) -> Any:
    """Sanitize one produced prepared variant before final parameter audit.

    This is intentionally later than gap detection, filler, and protonation. The
    sanitizer therefore tests the model FRUTON actually produced, not the raw
    representative monomer. It is still sanitize-only: it selects alternate
    locations, normalizes occupancies, drops non-polymer heterogens for the
    protein-only audit route, and writes a cleaned audit copy without repairing
    missing atoms or replacing non-standard residues.
    """

    output_pdb_path = _build_sanitized_prepared_variant_path(
        pdb_id=pdb_id,
        pdb_dir=pdb_dir,
        variant_label=variant_label,
    )
    log_path = _build_prepared_variant_sanitize_log_path(
        pdb_id=pdb_id,
        pdb_dir=pdb_dir,
        variant_label=variant_label,
    )

    sanitize_result = sanitize_pdb_for_gromacs(
        input_pdb_path=prepared_output_path,
        output_pdb_path=output_pdb_path,
        pdb_id=pdb_id,
        variant_label=f"prepared_{variant_label}",
        force_field=DEFAULT_GROMACS_FORCE_FIELD,
        log_path=log_path,
        fail_on_missing_heavy_atoms=False,
        fail_on_nonstandard_residues=False,
    )

    _append_fruton_log(
        "step_12:sanitize_prepared_variant:result",
        [
            f"pdb_id                            : {pdb_id}",
            f"variant_label                     : {variant_label}",
            f"prepared_output_path              : {prepared_output_path}",
            f"sanitized_output_path             : {getattr(sanitize_result, 'output_pdb_path', '')}",
            f"sanitize_log_path                 : {getattr(sanitize_result, 'log_path', '')}",
            f"status                            : {getattr(sanitize_result, 'status', '')}",
            f"can_run_gromacs                   : {getattr(sanitize_result, 'can_run_gromacs', '')}",
            f"atom_count_in                     : {getattr(sanitize_result, 'atom_count_in', '')}",
            f"atom_count_out                    : {getattr(sanitize_result, 'atom_count_out', '')}",
            f"selected_altloc_count             : {getattr(sanitize_result, 'selected_altloc_count', '')}",
            f"normalized_occupancy_count        : {getattr(sanitize_result, 'normalized_occupancy_count', '')}",
            f"missing_heavy_atom_count          : {getattr(sanitize_result, 'missing_heavy_atom_count', '')}",
            f"nonstandard_residue_names         : {', '.join(getattr(sanitize_result, 'nonstandard_residue_names', ())) or '(none)'}",
            f"issues                            : {_summarize_sanitize_issues(sanitize_result)}",
        ],
    )

    return sanitize_result


def _build_variant_parameter_audit_manifest_path(
    *,
    pdb_id: str,
    pdb_dir: Path,
    variant_label: str,
) -> Path:
    """Return the JSON manifest path for one prepared-variant audit."""

    safe_variant_label = re.sub(r"[^A-Za-z0-9_.-]+", "_", variant_label).strip("_")
    if not safe_variant_label:
        safe_variant_label = "variant"

    return pdb_dir / f"parameter_audit_{safe_variant_label}.json"


def _build_variant_parameter_audit_log_path(
    *,
    pdb_id: str,
    pdb_dir: Path,
    variant_label: str,
) -> Path:
    """Return the text log path for one prepared-variant audit."""

    safe_variant_label = re.sub(r"[^A-Za-z0-9_.-]+", "_", variant_label).strip("_")
    if not safe_variant_label:
        safe_variant_label = "variant"

    return pdb_dir.parent / "logs" / pdb_id / f"parameter_audit_{safe_variant_label}.log"


def _run_parameter_audit_for_variant(
    *,
    pdb_id: str,
    pdb_dir: Path,
    variant_label: str,
    input_pdb_path: Path,
) -> Any:
    """Run parameter audit on the sanitized protein-only copy of a final variant.

    This audit is the final end-model decision gate. Unlike the earlier
    representative-protein audit, it examines the exact produced model variant
    after FRUTON has completed gap handling, filler selection, protonation, and
    prepared-structure assembly.
    """

    audit_manifest_path = _build_variant_parameter_audit_manifest_path(
        pdb_id=pdb_id,
        pdb_dir=pdb_dir,
        variant_label=variant_label,
    )
    audit_log_path = _build_variant_parameter_audit_log_path(
        pdb_id=pdb_id,
        pdb_dir=pdb_dir,
        variant_label=variant_label,
    )
    audit_log_path.parent.mkdir(parents=True, exist_ok=True)

    return audit_protein_residue_parameters(
        input_pdb_path=input_pdb_path,
        pdb_id=pdb_id,
        force_field=DEFAULT_GROMACS_FORCE_FIELD,
        water_model=DEFAULT_GROMACS_WATER_MODEL,
        manifest_path=audit_manifest_path,
        log_path=audit_log_path,
    )


def _append_parameter_audit_variant_summary(
    *,
    pdb_id: str,
    variant_label: str,
    audit_result: Any,
    accepted: bool,
) -> None:
    """Write compact final variant-audit decision into the FRUTON log."""

    _append_fruton_log(
        f"step_13:parameter_audit_variant:{pdb_id}:{variant_label}",
        [
            f"pdb_id                       : {pdb_id}",
            f"variant_label                : {variant_label}",
            f"accepted_as_end_model        : {accepted}",
            f"status                       : {getattr(audit_result, 'status', '')}",
            f"input_pdb_path               : {getattr(audit_result, 'input_pdb_path', '')}",
            f"manifest_path                : {getattr(audit_result, 'manifest_path', '')}",
            f"log_path                     : {getattr(audit_result, 'log_path', '')}",
            f"gromacs_probe_success        : {getattr(audit_result, 'gromacs_probe_success', '')}",
            f"can_use_current_force_field  : {getattr(audit_result, 'can_use_current_force_field', '')}",
            f"requires_repair              : {getattr(audit_result, 'requires_repair', '')}",
            f"requires_external_parameters : {getattr(audit_result, 'requires_external_parameters', '')}",
            f"requires_qm_parameters       : {getattr(audit_result, 'requires_qm_parameters', '')}",
            f"requires_metal_parameters    : {getattr(audit_result, 'requires_metal_parameters', '')}",
            f"supported_nonstandard        : {', '.join(getattr(audit_result, 'supported_nonstandard_residue_names', ())) or '(none)'}",
            f"unsupported_residues         : {', '.join(getattr(audit_result, 'unsupported_residue_names', ())) or '(none)'}",
        ],
    )




def _run_metals_check_for_input_inventory(
    *,
    pdb_id: str,
    pdb_dir: Path,
    monomer_path: Path,
) -> Any:
    """Record metal identity from the component-split input structure.

    This inventory runs before any prepared variant exists. It is deliberately
    independent from final-model acceptance: it fills the workbook with metal
    identity and standard-ion policy even for structures that later fail
    protonation or parameter audit. Crystal waters and the separated metals
    component are included so coordination analysis has the same context that
    the prepared-variant check will later use.
    """

    result = run_metals_inventory_for_structure(
        pdb_id=pdb_id,
        protein_dir=pdb_dir,
        model_pdb_path=monomer_path,
        metals_pdb_path=pdb_dir / "components" / f"{pdb_id}_metals.pdb",
        water_pdb_path=pdb_dir / "components" / f"{pdb_id}_water.pdb",
        output_dir=pdb_dir / "components" / "metals_check" / "input_inventory",
        use_chimera=False,
    )
    _append_fruton_log(
        "metalls_check:input_inventory",
        [
            f"pdb_id                    : {pdb_id}",
            f"monomer_path              : {monomer_path}",
            f"status                    : {getattr(result, 'status', '')}",
            f"metal_class               : {getattr(result, 'metal_class', '')}",
            f"ion_type_text             : {getattr(result, 'ion_type_text', '')}",
            f"parameter_reference_text  : {getattr(result, 'parameter_reference_text', '')}",
            f"found_geometry_text       : {getattr(result, 'found_geometry_text', '')}",
            f"probable_geometry_text    : {getattr(result, 'probable_geometry_text', '')}",
            f"geometry_match            : {getattr(result, 'geometry_match', '')}",
            f"model_metal_ready         : {getattr(result, 'model_metal_ready', '')}",
            f"contacts_tsv_path         : {getattr(result, 'contacts_tsv_path', '')}",
            f"manifest_path             : {getattr(result, 'manifest_path', '')}",
            f"log_path                  : {getattr(result, 'log_path', '')}",
        ],
    )
    return result


def _run_metals_check_for_prepared_variant(
    *,
    pdb_id: str,
    pdb_dir: Path,
    variant_label: str,
    prepared_output_path: Path,
) -> Any:
    """Run metal classification on one produced prepared variant."""

    result = run_metals_check_for_model(
        pdb_id=pdb_id,
        protein_dir=pdb_dir,
        variant_label=variant_label,
        model_pdb_path=prepared_output_path,
        metals_pdb_path=pdb_dir / "components" / f"{pdb_id}_metals.pdb",
        water_pdb_path=pdb_dir / "components" / f"{pdb_id}_water.pdb",
        output_dir=pdb_dir / "components" / "metals_check" / variant_label,
        use_chimera=True,
    )
    _append_fruton_log(
        "metalls_check:prepared_variant",
        [
            f"pdb_id                    : {pdb_id}",
            f"variant_label             : {variant_label}",
            f"prepared_output_path      : {prepared_output_path}",
            f"status                    : {getattr(result, 'status', '')}",
            f"metal_class               : {getattr(result, 'metal_class', '')}",
            f"ion_type_text             : {getattr(result, 'ion_type_text', '')}",
            f"parameter_reference_text  : {getattr(result, 'parameter_reference_text', '')}",
            f"found_geometry_text       : {getattr(result, 'found_geometry_text', '')}",
            f"probable_geometry_text    : {getattr(result, 'probable_geometry_text', '')}",
            f"geometry_match            : {getattr(result, 'geometry_match', '')}",
            f"model_metal_ready         : {getattr(result, 'model_metal_ready', '')}",
            f"contacts_tsv_path         : {getattr(result, 'contacts_tsv_path', '')}",
            f"manifest_path             : {getattr(result, 'manifest_path', '')}",
            f"log_path                  : {getattr(result, 'log_path', '')}",
        ],
    )
    return result


def _store_metal_check_result_in_record(
    *,
    pipeline_record: dict[str, str],
    metal_check_result: Any,
) -> None:
    """Store the compact metal check summary in the flat pipeline state."""

    pipeline_record[METALS_CHECK_STATUS_COLUMN_NAME] = str(getattr(metal_check_result, "status", ""))
    pipeline_record[METALS_ION_TYPE_COLUMN_NAME] = str(getattr(metal_check_result, "ion_type_text", ""))
    pipeline_record[METALS_CLASS_COLUMN_NAME] = str(getattr(metal_check_result, "metal_class", ""))
    pipeline_record[METALS_PARAMETER_REFERENCE_COLUMN_NAME] = str(getattr(metal_check_result, "parameter_reference_text", ""))
    pipeline_record[METALS_GEOMETRY_FOUND_COLUMN_NAME] = str(getattr(metal_check_result, "found_geometry_text", ""))
    pipeline_record[METALS_GEOMETRY_PROBABLE_COLUMN_NAME] = str(getattr(metal_check_result, "probable_geometry_text", ""))
    pipeline_record[METALS_GEOMETRY_MATCH_COLUMN_NAME] = str(getattr(metal_check_result, "geometry_match", ""))
    pipeline_record[METALS_MODEL_READY_COLUMN_NAME] = "yes" if bool(getattr(metal_check_result, "model_metal_ready", False)) else "no"
    pipeline_record[METALS_CHECK_LOG_PATH_COLUMN_NAME] = str(getattr(metal_check_result, "log_path", ""))
    pipeline_record[METALS_CHECK_MANIFEST_PATH_COLUMN_NAME] = str(getattr(metal_check_result, "manifest_path", ""))


def _variant_metal_check_accepts_end_model(variant_result: dict[str, str]) -> bool:
    """Return whether the variant metal check allows final-model acceptance."""

    value = str(variant_result.get("metals.model_ready", "")).strip().lower()
    return value in {"", "yes", "true", "1"}

def _variant_audit_accepts_end_model(audit_result: Any) -> bool:
    """Return whether a sanitized prepared variant can remain an end model."""

    if not bool(getattr(audit_result, "can_use_current_force_field", False)):
        return False
    if bool(getattr(audit_result, "requires_repair", False)):
        return False
    if bool(getattr(audit_result, "requires_external_parameters", False)):
        return False
    if bool(getattr(audit_result, "requires_qm_parameters", False)):
        return False
    if bool(getattr(audit_result, "requires_metal_parameters", False)):
        return False
    return True


def _clear_prepared_variant_output_paths(pipeline_record: dict[str, str]) -> None:
    """Clear prepared-variant path fields before writing accepted end models."""

    pipeline_record[PREPARED_STRUCTURE_VARIANT_COLUMN_NAME] = ""
    pipeline_record[PREPARED_STRUCTURE_PROTEIN_INPUT_PATH_COLUMN_NAME] = ""
    pipeline_record[PREPARED_STRUCTURE_OUTPUT_PATH_COLUMN_NAME] = ""
    pipeline_record[PREPARED_GAPS_OUTPUT_PATH_COLUMN_NAME] = ""
    pipeline_record[PREPARED_MODELLER_OUTPUT_PATH_COLUMN_NAME] = ""
    pipeline_record[PREPARED_ALPHAFOLD_OUTPUT_PATH_COLUMN_NAME] = ""
    pipeline_record[AVAILABLE_MODELS_COLUMN_NAME] = ""


def _apply_accepted_variant_audit_decision(
    *,
    pipeline_record: dict[str, str],
    accepted_variant_result_list: list[dict[str, str]],
) -> None:
    """Update final prepared fields after variant-level parameter audit.

    Step 11 may produce several tentative prepared variants. Step 13 is the
    final end-model decision: only variants whose sanitized prepared model passes
    parameter audit remain visible as final prepared outputs in the pipeline
    state and workbook.
    """

    _clear_prepared_variant_output_paths(pipeline_record)

    if not accepted_variant_result_list:
        pipeline_record[PREPARED_STRUCTURE_STATUS_COLUMN_NAME] = STATUS_FAILED
        return

    representative_variant_result = _choose_representative_successful_variant(
        accepted_variant_result_list
    )

    accepted_label_list = [
        result["variant_label"] for result in accepted_variant_result_list
    ]
    pipeline_record[PREPARED_STRUCTURE_STATUS_COLUMN_NAME] = STATUS_SUCCESS
    pipeline_record[PREPARED_STRUCTURE_VARIANT_COLUMN_NAME] = ",".join(accepted_label_list)
    pipeline_record[PREPARED_STRUCTURE_PROTEIN_INPUT_PATH_COLUMN_NAME] = (
        representative_variant_result.get("final_protein_input_paths", "")
    )
    pipeline_record[PREPARED_STRUCTURE_OUTPUT_PATH_COLUMN_NAME] = (
        representative_variant_result.get("prepared_output_path", "")
    )

    for result in accepted_variant_result_list:
        variant_label = result.get("variant_label", "")
        model_source = result.get("model_source", "")
        prepared_output_path = result.get("prepared_output_path", "")

        if variant_label == "gaps":
            pipeline_record[PREPARED_GAPS_OUTPUT_PATH_COLUMN_NAME] = prepared_output_path
        elif variant_label in {"best_complete", "large_gap_complete"}:
            if model_source == "modeller":
                pipeline_record[PREPARED_MODELLER_OUTPUT_PATH_COLUMN_NAME] = prepared_output_path
            elif model_source == "alphafold":
                pipeline_record[PREPARED_ALPHAFOLD_OUTPUT_PATH_COLUMN_NAME] = prepared_output_path

    for metal_column_name in [
        METALS_CHECK_STATUS_COLUMN_NAME,
        METALS_ION_TYPE_COLUMN_NAME,
        METALS_CLASS_COLUMN_NAME,
        METALS_PARAMETER_REFERENCE_COLUMN_NAME,
        METALS_GEOMETRY_FOUND_COLUMN_NAME,
        METALS_GEOMETRY_PROBABLE_COLUMN_NAME,
        METALS_GEOMETRY_MATCH_COLUMN_NAME,
        METALS_MODEL_READY_COLUMN_NAME,
        METALS_CHECK_LOG_PATH_COLUMN_NAME,
        METALS_CHECK_MANIFEST_PATH_COLUMN_NAME,
    ]:
        if representative_variant_result.get(metal_column_name, ""):
            pipeline_record[metal_column_name] = representative_variant_result.get(metal_column_name, "")

    _update_available_models_field(pipeline_record)


def _recompute_prepared_summary_from_state(
    *,
    summary: dict[str, int],
    pipeline_record_list: list[dict[str, str]],
) -> None:
    """Recompute final prepared-output counters after variant audit decisions."""

    summary["single_best_available_written"] = 0
    summary["gaps_variant_written"] = 0
    summary["modeller_complete_written"] = 0
    summary["alphafold_complete_written"] = 0
    summary["proteins_with_prepared_output"] = 0
    summary["proteins_failed"] = 0

    for pipeline_record in pipeline_record_list:
        prepared_status = str(
            pipeline_record.get(PREPARED_STRUCTURE_STATUS_COLUMN_NAME, "")
        ).strip()
        variant_text = str(
            pipeline_record.get(PREPARED_STRUCTURE_VARIANT_COLUMN_NAME, "")
        ).strip()
        variant_set = {
            token.strip()
            for token in variant_text.split(",")
            if token.strip()
        }

        if prepared_status == STATUS_SUCCESS and variant_set:
            summary["proteins_with_prepared_output"] += 1
        else:
            summary["proteins_failed"] += 1
            continue

        if "single" in variant_set:
            summary["single_best_available_written"] += 1
        if "gaps" in variant_set:
            summary["gaps_variant_written"] += 1
        if str(pipeline_record.get(PREPARED_MODELLER_OUTPUT_PATH_COLUMN_NAME, "")).strip():
            summary["modeller_complete_written"] += 1
        if str(pipeline_record.get(PREPARED_ALPHAFOLD_OUTPUT_PATH_COLUMN_NAME, "")).strip():
            summary["alphafold_complete_written"] += 1
def _build_fragment_variant_label(
    base_variant_label: str,
    fragment_index: int,
) -> str:
    return f"{base_variant_label}_fragment_{fragment_index:02d}"


def _run_linear_protonation_route_for_variant(
    pdb_id: str,
    pdb_dir: Path,
    variant_label: str,
    model_source: str,
    model_path: Path | None,
) -> dict[str, str]:
    _append_fruton_log(
        "protonation_route_linear:start",
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
            "route              : monomer/model -> sanitize -> gmx pdb2gmx protonation -> prepared_structure",
        ],
    )

    selected_input_path, input_source_label = _choose_variant_structure_input_path(
        pdb_id=pdb_id,
        pdb_dir=pdb_dir,
        model_source=model_source,
        model_path=model_path,
    )

    try:
        sanitize_result = _run_sanitize_for_protonation_input(
            pdb_id=pdb_id,
            pdb_dir=pdb_dir,
            variant_label=variant_label,
            input_pdb=selected_input_path,
        )
    except Exception as error:
        _log_fruton_exception(
            f"protonation_route_linear:sanitize_exception:{pdb_id}:{variant_label}",
            error,
        )
        return {
            "status": STATUS_FAILED,
            "variant_label": variant_label,
            "structure_input_path": str(selected_input_path),
            "structure_input_source": input_source_label,
            "internal_capping_status": STATUS_SKIPPED,
            "internal_capping_input_path": "",
            "internal_capping_output_path": "",
            "protonation_input_source": "sanitize_exception",
            "protonation_input_path": str(selected_input_path),
            "protonation_output_path": "",
            "final_protein_input_paths": "",
        }

    sanitized_input_path = Path(str(getattr(sanitize_result, "output_pdb_path", "")))
    if not bool(getattr(sanitize_result, "can_run_gromacs", False)):
        _append_fruton_log(
            "protonation_route_linear:sanitize_failed",
            [
                f"pdb_id                 : {pdb_id}",
                f"variant_label          : {variant_label}",
                f"input_path             : {selected_input_path}",
                f"sanitized_input_path   : {sanitized_input_path}",
                f"sanitize_log_path      : {getattr(sanitize_result, 'log_path', '')}",
                f"status                 : {getattr(sanitize_result, 'status', '')}",
                f"issues                 : {_summarize_sanitize_issues(sanitize_result)}",
            ],
        )
        return {
            "status": STATUS_FAILED,
            "variant_label": variant_label,
            "structure_input_path": str(selected_input_path),
            "structure_input_source": input_source_label,
            "internal_capping_status": STATUS_SKIPPED,
            "internal_capping_input_path": "",
            "internal_capping_output_path": "",
            "protonation_input_source": "sanitize_failed",
            "protonation_input_path": str(sanitized_input_path),
            "protonation_output_path": "",
            "final_protein_input_paths": "",
        }

    protonation_result = _run_fruton_protonation_for_variant(
        pdb_id=pdb_id,
        pdb_dir=pdb_dir,
        variant_label=variant_label,
        input_pdb=sanitized_input_path,
        input_source=input_source_label,
    )

    if not protonation_result.get("protonation_success", False):
        _append_fruton_log(
            "protonation_route_linear:protonation_failed",
            [
                f"pdb_id                    : {pdb_id}",
                f"variant_label             : {variant_label}",
                f"input_path                : {selected_input_path}",
                f"sanitized_input_path      : {sanitized_input_path}",
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
            "internal_capping_status": STATUS_SKIPPED,
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
            "final_protein_input_paths": "",
        }

    protonation_output_path = Path(str(protonation_result["protonation_output_path"]))

    return {
        "status": STATUS_SUCCESS,
        "variant_label": variant_label,
        "structure_input_path": str(selected_input_path),
        "structure_input_source": input_source_label,
        "internal_capping_status": STATUS_SKIPPED,
        "internal_capping_input_path": "",
        "internal_capping_output_path": "",
        "protonation_input_source": str(protonation_result["protonation_input_source"]),
        "protonation_input_path": str(protonation_result["protonation_input_path"]),
        "protonation_output_path": str(protonation_output_path),
        "final_protein_input_paths": str(protonation_output_path),
    }


def _run_gaps_fragment_route(
    pdb_id: str,
    pdb_dir: Path,
    variant_label: str,
) -> dict[str, str]:
    """Run the retained gaps variant as temporary physical fragments.

    The conservative gaps variant should preserve the observed PDB fragments
    rather than reconstructing missing residues. For this route, FRUTON splits
    the representative monomer protein into connected fragments, lets GROMACS
    treat every fragment as its own molecule with proper termini, and then asks
    prepared_structure to merge the protonated fragment coordinate outputs. The
    raw split fragments are temporary implementation artefacts; the protonated
    fragment outputs remain visible until prepared assembly so failures stay
    inspectable.
    """
    if variant_label != "gaps":
        raise ValueError(
            f"_run_gaps_fragment_route only supports variant_label='gaps', got {variant_label!r}"
        )

    monomer_protein_input_path = _select_default_protein_input_path(
        pdb_id=pdb_id,
        pdb_dir=pdb_dir,
    )
    components_dir = pdb_dir / "components"
    temporary_fragment_dir = components_dir / ".tmp_gaps_fragments"

    if temporary_fragment_dir.exists():
        import shutil

        shutil.rmtree(temporary_fragment_dir)
    temporary_fragment_dir.mkdir(parents=True, exist_ok=True)

    _append_fruton_log(
        "protonation_route_gaps:start",
        [
            f"pdb_id                     : {pdb_id}",
            f"variant_label              : {variant_label}",
            f"monomer_protein_input_path : {monomer_protein_input_path}",
            f"temporary_fragment_dir     : {temporary_fragment_dir}",
            "internal_capping            : skipped",
            "internal_capping_backend    : none",
            f"protonation_method         : {FRUTON_PROTONATION_METHOD}",
            f"gromacs_ff                 : {DEFAULT_GROMACS_FORCE_FIELD}",
            f"gromacs_water              : {DEFAULT_GROMACS_WATER_MODEL}",
            "route                      : monomer protein -> temporary fragment split -> per-fragment sanitize -> per-fragment gmx pdb2gmx -> prepared_structure merge",
        ],
    )

    try:
        split_result = split_gap_variant_structure(
            pdb_id=pdb_id,
            protein_dir=pdb_dir,
            variant_label=variant_label,
            input_pdb_path=monomer_protein_input_path,
            output_dir=temporary_fragment_dir,
            verbose=False,
        )
    except Exception as error:
        _log_fruton_exception(
            f"protonation_route_gaps:fragment_split_failed:{pdb_id}",
            error,
        )
        return {
            "status": STATUS_FAILED,
            "variant_label": variant_label,
            "structure_input_path": str(monomer_protein_input_path),
            "structure_input_source": "protein_monomer",
            "internal_capping_status": STATUS_SKIPPED,
            "internal_capping_input_path": "",
            "internal_capping_output_path": "",
            "protonation_input_source": "",
            "protonation_input_path": "",
            "protonation_output_path": "",
            "final_protein_input_paths": "",
        }

    fragment_path_list = _extract_fragment_paths_from_split_result(split_result)
    if not fragment_path_list:
        _append_fruton_log(
            "protonation_route_gaps:fragment_split_empty",
            [
                f"pdb_id                 : {pdb_id}",
                f"variant_label          : {variant_label}",
                f"temporary_fragment_dir : {temporary_fragment_dir}",
                f"split_result           : {split_result}",
            ],
        )
        return {
            "status": STATUS_FAILED,
            "variant_label": variant_label,
            "structure_input_path": str(monomer_protein_input_path),
            "structure_input_source": "protein_monomer",
            "internal_capping_status": STATUS_SKIPPED,
            "internal_capping_input_path": "",
            "internal_capping_output_path": "",
            "protonation_input_source": "",
            "protonation_input_path": "",
            "protonation_output_path": "",
            "final_protein_input_paths": "",
        }

    protonated_fragment_path_list: list[Path] = []
    protonation_input_path_list: list[Path] = []
    protonation_output_path_list: list[Path] = []

    for fragment_index, fragment_path in enumerate(fragment_path_list, start=1):
        fragment_variant_label = f"{variant_label}_fragment_{fragment_index:02d}"

        try:
            sanitize_result = _run_sanitize_for_protonation_input(
                pdb_id=pdb_id,
                pdb_dir=pdb_dir,
                variant_label=fragment_variant_label,
                input_pdb=fragment_path,
            )
        except Exception as error:
            protonation_input_path_list.append(fragment_path)
            _log_fruton_exception(
                f"protonation_route_gaps:fragment_sanitize_exception:{pdb_id}:{fragment_variant_label}",
                error,
            )
            continue

        sanitized_fragment_path = Path(str(getattr(sanitize_result, "output_pdb_path", "")))
        protonation_input_path_list.append(sanitized_fragment_path)

        if not bool(getattr(sanitize_result, "can_run_gromacs", False)):
            _append_fruton_log(
                "protonation_route_gaps:fragment_sanitize_failed",
                [
                    f"pdb_id                  : {pdb_id}",
                    f"variant_label           : {variant_label}",
                    f"fragment_variant_label  : {fragment_variant_label}",
                    f"raw_fragment_path       : {fragment_path}",
                    f"sanitized_fragment_path : {sanitized_fragment_path}",
                    f"sanitize_log_path       : {getattr(sanitize_result, 'log_path', '')}",
                    f"status                  : {getattr(sanitize_result, 'status', '')}",
                    f"issues                  : {_summarize_sanitize_issues(sanitize_result)}",
                ],
            )
            continue

        try:
            protonation_result = _run_fruton_protonation_for_variant(
                pdb_id=pdb_id,
                pdb_dir=pdb_dir,
                variant_label=fragment_variant_label,
                input_pdb=sanitized_fragment_path,
                input_source="protein",
            )
        except Exception as error:
            _log_fruton_exception(
                f"protonation_route_gaps:fragment_protonation_exception:{pdb_id}:{fragment_variant_label}",
                error,
            )
            continue

        output_path_text = str(
            protonation_result.get("protonation_output_path", "")
        ).strip()
        output_path = Path(output_path_text) if output_path_text else None
        if (
            protonation_result.get("protonation_success", False)
            and output_path is not None
            and output_path.exists()
            and output_path.stat().st_size > 0
        ):
            protonated_fragment_path_list.append(output_path)
            protonation_output_path_list.append(output_path)
            continue

        _append_fruton_log(
            "protonation_route_gaps:fragment_protonation_failed",
            [
                f"pdb_id                    : {pdb_id}",
                f"variant_label             : {variant_label}",
                f"fragment_variant_label    : {fragment_variant_label}",
                f"raw_fragment_path         : {fragment_path}",
                f"sanitized_fragment_path   : {sanitized_fragment_path}",
                f"protonation_output_path   : {output_path_text}",
                f"protonation_error_message : {protonation_result.get('protonation_error_message', '')}",
            ],
        )

    if len(protonated_fragment_path_list) != len(fragment_path_list):
        return {
            "status": STATUS_FAILED,
            "variant_label": variant_label,
            "structure_input_path": str(monomer_protein_input_path),
            "structure_input_source": "protein_monomer",
            "internal_capping_status": STATUS_SKIPPED,
            "internal_capping_input_path": "",
            "internal_capping_output_path": "",
            "protonation_input_source": "temporary_fragments",
            "protonation_input_path": _path_list_to_text(protonation_input_path_list),
            "protonation_output_path": _path_list_to_text(protonation_output_path_list),
            "final_protein_input_paths": "",
        }

    _append_fruton_log(
        "protonation_route_gaps:success",
        [
            f"pdb_id                         : {pdb_id}",
            f"variant_label                  : {variant_label}",
            f"n_raw_fragments                : {len(fragment_path_list)}",
            f"n_protonated_fragments         : {len(protonated_fragment_path_list)}",
            f"raw_fragment_paths             : {_path_list_to_text(fragment_path_list)}",
            f"sanitized_fragment_paths       : {_path_list_to_text(protonation_input_path_list)}",
            f"protonated_fragment_paths      : {_path_list_to_text(protonated_fragment_path_list)}",
            "prepared_structure_input_policy : merge protonated fragment PDB files with TER separation",
        ],
    )

    return {
        "status": STATUS_SUCCESS,
        "variant_label": variant_label,
        "structure_input_path": str(monomer_protein_input_path),
        "structure_input_source": "protein_monomer",
        "internal_capping_status": STATUS_SKIPPED,
        "internal_capping_input_path": "",
        "internal_capping_output_path": "",
        "protonation_input_source": "temporary_fragments",
        "protonation_input_path": _path_list_to_text(protonation_input_path_list),
        "protonation_output_path": _path_list_to_text(protonated_fragment_path_list),
        "final_protein_input_paths": _path_list_to_text(protonated_fragment_path_list),
    }

def _run_protonation_route_for_variant(
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

    return _run_linear_protonation_route_for_variant(
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


def _choose_representative_successful_variant(
    successful_variant_result_list: list[dict[str, str]],
) -> dict[str, str]:
    """Choose which successful variant becomes the representative output.

    Gap-containing proteins can produce both a fragmented gaps variant and a
    complete filled variant. The fragmented gaps variant is useful to retain, but
    it should not become the main prepared-structure output when a complete model
    exists. This selector keeps the single variant as representative for
    proteins without gaps, prefers complete filled models for proteins with gaps,
    and only falls back to the gaps variant when no complete model was prepared.
    """

    priority = [
        "single",
        "best_complete",
        "large_gap_complete",
        "gaps",
    ]

    by_label = {
        result["variant_label"]: result
        for result in successful_variant_result_list
    }

    for variant_label in priority:
        if variant_label in by_label:
            return by_label[variant_label]

    raise RuntimeError("No successful prepared variant available.")


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

    _screen_step(14, "save_pipeline_json")
    _run_mutating_call(log_title="step_14:save_pipeline_json:captured_output", work_label="save_pipeline_json", screen_pdb_id="-", func=save_pipeline_table, protein_record_list=pipeline_record_list, json_path=pipeline_json_path)
    _screen_step(15, "write_pipeline_xlsx")
    _run_mutating_call(log_title="step_15:write_pipeline_xlsx:captured_output", work_label="write_pipeline_xlsx", screen_pdb_id="-", func=write_pipeline_to_xlsx, protein_record_list=pipeline_record_list, output_path=pipeline_xlsx_path)

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