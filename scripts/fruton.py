#!/usr/bin/env python3

from __future__ import annotations

import argparse
import contextlib
import io
import json
import os
import shutil
import sys
import random
import re
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

FRUTON_DEPENDENCY_REPORT = ensure_fruton_dependencies()

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
    FASTA_DIRECTORY_COLUMN_NAME,
    FASTA_FILES_DONE_COLUMN_NAME,
    FILLER_DIRECTORY_COLUMN_NAME,
    FILLER_MODEL_PATH_COLUMN_NAME,
    FILLER_MODEL_SOURCE_COLUMN_NAME,
    FILLER_STATUS_COLUMN_NAME,
    GAP_BOUNDARIES_COLUMN_NAME,
    GAP_SIZES_COLUMN_NAME,
    HAS_GAPS_COLUMN_NAME,
    HAS_COFACTORS_COLUMN_NAME,
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
    PREPARED_DIRECTORY_COLUMN_NAME,
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
    PROTONATION_PROPKA_HIS_COLUMN_NAME,
    PROTONATION_STATUS_COLUMN_NAME,
    RANGE_COLUMN_NAME,
    SEQUENCE_ALIGNMENT_DONE_COLUMN_NAME,
    METALL_PARAMS_DIRECTORY_COLUMN_NAME,
    METALL_PARAMS_MANIFEST_PATH_COLUMN_NAME,
    METALL_PARAMS_SITE_COUNT_COLUMN_NAME,
    METALL_PARAMS_STATUS_COLUMN_NAME,
    NONSTD_RESIDUE_PARAMS_MANIFEST_PATH_COLUMN_NAME,
    NONSTD_RESIDUE_PARAMS_N_RESIDUES_COLUMN_NAME,
    NONSTD_RESIDUE_PARAMS_STATUS_COLUMN_NAME,
    GAUSSIAN_PARAMS_STATUS_COLUMN_NAME,
    GAUSSIAN_PARAMS_MANIFEST_PATH_COLUMN_NAME,
    GAUSSIAN_PARAMS_MCPB_CHARGE_COLUMN_NAME,
    GAUSSIAN_PARAMS_MCPB_SPIN_COLUMN_NAME,
    GAUSSIAN_PARAMS_RESP_CHARGE_COLUMN_NAME,
    GAUSSIAN_PARAMS_FRCMOD_PATHS_COLUMN_NAME,
    MODEL_EVALUATION_STATUS_COLUMN_NAME,
    MODEL_EVALUATION_PCT_FAVORED_COLUMN_NAME,
    MODEL_EVALUATION_PCT_OUTLIER_COLUMN_NAME,
    MODEL_EVALUATION_CLASHSCORE_COLUMN_NAME,
    MODEL_EVALUATION_OVERALL_QUALITY_COLUMN_NAME,
    MODEL_EVALUATION_RAMA_PLOT_PATH_COLUMN_NAME,
    MODEL_EVALUATION_MANIFEST_PATH_COLUMN_NAME,
    REPORT_STATUS_COLUMN_NAME,
    REPORT_PDF_PATH_COLUMN_NAME,
    REPORT_FIGURE_PNG_PATH_COLUMN_NAME,
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
from stack_protein_preparation.pipeline_table import get_record_by_pdb_id
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
    _build_prepared_wat_structure_for_variant,
    _build_prepared_cofactor_structure_for_variant,
    _run_metall_params_for_protein,
    _run_nonstd_residue_params_for_protein,
    _prepare_crystal_water_component,
    _check_water_clashes_for_prepared_variant,
)
from stack_protein_preparation.model_evaluation import run_model_evaluation
from stack_protein_preparation.protein_report import generate_protein_report
from stack_protein_preparation.gaussian_hpc import (
    SlurmConfig,
    run_gaussian_parametrization_for_all_proteins,
)
from stack_protein_preparation.parameter_audit import (
    probe_amber_resp_parameters,
    RespParameterProbeResult,
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



_SCREEN_WIDTH = 116
_SCREEN_STEP_GAP_LINES = 2
_FRUTON_LOGO_FRU_WIDTH = 29

FRUTON_STRAWBERRY_LOGO = [
    "  ▗▄▖▗▄▖  ",
    "▗▄█████▄▖ ",
    "███ █ ███ ",
    "█████████ ",
    "▜██ █ ██▛ ",
    " ▜█████▛  ",
]

_C_GREEN = "\033[32m"
_C_RED = "\033[31m"
_C_YELLOW = "\033[33m"
_C_BLACK = "\033[30m"
_C_BLUE = "\033[34m"
_C_CYAN = "\033[36m"
_C_RESET = "\033[0m"
_C_DIM = "\033[2m"
_C_BOLD = "\033[1m"
_ANSI_RE = re.compile(r"\x1b\[[0-9;]*m")


def _terminal_uses_color() -> bool:
    if "NO_COLOR" in os.environ:
        return False
    if os.environ.get("TERM", "").strip().lower() == "dumb":
        return False
    return sys.stdout.isatty()


def _ansi(text: str, code: str) -> str:
    if not _terminal_uses_color():
        return text
    return f"{code}{text}{_C_RESET}"


def _green(s: str) -> str:
    return _ansi(s, _C_GREEN)


def _red(s: str) -> str:
    return _ansi(s, _C_RED)


def _yellow(s: str) -> str:
    return _ansi(s, _C_YELLOW)


def _black(s: str) -> str:
    return _ansi(s, _C_BLACK)


def _blue(s: str) -> str:
    return _ansi(s, _C_BLUE)


def _cyan(s: str) -> str:
    return _ansi(s, _C_CYAN)


def _dim(s: str) -> str:
    return _ansi(s, _C_DIM)


def _bold(s: str) -> str:
    return _ansi(s, _C_BOLD)


def _visible_len(text: str) -> int:
    return len(_ANSI_RE.sub("", text))


def _pad_visible(text: str, width: int) -> str:
    return text + " " * max(width - _visible_len(text), 0)


def _truncate_visible(text: str, width: int) -> str:
    plain = _ANSI_RE.sub("", text)
    if len(plain) <= width:
        return text
    return plain[: max(width - 1, 0)] + "…"


def _hr(char: str = "─", width: int = _SCREEN_WIDTH) -> str:
    return char * width


def _terminal_background_is_black_or_dark() -> bool:
    """Best-effort terminal background detection for the FRUTON logo.

    Background color detection is not portable across terminal emulators and
    HPC login shells. Some terminals expose ``COLORFGBG`` where the last field
    is the background color index; black is usually 0 and bright black/dark
    gray is usually 8. If detection is unavailable, FRUTON assumes a dark
    terminal because yellow remains readable while black can disappear. Users
    can override the decision with ``FRUTON_TERMINAL_BACKGROUND=light`` or
    ``FRUTON_TERMINAL_BACKGROUND=dark``.
    """

    override = os.environ.get("FRUTON_TERMINAL_BACKGROUND", "").strip().lower()
    if override in {"black", "dark"}:
        return True
    if override in {"white", "light"}:
        return False

    colorfgbg = os.environ.get("COLORFGBG", "").strip()
    if colorfgbg:
        background = colorfgbg.split(";")[-1]
        if background.isdigit():
            return int(background) in {0, 8}

    return True


def _fruton_ton_color_code() -> str:
    """Return the ANSI color for the TON half of the FRUTON logo."""

    override = os.environ.get("FRUTON_LOGO_TON_COLOR", "").strip().lower()
    if override == "black":
        return _C_BLACK
    if override == "yellow":
        return _C_YELLOW

    if _terminal_background_is_black_or_dark():
        return _C_YELLOW
    return _C_BLACK


def _format_fruton_logo_line(line: str, *, width: int | None = None) -> str:
    """Color one logo line as FRU/TON while preserving visible alignment."""

    rendered = line.center(width) if width is not None else line
    if not _terminal_uses_color():
        return rendered

    left_padding = 0
    if width is not None:
        left_padding = max((width - len(line)) // 2, 0)

    split_index = min(left_padding + _FRUTON_LOGO_FRU_WIDTH, len(rendered))
    fru = rendered[:split_index]
    ton = rendered[split_index:]
    return f"{_C_RED}{fru}{_C_RESET}{_fruton_ton_color_code()}{ton}{_C_RESET}"


def _format_fruton_logo_with_strawberry() -> str:
    """Return the colored FRUTON startup logo with a small strawberry mark.

    The strawberry is intentionally built from box-drawing characters rather
    than emoji so it remains monospaced and compatible with SSH terminals. The
    FRU/T0N split is kept as a column split rather than a string split because
    the Figlet-style letters have different visible widths. Interior spaces in
    the strawberry body are rendered as yellow blocks, which reads as seeds on
    a dark terminal. The plain fallback keeps the same geometry when colors are
    disabled.
    """

    logo_lines = FRUTON_LOGO.splitlines()
    logo_cols = max(len(line) for line in logo_lines)
    gap = " " * 5
    lines: list[str] = []

    for row, logo_line in enumerate(logo_lines):
        colored_logo = _format_fruton_logo_line(logo_line)
        berry_line = FRUTON_STRAWBERRY_LOGO[row]
        berry_chars: list[str] = []

        for col, char in enumerate(berry_line):
            if row < 2:
                berry_chars.append(_green(char) if char != " " else " ")
                continue

            if char == " ":
                left_has_body = any(c != " " for c in berry_line[:col])
                right_has_body = any(c != " " for c in berry_line[col + 1 :])
                berry_chars.append(_yellow("█") if left_has_body and right_has_body else " ")
                continue

            berry_chars.append(_red(char))

        lines.append(colored_logo + " " * (logo_cols - len(logo_line)) + gap + "".join(berry_chars))

    return "\n".join(lines)


def _print_box(title: str, rows: list[str], *, width: int = _SCREEN_WIDTH) -> None:
    """Print a rounded box with text rows padded by visible character width."""

    title_text = f" {title} "
    print("╭─" + title_text + "─" * max(width - len(title_text) - 3, 1) + "╮")
    inner_width = width - 4
    for row in rows:
        row = _truncate_visible(row, inner_width)
        print("│ " + _pad_visible(row, inner_width) + " │")
    print("╰" + "─" * (width - 2) + "╯")


def _print_abstract_separator(*, width: int = _SCREEN_WIDTH) -> None:
    """Print a subtle abstract dotted separator for the lower pipeline area."""

    motifs = ("○", "·", "·", "○", "·", "∘", "·", "○", "·", "·", "∘", "·")
    parts: list[str] = []
    for index in range(width):
        if index % 8 == 0:
            parts.append(_red("○") if index % 24 == 0 else _yellow("○"))
        elif index % 4 == 0:
            parts.append(_dim("·"))
        else:
            parts.append(" ")
    print("".join(parts).rstrip())


def _print_logo() -> None:
    print()
    print(_format_fruton_logo_with_strawberry())
    print(_hr("─", _SCREEN_WIDTH))
    _print_box(
        "FRUTON",
        [
            FRUTON_EXPANSION,
            "Protein preparation with explicit reconstruction, protonation, and topology checks.",
            "Screen output is compact; full subprocess stdout/stderr is written to fruton.log.",
        ],
    )
    print()


def _print_dependency_report() -> None:
    checks = FRUTON_DEPENDENCY_REPORT.checks
    n_ok = sum(1 for check in checks if check.available)
    n_fail = max(len(checks) - n_ok, 0)
    status = _green("OK") if n_fail == 0 else _red(f"{n_fail} missing")

    width = _SCREEN_WIDTH
    print("╭─ dependency check " + "─" * (width - 22) + "╮")
    print("│ " + "Available dependencies".ljust(39) + " │ " + _pad_visible(f"{n_ok}/{len(checks)}", width - 47) + " │")
    print("│ " + "Dependency status".ljust(39) + " │ " + _pad_visible(status, width - 47) + " │")
    print("├" + "─" * (width - 2) + "┤")

    col_width = 35
    n_cols = 3
    for index in range(0, len(checks), n_cols):
        row = checks[index : index + n_cols]
        parts: list[str] = []
        for check in row:
            mark = _green("✓") if check.available else _red("✗")
            label = f"{mark}  {check.name}"
            parts.append(_pad_visible(label, col_width))
        content = "".join(parts)
        print("│ " + _pad_visible(content, width - 4) + " │")

    print("╰" + "─" * (width - 2) + "╯")

    missing_with_hints = [check for check in checks if not check.available and check.hint.strip()]
    if missing_with_hints:
        print()
        _print_box(
            "missing dependency hints",
            [
                f"{check.name}: {line}"
                for check in missing_with_hints
                for line in check.hint.strip().splitlines()
            ],
        )
    print()


def _print_run_header(
    *,
    protein_data_dir: Path,
    n_records: int,
    target_pdb_id: str | None,
    from_step: int | None,
    force_field: str,
    water_model: str,
    ph: float,
) -> None:
    if "_PROGRESS" in globals():
        _PROGRESS.close_current_step()

    filter_text = f"filter {target_pdb_id}" if target_pdb_id else "off"
    rerun_text = f"from step {from_step}" if from_step is not None else "off"
    _print_box(
        "run context",
        [
            f"Data directory        {protein_data_dir}",
            f"Active proteins       {n_records}",
            f"PDB filter            {filter_text}",
            f"Force field           {force_field}",
            f"Water model           {water_model}",
            f"pH compatibility      {ph}",
            f"Rerun mode            {rerun_text}",
            f"FRUTON log            {FRUTON_LOG_PATH}",
        ],
    )
    _print_abstract_separator()
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
        lines.append(f"║{_format_fruton_logo_line(logo_line, width=content_width)}║")

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
    print(_hr("═"))

    rows = [
        ("PDB records processed", str(summary["total_records"])),
        ("Proteins with max gap > 5", str(summary["max_gap_gt_5"])),
        ("Proteins with max gap ≥ 8", str(summary["max_gap_ge_8"])),
        ("Single variant prepared", str(summary["single_best_available_written"])),
        ("Gaps variant prepared", str(summary["gaps_variant_written"])),
        ("MODELLER-complete variant prepared", str(summary["modeller_complete_written"])),
        ("AlphaFold-complete variant prepared", str(summary["alphafold_complete_written"])),
        ("Proteins with ≥ 1 prepared output", str(summary["proteins_with_prepared_output"])),
        ("Proteins with no prepared output", str(summary["proteins_failed"])),
    ]

    print(_build_summary_box(title="FRUTON  ·  run complete", rows=rows))
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
    """Render boxed, action-based FRUTON progress without a UI dependency.

    The progress display treats one call to ``_run_mutating_call`` as one real
    action. This matches the runner better than historical step numbers because
    FRUTON branches per protein, per retained variant, and per optional output.
    Step headers are printed as rounded panels with a red hollow dot separator,
    while the live progress line is updated only during real mutating actions.
    The class deliberately stays dependency-free so the runner remains usable on
    HPC login nodes, SSH sessions, and copied terminal logs.
    """

    def __init__(self, *, total_steps: int) -> None:
        self.total_steps = total_steps
        self.current_step = 0
        self.current_step_name = "startup"
        self.current_label = "initialising"
        self.current_pdb_id = "-"
        self.current_status = "idle"
        self.planned_actions = 0
        self.completed_actions = 0
        self.total_work = 0
        self.finalized = False
        self.started_at = time.monotonic()
        self.enabled = sys.stdout.isatty() and not _verbose_screen_enabled()
        self._has_started_step = False
        self._step_box_open = False
        self._last_rendered_progress = False

    def configure(self, *, total_work: int) -> None:
        """Set the estimated total action count so the bar scales correctly.

        The estimate is a best-effort count based on record count and typical
        pipeline branching. The bar uses max(estimate, planned_actions) as the
        denominator so it never exceeds 100% even if the estimate is too low,
        and scales across the full run instead of jumping to near-100% after
        a handful of actions.
        """

        self.total_work = total_work

    def start_step(self, *, step_number: int, step_name: str) -> None:
        self.close_current_step()

        self.current_step = step_number
        self.current_step_name = step_name
        self.current_label = step_name
        self.current_pdb_id = "-"
        self.current_status = "section"

        if self._has_started_step:
            print("\n" * _SCREEN_STEP_GAP_LINES, end="")
        else:
            print()
            self._has_started_step = True

        width = _SCREEN_WIDTH
        prefix = f" step {step_number:02d}/{self.total_steps:02d}  "
        suffix = f"  {step_name} "
        title_visible_len = len(prefix) + 1 + len(suffix)
        top = "╭─" + prefix + _red("○") + suffix + "─" * max(width - title_visible_len - 3, 1) + "╮"
        print(top)

        variants = _STEP_COMMENTARY.get(step_number)
        if variants:
            self._print_box_line(random.choice(variants), color="muted")

        self._step_box_open = True
        self._last_rendered_progress = False

    def close_current_step(self) -> None:
        if not self._step_box_open:
            return
        self.clear()
        print("╰" + "─" * (_SCREEN_WIDTH - 2) + "╯")
        self._step_box_open = False
        self._last_rendered_progress = False

    def register_action(self, *, label: str, pdb_id: str = "-") -> None:
        self.planned_actions += 1
        self.current_label = label
        self.current_pdb_id = pdb_id or "-"
        self.current_status = "running"
        self.render()

    def advance(self, *, label: str, pdb_id: str = "-", status: str = "done") -> None:
        self.completed_actions = min(self.completed_actions + 1, self.planned_actions)
        self.current_label = label
        self.current_pdb_id = pdb_id or "-"
        self.current_status = status
        self.render()

    def set_status(self, *, label: str, pdb_id: str = "-", status: str = "running") -> None:
        self.current_label = label
        self.current_pdb_id = pdb_id or "-"
        self.current_status = status
        self.render()

    def finish(self) -> None:
        self.finalized = True
        self.render()
        if self.enabled:
            print()
        self.close_current_step()

    def clear(self) -> None:
        if not self.enabled:
            return
        width = shutil.get_terminal_size((_SCREEN_WIDTH, 20)).columns
        print("\r" + " " * (width - 1) + "\r", end="", flush=True)

    def _print_box_line(self, text: str, *, color: str = "normal") -> None:
        content_width = _SCREEN_WIDTH - 4
        words = text.split(" ")
        lines: list[str] = []
        current = ""
        for word in words:
            candidate = (current + " " + word).lstrip() if current else word
            if len(candidate) <= content_width:
                current = candidate
            else:
                if current:
                    lines.append(current)
                current = word if len(word) <= content_width else word[:content_width]
        if current:
            lines.append(current)
        if not lines:
            lines = [""]
        for line in lines:
            rendered = _dim(line) if color == "muted" else line
            print("│ " + _pad_visible(rendered, content_width) + " │")

    def print_notice(self, message: str) -> None:
        if self._step_box_open:
            self._print_box_line(f"note  {message}")
        else:
            print(f"│  {_blue('note')}  {message}")

    def print_error(self, message: str) -> None:
        if self._step_box_open:
            self._print_box_line(f"{_red('FAIL')}  {message}")
        else:
            print(f"│  {_red('FAIL')}  {message}")

    def render(self) -> None:
        if not self.enabled:
            return

        width = shutil.get_terminal_size((_SCREEN_WIDTH, 20)).columns
        if self.finalized:
            denominator = max(self.planned_actions, self.total_work, 1)
        else:
            denominator = max(self.planned_actions + 1, self.total_work, 1)

        fraction = min(max(self.completed_actions / denominator, 0.0), 1.0)
        bar_width = 24
        filled = int(round(bar_width * fraction))
        bar = _green("█" * filled) + _dim("░" * (bar_width - filled))

        elapsed = int(time.monotonic() - self.started_at)
        elapsed_str = f"{elapsed // 60}m{elapsed % 60:02d}s" if elapsed >= 60 else f"{elapsed}s"
        status = _red("failed") if self.current_status == "failed" else self.current_status

        content = (
            f"[{bar}]  "
            f"action {self.completed_actions:03d}  "
            f"{self.current_pdb_id:<6}  "
            f"{self.current_label} · {status}  "
            f"{_dim(elapsed_str)}"
        )
        line = "│  " + content
        line = _truncate_visible(line, _SCREEN_WIDTH - 2)
        line = _pad_visible(line, _SCREEN_WIDTH - 1) + "│"

        if _visible_len(line) >= width:
            line = _truncate_visible(line, width - 1)
        print("\r" + line, end="", flush=True)
        self._last_rendered_progress = True


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


_PROGRESS = TerminalProgress(total_steps=19)

_STEP_COMMENTARY: dict[int, list[str]] = {
    1: [
        "Downloading PDB files. Crystallographers spent months at the synchrotron so you don't have to.",
        "Fetching coordinates from the PDB. Each file is the average position of an atom that was never actually still.",
        "Pulling structures from the Protein Data Bank. Someone solved these at 100 K so your simulation can fall apart at 300 K.",
        "Downloading PDB entries. The experimentalists froze it, cryo'd it, and shot X-rays at it. You get a text file.",
    ],
    4: [
        "SEQRES = what the biologist cloned. ATOM = what the X-rays actually saw. The difference is your problem now.",
        "Parsing SEQRES and ATOM records. One is the protein as intended; the other is the protein as expressed, crystallised, and partially observed.",
        "Cross-referencing sequence and structure. The gaps between SEQRES and ATOM are where biology refused to cooperate with crystallography.",
        "Reading PDB headers. SEQRES is the wish list. ATOM is what made it into the electron density. We need both.",
    ],
    5: [
        "Aligning sequences so FRUTON knows which residues the crystallographer quietly gave up on resolving.",
        "Running sequence alignment. Every dash in the output is a residue that was there in solution and nowhere to be found in the crystal.",
        "Mapping SEQRES onto ATOM coordinates. The alignment is how FRUTON decides which missing residues are a modeling problem versus a cloning artifact.",
        "Aligning full sequence to observed structure. Gaps are not errors — they are flexible regions politely declining to be photographed.",
    ],
    6: [
        "Renumbering insertion codes. Yes, someone once numbered residues 122, 122A, 122B, 123. We don't talk about it.",
        "Cleaning up insertion codes left by PDB depositors who ran out of numbers and started using the alphabet.",
        "Standardising residue numbering. Insertion codes exist because sequential integers cannot capture the full chaos of protein evolution.",
        "Fixing insertion code numbering. The PDB allows 122, 122A, 122B as valid residue identifiers. GROMACS does not share this flexibility.",
    ],
    7: [
        "Picking one representative chain. The other copies in the asymmetric unit are just expensive crystallographic ballast.",
        "Selecting the biological unit chain. Crystal packing produces multiple copies; only one of them is the protein your MD cares about.",
        "Choosing the representative chain from the asymmetric unit. The others are required for diffraction; they are not required for you.",
        "Identifying the target chain. Crystallographers deposit all chains out of obligation. We keep one out of necessity.",
    ],
    8: [
        "Counting gaps. Each one is a loop that refused to sit still long enough to be photographed at −173 °C.",
        "Enumerating missing residues. A gap in the structure is a loop that was too dynamic to resolve — and too important to ignore.",
        "Cataloguing unresolved regions. These are not missing data. They are regions where the protein was moving faster than the X-rays.",
        "Surveying structural gaps. Every missing loop was present in the crystal. It just would not stop moving long enough to leave coordinates.",
    ],
    9: [
        "Asking MODELLER to invent plausible coordinates for the residues that ghosted the crystallographer.",
        "Running loop modelling. MODELLER will generate coordinates for regions that the electron density refused to provide.",
        "Reconstructing missing loops with MODELLER. The output is physically reasonable and energetically defensible. It is not experimental data.",
        "Homology modelling of disordered regions. The loops were flexible in the crystal; MODELLER will make them look confident anyway.",
    ],
    10: [
        "Protonation determines charge; charge determines force; force determines trajectory. This step commits the physics. May the force be with you.",
        "Running gmx pdb2gmx. Somewhere in Groningen, a developer is quietly proud that this still works.",
        "Assigning protonation states. At physiological pH, histidine is neither fully protonated nor fully deprotonated, and it is absolutely your problem.",
        "Adding hydrogens. The PDB does not store them; the force field cannot live without them. This step bridges the gap.",
    ],
    11: [
        "Cleaning up after protonation. Alternate conformations collapsed, occupancies normalised, dignity partially restored.",
        "Sanitising the protonated structure. Alternate conformers selected, partial occupancies resolved, redundant atoms removed.",
        "Post-protonation cleanup. The structure that emerges from pdb2gmx is correct but messy. This step makes it presentable.",
        "Normalising the protonated PDB. Alternate side-chain conformations are physically real but computationally inconvenient — we pick one.",
    ],
    12: [
        "Metal parametrization. Because telling GROMACS 'it's just a zinc, don't worry about it' has never worked.",
        "Handling metal ions. Standard force fields treat zinc as a point charge. That is wrong, but it is consistently wrong, which is the best we can do.",
        "Parametrizing metal centres. Zinc, calcium, magnesium — each one is a coordination chemistry problem disguised as a topology entry.",
        "Metal ion setup. If the active site has a metal, the force field has an opinion about it. This step reconciles the two.",
    ],
    13: [
        "Deriving RESP charges for non-standard residues. The force field has no idea what a phosphotyrosine is. We do.",
        "Calculating electrostatic potential for ligand parametrization. Gaussian will spend real CPU time computing charges that AMBER cannot guess.",
        "Fitting RESP charges. The force field was trained on standard amino acids. Non-standard residues must earn their parameters.",
        "Generating partial charges for novel residues. RESP fitting ensures the electrostatics are consistent with the QM wavefunction, not invented.",
    ],
    14: [
        "Submitting Gaussian jobs to HPC. HF/6-31G* — expensive enough to be right, cheap enough to finish before your thesis.",
        "Running quantum chemistry calculations. HF/6-31G* is the level of theory at which RESP charges best reproduce experimental solvation free energies.",
        "Queuing Gaussian on the cluster. The HPC scheduler will decide when your science happens. FRUTON will wait.",
        "Hartree-Fock single-point calculations in progress. The wavefunction does not care about your deadline.",
    ],
    15: [
        "Auditing force field coverage. If even one atom has no parameters, GROMACS will tell you in the least helpful way possible.",
        "Checking that every atom in the system has a force field entry. GROMACS error messages for missing parameters are famously unhelpful.",
        "Verifying parameter completeness. A topology with a missing atom type will run fine right up until it segfaults at step one.",
        "Parameter audit. The force field covers standard residues perfectly and everything else optimistically. This step finds the gaps.",
    ],
    16: [
        "Ramachandran check. If more than 2% of residues are outliers, someone had a bad day at the refinement stage.",
        "Validating backbone geometry. Ramachandran outliers are not a simulation problem — they are a structure quality problem inherited from the PDB.",
        "Checking phi/psi angles. A well-refined structure should have >98% of residues in favoured regions. This step verifies the inheritance.",
        "Ramachandran analysis. Residues in disallowed regions were placed there by the refinement program. They are now your starting coordinates.",
    ],
    17: [
        "Writing the provenance report. Future you — or your reviewer — will be glad this exists.",
        "Generating the pipeline report. Every decision FRUTON made is documented here, in case a referee asks why.",
        "Exporting provenance metadata. The report answers: what was the input, what was changed, and why the structure looks the way it does.",
        "Building the audit trail. Science is reproducible only if you can remember what you did. FRUTON remembers for you.",
    ],
    18: [
        "Saving pipeline state to JSON. The one file format that will still be readable in 30 years.",
        "Serialising run state to JSON. If FRUTON crashes after this step, you can restart without losing anything.",
        "Writing pipeline checkpoint. JSON is not glamorous, but it has survived every file format war of the last two decades.",
        "Persisting run metadata. The JSON file is the single source of truth for what this pipeline produced and how.",
    ],
    19: [
        "Writing to Excel. Science ends; bureaucracy is eternal.",
        "Exporting results to xlsx. The PI asked for a spreadsheet. The spreadsheet will outlive the simulation.",
        "Generating the Excel summary. Not every stakeholder reads JSON. Some of them have pivot tables.",
        "Writing xlsx. The data will be copy-pasted into a PowerPoint within 48 hours. FRUTON has made peace with this.",
    ],
}


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

    fixed_units = 4 + 6  # pdb_sync + save_json + write_xlsx + headroom
    per_record_units = 15  # fasta/align/insertion/rep_unit/components/monomer/metals/gap/filler
    # + protonation/prepared/sanitize/parameter_audit + parametrization + stage_parameters
    return fixed_units + record_count * per_record_units


def _screen_step(step_number: int, step_name: str) -> None:
    _PROGRESS.start_step(step_number=step_number, step_name=step_name)


def _screen_item(message: str) -> None:
    _append_fruton_log("screen", [message])
    if _verbose_screen_enabled():
        _PROGRESS.clear()
        _PROGRESS.print_notice(message)


def _screen_notice(message: str) -> None:
    _PROGRESS.clear()
    _PROGRESS.print_notice(message)


def _screen_error(message: str) -> None:
    _PROGRESS.clear()
    _PROGRESS.print_error(message)


def _find_nonstd_resp_param_files(
    manifest_path: str,
    unsupported_residue_names: tuple[str, ...],
) -> tuple[dict[str, Path], dict[str, Path]]:
    """Return (frcmod_paths, mol2_paths) dicts for each unsupported residue.

    Reads the nonstd_params_manifest.json produced by step 13 and looks up
    step03_resp_params/{resname}.frcmod and {resname}.mol2 inside each
    residue's resp_dir.  Only residues whose resname matches an entry in
    unsupported_residue_names are returned.  Files that do not exist yet
    (Gaussian still pending) are simply absent from the returned dicts.
    """
    import json as _json

    frcmod_by_resname: dict[str, Path] = {}
    mol2_by_resname: dict[str, Path] = {}

    if not manifest_path:
        return frcmod_by_resname, mol2_by_resname

    manifest_file = Path(manifest_path)
    if not manifest_file.is_file():
        return frcmod_by_resname, mol2_by_resname

    try:
        manifest = _json.loads(manifest_file.read_text(encoding="utf-8"))
    except Exception:
        return frcmod_by_resname, mol2_by_resname

    nonstd_dir = manifest_file.parent
    upper_names = {n.upper() for n in unsupported_residue_names}

    for residue_entry in manifest.get("residues", []):
        resname = str(residue_entry.get("resname", "")).strip().upper()
        if resname not in upper_names:
            continue
        label = residue_entry.get("label", "")
        resp_dir = nonstd_dir / label / "resp"
        resp_params = resp_dir / "step03_resp_params"
        frcmod = resp_params / f"{resname}.frcmod"
        mol2 = resp_params / f"{resname}.mol2"
        # Use the first match found per residue name (handles multiple sites)
        if resname not in frcmod_by_resname and frcmod.is_file():
            frcmod_by_resname[resname] = frcmod
        if resname not in mol2_by_resname and mol2.is_file():
            mol2_by_resname[resname] = mol2

    return frcmod_by_resname, mol2_by_resname


# Elements that require MCPB (QM-optimised) force field — standard non-bonded
# parameters are NOT acceptable for these transition metals.
_TRANSITION_METALS_NEEDING_MCPB: frozenset[str] = frozenset(
    {"ZN", "CU", "FE", "MN", "CO", "NI", "MO", "CD", "HG", "PD", "PT"}
)


def _find_amber_dat_dir() -> "Path | None":
    """Return the AMBER dat/leap directory, or None if not found."""
    import subprocess as _sp

    amberhome = os.environ.get("AMBERHOME", "")
    if amberhome:
        p = Path(amberhome) / "dat" / "leap"
        if p.is_dir():
            return p
    try:
        r = _sp.run(
            ["find", "/opt/cesga/2020", "-maxdepth", "15", "-name", "leaprc.phosaa14SB"],
            capture_output=True, text=True, timeout=20,
        )
        for line in r.stdout.strip().splitlines():
            p = Path(line.strip()).parent.parent  # .../dat/leap/cmd/ → .../dat/leap/
            if (p / "lib" / "phosaa14SB.lib").is_file():
                return p
    except Exception:
        pass
    return None


def _protein_has_nonstd_needing_gaussian(pipeline_record: dict) -> bool:
    """Return True if at least one nonstd residue needs Gaussian (not found in a database)."""
    manifest_path = str(
        pipeline_record.get(NONSTD_RESIDUE_PARAMS_MANIFEST_PATH_COLUMN_NAME, "") or ""
    )
    if not manifest_path or not Path(manifest_path).is_file():
        return True  # conservative: assume yes if manifest missing
    try:
        manifest = json.loads(Path(manifest_path).read_text(encoding="utf-8"))
        return any(
            not r.get("in_database", False) and r.get("resp_scaffold_generated", False)
            for r in manifest.get("residues", [])
        )
    except Exception:
        return True


def _element_from_ion(ion_str: str) -> str:
    m = re.match(r"([A-Za-z]+)", ion_str.strip())
    return m.group(1).upper() if m else ion_str.strip().upper()


def _find_phosaa14sb_xml() -> "Path | None":
    """Locate the OpenMM XML for phosaa14SB (from openmmforcefields package)."""
    try:
        import openmmforcefields as _off
        p = Path(_off.__file__).parent / "ffxml" / "amber" / "phosaa14SB.xml"
        if p.is_file():
            return p
    except ImportError:
        pass
    # Fallback: known path in the gbsa-grub pixi environment on this cluster
    fallback = Path(
        "/mnt/netapp2/Store_uni/home/otras/hcx/lwa/repos/gbsa-grub"
        "/.pixi/envs/default/lib/python3.12/site-packages"
        "/openmmforcefields/ffxml/amber/phosaa14SB.xml"
    )
    if fallback.is_file():
        return fallback
    return None


def _stage_amber_parameters_for_protein(
    *,
    pdb_id: str,
    pdb_dir: Path,
    pipeline_record: dict,
    amber_dat_dir: "Path | None",
) -> dict:
    """Stage all AMBER parameter files for tleap into each prepared variant directory.

    Creates {prepared_variant_dir}/amber_params/ with copies of frcmod/mol2/lib files
    needed to build the AMBER system, plus a ready-to-use tleap_build.in script.
    Transition metals (Zn, Cu, Fe …) that require MCPB but whose parametrization is
    incomplete are flagged in INCOMPLETE.txt rather than substituted with non-bonded params.
    """
    prepared_output_path = str(pipeline_record.get("prepared_structure.output_path", "") or "")
    if not prepared_output_path or not Path(prepared_output_path).is_file():
        return {"status": "skipped", "message": "no prepared structure"}

    prepared_variant_dir = Path(prepared_output_path).parent
    params_dir = prepared_variant_dir / "amber_params"
    params_dir.mkdir(parents=True, exist_ok=True)

    # ── split prepared PDB into protein / cofactor / WAT variants ─────────────
    components_dir = pdb_dir / "components"

    def _resnames_in_pdb(path: Path) -> set:
        names: set = set()
        if not path.is_file():
            return names
        for ln in path.read_text(encoding="utf-8", errors="replace").splitlines():
            if ln[:6] in ("ATOM  ", "HETATM"):
                names.add(ln[17:21].strip().upper())
        return names

    ligand_resnames = _resnames_in_pdb(components_dir / f"{pdb_id}_ligand.pdb")
    cofactor_resnames = _resnames_in_pdb(components_dir / f"{pdb_id}_cofactor.pdb")
    metal_resnames = _resnames_in_pdb(components_dir / f"{pdb_id}_metals.pdb")

    src_pdb = Path(prepared_output_path)
    header_lines: list[str] = []
    protein_lines: list[str] = []
    cofactor_lines: list[str] = []
    metal_hetatm_lines: list[str] = []

    for raw_line in src_pdb.read_text(encoding="utf-8", errors="replace").splitlines():
        line = raw_line.rstrip()
        if not line:
            continue
        rec = line[:6]
        if rec in ("REMARK", "HEADER", "TITLE ", "CRYST1",
                   "ORIGX1", "ORIGX2", "ORIGX3", "SCALE1", "SCALE2", "SCALE3"):
            header_lines.append(line)
        elif rec == "ATOM  ":
            protein_lines.append(line)
        elif rec.strip() == "TER":
            protein_lines.append(line)
        elif rec == "HETATM":
            resname = line[17:21].strip().upper()
            if resname in ligand_resnames or resname in cofactor_resnames:
                cofactor_lines.append(line)
            elif resname in metal_resnames:
                metal_hetatm_lines.append(line)
            else:
                # non-standard amino acid (PTR, SEP, TPO, CSD …) — stays in protein
                protein_lines.append(line)

    water_lines: list[str] = []
    water_comp = components_dir / f"{pdb_id}_water.pdb"
    if water_comp.is_file():
        for raw in water_comp.read_text(encoding="utf-8", errors="replace").splitlines():
            ln = raw.rstrip()
            if ln and ln[:6] in ("ATOM  ", "HETATM"):
                water_lines.append(ln)

    # PDBID.pdb — pure protein (overwrites existing)
    src_pdb.write_text(
        "\n".join(header_lines + protein_lines + ["END"]) + "\n",
        encoding="utf-8",
    )
    # cofactor.pdb — HETATM ligands and cofactors (only written if present)
    cofactor_pdb_path = prepared_variant_dir / "cofactor.pdb"
    if cofactor_lines:
        cofactor_pdb_path.write_text(
            "\n".join(cofactor_lines + ["END"]) + "\n",
            encoding="utf-8",
        )
    # PDBID_WAT.pdb — protein + metal ions + crystal waters (overwrites existing)
    wat_pdb_path = prepared_variant_dir / f"{src_pdb.stem}_WAT.pdb"
    wat_pdb_path.write_text(
        "\n".join(
            header_lines + protein_lines + metal_hetatm_lines + water_lines + ["END"]
        ) + "\n",
        encoding="utf-8",
    )

    extra_leaprc: list[str] = []
    loadoff_lines: list[str] = []
    frcmod_lines: list[str] = []
    mol2_loads: list[str] = []
    incomplete_notes: list[str] = []
    has_phosaa14sb = False

    # ── nonstd residue parameters ─────────────────────────────────────────────
    nonstd_manifest_path = str(
        pipeline_record.get(NONSTD_RESIDUE_PARAMS_MANIFEST_PATH_COLUMN_NAME, "") or ""
    )
    if nonstd_manifest_path and Path(nonstd_manifest_path).is_file():
        try:
            nonstd_manifest = json.loads(Path(nonstd_manifest_path).read_text(encoding="utf-8"))
        except Exception:
            nonstd_manifest = {}
        nonstd_dir = Path(nonstd_manifest_path).parent
        for res in nonstd_manifest.get("residues", []):
            resname = str(res.get("resname", "")).strip().upper()
            if not resname:
                continue
            if res.get("in_database"):
                db = str(res.get("database", "")).strip()
                if db.lower() == "phosaa14sb":
                    has_phosaa14sb = True
                if amber_dat_dir is not None:
                    lib_src = amber_dat_dir / "lib" / f"{db}.lib"
                    if lib_src.is_file():
                        lib_dst = params_dir / f"{db}.lib"
                        if not lib_dst.exists():
                            shutil.copy2(lib_src, lib_dst)
                        if f"loadOff amber_params/{db}.lib" not in loadoff_lines:
                            loadoff_lines.append(f"loadOff amber_params/{db}.lib")
                    else:
                        lr = f"source leaprc.{db}"
                        if lr not in extra_leaprc:
                            extra_leaprc.append(lr)
                else:
                    lr = f"source leaprc.{db}"
                    if lr not in extra_leaprc:
                        extra_leaprc.append(lr)
            else:
                # RESP frcmod+mol2 expected in step03_resp_params
                label = str(res.get("label", "")).strip()
                resp_params = nonstd_dir / label / "resp" / "step03_resp_params"
                frcmod_src = resp_params / f"{resname}.frcmod"
                mol2_src = resp_params / f"{resname}.mol2"
                if frcmod_src.is_file() and mol2_src.is_file():
                    shutil.copy2(frcmod_src, params_dir / f"{resname}.frcmod")
                    shutil.copy2(mol2_src, params_dir / f"{resname}.mol2")
                    frcmod_lines.append(f"loadAmberParams amber_params/{resname}.frcmod")
                    mol2_loads.append(f"{resname} = loadMol2 amber_params/{resname}.mol2")
                else:
                    gauss_status = str(pipeline_record.get(GAUSSIAN_PARAMS_STATUS_COLUMN_NAME, "") or "")
                    if gauss_status == STATUS_FAILED:
                        incomplete_notes.append(
                            f"{resname}: RESP Gaussian failed — re-run Gaussian and "
                            f"run_after_gaussian.sh, then re-stage."
                        )
                    else:
                        incomplete_notes.append(
                            f"{resname}: RESP parametrization pending."
                        )

    # Copy OpenMM XML for phosaa14SB so gbsa-grub can load it via extra_ff_files
    if has_phosaa14sb:
        xml_src = _find_phosaa14sb_xml()
        if xml_src is not None:
            dst_xml = params_dir / "phosaa14SB.xml"
            if not dst_xml.exists():
                shutil.copy2(xml_src, dst_xml)

    # ── metal parameters ──────────────────────────────────────────────────────
    has_metals = str(pipeline_record.get(HAS_METALS_COLUMN_NAME, "")).strip().lower() == "yes"
    if has_metals:
        ion_type_str = str(pipeline_record.get("metals.ion_type", "") or "").strip()
        ion_tokens = [t.strip() for t in ion_type_str.replace(";", ",").split(",") if t.strip()]
        ion_elements = {_element_from_ion(t) for t in ion_tokens}
        transition_elements = ion_elements & _TRANSITION_METALS_NEEDING_MCPB

        # Collect MCPB frcmod+mol2 from gaussian_params manifest
        mcpb_files_staged: list[str] = []
        gauss_manifest_path = str(
            pipeline_record.get(GAUSSIAN_PARAMS_MANIFEST_PATH_COLUMN_NAME, "") or ""
        )
        if gauss_manifest_path and Path(gauss_manifest_path).is_file():
            try:
                gauss_manifest = json.loads(Path(gauss_manifest_path).read_text(encoding="utf-8"))
            except Exception:
                gauss_manifest = {}
            mcpb_result = gauss_manifest.get("mcpb_result", {})
            if str(mcpb_result.get("status", "")) == STATUS_SUCCESS:
                for site in mcpb_result.get("site_results", []):
                    if str(site.get("status", "")) != STATUS_SUCCESS:
                        continue
                    group_name = str(site.get("group_name", "")).strip()
                    if not group_name:
                        continue
                    amber_params_dir = (
                        pdb_dir / "metall_params" / group_name / "03_mcpb" / "step03_amber_params"
                    )
                    for frcmod_src in amber_params_dir.glob("*_mcpbpy.frcmod"):
                        dst = params_dir / frcmod_src.name
                        if not dst.exists():
                            shutil.copy2(frcmod_src, dst)
                        frcmod_lines.append(f"loadAmberParams amber_params/{frcmod_src.name}")
                        mcpb_files_staged.append(frcmod_src.stem)
                    for mol2_src in amber_params_dir.glob("*.mol2"):
                        dst = params_dir / mol2_src.name
                        if not dst.exists():
                            shutil.copy2(mol2_src, dst)
                        stem = mol2_src.stem.upper()
                        mol2_loads.append(f"{stem} = loadMol2 amber_params/{mol2_src.name}")

        # Transition metals without completed MCPB → flag as incomplete
        if transition_elements and not mcpb_files_staged:
            gauss_status = str(pipeline_record.get(GAUSSIAN_PARAMS_STATUS_COLUMN_NAME, "") or "")
            metall_status = str(pipeline_record.get(METALL_PARAMS_MANIFEST_PATH_COLUMN_NAME, "") or "")
            metals_class = str(pipeline_record.get("metals.class", "") or "")
            if gauss_status == STATUS_FAILED:
                incomplete_notes.append(
                    f"{', '.join(sorted(transition_elements))}: MCPB Gaussian failed — "
                    f"re-run Gaussian and run_after_gaussian.sh for the metal site."
                )
            elif metals_class in ("transition", "mixed"):
                incomplete_notes.append(
                    f"{', '.join(sorted(transition_elements))}: MCPB not completed — "
                    f"geometry mismatch or Gaussian not yet run. Manual review required."
                )
            else:
                incomplete_notes.append(
                    f"{', '.join(sorted(transition_elements))}: MCPB parameters missing."
                )

    # ── write tleap_build.in ─────────────────────────────────────────────────
    struct_name = src_pdb.name  # e.g. 1A5H.pdb

    tleap_lines: list[str] = []
    if incomplete_notes:
        tleap_lines += [
            "# WARNING: incomplete parameters — see INCOMPLETE.txt before running.",
            "",
        ]
    tleap_lines += [
        "# Run from the prepared variant directory:",
        f"#   cd {prepared_variant_dir}",
        "#   module load amber  # (or equivalent)",
        "#   tleap -f amber_params/tleap_build.in",
        "#",
        "# Available structure files in this directory:",
        f"#   {struct_name} — pure protein (ATOM + non-standard residues)",
        f"#   {src_pdb.stem}_WAT.pdb — protein + crystal waters + metal ions",
        "#   cofactor.pdb — HETATM ligands and cofactors (user selects what to dock)",
        "",
        "source leaprc.protein.ff14SB",
        "source leaprc.water.tip3p",
    ]
    for lr in sorted(set(extra_leaprc)):
        tleap_lines.append(lr)
    for lo in sorted(set(loadoff_lines)):
        tleap_lines.append(lo)
    for fm in frcmod_lines:
        tleap_lines.append(fm)
    for ml in mol2_loads:
        tleap_lines.append(ml)
    tleap_lines += [
        "",
        f"mol = loadPdb {struct_name}",
        "check mol",
        "",
        "# Solvation — adjust box size / neutralization to your needs:",
        "# solvatebox mol TIP3PBOX 12.0",
        "# addions mol Na+ 0",
        "",
        f"saveAmberParm mol amber_params/{pdb_id}.prmtop amber_params/{pdb_id}.inpcrd",
        "quit",
    ]
    (params_dir / "tleap_build.in").write_text("\n".join(tleap_lines) + "\n", encoding="utf-8")

    if incomplete_notes:
        content = (
            f"# {pdb_id} — incomplete AMBER parameters\n\n"
            + "\n".join(f"  - {n}" for n in incomplete_notes)
            + "\n"
        )
        (params_dir / "INCOMPLETE.txt").write_text(content, encoding="utf-8")

    status = "incomplete" if incomplete_notes else "success"
    msg = f"{len(incomplete_notes)} parameter(s) missing" if incomplete_notes else ""
    return {"status": status, "params_dir": str(params_dir), "message": msg}


def _run_mutating_call(
    *,
    log_title: str,
    work_label: str,
    screen_pdb_id: str,
    func: Any,
    **kwargs: Any,
) -> Any:
    """Run one file-changing function while keeping terminal output compact.

    Each call is counted as one visible FRUTON action. This keeps the progress
    display tied to executed work instead of guessed per-record weights or
    historical pipeline step numbers. Captured stdout and stderr are written to
    the FRUTON log on success and failure. Exceptions are re-raised so the
    existing caller-level pipeline semantics remain unchanged.
    """

    _PROGRESS.register_action(label=work_label, pdb_id=screen_pdb_id)

    captured = _CapturedOutput()
    try:
        with captured:
            result = func(**kwargs)
    except Exception:
        captured_lines = captured.captured_lines()
        if captured_lines:
            _append_fruton_log(log_title, captured_lines)
        _PROGRESS.advance(label=work_label, pdb_id=screen_pdb_id, status="failed")
        raise

    captured_lines = captured.captured_lines()
    if captured_lines:
        _append_fruton_log(log_title, captured_lines)

    _PROGRESS.advance(label=work_label, pdb_id=screen_pdb_id, status="done")
    return result



def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="FRUTON — protein preparation pipeline",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "--data-dir",
        type=Path,
        default=None,
        metavar="DIR",
        help=(
            "Directory containing pdb_ids.csv and pipeline data. "
            "Defaults to the current working directory if it contains pdb_ids.csv, "
            "otherwise <project_root>/data/proteins."
        ),
    )
    parser.add_argument(
        "--pdb-id",
        default=None,
        metavar="PDB_ID",
        help="Process only this PDB ID (case-insensitive).",
    )
    parser.add_argument(
        "--from-step",
        type=int,
        default=None,
        metavar="N",
        help=(
            "Resume pipeline from step N (4–17). Steps before N are skipped and "
            "existing pipeline.json state is loaded for those steps. "
            "Requires --pdb-id when used to track per-protein rerun history."
        ),
    )
    slurm = parser.add_argument_group("Slurm (Gaussian parametrization)")
    slurm.add_argument(
        "--slurm-partition", default="short", metavar="NAME",
        help="Slurm partition for Gaussian jobs (default: short)",
    )
    slurm.add_argument(
        "--slurm-time", default="06:00:00", metavar="HH:MM:SS",
        help="Slurm time limit per Gaussian job (default: 06:00:00; short partition max is 6 h)",
    )
    slurm.add_argument(
        "--slurm-mem", default="2000M", metavar="MEM",
        help="Memory per CPU per Gaussian job, passed as --mem-per-cpu (default: 2000M; 24 CPUs × 2000M = 48 GB)",
    )
    slurm.add_argument(
        "--slurm-cpus", type=int, default=24, metavar="N",
        help="CPUs per Gaussian job (default: 24)",
    )
    slurm.add_argument(
        "--slurm-gpus", type=int, default=0, metavar="N",
        help="GPUs per Gaussian job (default: 0; CPU-only Gaussian run)",
    )
    return parser.parse_args()


def _resolve_protein_data_dir(data_dir_arg: Path | None) -> Path:
    if data_dir_arg is not None:
        return data_dir_arg.resolve()
    cwd = Path.cwd()
    if (cwd / "pdb_ids.csv").is_file():
        return cwd
    raise SystemExit(
        f"\nERROR: pdb_ids.csv not found in the current directory:\n"
        f"  {cwd}\n\n"
        f"Run fruton.py from the directory that contains pdb_ids.csv "
        f"(e.g. your $LUSTRE working directory),\n"
        f"or pass --data-dir <path> explicitly.\n\n"
        f"All protein directories and intermediate files — including Gaussian inputs — "
        f"will be created as subdirectories of that location."
    )


def _write_rerun_markers(
    pipeline_record_list: list[dict[str, str]],
    from_step: int,
) -> None:
    """Append a rerun marker to each protein's rerun_markers.json file."""
    now = datetime.utcnow().isoformat(timespec="seconds") + "Z"
    for record in pipeline_record_list:
        pdb_dir = Path(record.get(PDB_DIRECTORY_COLUMN_NAME, ""))
        if not pdb_dir.name:
            continue
        markers_path = pdb_dir / "rerun_markers.json"
        try:
            existing: list[dict] = json.loads(markers_path.read_text()) if markers_path.is_file() else []
        except Exception:
            existing = []
        existing.append({"step": from_step, "timestamp": now})
        markers_path.write_text(json.dumps(existing, indent=2))


def run_pipeline() -> None:
    global FRUTON_LOG_DIR, FRUTON_LOG_PATH

    args = _parse_args()
    protein_data_dir = _resolve_protein_data_dir(args.data_dir)
    pdb_ids_csv_path = protein_data_dir / "pdb_ids.csv"
    pipeline_json_path = protein_data_dir / "pipeline.json"
    pipeline_xlsx_path = protein_data_dir / "pipeline.xlsx"

    target_pdb_id: str | None = args.pdb_id.strip().upper() if args.pdb_id else None
    from_step: int | None = args.from_step

    def _step_runs(n: int) -> bool:
        return from_step is None or n >= from_step

    FRUTON_LOG_DIR = protein_data_dir / "logs" / "fruton"
    FRUTON_LOG_PATH = protein_data_dir / "logs" / "fruton" / "fruton.log"

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
    _screen_notice(f"pipeline initialised  ·  {_timestamp()}")
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
    _PROGRESS.configure(total_work=_estimate_initial_work_units(len(pdb_record_list)))
    _screen_item(f"loaded {len(pdb_record_list)} input records")
    _append_fruton_log(
        "step_2:read_input_csv",
        [f"record_count: {len(pdb_record_list)}"],
    )

    _screen_step(3, "build_pipeline_state")
    fresh_record_list = build_pipeline_records_from_input_csv(
        pdb_record_list=pdb_record_list,
        protein_data_dir=protein_data_dir,
    )
    if target_pdb_id is not None:
        fresh_record_list = [r for r in fresh_record_list if r[PDB_ID_COLUMN_NAME] == target_pdb_id]
        if not fresh_record_list:
            _screen_error(f"--pdb-id {target_pdb_id!r} not found in {pdb_ids_csv_path}")
            return

    if from_step is not None:
        # Rerun mode: load existing pipeline state and overlay on fresh records so that
        # fields from steps < from_step are preserved from the previous run.
        existing_records = load_pipeline_table(pipeline_json_path)
        pipeline_record_list = []
        for fresh_rec in fresh_record_list:
            pdb_id = fresh_rec[PDB_ID_COLUMN_NAME]
            existing_rec = get_record_by_pdb_id(existing_records, pdb_id)
            if existing_rec is not None:
                merged = dict(existing_rec)
                # Always use current paths (machine may differ)
                for col in (
                    PDB_DIRECTORY_COLUMN_NAME, FASTA_DIRECTORY_COLUMN_NAME,
                    ALIGNMENT_DIRECTORY_COLUMN_NAME, COMPONENTS_DIRECTORY_COLUMN_NAME,
                    PREPARED_DIRECTORY_COLUMN_NAME,
                ):
                    merged[col] = fresh_rec[col]
                pipeline_record_list.append(merged)
            else:
                _screen_notice(f"no existing state for {pdb_id} — starting from scratch")
                pipeline_record_list.append(fresh_rec)
        _screen_notice(f"rerun mode: from step {from_step}" + (f", protein {target_pdb_id}" if target_pdb_id else ""))
    else:
        pipeline_record_list = fresh_record_list

    _screen_notice(f"active pipeline records: {len(pipeline_record_list)}")
    _print_run_header(
        protein_data_dir=protein_data_dir,
        n_records=len(pipeline_record_list),
        target_pdb_id=target_pdb_id,
        from_step=from_step,
        force_field=DEFAULT_GROMACS_FORCE_FIELD,
        water_model=DEFAULT_GROMACS_WATER_MODEL,
        ph=FRUTON_PROTONATION_PH,
    )
    _append_fruton_log(
        "step_3:build_pipeline_state",
        [f"record_count: {len(pipeline_record_list)}", f"from_step: {from_step}", f"target_pdb_id: {target_pdb_id}"],
    )

    # Write rerun marker for each protein being rerun so the report can show the red box.
    if from_step is not None:
        _write_rerun_markers(pipeline_record_list, from_step)

    _screen_step(4, "fasta_files")
    if _step_runs(4):
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
    else:
        _screen_notice(f"step 4 skipped (rerunning from step {from_step})")

    _screen_step(5, "sequence_alignment")
    if _step_runs(5):
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
    else:
        _screen_notice(f"step 5 skipped (rerunning from step {from_step})")

    _screen_step(6, "insertion_codes")
    if _step_runs(6):
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
    else:
        _screen_notice(f"step 6 skipped (rerunning from step {from_step})")

    _screen_step(7, "representative_unit and component_split")
    if _step_runs(7):
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

            if component_summary.get("n_water_lines", 0) > 0:
                try:
                    _run_mutating_call(
                        log_title=f"step_8:crystal_water_preparation:{pdb_id}:captured_output",
                        work_label="crystal_water_preparation",
                        screen_pdb_id=pdb_id,
                        func=_prepare_crystal_water_component,
                        pdb_id=pdb_id,
                        pdb_dir=pdb_dir,
                    )
                    _screen_item(f"crystal_water_preparation -> {pdb_id}: SOL (OW+HW1+HW2)")
                except Exception as error:
                    _screen_error(f"crystal_water_preparation failed for {pdb_id}: {error!r}")
                    _log_fruton_exception(f"step_8:crystal_water_preparation:{pdb_id}", error)

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
    else:
        _screen_notice(f"step 7 skipped (rerunning from step {from_step})")

    _screen_step(8, "gap_detection")
    if _step_runs(8):
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
            gap_windows = gap_summary.get("gaps", [])
            pipeline_record[GAP_BOUNDARIES_COLUMN_NAME] = (
                json.dumps(gap_windows) if gap_windows else ""
            )

            max_gap_size = max(gap_sizes) if gap_sizes else 0
            if max_gap_size > 5:
                summary["max_gap_gt_5"] += 1
            if max_gap_size >= 8:
                summary["max_gap_ge_8"] += 1
    else:
        _screen_notice(f"step 8 skipped (rerunning from step {from_step})")

    _screen_step(9, "filler")
    prepared_variant_results_by_pdb: dict[str, list[dict[str, str]]] = {}
    if _step_runs(9):
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
    else:
        _screen_notice(f"step 9 skipped (rerunning from step {from_step})")

    _screen_step(10, "GROMACS protonation + prepared variants")
    if _step_runs(10):
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
                pipeline_record[PROTONATION_PROPKA_HIS_COLUMN_NAME] = route_result.get(
                    "protonation_propka_his_assignments", ""
                )

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
                    if len(final_input_path_list) == 1:
                        wat_output_path = _run_mutating_call(
                            log_title=f"step_11:prepared_structure_wat:{pdb_id}:{variant_label}:captured_output",
                            work_label=f"prepared_wat:{variant_label}",
                            screen_pdb_id=pdb_id,
                            func=_build_prepared_wat_structure_for_variant,
                            pdb_id=pdb_id,
                            pdb_dir=pdb_dir,
                            variant_label=variant_label,
                            final_protein_input_path=final_input_path_list[0],
                            final_protein_input_paths=None,
                        )
                    else:
                        wat_output_path = _run_mutating_call(
                            log_title=f"step_11:prepared_structure_wat:{pdb_id}:{variant_label}:captured_output",
                            work_label=f"prepared_wat:{variant_label}",
                            screen_pdb_id=pdb_id,
                            func=_build_prepared_wat_structure_for_variant,
                            pdb_id=pdb_id,
                            pdb_dir=pdb_dir,
                            variant_label=variant_label,
                            final_protein_input_path=None,
                            final_protein_input_paths=final_input_path_list,
                        )
                    if wat_output_path is not None:
                        _screen_item(f"prepared_structure_wat -> {pdb_id} [{variant_label}]: {Path(wat_output_path).name}")
                    else:
                        _screen_item(f"prepared_structure_wat -> {pdb_id} [{variant_label}]: skipped (no crystal waters)")
                except Exception as error:
                    _screen_error(f"prepared_structure_wat failed for {pdb_id} [{variant_label}]: {error!r}")
                    _log_fruton_exception(
                        f"step_11:prepared_structure_wat:{pdb_id}:{variant_label}",
                        error,
                    )

                try:
                    if len(final_input_path_list) == 1:
                        cof_output_path = _run_mutating_call(
                            log_title=f"step_11:prepared_structure_cofactor:{pdb_id}:{variant_label}:captured_output",
                            work_label=f"prepared_cofactor:{variant_label}",
                            screen_pdb_id=pdb_id,
                            func=_build_prepared_cofactor_structure_for_variant,
                            pdb_id=pdb_id,
                            pdb_dir=pdb_dir,
                            variant_label=variant_label,
                            final_protein_input_path=final_input_path_list[0],
                            final_protein_input_paths=None,
                        )
                    else:
                        cof_output_path = _run_mutating_call(
                            log_title=f"step_11:prepared_structure_cofactor:{pdb_id}:{variant_label}:captured_output",
                            work_label=f"prepared_cofactor:{variant_label}",
                            screen_pdb_id=pdb_id,
                            func=_build_prepared_cofactor_structure_for_variant,
                            pdb_id=pdb_id,
                            pdb_dir=pdb_dir,
                            variant_label=variant_label,
                            final_protein_input_path=None,
                            final_protein_input_paths=final_input_path_list,
                        )
                    if cof_output_path is not None:
                        _screen_item(f"prepared_structure_cofactor -> {pdb_id} [{variant_label}]: {Path(cof_output_path).name}")
                    else:
                        _screen_item(f"prepared_structure_cofactor -> {pdb_id} [{variant_label}]: skipped (no cofactors/ligands)")
                except Exception as error:
                    _screen_error(f"prepared_structure_cofactor failed for {pdb_id} [{variant_label}]: {error!r}")
                    _log_fruton_exception(
                        f"step_11:prepared_structure_cofactor:{pdb_id}:{variant_label}",
                        error,
                    )

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

                if model_source in ("modeller", "alphafold"):
                    try:
                        clash_result = _check_water_clashes_for_prepared_variant(
                            pdb_id=pdb_id,
                            variant_label=variant_label,
                            prepared_output_path=Path(prepared_output_path),
                        )
                        n_hard = clash_result.get("n_hard_clashes", 0)
                        n_soft = clash_result.get("n_soft_clashes", 0)
                        if n_hard or n_soft:
                            _screen_item(
                                f"water_clash_check -> {pdb_id} [{variant_label}]: "
                                f"{n_hard} hard, {n_soft} soft"
                            )
                    except Exception as error:
                        _log_fruton_exception(
                            f"step_11:water_clash_check:{pdb_id}:{variant_label}",
                            error,
                        )

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
    else:
        _screen_notice(f"step 10 skipped (rerunning from step {from_step})")

    _screen_step(11, "sanitize variants")
    sanitized_variant_results_by_pdb: dict[str, list[dict[str, str]]] = {}
    if _step_runs(11):

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
    else:
        _screen_notice(f"step 11 skipped (rerunning from step {from_step})")

    _screen_step(12, "metall_params")
    if _step_runs(12):
        for pipeline_record in pipeline_record_list:
            pdb_id = pipeline_record[PDB_ID_COLUMN_NAME]
            pdb_dir = Path(pipeline_record[PDB_DIRECTORY_COLUMN_NAME])

            has_metals = str(pipeline_record.get(HAS_METALS_COLUMN_NAME, "")).strip().lower()
            if has_metals != "yes":
                _screen_item(f"metall_params skipped for {pdb_id}: no metals")
                continue

            try:
                metall_params_result = _run_mutating_call(
                    log_title=f"step_13:metall_params:{pdb_id}:captured_output",
                    work_label=f"metall_params:{pdb_id}",
                    screen_pdb_id=pdb_id,
                    func=_run_metall_params_for_protein,
                    pdb_id=pdb_id,
                    pdb_dir=pdb_dir,
                )
            except Exception as error:
                _screen_error(f"metall_params failed for {pdb_id}: {error!r}")
                _log_fruton_exception(f"step_13:metall_params:{pdb_id}", error)
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
            pipeline_record[METALL_PARAMS_DIRECTORY_COLUMN_NAME] = str(
                metall_params_result.get("metall_params_dir", "") or ""
            )
            site_count = metall_params_result.get("transition_metal_site_count", 0)
            status = metall_params_result.get("status", "")
            _screen_item(f"metall_params -> {pdb_id}: {status} ({site_count} transition-metal site(s))")
    else:
        _screen_notice(f"step 12 skipped (rerunning from step {from_step})")

    _screen_step(13, "nonstd_residue_params")
    if _step_runs(13):
        for pipeline_record in pipeline_record_list:
            pdb_id = pipeline_record[PDB_ID_COLUMN_NAME]
            pdb_dir = Path(pipeline_record[PDB_DIRECTORY_COLUMN_NAME])

            has_nonstd = str(pipeline_record.get(HAS_NONSTANDARD_RESIDUES_COLUMN_NAME, "")).strip().lower()
            if has_nonstd != "yes":
                _screen_item(f"nonstd_residue_params skipped for {pdb_id}: no non-standard residues")
                continue

            try:
                nonstd_result = _run_mutating_call(
                    log_title=f"step_13:nonstd_residue_params:{pdb_id}:captured_output",
                    work_label=f"nonstd_residue_params:{pdb_id}",
                    screen_pdb_id=pdb_id,
                    func=_run_nonstd_residue_params_for_protein,
                    pdb_id=pdb_id,
                    pdb_dir=pdb_dir,
                )
            except Exception as error:
                _screen_error(f"nonstd_residue_params failed for {pdb_id}: {error!r}")
                _log_fruton_exception(f"step_13:nonstd_residue_params:{pdb_id}", error)
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
    else:
        _screen_notice(f"step 13 skipped (rerunning from step {from_step})")

    _screen_step(14, "gaussian_parametrization (MCPB + RESP, HPC — all proteins at once)")
    if _step_runs(14):

        # Ensure RESP scaffold exists for any has_nonstd protein that is missing it
        # (e.g. when re-running from step 14 after deleting files and step 13 was skipped).
        for pipeline_record in pipeline_record_list:
            pdb_id = pipeline_record[PDB_ID_COLUMN_NAME]
            pdb_dir = Path(pipeline_record[PDB_DIRECTORY_COLUMN_NAME])
            has_nonstd = str(pipeline_record.get(HAS_NONSTANDARD_RESIDUES_COLUMN_NAME, "")).strip().lower() == "yes"
            if not has_nonstd:
                continue
            nonstd_dir = pdb_dir / "nonstandard_params"
            resp_dirs_exist = any(
                (label_dir / "resp" / "commands.sh").is_file()
                for label_dir in nonstd_dir.iterdir()
                if label_dir.is_dir()
            ) if nonstd_dir.is_dir() else False
            if not resp_dirs_exist:
                _screen_item(f"nonstd_residue_params re-run for {pdb_id}: resp scaffold missing")
                try:
                    _run_mutating_call(
                        log_title=f"step_14:nonstd_residue_params_recovery:{pdb_id}:captured_output",
                        work_label=f"nonstd_residue_params_recovery:{pdb_id}",
                        screen_pdb_id=pdb_id,
                        func=_run_nonstd_residue_params_for_protein,
                        pdb_id=pdb_id,
                        pdb_dir=pdb_dir,
                    )
                except Exception as error:
                    _screen_error(f"nonstd_residue_params recovery failed for {pdb_id}: {error!r}")

        # Collect all proteins that need Gaussian (metals or non-standard residues)
        gaussian_proteins: list[tuple[Path, str]] = []
        for pipeline_record in pipeline_record_list:
            pdb_id = pipeline_record[PDB_ID_COLUMN_NAME]
            pdb_dir = Path(pipeline_record[PDB_DIRECTORY_COLUMN_NAME])
            has_metals = str(pipeline_record.get(HAS_METALS_COLUMN_NAME, "")).strip().lower() == "yes"
            has_nonstd = str(pipeline_record.get(HAS_NONSTANDARD_RESIDUES_COLUMN_NAME, "")).strip().lower() == "yes"
            needs_gaussian = has_metals or (
                has_nonstd and _protein_has_nonstd_needing_gaussian(pipeline_record)
            )
            if needs_gaussian:
                gaussian_proteins.append((pdb_dir, pdb_id))
            else:
                _screen_item(
                    f"gaussian_parametrization skipped for {pdb_id}: "
                    f"no metals or non-standard residues requiring Gaussian"
                )

        if gaussian_proteins:
            slurm_config = SlurmConfig(
                partition=args.slurm_partition,
                time=args.slurm_time,
                mem=args.slurm_mem,
                cpus=args.slurm_cpus,
                gpus=args.slurm_gpus,
            )
            try:
                all_gauss_results = _run_mutating_call(
                    log_title="step_15:gaussian_parametrization:all:captured_output",
                    work_label="gaussian_parametrization:all_proteins",
                    screen_pdb_id="[all]",
                    func=run_gaussian_parametrization_for_all_proteins,
                    proteins=gaussian_proteins,
                    slurm_config=slurm_config,
                )
            except Exception as error:
                _screen_error(f"gaussian_parametrization failed: {error!r}")
                _log_fruton_exception("step_15:gaussian_parametrization:all", error)
                for _, pdb_id in gaussian_proteins:
                    for rec in pipeline_record_list:
                        if rec[PDB_ID_COLUMN_NAME] == pdb_id:
                            rec[GAUSSIAN_PARAMS_STATUS_COLUMN_NAME] = STATUS_FAILED
                all_gauss_results = None

            if all_gauss_results is not None:
                for pipeline_record in pipeline_record_list:
                    pdb_id = pipeline_record[PDB_ID_COLUMN_NAME]
                    if pdb_id not in {pid for _, pid in gaussian_proteins}:
                        continue
                    gauss_result = all_gauss_results.get(pdb_id, {})
                    status = str(gauss_result.get("status", "")).strip()
                    pipeline_record[GAUSSIAN_PARAMS_STATUS_COLUMN_NAME] = status
                    pipeline_record[GAUSSIAN_PARAMS_MANIFEST_PATH_COLUMN_NAME] = str(
                        gauss_result.get("manifest_path", "") or ""
                    )
                    pipeline_record[GAUSSIAN_PARAMS_MCPB_CHARGE_COLUMN_NAME] = str(
                        gauss_result.get("mcpb_metal_charge", "") or ""
                    )
                    pipeline_record[GAUSSIAN_PARAMS_MCPB_SPIN_COLUMN_NAME] = str(
                        gauss_result.get("mcpb_metal_spin", "") or ""
                    )
                    pipeline_record[GAUSSIAN_PARAMS_RESP_CHARGE_COLUMN_NAME] = str(
                        gauss_result.get("resp_charge", "") or ""
                    )
                    pipeline_record[GAUSSIAN_PARAMS_FRCMOD_PATHS_COLUMN_NAME] = str(
                        gauss_result.get("frcmod_paths", "") or ""
                    )
                    msg = gauss_result.get("message", "")
                    if status == STATUS_REQUIRED:
                        _screen_notice(f"gaussian_parametrization -> {pdb_id}: input files prepared — submit on CESGA to run Gaussian")
                    else:
                        _screen_item(f"gaussian_parametrization -> {pdb_id}: {status}" + (f" — {msg}" if msg else ""))
    else:
        _screen_notice(f"step 14 skipped (rerunning from step {from_step})")

    # --- Gaussian pending gate -------------------------------------------
    # If any protein still has STATUS_REQUIRED (Gaussian inputs prepared but
    # jobs not yet submitted / run_after_gaussian.sh not yet executed), halt
    # the pipeline here.  Save a checkpoint so no work is lost.  Steps 15-17
    # (parameter_audit, model_evaluation, protein_report) require finalized
    # parameters and must not run on incomplete data.
    #
    # Before halting, verify that frcmod+mol2 files do not already exist on
    # disk.  If run_after_gaussian.sh has already completed (e.g. when resuming
    # with --from-step 15), the status column is stale — update it to
    # STATUS_SUCCESS so the pipeline can continue without re-running step 14.
    for _r in pipeline_record_list:
        if str(_r.get(GAUSSIAN_PARAMS_STATUS_COLUMN_NAME, "")).strip() != STATUS_REQUIRED:
            continue
        _manifest_path = str(_r.get(NONSTD_RESIDUE_PARAMS_MANIFEST_PATH_COLUMN_NAME, "") or "")
        _unsupported = tuple(
            n.strip()
            for n in str(_r.get("parameter_audit.unsupported_residues", "") or "").split("|")
            if n.strip()
        )
        if _unsupported and _manifest_path:
            _frcmod_found, _mol2_found = _find_nonstd_resp_param_files(
                manifest_path=_manifest_path,
                unsupported_residue_names=_unsupported,
            )
            if len(_frcmod_found) == len(_unsupported) and len(_mol2_found) == len(_unsupported):
                _r[GAUSSIAN_PARAMS_STATUS_COLUMN_NAME] = STATUS_SUCCESS
                _screen_item(
                    f"gaussian_parametrization -> {_r[PDB_ID_COLUMN_NAME]}: "
                    f"frcmod+mol2 found on disk — marking complete"
                )

    _gaussian_pending_ids = [
        r[PDB_ID_COLUMN_NAME]
        for r in pipeline_record_list
        if str(r.get(GAUSSIAN_PARAMS_STATUS_COLUMN_NAME, "")).strip() == STATUS_REQUIRED
    ]
    if _gaussian_pending_ids:
        _screen_notice(
            "Gaussian jobs required for: " + ", ".join(_gaussian_pending_ids)
        )
        _screen_notice(
            "Pipeline halted — Gaussian must complete before the pipeline can continue.\n"
            "  1. Submit jobs:  run each protein's submit_gaussian.sh (RESP) or\n"
            "                   MCPB_submit.sh (metals) on the HPC cluster\n"
            "  2. After jobs finish: run run_after_gaussian.sh in each\n"
            "                   nonstandard_params/<label>/resp/ directory\n"
            "  3. Resume:       fruton --from-step 15  (or full re-run)"
        )
        _screen_step(19, "save_pipeline_json (checkpoint)")
        _records_to_save = pipeline_record_list
        if target_pdb_id is not None and pipeline_json_path.is_file():
            _existing = load_pipeline_table(pipeline_json_path)
            _updated_ids = {r[PDB_ID_COLUMN_NAME] for r in pipeline_record_list}
            _records_to_save = (
                [r for r in _existing if r[PDB_ID_COLUMN_NAME] not in _updated_ids]
                + list(pipeline_record_list)
            )
        _run_mutating_call(
            log_title="step_19:save_pipeline_json:captured_output",
            work_label="save_pipeline_json",
            screen_pdb_id="-",
            func=save_pipeline_table,
            protein_record_list=_records_to_save,
            json_path=pipeline_json_path,
        )
        _screen_step(20, "write_pipeline_xlsx (checkpoint)")
        _run_mutating_call(
            log_title="step_20:write_pipeline_xlsx:captured_output",
            work_label="write_pipeline_xlsx",
            screen_pdb_id="-",
            func=write_pipeline_to_xlsx,
            protein_record_list=pipeline_record_list,
            output_path=pipeline_xlsx_path,
        )
        _PROGRESS.finish()
        _screen_notice(f"checkpoint JSON written: {pipeline_json_path}")
        _screen_notice(f"checkpoint XLSX written: {pipeline_xlsx_path}")
        _append_fruton_log(
            "run_pipeline:halted_for_gaussian",
            [
                f"gaussian_pending   : {', '.join(_gaussian_pending_ids)}",
                f"pipeline_json_path : {pipeline_json_path}",
                f"pipeline_xlsx_path : {pipeline_xlsx_path}",
                f"summary            : {summary}",
            ],
        )
        _print_summary(summary)
        print(f"  {_timestamp()}  ·  halted — Gaussian pending")
        return
    # ---------------------------------------------------------------------

    _screen_step(15, "parameter_audit of variants")
    if _step_runs(15):

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
            metall_scaffold_ready = (
                str(pipeline_record.get(METALL_PARAMS_STATUS_COLUMN_NAME, "")).strip().lower()
                == STATUS_SUCCESS
            )
            nonstd_scaffold_ready = (
                str(pipeline_record.get(NONSTD_RESIDUE_PARAMS_STATUS_COLUMN_NAME, "")).strip().lower()
                == STATUS_SUCCESS
            )

            for variant_result in sanitized_variant_result_list:
                variant_label = variant_result["variant_label"]
                sanitized_input_path = Path(variant_result["sanitized_audit_input_path"])

                try:
                    parameter_audit_result = _run_mutating_call(
                        log_title=f"step_16:parameter_audit_variant:{pdb_id}:{variant_label}:captured_output",
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
                        f"step_16:parameter_audit_variant:{pdb_id}:{variant_label}",
                        error,
                    )
                    continue

                accepted = _variant_audit_accepts_end_model(parameter_audit_result)

                # Metal gate: if the metals_check vetoed this variant, check whether the
                # MCPB scaffold was generated in step 12.  If so, accept — Gaussian still
                # needs to run on HPC, but the preparation side is complete.
                if not _variant_metal_check_accepts_end_model(variant_result):
                    if metall_scaffold_ready:
                        accepted = True
                    else:
                        accepted = False

                # Non-standard residue gate (three tiers):
                #   1. frcmod+mol2 exist (Gaussian done) → tleap probe → verified accept
                #   2. scaffold ready, Gaussian pending   → provisional accept
                #   3. scaffold missing                  → reject
                resp_probe: RespParameterProbeResult | None = None
                if not accepted and getattr(parameter_audit_result, "requires_qm_parameters", False):
                    if nonstd_scaffold_ready:
                        unsupported_names = getattr(
                            parameter_audit_result, "unsupported_residue_names", ()
                        )
                        nonstd_manifest_path = str(
                            pipeline_record.get(NONSTD_RESIDUE_PARAMS_MANIFEST_PATH_COLUMN_NAME, "") or ""
                        )
                        frcmod_paths, mol2_paths = _find_nonstd_resp_param_files(
                            manifest_path=nonstd_manifest_path,
                            unsupported_residue_names=unsupported_names,
                        )
                        if unsupported_names:
                            try:
                                resp_probe = probe_amber_resp_parameters(
                                    unsupported_residue_names=unsupported_names,
                                    frcmod_paths=frcmod_paths,
                                    mol2_paths=mol2_paths,
                                )
                            except Exception as probe_exc:
                                _log_fruton_exception(
                                    f"step_15:resp_probe:{pdb_id}:{variant_label}",
                                    probe_exc,
                                )
                        accepted = True

                _append_parameter_audit_variant_summary(
                    pdb_id=pdb_id,
                    variant_label=variant_label,
                    audit_result=parameter_audit_result,
                    accepted=accepted,
                )

                if accepted:
                    accepted_variant_result_list.append(variant_result)
                    if metall_scaffold_ready and not _variant_metal_check_accepts_end_model(variant_result):
                        _screen_item(
                            f"parameter_audit variants -> {pdb_id} [{variant_label}]: accepted"
                            f" (MCPB scaffold generated — Gaussian pending)"
                        )
                    elif resp_probe is not None and getattr(parameter_audit_result, "requires_qm_parameters", False):
                        if resp_probe.parameters_complete and resp_probe.tleap_success:
                            _screen_item(
                                f"parameter_audit variants -> {pdb_id} [{variant_label}]: accepted"
                                f" (RESP parameters verified by tleap)"
                            )
                        elif resp_probe.parameters_complete and not resp_probe.tleap_success:
                            _screen_notice(
                                f"parameter_audit variants -> {pdb_id} [{variant_label}]: accepted"
                                f" (RESP files present but tleap flagged issues"
                                + (f" — {resp_probe.tleap_probe.error_message}" if resp_probe.tleap_probe else "")
                                + ")"
                            )
                        else:
                            _screen_item(
                                f"parameter_audit variants -> {pdb_id} [{variant_label}]: accepted"
                                f" (RESP scaffold generated — Gaussian pending)"
                            )
                    elif nonstd_scaffold_ready and getattr(parameter_audit_result, "requires_qm_parameters", False):
                        _screen_item(
                            f"parameter_audit variants -> {pdb_id} [{variant_label}]: accepted"
                            f" (RESP scaffold generated — Gaussian pending)"
                        )
                    else:
                        _screen_item(f"parameter_audit variants -> {pdb_id} [{variant_label}]: accepted")
                elif getattr(parameter_audit_result, "requires_metal_parameters", False):
                    _screen_notice(
                        f"parameter_audit variants -> {pdb_id} [{variant_label}]: metal parameters required"
                        f" — run metall_params first"
                    )
                elif getattr(parameter_audit_result, "requires_qm_parameters", False):
                    _screen_notice(
                        f"parameter_audit variants -> {pdb_id} [{variant_label}]: unsupported residue parameters required"
                        f" — run nonstd_residue_params first"
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
    else:
        _screen_notice(f"step 15 skipped (rerunning from step {from_step})")

    _screen_step(16, "stage_amber_parameters")
    if _step_runs(16):
        _amber_dat_dir = _find_amber_dat_dir()
        if _amber_dat_dir is None:
            _screen_notice(
                "stage_amber_parameters: AMBER dat/leap directory not found — skipping"
            )
        else:
            for pipeline_record in pipeline_record_list:
                pdb_id = pipeline_record[PDB_ID_COLUMN_NAME]
                pdb_dir = Path(pipeline_record[PDB_DIRECTORY_COLUMN_NAME])
                try:
                    stage_result = _run_mutating_call(
                        log_title=f"step_16:stage_amber_parameters:{pdb_id}:captured_output",
                        work_label=f"stage_amber_parameters:{pdb_id}",
                        screen_pdb_id=pdb_id,
                        func=_stage_amber_parameters_for_protein,
                        pdb_id=pdb_id,
                        pdb_dir=pdb_dir,
                        pipeline_record=pipeline_record,
                        amber_dat_dir=_amber_dat_dir,
                    )
                except Exception as error:
                    _screen_error(f"stage_amber_parameters failed for {pdb_id}: {error!r}")
                    _log_fruton_exception(f"step_16:stage_amber_parameters:{pdb_id}", error)
                    continue
                status = stage_result.get("status", "")
                msg = stage_result.get("message", "")
                _screen_item(
                    f"stage_amber_parameters -> {pdb_id}: {status}"
                    + (f" — {msg}" if msg else "")
                )
    else:
        _screen_notice(f"step 16 skipped (rerunning from step {from_step})")

    _screen_step(17, "model_evaluation (Modeller-filled only)")
    if _step_runs(17):
        for pipeline_record in pipeline_record_list:
            pdb_id = pipeline_record[PDB_ID_COLUMN_NAME]
            pdb_dir = Path(pipeline_record[PDB_DIRECTORY_COLUMN_NAME])
            try:
                eval_result = _run_mutating_call(
                    log_title=f"step_17:model_evaluation:{pdb_id}:captured_output",
                    work_label=f"model_evaluation:{pdb_id}",
                    screen_pdb_id=pdb_id,
                    func=run_model_evaluation,
                    protein_dir=pdb_dir,
                    pdb_id=pdb_id,
                    pipeline_record=pipeline_record,
                )
            except Exception as error:
                _screen_error(f"model_evaluation failed for {pdb_id}: {error!r}")
                _log_fruton_exception(f"step_17:model_evaluation:{pdb_id}", error)
                pipeline_record[MODEL_EVALUATION_STATUS_COLUMN_NAME] = STATUS_FAILED
                continue

            status = str(eval_result.get("status", "")).strip()
            pipeline_record[MODEL_EVALUATION_STATUS_COLUMN_NAME] = status
            pipeline_record[MODEL_EVALUATION_PCT_FAVORED_COLUMN_NAME] = (
                str(eval_result.get("pct_favored", "")) if eval_result.get("pct_favored") is not None else ""
            )
            pipeline_record[MODEL_EVALUATION_PCT_OUTLIER_COLUMN_NAME] = (
                str(eval_result.get("pct_outlier", "")) if eval_result.get("pct_outlier") is not None else ""
            )
            pipeline_record[MODEL_EVALUATION_CLASHSCORE_COLUMN_NAME] = (
                str(eval_result.get("clashscore", "")) if eval_result.get("clashscore") is not None else ""
            )
            pipeline_record[MODEL_EVALUATION_OVERALL_QUALITY_COLUMN_NAME] = str(eval_result.get("overall_quality", ""))
            pipeline_record[MODEL_EVALUATION_RAMA_PLOT_PATH_COLUMN_NAME] = str(eval_result.get("rama_plot_path", ""))
            pipeline_record[MODEL_EVALUATION_MANIFEST_PATH_COLUMN_NAME] = str(eval_result.get("manifest_path", ""))

            if status == "skipped":
                _screen_item(f"model_evaluation -> {pdb_id}: skipped (no Modeller gap-filling)")
            elif status == "success":
                pct_fav = eval_result.get("pct_favored")
                quality = eval_result.get("overall_quality", "")
                _screen_item(f"model_evaluation -> {pdb_id}: {quality} (favored {pct_fav}%)")
            else:
                _screen_item(f"model_evaluation -> {pdb_id}: {status} — {eval_result.get('message', '')}")

        global_bib_path = protein_data_dir / "references.bib"
    else:
        _screen_notice(f"step 17 skipped (rerunning from step {from_step})")

    _screen_step(18, "protein_report")
    if _step_runs(18):
        for pipeline_record in pipeline_record_list:
            pdb_id = pipeline_record[PDB_ID_COLUMN_NAME]
            pdb_dir = Path(pipeline_record[PDB_DIRECTORY_COLUMN_NAME])
            try:
                report_result = _run_mutating_call(
                    log_title=f"step_18:protein_report:{pdb_id}:captured_output",
                    work_label=f"protein_report:{pdb_id}",
                    screen_pdb_id=pdb_id,
                    func=generate_protein_report,
                    protein_dir=pdb_dir,
                    pdb_id=pdb_id,
                    pipeline_record=pipeline_record,
                    global_bib_path=global_bib_path,
                )
            except Exception as error:
                _screen_error(f"protein_report failed for {pdb_id}: {error!r}")
                _log_fruton_exception(f"step_18:protein_report:{pdb_id}", error)
                pipeline_record[REPORT_STATUS_COLUMN_NAME] = STATUS_FAILED
                continue

            status = str(report_result.get("status", "")).strip()
            pipeline_record[REPORT_STATUS_COLUMN_NAME] = status
            pipeline_record[REPORT_PDF_PATH_COLUMN_NAME] = str(report_result.get("report_pdf_path", ""))
            pipeline_record[REPORT_FIGURE_PNG_PATH_COLUMN_NAME] = str(report_result.get("figure_png_path", ""))

            msg = report_result.get("message", "")
            _screen_item(f"protein_report -> {pdb_id}: {status}" + (f" — {msg}" if msg else ""))
    else:
        _screen_notice(f"step 18 skipped (rerunning from step {from_step})")

    _screen_step(19, "save_pipeline_json")
    # When running with --pdb-id, merge the updated record(s) back into the
    # full pipeline table so the rest of the proteins are not overwritten.
    _records_to_save = pipeline_record_list
    if target_pdb_id is not None and pipeline_json_path.is_file():
        _existing = load_pipeline_table(pipeline_json_path)
        _updated_ids = {r[PDB_ID_COLUMN_NAME] for r in pipeline_record_list}
        _records_to_save = [r for r in _existing if r[PDB_ID_COLUMN_NAME] not in _updated_ids] + list(pipeline_record_list)
    _run_mutating_call(log_title="step_19:save_pipeline_json:captured_output", work_label="save_pipeline_json", screen_pdb_id="-", func=save_pipeline_table, protein_record_list=_records_to_save, json_path=pipeline_json_path)
    _screen_step(20, "write_pipeline_xlsx")
    _run_mutating_call(log_title="step_20:write_pipeline_xlsx:captured_output", work_label="write_pipeline_xlsx", screen_pdb_id="-", func=write_pipeline_to_xlsx, protein_record_list=pipeline_record_list, output_path=pipeline_xlsx_path)
    _PROGRESS.finish()

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
    print(f"  {_timestamp()}  ·  done")


def main() -> None:
    run_pipeline()


if __name__ == "__main__":
    main()
