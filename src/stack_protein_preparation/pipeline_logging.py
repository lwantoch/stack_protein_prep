# /home/grheco/repositorios/stack_protein_prep/src/stack_protein_preparation/pipeline_logging.py

"""
Shared logging and terminal-rendering utilities for the FRUTON pipeline.

Purpose
-------
- create one consistent logging system for all pipeline modules
- write per-protein / per-module / per-variant log files
- capture Python stdout/stderr into log files
- capture external subprocess stdout/stderr into log files
- provide structured terminal output with varied box-drawing styles

Design goals
------------
- screen output should stay compact and readable
- log files should be verbose and useful for debugging
- each module should be able to log independently
- no hidden global state
- safe to import from any pipeline module

Directory convention
--------------------
Logs are written under:

    data/proteins/logs/

Examples:
    data/proteins/logs/fruton/pipeline.log
    data/proteins/logs/1UOU/protonation.gaps.log
    data/proteins/logs/2P1T/filler.log

Recommended use
---------------
1. Build a log path for the current module call:
       log_paths = build_module_log_paths(...)

2. Write a header:
       write_log_header(...)

3. Print a compact terminal box:
       print_section_box(...)
       print_step_start_box(...)
       print_step_end_box(...)

4. Capture Python output:
       with capture_python_output_to_log(log_path):
           ...

5. Capture subprocess output:
       run_subprocess_logged(...)

Important
---------
This module does NOT know pipeline business logic.
It only handles logging and terminal presentation.
"""

from __future__ import annotations

import contextlib
import io
import os
import subprocess
import traceback
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import Iterable

# ---------------------------------------------------------------------------
# Data containers
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class ModuleLogPaths:
    """
    Resolved log paths for one module invocation.
    """

    logs_root_dir: Path
    protein_log_dir: Path
    module_log_path: Path


# ---------------------------------------------------------------------------
# Time helpers
# ---------------------------------------------------------------------------


def current_timestamp() -> str:
    """
    Return a human-readable local timestamp.
    """
    return datetime.now().strftime("%Y-%m-%d %H:%M:%S")


def current_timestamp_for_filename() -> str:
    """
    Return a compact timestamp useful for file naming if ever needed later.
    """
    return datetime.now().strftime("%Y%m%d_%H%M%S")


# ---------------------------------------------------------------------------
# Log path helpers
# ---------------------------------------------------------------------------


def ensure_logs_root_dir(protein_data_dir: Path) -> Path:
    """
    Ensure the central pipeline logs root directory exists.

    Expected result:
        data/proteins/logs
    """
    logs_root_dir = protein_data_dir / "logs"
    logs_root_dir.mkdir(parents=True, exist_ok=True)
    return logs_root_dir


def ensure_global_log_dir(protein_data_dir: Path, group_name: str) -> Path:
    """
    Ensure a named global log directory exists.

    Example:
        data/proteins/logs/fruton
    """
    logs_root_dir = ensure_logs_root_dir(protein_data_dir)
    global_log_dir = logs_root_dir / group_name
    global_log_dir.mkdir(parents=True, exist_ok=True)
    return global_log_dir


def build_module_log_paths(
    protein_data_dir: Path,
    pdb_id: str,
    module_name: str,
    variant_label: str | None = None,
) -> ModuleLogPaths:
    """
    Build per-protein log paths for one module invocation.

    Examples
    --------
    pdb_id='1UOU', module_name='protonation', variant_label='gaps'
        -> data/proteins/logs/1UOU/protonation.gaps.log

    pdb_id='2P1T', module_name='filler', variant_label=None
        -> data/proteins/logs/2P1T/filler.log
    """
    normalized_pdb_id = str(pdb_id).strip().upper()
    normalized_module_name = str(module_name).strip()

    if not normalized_pdb_id:
        raise ValueError("pdb_id must not be empty")
    if not normalized_module_name:
        raise ValueError("module_name must not be empty")

    logs_root_dir = ensure_logs_root_dir(protein_data_dir)
    protein_log_dir = logs_root_dir / normalized_pdb_id
    protein_log_dir.mkdir(parents=True, exist_ok=True)

    if variant_label is None or not str(variant_label).strip():
        log_filename = f"{normalized_module_name}.log"
    else:
        log_filename = f"{normalized_module_name}.{str(variant_label).strip()}.log"

    return ModuleLogPaths(
        logs_root_dir=logs_root_dir,
        protein_log_dir=protein_log_dir,
        module_log_path=protein_log_dir / log_filename,
    )


def build_global_log_path(
    protein_data_dir: Path,
    group_name: str,
    log_filename: str,
) -> Path:
    """
    Build a global log path under a named group directory.

    Example:
        data/proteins/logs/fruton/pipeline.log
    """
    global_log_dir = ensure_global_log_dir(protein_data_dir, group_name)
    return global_log_dir / log_filename


# ---------------------------------------------------------------------------
# Low-level log writing
# ---------------------------------------------------------------------------


def append_log_text(log_path: Path, text: str) -> None:
    """
    Append raw text to a log file.

    Behavior
    --------
    - creates parent directories if needed
    - writes nothing for empty text
    - always terminates with a newline
    """
    if not text:
        return

    log_path.parent.mkdir(parents=True, exist_ok=True)

    with log_path.open("a", encoding="utf-8") as handle:
        if text.endswith("\n"):
            handle.write(text)
        else:
            handle.write(text + "\n")


def append_log_lines(log_path: Path, line_list: Iterable[str]) -> None:
    """
    Append multiple lines to a log file.
    """
    for line in line_list:
        append_log_text(log_path, line)


def append_exception_traceback(log_path: Path, error: BaseException) -> None:
    """
    Append a full traceback for an exception to a log file.
    """
    traceback_text = "".join(
        traceback.format_exception(type(error), error, error.__traceback__)
    )
    append_log_text(log_path, "[EXCEPTION]")
    append_log_text(log_path, traceback_text.rstrip())


def write_log_header(
    log_path: Path,
    *,
    module_name: str,
    pdb_id: str | None = None,
    variant_label: str | None = None,
    extra_lines: Iterable[str] | None = None,
) -> None:
    """
    Append a standard section header to a log file.
    """
    header_line = "=" * 80
    descriptor_parts = [f"[{current_timestamp()}]", f"module={module_name}"]

    if pdb_id is not None and str(pdb_id).strip():
        descriptor_parts.append(f"pdb_id={str(pdb_id).strip().upper()}")

    if variant_label is not None and str(variant_label).strip():
        descriptor_parts.append(f"variant={str(variant_label).strip()}")

    append_log_lines(
        log_path,
        [
            header_line,
            " ".join(descriptor_parts),
            header_line,
        ],
    )

    if extra_lines:
        append_log_lines(log_path, extra_lines)


def append_key_value_block(
    log_path: Path,
    title: str,
    key_value_pairs: list[tuple[str, str]],
) -> None:
    """
    Append a small titled key/value block to a log file.
    """
    append_log_text(log_path, f"[{title}]")
    for key, value in key_value_pairs:
        append_log_text(log_path, f"{key}: {value}")
    append_log_text(log_path, "")


def append_separator(log_path: Path, character: str = "-", width: int = 80) -> None:
    """
    Append a simple separator line to a log file.
    """
    append_log_text(log_path, character * width)


# ---------------------------------------------------------------------------
# Python stdout/stderr capture
# ---------------------------------------------------------------------------


@contextlib.contextmanager
def capture_python_output_to_log(log_path: Path):
    """
    Capture Python stdout/stderr during a block and append both to a log file.

    Example
    -------
        with capture_python_output_to_log(log_path):
            run_alignments_for_pdb_directory(...)
    """
    stdout_buffer = io.StringIO()
    stderr_buffer = io.StringIO()

    with (
        contextlib.redirect_stdout(stdout_buffer),
        contextlib.redirect_stderr(stderr_buffer),
    ):
        yield

    stdout_text = stdout_buffer.getvalue().strip()
    stderr_text = stderr_buffer.getvalue().strip()

    if stdout_text:
        append_log_text(log_path, "[PYTHON STDOUT]")
        append_log_text(log_path, stdout_text)

    if stderr_text:
        append_log_text(log_path, "[PYTHON STDERR]")
        append_log_text(log_path, stderr_text)


# ---------------------------------------------------------------------------
# Subprocess capture
# ---------------------------------------------------------------------------


def run_subprocess_logged(
    command: list[str],
    *,
    log_path: Path,
    cwd: Path | None = None,
    env: dict[str, str] | None = None,
    check: bool = False,
) -> subprocess.CompletedProcess[str]:
    """
    Run a subprocess and append full execution details to the log file.

    Logged data
    -----------
    - exact command
    - working directory
    - selected environment overrides
    - return code
    - stdout
    - stderr

    Important
    ---------
    This function does not hide stderr.
    It is intended for maximal debug visibility.
    """
    if not command:
        raise ValueError("command must not be empty")

    append_log_text(log_path, "[SUBPROCESS COMMAND]")
    append_log_text(log_path, " ".join(command))

    if cwd is not None:
        append_log_text(log_path, f"[SUBPROCESS CWD] {cwd}")

    if env:
        append_log_text(log_path, "[SUBPROCESS ENV OVERRIDES]")
        for key in sorted(env):
            append_log_text(log_path, f"{key}={env[key]}")

    completed = subprocess.run(
        command,
        cwd=str(cwd) if cwd is not None else None,
        env=env,
        text=True,
        capture_output=True,
        check=False,
    )

    append_log_text(log_path, f"[SUBPROCESS RETURN CODE] {completed.returncode}")

    if completed.stdout:
        append_log_text(log_path, "[SUBPROCESS STDOUT]")
        append_log_text(log_path, completed.stdout.rstrip())

    if completed.stderr:
        append_log_text(log_path, "[SUBPROCESS STDERR]")
        append_log_text(log_path, completed.stderr.rstrip())

    if check and completed.returncode != 0:
        raise subprocess.CalledProcessError(
            returncode=completed.returncode,
            cmd=command,
            output=completed.stdout,
            stderr=completed.stderr,
        )

    return completed


# ---------------------------------------------------------------------------
# Terminal box styles
# ---------------------------------------------------------------------------


def _normalize_box_width(
    title: str,
    line_list: list[str],
    requested_inner_width: int | None = None,
) -> int:
    content_width = max([len(title)] + [len(line) for line in line_list] + [50])
    if requested_inner_width is None:
        return content_width + 2
    return max(requested_inner_width, content_width + 2)


def _content_line(
    text: str,
    inner_width: int,
    left_border: str,
    right_border: str,
) -> str:
    return f"{left_border} {text.ljust(inner_width - 2)} {right_border}"


def print_double_header_box(
    title: str,
    line_list: list[str],
    *,
    inner_width: int | None = None,
) -> None:
    """
    Print a box with a heavy outer frame and heavy title separator.
    """
    inner_width = _normalize_box_width(title, line_list, inner_width)

    print(f"╔{'═' * (inner_width + 2)}╗")
    print(_content_line(title, inner_width, "║", "║"))
    print(f"╠{'═' * (inner_width + 2)}╣")
    for line in line_list:
        print(_content_line(line, inner_width, "║", "║"))
    print(f"╚{'═' * (inner_width + 2)}╝")


def print_mixed_box(
    title: str,
    line_list: list[str],
    *,
    inner_width: int | None = None,
) -> None:
    """
    Print a box with a heavy top/header and lighter body separators.
    """
    inner_width = _normalize_box_width(title, line_list, inner_width)

    print(f"╔{'═' * (inner_width + 2)}╗")
    print(_content_line(title, inner_width, "║", "║"))
    print(f"╟{'─' * (inner_width + 2)}╢")
    for line in line_list:
        print(_content_line(line, inner_width, "│", "│"))
    print(f"╚{'═' * (inner_width + 2)}╝")


def print_light_table_box(
    title: str,
    row_list: list[tuple[str, str]],
    *,
    left_header: str = "Field",
    right_header: str = "Value",
) -> None:
    """
    Print a compact two-column table with lighter body grid.
    """
    left_width = max(len(left_header), *(len(left) for left, _ in row_list))
    right_width = max(len(right_header), *(len(right) for _, right in row_list))

    top_width = left_width + right_width + 7
    print(f"╔{'═' * top_width}╗")
    print(f"║ {title.ljust(top_width - 2)} ║")
    print(f"╠{'═' * (left_width + 2)}╤{'═' * (right_width + 2)}╣")
    print(f"║ {left_header.ljust(left_width)} │ {right_header.ljust(right_width)} ║")
    print(f"╟{'─' * (left_width + 2)}┼{'─' * (right_width + 2)}╢")

    for index, (left, right) in enumerate(row_list):
        print(f"║ {left.ljust(left_width)} │ {right.ljust(right_width)} ║")
        if index != len(row_list) - 1:
            print(f"╟{'─' * (left_width + 2)}┼{'─' * (right_width + 2)}╢")

    print(f"╚{'═' * (left_width + 2)}╧{'═' * (right_width + 2)}╝")


def print_round_box(
    title: str,
    line_list: list[str],
    *,
    inner_width: int | None = None,
) -> None:
    """
    Print a softer rounded box style for visual variety.
    """
    inner_width = _normalize_box_width(title, line_list, inner_width)

    print(f"╭{'─' * (inner_width + 2)}╮")
    print(_content_line(title, inner_width, "│", "│"))
    print(f"├{'─' * (inner_width + 2)}┤")
    for line in line_list:
        print(_content_line(line, inner_width, "│", "│"))
    print(f"╰{'─' * (inner_width + 2)}╯")


def print_banner(title: str) -> None:
    """
    Print a simple banner separator for major pipeline sections.
    """
    width = max(60, len(title) + 8)
    print(f"┏{'━' * width}┓")
    print(f"┃ {title.center(width - 2)} ┃")
    print(f"┗{'━' * width}┛")


def print_section_box(
    title: str,
    line_list: list[str] | None = None,
) -> None:
    """
    Generic section box expected by the module docstring examples.
    """
    print_round_box(title=title, line_list=line_list or [])


def print_step_start_box(
    *,
    module_name: str,
    pdb_id: str,
    variant_label: str | None = None,
    log_path: Path | None = None,
) -> None:
    """
    Print a compact START box for one module invocation.
    """
    line_list = [f"PDB ID   : {pdb_id.upper()}"]
    if variant_label:
        line_list.append(f"Variant  : {variant_label}")
    if log_path is not None:
        line_list.append(f"Log file : {log_path}")

    print_mixed_box(
        title=f"{module_name} · START",
        line_list=line_list,
    )


def print_step_end_box(
    *,
    module_name: str,
    pdb_id: str,
    status: str,
    detail_lines: list[str] | None = None,
) -> None:
    """
    Print a compact END box for one module invocation.
    """
    line_list = [
        f"PDB ID  : {pdb_id.upper()}",
        f"Status  : {status}",
    ]

    if detail_lines:
        line_list.extend(detail_lines)

    print_double_header_box(
        title=f"{module_name} · END",
        line_list=line_list,
    )


def print_status_table(
    title: str,
    row_list: list[tuple[str, str]],
) -> None:
    """
    Print a two-column status table.
    """
    print_light_table_box(
        title=title,
        row_list=row_list,
        left_header="Item",
        right_header="Status",
    )


# ---------------------------------------------------------------------------
# Convenience high-level helpers
# ---------------------------------------------------------------------------


def log_and_print_start(
    *,
    protein_data_dir: Path,
    pdb_id: str,
    module_name: str,
    variant_label: str | None = None,
    extra_log_lines: Iterable[str] | None = None,
) -> Path:
    """
    Build the log path, write a standard header, and print a START box.

    Returns
    -------
    Path
        The concrete module log path.
    """
    log_paths = build_module_log_paths(
        protein_data_dir=protein_data_dir,
        pdb_id=pdb_id,
        module_name=module_name,
        variant_label=variant_label,
    )

    write_log_header(
        log_path=log_paths.module_log_path,
        module_name=module_name,
        pdb_id=pdb_id,
        variant_label=variant_label,
        extra_lines=extra_log_lines,
    )

    print_step_start_box(
        module_name=module_name,
        pdb_id=pdb_id,
        variant_label=variant_label,
        log_path=log_paths.module_log_path,
    )

    return log_paths.module_log_path


def log_success(
    log_path: Path,
    message: str = "success",
) -> None:
    """
    Append a standardized success marker to a log file.
    """
    append_log_text(log_path, f"[RESULT] {message}")


def log_failure(
    log_path: Path,
    message: str,
    error: BaseException | None = None,
) -> None:
    """
    Append a standardized failure marker to a log file.
    """
    append_log_text(log_path, f"[RESULT] FAILED: {message}")
    if error is not None:
        append_exception_traceback(log_path, error)


def environment_snapshot_text() -> list[str]:
    """
    Return a small environment snapshot that can be useful in logs.

    Intended for optional use in selected modules.
    """
    return [
        f"timestamp={current_timestamp()}",
        f"cwd={Path.cwd()}",
        f"python_executable={os.sys.executable}",
    ]
