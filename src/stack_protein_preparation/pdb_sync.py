# /home/grheco/repositorios/stack_protein_prep/src/stack_protein_preparation/pdb_sync.py

"""
Synchronize a CSV file of PDB IDs with a protein data directory and download
range-filtered PDB files.

CSV convention
--------------
The CSV must contain a 'pdb_id' column.
It may optionally contain a 'range' column.

Example:

    pdb_id,range
    1ABC,10-280
    2XYZ,
    3DEF,5-190

Meaning of 'range'
------------------
The optional 'range' column describes the residue-number interval that should be
kept for polymer ATOM records after downloading the PDB file.

Current trimming behavior
-------------------------
- If range is empty:
    keep the downloaded PDB unchanged.
- If range is present, for example '20-280':
    keep only ATOM residues whose residue number is within the inclusive range.
- Waters are always kept, even if their residue number lies outside the range.
- Other HETATM records are currently kept unchanged.

Important clarification
-----------------------
A requested range such as '20-280' does NOT guarantee that the final local PDB
will start with residue 20.

Why
---
PDB coordinate files often do not contain all residues in the biological
sequence. For example, residues 20-32 may be missing in the experimental
structure, so the first observed ATOM residue inside the requested range may be
33.

Therefore this module now also reports:
- requested range start/end
- observed ATOM residue start/end after filtering
- whether residues appear to be missing at the requested start/end

Design boundaries
-----------------
This module is responsible for:
- synchronizing PDB IDs between CSV and directories
- preserving an optional 'range' column from the CSV
- downloading missing PDB files
- trimming downloaded polymer ATOM records to the requested residue range
- preserving waters during trimming
- reporting requested vs observed residue bounds

This module is NOT responsible for:
- validating whether the range is biologically correct
- checking sequence consistency against UniProt
- gap detection inside the retained structure
- ligand chemistry decisions
"""

from __future__ import annotations

import csv
import traceback
import urllib.request
from datetime import datetime
from pathlib import Path

PDB_ID_COLUMN_NAME = "pdb_id"
RANGE_COLUMN_NAME = "range"
DEFAULT_PDB_ID_CSV_FILENAME = "pdb_ids.csv"
RCSB_PDB_DOWNLOAD_URL_TEMPLATE = "https://files.rcsb.org/download/{pdb_id}.pdb"

WATER_RESIDUE_NAMES = {
    "HOH",
    "WAT",
    "H2O",
    "TIP",
    "TIP3",
    "TIP3P",
    "SOL",
}

MODULE_NAME = "pdb_sync"
IGNORED_SUBDIRECTORY_NAMES = {
    "logs",
    "__pycache__",
}


# ---------------------------------------------------------------------------
# Terminal output helpers
# ---------------------------------------------------------------------------


def _screen_header(title: str, lines: list[str]) -> None:
    """
    Render a compact header block for terminal output.
    """
    all_lines = [title, *lines]
    width = max(len(line) for line in all_lines) if all_lines else len(title)

    print(f"┏{'━' * (width + 2)}┓")
    print(f"┃ {title.ljust(width)} ┃")
    if lines:
        print(f"┣{'━' * (width + 2)}┫")
        for line in lines:
            print(f"┃ {line.ljust(width)} ┃")
    print(f"┗{'━' * (width + 2)}┛")


def _screen_step(step_name: str) -> None:
    print(f"\n╔═ {step_name}")


def _screen_sub(message: str) -> None:
    print(f"╟─ {message}")


def _screen_result(message: str) -> None:
    print(f"╚═ {message}")


# ---------------------------------------------------------------------------
# Logging helpers
# ---------------------------------------------------------------------------


def _timestamp() -> str:
    return datetime.now().strftime("%Y-%m-%d %H:%M:%S")


def _get_logs_root_from_protein_data_dir(protein_data_dir: Path) -> Path:
    """
    Standard logs root for this module:

        data/proteins/logs/pdb_sync/
    """
    logs_dir = protein_data_dir / "logs" / MODULE_NAME
    logs_dir.mkdir(parents=True, exist_ok=True)
    return logs_dir


def _get_logs_root_from_pdb_file_path(pdb_file_path: Path) -> Path:
    """
    Infer the logs root from:
        data/proteins/<PDB_ID>/<PDB_ID>.pdb
    """
    protein_data_dir = pdb_file_path.parent.parent
    return _get_logs_root_from_protein_data_dir(protein_data_dir)


def _protein_log_path_from_pdb_id(
    protein_data_dir: Path,
    pdb_id: str,
) -> Path:
    """
    One log file per protein/PDB ID.
    """
    logs_root = _get_logs_root_from_protein_data_dir(protein_data_dir)
    return logs_root / f"{normalize_pdb_id(pdb_id)}.log"


def _session_log_path(protein_data_dir: Path) -> Path:
    """
    Session/global module log for the overall sync run.
    """
    logs_root = _get_logs_root_from_protein_data_dir(protein_data_dir)
    return logs_root / "_session.log"


def _append_log(log_path: Path, title: str, lines: list[str]) -> None:
    """
    Append a structured log block.
    """
    log_path.parent.mkdir(parents=True, exist_ok=True)

    with log_path.open("a", encoding="utf-8") as handle:
        handle.write("═" * 100 + "\n")
        handle.write(f"[{_timestamp()}] {title}\n")
        handle.write("─" * 100 + "\n")
        for line in lines:
            handle.write(f"{line}\n")
        handle.write("\n")


def _append_dual_log(
    session_log_path: Path,
    protein_log_path: Path | None,
    title: str,
    lines: list[str],
) -> None:
    """
    Write the same block to the session log and optionally to a per-protein log.
    """
    _append_log(session_log_path, title, lines)
    if protein_log_path is not None:
        _append_log(protein_log_path, title, lines)


def _log_exception(log_path: Path, title: str, exc: Exception) -> None:
    _append_log(
        log_path=log_path,
        title=title,
        lines=[
            f"exception_type: {type(exc).__name__}",
            f"exception_repr: {exc!r}",
            "traceback:",
            traceback.format_exc().rstrip(),
        ],
    )


# ---------------------------------------------------------------------------
# PDB/range normalization
# ---------------------------------------------------------------------------


def normalize_pdb_id(raw_pdb_id: str) -> str:
    """
    Normalize a PDB ID by stripping whitespace and converting to uppercase.
    """
    return raw_pdb_id.strip().upper()


def _is_valid_pdb_id_text(text: str) -> bool:
    """
    Return True only for valid 4-character PDB-style IDs.

    Rule:
    - exactly 4 characters
    - alphanumeric
    - first character must be a digit

    Examples
    --------
    valid:
        1ABC
        3RM2

    invalid:
        LOGS
        logs
        prepared
    """
    normalized = normalize_pdb_id(text)

    return len(normalized) == 4 and normalized.isalnum() and normalized[0].isdigit()


def normalize_range_value(raw_range_value: str | None) -> str:
    """
    Normalize a residue range value from the CSV.
    """
    if raw_range_value is None:
        return ""
    return raw_range_value.strip()


def parse_residue_range(range_value: str) -> tuple[int, int] | None:
    """
    Parse a simple inclusive residue range string such as '10-280'.
    """
    normalized_range_value = normalize_range_value(range_value)

    if not normalized_range_value:
        return None

    parts = normalized_range_value.split("-")
    if len(parts) != 2:
        raise ValueError(
            f"Invalid range value {normalized_range_value!r}. Expected format like '10-280'."
        )

    start_text, end_text = parts[0].strip(), parts[1].strip()

    try:
        start_residue_number = int(start_text)
        end_residue_number = int(end_text)
    except ValueError as exc:
        raise ValueError(
            f"Invalid range value {normalized_range_value!r}. "
            "Range boundaries must be integers."
        ) from exc

    if start_residue_number > end_residue_number:
        raise ValueError(
            f"Invalid range value {normalized_range_value!r}. "
            "Range start must be <= range end."
        )

    return (start_residue_number, end_residue_number)


# ---------------------------------------------------------------------------
# CSV handling
# ---------------------------------------------------------------------------


def read_pdb_records_from_csv(pdb_id_csv_path: Path) -> list[dict[str, str]]:
    """
    Read PDB records from a CSV file.
    """
    pdb_record_list: list[dict[str, str]] = []

    with pdb_id_csv_path.open("r", encoding="utf-8", newline="") as csv_handle:
        csv_reader = csv.DictReader(csv_handle)
        csv_fieldnames = csv_reader.fieldnames

        if csv_fieldnames is None:
            raise ValueError(
                f"CSV file {pdb_id_csv_path} is empty or missing a header row."
            )

        if PDB_ID_COLUMN_NAME not in csv_fieldnames:
            raise ValueError(
                f"CSV file {pdb_id_csv_path} must contain a '{PDB_ID_COLUMN_NAME}' "
                f"column. Found columns: {csv_fieldnames}"
            )

        csv_has_range_column = RANGE_COLUMN_NAME in csv_fieldnames

        for csv_row in csv_reader:
            raw_pdb_id = csv_row.get(PDB_ID_COLUMN_NAME, "")
            normalized_pdb_id = normalize_pdb_id(raw_pdb_id)

            if not normalized_pdb_id:
                continue

            if not _is_valid_pdb_id_text(normalized_pdb_id):
                continue

            raw_range_value = ""
            if csv_has_range_column:
                raw_range_value = csv_row.get(RANGE_COLUMN_NAME, "")

            normalized_range_value = normalize_range_value(raw_range_value)

            pdb_record = {
                PDB_ID_COLUMN_NAME: normalized_pdb_id,
                RANGE_COLUMN_NAME: normalized_range_value,
            }
            pdb_record_list.append(pdb_record)

    return pdb_record_list


def write_pdb_records_to_csv(
    pdb_record_list: list[dict[str, str]],
    pdb_id_csv_path: Path,
) -> None:
    """
    Write PDB records to a CSV file with columns:
    - pdb_id
    - range
    """
    pdb_id_to_record: dict[str, dict[str, str]] = {}

    for pdb_record in pdb_record_list:
        normalized_pdb_id = normalize_pdb_id(pdb_record.get(PDB_ID_COLUMN_NAME, ""))

        if not normalized_pdb_id:
            continue

        if not _is_valid_pdb_id_text(normalized_pdb_id):
            continue

        normalized_range_value = normalize_range_value(
            pdb_record.get(RANGE_COLUMN_NAME, "")
        )

        cleaned_record = {
            PDB_ID_COLUMN_NAME: normalized_pdb_id,
            RANGE_COLUMN_NAME: normalized_range_value,
        }

        pdb_id_to_record[normalized_pdb_id] = cleaned_record

    sorted_pdb_id_list = sorted(pdb_id_to_record.keys())

    with pdb_id_csv_path.open("w", encoding="utf-8", newline="") as csv_handle:
        csv_writer = csv.DictWriter(
            csv_handle,
            fieldnames=[PDB_ID_COLUMN_NAME, RANGE_COLUMN_NAME],
        )
        csv_writer.writeheader()

        for pdb_id in sorted_pdb_id_list:
            csv_writer.writerow(pdb_id_to_record[pdb_id])


# ---------------------------------------------------------------------------
# Directory inference helpers
# ---------------------------------------------------------------------------


def _is_valid_pdb_directory_name(directory_name: str) -> bool:
    """
    Return True only for real protein/PDB directory names.
    """
    normalized_name = normalize_pdb_id(directory_name)

    if normalized_name.lower() in {name.lower() for name in IGNORED_SUBDIRECTORY_NAMES}:
        return False

    return _is_valid_pdb_id_text(normalized_name)


def _iter_valid_protein_subdirectories(protein_data_dir: Path) -> list[Path]:
    """
    Return only valid immediate protein/PDB subdirectories inside data/proteins.
    """
    valid_subdirectory_list: list[Path] = []

    for filesystem_entry in protein_data_dir.iterdir():
        if not filesystem_entry.is_dir():
            continue

        if not _is_valid_pdb_directory_name(filesystem_entry.name):
            continue

        valid_subdirectory_list.append(filesystem_entry)

    return sorted(valid_subdirectory_list)


def get_pdb_records_from_subdirectories(protein_data_dir: Path) -> list[dict[str, str]]:
    """
    Infer PDB records from immediate subdirectory names in data/proteins.
    """
    pdb_record_list: list[dict[str, str]] = []
    seen_pdb_id_set: set[str] = set()

    for filesystem_entry in _iter_valid_protein_subdirectories(protein_data_dir):
        normalized_pdb_id = normalize_pdb_id(filesystem_entry.name)

        if normalized_pdb_id in seen_pdb_id_set:
            continue

        pdb_record = {
            PDB_ID_COLUMN_NAME: normalized_pdb_id,
            RANGE_COLUMN_NAME: "",
        }
        pdb_record_list.append(pdb_record)
        seen_pdb_id_set.add(normalized_pdb_id)

    pdb_record_list.sort(key=lambda pdb_record: pdb_record[PDB_ID_COLUMN_NAME])
    return pdb_record_list


def extract_pdb_id_list_from_pdb_record_list(
    pdb_record_list: list[dict[str, str]],
) -> list[str]:
    """
    Extract only the valid PDB IDs from a list of PDB records.
    """
    pdb_id_list: list[str] = []
    seen_pdb_id_set: set[str] = set()

    for pdb_record in pdb_record_list:
        normalized_pdb_id = normalize_pdb_id(pdb_record.get(PDB_ID_COLUMN_NAME, ""))

        if not normalized_pdb_id:
            continue

        if not _is_valid_pdb_id_text(normalized_pdb_id):
            continue

        if normalized_pdb_id in seen_pdb_id_set:
            continue

        pdb_id_list.append(normalized_pdb_id)
        seen_pdb_id_set.add(normalized_pdb_id)

    pdb_id_list.sort()
    return pdb_id_list


def create_missing_subdirectories(
    protein_data_dir: Path,
    pdb_id_list: list[str],
) -> list[Path]:
    """
    Create one subdirectory per valid PDB ID if it does not already exist.
    """
    created_subdirectory_path_list: list[Path] = []

    for pdb_id in sorted(set(pdb_id_list)):
        normalized_pdb_id = normalize_pdb_id(pdb_id)

        if not _is_valid_pdb_id_text(normalized_pdb_id):
            continue

        protein_subdirectory_path = protein_data_dir / normalized_pdb_id

        if not protein_subdirectory_path.exists():
            protein_subdirectory_path.mkdir(parents=True, exist_ok=True)
            created_subdirectory_path_list.append(protein_subdirectory_path)

    return created_subdirectory_path_list


# ---------------------------------------------------------------------------
# PDB line helpers
# ---------------------------------------------------------------------------


def _record_type(line: str) -> str:
    if line.startswith("ATOM"):
        return "ATOM"
    if line.startswith("HETATM"):
        return "HETATM"
    return ""


def _resname(line: str) -> str:
    return line[17:20].strip().upper()


def _resseq(line: str) -> int | None:
    """
    Parse the residue sequence number from a PDB ATOM/HETATM line.
    """
    raw_resseq = line[22:26].strip()

    if not raw_resseq:
        return None

    try:
        return int(raw_resseq)
    except ValueError:
        return None


def _is_water_line(line: str) -> bool:
    return _record_type(line) == "HETATM" and _resname(line) in WATER_RESIDUE_NAMES


# ---------------------------------------------------------------------------
# Range summary helpers
# ---------------------------------------------------------------------------


def get_observed_polymer_residue_bounds(
    pdb_path: Path,
) -> tuple[int | None, int | None]:
    """
    Return the minimum and maximum observed ATOM residue numbers in a PDB file.
    """
    observed_residue_numbers: set[int] = set()

    with pdb_path.open("r", encoding="utf-8") as handle:
        for line in handle:
            if _record_type(line) != "ATOM":
                continue

            residue_number = _resseq(line)
            if residue_number is None:
                continue

            observed_residue_numbers.add(residue_number)

    if not observed_residue_numbers:
        return (None, None)

    return (min(observed_residue_numbers), max(observed_residue_numbers))


def summarize_requested_vs_observed_range(
    pdb_path: Path,
    residue_range: str,
) -> dict[str, int | bool | None]:
    """
    Compare the requested CSV range against the actually observed ATOM residue bounds.
    """
    parsed_range = parse_residue_range(residue_range)

    if parsed_range is None:
        observed_start, observed_end = get_observed_polymer_residue_bounds(pdb_path)
        return {
            "requested_start": None,
            "requested_end": None,
            "observed_start": observed_start,
            "observed_end": observed_end,
            "start_missing": False,
            "end_missing": False,
        }

    requested_start, requested_end = parsed_range
    observed_start, observed_end = get_observed_polymer_residue_bounds(pdb_path)

    return {
        "requested_start": requested_start,
        "requested_end": requested_end,
        "observed_start": observed_start,
        "observed_end": observed_end,
        "start_missing": (
            observed_start is not None and observed_start > requested_start
        ),
        "end_missing": (observed_end is not None and observed_end < requested_end),
    }


# ---------------------------------------------------------------------------
# Trimming / download
# ---------------------------------------------------------------------------


def trim_pdb_to_residue_range(
    input_pdb_path: Path,
    output_pdb_path: Path,
    residue_range: str,
) -> dict[str, int | bool | None]:
    """
    Trim polymer ATOM records to the requested residue range.
    """
    protein_data_dir = output_pdb_path.parent.parent
    session_log_path = _session_log_path(protein_data_dir)
    protein_log_path = _protein_log_path_from_pdb_id(
        protein_data_dir,
        output_pdb_path.stem,
    )

    parsed_range = parse_residue_range(residue_range)

    if parsed_range is None:
        output_pdb_path.write_text(
            input_pdb_path.read_text(encoding="utf-8"),
            encoding="utf-8",
        )

        range_summary = summarize_requested_vs_observed_range(
            output_pdb_path,
            residue_range,
        )

        _append_dual_log(
            session_log_path=session_log_path,
            protein_log_path=protein_log_path,
            title="trim_pdb_to_residue_range",
            lines=[
                f"input_pdb_path  : {input_pdb_path}",
                f"output_pdb_path : {output_pdb_path}",
                f"residue_range   : {residue_range!r}",
                "mode            : copy_unchanged",
                f"requested_start : {range_summary['requested_start']}",
                f"requested_end   : {range_summary['requested_end']}",
                f"observed_start  : {range_summary['observed_start']}",
                f"observed_end    : {range_summary['observed_end']}",
                f"start_missing   : {range_summary['start_missing']}",
                f"end_missing     : {range_summary['end_missing']}",
            ],
        )
        return range_summary

    start_residue_number, end_residue_number = parsed_range

    kept_lines: list[str] = []
    total_atom_lines = 0
    kept_atom_lines = 0
    total_water_hetatm_lines = 0
    kept_water_hetatm_lines = 0
    total_other_hetatm_lines = 0
    kept_other_hetatm_lines = 0

    with input_pdb_path.open("r", encoding="utf-8") as handle:
        for line in handle:
            record_type = _record_type(line)

            if record_type == "ATOM":
                total_atom_lines += 1
                residue_number = _resseq(line)
                if residue_number is None:
                    continue

                if start_residue_number <= residue_number <= end_residue_number:
                    kept_lines.append(line)
                    kept_atom_lines += 1
                continue

            if _is_water_line(line):
                total_water_hetatm_lines += 1
                kept_lines.append(line)
                kept_water_hetatm_lines += 1
                continue

            if record_type == "HETATM":
                total_other_hetatm_lines += 1
                kept_other_hetatm_lines += 1

            kept_lines.append(line)

    with output_pdb_path.open("w", encoding="utf-8") as handle:
        handle.writelines(kept_lines)

    range_summary = summarize_requested_vs_observed_range(
        output_pdb_path,
        residue_range,
    )

    _append_dual_log(
        session_log_path=session_log_path,
        protein_log_path=protein_log_path,
        title="trim_pdb_to_residue_range",
        lines=[
            f"input_pdb_path        : {input_pdb_path}",
            f"output_pdb_path       : {output_pdb_path}",
            f"residue_range         : {residue_range!r}",
            f"parsed_start          : {start_residue_number}",
            f"parsed_end            : {end_residue_number}",
            f"total_atom_lines      : {total_atom_lines}",
            f"kept_atom_lines       : {kept_atom_lines}",
            f"total_water_hetatm    : {total_water_hetatm_lines}",
            f"kept_water_hetatm     : {kept_water_hetatm_lines}",
            f"total_other_hetatm    : {total_other_hetatm_lines}",
            f"kept_other_hetatm     : {kept_other_hetatm_lines}",
            f"requested_start       : {range_summary['requested_start']}",
            f"requested_end         : {range_summary['requested_end']}",
            f"observed_start        : {range_summary['observed_start']}",
            f"observed_end          : {range_summary['observed_end']}",
            f"start_missing         : {range_summary['start_missing']}",
            f"end_missing           : {range_summary['end_missing']}",
        ],
    )

    return range_summary


def download_raw_pdb_file(pdb_id: str, target_pdb_file_path: Path) -> None:
    """
    Download a raw PDB file from the RCSB archive.
    """
    normalized_pdb_id = normalize_pdb_id(pdb_id)
    rcsb_pdb_download_url = RCSB_PDB_DOWNLOAD_URL_TEMPLATE.format(
        pdb_id=normalized_pdb_id
    )

    target_pdb_file_path.parent.mkdir(parents=True, exist_ok=True)
    urllib.request.urlretrieve(rcsb_pdb_download_url, target_pdb_file_path)


def download_and_prepare_pdb_file(
    pdb_id: str,
    residue_range: str,
    target_pdb_file_path: Path,
) -> dict[str, int | bool | None]:
    """
    Download a PDB file and apply optional range trimming.
    """
    normalized_pdb_id = normalize_pdb_id(pdb_id)
    temporary_raw_pdb_file_path = target_pdb_file_path.with_suffix(".raw.pdb")

    protein_data_dir = target_pdb_file_path.parent.parent
    session_log_path = _session_log_path(protein_data_dir)
    protein_log_path = _protein_log_path_from_pdb_id(
        protein_data_dir,
        normalized_pdb_id,
    )

    try:
        _append_dual_log(
            session_log_path=session_log_path,
            protein_log_path=protein_log_path,
            title="download_and_prepare_pdb_file:start",
            lines=[
                f"pdb_id                 : {normalized_pdb_id}",
                f"residue_range          : {residue_range!r}",
                f"temporary_raw_pdb_path : {temporary_raw_pdb_file_path}",
                f"target_pdb_file_path   : {target_pdb_file_path}",
                f"download_url           : {RCSB_PDB_DOWNLOAD_URL_TEMPLATE.format(pdb_id=normalized_pdb_id)}",
            ],
        )

        download_raw_pdb_file(normalized_pdb_id, temporary_raw_pdb_file_path)

        _append_dual_log(
            session_log_path=session_log_path,
            protein_log_path=protein_log_path,
            title="download_and_prepare_pdb_file:downloaded_raw",
            lines=[
                f"temporary_raw_pdb_path : {temporary_raw_pdb_file_path}",
                f"exists                 : {temporary_raw_pdb_file_path.exists()}",
                f"size_bytes             : {temporary_raw_pdb_file_path.stat().st_size if temporary_raw_pdb_file_path.exists() else 'n/a'}",
            ],
        )

        range_summary = trim_pdb_to_residue_range(
            input_pdb_path=temporary_raw_pdb_file_path,
            output_pdb_path=target_pdb_file_path,
            residue_range=residue_range,
        )

        _append_dual_log(
            session_log_path=session_log_path,
            protein_log_path=protein_log_path,
            title="download_and_prepare_pdb_file:done",
            lines=[
                f"pdb_id          : {normalized_pdb_id}",
                f"residue_range   : {residue_range!r}",
                f"requested_start : {range_summary['requested_start']}",
                f"requested_end   : {range_summary['requested_end']}",
                f"observed_start  : {range_summary['observed_start']}",
                f"observed_end    : {range_summary['observed_end']}",
                f"start_missing   : {range_summary['start_missing']}",
                f"end_missing     : {range_summary['end_missing']}",
            ],
        )

        return range_summary

    except Exception as exc:
        _log_exception(session_log_path, "download_and_prepare_pdb_file:exception", exc)
        _log_exception(protein_log_path, "download_and_prepare_pdb_file:exception", exc)
        raise
    finally:
        if temporary_raw_pdb_file_path.exists():
            temporary_raw_pdb_file_path.unlink()

            _append_dual_log(
                session_log_path=session_log_path,
                protein_log_path=protein_log_path,
                title="download_and_prepare_pdb_file:cleanup",
                lines=[
                    f"removed_temporary_raw_pdb_path : {temporary_raw_pdb_file_path}",
                ],
            )


def download_missing_pdb_files(
    protein_data_dir: Path,
    pdb_record_list: list[dict[str, str]],
) -> list[Path]:
    """
    Download missing PDB files for a list of PDB records.
    """
    downloaded_pdb_file_path_list: list[Path] = []
    session_log_path = _session_log_path(protein_data_dir)

    record_by_pdb_id: dict[str, dict[str, str]] = {}

    for pdb_record in pdb_record_list:
        normalized_pdb_id = normalize_pdb_id(pdb_record.get(PDB_ID_COLUMN_NAME, ""))
        if not normalized_pdb_id:
            continue

        if not _is_valid_pdb_id_text(normalized_pdb_id):
            continue

        record_by_pdb_id[normalized_pdb_id] = {
            PDB_ID_COLUMN_NAME: normalized_pdb_id,
            RANGE_COLUMN_NAME: normalize_range_value(
                pdb_record.get(RANGE_COLUMN_NAME, "")
            ),
        }

    _append_log(
        session_log_path,
        "download_missing_pdb_files:start",
        [
            f"record_count      : {len(record_by_pdb_id)}",
            f"protein_data_dir  : {protein_data_dir}",
        ],
    )

    for pdb_id in sorted(record_by_pdb_id):
        if not _is_valid_pdb_id_text(pdb_id):
            continue

        pdb_record = record_by_pdb_id[pdb_id]
        residue_range = pdb_record[RANGE_COLUMN_NAME]

        protein_subdirectory_path = protein_data_dir / pdb_id
        pdb_file_path = protein_subdirectory_path / f"{pdb_id}.pdb"
        protein_log_path = _protein_log_path_from_pdb_id(protein_data_dir, pdb_id)

        protein_subdirectory_path.mkdir(parents=True, exist_ok=True)

        if not pdb_file_path.exists():
            _screen_sub(f"{pdb_id} | download")
            try:
                range_summary = download_and_prepare_pdb_file(
                    pdb_id=pdb_id,
                    residue_range=residue_range,
                    target_pdb_file_path=pdb_file_path,
                )
            except Exception as exc:
                _append_dual_log(
                    session_log_path=session_log_path,
                    protein_log_path=protein_log_path,
                    title="download_missing_pdb_files:failed",
                    lines=[
                        f"pdb_id        : {pdb_id}",
                        f"residue_range : {residue_range!r}",
                        f"pdb_file_path : {pdb_file_path}",
                        f"error         : {exc!r}",
                    ],
                )
                raise

            _append_dual_log(
                session_log_path=session_log_path,
                protein_log_path=protein_log_path,
                title="download_missing_pdb_files:summary",
                lines=[
                    f"pdb_id          : {pdb_id}",
                    f"residue_range   : {residue_range!r}",
                    f"pdb_file_path   : {pdb_file_path}",
                    f"requested_start : {range_summary['requested_start']}",
                    f"requested_end   : {range_summary['requested_end']}",
                    f"observed_start  : {range_summary['observed_start']}",
                    f"observed_end    : {range_summary['observed_end']}",
                    f"start_missing   : {range_summary['start_missing']}",
                    f"end_missing     : {range_summary['end_missing']}",
                ],
            )

            downloaded_pdb_file_path_list.append(pdb_file_path)
        else:
            _screen_sub(f"{pdb_id} | already present")
            _append_dual_log(
                session_log_path=session_log_path,
                protein_log_path=protein_log_path,
                title="download_missing_pdb_files:skip_existing",
                lines=[
                    f"pdb_id        : {pdb_id}",
                    f"residue_range : {residue_range!r}",
                    f"pdb_file_path : {pdb_file_path}",
                ],
            )

    _append_log(
        session_log_path,
        "download_missing_pdb_files:done",
        [
            f"downloaded_count : {len(downloaded_pdb_file_path_list)}",
        ],
    )

    return downloaded_pdb_file_path_list


# ---------------------------------------------------------------------------
# Merge helpers
# ---------------------------------------------------------------------------


def merge_csv_records_with_directory_records(
    csv_pdb_record_list: list[dict[str, str]],
    directory_pdb_record_list: list[dict[str, str]],
) -> list[dict[str, str]]:
    """
    Merge records from CSV and records inferred from directories.
    """
    merged_pdb_id_to_record: dict[str, dict[str, str]] = {}

    for directory_pdb_record in directory_pdb_record_list:
        directory_pdb_id = normalize_pdb_id(
            directory_pdb_record.get(PDB_ID_COLUMN_NAME, "")
        )

        if not directory_pdb_id:
            continue

        if not _is_valid_pdb_id_text(directory_pdb_id):
            continue

        merged_pdb_id_to_record[directory_pdb_id] = {
            PDB_ID_COLUMN_NAME: directory_pdb_id,
            RANGE_COLUMN_NAME: "",
        }

    for csv_pdb_record in csv_pdb_record_list:
        csv_pdb_id = normalize_pdb_id(csv_pdb_record.get(PDB_ID_COLUMN_NAME, ""))

        if not csv_pdb_id:
            continue

        if not _is_valid_pdb_id_text(csv_pdb_id):
            continue

        csv_range_value = normalize_range_value(
            csv_pdb_record.get(RANGE_COLUMN_NAME, "")
        )

        merged_pdb_id_to_record[csv_pdb_id] = {
            PDB_ID_COLUMN_NAME: csv_pdb_id,
            RANGE_COLUMN_NAME: csv_range_value,
        }

    merged_pdb_record_list = list(merged_pdb_id_to_record.values())
    merged_pdb_record_list.sort(key=lambda pdb_record: pdb_record[PDB_ID_COLUMN_NAME])

    return merged_pdb_record_list


# ---------------------------------------------------------------------------
# Main sync entry point
# ---------------------------------------------------------------------------


def sync_pdb_csv_and_directories(
    protein_data_dir: Path,
    pdb_id_csv_filename: str = DEFAULT_PDB_ID_CSV_FILENAME,
) -> None:
    """
    Synchronize the PDB CSV file and the protein subdirectory structure.

    Cases
    -----
    Case 1:
        No CSV exists, but protein subdirectories exist
        -> create CSV from subdirectory names
        -> range values will be empty

    Case 2:
        CSV exists, but no protein subdirectories exist
        -> create subdirectories and download range-aware PDB files

    Case 3:
        Both CSV and protein subdirectories exist
        -> merge both sides
        -> preserve CSV range values
        -> add missing IDs from either side
        -> download missing PDB files
    """
    protein_data_dir.mkdir(parents=True, exist_ok=True)
    pdb_id_csv_path = protein_data_dir / pdb_id_csv_filename

    valid_protein_subdirectory_list = _iter_valid_protein_subdirectories(
        protein_data_dir
    )
    any_protein_subdirectory_exists = len(valid_protein_subdirectory_list) > 0

    logs_root = _get_logs_root_from_protein_data_dir(protein_data_dir)
    session_log_path = _session_log_path(protein_data_dir)

    csv_file_exists = pdb_id_csv_path.exists()

    _screen_header(
        "PDB SYNC",
        [
            f"data_dir   : {protein_data_dir}",
            f"csv_path   : {pdb_id_csv_path}",
            f"logs_dir   : {logs_root}",
            f"csv_exists : {csv_file_exists}",
            f"has_dirs   : {any_protein_subdirectory_exists}",
        ],
    )

    _append_log(
        session_log_path,
        "sync_pdb_csv_and_directories:start",
        [
            f"protein_data_dir            : {protein_data_dir}",
            f"pdb_id_csv_path             : {pdb_id_csv_path}",
            f"logs_root                   : {logs_root}",
            f"csv_file_exists             : {csv_file_exists}",
            f"any_protein_subdirs_exist   : {any_protein_subdirectory_exists}",
            f"valid_protein_subdir_names  : {[path.name for path in valid_protein_subdirectory_list]}",
        ],
    )

    if not csv_file_exists and any_protein_subdirectory_exists:
        _screen_step("case 1 | no csv / dirs exist")

        directory_pdb_record_list = get_pdb_records_from_subdirectories(
            protein_data_dir
        )
        write_pdb_records_to_csv(directory_pdb_record_list, pdb_id_csv_path)

        _append_log(
            session_log_path,
            "sync_pdb_csv_and_directories:case1",
            [
                f"directory_record_count : {len(directory_pdb_record_list)}",
                f"csv_written            : {pdb_id_csv_path}",
            ],
        )

        _screen_result(f"created CSV from {len(directory_pdb_record_list)} directories")
        return

    if csv_file_exists and not any_protein_subdirectory_exists:
        _screen_step("case 2 | csv exists / no dirs")

        csv_pdb_record_list = read_pdb_records_from_csv(pdb_id_csv_path)
        csv_pdb_id_list = extract_pdb_id_list_from_pdb_record_list(csv_pdb_record_list)

        created_dirs = create_missing_subdirectories(protein_data_dir, csv_pdb_id_list)
        downloaded_files = download_missing_pdb_files(
            protein_data_dir,
            csv_pdb_record_list,
        )

        _append_log(
            session_log_path,
            "sync_pdb_csv_and_directories:case2",
            [
                f"csv_record_count      : {len(csv_pdb_record_list)}",
                f"created_dir_count     : {len(created_dirs)}",
                f"downloaded_file_count : {len(downloaded_files)}",
            ],
        )

        _screen_sub(f"created {len(created_dirs)} directories")
        _screen_result("download phase completed")
        return

    if not csv_file_exists and not any_protein_subdirectory_exists:
        _screen_step("case 0 | nothing present")
        _append_log(
            session_log_path,
            "sync_pdb_csv_and_directories:case0",
            ["no CSV and no protein subdirectories"],
        )
        _screen_result("no CSV and no protein subdirectories")
        return

    _screen_step("case 3 | csv and dirs both exist")

    csv_pdb_record_list = read_pdb_records_from_csv(pdb_id_csv_path)
    directory_pdb_record_list = get_pdb_records_from_subdirectories(protein_data_dir)

    csv_pdb_id_set = set(extract_pdb_id_list_from_pdb_record_list(csv_pdb_record_list))
    directory_pdb_id_set = set(
        extract_pdb_id_list_from_pdb_record_list(directory_pdb_record_list)
    )

    missing_subdirectory_pdb_id_list = sorted(csv_pdb_id_set - directory_pdb_id_set)
    missing_csv_pdb_id_list = sorted(directory_pdb_id_set - csv_pdb_id_set)

    if missing_subdirectory_pdb_id_list:
        create_missing_subdirectories(
            protein_data_dir,
            missing_subdirectory_pdb_id_list,
        )
        _screen_sub(
            f"dirs created from CSV only: {', '.join(missing_subdirectory_pdb_id_list)}"
        )

    if missing_csv_pdb_id_list:
        _screen_sub(f"dirs only / added to CSV: {', '.join(missing_csv_pdb_id_list)}")

    merged_pdb_record_list = merge_csv_records_with_directory_records(
        csv_pdb_record_list,
        directory_pdb_record_list,
    )

    write_pdb_records_to_csv(merged_pdb_record_list, pdb_id_csv_path)
    _screen_sub(f"CSV synchronized ({len(merged_pdb_record_list)} records)")

    downloaded_files = download_missing_pdb_files(
        protein_data_dir,
        merged_pdb_record_list,
    )

    _append_log(
        session_log_path,
        "sync_pdb_csv_and_directories:case3",
        [
            f"csv_record_count             : {len(csv_pdb_record_list)}",
            f"directory_record_count       : {len(directory_pdb_record_list)}",
            f"missing_subdirectory_pdb_ids : {missing_subdirectory_pdb_id_list}",
            f"missing_csv_pdb_ids          : {missing_csv_pdb_id_list}",
            f"merged_record_count          : {len(merged_pdb_record_list)}",
            f"downloaded_file_count        : {len(downloaded_files)}",
        ],
    )

    _screen_result("download phase completed")
