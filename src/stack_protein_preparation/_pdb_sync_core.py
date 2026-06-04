"""Core helpers for pdb_sync: constants, logging, parsing, CSV I/O, PDB line utilities."""
from __future__ import annotations

import csv
import re
import traceback
import urllib.request
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path

from stack_protein_preparation.pdb_components import WATER_NAMES as WATER_RESIDUE_NAMES

PDB_ID_COLUMN_NAME = "pdb_id"
RANGE_COLUMN_NAME = "range"
DEFAULT_PDB_ID_CSV_FILENAME = "pdb_ids.csv"
RCSB_PDB_DOWNLOAD_URL_TEMPLATE = "https://files.rcsb.org/download/{pdb_id}.pdb"

MODULE_NAME = "pdb_sync"
IGNORED_SUBDIRECTORY_NAMES = {
    "logs",
    "__pycache__",
}


# ---------------------------------------------------------------------------
# Terminal output helpers
# ---------------------------------------------------------------------------


def _screen_header(title: str, lines: list[str]) -> None:
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
    logs_dir = protein_data_dir / "logs" / MODULE_NAME
    logs_dir.mkdir(parents=True, exist_ok=True)
    return logs_dir


def _get_logs_root_from_pdb_file_path(pdb_file_path: Path) -> Path:
    protein_data_dir = pdb_file_path.parent.parent
    return _get_logs_root_from_protein_data_dir(protein_data_dir)


def _protein_log_path_from_pdb_id(
    protein_data_dir: Path,
    pdb_id: str,
) -> Path:
    logs_root = _get_logs_root_from_protein_data_dir(protein_data_dir)
    return logs_root / f"{normalize_pdb_id(pdb_id)}.log"


def _session_log_path(protein_data_dir: Path) -> Path:
    logs_root = _get_logs_root_from_protein_data_dir(protein_data_dir)
    return logs_root / "_session.log"


def _append_log(log_path: Path, title: str, lines: list[str]) -> None:
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
    """Normalize a PDB ID by stripping whitespace and converting to uppercase."""
    return raw_pdb_id.strip().upper()


def _is_valid_pdb_id_text(text: str) -> bool:
    """Return True only for valid 4-character PDB-style IDs (digit-first, alphanumeric)."""
    normalized = normalize_pdb_id(text)
    return len(normalized) == 4 and normalized.isalnum() and normalized[0].isdigit()


def normalize_range_value(raw_range_value: str | None) -> str:
    """Normalize a residue range value from the CSV."""
    if raw_range_value is None:
        return ""
    return raw_range_value.strip()


@dataclass(frozen=True)
class ParsedRange:
    """Structured representation of a residue range, optionally chain-aware.

    For a plain integer range ('25-200'): start_chain and end_chain are None.
    For a chain-aware range ('A25-B15'): both chain fields are set.
    """

    start_resnum: int
    end_resnum: int
    start_chain: str | None = None
    end_chain: str | None = None

    @property
    def is_chain_aware(self) -> bool:
        return self.start_chain is not None

    @property
    def is_cross_chain(self) -> bool:
        return (
            self.start_chain is not None
            and self.end_chain is not None
            and self.start_chain != self.end_chain
        )


_CHAIN_AWARE_RANGE_RE = re.compile(r"([A-Za-z])(-?\d+)-([A-Za-z])(-?\d+)")
_LEGACY_RANGE_RE = re.compile(r"(-?\d+)-(-?\d+)")


def _parse_range_internal(normalized: str) -> ParsedRange:
    """Parse a non-empty, stripped range string into a ParsedRange."""
    m = _CHAIN_AWARE_RANGE_RE.fullmatch(normalized)
    if m:
        start_chain = m.group(1).upper()
        start_resnum = int(m.group(2))
        end_chain = m.group(3).upper()
        end_resnum = int(m.group(4))
        return ParsedRange(
            start_resnum=start_resnum,
            end_resnum=end_resnum,
            start_chain=start_chain,
            end_chain=end_chain,
        )
    m = _LEGACY_RANGE_RE.fullmatch(normalized)
    if m:
        start_resnum = int(m.group(1))
        end_resnum = int(m.group(2))
        if start_resnum > end_resnum:
            raise ValueError(
                f"Invalid range value {normalized!r}. Range start must be <= range end."
            )
        return ParsedRange(start_resnum=start_resnum, end_resnum=end_resnum)
    raise ValueError(
        f"Invalid range value {normalized!r}. "
        "Expected '10-280' (legacy) or 'A10-B280' (chain-aware)."
    )


def parse_residue_range(range_value: str) -> tuple[int, int] | None:
    """Parse a residue range string; return (start, end) integers or None.

    Accepts both the legacy '10-280' format and the new chain-aware 'A10-B280'
    format.  For chain-aware ranges the chain identifiers are discarded here —
    use parse_residue_range_chain_aware when chain information is needed.
    """
    normalized = normalize_range_value(range_value)
    if not normalized:
        return None
    parsed = _parse_range_internal(normalized)
    return (parsed.start_resnum, parsed.end_resnum)


def parse_residue_range_chain_aware(range_value: str) -> ParsedRange | None:
    """Parse a residue range string and return a full ParsedRange, or None if empty."""
    normalized = normalize_range_value(range_value)
    if not normalized:
        return None
    return _parse_range_internal(normalized)


# ---------------------------------------------------------------------------
# CSV handling
# ---------------------------------------------------------------------------


def read_pdb_records_from_csv(pdb_id_csv_path: Path) -> list[dict[str, str]]:
    """Read PDB records from a CSV file."""
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
    """Write PDB records to a CSV file with columns: pdb_id, range."""
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
    """Return True only for real protein/PDB directory names."""
    normalized_name = normalize_pdb_id(directory_name)

    if normalized_name.lower() in {name.lower() for name in IGNORED_SUBDIRECTORY_NAMES}:
        return False

    return _is_valid_pdb_id_text(normalized_name)


def _iter_valid_protein_subdirectories(protein_data_dir: Path) -> list[Path]:
    """Return only valid immediate protein/PDB subdirectories inside data/proteins."""
    valid_subdirectory_list: list[Path] = []

    for filesystem_entry in protein_data_dir.iterdir():
        if not filesystem_entry.is_dir():
            continue

        if not _is_valid_pdb_directory_name(filesystem_entry.name):
            continue

        valid_subdirectory_list.append(filesystem_entry)

    return sorted(valid_subdirectory_list)


def get_pdb_records_from_subdirectories(protein_data_dir: Path) -> list[dict[str, str]]:
    """Infer PDB records from immediate subdirectory names in data/proteins."""
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
    """Extract only the valid PDB IDs from a list of PDB records."""
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
    """Create one subdirectory per valid PDB ID if it does not already exist."""
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
    """Parse the residue sequence number from a PDB ATOM/HETATM line."""
    raw_resseq = line[22:26].strip()

    if not raw_resseq:
        return None

    try:
        return int(raw_resseq)
    except ValueError:
        return None


def _chain_id(line: str) -> str:
    """Return the chain identifier from a PDB ATOM/HETATM line (column 22, 0-indexed 21)."""
    if len(line) > 21:
        return line[21]
    return ""


def _is_water_line(line: str) -> bool:
    return _record_type(line) == "HETATM" and _resname(line) in WATER_RESIDUE_NAMES


def _atom_line_passes_range(line: str, residue_number: int, parsed_range: ParsedRange) -> bool:
    """Return True if this ATOM line should be kept given the parsed range."""
    if not parsed_range.is_chain_aware:
        return parsed_range.start_resnum <= residue_number <= parsed_range.end_resnum

    chain = _chain_id(line)
    start_chain = parsed_range.start_chain
    end_chain = parsed_range.end_chain

    if chain == start_chain and chain == end_chain:
        return parsed_range.start_resnum <= residue_number <= parsed_range.end_resnum
    if chain == start_chain:
        return residue_number >= parsed_range.start_resnum
    if chain == end_chain:
        return residue_number <= parsed_range.end_resnum
    return False


def get_observed_polymer_residue_bounds(
    pdb_path: Path,
) -> tuple[int | None, int | None]:
    """Return the minimum and maximum observed ATOM residue numbers in a PDB file."""
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
    """Compare the requested CSV range against the actually observed ATOM residue bounds."""
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


def trim_pdb_to_residue_range(
    input_pdb_path: Path,
    output_pdb_path: Path,
    residue_range: str,
) -> dict[str, int | bool | None]:
    """Trim polymer ATOM records to the requested residue range."""
    protein_data_dir = output_pdb_path.parent.parent
    session_log_path = _session_log_path(protein_data_dir)
    protein_log_path = _protein_log_path_from_pdb_id(
        protein_data_dir,
        output_pdb_path.stem,
    )

    parsed_range = parse_residue_range_chain_aware(residue_range)

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

                if _atom_line_passes_range(line, residue_number, parsed_range):
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
            f"parsed_start          : {parsed_range.start_resnum}",
            f"parsed_end            : {parsed_range.end_resnum}",
            f"start_chain           : {parsed_range.start_chain}",
            f"end_chain             : {parsed_range.end_chain}",
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
