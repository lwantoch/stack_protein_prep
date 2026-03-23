"""
/home/grheco/repositorios/stack_protein_prep/src/stack_protein_preparation/pdb_sync.py

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
import shutil
import urllib.request
from pathlib import Path
from typing import Any

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


def normalize_pdb_id(raw_pdb_id: str) -> str:
    """
    Normalize a PDB ID by stripping whitespace and converting to uppercase.
    """
    return raw_pdb_id.strip().upper()


def normalize_range_value(raw_range_value: str | None) -> str:
    """
    Normalize a residue range value from the CSV.

    Rules
    -----
    - If the value is None, return an empty string.
    - Otherwise strip surrounding whitespace.
    - No further validation is performed here.
    """
    if raw_range_value is None:
        return ""

    return raw_range_value.strip()


def parse_residue_range(range_value: str) -> tuple[int, int] | None:
    """
    Parse a simple inclusive residue range string such as '10-280'.

    Parameters
    ----------
    range_value
        Range string from the CSV.

    Returns
    -------
    tuple[int, int] | None
        (start_residue_number, end_residue_number), or None if the input is empty.

    Raises
    ------
    ValueError
        If the range string is non-empty but not in the expected simple form.
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


def read_pdb_records_from_csv(pdb_id_csv_path: Path) -> list[dict[str, str]]:
    """
    Read PDB records from a CSV file.

    Expected behavior
    -----------------
    - 'pdb_id' column is required
    - 'range' column is optional
    - empty PDB ID rows are ignored
    - PDB IDs are normalized to uppercase
    - range values are preserved as plain text

    Returns
    -------
    list[dict[str, str]]
        Example:
        [
            {"pdb_id": "1ABC", "range": "10-280"},
            {"pdb_id": "2XYZ", "range": ""},
        ]
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

    Behavior
    --------
    - keeps exactly one row per PDB ID
    - if the same PDB ID appears multiple times, the last occurrence wins
    - output is sorted by PDB ID
    - range is written as empty string if missing
    """
    pdb_id_to_record: dict[str, dict[str, str]] = {}

    for pdb_record in pdb_record_list:
        normalized_pdb_id = normalize_pdb_id(pdb_record.get(PDB_ID_COLUMN_NAME, ""))

        if not normalized_pdb_id:
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


def get_pdb_records_from_subdirectories(protein_data_dir: Path) -> list[dict[str, str]]:
    """
    Infer PDB records from subdirectory names.

    Important
    ---------
    The filesystem only gives us the PDB ID.
    It does NOT tell us the residue range.

    Therefore, every record created from directories gets:
        range = ""

    Returns
    -------
    list[dict[str, str]]
        Example:
        [
            {"pdb_id": "1ABC", "range": ""},
            {"pdb_id": "2XYZ", "range": ""},
        ]
    """
    pdb_record_list: list[dict[str, str]] = []
    seen_pdb_id_set: set[str] = set()

    for filesystem_entry in protein_data_dir.iterdir():
        if not filesystem_entry.is_dir():
            continue

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
    Extract only the PDB IDs from a list of PDB records.

    This helper is useful because directory creation only needs the PDB ID.
    """
    pdb_id_list: list[str] = []
    seen_pdb_id_set: set[str] = set()

    for pdb_record in pdb_record_list:
        normalized_pdb_id = normalize_pdb_id(pdb_record.get(PDB_ID_COLUMN_NAME, ""))

        if not normalized_pdb_id:
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
    Create one subdirectory per PDB ID if it does not already exist.
    """
    created_subdirectory_path_list: list[Path] = []

    for pdb_id in sorted(set(pdb_id_list)):
        protein_subdirectory_path = protein_data_dir / pdb_id

        if not protein_subdirectory_path.exists():
            protein_subdirectory_path.mkdir(parents=True, exist_ok=True)
            created_subdirectory_path_list.append(protein_subdirectory_path)

    return created_subdirectory_path_list


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

    Returns None if parsing fails.
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


def get_observed_polymer_residue_bounds(
    pdb_path: Path,
) -> tuple[int | None, int | None]:
    """
    Return the minimum and maximum observed ATOM residue numbers in a PDB file.

    Notes
    -----
    - only ATOM records are considered
    - this reports observed coordinates, not requested or biological sequence bounds
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

    Returns
    -------
    dict[str, int | bool | None]
        Keys:
        - requested_start
        - requested_end
        - observed_start
        - observed_end
        - start_missing
        - end_missing
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


def trim_pdb_to_residue_range(
    input_pdb_path: Path,
    output_pdb_path: Path,
    residue_range: str,
) -> dict[str, int | bool | None]:
    """
    Trim polymer ATOM records to the requested residue range.

    Behavior
    --------
    - empty range: copy file unchanged
    - ATOM lines: keep only residues whose residue number is within the range
    - HETATM water lines: always keep
    - other HETATM lines: keep unchanged
    - other record types: keep unchanged

    Returns
    -------
    dict[str, int | bool | None]
        Summary of requested vs observed bounds for the written output file.
    """
    parsed_range = parse_residue_range(residue_range)

    if parsed_range is None:
        shutil.copyfile(input_pdb_path, output_pdb_path)
        return summarize_requested_vs_observed_range(output_pdb_path, residue_range)

    start_residue_number, end_residue_number = parsed_range

    kept_lines: list[str] = []

    with input_pdb_path.open("r", encoding="utf-8") as handle:
        for line in handle:
            record_type = _record_type(line)

            if record_type == "ATOM":
                residue_number = _resseq(line)
                if residue_number is None:
                    continue

                if start_residue_number <= residue_number <= end_residue_number:
                    kept_lines.append(line)
                continue

            if _is_water_line(line):
                kept_lines.append(line)
                continue

            kept_lines.append(line)

    with output_pdb_path.open("w", encoding="utf-8") as handle:
        handle.writelines(kept_lines)

    return summarize_requested_vs_observed_range(output_pdb_path, residue_range)


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

    Implementation detail
    ---------------------
    The raw download is first written to a temporary file in the same directory,
    then trimmed into the final target path.

    Returns
    -------
    dict[str, int | bool | None]
        Requested vs observed residue-bound summary.
    """
    normalized_pdb_id = normalize_pdb_id(pdb_id)
    temporary_raw_pdb_file_path = target_pdb_file_path.with_suffix(".raw.pdb")

    try:
        download_raw_pdb_file(normalized_pdb_id, temporary_raw_pdb_file_path)
        return trim_pdb_to_residue_range(
            input_pdb_path=temporary_raw_pdb_file_path,
            output_pdb_path=target_pdb_file_path,
            residue_range=residue_range,
        )
    finally:
        if temporary_raw_pdb_file_path.exists():
            temporary_raw_pdb_file_path.unlink()


def download_missing_pdb_files(
    protein_data_dir: Path,
    pdb_record_list: list[dict[str, str]],
) -> list[Path]:
    """
    Download missing PDB files for a list of PDB records.

    Range-aware behavior
    --------------------
    - the local PDB file is created from the CSV record
    - if a range is present, the downloaded file is trimmed accordingly
    - existing local PDB files are left untouched
    """
    downloaded_pdb_file_path_list: list[Path] = []

    record_by_pdb_id: dict[str, dict[str, str]] = {}

    for pdb_record in pdb_record_list:
        normalized_pdb_id = normalize_pdb_id(pdb_record.get(PDB_ID_COLUMN_NAME, ""))
        if not normalized_pdb_id:
            continue

        record_by_pdb_id[normalized_pdb_id] = {
            PDB_ID_COLUMN_NAME: normalized_pdb_id,
            RANGE_COLUMN_NAME: normalize_range_value(
                pdb_record.get(RANGE_COLUMN_NAME, "")
            ),
        }

    for pdb_id in sorted(record_by_pdb_id):
        pdb_record = record_by_pdb_id[pdb_id]
        residue_range = pdb_record[RANGE_COLUMN_NAME]

        protein_subdirectory_path = protein_data_dir / pdb_id
        pdb_file_path = protein_subdirectory_path / f"{pdb_id}.pdb"

        protein_subdirectory_path.mkdir(parents=True, exist_ok=True)

        if not pdb_file_path.exists():
            print(
                f"[INFO] Downloading PDB for {pdb_id} "
                f"(range={residue_range!r}) -> {pdb_file_path}"
            )
            range_summary = download_and_prepare_pdb_file(
                pdb_id=pdb_id,
                residue_range=residue_range,
                target_pdb_file_path=pdb_file_path,
            )

            print(
                f"[INFO] Prepared PDB for {pdb_id}: "
                f"requested=({range_summary['requested_start']}, "
                f"{range_summary['requested_end']}), "
                f"observed=({range_summary['observed_start']}, "
                f"{range_summary['observed_end']}), "
                f"start_missing={range_summary['start_missing']}, "
                f"end_missing={range_summary['end_missing']}"
            )

            downloaded_pdb_file_path_list.append(pdb_file_path)
        else:
            print(f"[INFO] PDB already exists for {pdb_id}: {pdb_file_path}")

    return downloaded_pdb_file_path_list


def merge_csv_records_with_directory_records(
    csv_pdb_record_list: list[dict[str, str]],
    directory_pdb_record_list: list[dict[str, str]],
) -> list[dict[str, str]]:
    """
    Merge records from CSV and records inferred from directories.

    Rule
    ----
    CSV data has priority for the 'range' field, because only the CSV can store
    that metadata.

    Practical meaning
    -----------------
    - if a PDB ID exists in CSV, keep its CSV range
    - if a PDB ID exists only as a directory, create a record with empty range
    """
    merged_pdb_id_to_record: dict[str, dict[str, str]] = {}

    for directory_pdb_record in directory_pdb_record_list:
        directory_pdb_id = normalize_pdb_id(
            directory_pdb_record.get(PDB_ID_COLUMN_NAME, "")
        )

        if not directory_pdb_id:
            continue

        merged_pdb_id_to_record[directory_pdb_id] = {
            PDB_ID_COLUMN_NAME: directory_pdb_id,
            RANGE_COLUMN_NAME: "",
        }

    for csv_pdb_record in csv_pdb_record_list:
        csv_pdb_id = normalize_pdb_id(csv_pdb_record.get(PDB_ID_COLUMN_NAME, ""))

        if not csv_pdb_id:
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


def sync_pdb_csv_and_directories(
    protein_data_dir: Path,
    pdb_id_csv_filename: str = DEFAULT_PDB_ID_CSV_FILENAME,
) -> None:
    """
    Synchronize the PDB CSV file and the protein subdirectory structure.

    Cases
    -----
    Case 1:
        No CSV exists, but subdirectories exist
        -> create CSV from subdirectory names
        -> range values will be empty

    Case 2:
        CSV exists, but no subdirectories exist
        -> create subdirectories and download range-aware PDB files

    Case 3:
        Both CSV and subdirectories exist
        -> merge both sides
        -> preserve CSV range values
        -> add missing IDs from either side
        -> download missing PDB files
    """
    protein_data_dir.mkdir(parents=True, exist_ok=True)
    pdb_id_csv_path = protein_data_dir / pdb_id_csv_filename

    csv_file_exists = pdb_id_csv_path.exists()
    any_subdirectory_exists = any(
        filesystem_entry.is_dir() for filesystem_entry in protein_data_dir.iterdir()
    )

    print(f"[INFO] Working in data directory: {protein_data_dir}")
    print(f"[INFO] CSV path: {pdb_id_csv_path}")
    print(f"[INFO] CSV exists: {csv_file_exists}")
    print(f"[INFO] Any subdirectories exist: {any_subdirectory_exists}")

    if not csv_file_exists and any_subdirectory_exists:
        directory_pdb_record_list = get_pdb_records_from_subdirectories(
            protein_data_dir
        )

        print("[INFO] Case 1 detected: no CSV, but subdirectories exist")
        print(
            "[INFO] Creating CSV from subdirectories. All range values will be empty."
        )

        write_pdb_records_to_csv(directory_pdb_record_list, pdb_id_csv_path)

        print(f"[INFO] Created CSV from subdirectories: {pdb_id_csv_path}")
        return

    if csv_file_exists and not any_subdirectory_exists:
        csv_pdb_record_list = read_pdb_records_from_csv(pdb_id_csv_path)
        csv_pdb_id_list = extract_pdb_id_list_from_pdb_record_list(csv_pdb_record_list)

        print("[INFO] Case 2 detected: CSV exists, but no subdirectories")
        print(f"[INFO] CSV records: {csv_pdb_record_list}")

        create_missing_subdirectories(protein_data_dir, csv_pdb_id_list)
        download_missing_pdb_files(protein_data_dir, csv_pdb_record_list)

        print("[INFO] Created subdirectories and downloaded missing PDB files")
        return

    if not csv_file_exists and not any_subdirectory_exists:
        print("[WARNING] Neither CSV nor subdirectories exist.")
        print("[WARNING] Nothing to sync yet.")
        print(
            f"[WARNING] Create either a CSV file with column '{PDB_ID_COLUMN_NAME}' "
            "or create subdirectories first."
        )
        return

    csv_pdb_record_list = read_pdb_records_from_csv(pdb_id_csv_path)
    directory_pdb_record_list = get_pdb_records_from_subdirectories(protein_data_dir)

    csv_pdb_id_set = set(extract_pdb_id_list_from_pdb_record_list(csv_pdb_record_list))
    directory_pdb_id_set = set(
        extract_pdb_id_list_from_pdb_record_list(directory_pdb_record_list)
    )

    print("[INFO] Case 3 detected: CSV and subdirectories both exist")
    print(f"[INFO] CSV IDs: {sorted(csv_pdb_id_set)}")
    print(f"[INFO] Subdirectory IDs: {sorted(directory_pdb_id_set)}")

    missing_subdirectory_pdb_id_list = sorted(csv_pdb_id_set - directory_pdb_id_set)
    missing_csv_pdb_id_list = sorted(directory_pdb_id_set - csv_pdb_id_set)

    if missing_subdirectory_pdb_id_list:
        print(
            "[INFO] These IDs exist in CSV but are missing as subdirectories: "
            f"{missing_subdirectory_pdb_id_list}"
        )
        create_missing_subdirectories(
            protein_data_dir,
            missing_subdirectory_pdb_id_list,
        )

    if missing_csv_pdb_id_list:
        print(
            "[INFO] These IDs exist as subdirectories but are missing in CSV: "
            f"{missing_csv_pdb_id_list}"
        )

    merged_pdb_record_list = merge_csv_records_with_directory_records(
        csv_pdb_record_list,
        directory_pdb_record_list,
    )

    write_pdb_records_to_csv(merged_pdb_record_list, pdb_id_csv_path)
    print(f"[INFO] Updated CSV with all known IDs: {pdb_id_csv_path}")

    download_missing_pdb_files(protein_data_dir, merged_pdb_record_list)

    print("[INFO] Synchronization finished successfully")
