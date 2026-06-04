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
- If range is present as a plain integer range, for example '20-280':
    keep only ATOM residues whose residue number is within the inclusive range
    (applied to all chains).
- If range is present as a chain-aware range, for example 'A25-B15':
    chain A: keep residues whose residue number >= 25 (start endpoint)
    chain B: keep residues whose residue number <= 15 (end endpoint)
    chains not named in the range are dropped entirely.
    For a same-chain range like 'A25-A200', only chain A is kept (25–200).
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

import urllib.request
from pathlib import Path

from stack_protein_preparation._pdb_sync_core import (  # noqa: F401  (re-exported)
    DEFAULT_PDB_ID_CSV_FILENAME,
    IGNORED_SUBDIRECTORY_NAMES,
    MODULE_NAME,
    PDB_ID_COLUMN_NAME,
    RCSB_PDB_DOWNLOAD_URL_TEMPLATE,
    RANGE_COLUMN_NAME,
    ParsedRange,
    _append_dual_log,
    _append_log,
    _atom_line_passes_range,
    _chain_id,
    _get_logs_root_from_pdb_file_path,
    _get_logs_root_from_protein_data_dir,
    _is_valid_pdb_directory_name,
    _is_valid_pdb_id_text,
    _is_water_line,
    _iter_valid_protein_subdirectories,
    _log_exception,
    _parse_range_internal,
    _protein_log_path_from_pdb_id,
    _record_type,
    _resname,
    _resseq,
    _screen_header,
    _screen_result,
    _screen_step,
    _screen_sub,
    _session_log_path,
    _timestamp,
    create_missing_subdirectories,
    extract_pdb_id_list_from_pdb_record_list,
    get_observed_polymer_residue_bounds,
    get_pdb_records_from_subdirectories,
    normalize_pdb_id,
    normalize_range_value,
    parse_residue_range,
    parse_residue_range_chain_aware,
    read_pdb_records_from_csv,
    summarize_requested_vs_observed_range,
    trim_pdb_to_residue_range,
    write_pdb_records_to_csv,
)


def download_raw_pdb_file(pdb_id: str, target_pdb_file_path: Path) -> None:
    """Download a raw PDB file from the RCSB archive."""
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
    """Download a PDB file and apply optional range trimming."""
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
    """Download missing PDB files for a list of PDB records."""
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


def merge_csv_records_with_directory_records(
    csv_pdb_record_list: list[dict[str, str]],
    directory_pdb_record_list: list[dict[str, str]],
) -> list[dict[str, str]]:
    """Merge records from CSV and records inferred from directories."""
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
