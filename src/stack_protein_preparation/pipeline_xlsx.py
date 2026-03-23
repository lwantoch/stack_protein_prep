"""
pipeline_xlsx.py

Excel export for pipeline state.

Purpose
-------
- create pipeline records
- define Excel column order
- write pipeline records to an XLSX file
- color status cells using the agreed color logic

Color convention
----------------
- "success"  -> green
- "warning"  -> yellow
- "required" -> red

Important
---------
- pdb_id and range are not colored
- only status columns are colored
"""

from pathlib import Path

from openpyxl import Workbook
from openpyxl.styles import PatternFill

PDB_ID_COLUMN_NAME = "pdb_id"
RANGE_COLUMN_NAME = "range"

PDB_SYNC_DONE_COLUMN_NAME = "pdb_sync_done"
SEQUENCE_ALIGNMENT_DONE_COLUMN_NAME = "sequence_alignment_done"
NUMBERING_CHECK_DONE_COLUMN_NAME = "numbering_check_done"
MODELLER_DONE_COLUMN_NAME = "modeller_done"
PROTONATION_DONE_COLUMN_NAME = "protonation_done"


GREEN_FILL = PatternFill(
    start_color="C6EFCE",
    end_color="C6EFCE",
    fill_type="solid",
)

YELLOW_FILL = PatternFill(
    start_color="FFEB9C",
    end_color="FFEB9C",
    fill_type="solid",
)

RED_FILL = PatternFill(
    start_color="FFC7CE",
    end_color="FFC7CE",
    fill_type="solid",
)


def get_column_order() -> list[str]:
    """
    Return the fixed column order for the pipeline Excel file.
    """
    column_order = [
        PDB_ID_COLUMN_NAME,
        RANGE_COLUMN_NAME,
        PDB_SYNC_DONE_COLUMN_NAME,
        SEQUENCE_ALIGNMENT_DONE_COLUMN_NAME,
        NUMBERING_CHECK_DONE_COLUMN_NAME,
        MODELLER_DONE_COLUMN_NAME,
        PROTONATION_DONE_COLUMN_NAME,
    ]

    return column_order


def create_empty_record() -> dict[str, str]:
    """
    Create one empty pipeline record with all expected columns.
    """
    empty_record = {}

    for column_name in get_column_order():
        empty_record[column_name] = ""

    return empty_record


def create_pipeline_record_from_input_values(
    pdb_id: str,
    residue_range: str = "",
) -> dict[str, str]:
    """
    Create one pipeline record from the basic input fields.

    Parameters
    ----------
    pdb_id
        Protein Data Bank identifier.
    residue_range
        Optional residue range string.

    Returns
    -------
    dict[str, str]
        A complete pipeline record with empty status fields.
    """
    pipeline_record = create_empty_record()

    normalized_pdb_id = str(pdb_id).strip().upper()
    normalized_residue_range = str(residue_range).strip()

    pipeline_record[PDB_ID_COLUMN_NAME] = normalized_pdb_id
    pipeline_record[RANGE_COLUMN_NAME] = normalized_residue_range

    return pipeline_record


def ensure_all_records_have_all_columns(
    pipeline_record_list: list[dict[str, str]],
) -> None:
    """
    Ensure that every record has all expected columns.
    """
    expected_column_name_list = get_column_order()

    for pipeline_record in pipeline_record_list:
        for column_name in expected_column_name_list:
            if column_name not in pipeline_record:
                pipeline_record[column_name] = ""


def normalize_records(
    pipeline_record_list: list[dict[str, str]],
) -> None:
    """
    Normalize record values before writing.

    Rules
    -----
    - pdb_id is uppercase
    - all values are strings without surrounding whitespace
    """
    for pipeline_record in pipeline_record_list:
        raw_pdb_id = pipeline_record.get(PDB_ID_COLUMN_NAME, "")
        normalized_pdb_id = str(raw_pdb_id).strip().upper()
        pipeline_record[PDB_ID_COLUMN_NAME] = normalized_pdb_id

        for column_name in pipeline_record:
            if column_name == PDB_ID_COLUMN_NAME:
                continue

            raw_value = pipeline_record.get(column_name, "")
            normalized_value = str(raw_value).strip()
            pipeline_record[column_name] = normalized_value


def sort_records_by_pdb_id(
    pipeline_record_list: list[dict[str, str]],
) -> list[dict[str, str]]:
    """
    Return a new list of records sorted by pdb_id.

    Records with empty pdb_id are ignored.
    """
    sortable_record_list: list[dict[str, str]] = []

    for pipeline_record in pipeline_record_list:
        pdb_id = pipeline_record.get(PDB_ID_COLUMN_NAME, "").strip().upper()

        if not pdb_id:
            continue

        sortable_record_list.append(pipeline_record)

    sorted_record_list = sorted(
        sortable_record_list,
        key=lambda pipeline_record: pipeline_record[PDB_ID_COLUMN_NAME],
    )

    return sorted_record_list


def create_unique_record_list(
    pipeline_record_list: list[dict[str, str]],
) -> list[dict[str, str]]:
    """
    Keep only one record per pdb_id.

    Rule
    ----
    If the same pdb_id appears multiple times, the last record wins.
    """
    pdb_id_to_record_dict: dict[str, dict[str, str]] = {}

    for pipeline_record in pipeline_record_list:
        pdb_id = pipeline_record.get(PDB_ID_COLUMN_NAME, "").strip().upper()

        if not pdb_id:
            continue

        copied_record = dict(pipeline_record)
        copied_record[PDB_ID_COLUMN_NAME] = pdb_id

        pdb_id_to_record_dict[pdb_id] = copied_record

    unique_record_list = list(pdb_id_to_record_dict.values())
    unique_record_list = sort_records_by_pdb_id(unique_record_list)

    return unique_record_list


def apply_cell_color(cell, value: str) -> None:
    """
    Apply the agreed status color to one Excel cell.
    """
    if value == "success":
        cell.fill = GREEN_FILL
    elif value == "warning":
        cell.fill = YELLOW_FILL
    elif value == "required":
        cell.fill = RED_FILL


def write_pipeline_to_xlsx(
    pipeline_record_list: list[dict[str, str]],
    output_path: Path,
) -> None:
    """
    Write pipeline records to an XLSX file.

    Parameters
    ----------
    pipeline_record_list
        List of pipeline record dictionaries.
    output_path
        Target Excel file path.
    """
    working_record_list: list[dict[str, str]] = []

    for pipeline_record in pipeline_record_list:
        copied_record = dict(pipeline_record)
        working_record_list.append(copied_record)

    ensure_all_records_have_all_columns(working_record_list)
    normalize_records(working_record_list)
    unique_sorted_record_list = create_unique_record_list(working_record_list)

    column_order = get_column_order()

    workbook = Workbook()
    worksheet = workbook.active
    worksheet.title = "pipeline"

    worksheet.append(column_order)

    for pipeline_record in unique_sorted_record_list:
        row_value_list: list[str] = []

        for column_name in column_order:
            cell_value = pipeline_record.get(column_name, "")
            row_value_list.append(cell_value)

        worksheet.append(row_value_list)

        current_row_index = worksheet.max_row

        for column_index, column_name in enumerate(column_order, start=1):
            if column_name in [PDB_ID_COLUMN_NAME, RANGE_COLUMN_NAME]:
                continue

            cell = worksheet.cell(row=current_row_index, column=column_index)
            cell_value = str(cell.value)

            apply_cell_color(cell, cell_value)

    output_path.parent.mkdir(parents=True, exist_ok=True)
    workbook.save(output_path)
