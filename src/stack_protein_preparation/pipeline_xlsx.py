"""
/home/grheco/repositorios/stack_protein_prep/src/stack_protein_preparation/pipeline_xlsx.py

Excel export for pipeline state.

Purpose
-------
- define the exported Excel column order
- avoid writing completely empty columns
- write pipeline records to an XLSX file
- apply clear, consistent status colors

Color convention
----------------
- "success"  -> green
- "warning"  -> yellow
- "required" -> red
- "skipped"  -> grey
- "failed"   -> orange-red
- empty      -> no fill

Important
---------
- only real status columns are color-filled
- metadata columns are not color-filled
- completely empty columns are removed from the export
- pdb_id is always kept, even if a user somehow gives empty records
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

from openpyxl import Workbook
from openpyxl.styles import Alignment, Font, PatternFill
from openpyxl.utils import get_column_letter

from stack_protein_preparation.pipeline_state import (
    PDB_ID_COLUMN_NAME,
    STATE_COLUMN_NAME_LIST,
    STATUS_FAILED,
    STATUS_REQUIRED,
    STATUS_SKIPPED,
    STATUS_SUCCESS,
    STATUS_WARNING,
    STEP_STATUS_COLUMN_NAME_LIST,
    create_empty_protein_record,
)

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

GREY_FILL = PatternFill(
    start_color="D9D9D9",
    end_color="D9D9D9",
    fill_type="solid",
)

ORANGE_RED_FILL = PatternFill(
    start_color="F4B084",
    end_color="F4B084",
    fill_type="solid",
)

HEADER_FILL = PatternFill(
    start_color="D9E2F3",
    end_color="D9E2F3",
    fill_type="solid",
)

HEADER_FONT = Font(bold=True)
HEADER_ALIGNMENT = Alignment(horizontal="center", vertical="center")
DATA_ALIGNMENT = Alignment(vertical="top")


def get_excel_column_order() -> list[str]:
    """
    Return the base column order for Excel export.
    """
    return list(STATE_COLUMN_NAME_LIST)


def ensure_all_records_have_all_columns(
    protein_record_list: list[dict[str, Any]],
) -> None:
    """
    Ensure that every record has all expected pipeline columns.
    """
    empty_record = create_empty_protein_record()

    for protein_record in protein_record_list:
        for column_name in empty_record:
            if column_name not in protein_record:
                protein_record[column_name] = ""


def normalize_record_values(
    protein_record_list: list[dict[str, Any]],
) -> list[dict[str, str]]:
    """
    Return normalized copies of protein records.

    Rules
    -----
    - pdb_id is uppercase
    - all values are converted to stripped strings
    - records with empty pdb_id are skipped
    """
    normalized_record_list: list[dict[str, str]] = []
    column_order = get_excel_column_order()

    for protein_record in protein_record_list:
        normalized_record: dict[str, str] = {}

        for column_name in column_order:
            raw_value = protein_record.get(column_name, "")
            normalized_value = str(raw_value).strip()

            if column_name == PDB_ID_COLUMN_NAME:
                normalized_value = normalized_value.upper()

            normalized_record[column_name] = normalized_value

        if not normalized_record[PDB_ID_COLUMN_NAME]:
            continue

        normalized_record_list.append(normalized_record)

    return normalized_record_list


def create_unique_sorted_record_list(
    protein_record_list: list[dict[str, str]],
) -> list[dict[str, str]]:
    """
    Return unique records sorted by pdb_id.

    If the same pdb_id appears multiple times, the last record wins.
    """
    record_by_pdb_id: dict[str, dict[str, str]] = {}

    for protein_record in protein_record_list:
        pdb_id = str(protein_record.get(PDB_ID_COLUMN_NAME, "")).strip().upper()
        if not pdb_id:
            continue

        copied_record = dict(protein_record)
        copied_record[PDB_ID_COLUMN_NAME] = pdb_id
        record_by_pdb_id[pdb_id] = copied_record

    return [record_by_pdb_id[pdb_id] for pdb_id in sorted(record_by_pdb_id.keys())]


def get_nonempty_column_order(
    protein_record_list: list[dict[str, str]],
) -> list[str]:
    """
    Return only columns that are actually needed in the export.

    Rules
    -----
    - always keep pdb_id
    - keep a column if at least one record has a non-empty value there
    """
    base_column_order = get_excel_column_order()
    filtered_column_order: list[str] = []

    for column_name in base_column_order:
        if column_name == PDB_ID_COLUMN_NAME:
            filtered_column_order.append(column_name)
            continue

        has_nonempty_value = any(
            str(protein_record.get(column_name, "")).strip()
            for protein_record in protein_record_list
        )

        if has_nonempty_value:
            filtered_column_order.append(column_name)

    return filtered_column_order


def apply_status_cell_color(cell, value: str) -> None:
    """
    Apply the agreed status color to one Excel cell.
    """
    normalized_value = str(value).strip()

    if normalized_value == STATUS_SUCCESS:
        cell.fill = GREEN_FILL
    elif normalized_value == STATUS_WARNING:
        cell.fill = YELLOW_FILL
    elif normalized_value == STATUS_REQUIRED:
        cell.fill = RED_FILL
    elif normalized_value == STATUS_SKIPPED:
        cell.fill = GREY_FILL
    elif normalized_value == STATUS_FAILED:
        cell.fill = ORANGE_RED_FILL


def style_header_row(worksheet, column_order: list[str]) -> None:
    """
    Apply styling to the header row.
    """
    for column_index, _column_name in enumerate(column_order, start=1):
        cell = worksheet.cell(row=1, column=column_index)
        cell.fill = HEADER_FILL
        cell.font = HEADER_FONT
        cell.alignment = HEADER_ALIGNMENT


def style_data_cells(worksheet) -> None:
    """
    Apply base alignment to non-header cells.
    """
    for row in worksheet.iter_rows(min_row=2):
        for cell in row:
            cell.alignment = DATA_ALIGNMENT


def autosize_worksheet_columns(worksheet) -> None:
    """
    Set worksheet column widths based on content length.
    """
    for column_index in range(1, worksheet.max_column + 1):
        column_letter = get_column_letter(column_index)
        max_length = 0

        for row_index in range(1, worksheet.max_row + 1):
            cell = worksheet.cell(row=row_index, column=column_index)
            cell_value = "" if cell.value is None else str(cell.value)
            max_length = max(max_length, len(cell_value))

        adjusted_width = min(max(max_length + 2, 12), 60)
        worksheet.column_dimensions[column_letter].width = adjusted_width


def write_pipeline_to_xlsx(
    protein_record_list: list[dict[str, Any]],
    output_path: Path,
) -> None:
    """
    Write pipeline records to an XLSX file.
    """
    working_record_list = [
        dict(protein_record) for protein_record in protein_record_list
    ]

    ensure_all_records_have_all_columns(working_record_list)
    normalized_record_list = normalize_record_values(working_record_list)
    unique_sorted_record_list = create_unique_sorted_record_list(normalized_record_list)

    column_order = get_nonempty_column_order(unique_sorted_record_list)

    workbook = Workbook()
    worksheet = workbook.active
    worksheet.title = "pipeline"
    worksheet.freeze_panes = "A2"

    worksheet.append(column_order)
    style_header_row(worksheet, column_order)

    for protein_record in unique_sorted_record_list:
        row_value_list = [
            protein_record.get(column_name, "") for column_name in column_order
        ]
        worksheet.append(row_value_list)

        current_row_index = worksheet.max_row

        for column_index, column_name in enumerate(column_order, start=1):
            if column_name not in STEP_STATUS_COLUMN_NAME_LIST:
                continue

            cell = worksheet.cell(row=current_row_index, column=column_index)
            cell_value = "" if cell.value is None else str(cell.value).strip()
            apply_status_cell_color(cell, cell_value)

    style_data_cells(worksheet)
    autosize_worksheet_columns(worksheet)

    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    workbook.save(output_path)
