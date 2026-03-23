# /home/grheco/repositorios/stack_protein_prep/src/stack_protein_preparation/pipeline_xlsx.py

"""
pipeline_xlsx.py

Excel export for pipeline state.

Purpose
-------
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
- non-status metadata columns are not colored
- only status columns are colored
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

from openpyxl import Workbook
from openpyxl.styles import PatternFill

from stack_protein_preparation.pipeline_state import (
    PDB_ID_COLUMN_NAME,
    STATE_COLUMN_NAME_LIST,
    STATUS_REQUIRED,
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


def get_excel_column_order() -> list[str]:
    """
    Return the fixed column order for the pipeline Excel file.

    Returns
    -------
    list[str]
        Ordered list of column names to write.
    """
    return list(STATE_COLUMN_NAME_LIST)


def ensure_all_records_have_all_columns(
    protein_record_list: list[dict[str, str]],
) -> None:
    """
    Ensure that every record has all expected pipeline columns.

    Parameters
    ----------
    protein_record_list
        List of protein records to repair in memory.
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

    Parameters
    ----------
    protein_record_list
        Input records.

    Returns
    -------
    list[dict[str, str]]
        Normalized record copies.
    """
    normalized_record_list: list[dict[str, str]] = []

    for protein_record in protein_record_list:
        normalized_record: dict[str, str] = {}

        for column_name in get_excel_column_order():
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

    Rule
    ----
    If the same pdb_id appears multiple times, the last record wins.

    Parameters
    ----------
    protein_record_list
        Input records.

    Returns
    -------
    list[dict[str, str]]
        Unique sorted records.
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


def apply_status_cell_color(cell, value: str) -> None:
    """
    Apply the agreed status color to one Excel cell.

    Parameters
    ----------
    cell
        Openpyxl worksheet cell.
    value
        Cell string value.
    """
    if value == STATUS_SUCCESS:
        cell.fill = GREEN_FILL
    elif value == STATUS_WARNING:
        cell.fill = YELLOW_FILL
    elif value == STATUS_REQUIRED:
        cell.fill = RED_FILL


def autosize_worksheet_columns(worksheet) -> None:
    """
    Set worksheet column widths based on content length.

    Parameters
    ----------
    worksheet
        Openpyxl worksheet object.
    """
    for column_cells in worksheet.columns:
        max_length = 0
        column_letter = column_cells[0].column_letter

        for cell in column_cells:
            cell_value = "" if cell.value is None else str(cell.value)
            max_length = max(max_length, len(cell_value))

        worksheet.column_dimensions[column_letter].width = max_length + 2


def write_pipeline_to_xlsx(
    protein_record_list: list[dict[str, str]],
    output_path: Path,
) -> None:
    """
    Write pipeline records to an XLSX file.

    Parameters
    ----------
    protein_record_list
        List of protein record dictionaries.
    output_path
        Target Excel file path.
    """
    working_record_list = [
        dict(protein_record) for protein_record in protein_record_list
    ]

    ensure_all_records_have_all_columns(working_record_list)
    normalized_record_list = normalize_record_values(working_record_list)
    unique_sorted_record_list = create_unique_sorted_record_list(normalized_record_list)

    column_order = get_excel_column_order()

    workbook = Workbook()
    worksheet = workbook.active
    worksheet.title = "pipeline"

    worksheet.append(column_order)

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

    autosize_worksheet_columns(worksheet)

    output_path.parent.mkdir(parents=True, exist_ok=True)
    workbook.save(output_path)
