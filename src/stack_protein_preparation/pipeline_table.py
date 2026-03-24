# /home/grheco/repositorios/stack_protein_prep/src/stack_protein_preparation/pipeline_table.py

"""
/home/grheco/repositorios/stack_protein_prep/src/stack_protein_preparation/pipeline_table.py

Central table management for the pipeline.

Purpose
-------
This module manages the main pipeline table, which stores one record per protein.

Design
------
- the table is stored as a JSON file
- the JSON content is a list of dictionaries
- each dictionary represents one protein
- data is handled in memory as a list of dictionaries

Important
---------
This module does NOT perform any scientific logic.
It only reads, writes, and updates table data.
"""

from __future__ import annotations

import json
from pathlib import Path

from stack_protein_preparation.pipeline_state import (
    PDB_ID_COLUMN_NAME,
    RANGE_COLUMN_NAME,
    STATE_COLUMN_NAME_LIST,
    create_protein_record,
    ensure_all_state_columns_exist,
)


def load_pipeline_table(json_path: Path) -> list[dict[str, str]]:
    """
    Load the pipeline table from JSON.

    Expected JSON structure
    -----------------------
    A list of dictionaries, for example:

    [
        {
            "pdb_id": "1ABC",
            "range": "5-280",
            "pdb_sync_done": "success"
        },
        {
            "pdb_id": "2XYZ",
            "range": "",
            "sequence_alignment_done": "required"
        }
    ]

    Behavior
    --------
    - returns an empty list if the file does not exist
    - returns an empty list if the file is empty
    - normalizes pdb_id to uppercase
    - preserves all known pipeline state columns
    - repairs missing columns in older files

    Parameters
    ----------
    json_path
        Path to the pipeline JSON file.

    Returns
    -------
    list[dict[str, str]]
        List of protein records.

    Raises
    ------
    ValueError
        If the JSON root object is not a list.
    """
    json_path = Path(json_path)
    protein_record_list: list[dict[str, str]] = []

    if not json_path.exists():
        return protein_record_list

    raw_text = json_path.read_text(encoding="utf-8").strip()

    if not raw_text:
        return protein_record_list

    loaded_data = json.loads(raw_text)

    if not isinstance(loaded_data, list):
        raise ValueError(
            f"Pipeline JSON must contain a list of records, got {type(loaded_data).__name__}"
        )

    for item in loaded_data:
        if not isinstance(item, dict):
            continue

        pdb_id = str(item.get(PDB_ID_COLUMN_NAME, "")).strip().upper()
        residue_range = str(item.get(RANGE_COLUMN_NAME, "")).strip()

        if not pdb_id:
            continue

        protein_record = create_protein_record(
            pdb_id=pdb_id,
            residue_range=residue_range,
        )

        for column_name in STATE_COLUMN_NAME_LIST:
            if column_name not in item:
                continue

            raw_value = item.get(column_name, "")
            protein_record[column_name] = str(raw_value).strip()

        protein_record_list.append(protein_record)

    ensure_all_state_columns_exist(protein_record_list)
    return protein_record_list


def save_pipeline_table(
    protein_record_list: list[dict[str, str]],
    json_path: Path,
) -> None:
    """
    Save protein records to JSON.

    Behavior
    --------
    - overwrites the file
    - keeps only unique pdb_id entries
    - sorts records by pdb_id
    - writes all known pipeline state columns
    - creates parent directory if needed

    Parameters
    ----------
    protein_record_list
        List of protein records to save.
    json_path
        Path to the output JSON file.
    """
    json_path = Path(json_path)

    working_record_list = [
        dict(protein_record) for protein_record in protein_record_list
    ]
    ensure_all_state_columns_exist(working_record_list)

    unique_record_by_pdb_id: dict[str, dict[str, str]] = {}

    for protein_record in working_record_list:
        pdb_id = str(protein_record.get(PDB_ID_COLUMN_NAME, "")).strip().upper()

        if not pdb_id:
            continue

        normalized_record = create_protein_record(
            pdb_id=pdb_id,
            residue_range=str(protein_record.get(RANGE_COLUMN_NAME, "")).strip(),
        )

        for column_name in STATE_COLUMN_NAME_LIST:
            raw_value = protein_record.get(column_name, "")
            normalized_record[column_name] = str(raw_value).strip()

        unique_record_by_pdb_id[pdb_id] = normalized_record

    sorted_pdb_id_list = sorted(unique_record_by_pdb_id.keys())

    output_record_list = [
        unique_record_by_pdb_id[pdb_id] for pdb_id in sorted_pdb_id_list
    ]

    json_path.parent.mkdir(parents=True, exist_ok=True)
    json_path.write_text(
        json.dumps(output_record_list, indent=4, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )


def get_record_by_pdb_id(
    protein_record_list: list[dict[str, str]],
    pdb_id: str,
) -> dict[str, str] | None:
    """
    Find one protein record by PDB ID.

    Parameters
    ----------
    protein_record_list
        List of protein records.
    pdb_id
        PDB identifier to search for.

    Returns
    -------
    dict[str, str] | None
        Matching protein record, or None if not found.
    """
    normalized_pdb_id = str(pdb_id).strip().upper()

    for protein_record in protein_record_list:
        current_pdb_id = str(protein_record.get(PDB_ID_COLUMN_NAME, "")).strip().upper()

        if current_pdb_id == normalized_pdb_id:
            return protein_record

    return None


def update_record(
    protein_record_list: list[dict[str, str]],
    pdb_id: str,
    column_name: str,
    value: str,
) -> None:
    """
    Update one value in one protein record.

    Parameters
    ----------
    protein_record_list
        List of protein records.
    pdb_id
        PDB identifier of the record to update.
    column_name
        Column name to update.
    value
        New value to store.

    Raises
    ------
    KeyError
        If the protein record does not exist.
    ValueError
        If the column name is not a known pipeline state column.
    """
    if column_name not in STATE_COLUMN_NAME_LIST:
        allowed_column_name_string = ", ".join(STATE_COLUMN_NAME_LIST)
        raise ValueError(
            f"Unknown pipeline state column: '{column_name}'. "
            f"Allowed columns: {allowed_column_name_string}"
        )

    protein_record = get_record_by_pdb_id(protein_record_list, pdb_id)

    if protein_record is None:
        raise KeyError(f"No protein record found for pdb_id='{pdb_id}'")

    protein_record[column_name] = str(value).strip()
