"""
pipeline_table.py

Central table management for the pipeline.

Purpose
-------
This module manages the main pipeline table, which stores:

- pdb_id
- optional residue range
- future pipeline status fields

Design
------
- The table is stored as a CSV file
- Each row represents one protein
- Data is handled as a list of dictionaries

Important
---------
This module does NOT perform any scientific logic.
It only reads, writes, and updates table data.
"""

import csv
from pathlib import Path

PDB_ID_COLUMN_NAME = "pdb_id"
RANGE_COLUMN_NAME = "range"


def load_pipeline_table(csv_path: Path) -> list[dict[str, str]]:
    """
    Load the pipeline table from CSV.

    Returns
    -------
    list of records (dict)
    """
    records: list[dict[str, str]] = []

    if not csv_path.exists():
        return records

    with csv_path.open("r", encoding="utf-8", newline="") as file:
        reader = csv.DictReader(file)

        for row in reader:
            pdb_id = row.get(PDB_ID_COLUMN_NAME, "").strip().upper()
            range_value = row.get(RANGE_COLUMN_NAME, "").strip()

            if not pdb_id:
                continue

            record = {
                PDB_ID_COLUMN_NAME: pdb_id,
                RANGE_COLUMN_NAME: range_value,
            }

            records.append(record)

    return records


def save_pipeline_table(records: list[dict[str, str]], csv_path: Path) -> None:
    """
    Save records to CSV.

    Behavior
    --------
    - overwrites file
    - keeps only unique pdb_id entries
    - sorts by pdb_id
    """
    unique_records: dict[str, dict[str, str]] = {}

    for record in records:
        pdb_id = record.get(PDB_ID_COLUMN_NAME, "").strip().upper()

        if not pdb_id:
            continue

        range_value = record.get(RANGE_COLUMN_NAME, "").strip()

        unique_records[pdb_id] = {
            PDB_ID_COLUMN_NAME: pdb_id,
            RANGE_COLUMN_NAME: range_value,
        }

    sorted_ids = sorted(unique_records.keys())

    with csv_path.open("w", encoding="utf-8", newline="") as file:
        writer = csv.DictWriter(
            file,
            fieldnames=[PDB_ID_COLUMN_NAME, RANGE_COLUMN_NAME],
        )

        writer.writeheader()

        for pdb_id in sorted_ids:
            writer.writerow(unique_records[pdb_id])


def add_column_if_missing(
    records: list[dict[str, str]],
    column_name: str,
    default_value: str = "",
) -> None:
    """
    Ensure all records contain a specific column.

    Example:
    - later for fields like 'alignment_done'
    """
    for record in records:
        if column_name not in record:
            record[column_name] = default_value


def update_record(
    records: list[dict[str, str]],
    pdb_id: str,
    column_name: str,
    value: str,
) -> None:
    """
    Update a single value in the table.
    """
    normalized_pdb_id = pdb_id.strip().upper()

    for record in records:
        if record.get(PDB_ID_COLUMN_NAME) == normalized_pdb_id:
            record[column_name] = value
