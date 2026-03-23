# /home/grheco/repositorios/stack_protein_prep/src/stack_protein_preparation/pipeline_state.py

"""
Central state definition for the protein preparation pipeline.

Purpose
-------
This module defines which fields belong to one protein record in the pipeline.

Important idea
--------------
At this stage, the pipeline state is represented as a normal Python dictionary.

Why this design
---------------
- easy to understand
- easy to debug
- easy to later map to JSON or XLSX columns

This module does NOT read or write files.
It only defines and prepares record structure.
"""

from __future__ import annotations

PDB_ID_COLUMN_NAME = "pdb_id"
RANGE_COLUMN_NAME = "range"

PDB_DIRECTORY_COLUMN_NAME = "pdb_directory"
FASTA_DIRECTORY_COLUMN_NAME = "fasta_directory"
ALIGNMENT_DIRECTORY_COLUMN_NAME = "alignment_directory"

PDB_SYNC_DONE_COLUMN_NAME = "pdb_sync_done"
FASTA_FILES_DONE_COLUMN_NAME = "fasta_files_done"
SEQUENCE_ALIGNMENT_DONE_COLUMN_NAME = "sequence_alignment_done"
NUMBERING_CHECK_DONE_COLUMN_NAME = "numbering_check_done"
MODELLER_DONE_COLUMN_NAME = "modeller_done"
PROTONATION_DONE_COLUMN_NAME = "protonation_done"

STATE_COLUMN_NAME_LIST = [
    PDB_ID_COLUMN_NAME,
    RANGE_COLUMN_NAME,
    PDB_DIRECTORY_COLUMN_NAME,
    FASTA_DIRECTORY_COLUMN_NAME,
    ALIGNMENT_DIRECTORY_COLUMN_NAME,
    PDB_SYNC_DONE_COLUMN_NAME,
    FASTA_FILES_DONE_COLUMN_NAME,
    SEQUENCE_ALIGNMENT_DONE_COLUMN_NAME,
    NUMBERING_CHECK_DONE_COLUMN_NAME,
    MODELLER_DONE_COLUMN_NAME,
    PROTONATION_DONE_COLUMN_NAME,
]

STEP_STATUS_COLUMN_NAME_LIST = [
    PDB_SYNC_DONE_COLUMN_NAME,
    FASTA_FILES_DONE_COLUMN_NAME,
    SEQUENCE_ALIGNMENT_DONE_COLUMN_NAME,
    NUMBERING_CHECK_DONE_COLUMN_NAME,
    MODELLER_DONE_COLUMN_NAME,
    PROTONATION_DONE_COLUMN_NAME,
]

STATUS_EMPTY = ""
STATUS_SUCCESS = "success"
STATUS_WARNING = "warning"
STATUS_REQUIRED = "required"

ALLOWED_STATUS_VALUE_LIST = [
    STATUS_EMPTY,
    STATUS_SUCCESS,
    STATUS_WARNING,
    STATUS_REQUIRED,
]


def create_empty_protein_record() -> dict[str, str]:
    """
    Create one empty protein record with all expected columns.

    Returns
    -------
    dict[str, str]
        A record with default empty values.
    """
    protein_record = {column_name: "" for column_name in STATE_COLUMN_NAME_LIST}
    return protein_record


def create_protein_record(
    pdb_id: str,
    residue_range: str = "",
) -> dict[str, str]:
    """
    Create one protein record with required initial values.

    Parameters
    ----------
    pdb_id
        Protein Data Bank identifier.
    residue_range
        Optional residue range string from the input table.

    Returns
    -------
    dict[str, str]
        One initialized protein record.
    """
    protein_record = create_empty_protein_record()

    normalized_pdb_id = str(pdb_id).strip().upper()
    normalized_residue_range = str(residue_range).strip()

    protein_record[PDB_ID_COLUMN_NAME] = normalized_pdb_id
    protein_record[RANGE_COLUMN_NAME] = normalized_residue_range

    return protein_record


def ensure_all_state_columns_exist(
    protein_record_list: list[dict[str, str]],
) -> None:
    """
    Ensure that every record contains all expected state columns.

    Why this is useful
    ------------------
    Older JSON/XLSX files may later miss newly added columns.
    This function can repair that in memory by adding missing keys.

    Parameters
    ----------
    protein_record_list
        List of protein records that should all follow the same structure.
    """
    for protein_record in protein_record_list:
        for column_name in STATE_COLUMN_NAME_LIST:
            if column_name not in protein_record:
                protein_record[column_name] = ""


def validate_step_status_column_name(step_status_column_name: str) -> None:
    """
    Validate that a given column name is a known pipeline step status column.

    Parameters
    ----------
    step_status_column_name
        Column name to validate.

    Raises
    ------
    ValueError
        If the provided column name is not a known step status column.
    """
    if step_status_column_name not in STEP_STATUS_COLUMN_NAME_LIST:
        allowed_column_name_string = ", ".join(STEP_STATUS_COLUMN_NAME_LIST)
        raise ValueError(
            f"Unknown step status column: '{step_status_column_name}'. "
            f"Allowed columns: {allowed_column_name_string}"
        )


def validate_status_value(status_value: str) -> None:
    """
    Validate that a given status value is allowed.

    Parameters
    ----------
    status_value
        Status value to validate.

    Raises
    ------
    ValueError
        If the provided status value is not allowed.
    """
    if status_value not in ALLOWED_STATUS_VALUE_LIST:
        allowed_status_value_string = ", ".join(
            repr(value) for value in ALLOWED_STATUS_VALUE_LIST
        )
        raise ValueError(
            f"Unknown status value: {status_value!r}. "
            f"Allowed values: {allowed_status_value_string}"
        )


def set_step_status(
    protein_record: dict[str, str],
    step_status_column_name: str,
    status_value: str,
) -> None:
    """
    Set one pipeline step status for a single protein record.

    Parameters
    ----------
    protein_record
        One protein record dictionary.
    step_status_column_name
        Name of the status column to update.
    status_value
        New status value.
    """
    validate_step_status_column_name(step_status_column_name)
    validate_status_value(status_value)
    protein_record[step_status_column_name] = status_value


def get_step_status(
    protein_record: dict[str, str],
    step_status_column_name: str,
) -> str:
    """
    Return the current status value of one pipeline step.

    Parameters
    ----------
    protein_record
        One protein record dictionary.
    step_status_column_name
        Name of the status column to inspect.

    Returns
    -------
    str
        Current status value.
    """
    validate_step_status_column_name(step_status_column_name)
    return str(protein_record.get(step_status_column_name, "")).strip()
