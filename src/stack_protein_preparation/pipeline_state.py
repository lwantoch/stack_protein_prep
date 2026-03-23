"""
pipeline_state.py

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
- easy to later map to CSV or XLSX columns

This module does NOT read or write files.
It only defines and prepares record structure.
"""

PDB_ID_COLUMN_NAME = "pdb_id"
RANGE_COLUMN_NAME = "range"

PDB_SYNC_DONE_COLUMN_NAME = "pdb_sync_done"
SEQUENCE_ALIGNMENT_DONE_COLUMN_NAME = "sequence_alignment_done"
NUMBERING_CHECK_DONE_COLUMN_NAME = "numbering_check_done"
MODELLER_DONE_COLUMN_NAME = "modeller_done"
PROTONATION_DONE_COLUMN_NAME = "protonation_done"


def create_empty_protein_record() -> dict[str, str]:
    """
    Create one empty protein record with all expected columns.

    Returns
    -------
    dict[str, str]
        A record with default empty values.
    """
    protein_record = {
        PDB_ID_COLUMN_NAME: "",
        RANGE_COLUMN_NAME: "",
        PDB_SYNC_DONE_COLUMN_NAME: "",
        SEQUENCE_ALIGNMENT_DONE_COLUMN_NAME: "",
        NUMBERING_CHECK_DONE_COLUMN_NAME: "",
        MODELLER_DONE_COLUMN_NAME: "",
        PROTONATION_DONE_COLUMN_NAME: "",
    }

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

    normalized_pdb_id = pdb_id.strip().upper()
    normalized_residue_range = residue_range.strip()

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
    Older CSV/XLSX files may later miss newly added columns.
    This function can repair that in memory by adding missing keys.

    Parameters
    ----------
    protein_record_list
        List of protein records that should all follow the same structure.
    """
    empty_record = create_empty_protein_record()
    expected_column_name_list = list(empty_record.keys())

    for protein_record in protein_record_list:
        for column_name in expected_column_name_list:
            if column_name not in protein_record:
                protein_record[column_name] = ""


def mark_step_as_done(
    protein_record: dict[str, str],
    step_done_column_name: str,
) -> None:
    """
    Mark one pipeline step as done for a single protein record.

    Parameters
    ----------
    protein_record
        One protein record dictionary.
    step_done_column_name
        Name of the status column to update.
    """
    protein_record[step_done_column_name] = "yes"


def mark_step_as_not_done(
    protein_record: dict[str, str],
    step_done_column_name: str,
) -> None:
    """
    Mark one pipeline step as not done for a single protein record.

    Parameters
    ----------
    protein_record
        One protein record dictionary.
    step_done_column_name
        Name of the status column to update.
    """
    protein_record[step_done_column_name] = ""
