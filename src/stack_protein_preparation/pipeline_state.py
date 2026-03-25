"""
/home/grheco/repositorios/stack_protein_prep/src/stack_protein_preparation/pipeline_state.py

Central state definition for the protein preparation pipeline.

Purpose
-------
This module defines which fields belong to one protein record in the pipeline.

Design choice
-------------
The state should stay small and practical.
Only fields that are actually useful for orchestration, debugging, JSON export,
and XLSX export should live here.

This module does NOT read or write files.
It only defines and prepares record structure.
"""

from __future__ import annotations

PDB_ID_COLUMN_NAME = "pdb_id"
RANGE_COLUMN_NAME = "range"

PDB_DIRECTORY_COLUMN_NAME = "pdb_directory"
FASTA_DIRECTORY_COLUMN_NAME = "fasta_directory"
ALIGNMENT_DIRECTORY_COLUMN_NAME = "alignment_directory"
COMPONENTS_DIRECTORY_COLUMN_NAME = "components_directory"
FILLER_DIRECTORY_COLUMN_NAME = "filler_directory"
PREPARED_DIRECTORY_COLUMN_NAME = "prepared_directory"

UNIPROT_ID_COLUMN_NAME = "uniprot_id"

N_GAPS_COLUMN_NAME = "n_gaps"
GAP_SIZES_COLUMN_NAME = "gap_sizes"

PDB_SYNC_DONE_COLUMN_NAME = "pdb_sync_done"
FASTA_FILES_DONE_COLUMN_NAME = "fasta_files_done"
SEQUENCE_ALIGNMENT_DONE_COLUMN_NAME = "sequence_alignment_done"
INSERTION_CODES_DONE_COLUMN_NAME = "insertion_codes_done"
NUMBERING_CHECK_DONE_COLUMN_NAME = "numbering_check_done"

HAS_GAPS_COLUMN_NAME = "has_gaps"
HAS_METALS_COLUMN_NAME = "has_metals"
HAS_LIGANDS_COLUMN_NAME = "has_ligands"
HAS_NONSTANDARD_RESIDUES_COLUMN_NAME = "has_nonstandard_residues"

FILLER_STATUS_COLUMN_NAME = "filler.status"
FILLER_MODEL_PATH_COLUMN_NAME = "filler.model_path"

PROTONATION_STATUS_COLUMN_NAME = "protonation.status"
PROTONATION_INPUT_SOURCE_COLUMN_NAME = "protonation.input_source"
PROTONATION_INPUT_PATH_COLUMN_NAME = "protonation.input_path"
PROTONATION_OUTPUT_PATH_COLUMN_NAME = "protonation.output_path"

AMBER_RENAMING_STATUS_COLUMN_NAME = "amber_renaming.status"
AMBER_INPUT_PATH_COLUMN_NAME = "amber_renaming.input_path"
AMBER_OUTPUT_PATH_COLUMN_NAME = "amber_renaming.output_path"

AMBER_TERMINI_STATUS_COLUMN_NAME = "amber_termini.status"
AMBER_TERMINI_INPUT_PATH_COLUMN_NAME = "amber_termini.input_path"
AMBER_TERMINI_OUTPUT_PATH_COLUMN_NAME = "amber_termini.output_path"

INTERNAL_CAPPING_STATUS_COLUMN_NAME = "internal_capping.status"
INTERNAL_CAPPING_INPUT_PATH_COLUMN_NAME = "internal_capping.input_path"
INTERNAL_CAPPING_OUTPUT_PATH_COLUMN_NAME = "internal_capping.output_path"

NUMBERING_RESTORE_STATUS_COLUMN_NAME = "numbering_restore.status"
NUMBERING_RESTORE_INPUT_PATH_COLUMN_NAME = "numbering_restore.input_path"
NUMBERING_RESTORE_OUTPUT_PATH_COLUMN_NAME = "numbering_restore.output_path"
NUMBERING_RESTORE_MAPPING_PATH_COLUMN_NAME = "numbering_restore.mapping_path"
NUMBERING_RESTORE_SOURCE_COLUMN_NAME = "numbering_restore.source"

PREPARED_STRUCTURE_STATUS_COLUMN_NAME = "prepared_structure.status"
PREPARED_STRUCTURE_VARIANT_COLUMN_NAME = "prepared_structure.variant"
PREPARED_STRUCTURE_PROTEIN_INPUT_PATH_COLUMN_NAME = (
    "prepared_structure.protein_input_path"
)
PREPARED_STRUCTURE_OUTPUT_PATH_COLUMN_NAME = "prepared_structure.output_path"

STATE_COLUMN_NAME_LIST = [
    PDB_ID_COLUMN_NAME,
    RANGE_COLUMN_NAME,
    PDB_DIRECTORY_COLUMN_NAME,
    FASTA_DIRECTORY_COLUMN_NAME,
    ALIGNMENT_DIRECTORY_COLUMN_NAME,
    COMPONENTS_DIRECTORY_COLUMN_NAME,
    FILLER_DIRECTORY_COLUMN_NAME,
    PREPARED_DIRECTORY_COLUMN_NAME,
    UNIPROT_ID_COLUMN_NAME,
    N_GAPS_COLUMN_NAME,
    GAP_SIZES_COLUMN_NAME,
    HAS_GAPS_COLUMN_NAME,
    PDB_SYNC_DONE_COLUMN_NAME,
    FASTA_FILES_DONE_COLUMN_NAME,
    SEQUENCE_ALIGNMENT_DONE_COLUMN_NAME,
    INSERTION_CODES_DONE_COLUMN_NAME,
    NUMBERING_CHECK_DONE_COLUMN_NAME,
    HAS_METALS_COLUMN_NAME,
    HAS_LIGANDS_COLUMN_NAME,
    HAS_NONSTANDARD_RESIDUES_COLUMN_NAME,
    FILLER_STATUS_COLUMN_NAME,
    FILLER_MODEL_PATH_COLUMN_NAME,
    PROTONATION_STATUS_COLUMN_NAME,
    PROTONATION_INPUT_SOURCE_COLUMN_NAME,
    PROTONATION_INPUT_PATH_COLUMN_NAME,
    PROTONATION_OUTPUT_PATH_COLUMN_NAME,
    AMBER_RENAMING_STATUS_COLUMN_NAME,
    AMBER_INPUT_PATH_COLUMN_NAME,
    AMBER_OUTPUT_PATH_COLUMN_NAME,
    AMBER_TERMINI_STATUS_COLUMN_NAME,
    AMBER_TERMINI_INPUT_PATH_COLUMN_NAME,
    AMBER_TERMINI_OUTPUT_PATH_COLUMN_NAME,
    INTERNAL_CAPPING_STATUS_COLUMN_NAME,
    INTERNAL_CAPPING_INPUT_PATH_COLUMN_NAME,
    INTERNAL_CAPPING_OUTPUT_PATH_COLUMN_NAME,
    NUMBERING_RESTORE_STATUS_COLUMN_NAME,
    NUMBERING_RESTORE_INPUT_PATH_COLUMN_NAME,
    NUMBERING_RESTORE_OUTPUT_PATH_COLUMN_NAME,
    NUMBERING_RESTORE_MAPPING_PATH_COLUMN_NAME,
    NUMBERING_RESTORE_SOURCE_COLUMN_NAME,
    PREPARED_STRUCTURE_STATUS_COLUMN_NAME,
    PREPARED_STRUCTURE_VARIANT_COLUMN_NAME,
    PREPARED_STRUCTURE_PROTEIN_INPUT_PATH_COLUMN_NAME,
    PREPARED_STRUCTURE_OUTPUT_PATH_COLUMN_NAME,
]

STEP_STATUS_COLUMN_NAME_LIST = [
    PDB_SYNC_DONE_COLUMN_NAME,
    FASTA_FILES_DONE_COLUMN_NAME,
    SEQUENCE_ALIGNMENT_DONE_COLUMN_NAME,
    INSERTION_CODES_DONE_COLUMN_NAME,
    NUMBERING_CHECK_DONE_COLUMN_NAME,
    FILLER_STATUS_COLUMN_NAME,
    PROTONATION_STATUS_COLUMN_NAME,
    AMBER_RENAMING_STATUS_COLUMN_NAME,
    AMBER_TERMINI_STATUS_COLUMN_NAME,
    INTERNAL_CAPPING_STATUS_COLUMN_NAME,
    NUMBERING_RESTORE_STATUS_COLUMN_NAME,
    PREPARED_STRUCTURE_STATUS_COLUMN_NAME,
]

STATUS_EMPTY = ""
STATUS_SUCCESS = "success"
STATUS_WARNING = "warning"
STATUS_REQUIRED = "required"
STATUS_SKIPPED = "skipped"
STATUS_FAILED = "failed"

ALLOWED_STATUS_VALUE_LIST = [
    STATUS_EMPTY,
    STATUS_SUCCESS,
    STATUS_WARNING,
    STATUS_REQUIRED,
    STATUS_SKIPPED,
    STATUS_FAILED,
]


def create_empty_protein_record() -> dict[str, str]:
    """
    Create one empty protein record with all expected columns.
    """
    return {column_name: "" for column_name in STATE_COLUMN_NAME_LIST}


def create_protein_record(
    pdb_id: str,
    residue_range: str = "",
) -> dict[str, str]:
    """
    Create one protein record with required initial values.
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
    """
    for protein_record in protein_record_list:
        for column_name in STATE_COLUMN_NAME_LIST:
            if column_name not in protein_record:
                protein_record[column_name] = ""


def validate_step_status_column_name(step_status_column_name: str) -> None:
    """
    Validate that a given column name is a known pipeline step status column.
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
    """
    validate_step_status_column_name(step_status_column_name)
    return str(protein_record.get(step_status_column_name, "")).strip()
