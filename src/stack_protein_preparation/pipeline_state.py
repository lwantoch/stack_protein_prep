"""
/home/grheco/repositorios/stack_protein_prep/src/stack_protein_preparation/pipeline_state.py

Central state definition for the protein preparation pipeline.

Purpose
-------
This module defines which fields belong to one protein record in the pipeline.

Design choice
-------------
The state should stay practical and flat enough for:
- orchestration
- JSON persistence
- XLSX export
- debugging
- downstream branching

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
METALL_PARAMS_DIRECTORY_COLUMN_NAME = "metall_params_directory"

UNIPROT_ID_COLUMN_NAME = "uniprot_id"

N_GAPS_COLUMN_NAME = "n_gaps"
GAP_SIZES_COLUMN_NAME = "gap_sizes"
HAS_GAPS_COLUMN_NAME = "has_gaps"
GAP_BOUNDARIES_COLUMN_NAME = "gap_boundaries"

PDB_SYNC_DONE_COLUMN_NAME = "pdb_sync_done"
FASTA_FILES_DONE_COLUMN_NAME = "fasta_files_done"
SEQUENCE_ALIGNMENT_DONE_COLUMN_NAME = "sequence_alignment_done"
INSERTION_CODES_DONE_COLUMN_NAME = "insertion_codes_done"

HAS_METALS_COLUMN_NAME = "has_metals"
HAS_LIGANDS_COLUMN_NAME = "has_ligands"
HAS_NONSTANDARD_RESIDUES_COLUMN_NAME = "has_nonstandard_residues"

FILLER_STATUS_COLUMN_NAME = "filler.status"
FILLER_MODEL_PATH_COLUMN_NAME = "filler.model_path"
FILLER_MODEL_SOURCE_COLUMN_NAME = "filler.model_source"

PROTONATION_STATUS_COLUMN_NAME = "protonation.status"
PROTONATION_INPUT_SOURCE_COLUMN_NAME = "protonation.input_source"
PROTONATION_INPUT_PATH_COLUMN_NAME = "protonation.input_path"
PROTONATION_OUTPUT_PATH_COLUMN_NAME = "protonation.output_path"
PROTONATION_PROPKA_HIS_COLUMN_NAME = "protonation.propka_his_assignments"

INTERNAL_CAPPING_STATUS_COLUMN_NAME = "internal_capping.status"
INTERNAL_CAPPING_INPUT_PATH_COLUMN_NAME = "internal_capping.input_path"
INTERNAL_CAPPING_OUTPUT_PATH_COLUMN_NAME = "internal_capping.output_path"

PREPARED_STRUCTURE_STATUS_COLUMN_NAME = "prepared_structure.status"
PREPARED_STRUCTURE_VARIANT_COLUMN_NAME = "prepared_structure.variant"
PREPARED_STRUCTURE_PROTEIN_INPUT_PATH_COLUMN_NAME = (
    "prepared_structure.protein_input_path"
)
PREPARED_STRUCTURE_OUTPUT_PATH_COLUMN_NAME = "prepared_structure.output_path"

# Variant logic summary
VARIANT_POLICY_COLUMN_NAME = "variant.policy"
RETAIN_GAPS_VARIANT_COLUMN_NAME = "variant.retain_gaps"
RETAIN_MODELLER_VARIANT_COLUMN_NAME = "variant.retain_modeller"
RETAIN_ALPHAFOLD_VARIANT_COLUMN_NAME = "variant.retain_alphafold"

# Optional future-ready per-variant prepared outputs
PREPARED_GAPS_OUTPUT_PATH_COLUMN_NAME = "prepared.gaps.output_path"
PREPARED_MODELLER_OUTPUT_PATH_COLUMN_NAME = "prepared.modeller.output_path"
PREPARED_ALPHAFOLD_OUTPUT_PATH_COLUMN_NAME = "prepared.alphafold.output_path"

# Human-readable final summary of which model families were actually produced.
# Example values:
# - "initial"
# - "gaps | modeller"
# - "gaps | alphafold"
# - "initial | gaps | modeller"
AVAILABLE_MODELS_COLUMN_NAME = "available_models"


# Metal check / geometry summary
METALS_CHECK_STATUS_COLUMN_NAME = "metals_check.status"
METALS_ION_TYPE_COLUMN_NAME = "metals.ion_type"
METALS_CLASS_COLUMN_NAME = "metals.class"
METALS_PARAMETER_REFERENCE_COLUMN_NAME = "metals.parameter_reference"
METALS_GEOMETRY_FOUND_COLUMN_NAME = "metals.geometry_found"
METALS_GEOMETRY_PROBABLE_COLUMN_NAME = "metals.geometry_probable"
METALS_GEOMETRY_MATCH_COLUMN_NAME = "metals.geometry_match"
METALS_MODEL_READY_COLUMN_NAME = "metals.model_ready"
METALS_CHECK_LOG_PATH_COLUMN_NAME = "metals.check_log_path"
METALS_CHECK_MANIFEST_PATH_COLUMN_NAME = "metals.check_manifest_path"

# Metal parametrization branch
METALL_PARAMS_STATUS_COLUMN_NAME = "metall_params.status"
METALL_PARAMS_SITE_COUNT_COLUMN_NAME = "metall_params.site_count"
METALL_PARAMS_MANIFEST_PATH_COLUMN_NAME = "metall_params.manifest_path"

# Non-standard residue parametrization branch
NONSTD_RESIDUE_PARAMS_STATUS_COLUMN_NAME = "nonstd_residue_params.status"
NONSTD_RESIDUE_PARAMS_N_RESIDUES_COLUMN_NAME = "nonstd_residue_params.n_residues"
NONSTD_RESIDUE_PARAMS_MANIFEST_PATH_COLUMN_NAME = "nonstd_residue_params.manifest_path"

# Model evaluation (Ramachandran, clashscore) — only for Modeller-filled proteins
MODEL_EVALUATION_STATUS_COLUMN_NAME = "model_eval.status"
MODEL_EVALUATION_PCT_FAVORED_COLUMN_NAME = "model_eval.rama_pct_favored"
MODEL_EVALUATION_PCT_OUTLIER_COLUMN_NAME = "model_eval.rama_pct_outlier"
MODEL_EVALUATION_CLASHSCORE_COLUMN_NAME = "model_eval.clashscore"
MODEL_EVALUATION_OVERALL_QUALITY_COLUMN_NAME = "model_eval.overall_quality"
MODEL_EVALUATION_RAMA_PLOT_PATH_COLUMN_NAME = "model_eval.rama_plot_path"
MODEL_EVALUATION_MANIFEST_PATH_COLUMN_NAME = "model_eval.manifest_path"

# Per-protein PDF report
REPORT_STATUS_COLUMN_NAME = "report.status"
REPORT_PDF_PATH_COLUMN_NAME = "report.pdf_path"
REPORT_FIGURE_PNG_PATH_COLUMN_NAME = "report.figure_png_path"

STATE_COLUMN_NAME_LIST = [
    PDB_ID_COLUMN_NAME,
    RANGE_COLUMN_NAME,
    PDB_DIRECTORY_COLUMN_NAME,
    FASTA_DIRECTORY_COLUMN_NAME,
    ALIGNMENT_DIRECTORY_COLUMN_NAME,
    COMPONENTS_DIRECTORY_COLUMN_NAME,
    FILLER_DIRECTORY_COLUMN_NAME,
    PREPARED_DIRECTORY_COLUMN_NAME,
    METALL_PARAMS_DIRECTORY_COLUMN_NAME,
    UNIPROT_ID_COLUMN_NAME,
    N_GAPS_COLUMN_NAME,
    GAP_SIZES_COLUMN_NAME,
    HAS_GAPS_COLUMN_NAME,
    GAP_BOUNDARIES_COLUMN_NAME,
    PDB_SYNC_DONE_COLUMN_NAME,
    FASTA_FILES_DONE_COLUMN_NAME,
    SEQUENCE_ALIGNMENT_DONE_COLUMN_NAME,
    INSERTION_CODES_DONE_COLUMN_NAME,
    HAS_METALS_COLUMN_NAME,
    HAS_LIGANDS_COLUMN_NAME,
    HAS_NONSTANDARD_RESIDUES_COLUMN_NAME,
    FILLER_STATUS_COLUMN_NAME,
    FILLER_MODEL_PATH_COLUMN_NAME,
    FILLER_MODEL_SOURCE_COLUMN_NAME,
    PROTONATION_STATUS_COLUMN_NAME,
    PROTONATION_INPUT_SOURCE_COLUMN_NAME,
    PROTONATION_INPUT_PATH_COLUMN_NAME,
    PROTONATION_OUTPUT_PATH_COLUMN_NAME,
    PROTONATION_PROPKA_HIS_COLUMN_NAME,
    INTERNAL_CAPPING_STATUS_COLUMN_NAME,
    INTERNAL_CAPPING_INPUT_PATH_COLUMN_NAME,
    INTERNAL_CAPPING_OUTPUT_PATH_COLUMN_NAME,
    PREPARED_STRUCTURE_STATUS_COLUMN_NAME,
    PREPARED_STRUCTURE_VARIANT_COLUMN_NAME,
    PREPARED_STRUCTURE_PROTEIN_INPUT_PATH_COLUMN_NAME,
    PREPARED_STRUCTURE_OUTPUT_PATH_COLUMN_NAME,
    VARIANT_POLICY_COLUMN_NAME,
    RETAIN_GAPS_VARIANT_COLUMN_NAME,
    RETAIN_MODELLER_VARIANT_COLUMN_NAME,
    RETAIN_ALPHAFOLD_VARIANT_COLUMN_NAME,
    PREPARED_GAPS_OUTPUT_PATH_COLUMN_NAME,
    PREPARED_MODELLER_OUTPUT_PATH_COLUMN_NAME,
    PREPARED_ALPHAFOLD_OUTPUT_PATH_COLUMN_NAME,
    AVAILABLE_MODELS_COLUMN_NAME,
    METALS_CHECK_STATUS_COLUMN_NAME,
    METALS_ION_TYPE_COLUMN_NAME,
    METALS_CLASS_COLUMN_NAME,
    METALS_PARAMETER_REFERENCE_COLUMN_NAME,
    METALS_GEOMETRY_FOUND_COLUMN_NAME,
    METALS_GEOMETRY_PROBABLE_COLUMN_NAME,
    METALS_GEOMETRY_MATCH_COLUMN_NAME,
    METALS_MODEL_READY_COLUMN_NAME,
    METALS_CHECK_LOG_PATH_COLUMN_NAME,
    METALS_CHECK_MANIFEST_PATH_COLUMN_NAME,
    METALL_PARAMS_STATUS_COLUMN_NAME,
    METALL_PARAMS_SITE_COUNT_COLUMN_NAME,
    METALL_PARAMS_MANIFEST_PATH_COLUMN_NAME,
    NONSTD_RESIDUE_PARAMS_STATUS_COLUMN_NAME,
    NONSTD_RESIDUE_PARAMS_N_RESIDUES_COLUMN_NAME,
    NONSTD_RESIDUE_PARAMS_MANIFEST_PATH_COLUMN_NAME,
    MODEL_EVALUATION_STATUS_COLUMN_NAME,
    MODEL_EVALUATION_PCT_FAVORED_COLUMN_NAME,
    MODEL_EVALUATION_PCT_OUTLIER_COLUMN_NAME,
    MODEL_EVALUATION_CLASHSCORE_COLUMN_NAME,
    MODEL_EVALUATION_OVERALL_QUALITY_COLUMN_NAME,
    MODEL_EVALUATION_RAMA_PLOT_PATH_COLUMN_NAME,
    MODEL_EVALUATION_MANIFEST_PATH_COLUMN_NAME,
    REPORT_STATUS_COLUMN_NAME,
    REPORT_PDF_PATH_COLUMN_NAME,
    REPORT_FIGURE_PNG_PATH_COLUMN_NAME,
]

STEP_STATUS_COLUMN_NAME_LIST = [
    PDB_SYNC_DONE_COLUMN_NAME,
    FASTA_FILES_DONE_COLUMN_NAME,
    SEQUENCE_ALIGNMENT_DONE_COLUMN_NAME,
    INSERTION_CODES_DONE_COLUMN_NAME,
    FILLER_STATUS_COLUMN_NAME,
    PROTONATION_STATUS_COLUMN_NAME,
    INTERNAL_CAPPING_STATUS_COLUMN_NAME,
    PREPARED_STRUCTURE_STATUS_COLUMN_NAME,
    METALS_CHECK_STATUS_COLUMN_NAME,
    METALL_PARAMS_STATUS_COLUMN_NAME,
    NONSTD_RESIDUE_PARAMS_STATUS_COLUMN_NAME,
    MODEL_EVALUATION_STATUS_COLUMN_NAME,
    REPORT_STATUS_COLUMN_NAME,
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
    protein_record[PDB_ID_COLUMN_NAME] = str(pdb_id).strip().upper()
    protein_record[RANGE_COLUMN_NAME] = str(residue_range).strip()
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