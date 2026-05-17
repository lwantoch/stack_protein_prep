"""Professional XLSX export for FRUTON pipeline state.

The exporter writes the flat FRUTON pipeline records into a compact multi-sheet
workbook. The first worksheet is a decision dashboard: it starts with a compact
logo header, a small KPI strip, and then shows only the columns needed to decide
whether a final model exists, which chemistry flags are relevant, which model
family is available, and which stage failed. Detailed paths and full-state data remain available in later
worksheets for debugging, but they do not dominate the first view.

The visual style uses an IDIS-aligned FRUTON palette rather than default office colors.
Navy (#2F3761) is used for structure, headers, and final-model text accents. Gold (#8A730F) marks review/warning states. A pre-blended soft navy
(#D5D7DF) approximates #2F3761 at 80% transparency over white and is used as the
main table body shade. Gold is used only for the accent row and review signals.
Red (#C00000) is reserved for hard failures and PDB IDs where no final model was
created.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

from openpyxl import Workbook
from openpyxl.drawing.image import Image as OpenpyxlImage
from openpyxl.drawing.spreadsheet_drawing import AnchorMarker, OneCellAnchor
from openpyxl.drawing.xdr import XDRPositiveSize2D
from openpyxl.utils.units import inch_to_EMU
from openpyxl.styles import Alignment, Border, Font, PatternFill, Side
from openpyxl.utils import get_column_letter

from stack_protein_preparation.pipeline_state import (
    ALIGNMENT_DIRECTORY_COLUMN_NAME,
    AMBER_INPUT_PATH_COLUMN_NAME,
    AMBER_OUTPUT_PATH_COLUMN_NAME,
    AMBER_RENAMING_STATUS_COLUMN_NAME,
    AMBER_TERMINI_INPUT_PATH_COLUMN_NAME,
    AMBER_TERMINI_OUTPUT_PATH_COLUMN_NAME,
    AMBER_TERMINI_STATUS_COLUMN_NAME,
    AVAILABLE_MODELS_COLUMN_NAME,
    COMPONENTS_DIRECTORY_COLUMN_NAME,
    FASTA_DIRECTORY_COLUMN_NAME,
    FASTA_FILES_DONE_COLUMN_NAME,
    FILLER_DIRECTORY_COLUMN_NAME,
    FILLER_MODEL_PATH_COLUMN_NAME,
    FILLER_MODEL_SOURCE_COLUMN_NAME,
    FILLER_STATUS_COLUMN_NAME,
    GAP_SIZES_COLUMN_NAME,
    HAS_GAPS_COLUMN_NAME,
    HAS_LIGANDS_COLUMN_NAME,
    HAS_METALS_COLUMN_NAME,
    METALS_ION_TYPE_COLUMN_NAME,
    METALS_CLASS_COLUMN_NAME,
    METALS_MODEL_READY_COLUMN_NAME,
    HAS_NONSTANDARD_RESIDUES_COLUMN_NAME,
    INSERTION_CODES_DONE_COLUMN_NAME,
    INTERNAL_CAPPING_INPUT_PATH_COLUMN_NAME,
    INTERNAL_CAPPING_OUTPUT_PATH_COLUMN_NAME,
    INTERNAL_CAPPING_STATUS_COLUMN_NAME,
    METALL_PARAMS_DIRECTORY_COLUMN_NAME,
    METALL_PARAMS_MANIFEST_PATH_COLUMN_NAME,
    METALL_PARAMS_SITE_COUNT_COLUMN_NAME,
    METALL_PARAMS_STATUS_COLUMN_NAME,
    NONSTD_RESIDUE_PARAMS_MANIFEST_PATH_COLUMN_NAME,
    NONSTD_RESIDUE_PARAMS_N_RESIDUES_COLUMN_NAME,
    NONSTD_RESIDUE_PARAMS_STATUS_COLUMN_NAME,
    N_GAPS_COLUMN_NAME,
    PDB_DIRECTORY_COLUMN_NAME,
    PDB_ID_COLUMN_NAME,
    PDB_SYNC_DONE_COLUMN_NAME,
    PREPARED_ALPHAFOLD_OUTPUT_PATH_COLUMN_NAME,
    PREPARED_DIRECTORY_COLUMN_NAME,
    PREPARED_GAPS_OUTPUT_PATH_COLUMN_NAME,
    PREPARED_MODELLER_OUTPUT_PATH_COLUMN_NAME,
    PREPARED_STRUCTURE_OUTPUT_PATH_COLUMN_NAME,
    PREPARED_STRUCTURE_PROTEIN_INPUT_PATH_COLUMN_NAME,
    PREPARED_STRUCTURE_STATUS_COLUMN_NAME,
    PREPARED_STRUCTURE_VARIANT_COLUMN_NAME,
    PROTONATION_INPUT_PATH_COLUMN_NAME,
    PROTONATION_INPUT_SOURCE_COLUMN_NAME,
    PROTONATION_OUTPUT_PATH_COLUMN_NAME,
    PROTONATION_STATUS_COLUMN_NAME,
    RANGE_COLUMN_NAME,
    RETAIN_ALPHAFOLD_VARIANT_COLUMN_NAME,
    RETAIN_GAPS_VARIANT_COLUMN_NAME,
    RETAIN_MODELLER_VARIANT_COLUMN_NAME,
    SEQUENCE_ALIGNMENT_DONE_COLUMN_NAME,
    STATE_COLUMN_NAME_LIST,
    STATUS_FAILED,
    STATUS_REQUIRED,
    STATUS_SKIPPED,
    STATUS_SUCCESS,
    STATUS_WARNING,
    STEP_STATUS_COLUMN_NAME_LIST,
    UNIPROT_ID_COLUMN_NAME,
    VARIANT_POLICY_COLUMN_NAME,
    create_empty_protein_record,
)

# ---------------------------------------------------------------------------
# XLSX-only derived columns
# ---------------------------------------------------------------------------

FINAL_MODEL_CREATED_COLUMN_NAME = "final_model"
FINAL_MODEL_TYPE_COLUMN_NAME = "final_model_type"
FINAL_MODEL_PATH_COLUMN_NAME = "final_model_path"

PARAMETER_AUDIT_STATUS_COLUMN_NAME = "parameter_audit.status"
PARAMETER_AUDIT_REQUIRES_REPAIR_COLUMN_NAME = "parameter_audit.requires_repair"
PARAMETER_AUDIT_REQUIRES_QM_COLUMN_NAME = "parameter_audit.requires_qm_parameters"
PARAMETER_AUDIT_REQUIRES_METAL_COLUMN_NAME = "parameter_audit.requires_metal_parameters"
PARAMETER_AUDIT_JSON_PATH_COLUMN_NAME = "parameter_audit.json_path"
PARAMETER_AUDIT_LOG_PATH_COLUMN_NAME = "parameter_audit.log_path"

# Metal check / geometry summary columns.
# These names are kept local in the XLSX exporter so the workbook can render
# metal diagnostics as soon as FRUTON writes them, even if older JSON states do
# not yet carry all fields. They match the names used by ``metalls_check.py``
# and the updated ``pipeline_state.py``.
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

DERIVED_COLUMN_NAME_LIST = [
    FINAL_MODEL_CREATED_COLUMN_NAME,
    FINAL_MODEL_TYPE_COLUMN_NAME,
    FINAL_MODEL_PATH_COLUMN_NAME,
]

OPTIONAL_KNOWN_EXTRA_COLUMN_NAME_LIST = [
    PARAMETER_AUDIT_STATUS_COLUMN_NAME,
    PARAMETER_AUDIT_REQUIRES_REPAIR_COLUMN_NAME,
    PARAMETER_AUDIT_REQUIRES_QM_COLUMN_NAME,
    PARAMETER_AUDIT_REQUIRES_METAL_COLUMN_NAME,
    PARAMETER_AUDIT_JSON_PATH_COLUMN_NAME,
    PARAMETER_AUDIT_LOG_PATH_COLUMN_NAME,
    METALS_CHECK_STATUS_COLUMN_NAME,
    METALS_PARAMETER_REFERENCE_COLUMN_NAME,
    METALS_GEOMETRY_FOUND_COLUMN_NAME,
    METALS_GEOMETRY_PROBABLE_COLUMN_NAME,
    METALS_GEOMETRY_MATCH_COLUMN_NAME,
    METALS_CHECK_LOG_PATH_COLUMN_NAME,
    METALS_CHECK_MANIFEST_PATH_COLUMN_NAME,
]

LOGO_FILENAME_CANDIDATE_LIST = [
    "logo_xlsx.png",
    "logo.png",
    "fruton_logo.png",
    "FRUTON_logo.png",
]

# LibreOffice position/size equivalent for the overview logo.
# The values match the manually tuned screenshot: X=0.04 in, Y=0.00 in,
# W=0.74 in, H=0.84 in. OpenXML stores drawing positions in EMUs.
LOGO_OFFSET_X_IN = 0.04
LOGO_OFFSET_Y_IN = 0.00
LOGO_WIDTH_IN = 0.75
LOGO_HEIGHT_IN = 0.80


# ---------------------------------------------------------------------------
# Palette and style constants
# ---------------------------------------------------------------------------

FONT_NAME = "Liberation Sans"
FONT_SIZE = 10
HEADER_FONT_SIZE = 10
TITLE_FONT_SIZE = 14

NAVY_HEX = "2F3761"
GOLD_HEX = "8A730F"
SOFT_NAVY_HEX = "D5D7DF"
RED_HEX = "C00000"
WHITE_HEX = "FFFFFF"

GREEN_HEX = "70AD47"
SOFT_GREEN_HEX = "D9EAD3"
SOFT_GOLD_HEX = "C5B987"
SOFT_RED_HEX = "F4CCCC"
ORANGE_HEX = "F4B183"
GREY_HEX = "E7E6E6"
STRIPE_HEX = "E4E6ED"
GRID_HEX = "B7B7B7"
LIGHT_GRID_HEX = "D9D9D9"

TITLE_FILL = PatternFill(start_color=NAVY_HEX, end_color=NAVY_HEX, fill_type="solid")
HEADER_FILL = PatternFill(start_color=NAVY_HEX, end_color=NAVY_HEX, fill_type="solid")
KPI_LABEL_FILL = PatternFill(start_color=SOFT_NAVY_HEX, end_color=SOFT_NAVY_HEX, fill_type="solid")
KPI_VALUE_FILL = PatternFill(start_color=WHITE_HEX, end_color=WHITE_HEX, fill_type="solid")
BODY_FILL = PatternFill(start_color=SOFT_NAVY_HEX, end_color=SOFT_NAVY_HEX, fill_type="solid")
ROW_STRIPE_FILL = PatternFill(start_color=STRIPE_HEX, end_color=STRIPE_HEX, fill_type="solid")

FINAL_YES_FILL = PatternFill(start_color=WHITE_HEX, end_color=WHITE_HEX, fill_type="solid")
FINAL_NO_FILL = PatternFill(start_color=RED_HEX, end_color=RED_HEX, fill_type="solid")
SUCCESS_FILL = PatternFill(start_color="E8E3CF", end_color="E8E3CF", fill_type="solid")
REVIEW_FILL = PatternFill(start_color=SOFT_GOLD_HEX, end_color=SOFT_GOLD_HEX, fill_type="solid")
STRONG_REVIEW_FILL = PatternFill(start_color=GOLD_HEX, end_color=GOLD_HEX, fill_type="solid")
FAIL_FILL = PatternFill(start_color=SOFT_RED_HEX, end_color=SOFT_RED_HEX, fill_type="solid")
SOFT_NAVY_FILL = PatternFill(start_color=SOFT_NAVY_HEX, end_color=SOFT_NAVY_HEX, fill_type="solid")
SKIP_FILL = PatternFill(start_color=SOFT_NAVY_HEX, end_color=SOFT_NAVY_HEX, fill_type="solid")
NEUTRAL_FILL = PatternFill(start_color=GREY_HEX, end_color=GREY_HEX, fill_type="solid")
ORANGE_FILL = PatternFill(start_color=ORANGE_HEX, end_color=ORANGE_HEX, fill_type="solid")

MODEL_GROUP_FILL = PatternFill(start_color="F4F6FB", end_color="F4F6FB", fill_type="solid")
CHEMISTRY_GROUP_FILL = PatternFill(start_color="FBF8EA", end_color="FBF8EA", fill_type="solid")
GAP_GROUP_FILL = PatternFill(start_color="FBF8EA", end_color="FBF8EA", fill_type="solid")
STATUS_GROUP_FILL = PatternFill(start_color="F0F2F7", end_color="F0F2F7", fill_type="solid")

TITLE_FONT = Font(name=FONT_NAME, size=TITLE_FONT_SIZE, bold=True, color=WHITE_HEX)
HEADER_FONT = Font(name=FONT_NAME, size=HEADER_FONT_SIZE, bold=True, color=WHITE_HEX)
KPI_LABEL_FONT = Font(name=FONT_NAME, size=FONT_SIZE, bold=True, color=NAVY_HEX)
KPI_VALUE_FONT = Font(name=FONT_NAME, size=FONT_SIZE, bold=True, color="111827")
PDB_FONT = Font(name=FONT_NAME, size=FONT_SIZE, bold=True, color=NAVY_HEX)
PDB_MISSING_FONT = Font(name=FONT_NAME, size=FONT_SIZE, bold=True, color=RED_HEX)
FINAL_MODEL_YES_FONT = Font(name=FONT_NAME, size=FONT_SIZE, bold=True, color=NAVY_HEX)
DATA_FONT = Font(name=FONT_NAME, size=FONT_SIZE, color="111827")
MUTED_FONT = Font(name=FONT_NAME, size=FONT_SIZE, color="666666")
WHITE_BOLD_FONT = Font(name=FONT_NAME, size=FONT_SIZE, bold=True, color=WHITE_HEX)
LINK_FONT = Font(name=FONT_NAME, size=FONT_SIZE, color=NAVY_HEX, underline="single")

TITLE_ALIGNMENT = Alignment(horizontal="left", vertical="center")
HEADER_ALIGNMENT = Alignment(horizontal="center", vertical="center", wrap_text=True)
DATA_ALIGNMENT = Alignment(vertical="top", wrap_text=True)
CENTER_ALIGNMENT = Alignment(horizontal="center", vertical="center", wrap_text=True)

GRID_SIDE = Side(style="thin", color=GRID_HEX)
LIGHT_GRID_SIDE = Side(style="thin", color=LIGHT_GRID_HEX)
NAVY_SIDE = Side(style="thin", color=NAVY_HEX)
RED_SIDE = Side(style="medium", color=RED_HEX)

GRID_BORDER = Border(
    left=GRID_SIDE,
    right=GRID_SIDE,
    top=LIGHT_GRID_SIDE,
    bottom=LIGHT_GRID_SIDE,
)

HEADER_BORDER = Border(
    left=NAVY_SIDE,
    right=NAVY_SIDE,
    top=NAVY_SIDE,
    bottom=NAVY_SIDE,
)

MISSING_ROW_BORDER = Border(
    left=RED_SIDE,
    right=RED_SIDE,
    top=RED_SIDE,
    bottom=RED_SIDE,
)

# ---------------------------------------------------------------------------
# Column groups
# ---------------------------------------------------------------------------

HYPERLINK_COLUMN_NAME_SET = {
    PDB_DIRECTORY_COLUMN_NAME,
    FASTA_DIRECTORY_COLUMN_NAME,
    ALIGNMENT_DIRECTORY_COLUMN_NAME,
    COMPONENTS_DIRECTORY_COLUMN_NAME,
    FILLER_DIRECTORY_COLUMN_NAME,
    PREPARED_DIRECTORY_COLUMN_NAME,
    METALL_PARAMS_DIRECTORY_COLUMN_NAME,
    FILLER_MODEL_PATH_COLUMN_NAME,
    PROTONATION_INPUT_PATH_COLUMN_NAME,
    PROTONATION_OUTPUT_PATH_COLUMN_NAME,
    AMBER_INPUT_PATH_COLUMN_NAME,
    AMBER_OUTPUT_PATH_COLUMN_NAME,
    AMBER_TERMINI_INPUT_PATH_COLUMN_NAME,
    AMBER_TERMINI_OUTPUT_PATH_COLUMN_NAME,
    INTERNAL_CAPPING_INPUT_PATH_COLUMN_NAME,
    INTERNAL_CAPPING_OUTPUT_PATH_COLUMN_NAME,
    PREPARED_STRUCTURE_PROTEIN_INPUT_PATH_COLUMN_NAME,
    PREPARED_STRUCTURE_OUTPUT_PATH_COLUMN_NAME,
    PREPARED_GAPS_OUTPUT_PATH_COLUMN_NAME,
    PREPARED_MODELLER_OUTPUT_PATH_COLUMN_NAME,
    PREPARED_ALPHAFOLD_OUTPUT_PATH_COLUMN_NAME,
    METALL_PARAMS_SITE_COUNT_COLUMN_NAME,
    METALL_PARAMS_MANIFEST_PATH_COLUMN_NAME,
    NONSTD_RESIDUE_PARAMS_MANIFEST_PATH_COLUMN_NAME,
    FINAL_MODEL_PATH_COLUMN_NAME,
    PARAMETER_AUDIT_JSON_PATH_COLUMN_NAME,
    PARAMETER_AUDIT_LOG_PATH_COLUMN_NAME,
    METALS_CHECK_LOG_PATH_COLUMN_NAME,
    METALS_CHECK_MANIFEST_PATH_COLUMN_NAME,
}

STATUS_COLOR_COLUMN_NAME_SET = set(STEP_STATUS_COLUMN_NAME_LIST) | {
    PARAMETER_AUDIT_STATUS_COLUMN_NAME,
    METALS_CHECK_STATUS_COLUMN_NAME,
}

SEMANTIC_COLOR_COLUMN_NAME_SET = {
    FINAL_MODEL_CREATED_COLUMN_NAME,
    HAS_METALS_COLUMN_NAME,
    HAS_LIGANDS_COLUMN_NAME,
    HAS_NONSTANDARD_RESIDUES_COLUMN_NAME,
    HAS_GAPS_COLUMN_NAME,
    RETAIN_GAPS_VARIANT_COLUMN_NAME,
    RETAIN_MODELLER_VARIANT_COLUMN_NAME,
    RETAIN_ALPHAFOLD_VARIANT_COLUMN_NAME,
    PARAMETER_AUDIT_REQUIRES_REPAIR_COLUMN_NAME,
    PARAMETER_AUDIT_REQUIRES_QM_COLUMN_NAME,
    PARAMETER_AUDIT_REQUIRES_METAL_COLUMN_NAME,
    METALS_GEOMETRY_MATCH_COLUMN_NAME,
}

MODEL_GROUP_COLUMN_NAME_SET = {
    FINAL_MODEL_TYPE_COLUMN_NAME,
    AVAILABLE_MODELS_COLUMN_NAME,
}

CHEMISTRY_GROUP_COLUMN_NAME_SET = {
    HAS_METALS_COLUMN_NAME,
    HAS_NONSTANDARD_RESIDUES_COLUMN_NAME,
    HAS_LIGANDS_COLUMN_NAME,
    METALS_PARAMETER_REFERENCE_COLUMN_NAME,
    METALS_GEOMETRY_FOUND_COLUMN_NAME,
    METALS_GEOMETRY_PROBABLE_COLUMN_NAME,
    METALS_GEOMETRY_MATCH_COLUMN_NAME,
}

GAP_GROUP_COLUMN_NAME_SET = {
    HAS_GAPS_COLUMN_NAME,
    N_GAPS_COLUMN_NAME,
    GAP_SIZES_COLUMN_NAME,
}

STATUS_GROUP_COLUMN_NAME_SET = {
    FILLER_STATUS_COLUMN_NAME,
    PREPARED_STRUCTURE_STATUS_COLUMN_NAME,
    PREPARED_STRUCTURE_VARIANT_COLUMN_NAME,
}

OVERVIEW_COLUMN_NAME_LIST = [
    PDB_ID_COLUMN_NAME,
    UNIPROT_ID_COLUMN_NAME,
    RANGE_COLUMN_NAME,
    FINAL_MODEL_CREATED_COLUMN_NAME,
    FINAL_MODEL_TYPE_COLUMN_NAME,
    AVAILABLE_MODELS_COLUMN_NAME,
    HAS_METALS_COLUMN_NAME,
    METALL_PARAMS_STATUS_COLUMN_NAME,
    HAS_NONSTANDARD_RESIDUES_COLUMN_NAME,
    NONSTD_RESIDUE_PARAMS_STATUS_COLUMN_NAME,
    NONSTD_RESIDUE_PARAMS_N_RESIDUES_COLUMN_NAME,
    HAS_LIGANDS_COLUMN_NAME,
    HAS_GAPS_COLUMN_NAME,
    N_GAPS_COLUMN_NAME,
    GAP_SIZES_COLUMN_NAME,
    PARAMETER_AUDIT_STATUS_COLUMN_NAME,
    PARAMETER_AUDIT_REQUIRES_QM_COLUMN_NAME,
    FILLER_STATUS_COLUMN_NAME,
    PREPARED_STRUCTURE_STATUS_COLUMN_NAME,
    PREPARED_STRUCTURE_VARIANT_COLUMN_NAME,
]

PREPARATION_COLUMN_NAME_LIST = [
    PDB_ID_COLUMN_NAME,
    PDB_SYNC_DONE_COLUMN_NAME,
    FASTA_FILES_DONE_COLUMN_NAME,
    SEQUENCE_ALIGNMENT_DONE_COLUMN_NAME,
    INSERTION_CODES_DONE_COLUMN_NAME,
    PARAMETER_AUDIT_STATUS_COLUMN_NAME,
    PARAMETER_AUDIT_REQUIRES_REPAIR_COLUMN_NAME,
    PARAMETER_AUDIT_REQUIRES_QM_COLUMN_NAME,
    PARAMETER_AUDIT_REQUIRES_METAL_COLUMN_NAME,
    FILLER_STATUS_COLUMN_NAME,
    FILLER_MODEL_SOURCE_COLUMN_NAME,
]

COMPONENTS_COLUMN_NAME_LIST = [
    PDB_ID_COLUMN_NAME,
    HAS_METALS_COLUMN_NAME,
    METALS_CHECK_STATUS_COLUMN_NAME,
    METALS_ION_TYPE_COLUMN_NAME,
    METALS_CLASS_COLUMN_NAME,
    METALS_MODEL_READY_COLUMN_NAME,
    METALS_PARAMETER_REFERENCE_COLUMN_NAME,
    HAS_LIGANDS_COLUMN_NAME,
    HAS_NONSTANDARD_RESIDUES_COLUMN_NAME,
    NONSTD_RESIDUE_PARAMS_STATUS_COLUMN_NAME,
    NONSTD_RESIDUE_PARAMS_N_RESIDUES_COLUMN_NAME,
]

GAPS_AND_VARIANTS_COLUMN_NAME_LIST = [
    PDB_ID_COLUMN_NAME,
    N_GAPS_COLUMN_NAME,
    GAP_SIZES_COLUMN_NAME,
    HAS_GAPS_COLUMN_NAME,
    FILLER_STATUS_COLUMN_NAME,
    FILLER_MODEL_SOURCE_COLUMN_NAME,
    VARIANT_POLICY_COLUMN_NAME,
    RETAIN_GAPS_VARIANT_COLUMN_NAME,
    RETAIN_MODELLER_VARIANT_COLUMN_NAME,
    RETAIN_ALPHAFOLD_VARIANT_COLUMN_NAME,
    PREPARED_STRUCTURE_VARIANT_COLUMN_NAME,
    AVAILABLE_MODELS_COLUMN_NAME,
]

CHEMISTRY_COLUMN_NAME_LIST = [
    PDB_ID_COLUMN_NAME,
    PROTONATION_STATUS_COLUMN_NAME,
    PROTONATION_INPUT_SOURCE_COLUMN_NAME,
    INTERNAL_CAPPING_STATUS_COLUMN_NAME,
    METALS_CHECK_STATUS_COLUMN_NAME,
    METALS_ION_TYPE_COLUMN_NAME,
    METALS_CLASS_COLUMN_NAME,
    METALS_MODEL_READY_COLUMN_NAME,
    METALS_PARAMETER_REFERENCE_COLUMN_NAME,
    METALS_GEOMETRY_FOUND_COLUMN_NAME,
    METALS_GEOMETRY_PROBABLE_COLUMN_NAME,
    METALS_GEOMETRY_MATCH_COLUMN_NAME,
    METALL_PARAMS_STATUS_COLUMN_NAME,
    NONSTD_RESIDUE_PARAMS_STATUS_COLUMN_NAME,
    NONSTD_RESIDUE_PARAMS_N_RESIDUES_COLUMN_NAME,
    PARAMETER_AUDIT_STATUS_COLUMN_NAME,
    PARAMETER_AUDIT_REQUIRES_REPAIR_COLUMN_NAME,
    PARAMETER_AUDIT_REQUIRES_QM_COLUMN_NAME,
    PARAMETER_AUDIT_REQUIRES_METAL_COLUMN_NAME,
]

PREPARED_OUTPUTS_COLUMN_NAME_LIST = [
    PDB_ID_COLUMN_NAME,
    FINAL_MODEL_CREATED_COLUMN_NAME,
    FINAL_MODEL_TYPE_COLUMN_NAME,
    AVAILABLE_MODELS_COLUMN_NAME,
    METALS_MODEL_READY_COLUMN_NAME,
    METALS_CHECK_STATUS_COLUMN_NAME,
    PREPARED_STRUCTURE_STATUS_COLUMN_NAME,
    PREPARED_STRUCTURE_VARIANT_COLUMN_NAME,
    FINAL_MODEL_PATH_COLUMN_NAME,
    PREPARED_STRUCTURE_OUTPUT_PATH_COLUMN_NAME,
    PREPARED_GAPS_OUTPUT_PATH_COLUMN_NAME,
    PREPARED_MODELLER_OUTPUT_PATH_COLUMN_NAME,
    PREPARED_ALPHAFOLD_OUTPUT_PATH_COLUMN_NAME,
]

PATHS_COLUMN_NAME_LIST = [
    PDB_ID_COLUMN_NAME,
    PDB_DIRECTORY_COLUMN_NAME,
    COMPONENTS_DIRECTORY_COLUMN_NAME,
    ALIGNMENT_DIRECTORY_COLUMN_NAME,
    FILLER_DIRECTORY_COLUMN_NAME,
    PREPARED_DIRECTORY_COLUMN_NAME,
    FINAL_MODEL_PATH_COLUMN_NAME,
    FILLER_MODEL_PATH_COLUMN_NAME,
    PARAMETER_AUDIT_JSON_PATH_COLUMN_NAME,
    PARAMETER_AUDIT_LOG_PATH_COLUMN_NAME,
    METALS_CHECK_LOG_PATH_COLUMN_NAME,
    METALS_CHECK_MANIFEST_PATH_COLUMN_NAME,
    NONSTD_RESIDUE_PARAMS_MANIFEST_PATH_COLUMN_NAME,
    PROTONATION_INPUT_PATH_COLUMN_NAME,
    PROTONATION_OUTPUT_PATH_COLUMN_NAME,
    INTERNAL_CAPPING_INPUT_PATH_COLUMN_NAME,
    INTERNAL_CAPPING_OUTPUT_PATH_COLUMN_NAME,
    PREPARED_STRUCTURE_PROTEIN_INPUT_PATH_COLUMN_NAME,
    PREPARED_STRUCTURE_OUTPUT_PATH_COLUMN_NAME,
    PREPARED_GAPS_OUTPUT_PATH_COLUMN_NAME,
    PREPARED_MODELLER_OUTPUT_PATH_COLUMN_NAME,
    PREPARED_ALPHAFOLD_OUTPUT_PATH_COLUMN_NAME,
]

WORKSHEET_LAYOUTS: list[tuple[str, list[str]]] = [
    ("overview", OVERVIEW_COLUMN_NAME_LIST),
    ("preparation", PREPARATION_COLUMN_NAME_LIST),
    ("components", COMPONENTS_COLUMN_NAME_LIST),
    ("gaps_and_variants", GAPS_AND_VARIANTS_COLUMN_NAME_LIST),
    ("chemistry", CHEMISTRY_COLUMN_NAME_LIST),
    ("prepared_outputs", PREPARED_OUTPUTS_COLUMN_NAME_LIST),
    ("paths", PATHS_COLUMN_NAME_LIST),
    ("full_state", []),
]


# ---------------------------------------------------------------------------
# Record normalization and derived values
# ---------------------------------------------------------------------------


def get_excel_column_order(
    protein_record_list: list[dict[str, Any]] | None = None,
) -> list[str]:
    """Return the full export column order including known extra fields.

    The JSON state is defined centrally, but some recent FRUTON stages may add
    extra diagnostic keys before ``pipeline_state.py`` is extended. The XLSX
    exporter keeps those values visible in ``full_state`` instead of discarding
    them. Derived display columns are appended first so they can be reused by
    compact worksheets without changing the persisted JSON schema.
    """

    column_order = list(STATE_COLUMN_NAME_LIST)

    for column_name in OPTIONAL_KNOWN_EXTRA_COLUMN_NAME_LIST:
        if column_name not in column_order:
            column_order.append(column_name)

    if protein_record_list:
        extra_column_names = sorted(
            {
                str(column_name)
                for protein_record in protein_record_list
                for column_name in protein_record
                if str(column_name) not in column_order
            }
        )
        column_order.extend(extra_column_names)

    for column_name in DERIVED_COLUMN_NAME_LIST:
        if column_name not in column_order:
            column_order.append(column_name)

    return column_order


def ensure_all_records_have_all_columns(
    protein_record_list: list[dict[str, Any]],
) -> None:
    """Ensure that every record has at least the core pipeline columns."""

    empty_record = create_empty_protein_record()

    for protein_record in protein_record_list:
        for column_name in empty_record:
            if column_name not in protein_record:
                protein_record[column_name] = ""


def normalize_record_values(
    protein_record_list: list[dict[str, Any]],
) -> list[dict[str, str]]:
    """Return normalized copies of protein records for export.

    The function preserves known state columns, optional diagnostic columns, and
    any extra keys present in the incoming record dictionaries. It then enriches
    each record with XLSX-only derived fields that summarize whether a prepared
    final model exists and which path should be treated as the main output. This
    keeps the workbook useful even when the JSON state evolves faster than the
    spreadsheet module.
    """

    normalized_record_list: list[dict[str, str]] = []
    column_order = get_excel_column_order(protein_record_list)

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

        _add_derived_final_model_fields(normalized_record)
        normalized_record_list.append(normalized_record)

    return normalized_record_list


def create_unique_sorted_record_list(
    protein_record_list: list[dict[str, str]],
) -> list[dict[str, str]]:
    """Return unique records sorted by PDB ID."""

    record_by_pdb_id: dict[str, dict[str, str]] = {}

    for protein_record in protein_record_list:
        pdb_id = str(protein_record.get(PDB_ID_COLUMN_NAME, "")).strip().upper()
        if not pdb_id:
            continue

        copied_record = dict(protein_record)
        copied_record[PDB_ID_COLUMN_NAME] = pdb_id
        record_by_pdb_id[pdb_id] = copied_record

    return [record_by_pdb_id[pdb_id] for pdb_id in sorted(record_by_pdb_id.keys())]


def _add_derived_final_model_fields(protein_record: dict[str, str]) -> None:
    """Add spreadsheet-only final-model summary fields to one record."""

    final_path = _choose_final_model_path(protein_record)
    final_type = _choose_final_model_type(protein_record, final_path)

    protein_record[FINAL_MODEL_CREATED_COLUMN_NAME] = "yes" if final_path else "no"
    protein_record[FINAL_MODEL_TYPE_COLUMN_NAME] = final_type
    protein_record[FINAL_MODEL_PATH_COLUMN_NAME] = final_path


def _choose_final_model_path(protein_record: dict[str, str]) -> str:
    """Return the preferred final prepared model path, or an empty string."""

    for column_name in (
        PREPARED_STRUCTURE_OUTPUT_PATH_COLUMN_NAME,
        PREPARED_MODELLER_OUTPUT_PATH_COLUMN_NAME,
        PREPARED_ALPHAFOLD_OUTPUT_PATH_COLUMN_NAME,
        PREPARED_GAPS_OUTPUT_PATH_COLUMN_NAME,
    ):
        value = str(protein_record.get(column_name, "")).strip()
        if value:
            return value

    return ""


def _choose_final_model_type(
    protein_record: dict[str, str],
    final_path: str,
) -> str:
    """Return a compact label for the final model type."""

    if not final_path:
        return ""

    prepared_variant = str(
        protein_record.get(PREPARED_STRUCTURE_VARIANT_COLUMN_NAME, "")
    ).strip()
    if prepared_variant:
        return prepared_variant

    available_models = str(protein_record.get(AVAILABLE_MODELS_COLUMN_NAME, "")).strip()
    if available_models:
        return available_models

    return "prepared"


def _record_has_final_model(protein_record: dict[str, str]) -> bool:
    return str(protein_record.get(FINAL_MODEL_CREATED_COLUMN_NAME, "")).lower() == "yes"


# ---------------------------------------------------------------------------
# Column filtering
# ---------------------------------------------------------------------------


def get_nonempty_column_order_for_subset(
    protein_record_list: list[dict[str, str]],
    base_column_order: list[str],
    always_keep_column_name_set: set[str] | None = None,
) -> list[str]:
    """Return only useful columns from one requested worksheet layout."""

    always_keep_column_name_set = always_keep_column_name_set or set()
    filtered_column_order: list[str] = []

    for column_name in base_column_order:
        if column_name in always_keep_column_name_set:
            filtered_column_order.append(column_name)
            continue

        has_nonempty_value = any(
            str(protein_record.get(column_name, "")).strip()
            for protein_record in protein_record_list
        )

        if has_nonempty_value:
            filtered_column_order.append(column_name)

    return filtered_column_order



# ---------------------------------------------------------------------------
# Logo helpers
# ---------------------------------------------------------------------------


def _resolve_logo_path(
    *,
    output_path: Path,
    logo_path: str | Path | None,
) -> Path | None:
    """Return the first available logo path for the overview header.

    The function keeps the public writer convenient for the FRUTON runner by
    accepting an optional explicit path while also checking sensible project
    locations. The most useful default is ``data/logo.png`` at the repository
    root, because the workbook is usually written to ``data/proteins``. Missing
    logos never fail the export; the sheet falls back to a clean text-only
    header so automated pipeline runs remain robust.
    """

    candidate_path_list: list[Path] = []

    if logo_path is not None:
        candidate_path_list.append(Path(logo_path).expanduser())

    output_path = Path(output_path).expanduser()
    for filename in LOGO_FILENAME_CANDIDATE_LIST:
        candidate_path_list.extend(
            [
                output_path.parent / filename,
                output_path.parent.parent / filename,
                Path.cwd() / "data" / filename,
                Path.cwd() / filename,
                Path("/mnt/data") / filename,
            ]
        )

    for candidate_path in candidate_path_list:
        try:
            resolved_path = candidate_path.resolve()
        except Exception:
            resolved_path = candidate_path

        if resolved_path.exists() and resolved_path.is_file():
            return resolved_path

    return None


def _add_logo_to_worksheet(
    worksheet,
    *,
    logo_path: Path | None,
) -> bool:
    """Place the FRUTON logo using the manually tuned Calc geometry.

    LibreOffice reports object position and size in inches, while XLSX stores
    drawing anchors in EMUs. A one-cell anchor with explicit EMU offsets keeps
    the logo in the top-left dashboard band without depending on column widths
    or row heights. The size intentionally follows the inspected Calc values:
    X=0.04 in, Y=0.00 in, W=0.74 in, H=0.84 in.
    """

    if logo_path is None:
        return False

    try:
        logo_image = OpenpyxlImage(str(logo_path))
    except Exception:
        return False

    marker = AnchorMarker(
        col=0,
        colOff=inch_to_EMU(LOGO_OFFSET_X_IN),
        row=0,
        rowOff=inch_to_EMU(LOGO_OFFSET_Y_IN),
    )
    size = XDRPositiveSize2D(
        cx=inch_to_EMU(LOGO_WIDTH_IN),
        cy=inch_to_EMU(LOGO_HEIGHT_IN),
    )
    logo_image.anchor = OneCellAnchor(_from=marker, ext=size)
    worksheet.add_image(logo_image)
    return True

# ---------------------------------------------------------------------------
# KPI helpers
# ---------------------------------------------------------------------------


def _count_records_where(
    protein_record_list: list[dict[str, str]],
    column_name: str,
    value: str,
) -> int:
    return sum(
        1
        for protein_record in protein_record_list
        if protein_record.get(column_name, "").strip().lower() == value
    )


def _build_overview_kpi_list(
    protein_record_list: list[dict[str, str]],
) -> list[tuple[str, int]]:
    """Return compact top-row metrics for the overview worksheet."""

    total = len(protein_record_list)
    final_models = _count_records_where(
        protein_record_list,
        FINAL_MODEL_CREATED_COLUMN_NAME,
        "yes",
    )

    return [
        ("records", total),
        ("final models", final_models),
        ("missing final", total - final_models),
        (
            "metals",
            _count_records_where(protein_record_list, HAS_METALS_COLUMN_NAME, "yes"),
        ),
        (
            "metal ready",
            _count_records_where(protein_record_list, METALS_MODEL_READY_COLUMN_NAME, "yes"),
        ),
        (
            "nonstandard",
            _count_records_where(
                protein_record_list,
                HAS_NONSTANDARD_RESIDUES_COLUMN_NAME,
                "yes",
            ),
        ),
        (
            "gap cases",
            _count_records_where(protein_record_list, HAS_GAPS_COLUMN_NAME, "yes"),
        ),
    ]


# ---------------------------------------------------------------------------
# Styling helpers
# ---------------------------------------------------------------------------


def _set_white_bold(cell) -> None:
    cell.font = WHITE_BOLD_FONT
    cell.alignment = CENTER_ALIGNMENT


def apply_status_cell_color(cell, value: str) -> None:
    """Apply status coloring to one cell."""

    normalized_value = str(value).strip().lower()

    if normalized_value == STATUS_SUCCESS:
        cell.fill = SUCCESS_FILL
    elif normalized_value == STATUS_WARNING:
        cell.fill = STRONG_REVIEW_FILL
        _set_white_bold(cell)
    elif normalized_value == STATUS_REQUIRED:
        cell.fill = REVIEW_FILL
    elif normalized_value == STATUS_SKIPPED:
        cell.fill = SKIP_FILL
        cell.font = MUTED_FONT
    elif normalized_value == STATUS_FAILED:
        cell.fill = FAIL_FILL


def apply_semantic_yes_no_color(
    cell,
    column_name: str,
    value: str,
) -> None:
    """Apply value-based colors with one meaning per color.

    Normal yes/no values stay mostly unpainted. Strong fill colors are reserved
    for values that change interpretation: metals are blockers, non-standard
    residues and gaps need review, retained variants are available, and missing
    final models are marked by red text. This avoids the ambiguous situation
    where pale yellow or red appears only because a column belongs to a visual
    group rather than because the value matters.
    """

    normalized_value = str(value).strip().lower()

    if column_name == FINAL_MODEL_CREATED_COLUMN_NAME:
        cell.alignment = CENTER_ALIGNMENT
        if normalized_value == "yes":
            cell.font = FINAL_MODEL_YES_FONT
        elif normalized_value == "no":
            cell.font = PDB_MISSING_FONT
        return

    if column_name == HAS_LIGANDS_COLUMN_NAME:
        return

    if column_name == HAS_METALS_COLUMN_NAME:
        if normalized_value == "yes":
            cell.fill = REVIEW_FILL
        return

    if column_name == HAS_NONSTANDARD_RESIDUES_COLUMN_NAME:
        if normalized_value == "yes":
            cell.fill = REVIEW_FILL
        return

    if column_name == HAS_GAPS_COLUMN_NAME:
        if normalized_value == "yes":
            cell.fill = REVIEW_FILL
        return

    if column_name in {
        RETAIN_GAPS_VARIANT_COLUMN_NAME,
        RETAIN_MODELLER_VARIANT_COLUMN_NAME,
        RETAIN_ALPHAFOLD_VARIANT_COLUMN_NAME,
    }:
        if normalized_value == "yes":
            cell.fill = SUCCESS_FILL
        return

    if column_name == METALS_MODEL_READY_COLUMN_NAME:
        cell.alignment = CENTER_ALIGNMENT
        if normalized_value == "yes":
            cell.fill = SUCCESS_FILL
        elif normalized_value == "no":
            cell.fill = FAIL_FILL
        elif normalized_value in {"review", "warning", "required"}:
            cell.fill = REVIEW_FILL
        return

    if column_name == METALS_GEOMETRY_MATCH_COLUMN_NAME:
        cell.alignment = CENTER_ALIGNMENT
        if normalized_value == "yes":
            cell.font = Font(name=FONT_NAME, size=FONT_SIZE, bold=True, color=GOLD_HEX)
        elif normalized_value == "no":
            cell.font = PDB_MISSING_FONT
        return

    if column_name in {
        PARAMETER_AUDIT_REQUIRES_REPAIR_COLUMN_NAME,
        PARAMETER_AUDIT_REQUIRES_QM_COLUMN_NAME,
        PARAMETER_AUDIT_REQUIRES_METAL_COLUMN_NAME,
    }:
        if normalized_value == "true":
            cell.fill = FAIL_FILL
        return




def apply_metal_geometry_pair_color(
    worksheet,
    protein_record_list: list[dict[str, str]],
    column_order: list[str],
    *,
    first_data_row_index: int,
) -> None:
    """Color found/probable metal geometry text from the match decision.

    Transition-metal rows need two visually linked geometry cells. If the
    predicted and observed geometry agree, both geometry texts are gold. If
    they disagree, both geometry texts are red. The cells keep the shared body
    fill, so the color works as a precise diagnostic signal instead of adding
    another large block of background shading.
    """

    try:
        found_column_index = column_order.index(METALS_GEOMETRY_FOUND_COLUMN_NAME) + 1
        probable_column_index = column_order.index(METALS_GEOMETRY_PROBABLE_COLUMN_NAME) + 1
    except ValueError:
        return

    for offset, protein_record in enumerate(protein_record_list):
        geometry_match = str(
            protein_record.get(METALS_GEOMETRY_MATCH_COLUMN_NAME, "")
        ).strip().lower()

        if geometry_match == "yes":
            geometry_font = Font(name=FONT_NAME, size=FONT_SIZE, bold=True, color=GOLD_HEX)
        elif geometry_match == "no":
            geometry_font = PDB_MISSING_FONT
        else:
            continue

        row_index = first_data_row_index + offset
        for column_index in (found_column_index, probable_column_index):
            cell = worksheet.cell(row=row_index, column=column_index)
            if cell.value not in (None, ""):
                cell.font = geometry_font


def style_header_row(
    worksheet,
    column_order: list[str],
    header_row_index: int,
) -> None:
    """Apply professional header styling."""

    for column_index, _column_name in enumerate(column_order, start=1):
        cell = worksheet.cell(row=header_row_index, column=column_index)
        cell.fill = HEADER_FILL
        cell.font = HEADER_FONT
        cell.alignment = HEADER_ALIGNMENT
        cell.border = HEADER_BORDER


def apply_base_data_style(
    worksheet,
    *,
    first_data_row_index: int,
) -> None:
    """Apply the consistent body base before semantic colors.

    The overview should read as one designed table rather than as many
    independent colored blocks. Every data row therefore starts with the same
    soft navy base, with a slightly lighter stripe on alternating rows. Semantic
    colors are applied afterwards and only override cells whose value actually
    matters.
    """

    for row_index in range(first_data_row_index, worksheet.max_row + 1):
        base_fill = (
            ROW_STRIPE_FILL
            if (row_index - first_data_row_index) % 2 == 1
            else BODY_FILL
        )

        for cell in worksheet[row_index]:
            cell.font = DATA_FONT
            cell.alignment = DATA_ALIGNMENT
            cell.border = GRID_BORDER
            cell.fill = base_fill


def apply_group_body_shading(
    worksheet,
    column_order: list[str],
    *,
    first_data_row_index: int,
) -> None:
    """Keep body shading value-based instead of group-based.

    The earlier column-group background colors made the overview harder to
    interpret because pale gold could mean either a chemistry group, a gap group,
    or a real warning. The body now uses only subtle row striping plus semantic
    value colors. This gives a stricter grammar: red means missing/failure,
    gold means review/gap/non-standard, green means success, and grey means
    skipped.
    """

    return


def apply_row_striping(
    worksheet,
    *,
    first_data_row_index: int,
) -> None:
    """No-op: striping is applied in ``apply_base_data_style``.

    Keeping this function preserves the writer call order while preventing a
    later styling pass from touching semantically colored cells.
    """

    return


def mark_missing_model_pdb_ids(
    worksheet,
    protein_record_list: list[dict[str, str]],
    column_order: list[str],
    *,
    first_data_row_index: int,
) -> None:
    """Mark missing final models by coloring the PDB ID red.

    A full red row border makes the overview visually noisy because it competes
    with every status and chemistry color. The PDB ID is the natural row anchor,
    so missing final models are marked there instead. This keeps the warning
    visible even when the user scrolls horizontally, without turning the whole
    sheet into a red grid.
    """

    try:
        pdb_column_index = column_order.index(PDB_ID_COLUMN_NAME) + 1
    except ValueError:
        return

    for offset, protein_record in enumerate(protein_record_list):
        if _record_has_final_model(protein_record):
            continue

        row_index = first_data_row_index + offset
        cell = worksheet.cell(row=row_index, column=pdb_column_index)
        cell.font = PDB_MISSING_FONT
        cell.alignment = CENTER_ALIGNMENT


def autosize_worksheet_columns(worksheet) -> None:
    """Set readable bounded widths for all worksheet columns."""

    for column_index in range(1, worksheet.max_column + 1):
        column_letter = get_column_letter(column_index)
        max_length = 0

        for row_index in range(1, worksheet.max_row + 1):
            cell = worksheet.cell(row=row_index, column=column_index)
            cell_value = "" if cell.value is None else str(cell.value)
            max_length = max(max_length, len(cell_value))

        header_value = str(worksheet.cell(row=1, column=column_index).value or "")

        if header_value.endswith("_path") or "directory" in header_value:
            adjusted_width = 12
        elif column_index == 1:
            adjusted_width = 11
        else:
            adjusted_width = min(max(max_length + 2, 10), 28)

        worksheet.column_dimensions[column_letter].width = adjusted_width


def set_worksheet_view_options(
    worksheet,
    *,
    freeze_cell: str,
    auto_filter_range: str,
) -> None:
    """Apply stable worksheet usability defaults."""

    worksheet.freeze_panes = freeze_cell
    worksheet.auto_filter.ref = auto_filter_range
    worksheet.sheet_view.showGridLines = False


def _set_click_hyperlink(cell, target_path: str) -> None:
    """Store the path as hyperlink target, but show only ``open``."""

    if not target_path:
        return

    try:
        target = Path(target_path).expanduser().resolve()
        cell.value = "open"
        cell.hyperlink = target.as_uri()
        cell.font = LINK_FONT
        cell.alignment = CENTER_ALIGNMENT
    except Exception:
        cell.value = target_path


def _set_row_height(worksheet, row_index: int, height: float) -> None:
    worksheet.row_dimensions[row_index].height = height


# ---------------------------------------------------------------------------
# Worksheet writing
# ---------------------------------------------------------------------------


def _write_overview_title_and_kpis(
    worksheet,
    protein_record_list: list[dict[str, str]],
    column_count: int,
    *,
    logo_path: Path | None,
) -> None:
    """Write the overview header, gold band, and KPI strip.

    The merge pattern follows the manually adjusted Calc layout: ``A1:A2`` is
    the logo cell, ``B1:last`` is the title, ``B2:last`` is the subtitle, and
    ``A3:last`` is one continuous gold accent row. The title is ``SUMMARY``
    because the workbook summarizes a framework run, not a generic pipeline.
    """

    max_column = max(column_count, 14)
    last_column_letter = get_column_letter(max_column)

    for row_index in range(1, 3):
        for column_index in range(1, max_column + 1):
            cell = worksheet.cell(row=row_index, column=column_index)
            cell.fill = TITLE_FILL
            cell.border = HEADER_BORDER

    worksheet.merge_cells("A1:A2")
    logo_added = _add_logo_to_worksheet(worksheet, logo_path=logo_path)

    logo_cell = worksheet.cell(row=1, column=1)
    if not logo_added:
        logo_cell.value = "FRUTON"
        logo_cell.font = TITLE_FONT
        logo_cell.alignment = CENTER_ALIGNMENT

    worksheet.merge_cells(f"B1:{last_column_letter}1")
    title_cell = worksheet.cell(row=1, column=2)
    title_cell.value = "SUMMARY"
    title_cell.fill = TITLE_FILL
    title_cell.font = TITLE_FONT
    title_cell.alignment = TITLE_ALIGNMENT
    title_cell.border = HEADER_BORDER

    worksheet.merge_cells(f"B2:{last_column_letter}2")
    subtitle_cell = worksheet.cell(row=2, column=2)
    subtitle_cell.value = "model availability · gaps · chemistry · preparation status"
    subtitle_cell.fill = TITLE_FILL
    subtitle_cell.font = Font(name=FONT_NAME, size=10, color=WHITE_HEX)
    subtitle_cell.alignment = TITLE_ALIGNMENT
    subtitle_cell.border = HEADER_BORDER

    worksheet.merge_cells(f"A3:{last_column_letter}3")
    accent_cell = worksheet.cell(row=3, column=1)
    accent_cell.value = ""
    accent_cell.fill = STRONG_REVIEW_FILL
    accent_cell.border = HEADER_BORDER

    for column_index in range(1, max_column + 1):
        cell = worksheet.cell(row=4, column=column_index)
        cell.fill = KPI_VALUE_FILL
        cell.border = GRID_BORDER

    _set_row_height(worksheet, 1, 31.0)
    _set_row_height(worksheet, 2, 20.25)
    _set_row_height(worksheet, 3, 10.1)
    _set_row_height(worksheet, 4, 22)

    kpi_list = _build_overview_kpi_list(protein_record_list)
    start_column = 1

    for label, value in kpi_list:
        label_cell = worksheet.cell(row=4, column=start_column)
        value_cell = worksheet.cell(row=4, column=start_column + 1)

        label_cell.value = label
        value_cell.value = value

        label_cell.fill = KPI_LABEL_FILL
        label_cell.font = KPI_LABEL_FONT
        label_cell.alignment = CENTER_ALIGNMENT
        label_cell.border = GRID_BORDER

        value_cell.fill = KPI_VALUE_FILL
        value_cell.font = KPI_VALUE_FONT
        value_cell.alignment = CENTER_ALIGNMENT
        value_cell.border = GRID_BORDER

        if label == "missing final" and value:
            value_cell.font = PDB_MISSING_FONT

        start_column += 2


def _append_record_row(
    worksheet,
    protein_record: dict[str, str],
    column_order: list[str],
) -> None:
    row_value_list = [
        protein_record.get(column_name, "") for column_name in column_order
    ]
    worksheet.append(row_value_list)


def _style_record_rows(
    worksheet,
    protein_record_list: list[dict[str, str]],
    column_order: list[str],
    *,
    first_data_row_index: int,
) -> None:
    for offset, _protein_record in enumerate(protein_record_list):
        row_index = first_data_row_index + offset

        for column_index, column_name in enumerate(column_order, start=1):
            cell = worksheet.cell(row=row_index, column=column_index)
            cell_value = "" if cell.value is None else str(cell.value).strip()

            if column_name == PDB_ID_COLUMN_NAME:
                cell.font = PDB_FONT
                cell.alignment = CENTER_ALIGNMENT

            if column_name in STATUS_COLOR_COLUMN_NAME_SET:
                apply_status_cell_color(cell, cell_value)

            if column_name in SEMANTIC_COLOR_COLUMN_NAME_SET:
                apply_semantic_yes_no_color(
                    cell=cell,
                    column_name=column_name,
                    value=cell_value,
                )
            if column_name in {N_GAPS_COLUMN_NAME, GAP_SIZES_COLUMN_NAME}:
                if cell_value and cell_value.lower() not in {"0", "none"}:
                    cell.fill = REVIEW_FILL


            if column_name in HYPERLINK_COLUMN_NAME_SET and cell_value:
                _set_click_hyperlink(cell, cell_value)


def _write_one_worksheet(
    workbook: Workbook,
    worksheet_name: str,
    protein_record_list: list[dict[str, str]],
    base_column_order: list[str],
    *,
    always_keep_column_name_set: set[str] | None = None,
    logo_path: Path | None = None,
) -> None:
    column_order = get_nonempty_column_order_for_subset(
        protein_record_list=protein_record_list,
        base_column_order=base_column_order,
        always_keep_column_name_set=always_keep_column_name_set,
    )

    worksheet = workbook.create_sheet(title=worksheet_name)

    if worksheet_name == "overview":
        header_row_index = 6
        first_data_row_index = 7
        _write_overview_title_and_kpis(
            worksheet=worksheet,
            protein_record_list=protein_record_list,
            column_count=len(column_order),
            logo_path=logo_path,
        )
        worksheet.append([])
        worksheet.append(column_order)
    else:
        header_row_index = 1
        first_data_row_index = 2
        worksheet.append(column_order)

    style_header_row(
        worksheet=worksheet,
        column_order=column_order,
        header_row_index=header_row_index,
    )

    for protein_record in protein_record_list:
        _append_record_row(
            worksheet=worksheet,
            protein_record=protein_record,
            column_order=column_order,
        )

    apply_base_data_style(
        worksheet=worksheet,
        first_data_row_index=first_data_row_index,
    )
    apply_group_body_shading(
        worksheet=worksheet,
        column_order=column_order,
        first_data_row_index=first_data_row_index,
    )
    _style_record_rows(
        worksheet=worksheet,
        protein_record_list=protein_record_list,
        column_order=column_order,
        first_data_row_index=first_data_row_index,
    )
    apply_metal_geometry_pair_color(
        worksheet=worksheet,
        protein_record_list=protein_record_list,
        column_order=column_order,
        first_data_row_index=first_data_row_index,
    )
    apply_row_striping(
        worksheet=worksheet,
        first_data_row_index=first_data_row_index,
    )
    mark_missing_model_pdb_ids(
        worksheet=worksheet,
        protein_record_list=protein_record_list,
        column_order=column_order,
        first_data_row_index=first_data_row_index,
    )

    for row_index in range(first_data_row_index, worksheet.max_row + 1):
        _set_row_height(worksheet, row_index, 18)

    autosize_worksheet_columns(worksheet)

    last_column_letter = get_column_letter(max(1, len(column_order)))
    auto_filter_range = (
        f"A{header_row_index}:{last_column_letter}{worksheet.max_row}"
    )
    set_worksheet_view_options(
        worksheet=worksheet,
        freeze_cell=f"A{first_data_row_index}",
        auto_filter_range=auto_filter_range,
    )


# ---------------------------------------------------------------------------
# Public writer
# ---------------------------------------------------------------------------


def write_pipeline_to_xlsx(
    protein_record_list: list[dict[str, Any]],
    output_path: Path,
    logo_path: str | Path | None = None,
) -> None:
    """Write pipeline records to a polished multi-sheet XLSX file."""

    working_record_list = [
        dict(protein_record) for protein_record in protein_record_list
    ]

    ensure_all_records_have_all_columns(working_record_list)
    normalized_record_list = normalize_record_values(working_record_list)
    unique_sorted_record_list = create_unique_sorted_record_list(normalized_record_list)

    output_path = Path(output_path)
    resolved_logo_path = _resolve_logo_path(
        output_path=output_path,
        logo_path=logo_path,
    )

    workbook = Workbook()

    default_sheet = workbook.active
    workbook.remove(default_sheet)

    for worksheet_name, base_column_order in WORKSHEET_LAYOUTS:
        if worksheet_name == "full_state":
            base_column_order = get_excel_column_order(unique_sorted_record_list)

        always_keep = {PDB_ID_COLUMN_NAME}

        if worksheet_name == "overview":
            always_keep |= {
                UNIPROT_ID_COLUMN_NAME,
                RANGE_COLUMN_NAME,
                FINAL_MODEL_CREATED_COLUMN_NAME,
                FINAL_MODEL_TYPE_COLUMN_NAME,
                AVAILABLE_MODELS_COLUMN_NAME,
            }

        if worksheet_name in {"gaps_and_variants", "prepared_outputs"}:
            always_keep |= {
                FINAL_MODEL_CREATED_COLUMN_NAME,
                AVAILABLE_MODELS_COLUMN_NAME,
            }

        _write_one_worksheet(
            workbook=workbook,
            worksheet_name=worksheet_name,
            protein_record_list=unique_sorted_record_list,
            base_column_order=base_column_order,
            always_keep_column_name_set=always_keep,
            logo_path=resolved_logo_path if worksheet_name == "overview" else None,
        )

    output_path.parent.mkdir(parents=True, exist_ok=True)
    workbook.save(output_path)