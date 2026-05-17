"""Pure record-management helpers for the FRUTON pipeline.

All functions in this module are pure: they only depend on pipeline_state
constants and the pipeline_table helper.  No fruton globals are referenced
here.
"""

from __future__ import annotations

from pathlib import Path

from stack_protein_preparation.pipeline_state import (
    ALIGNMENT_DIRECTORY_COLUMN_NAME,
    AVAILABLE_MODELS_COLUMN_NAME,
    COMPONENTS_DIRECTORY_COLUMN_NAME,
    FASTA_DIRECTORY_COLUMN_NAME,
    FILLER_DIRECTORY_COLUMN_NAME,
    FILLER_MODEL_PATH_COLUMN_NAME,
    FILLER_MODEL_SOURCE_COLUMN_NAME,
    FILLER_STATUS_COLUMN_NAME,
    HAS_GAPS_COLUMN_NAME,
    HAS_LIGANDS_COLUMN_NAME,
    HAS_METALS_COLUMN_NAME,
    HAS_NONSTANDARD_RESIDUES_COLUMN_NAME,
    INTERNAL_CAPPING_INPUT_PATH_COLUMN_NAME,
    INTERNAL_CAPPING_OUTPUT_PATH_COLUMN_NAME,
    INTERNAL_CAPPING_STATUS_COLUMN_NAME,
    METALS_CHECK_LOG_PATH_COLUMN_NAME,
    METALS_CHECK_MANIFEST_PATH_COLUMN_NAME,
    METALS_CHECK_STATUS_COLUMN_NAME,
    METALS_CLASS_COLUMN_NAME,
    METALS_GEOMETRY_FOUND_COLUMN_NAME,
    METALS_GEOMETRY_MATCH_COLUMN_NAME,
    METALS_GEOMETRY_PROBABLE_COLUMN_NAME,
    METALS_ION_TYPE_COLUMN_NAME,
    METALS_MODEL_READY_COLUMN_NAME,
    METALS_PARAMETER_REFERENCE_COLUMN_NAME,
    N_GAPS_COLUMN_NAME,
    GAP_SIZES_COLUMN_NAME,
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
    STATUS_SUCCESS,
    UNIPROT_ID_COLUMN_NAME,
    VARIANT_POLICY_COLUMN_NAME,
    create_protein_record,
)
from stack_protein_preparation.pipeline_table import get_record_by_pdb_id


def build_pipeline_records_from_input_csv(
    pdb_record_list: list[dict[str, str]],
    protein_data_dir: Path,
) -> list[dict[str, str]]:
    pipeline_record_list: list[dict[str, str]] = []

    for pdb_record in pdb_record_list:
        pdb_id = str(pdb_record.get(PDB_ID_COLUMN_NAME, "")).strip().upper()
        residue_range = str(pdb_record.get(RANGE_COLUMN_NAME, "")).strip()

        if not pdb_id:
            continue

        pdb_directory = protein_data_dir / pdb_id
        fasta_directory = pdb_directory / "fasta"
        alignment_directory = fasta_directory / "alignments"
        components_directory = pdb_directory / "components"
        prepared_directory = pdb_directory / "prepared"

        pipeline_record = create_protein_record(
            pdb_id=pdb_id,
            residue_range=residue_range,
        )
        pipeline_record[PDB_DIRECTORY_COLUMN_NAME] = str(pdb_directory)
        pipeline_record[FASTA_DIRECTORY_COLUMN_NAME] = str(fasta_directory)
        pipeline_record[ALIGNMENT_DIRECTORY_COLUMN_NAME] = str(alignment_directory)
        pipeline_record[COMPONENTS_DIRECTORY_COLUMN_NAME] = str(components_directory)
        pipeline_record[PREPARED_DIRECTORY_COLUMN_NAME] = str(prepared_directory)
        pipeline_record[PDB_SYNC_DONE_COLUMN_NAME] = STATUS_SUCCESS

        pipeline_record_list.append(pipeline_record)

    return pipeline_record_list


def merge_existing_and_new_pipeline_records(
    existing_pipeline_record_list: list[dict[str, str]],
    new_pipeline_record_list: list[dict[str, str]],
) -> list[dict[str, str]]:
    merged_record_list: list[dict[str, str]] = []

    for new_record in new_pipeline_record_list:
        pdb_id = new_record[PDB_ID_COLUMN_NAME]
        existing_record = get_record_by_pdb_id(existing_pipeline_record_list, pdb_id)

        if existing_record is None:
            merged_record_list.append(new_record)
            continue

        merged_record = dict(existing_record)
        merged_record[PDB_ID_COLUMN_NAME] = new_record[PDB_ID_COLUMN_NAME]
        merged_record[RANGE_COLUMN_NAME] = new_record[RANGE_COLUMN_NAME]
        merged_record[PDB_DIRECTORY_COLUMN_NAME] = new_record[PDB_DIRECTORY_COLUMN_NAME]
        merged_record[FASTA_DIRECTORY_COLUMN_NAME] = new_record[
            FASTA_DIRECTORY_COLUMN_NAME
        ]
        merged_record[ALIGNMENT_DIRECTORY_COLUMN_NAME] = new_record[
            ALIGNMENT_DIRECTORY_COLUMN_NAME
        ]
        merged_record[COMPONENTS_DIRECTORY_COLUMN_NAME] = new_record[
            COMPONENTS_DIRECTORY_COLUMN_NAME
        ]
        merged_record[PREPARED_DIRECTORY_COLUMN_NAME] = new_record[
            PREPARED_DIRECTORY_COLUMN_NAME
        ]
        merged_record[PDB_SYNC_DONE_COLUMN_NAME] = STATUS_SUCCESS

        merged_record_list.append(merged_record)

    return merged_record_list


def _clear_component_fields(pipeline_record: dict[str, str]) -> None:
    pipeline_record[HAS_METALS_COLUMN_NAME] = ""
    pipeline_record[HAS_LIGANDS_COLUMN_NAME] = ""
    pipeline_record[HAS_NONSTANDARD_RESIDUES_COLUMN_NAME] = ""


def _clear_gap_fields(pipeline_record: dict[str, str]) -> None:
    pipeline_record[N_GAPS_COLUMN_NAME] = ""
    pipeline_record[GAP_SIZES_COLUMN_NAME] = ""
    pipeline_record[HAS_GAPS_COLUMN_NAME] = ""
    pipeline_record[VARIANT_POLICY_COLUMN_NAME] = ""
    pipeline_record[RETAIN_GAPS_VARIANT_COLUMN_NAME] = ""
    pipeline_record[RETAIN_MODELLER_VARIANT_COLUMN_NAME] = ""
    pipeline_record[RETAIN_ALPHAFOLD_VARIANT_COLUMN_NAME] = ""


def _clear_filler_fields(pipeline_record: dict[str, str]) -> None:
    pipeline_record[FILLER_DIRECTORY_COLUMN_NAME] = ""
    pipeline_record[FILLER_MODEL_PATH_COLUMN_NAME] = ""
    pipeline_record[FILLER_MODEL_SOURCE_COLUMN_NAME] = ""
    pipeline_record[FILLER_STATUS_COLUMN_NAME] = ""


def _clear_protonation_fields(pipeline_record: dict[str, str]) -> None:
    pipeline_record[PROTONATION_STATUS_COLUMN_NAME] = ""
    pipeline_record[PROTONATION_INPUT_SOURCE_COLUMN_NAME] = ""
    pipeline_record[PROTONATION_INPUT_PATH_COLUMN_NAME] = ""
    pipeline_record[PROTONATION_OUTPUT_PATH_COLUMN_NAME] = ""


def _clear_internal_capping_fields(pipeline_record: dict[str, str]) -> None:
    pipeline_record[INTERNAL_CAPPING_STATUS_COLUMN_NAME] = ""
    pipeline_record[INTERNAL_CAPPING_INPUT_PATH_COLUMN_NAME] = ""
    pipeline_record[INTERNAL_CAPPING_OUTPUT_PATH_COLUMN_NAME] = ""


def _clear_prepared_structure_fields(pipeline_record: dict[str, str]) -> None:
    pipeline_record[PREPARED_STRUCTURE_STATUS_COLUMN_NAME] = ""
    pipeline_record[PREPARED_STRUCTURE_VARIANT_COLUMN_NAME] = ""
    pipeline_record[PREPARED_STRUCTURE_PROTEIN_INPUT_PATH_COLUMN_NAME] = ""
    pipeline_record[PREPARED_STRUCTURE_OUTPUT_PATH_COLUMN_NAME] = ""
    pipeline_record[PREPARED_GAPS_OUTPUT_PATH_COLUMN_NAME] = ""
    pipeline_record[PREPARED_MODELLER_OUTPUT_PATH_COLUMN_NAME] = ""
    pipeline_record[PREPARED_ALPHAFOLD_OUTPUT_PATH_COLUMN_NAME] = ""
    pipeline_record[AVAILABLE_MODELS_COLUMN_NAME] = ""


def _clear_downstream_after_gap_stage(pipeline_record: dict[str, str]) -> None:
    _clear_filler_fields(pipeline_record)
    _clear_protonation_fields(pipeline_record)
    _clear_internal_capping_fields(pipeline_record)
    _clear_prepared_structure_fields(pipeline_record)


def _set_component_flags(
    pipeline_record: dict[str, str],
    component_summary: dict[str, object],
) -> None:
    pipeline_record[HAS_METALS_COLUMN_NAME] = (
        "yes" if bool(component_summary.get("has_metals", False)) else "no"
    )
    pipeline_record[HAS_LIGANDS_COLUMN_NAME] = (
        "yes" if bool(component_summary.get("has_ligands", False)) else "no"
    )
    pipeline_record[HAS_NONSTANDARD_RESIDUES_COLUMN_NAME] = (
        "yes"
        if bool(component_summary.get("has_nonstandard_residues", False))
        else "no"
    )


def _clear_prepared_variant_output_paths(pipeline_record: dict[str, str]) -> None:
    """Clear prepared-variant path fields before writing accepted end models."""

    pipeline_record[PREPARED_STRUCTURE_VARIANT_COLUMN_NAME] = ""
    pipeline_record[PREPARED_STRUCTURE_PROTEIN_INPUT_PATH_COLUMN_NAME] = ""
    pipeline_record[PREPARED_STRUCTURE_OUTPUT_PATH_COLUMN_NAME] = ""
    pipeline_record[PREPARED_GAPS_OUTPUT_PATH_COLUMN_NAME] = ""
    pipeline_record[PREPARED_MODELLER_OUTPUT_PATH_COLUMN_NAME] = ""
    pipeline_record[PREPARED_ALPHAFOLD_OUTPUT_PATH_COLUMN_NAME] = ""
    pipeline_record[AVAILABLE_MODELS_COLUMN_NAME] = ""


def _recompute_prepared_summary_from_state(
    *,
    summary: dict[str, int],
    pipeline_record_list: list[dict[str, str]],
) -> None:
    """Recompute final prepared-output counters after variant audit decisions."""

    summary["single_best_available_written"] = 0
    summary["gaps_variant_written"] = 0
    summary["modeller_complete_written"] = 0
    summary["alphafold_complete_written"] = 0
    summary["proteins_with_prepared_output"] = 0
    summary["proteins_failed"] = 0

    for pipeline_record in pipeline_record_list:
        prepared_status = str(
            pipeline_record.get(PREPARED_STRUCTURE_STATUS_COLUMN_NAME, "")
        ).strip()
        variant_text = str(
            pipeline_record.get(PREPARED_STRUCTURE_VARIANT_COLUMN_NAME, "")
        ).strip()
        variant_set = {
            token.strip()
            for token in variant_text.split(",")
            if token.strip()
        }

        if prepared_status == STATUS_SUCCESS and variant_set:
            summary["proteins_with_prepared_output"] += 1
        else:
            summary["proteins_failed"] += 1
            continue

        if "single" in variant_set:
            summary["single_best_available_written"] += 1
        if "gaps" in variant_set:
            summary["gaps_variant_written"] += 1
        if str(pipeline_record.get(PREPARED_MODELLER_OUTPUT_PATH_COLUMN_NAME, "")).strip():
            summary["modeller_complete_written"] += 1
        if str(pipeline_record.get(PREPARED_ALPHAFOLD_OUTPUT_PATH_COLUMN_NAME, "")).strip():
            summary["alphafold_complete_written"] += 1
