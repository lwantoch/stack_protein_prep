"""
/home/grheco/repositorios/stack_protein_prep/scripts/fruton.py

Main entry point for the FRUTON protein preparation pipeline.

Responsibilities
----------------
- define working directories
- make the src directory importable
- synchronize input CSV and protein directories
- build and merge pipeline records
- run the implemented pipeline modules in order
- save pipeline state to JSON
- export pipeline overview to XLSX

Implemented steps
-----------------
1. pdb_sync
2. fasta_files
3. sequence_alignment
4. insertion_codes
5. component_split
6. gap_detection
7. filler
8. protonation
9. amber_renaming
10. amber_termini
11. internal_capping
12. prepared_structure

Important
---------
- pipeline state is persisted in pipeline.json
- pipeline overview is exported to pipeline.xlsx
- this script keeps the state flat and practical
- downstream fields are cleared whenever an upstream dependency fails
- gap-aware prepared output is handled as:
    - prepared/<PDBID>.pdb                   if no gaps
    - prepared/gaps/<PDBID>.pdb             if gaps remain
    - prepared/complete/<PDBID>.pdb         if a filled model was used

Current limitation
------------------
This script builds the currently available prepared variant cleanly.
To build BOTH gaps and complete variants side-by-side for the same protein,
the protonation/renaming/termini/capping modules would need fully separate
variant-specific output names all the way through the pipeline.
"""

from __future__ import annotations

import re
import sys
import traceback
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_ROOT_DIR = SCRIPT_DIR.parent
SRC_DIR = PROJECT_ROOT_DIR / "src"

if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

from stack_protein_preparation.amber_renaming import amber_rename_protein_structure
from stack_protein_preparation.cap import cap_internal_gaps_for_pdb_directory
from stack_protein_preparation.fasta_files import create_fasta_files_for_pdb_directory
from stack_protein_preparation.filler import run_filler_for_chain
from stack_protein_preparation.gaps import summarize_gaps
from stack_protein_preparation.insertion_codes import (
    STATUS_FAILED as INSERTION_STATUS_FAILED,
)
from stack_protein_preparation.insertion_codes import (
    STATUS_NONE as INSERTION_STATUS_NONE,
)
from stack_protein_preparation.insertion_codes import (
    STATUS_SUCCESS as INSERTION_STATUS_SUCCESS,
)
from stack_protein_preparation.insertion_codes import (
    find_input_pdb_for_protein,
    process_pdb_for_delinsertion,
)
from stack_protein_preparation.pdb_components import split_pdb_components
from stack_protein_preparation.pdb_sync import (
    read_pdb_records_from_csv,
    sync_pdb_csv_and_directories,
)
from stack_protein_preparation.pipeline_state import (
    ALIGNMENT_DIRECTORY_COLUMN_NAME,
    AMBER_INPUT_PATH_COLUMN_NAME,
    AMBER_OUTPUT_PATH_COLUMN_NAME,
    AMBER_RENAMING_STATUS_COLUMN_NAME,
    AMBER_TERMINI_INPUT_PATH_COLUMN_NAME,
    AMBER_TERMINI_OUTPUT_PATH_COLUMN_NAME,
    AMBER_TERMINI_STATUS_COLUMN_NAME,
    COMPONENTS_DIRECTORY_COLUMN_NAME,
    FASTA_DIRECTORY_COLUMN_NAME,
    FASTA_FILES_DONE_COLUMN_NAME,
    FILLER_DIRECTORY_COLUMN_NAME,
    FILLER_MODEL_PATH_COLUMN_NAME,
    FILLER_STATUS_COLUMN_NAME,
    GAP_SIZES_COLUMN_NAME,
    HAS_GAPS_COLUMN_NAME,
    HAS_LIGANDS_COLUMN_NAME,
    HAS_METALS_COLUMN_NAME,
    HAS_NONSTANDARD_RESIDUES_COLUMN_NAME,
    INSERTION_CODES_DONE_COLUMN_NAME,
    INTERNAL_CAPPING_INPUT_PATH_COLUMN_NAME,
    INTERNAL_CAPPING_OUTPUT_PATH_COLUMN_NAME,
    INTERNAL_CAPPING_STATUS_COLUMN_NAME,
    N_GAPS_COLUMN_NAME,
    PDB_DIRECTORY_COLUMN_NAME,
    PDB_ID_COLUMN_NAME,
    PDB_SYNC_DONE_COLUMN_NAME,
    PREPARED_DIRECTORY_COLUMN_NAME,
    PREPARED_STRUCTURE_OUTPUT_PATH_COLUMN_NAME,
    PREPARED_STRUCTURE_PROTEIN_INPUT_PATH_COLUMN_NAME,
    PREPARED_STRUCTURE_STATUS_COLUMN_NAME,
    PREPARED_STRUCTURE_VARIANT_COLUMN_NAME,
    PROTONATION_INPUT_PATH_COLUMN_NAME,
    PROTONATION_INPUT_SOURCE_COLUMN_NAME,
    PROTONATION_OUTPUT_PATH_COLUMN_NAME,
    PROTONATION_STATUS_COLUMN_NAME,
    RANGE_COLUMN_NAME,
    SEQUENCE_ALIGNMENT_DONE_COLUMN_NAME,
    STATUS_FAILED,
    STATUS_REQUIRED,
    STATUS_SKIPPED,
    STATUS_SUCCESS,
    STATUS_WARNING,
    UNIPROT_ID_COLUMN_NAME,
    create_protein_record,
)
from stack_protein_preparation.pipeline_table import (
    get_record_by_pdb_id,
    load_pipeline_table,
    save_pipeline_table,
)
from stack_protein_preparation.pipeline_xlsx import write_pipeline_to_xlsx
from stack_protein_preparation.prepared_structure import (
    build_prepared_structure_for_pdb_directory,
)
from stack_protein_preparation.protonation import protonate_protein_structure
from stack_protein_preparation.sequence_alignment import (
    run_alignments_for_pdb_directory,
)
from stack_protein_preparation.terminus import (
    convert_protein_termini_for_pdb_directory,
)


def build_pipeline_records_from_input_csv(
    pdb_record_list: list[dict[str, str]],
    protein_data_dir: Path,
) -> list[dict[str, str]]:
    """
    Build initial pipeline records from the input PDB CSV records.
    """
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
    """
    Merge previously saved pipeline state with current input records.
    """
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


def _find_uniprot_id_for_protein(pdb_dir: Path) -> str:
    """
    Find UniProt ID from FASTA files in the protein directory.
    """
    fasta_dir = pdb_dir / "fasta"

    if not fasta_dir.exists():
        return ""

    for fasta_path in sorted(fasta_dir.glob("UniProt_*.fasta")):
        match = re.match(r"^UniProt_([A-Z0-9]+)\.fasta$", fasta_path.name)
        if match:
            return match.group(1)

    return ""


def _find_template_pdb_for_filler(
    pdb_id: str,
    pdb_dir: Path,
    components_dir: Path,
) -> Path | None:
    """
    Find the template PDB to be used by the filler step.
    """
    candidate_path_list = [
        components_dir / f"{pdb_id}_protein.pdb",
        pdb_dir / f"{pdb_id}_delins.pdb",
        pdb_dir / f"{pdb_id}.pdb",
    ]

    for candidate_path in candidate_path_list:
        if candidate_path.exists():
            return candidate_path

    return None


def _clear_component_fields(pipeline_record: dict[str, str]) -> None:
    pipeline_record[HAS_METALS_COLUMN_NAME] = ""
    pipeline_record[HAS_LIGANDS_COLUMN_NAME] = ""
    pipeline_record[HAS_NONSTANDARD_RESIDUES_COLUMN_NAME] = ""


def _clear_gap_fields(pipeline_record: dict[str, str]) -> None:
    pipeline_record[N_GAPS_COLUMN_NAME] = ""
    pipeline_record[GAP_SIZES_COLUMN_NAME] = ""
    pipeline_record[HAS_GAPS_COLUMN_NAME] = ""


def _clear_filler_fields(pipeline_record: dict[str, str]) -> None:
    pipeline_record[FILLER_DIRECTORY_COLUMN_NAME] = ""
    pipeline_record[FILLER_MODEL_PATH_COLUMN_NAME] = ""
    pipeline_record[FILLER_STATUS_COLUMN_NAME] = ""


def _clear_protonation_fields(pipeline_record: dict[str, str]) -> None:
    pipeline_record[PROTONATION_STATUS_COLUMN_NAME] = ""
    pipeline_record[PROTONATION_INPUT_SOURCE_COLUMN_NAME] = ""
    pipeline_record[PROTONATION_INPUT_PATH_COLUMN_NAME] = ""
    pipeline_record[PROTONATION_OUTPUT_PATH_COLUMN_NAME] = ""


def _clear_amber_renaming_fields(pipeline_record: dict[str, str]) -> None:
    pipeline_record[AMBER_RENAMING_STATUS_COLUMN_NAME] = ""
    pipeline_record[AMBER_INPUT_PATH_COLUMN_NAME] = ""
    pipeline_record[AMBER_OUTPUT_PATH_COLUMN_NAME] = ""


def _clear_amber_termini_fields(pipeline_record: dict[str, str]) -> None:
    pipeline_record[AMBER_TERMINI_STATUS_COLUMN_NAME] = ""
    pipeline_record[AMBER_TERMINI_INPUT_PATH_COLUMN_NAME] = ""
    pipeline_record[AMBER_TERMINI_OUTPUT_PATH_COLUMN_NAME] = ""


def _clear_internal_capping_fields(pipeline_record: dict[str, str]) -> None:
    pipeline_record[INTERNAL_CAPPING_STATUS_COLUMN_NAME] = ""
    pipeline_record[INTERNAL_CAPPING_INPUT_PATH_COLUMN_NAME] = ""
    pipeline_record[INTERNAL_CAPPING_OUTPUT_PATH_COLUMN_NAME] = ""


def _clear_prepared_structure_fields(pipeline_record: dict[str, str]) -> None:
    pipeline_record[PREPARED_STRUCTURE_STATUS_COLUMN_NAME] = ""
    pipeline_record[PREPARED_STRUCTURE_VARIANT_COLUMN_NAME] = ""
    pipeline_record[PREPARED_STRUCTURE_PROTEIN_INPUT_PATH_COLUMN_NAME] = ""
    pipeline_record[PREPARED_STRUCTURE_OUTPUT_PATH_COLUMN_NAME] = ""


def _clear_downstream_from_fasta_failure(pipeline_record: dict[str, str]) -> None:
    pipeline_record[SEQUENCE_ALIGNMENT_DONE_COLUMN_NAME] = STATUS_REQUIRED
    pipeline_record[INSERTION_CODES_DONE_COLUMN_NAME] = STATUS_REQUIRED
    _clear_component_fields(pipeline_record)
    _clear_gap_fields(pipeline_record)
    _clear_filler_fields(pipeline_record)
    _clear_protonation_fields(pipeline_record)
    _clear_amber_renaming_fields(pipeline_record)
    _clear_amber_termini_fields(pipeline_record)
    _clear_internal_capping_fields(pipeline_record)
    _clear_prepared_structure_fields(pipeline_record)


def _clear_downstream_from_alignment_failure(pipeline_record: dict[str, str]) -> None:
    pipeline_record[INSERTION_CODES_DONE_COLUMN_NAME] = STATUS_REQUIRED
    _clear_component_fields(pipeline_record)
    _clear_gap_fields(pipeline_record)
    _clear_filler_fields(pipeline_record)
    _clear_protonation_fields(pipeline_record)
    _clear_amber_renaming_fields(pipeline_record)
    _clear_amber_termini_fields(pipeline_record)
    _clear_internal_capping_fields(pipeline_record)
    _clear_prepared_structure_fields(pipeline_record)


def _clear_downstream_from_insertion_failure(pipeline_record: dict[str, str]) -> None:
    _clear_component_fields(pipeline_record)
    _clear_gap_fields(pipeline_record)
    _clear_filler_fields(pipeline_record)
    _clear_protonation_fields(pipeline_record)
    _clear_amber_renaming_fields(pipeline_record)
    _clear_amber_termini_fields(pipeline_record)
    _clear_internal_capping_fields(pipeline_record)
    _clear_prepared_structure_fields(pipeline_record)


def _clear_downstream_from_component_failure(pipeline_record: dict[str, str]) -> None:
    _clear_gap_fields(pipeline_record)
    _clear_filler_fields(pipeline_record)
    _clear_protonation_fields(pipeline_record)
    _clear_amber_renaming_fields(pipeline_record)
    _clear_amber_termini_fields(pipeline_record)
    _clear_internal_capping_fields(pipeline_record)
    _clear_prepared_structure_fields(pipeline_record)


def _clear_downstream_from_gap_failure(pipeline_record: dict[str, str]) -> None:
    _clear_filler_fields(pipeline_record)
    _clear_protonation_fields(pipeline_record)
    _clear_amber_renaming_fields(pipeline_record)
    _clear_amber_termini_fields(pipeline_record)
    _clear_internal_capping_fields(pipeline_record)
    _clear_prepared_structure_fields(pipeline_record)


def _clear_downstream_from_filler_failure(pipeline_record: dict[str, str]) -> None:
    _clear_protonation_fields(pipeline_record)
    _clear_amber_renaming_fields(pipeline_record)
    _clear_amber_termini_fields(pipeline_record)
    _clear_internal_capping_fields(pipeline_record)
    _clear_prepared_structure_fields(pipeline_record)


def _clear_downstream_from_protonation_failure(pipeline_record: dict[str, str]) -> None:
    _clear_amber_renaming_fields(pipeline_record)
    _clear_amber_termini_fields(pipeline_record)
    _clear_internal_capping_fields(pipeline_record)
    _clear_prepared_structure_fields(pipeline_record)


def _clear_downstream_from_amber_failure(pipeline_record: dict[str, str]) -> None:
    _clear_amber_termini_fields(pipeline_record)
    _clear_internal_capping_fields(pipeline_record)
    _clear_prepared_structure_fields(pipeline_record)


def _clear_downstream_from_termini_failure(pipeline_record: dict[str, str]) -> None:
    _clear_internal_capping_fields(pipeline_record)
    _clear_prepared_structure_fields(pipeline_record)


def _clear_downstream_from_internal_capping_failure(
    pipeline_record: dict[str, str],
) -> None:
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


def _run_protonation_for_protein(
    pdb_id: str,
    protein_dir: Path,
    modeller_model_path: Path | None,
    alphafold_model_path: Path | None,
) -> dict[str, str]:
    """
    Run protonation for one protein and return a flat result dictionary.
    """
    try:
        result = protonate_protein_structure(
            pdb_id=pdb_id,
            protein_dir=protein_dir,
            modeller_model_path=modeller_model_path,
            alphafold_model_path=alphafold_model_path,
            ph=7.4,
        )

        return {
            "status": (
                STATUS_SUCCESS
                if result.get("protonation_success", False)
                else STATUS_FAILED
            ),
            "input_source": str(result.get("protonation_input_source", "")).strip(),
            "input_path": str(result.get("protonation_input_path", "")).strip(),
            "output_path": str(result.get("protonation_output_path", "")).strip(),
        }

    except Exception as error:
        return {
            "status": STATUS_FAILED,
            "input_source": "",
            "input_path": "",
            "output_path": "",
            "error": str(error),
        }


def _run_amber_renaming_for_protein(
    pdb_id: str,
    protein_dir: Path,
) -> dict[str, str]:
    """
    Run AMBER residue renaming for one protein and return a flat result dictionary.
    """
    try:
        result = amber_rename_protein_structure(
            pdb_id=pdb_id,
            protein_dir=protein_dir,
            strict_his=False,
            disulf_min=1.8,
            disulf_max=2.2,
        )

        return {
            "status": (
                STATUS_SUCCESS
                if result.get("amber_renaming_success", False)
                else STATUS_FAILED
            ),
            "input_path": str(result.get("amber_input_path", "")).strip(),
            "output_path": str(result.get("amber_output_path", "")).strip(),
        }

    except Exception as error:
        return {
            "status": STATUS_FAILED,
            "input_path": "",
            "output_path": "",
            "error": str(error),
        }


def _determine_prepared_variant(pipeline_record: dict[str, str]) -> str | None:
    """
    Determine which prepared variant should be written.

    Rules
    -----
    - no gaps -> None  (plain prepared/<PDBID>.pdb)
    - gaps + filler model present -> "complete"
    - gaps without filler model -> "gaps"
    """
    has_gaps_text = str(pipeline_record.get(HAS_GAPS_COLUMN_NAME, "")).strip().lower()

    if has_gaps_text != "yes":
        return None

    filler_model_path = str(
        pipeline_record.get(FILLER_MODEL_PATH_COLUMN_NAME, "")
    ).strip()
    if filler_model_path:
        return "complete"

    return "gaps"


def run_pipeline() -> None:
    """
    Run the FRUTON protein preparation pipeline.
    """
    protein_data_dir = PROJECT_ROOT_DIR / "data" / "proteins"
    pdb_ids_csv_path = protein_data_dir / "pdb_ids.csv"
    pipeline_json_path = protein_data_dir / "pipeline.json"
    pipeline_xlsx_path = protein_data_dir / "pipeline.xlsx"

    print("[FRUTON] Starting")

    print("[FRUTON] Step 1: pdb_sync")
    sync_pdb_csv_and_directories(protein_data_dir)

    print("[FRUTON] Step 2: read_input_csv")
    pdb_record_list = read_pdb_records_from_csv(pdb_ids_csv_path)
    print(f"[FRUTON] Loaded {len(pdb_record_list)} input records")

    print("[FRUTON] Step 3: load_existing_pipeline_state")
    existing_pipeline_record_list = load_pipeline_table(pipeline_json_path)
    print(
        f"[FRUTON] Loaded {len(existing_pipeline_record_list)} existing pipeline records"
    )

    print("[FRUTON] Step 4: build_and_merge_pipeline_state")
    new_pipeline_record_list = build_pipeline_records_from_input_csv(
        pdb_record_list=pdb_record_list,
        protein_data_dir=protein_data_dir,
    )
    pipeline_record_list = merge_existing_and_new_pipeline_records(
        existing_pipeline_record_list=existing_pipeline_record_list,
        new_pipeline_record_list=new_pipeline_record_list,
    )
    print(f"[FRUTON] Active pipeline records: {len(pipeline_record_list)}")

    print("[FRUTON] Step 5: fasta_files")
    for pipeline_record in pipeline_record_list:
        pdb_id = pipeline_record[PDB_ID_COLUMN_NAME]
        pdb_dir = Path(pipeline_record[PDB_DIRECTORY_COLUMN_NAME])

        print(f"[FRUTON] fasta_files -> {pdb_id}")

        try:
            create_fasta_files_for_pdb_directory(pdb_dir)
            pipeline_record[FASTA_FILES_DONE_COLUMN_NAME] = STATUS_SUCCESS
            pipeline_record[UNIPROT_ID_COLUMN_NAME] = _find_uniprot_id_for_protein(
                pdb_dir
            )
        except Exception as error:
            print(f"[FRUTON] fasta_files failed for {pdb_id}: {error!r}")
            pipeline_record[FASTA_FILES_DONE_COLUMN_NAME] = STATUS_FAILED
            pipeline_record[UNIPROT_ID_COLUMN_NAME] = ""
            _clear_downstream_from_fasta_failure(pipeline_record)

    print("[FRUTON] Step 6: sequence_alignment")
    for pipeline_record in pipeline_record_list:
        pdb_id = pipeline_record[PDB_ID_COLUMN_NAME]
        pdb_dir = Path(pipeline_record[PDB_DIRECTORY_COLUMN_NAME])

        if pipeline_record.get(FASTA_FILES_DONE_COLUMN_NAME, "") != STATUS_SUCCESS:
            print(f"[FRUTON] sequence_alignment skipped for {pdb_id}")
            pipeline_record[SEQUENCE_ALIGNMENT_DONE_COLUMN_NAME] = STATUS_SKIPPED
            _clear_downstream_from_alignment_failure(pipeline_record)
            continue

        print(f"[FRUTON] sequence_alignment -> {pdb_id}")

        try:
            run_alignments_for_pdb_directory(pdb_dir)
            pipeline_record[SEQUENCE_ALIGNMENT_DONE_COLUMN_NAME] = STATUS_SUCCESS

            if not str(pipeline_record.get(UNIPROT_ID_COLUMN_NAME, "")).strip():
                pipeline_record[UNIPROT_ID_COLUMN_NAME] = _find_uniprot_id_for_protein(
                    pdb_dir
                )
        except Exception as error:
            print(f"[FRUTON] sequence_alignment failed for {pdb_id}: {error!r}")
            pipeline_record[SEQUENCE_ALIGNMENT_DONE_COLUMN_NAME] = STATUS_FAILED
            _clear_downstream_from_alignment_failure(pipeline_record)

    print("[FRUTON] Step 7: insertion_codes")
    for pipeline_record in pipeline_record_list:
        pdb_id = pipeline_record[PDB_ID_COLUMN_NAME]
        pdb_dir = Path(pipeline_record[PDB_DIRECTORY_COLUMN_NAME])

        if (
            pipeline_record.get(SEQUENCE_ALIGNMENT_DONE_COLUMN_NAME, "")
            != STATUS_SUCCESS
        ):
            print(f"[FRUTON] insertion_codes skipped for {pdb_id}")
            pipeline_record[INSERTION_CODES_DONE_COLUMN_NAME] = STATUS_SKIPPED
            _clear_downstream_from_insertion_failure(pipeline_record)
            continue

        print(f"[FRUTON] insertion_codes -> {pdb_id}")

        input_pdb_path = find_input_pdb_for_protein(pdb_dir)
        if input_pdb_path is None:
            print(f"[FRUTON] insertion_codes failed for {pdb_id}: input PDB not found")
            pipeline_record[INSERTION_CODES_DONE_COLUMN_NAME] = STATUS_FAILED
            _clear_downstream_from_insertion_failure(pipeline_record)
            continue

        output_pdb_path = pdb_dir / f"{pdb_id}_delins.pdb"

        try:
            insertion_result = process_pdb_for_delinsertion(
                input_pdb_path=input_pdb_path,
                output_pdb_path=output_pdb_path,
            )
        except Exception as error:
            print(f"[FRUTON] insertion_codes failed for {pdb_id}: {error!r}")
            pipeline_record[INSERTION_CODES_DONE_COLUMN_NAME] = STATUS_FAILED
            _clear_downstream_from_insertion_failure(pipeline_record)
            continue

        insertion_status = str(insertion_result.get("status", "")).strip().lower()

        if insertion_status in {INSERTION_STATUS_NONE, INSERTION_STATUS_SUCCESS}:
            pipeline_record[INSERTION_CODES_DONE_COLUMN_NAME] = STATUS_SUCCESS
        elif insertion_status == INSERTION_STATUS_FAILED:
            pipeline_record[INSERTION_CODES_DONE_COLUMN_NAME] = STATUS_FAILED
            _clear_downstream_from_insertion_failure(pipeline_record)
        else:
            pipeline_record[INSERTION_CODES_DONE_COLUMN_NAME] = STATUS_WARNING
            _clear_downstream_from_insertion_failure(pipeline_record)

    print("[FRUTON] Step 8: component_split")
    for pipeline_record in pipeline_record_list:
        pdb_id = pipeline_record[PDB_ID_COLUMN_NAME]
        pdb_dir = Path(pipeline_record[PDB_DIRECTORY_COLUMN_NAME])
        components_dir = Path(pipeline_record[COMPONENTS_DIRECTORY_COLUMN_NAME])

        if pipeline_record.get(INSERTION_CODES_DONE_COLUMN_NAME, "") != STATUS_SUCCESS:
            print(f"[FRUTON] component_split skipped for {pdb_id}")
            _clear_component_fields(pipeline_record)
            _clear_downstream_from_component_failure(pipeline_record)
            continue

        component_input_pdb_path = pdb_dir / f"{pdb_id}_delins.pdb"

        if not component_input_pdb_path.exists():
            print(f"[FRUTON] component_split skipped for {pdb_id}: missing input")
            _clear_component_fields(pipeline_record)
            _clear_downstream_from_component_failure(pipeline_record)
            continue

        print(f"[FRUTON] component_split -> {pdb_id}")

        try:
            component_summary = split_pdb_components(
                pdb_path=component_input_pdb_path,
                output_dir=components_dir,
                protein_stem=pdb_id,
            )
        except Exception as error:
            print(f"[FRUTON] component_split failed for {pdb_id}: {error!r}")
            _clear_component_fields(pipeline_record)
            _clear_downstream_from_component_failure(pipeline_record)
            continue

        pipeline_record[COMPONENTS_DIRECTORY_COLUMN_NAME] = str(
            component_summary["output_dir"]
        )
        _set_component_flags(pipeline_record, component_summary)

    print("[FRUTON] Step 9: gap_detection")
    for pipeline_record in pipeline_record_list:
        pdb_id = pipeline_record[PDB_ID_COLUMN_NAME]
        components_dir = Path(pipeline_record[COMPONENTS_DIRECTORY_COLUMN_NAME])

        if pipeline_record.get(INSERTION_CODES_DONE_COLUMN_NAME, "") != STATUS_SUCCESS:
            print(f"[FRUTON] gap_detection skipped for {pdb_id}")
            _clear_gap_fields(pipeline_record)
            _clear_downstream_from_gap_failure(pipeline_record)
            continue

        gap_input_pdb_path = components_dir / f"{pdb_id}_protein.pdb"

        if not gap_input_pdb_path.exists():
            print(f"[FRUTON] gap_detection skipped for {pdb_id}: missing protein file")
            _clear_gap_fields(pipeline_record)
            _clear_downstream_from_gap_failure(pipeline_record)
            continue

        print(f"[FRUTON] gap_detection -> {pdb_id}")

        try:
            gap_summary = summarize_gaps(gap_input_pdb_path)
        except Exception as error:
            print(f"[FRUTON] gap_detection failed for {pdb_id}: {error!r}")
            _clear_gap_fields(pipeline_record)
            _clear_downstream_from_gap_failure(pipeline_record)
            continue

        n_gaps = int(gap_summary.get("n_gaps", 0))
        gap_sizes = gap_summary.get("gap_sizes", [])

        pipeline_record[N_GAPS_COLUMN_NAME] = str(n_gaps)
        pipeline_record[GAP_SIZES_COLUMN_NAME] = (
            "none" if n_gaps == 0 else "|".join(str(gap_size) for gap_size in gap_sizes)
        )
        pipeline_record[HAS_GAPS_COLUMN_NAME] = "yes" if n_gaps > 0 else "no"

    print("[FRUTON] Step 10: filler")
    for pipeline_record in pipeline_record_list:
        pdb_id = pipeline_record[PDB_ID_COLUMN_NAME]
        pdb_dir = Path(pipeline_record[PDB_DIRECTORY_COLUMN_NAME])
        alignment_dir = Path(pipeline_record[ALIGNMENT_DIRECTORY_COLUMN_NAME])
        components_dir = Path(pipeline_record[COMPONENTS_DIRECTORY_COLUMN_NAME])

        if pipeline_record.get(HAS_GAPS_COLUMN_NAME, "") != "yes":
            print(f"[FRUTON] filler skipped for {pdb_id}: no gaps")
            _clear_filler_fields(pipeline_record)
            pipeline_record[FILLER_STATUS_COLUMN_NAME] = STATUS_SKIPPED
            continue

        if pipeline_record.get(INSERTION_CODES_DONE_COLUMN_NAME, "") != STATUS_SUCCESS:
            print(f"[FRUTON] filler skipped for {pdb_id}: upstream failure")
            _clear_filler_fields(pipeline_record)
            pipeline_record[FILLER_STATUS_COLUMN_NAME] = STATUS_SKIPPED
            _clear_downstream_from_filler_failure(pipeline_record)
            continue

        template_pdb_path = _find_template_pdb_for_filler(
            pdb_id=pdb_id,
            pdb_dir=pdb_dir,
            components_dir=components_dir,
        )

        if template_pdb_path is None or not alignment_dir.exists():
            print(f"[FRUTON] filler skipped for {pdb_id}: missing inputs")
            _clear_filler_fields(pipeline_record)
            pipeline_record[FILLER_STATUS_COLUMN_NAME] = STATUS_WARNING
            _clear_downstream_from_filler_failure(pipeline_record)
            continue

        print(f"[FRUTON] filler -> {pdb_id}")

        chain_id = "A"
        template_id = template_pdb_path.stem
        target_id = f"{pdb_id}_{chain_id}"
        filler_output_dir = alignment_dir / "filler" / chain_id
        uniprot_id = str(pipeline_record.get(UNIPROT_ID_COLUMN_NAME, "")).strip()

        try:
            filler_result = run_filler_for_chain(
                alignment_directory=alignment_dir,
                template_pdb_path=template_pdb_path,
                output_dir=filler_output_dir,
                template_id=template_id,
                target_id=target_id,
                chain_id=chain_id,
                final_model_name=f"{pdb_id}_protein_mod.pdb",
                starting_model=1,
                ending_model=20,
                uniprot_id=uniprot_id,
                residue_range=pipeline_record[RANGE_COLUMN_NAME],
            )
        except Exception as error:
            print(f"[FRUTON] filler failed for {pdb_id}: {error!r}")
            traceback.print_exc()
            _clear_filler_fields(pipeline_record)
            pipeline_record[FILLER_STATUS_COLUMN_NAME] = STATUS_FAILED
            _clear_downstream_from_filler_failure(pipeline_record)
            continue

        pipeline_record[FILLER_DIRECTORY_COLUMN_NAME] = str(filler_result.output_dir)
        pipeline_record[FILLER_MODEL_PATH_COLUMN_NAME] = (
            str(filler_result.final_model_path)
            if filler_result.final_model_path is not None
            else ""
        )

        if filler_result.final_model_path is not None:
            pipeline_record[FILLER_STATUS_COLUMN_NAME] = STATUS_SUCCESS
        elif filler_result.skipped:
            pipeline_record[FILLER_STATUS_COLUMN_NAME] = STATUS_WARNING
        else:
            pipeline_record[FILLER_STATUS_COLUMN_NAME] = STATUS_WARNING

    print("[FRUTON] Step 11: protonation")
    for pipeline_record in pipeline_record_list:
        pdb_id = pipeline_record[PDB_ID_COLUMN_NAME]
        protein_dir = Path(pipeline_record[PDB_DIRECTORY_COLUMN_NAME])
        components_dir = Path(pipeline_record[COMPONENTS_DIRECTORY_COLUMN_NAME])

        if pipeline_record.get(INSERTION_CODES_DONE_COLUMN_NAME, "") != STATUS_SUCCESS:
            print(f"[FRUTON] protonation skipped for {pdb_id}")
            _clear_protonation_fields(pipeline_record)
            pipeline_record[PROTONATION_STATUS_COLUMN_NAME] = STATUS_SKIPPED
            _clear_downstream_from_protonation_failure(pipeline_record)
            continue

        default_protein_path = components_dir / f"{pdb_id}_protein.pdb"
        if not default_protein_path.exists():
            print(f"[FRUTON] protonation skipped for {pdb_id}: missing protein file")
            _clear_protonation_fields(pipeline_record)
            pipeline_record[PROTONATION_STATUS_COLUMN_NAME] = STATUS_FAILED
            _clear_downstream_from_protonation_failure(pipeline_record)
            continue

        filler_model_path_text = str(
            pipeline_record.get(FILLER_MODEL_PATH_COLUMN_NAME, "")
        ).strip()

        modeller_model_path: Path | None = None
        alphafold_model_path: Path | None = None

        if filler_model_path_text:
            filler_model_path = Path(filler_model_path_text)
            if filler_model_path.exists():
                lower_name = filler_model_path.name.lower()
                if "alphafold" in lower_name:
                    alphafold_model_path = filler_model_path
                else:
                    modeller_model_path = filler_model_path

        print(f"[FRUTON] protonation -> {pdb_id}")

        protonation_result = _run_protonation_for_protein(
            pdb_id=pdb_id,
            protein_dir=protein_dir,
            modeller_model_path=modeller_model_path,
            alphafold_model_path=alphafold_model_path,
        )

        pipeline_record[PROTONATION_STATUS_COLUMN_NAME] = protonation_result["status"]
        pipeline_record[PROTONATION_INPUT_SOURCE_COLUMN_NAME] = protonation_result[
            "input_source"
        ]
        pipeline_record[PROTONATION_INPUT_PATH_COLUMN_NAME] = protonation_result[
            "input_path"
        ]
        pipeline_record[PROTONATION_OUTPUT_PATH_COLUMN_NAME] = protonation_result[
            "output_path"
        ]

        if pipeline_record[PROTONATION_STATUS_COLUMN_NAME] != STATUS_SUCCESS:
            _clear_downstream_from_protonation_failure(pipeline_record)

    print("[FRUTON] Step 12: amber_renaming")
    for pipeline_record in pipeline_record_list:
        pdb_id = pipeline_record[PDB_ID_COLUMN_NAME]
        protein_dir = Path(pipeline_record[PDB_DIRECTORY_COLUMN_NAME])

        if pipeline_record.get(PROTONATION_STATUS_COLUMN_NAME, "") != STATUS_SUCCESS:
            print(f"[FRUTON] amber_renaming skipped for {pdb_id}")
            _clear_amber_renaming_fields(pipeline_record)
            pipeline_record[AMBER_RENAMING_STATUS_COLUMN_NAME] = STATUS_SKIPPED
            _clear_downstream_from_amber_failure(pipeline_record)
            continue

        print(f"[FRUTON] amber_renaming -> {pdb_id}")

        amber_result = _run_amber_renaming_for_protein(
            pdb_id=pdb_id,
            protein_dir=protein_dir,
        )

        pipeline_record[AMBER_RENAMING_STATUS_COLUMN_NAME] = amber_result["status"]
        pipeline_record[AMBER_INPUT_PATH_COLUMN_NAME] = amber_result["input_path"]
        pipeline_record[AMBER_OUTPUT_PATH_COLUMN_NAME] = amber_result["output_path"]

        if pipeline_record[AMBER_RENAMING_STATUS_COLUMN_NAME] != STATUS_SUCCESS:
            _clear_downstream_from_amber_failure(pipeline_record)

    print("[FRUTON] Step 13: amber_termini")
    for pipeline_record in pipeline_record_list:
        pdb_id = pipeline_record[PDB_ID_COLUMN_NAME]
        pdb_dir = Path(pipeline_record[PDB_DIRECTORY_COLUMN_NAME])

        if pipeline_record.get(AMBER_RENAMING_STATUS_COLUMN_NAME, "") != STATUS_SUCCESS:
            print(f"[FRUTON] amber_termini skipped for {pdb_id}")
            _clear_amber_termini_fields(pipeline_record)
            pipeline_record[AMBER_TERMINI_STATUS_COLUMN_NAME] = STATUS_SKIPPED
            _clear_downstream_from_termini_failure(pipeline_record)
            continue

        amber_input_path_text = str(
            pipeline_record.get(AMBER_OUTPUT_PATH_COLUMN_NAME, "")
        ).strip()

        if not amber_input_path_text or not Path(amber_input_path_text).exists():
            print(f"[FRUTON] amber_termini skipped for {pdb_id}: missing amber input")
            _clear_amber_termini_fields(pipeline_record)
            pipeline_record[AMBER_TERMINI_STATUS_COLUMN_NAME] = STATUS_FAILED
            _clear_downstream_from_termini_failure(pipeline_record)
            continue

        print(f"[FRUTON] amber_termini -> {pdb_id}")

        try:
            amber_termini_output_path, _termini_result_list = (
                convert_protein_termini_for_pdb_directory(
                    pdb_directory=pdb_dir,
                    pdb_id=pdb_id,
                    input_pdb_path=amber_input_path_text,
                )
            )
            pipeline_record[AMBER_TERMINI_STATUS_COLUMN_NAME] = STATUS_SUCCESS
            pipeline_record[AMBER_TERMINI_INPUT_PATH_COLUMN_NAME] = (
                amber_input_path_text
            )
            pipeline_record[AMBER_TERMINI_OUTPUT_PATH_COLUMN_NAME] = str(
                amber_termini_output_path
            )
        except Exception as error:
            print(f"[FRUTON] amber_termini failed for {pdb_id}: {error!r}")
            _clear_amber_termini_fields(pipeline_record)
            pipeline_record[AMBER_TERMINI_STATUS_COLUMN_NAME] = STATUS_FAILED
            _clear_downstream_from_termini_failure(pipeline_record)

    print("[FRUTON] Step 14: internal_capping")
    for pipeline_record in pipeline_record_list:
        pdb_id = pipeline_record[PDB_ID_COLUMN_NAME]
        pdb_dir = Path(pipeline_record[PDB_DIRECTORY_COLUMN_NAME])

        if pipeline_record.get(AMBER_TERMINI_STATUS_COLUMN_NAME, "") != STATUS_SUCCESS:
            print(f"[FRUTON] internal_capping skipped for {pdb_id}")
            _clear_internal_capping_fields(pipeline_record)
            pipeline_record[INTERNAL_CAPPING_STATUS_COLUMN_NAME] = STATUS_SKIPPED
            _clear_downstream_from_internal_capping_failure(pipeline_record)
            continue

        amber_termini_input_path = str(
            pipeline_record.get(AMBER_TERMINI_OUTPUT_PATH_COLUMN_NAME, "")
        ).strip()

        if not amber_termini_input_path or not Path(amber_termini_input_path).exists():
            print(f"[FRUTON] internal_capping skipped for {pdb_id}: missing input")
            _clear_internal_capping_fields(pipeline_record)
            pipeline_record[INTERNAL_CAPPING_STATUS_COLUMN_NAME] = STATUS_FAILED
            _clear_downstream_from_internal_capping_failure(pipeline_record)
            continue

        print(f"[FRUTON] internal_capping -> {pdb_id}")

        try:
            internal_capped_output_path, _summary_list = (
                cap_internal_gaps_for_pdb_directory(
                    pdb_directory=pdb_dir,
                    pdb_id=pdb_id,
                    input_pdb_path=amber_termini_input_path,
                )
            )
            pipeline_record[INTERNAL_CAPPING_STATUS_COLUMN_NAME] = STATUS_SUCCESS
            pipeline_record[INTERNAL_CAPPING_INPUT_PATH_COLUMN_NAME] = (
                amber_termini_input_path
            )
            pipeline_record[INTERNAL_CAPPING_OUTPUT_PATH_COLUMN_NAME] = str(
                internal_capped_output_path
            )
        except Exception as error:
            print(f"[FRUTON] internal_capping failed for {pdb_id}: {error!r}")
            _clear_internal_capping_fields(pipeline_record)
            pipeline_record[INTERNAL_CAPPING_STATUS_COLUMN_NAME] = STATUS_FAILED
            _clear_downstream_from_internal_capping_failure(pipeline_record)

    print("[FRUTON] Step 15: prepared_structure")
    for pipeline_record in pipeline_record_list:
        pdb_id = pipeline_record[PDB_ID_COLUMN_NAME]
        pdb_dir = Path(pipeline_record[PDB_DIRECTORY_COLUMN_NAME])

        if (
            pipeline_record.get(INTERNAL_CAPPING_STATUS_COLUMN_NAME, "")
            != STATUS_SUCCESS
        ):
            print(f"[FRUTON] prepared_structure skipped for {pdb_id}")
            _clear_prepared_structure_fields(pipeline_record)
            pipeline_record[PREPARED_STRUCTURE_STATUS_COLUMN_NAME] = STATUS_SKIPPED
            continue

        internal_capped_input_path = str(
            pipeline_record.get(INTERNAL_CAPPING_OUTPUT_PATH_COLUMN_NAME, "")
        ).strip()

        if (
            not internal_capped_input_path
            or not Path(internal_capped_input_path).exists()
        ):
            print(f"[FRUTON] prepared_structure skipped for {pdb_id}: missing input")
            _clear_prepared_structure_fields(pipeline_record)
            pipeline_record[PREPARED_STRUCTURE_STATUS_COLUMN_NAME] = STATUS_FAILED
            continue

        structure_variant = _determine_prepared_variant(pipeline_record)
        had_gaps = (
            str(pipeline_record.get(HAS_GAPS_COLUMN_NAME, "")).strip().lower() == "yes"
        )

        print(f"[FRUTON] prepared_structure -> {pdb_id}")

        try:
            prepared_summary = build_prepared_structure_for_pdb_directory(
                pdb_directory=pdb_dir,
                pdb_id=pdb_id,
                had_gaps=had_gaps,
                structure_variant=structure_variant,
                protein_input_path=internal_capped_input_path,
            )

            pipeline_record[PREPARED_DIRECTORY_COLUMN_NAME] = str(
                Path(prepared_summary.output_pdb_path).parent
            )
            pipeline_record[PREPARED_STRUCTURE_STATUS_COLUMN_NAME] = STATUS_SUCCESS
            pipeline_record[PREPARED_STRUCTURE_VARIANT_COLUMN_NAME] = (
                structure_variant if structure_variant is not None else "none"
            )
            pipeline_record[PREPARED_STRUCTURE_PROTEIN_INPUT_PATH_COLUMN_NAME] = (
                internal_capped_input_path
            )
            pipeline_record[PREPARED_STRUCTURE_OUTPUT_PATH_COLUMN_NAME] = str(
                prepared_summary.output_pdb_path
            )
        except Exception as error:
            print(f"[FRUTON] prepared_structure failed for {pdb_id}: {error!r}")
            _clear_prepared_structure_fields(pipeline_record)
            pipeline_record[PREPARED_STRUCTURE_STATUS_COLUMN_NAME] = STATUS_FAILED

    print("[FRUTON] Step 16: save_pipeline_json")
    save_pipeline_table(
        pipeline_record_list,
        pipeline_json_path,
    )

    print("[FRUTON] Step 17: write_pipeline_xlsx")
    write_pipeline_to_xlsx(
        pipeline_record_list,
        pipeline_xlsx_path,
    )

    print(f"[FRUTON] Pipeline JSON written to: {pipeline_json_path}")
    print(f"[FRUTON] Pipeline XLSX written to: {pipeline_xlsx_path}")
    print("[FRUTON] Finished")


def main() -> None:
    """
    Command line entry point.
    """
    run_pipeline()


if __name__ == "__main__":
    main()
