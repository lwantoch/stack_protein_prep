"""
/home/grheco/repositorios/stack_protein_prep/scripts/run_pipeline.py

Main entry point for the protein preparation pipeline.

Responsibilities
----------------
- define working directories
- make the src directory importable
- synchronize input CSV and protein directories
- build and merge pipeline records
- run implemented pipeline modules in order
- save pipeline state to JSON
- export pipeline overview to XLSX

Current implemented steps
-------------------------
1. pdb_sync
2. fasta_files
3. sequence_alignment
4. insertion_codes
5. component_split
6. gap_detection
7. filler

Important
---------
- pipeline state is persisted in pipeline.json
- pipeline overview is exported to pipeline.xlsx
- each module updates its corresponding pipeline status column
- component summary is stored as:
    - components_directory
    - has_metals
    - has_ligands
    - has_nonstandard_residues
- gap summary is stored as:
    - n_gaps
    - gap_sizes
- filler summary is stored as:
    - filler_directory
    - filler_model_path
    - filler_status
- UniProt ID is stored as:
    - uniprot_id
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
    COMPONENTS_DIRECTORY_COLUMN_NAME,
    FASTA_DIRECTORY_COLUMN_NAME,
    FASTA_FILES_DONE_COLUMN_NAME,
    FILLER_DIRECTORY_COLUMN_NAME,
    FILLER_MODEL_PATH_COLUMN_NAME,
    FILLER_STATUS_COLUMN_NAME,
    GAP_SIZES_COLUMN_NAME,
    HAS_LIGANDS_COLUMN_NAME,
    HAS_METALS_COLUMN_NAME,
    HAS_NONSTANDARD_RESIDUES_COLUMN_NAME,
    INSERTION_CODES_DONE_COLUMN_NAME,
    N_GAPS_COLUMN_NAME,
    PDB_DIRECTORY_COLUMN_NAME,
    PDB_ID_COLUMN_NAME,
    PDB_SYNC_DONE_COLUMN_NAME,
    RANGE_COLUMN_NAME,
    SEQUENCE_ALIGNMENT_DONE_COLUMN_NAME,
    STATUS_REQUIRED,
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
from stack_protein_preparation.sequence_alignment import (
    run_alignments_for_pdb_directory,
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

        pipeline_record = create_protein_record(
            pdb_id=pdb_id,
            residue_range=residue_range,
        )
        pipeline_record[PDB_DIRECTORY_COLUMN_NAME] = str(pdb_directory)
        pipeline_record[FASTA_DIRECTORY_COLUMN_NAME] = str(fasta_directory)
        pipeline_record[ALIGNMENT_DIRECTORY_COLUMN_NAME] = str(alignment_directory)
        pipeline_record[COMPONENTS_DIRECTORY_COLUMN_NAME] = str(components_directory)
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
        merged_record[PDB_SYNC_DONE_COLUMN_NAME] = STATUS_SUCCESS

        merged_record_list.append(merged_record)

    return merged_record_list


def _clear_component_fields(pipeline_record: dict[str, str]) -> None:
    """
    Clear component-related output fields.
    """
    pipeline_record[HAS_METALS_COLUMN_NAME] = ""
    pipeline_record[HAS_LIGANDS_COLUMN_NAME] = ""
    pipeline_record[HAS_NONSTANDARD_RESIDUES_COLUMN_NAME] = ""


def _clear_gap_fields(pipeline_record: dict[str, str]) -> None:
    """
    Clear gap-related output fields.
    """
    pipeline_record[N_GAPS_COLUMN_NAME] = ""
    pipeline_record[GAP_SIZES_COLUMN_NAME] = ""


def _clear_filler_fields(pipeline_record: dict[str, str]) -> None:
    """
    Clear filler-related output fields.
    """
    pipeline_record[FILLER_DIRECTORY_COLUMN_NAME] = ""
    pipeline_record[FILLER_MODEL_PATH_COLUMN_NAME] = ""
    pipeline_record[FILLER_STATUS_COLUMN_NAME] = ""


def _set_component_status_fields(
    pipeline_record: dict[str, str],
    component_summary: dict[str, object],
) -> None:
    """
    Translate boolean component flags into pipeline warning/success statuses.
    """
    pipeline_record[HAS_METALS_COLUMN_NAME] = (
        STATUS_WARNING
        if bool(component_summary.get("has_metals", False))
        else STATUS_SUCCESS
    )
    pipeline_record[HAS_LIGANDS_COLUMN_NAME] = (
        STATUS_WARNING
        if bool(component_summary.get("has_ligands", False))
        else STATUS_SUCCESS
    )
    pipeline_record[HAS_NONSTANDARD_RESIDUES_COLUMN_NAME] = (
        STATUS_WARNING
        if bool(component_summary.get("has_nonstandard_residues", False))
        else STATUS_SUCCESS
    )


def _find_template_pdb_for_filler(
    pdb_id: str,
    pdb_dir: Path,
    components_dir: Path,
) -> Path | None:
    """
    Find the template PDB to be used by the filler step.

    Current best-effort search order
    --------------------------------
    1. <components_dir>/<PDB_ID>_protein.pdb
    2. <pdb_dir>/<PDB_ID>_delins.pdb
    3. <pdb_dir>/<PDB_ID>.pdb
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


def _find_uniprot_id_for_protein(pdb_dir: Path) -> str:
    """
    Find UniProt ID from FASTA files in the protein directory.

    Expected pattern
    ----------------
    UniProt_<UNIPROT_ID>.fasta
    """
    fasta_dir = pdb_dir / "fasta"

    if not fasta_dir.exists():
        return ""

    for fasta_path in sorted(fasta_dir.glob("UniProt_*.fasta")):
        match = re.match(r"^UniProt_([A-Z0-9]+)\.fasta$", fasta_path.name)
        if match:
            return match.group(1)

    return ""


def run_pipeline() -> None:
    """
    Run the currently implemented protein preparation pipeline.
    """
    protein_data_dir = PROJECT_ROOT_DIR / "data" / "proteins"
    pdb_ids_csv_path = protein_data_dir / "pdb_ids.csv"
    pipeline_json_path = protein_data_dir / "pipeline.json"
    pipeline_xlsx_path = protein_data_dir / "pipeline.xlsx"

    print("[PIPELINE] Starting")

    print("[PIPELINE] Step 1: PDB sync")
    sync_pdb_csv_and_directories(protein_data_dir)

    print("[PIPELINE] Step 2: Read input CSV")
    pdb_record_list = read_pdb_records_from_csv(pdb_ids_csv_path)
    print(f"[PIPELINE] Loaded {len(pdb_record_list)} input records")

    print("[PIPELINE] Step 3: Load existing pipeline state")
    existing_pipeline_record_list = load_pipeline_table(pipeline_json_path)
    print(
        f"[PIPELINE] Loaded {len(existing_pipeline_record_list)} existing pipeline records"
    )

    print("[PIPELINE] Step 4: Build and merge pipeline state")
    new_pipeline_record_list = build_pipeline_records_from_input_csv(
        pdb_record_list=pdb_record_list,
        protein_data_dir=protein_data_dir,
    )
    pipeline_record_list = merge_existing_and_new_pipeline_records(
        existing_pipeline_record_list=existing_pipeline_record_list,
        new_pipeline_record_list=new_pipeline_record_list,
    )
    print(f"[PIPELINE] Active pipeline records: {len(pipeline_record_list)}")

    print("[PIPELINE] Step 5: Create FASTA files")
    for pipeline_record in pipeline_record_list:
        pdb_id = pipeline_record[PDB_ID_COLUMN_NAME]
        pdb_dir = Path(pipeline_record[PDB_DIRECTORY_COLUMN_NAME])

        print(f"[PIPELINE] fasta_files -> {pdb_id}")

        try:
            create_fasta_files_for_pdb_directory(pdb_dir)
            pipeline_record[FASTA_FILES_DONE_COLUMN_NAME] = STATUS_SUCCESS
            pipeline_record[UNIPROT_ID_COLUMN_NAME] = _find_uniprot_id_for_protein(
                pdb_dir
            )
        except Exception as error:
            print(f"[PIPELINE] fasta_files failed for {pdb_id}: {error}")
            pipeline_record[FASTA_FILES_DONE_COLUMN_NAME] = STATUS_REQUIRED
            pipeline_record[SEQUENCE_ALIGNMENT_DONE_COLUMN_NAME] = STATUS_REQUIRED
            pipeline_record[INSERTION_CODES_DONE_COLUMN_NAME] = STATUS_REQUIRED
            pipeline_record[UNIPROT_ID_COLUMN_NAME] = ""
            _clear_component_fields(pipeline_record)
            _clear_gap_fields(pipeline_record)
            _clear_filler_fields(pipeline_record)

    print("[PIPELINE] Step 6: Run sequence alignment")
    for pipeline_record in pipeline_record_list:
        pdb_id = pipeline_record[PDB_ID_COLUMN_NAME]
        pdb_dir = Path(pipeline_record[PDB_DIRECTORY_COLUMN_NAME])

        if pipeline_record.get(FASTA_FILES_DONE_COLUMN_NAME, "") != STATUS_SUCCESS:
            print(
                f"[PIPELINE] sequence_alignment skipped for {pdb_id} "
                "(FASTA step failed)"
            )
            pipeline_record[SEQUENCE_ALIGNMENT_DONE_COLUMN_NAME] = STATUS_REQUIRED
            pipeline_record[INSERTION_CODES_DONE_COLUMN_NAME] = STATUS_REQUIRED
            _clear_component_fields(pipeline_record)
            _clear_gap_fields(pipeline_record)
            _clear_filler_fields(pipeline_record)
            continue

        print(f"[PIPELINE] sequence_alignment -> {pdb_id}")

        try:
            run_alignments_for_pdb_directory(pdb_dir)
            pipeline_record[SEQUENCE_ALIGNMENT_DONE_COLUMN_NAME] = STATUS_SUCCESS

            if not str(pipeline_record.get(UNIPROT_ID_COLUMN_NAME, "")).strip():
                pipeline_record[UNIPROT_ID_COLUMN_NAME] = _find_uniprot_id_for_protein(
                    pdb_dir
                )
        except Exception as error:
            print(f"[PIPELINE] sequence_alignment failed for {pdb_id}: {error}")
            pipeline_record[SEQUENCE_ALIGNMENT_DONE_COLUMN_NAME] = STATUS_REQUIRED
            pipeline_record[INSERTION_CODES_DONE_COLUMN_NAME] = STATUS_REQUIRED
            _clear_component_fields(pipeline_record)
            _clear_gap_fields(pipeline_record)
            _clear_filler_fields(pipeline_record)

    print("[PIPELINE] Step 7: Handle insertion codes")
    for pipeline_record in pipeline_record_list:
        pdb_id = pipeline_record[PDB_ID_COLUMN_NAME]
        pdb_dir = Path(pipeline_record[PDB_DIRECTORY_COLUMN_NAME])

        if (
            pipeline_record.get(SEQUENCE_ALIGNMENT_DONE_COLUMN_NAME, "")
            != STATUS_SUCCESS
        ):
            print(
                f"[PIPELINE] insertion_codes skipped for {pdb_id} "
                "(sequence alignment step failed)"
            )
            pipeline_record[INSERTION_CODES_DONE_COLUMN_NAME] = STATUS_REQUIRED
            _clear_component_fields(pipeline_record)
            _clear_gap_fields(pipeline_record)
            _clear_filler_fields(pipeline_record)
            continue

        print(f"[PIPELINE] insertion_codes -> {pdb_id}")

        input_pdb_path = find_input_pdb_for_protein(pdb_dir)
        if input_pdb_path is None:
            print(
                f"[PIPELINE] insertion_codes failed for {pdb_id}: input PDB not found"
            )
            pipeline_record[INSERTION_CODES_DONE_COLUMN_NAME] = STATUS_REQUIRED
            _clear_component_fields(pipeline_record)
            _clear_gap_fields(pipeline_record)
            _clear_filler_fields(pipeline_record)
            continue

        output_pdb_path = pdb_dir / f"{pdb_id}_delins.pdb"

        try:
            insertion_result = process_pdb_for_delinsertion(
                input_pdb_path=input_pdb_path,
                output_pdb_path=output_pdb_path,
            )
        except Exception as error:
            print(f"[PIPELINE] insertion_codes failed for {pdb_id}: {error}")
            pipeline_record[INSERTION_CODES_DONE_COLUMN_NAME] = STATUS_REQUIRED
            _clear_component_fields(pipeline_record)
            _clear_gap_fields(pipeline_record)
            _clear_filler_fields(pipeline_record)
            continue

        insertion_status = str(insertion_result.get("status", "")).strip().lower()

        if insertion_status in {INSERTION_STATUS_NONE, INSERTION_STATUS_SUCCESS}:
            pipeline_record[INSERTION_CODES_DONE_COLUMN_NAME] = STATUS_SUCCESS
        elif insertion_status == INSERTION_STATUS_FAILED:
            pipeline_record[INSERTION_CODES_DONE_COLUMN_NAME] = STATUS_REQUIRED
            _clear_component_fields(pipeline_record)
            _clear_gap_fields(pipeline_record)
            _clear_filler_fields(pipeline_record)
        else:
            pipeline_record[INSERTION_CODES_DONE_COLUMN_NAME] = STATUS_REQUIRED
            _clear_component_fields(pipeline_record)
            _clear_gap_fields(pipeline_record)
            _clear_filler_fields(pipeline_record)

    print("[PIPELINE] Step 8: Split components")
    for pipeline_record in pipeline_record_list:
        pdb_id = pipeline_record[PDB_ID_COLUMN_NAME]
        pdb_dir = Path(pipeline_record[PDB_DIRECTORY_COLUMN_NAME])
        components_dir = Path(pipeline_record[COMPONENTS_DIRECTORY_COLUMN_NAME])

        if pipeline_record.get(INSERTION_CODES_DONE_COLUMN_NAME, "") != STATUS_SUCCESS:
            print(
                f"[PIPELINE] component_split skipped for {pdb_id} "
                "(insertion_codes step failed)"
            )
            _clear_component_fields(pipeline_record)
            _clear_gap_fields(pipeline_record)
            _clear_filler_fields(pipeline_record)
            continue

        component_input_pdb_path = pdb_dir / f"{pdb_id}_delins.pdb"

        if not component_input_pdb_path.exists():
            print(
                f"[PIPELINE] component_split skipped for {pdb_id}: "
                f"missing file {component_input_pdb_path}"
            )
            _clear_component_fields(pipeline_record)
            _clear_gap_fields(pipeline_record)
            _clear_filler_fields(pipeline_record)
            continue

        print(f"[PIPELINE] component_split -> {pdb_id}")

        try:
            component_summary = split_pdb_components(
                pdb_path=component_input_pdb_path,
                output_dir=components_dir,
                protein_stem=pdb_id,
            )
        except Exception as error:
            print(f"[PIPELINE] component_split failed for {pdb_id}: {error}")
            _clear_component_fields(pipeline_record)
            _clear_gap_fields(pipeline_record)
            _clear_filler_fields(pipeline_record)
            continue

        pipeline_record[COMPONENTS_DIRECTORY_COLUMN_NAME] = str(
            component_summary["output_dir"]
        )
        _set_component_status_fields(pipeline_record, component_summary)

        print(
            f"[PIPELINE] component_split result for {pdb_id}: "
            f"metals={pipeline_record[HAS_METALS_COLUMN_NAME]}, "
            f"ligands={pipeline_record[HAS_LIGANDS_COLUMN_NAME]}, "
            f"nonstandard={pipeline_record[HAS_NONSTANDARD_RESIDUES_COLUMN_NAME]}, "
            f"components_dir={pipeline_record[COMPONENTS_DIRECTORY_COLUMN_NAME]}"
        )

    print("[PIPELINE] Step 9: Detect gaps")
    for pipeline_record in pipeline_record_list:
        pdb_id = pipeline_record[PDB_ID_COLUMN_NAME]
        components_dir = Path(pipeline_record[COMPONENTS_DIRECTORY_COLUMN_NAME])

        if pipeline_record.get(INSERTION_CODES_DONE_COLUMN_NAME, "") != STATUS_SUCCESS:
            print(
                f"[PIPELINE] gap_detection skipped for {pdb_id} "
                "(insertion_codes step failed)"
            )
            _clear_gap_fields(pipeline_record)
            _clear_filler_fields(pipeline_record)
            continue

        gap_input_pdb_path = components_dir / f"{pdb_id}_protein.pdb"

        if not gap_input_pdb_path.exists():
            print(
                f"[PIPELINE] gap_detection skipped for {pdb_id}: "
                f"missing file {gap_input_pdb_path}"
            )
            _clear_gap_fields(pipeline_record)
            _clear_filler_fields(pipeline_record)
            continue

        print(f"[PIPELINE] gap_detection -> {pdb_id}")

        try:
            gap_summary = summarize_gaps(gap_input_pdb_path)
        except Exception as error:
            print(f"[PIPELINE] gap_detection failed for {pdb_id}: {error}")
            _clear_gap_fields(pipeline_record)
            _clear_filler_fields(pipeline_record)
            continue

        n_gaps = int(gap_summary.get("n_gaps", 0))
        gap_sizes = gap_summary.get("gap_sizes", [])

        pipeline_record[N_GAPS_COLUMN_NAME] = str(n_gaps)

        if n_gaps == 0:
            pipeline_record[GAP_SIZES_COLUMN_NAME] = "none"
        else:
            pipeline_record[GAP_SIZES_COLUMN_NAME] = "|".join(
                str(gap_size) for gap_size in gap_sizes
            )

        print(
            f"[PIPELINE] gap_detection result for {pdb_id}: "
            f"n_gaps={pipeline_record[N_GAPS_COLUMN_NAME]}, "
            f"gap_sizes={pipeline_record[GAP_SIZES_COLUMN_NAME]!r}"
        )

    print("[PIPELINE] Step 10: Run filler")
    for pipeline_record in pipeline_record_list:
        pdb_id = pipeline_record[PDB_ID_COLUMN_NAME]
        pdb_dir = Path(pipeline_record[PDB_DIRECTORY_COLUMN_NAME])
        alignment_dir = Path(pipeline_record[ALIGNMENT_DIRECTORY_COLUMN_NAME])
        components_dir = Path(pipeline_record[COMPONENTS_DIRECTORY_COLUMN_NAME])

        if pipeline_record.get(INSERTION_CODES_DONE_COLUMN_NAME, "") != STATUS_SUCCESS:
            print(
                f"[PIPELINE] filler skipped for {pdb_id} (insertion_codes step failed)"
            )
            _clear_filler_fields(pipeline_record)
            continue

        if not alignment_dir.exists():
            print(
                f"[PIPELINE] filler skipped for {pdb_id}: "
                f"alignment directory not found: {alignment_dir}"
            )
            _clear_filler_fields(pipeline_record)
            continue

        template_pdb_path = _find_template_pdb_for_filler(
            pdb_id=pdb_id,
            pdb_dir=pdb_dir,
            components_dir=components_dir,
        )

        if template_pdb_path is None:
            print(f"[PIPELINE] filler skipped for {pdb_id}: no template PDB found")
            _clear_filler_fields(pipeline_record)
            continue

        n_gaps_text = str(pipeline_record.get(N_GAPS_COLUMN_NAME, "")).strip()

        if not n_gaps_text:
            print(
                f"[PIPELINE] filler skipped for {pdb_id} (gap_detection result missing)"
            )
            _clear_filler_fields(pipeline_record)
            continue

        if n_gaps_text == "0":
            print(
                f"[PIPELINE] filler skipped for {pdb_id} (no structural gaps detected)"
            )
            _clear_filler_fields(pipeline_record)
            pipeline_record[FILLER_STATUS_COLUMN_NAME] = STATUS_SUCCESS
            continue

        print(f"[PIPELINE] filler -> {pdb_id}")

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
            print(f"[PIPELINE] filler failed for {pdb_id}: {error!r}")
            traceback.print_exc()
            _clear_filler_fields(pipeline_record)
            continue

        pipeline_record[FILLER_DIRECTORY_COLUMN_NAME] = str(filler_result.output_dir)

        if filler_result.final_model_path is not None:
            pipeline_record[FILLER_MODEL_PATH_COLUMN_NAME] = str(
                filler_result.final_model_path
            )
        else:
            pipeline_record[FILLER_MODEL_PATH_COLUMN_NAME] = ""

        if filler_result.skipped:
            if filler_result.fill_decision.alphafold_candidate:
                pipeline_record[FILLER_STATUS_COLUMN_NAME] = STATUS_WARNING
            else:
                pipeline_record[FILLER_STATUS_COLUMN_NAME] = STATUS_SUCCESS
        else:
            if filler_result.fill_decision.overall_classification == "yellow":
                pipeline_record[FILLER_STATUS_COLUMN_NAME] = STATUS_WARNING
            else:
                pipeline_record[FILLER_STATUS_COLUMN_NAME] = STATUS_SUCCESS

        print(
            f"[PIPELINE] filler result for {pdb_id}: "
            f"filler_directory={pipeline_record[FILLER_DIRECTORY_COLUMN_NAME]!r}, "
            f"filler_model_path={pipeline_record[FILLER_MODEL_PATH_COLUMN_NAME]!r}, "
            f"filler_status={pipeline_record[FILLER_STATUS_COLUMN_NAME]!r}"
        )

    print("[PIPELINE] Step 11: Save pipeline JSON")
    save_pipeline_table(
        pipeline_record_list,
        pipeline_json_path,
    )

    print("[PIPELINE] Step 12: Write pipeline XLSX")
    write_pipeline_to_xlsx(
        pipeline_record_list,
        pipeline_xlsx_path,
    )

    print(f"[PIPELINE] Pipeline JSON written to: {pipeline_json_path}")
    print(f"[PIPELINE] Pipeline XLSX written to: {pipeline_xlsx_path}")
    print("[PIPELINE] Finished")


def main() -> None:
    """
    Command line entry point.
    """
    run_pipeline()


if __name__ == "__main__":
    main()
