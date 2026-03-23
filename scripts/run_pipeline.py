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
5. gap_detection

Important
---------
- pipeline state is persisted in pipeline.json
- pipeline overview is exported to pipeline.xlsx
- each module updates its corresponding pipeline status column
- gap information is stored as:
    - n_gaps
    - gap_sizes
"""

from __future__ import annotations

import sys
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_ROOT_DIR = SCRIPT_DIR.parent
SRC_DIR = PROJECT_ROOT_DIR / "src"

if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

from stack_protein_preparation.fasta_files import (
    create_fasta_files_for_pdb_directory,
)
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
from stack_protein_preparation.pdb_sync import (
    read_pdb_records_from_csv,
    sync_pdb_csv_and_directories,
)
from stack_protein_preparation.pipeline_state import (
    ALIGNMENT_DIRECTORY_COLUMN_NAME,
    FASTA_DIRECTORY_COLUMN_NAME,
    FASTA_FILES_DONE_COLUMN_NAME,
    GAP_SIZES_COLUMN_NAME,
    INSERTION_CODES_DONE_COLUMN_NAME,
    N_GAPS_COLUMN_NAME,
    PDB_DIRECTORY_COLUMN_NAME,
    PDB_ID_COLUMN_NAME,
    PDB_SYNC_DONE_COLUMN_NAME,
    RANGE_COLUMN_NAME,
    SEQUENCE_ALIGNMENT_DONE_COLUMN_NAME,
    STATUS_REQUIRED,
    STATUS_SUCCESS,
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

        pipeline_record = create_protein_record(
            pdb_id=pdb_id,
            residue_range=residue_range,
        )
        pipeline_record[PDB_DIRECTORY_COLUMN_NAME] = str(pdb_directory)
        pipeline_record[FASTA_DIRECTORY_COLUMN_NAME] = str(fasta_directory)
        pipeline_record[ALIGNMENT_DIRECTORY_COLUMN_NAME] = str(alignment_directory)
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
        merged_record[PDB_SYNC_DONE_COLUMN_NAME] = STATUS_SUCCESS

        merged_record_list.append(merged_record)

    return merged_record_list


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
        except Exception as error:
            print(f"[PIPELINE] fasta_files failed for {pdb_id}: {error}")
            pipeline_record[FASTA_FILES_DONE_COLUMN_NAME] = STATUS_REQUIRED
            pipeline_record[SEQUENCE_ALIGNMENT_DONE_COLUMN_NAME] = STATUS_REQUIRED
            pipeline_record[INSERTION_CODES_DONE_COLUMN_NAME] = STATUS_REQUIRED
            pipeline_record[N_GAPS_COLUMN_NAME] = ""
            pipeline_record[GAP_SIZES_COLUMN_NAME] = ""

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
            pipeline_record[N_GAPS_COLUMN_NAME] = ""
            pipeline_record[GAP_SIZES_COLUMN_NAME] = ""
            continue

        print(f"[PIPELINE] sequence_alignment -> {pdb_id}")

        try:
            run_alignments_for_pdb_directory(pdb_dir)
            pipeline_record[SEQUENCE_ALIGNMENT_DONE_COLUMN_NAME] = STATUS_SUCCESS
        except Exception as error:
            print(f"[PIPELINE] sequence_alignment failed for {pdb_id}: {error}")
            pipeline_record[SEQUENCE_ALIGNMENT_DONE_COLUMN_NAME] = STATUS_REQUIRED
            pipeline_record[INSERTION_CODES_DONE_COLUMN_NAME] = STATUS_REQUIRED
            pipeline_record[N_GAPS_COLUMN_NAME] = ""
            pipeline_record[GAP_SIZES_COLUMN_NAME] = ""

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
            pipeline_record[N_GAPS_COLUMN_NAME] = ""
            pipeline_record[GAP_SIZES_COLUMN_NAME] = ""
            continue

        print(f"[PIPELINE] insertion_codes -> {pdb_id}")

        input_pdb_path = find_input_pdb_for_protein(pdb_dir)
        if input_pdb_path is None:
            print(
                f"[PIPELINE] insertion_codes failed for {pdb_id}: input PDB not found"
            )
            pipeline_record[INSERTION_CODES_DONE_COLUMN_NAME] = STATUS_REQUIRED
            pipeline_record[N_GAPS_COLUMN_NAME] = ""
            pipeline_record[GAP_SIZES_COLUMN_NAME] = ""
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
            pipeline_record[N_GAPS_COLUMN_NAME] = ""
            pipeline_record[GAP_SIZES_COLUMN_NAME] = ""
            continue

        insertion_status = str(insertion_result.get("status", "")).strip().lower()

        if insertion_status in {INSERTION_STATUS_NONE, INSERTION_STATUS_SUCCESS}:
            pipeline_record[INSERTION_CODES_DONE_COLUMN_NAME] = STATUS_SUCCESS
        elif insertion_status == INSERTION_STATUS_FAILED:
            pipeline_record[INSERTION_CODES_DONE_COLUMN_NAME] = STATUS_REQUIRED
            pipeline_record[N_GAPS_COLUMN_NAME] = ""
            pipeline_record[GAP_SIZES_COLUMN_NAME] = ""
        else:
            pipeline_record[INSERTION_CODES_DONE_COLUMN_NAME] = STATUS_REQUIRED
            pipeline_record[N_GAPS_COLUMN_NAME] = ""
            pipeline_record[GAP_SIZES_COLUMN_NAME] = ""

    print("[PIPELINE] Step 8: Detect gaps")
    for pipeline_record in pipeline_record_list:
        pdb_id = pipeline_record[PDB_ID_COLUMN_NAME]
        pdb_dir = Path(pipeline_record[PDB_DIRECTORY_COLUMN_NAME])

        if pipeline_record.get(INSERTION_CODES_DONE_COLUMN_NAME, "") != STATUS_SUCCESS:
            print(
                f"[PIPELINE] gap_detection skipped for {pdb_id} "
                "(insertion_codes step failed)"
            )
            pipeline_record[N_GAPS_COLUMN_NAME] = ""
            pipeline_record[GAP_SIZES_COLUMN_NAME] = ""
            continue

        gap_input_pdb_path = pdb_dir / f"{pdb_id}_delins.pdb"

        if not gap_input_pdb_path.exists():
            print(
                f"[PIPELINE] gap_detection skipped for {pdb_id}: "
                f"missing file {gap_input_pdb_path}"
            )
            pipeline_record[N_GAPS_COLUMN_NAME] = ""
            pipeline_record[GAP_SIZES_COLUMN_NAME] = ""
            continue

        print(f"[PIPELINE] gap_detection -> {pdb_id}")

        try:
            gap_summary = summarize_gaps(gap_input_pdb_path)
        except Exception as error:
            print(f"[PIPELINE] gap_detection failed for {pdb_id}: {error}")
            pipeline_record[N_GAPS_COLUMN_NAME] = ""
            pipeline_record[GAP_SIZES_COLUMN_NAME] = ""
            continue

        n_gaps = int(gap_summary.get("n_gaps", 0))
        gap_sizes = gap_summary.get("gap_sizes", [])

        pipeline_record[N_GAPS_COLUMN_NAME] = str(n_gaps)
        pipeline_record[GAP_SIZES_COLUMN_NAME] = "|".join(
            str(gap_size) for gap_size in gap_sizes
        )

        print(
            f"[PIPELINE] gap_detection result for {pdb_id}: "
            f"n_gaps={pipeline_record[N_GAPS_COLUMN_NAME]}, "
            f"gap_sizes={pipeline_record[GAP_SIZES_COLUMN_NAME]!r}"
        )

    print("[PIPELINE] Step 9: Save pipeline JSON")
    save_pipeline_table(
        pipeline_record_list,
        pipeline_json_path,
    )

    print("[PIPELINE] Step 10: Write pipeline XLSX")
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
