"""
run_pipeline.py

Main entry point for the protein preparation pipeline.

Responsibilities
----------------
- define working directory
- make the src directory importable
- run pipeline steps in order
- build pipeline state
- write final XLSX file
"""

import sys
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_ROOT_DIR = SCRIPT_DIR.parent
SRC_DIR = PROJECT_ROOT_DIR / "src"

if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))


from stack_protein_preparation.pdb_sync import (
    read_pdb_records_from_csv,
    sync_pdb_csv_and_directories,
)
from stack_protein_preparation.pipeline_xlsx import (
    create_pipeline_record_from_input_values,
    write_pipeline_to_xlsx,
)


def run_pipeline() -> None:
    """
    Run the current minimal pipeline.
    """
    protein_data_dir = PROJECT_ROOT_DIR / "data" / "proteins"
    csv_path = protein_data_dir / "pdb_ids.csv"
    xlsx_output_path = protein_data_dir / "pipeline.xlsx"

    print("[PIPELINE] Starting")

    print("[PIPELINE] Step 1: PDB sync")
    sync_pdb_csv_and_directories(protein_data_dir)

    print("[PIPELINE] Step 2: Read CSV")
    pdb_record_list = read_pdb_records_from_csv(csv_path)

    print(f"[PIPELINE] Loaded {len(pdb_record_list)} records")

    print("[PIPELINE] Step 3: Build pipeline state")
    pipeline_record_list: list[dict[str, str]] = []

    for pdb_record in pdb_record_list:
        pdb_id = pdb_record.get("pdb_id", "")
        residue_range = pdb_record.get("range", "")

        pipeline_record = create_pipeline_record_from_input_values(
            pdb_id,
            residue_range,
        )

        pipeline_record["pdb_sync_done"] = "success"
        pipeline_record_list.append(pipeline_record)

    print("[PIPELINE] Step 4: Write XLSX")
    write_pipeline_to_xlsx(
        pipeline_record_list,
        xlsx_output_path,
    )

    print(f"[PIPELINE] XLSX written to: {xlsx_output_path}")
    print("[PIPELINE] Finished")


def main() -> None:
    """
    Command line entry point.
    """
    run_pipeline()


if __name__ == "__main__":
    main()
