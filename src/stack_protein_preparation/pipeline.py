"""
pipeline.py

Central pipeline execution module.

This file defines the main order of pipeline steps.

Current design
--------------
- The pipeline is executed step by step in a fixed order.
- Each step is implemented in a separate module inside 'steps/'.
- This file only coordinates those steps.

Important
---------
This module contains NO command line interface logic.
It should be importable and testable.
"""

from pathlib import Path

from stack_protein_preparation.steps.pdb_sync import (
    sync_pdb_csv_and_directories,
)


def run_pipeline(protein_data_dir: Path) -> None:
    """
    Run the full pipeline.

    Parameters
    ----------
    protein_data_dir
        Directory that contains:
        - the CSV file with PDB IDs
        - one subdirectory per protein
    """
    print("[PIPELINE] Starting pipeline")

    print(f"[PIPELINE] Using data directory: {protein_data_dir}")

    # --- Step 1: PDB synchronization ---
    print("[PIPELINE] Step 1: PDB synchronization")

    sync_pdb_csv_and_directories(protein_data_dir)

    print("[PIPELINE] Pipeline finished successfully")
