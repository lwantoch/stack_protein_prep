# /home/grheco/repositorios/stack_protein_prep/src/stack_protein_preparation/insertion_codes.py

"""
insertion_codes.py

Handle insertion codes in PDB files.

Responsibilities
----------------
- detect insertion codes
- remove them using pdb_delinsertion
- provide standardized status output

Design
------
- no structural modification beyond insertion handling
- returns structured result dict for pipeline integration

Important
---------
- does NOT perform renumbering
- does NOT depend on alignment/mapping
"""

from __future__ import annotations

import shutil
import subprocess
from pathlib import Path

STATUS_NONE = "none"
STATUS_SUCCESS = "success"
STATUS_FAILED = "failed"

COLOR_GREEN = "green"
COLOR_RED = "red"


def pdb_has_insertion_codes(pdb_path: Path) -> bool:
    """
    Detect insertion codes directly from a PDB file.

    In the PDB format, the insertion code is column 27 (1-based),
    i.e. line[26] in Python 0-based indexing.
    """
    with pdb_path.open("r", encoding="utf-8") as handle:
        for line in handle:
            if not (line.startswith("ATOM") or line.startswith("HETATM")):
                continue

            if len(line) < 27:
                continue

            insertion_code = line[26].strip()
            if insertion_code:
                return True

    return False


def run_pdb_delinsertion(
    input_pdb_path: Path,
    output_pdb_path: Path,
) -> tuple[bool, str]:
    """
    Run pdb_delinsertion on one file and write output to output_pdb_path.

    Returns
    -------
    tuple[bool, str]
        (success, message)
    """
    try:
        result = subprocess.run(
            ["pdb_delinsertion", str(input_pdb_path)],
            capture_output=True,
            text=True,
            check=False,
        )
    except FileNotFoundError:
        return False, "pdb_delinsertion executable not found in PATH"

    if result.returncode != 0:
        return False, result.stderr.strip() or "pdb_delinsertion failed"

    stdout = result.stdout.strip()
    if not stdout:
        return False, "pdb_delinsertion returned empty stdout"

    output_pdb_path.parent.mkdir(parents=True, exist_ok=True)
    output_pdb_path.write_text(result.stdout, encoding="utf-8")
    return True, "pdb_delinsertion completed"


def copy_pdb_to_output(
    input_pdb_path: Path,
    output_pdb_path: Path,
) -> tuple[bool, str]:
    """
    Copy input PDB to output location.
    """
    output_pdb_path.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(input_pdb_path, output_pdb_path)
    return True, "copied original pdb (no insertion codes)"


def process_pdb_for_delinsertion(
    input_pdb_path: Path,
    output_pdb_path: Path,
) -> dict[str, str | bool]:
    """
    Process one PDB:
    - if insertion codes exist -> run pdb_delinsertion
    - else -> copy file unchanged

    Status convention
    -----------------
    - none    -> green
    - success -> green
    - failed  -> red

    Returns
    -------
    dict[str, str | bool]
        Structured result dictionary.
    """
    had_insertion_codes = pdb_has_insertion_codes(input_pdb_path)

    if had_insertion_codes:
        success, message = run_pdb_delinsertion(
            input_pdb_path=input_pdb_path,
            output_pdb_path=output_pdb_path,
        )
        if success:
            status = STATUS_SUCCESS
            color = COLOR_GREEN
        else:
            status = STATUS_FAILED
            color = COLOR_RED
    else:
        success, message = copy_pdb_to_output(
            input_pdb_path=input_pdb_path,
            output_pdb_path=output_pdb_path,
        )
        status = STATUS_NONE
        color = COLOR_GREEN

    return {
        "input_pdb": str(input_pdb_path),
        "output_pdb": str(output_pdb_path),
        "had_insertion_codes": had_insertion_codes,
        "status": status,
        "color": color,
        "success": success,
        "message": message,
    }


def find_input_pdb_for_protein(protein_dir: Path) -> Path | None:
    """
    Find the input PDB for one protein directory.

    Search order:
    1. PDB directly inside data/proteins/<PDB_ID>/
    2. If none found, first PDB inside the first subdirectory

    This matches the current project layout where files like:
    data/proteins/1A5H/1A5H.pdb
    exist directly in the protein directory.
    """
    if not protein_dir.exists():
        return None

    direct_pdb_files = sorted(protein_dir.glob("*.pdb"))
    if direct_pdb_files:
        return direct_pdb_files[0]

    subdirs = sorted(path for path in protein_dir.iterdir() if path.is_dir())
    for subdir in subdirs:
        pdb_files = sorted(subdir.glob("*.pdb"))
        if pdb_files:
            return pdb_files[0]

    return None
