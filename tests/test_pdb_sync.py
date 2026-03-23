"""
test_pdb_sync.py

Tests for the pdb_sync module.

Design
------
- CSV files are stored in tests/testdata/
- directories are created inside tmp_path
- imports follow the CURRENT project structure
"""

import shutil
from pathlib import Path

from stack_protein_preparation.pdb_sync import (
    read_pdb_records_from_csv,
    sync_pdb_csv_and_directories,
)


def test_read_csv_with_range() -> None:
    """
    Test reading a CSV file with a range column.
    """
    csv_path = Path("tests/testdata/pdb_ids_with_range.csv")

    records = read_pdb_records_from_csv(csv_path)

    assert len(records) == 3

    assert records[0]["pdb_id"] == "2AFX"
    assert records[0]["range"] == "10-200"

    assert records[1]["pdb_id"] == "1B8O"
    assert records[1]["range"] == ""

    assert records[2]["pdb_id"] == "2W8Y"
    assert records[2]["range"] == "5-150"


def test_sync_csv_only(tmp_path: Path) -> None:
    """
    Test case:
    - CSV exists
    - no protein directories exist yet

    Expected behavior:
    - protein directories are created
    """
    source_csv_path = Path("tests/testdata/pdb_ids_with_range.csv")
    target_protein_data_dir = tmp_path / "proteins"

    target_protein_data_dir.mkdir()
    shutil.copy(source_csv_path, target_protein_data_dir / "pdb_ids.csv")

    sync_pdb_csv_and_directories(target_protein_data_dir)

    assert (target_protein_data_dir / "2AFX").exists()
    assert (target_protein_data_dir / "1B8O").exists()
    assert (target_protein_data_dir / "2W8Y").exists()


def test_sync_dirs_only(tmp_path: Path) -> None:
    """
    Test case:
    - protein directories exist
    - no CSV exists

    Expected behavior:
    - CSV is created
    """
    target_protein_data_dir = tmp_path / "proteins"
    target_protein_data_dir.mkdir()

    (target_protein_data_dir / "2AFX").mkdir()
    (target_protein_data_dir / "1B8O").mkdir()

    sync_pdb_csv_and_directories(target_protein_data_dir)

    csv_path = target_protein_data_dir / "pdb_ids.csv"

    assert csv_path.exists()

    csv_content = csv_path.read_text()

    assert "2AFX" in csv_content
    assert "1B8O" in csv_content


def test_sync_mismatch(tmp_path: Path) -> None:
    """
    Test case:
    - CSV and directories differ

    Expected behavior:
    - final CSV contains IDs from both sources
    """
    source_csv_path = Path("tests/testdata/pdb_ids_with_range.csv")
    target_protein_data_dir = tmp_path / "proteins"

    target_protein_data_dir.mkdir()
    shutil.copy(source_csv_path, target_protein_data_dir / "pdb_ids.csv")

    (target_protein_data_dir / "2AFX").mkdir()
    (target_protein_data_dir / "3DEF").mkdir()

    sync_pdb_csv_and_directories(target_protein_data_dir)

    final_csv_path = target_protein_data_dir / "pdb_ids.csv"
    final_csv_content = final_csv_path.read_text()

    assert "2AFX" in final_csv_content
    assert "1B8O" in final_csv_content
    assert "2W8Y" in final_csv_content
    assert "3DEF" in final_csv_content
