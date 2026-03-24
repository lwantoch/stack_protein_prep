# /home/grheco/repositorios/stack_protein_prep/tests/test_amber_renaming.py

from __future__ import annotations

from pathlib import Path

from stack_protein_preparation.amber_renaming import (
    amber_rename_protein_structure,
)


def _write_pdb(path: Path, content: str) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(content, encoding="utf-8")
    return path


def test_his_renaming_hid(tmp_path: Path) -> None:
    pdb_id = "TEST"
    protein_dir = tmp_path / pdb_id
    components_dir = protein_dir / "components"

    # HIS with only HD1 → HID
    pdb_content = """\
ATOM      1  N   HIS A   1      0.000   0.000   0.000
ATOM      2  CA  HIS A   1      0.000   0.000   0.000
ATOM      3  HD1 HIS A   1      0.000   0.000   0.000
END
"""
    input_path = _write_pdb(
        components_dir / f"{pdb_id}_proteinH.pdb",
        pdb_content,
    )

    result = amber_rename_protein_structure(
        pdb_id=pdb_id,
        protein_dir=protein_dir,
    )

    output_path = Path(result["amber_output_path"])

    assert result["amber_renaming_success"] is True
    assert output_path.exists()

    out_text = output_path.read_text()
    assert "HID" in out_text
    assert result["his_to_hid"] == 1
    assert result["his_to_hie"] == 0
    assert result["his_to_hip"] == 0


def test_his_renaming_hie(tmp_path: Path) -> None:
    pdb_id = "TEST2"
    protein_dir = tmp_path / pdb_id
    components_dir = protein_dir / "components"

    # HIS with only HE2 → HIE
    pdb_content = """\
ATOM      1  N   HIS A   1      0.000   0.000   0.000
ATOM      2  CA  HIS A   1      0.000   0.000   0.000
ATOM      3  HE2 HIS A   1      0.000   0.000   0.000
END
"""
    _write_pdb(
        components_dir / f"{pdb_id}_proteinH.pdb",
        pdb_content,
    )

    result = amber_rename_protein_structure(
        pdb_id=pdb_id,
        protein_dir=protein_dir,
    )

    output_path = Path(result["amber_output_path"])

    out_text = output_path.read_text()
    assert "HIE" in out_text
    assert result["his_to_hie"] == 1


def test_his_renaming_hip(tmp_path: Path) -> None:
    pdb_id = "TEST3"
    protein_dir = tmp_path / pdb_id
    components_dir = protein_dir / "components"

    # HIS with both → HIP
    pdb_content = """\
ATOM      1  N   HIS A   1      0.000   0.000   0.000
ATOM      2  CA  HIS A   1      0.000   0.000   0.000
ATOM      3  HD1 HIS A   1      0.000   0.000   0.000
ATOM      4  HE2 HIS A   1      0.000   0.000   0.000
END
"""
    _write_pdb(
        components_dir / f"{pdb_id}_proteinH.pdb",
        pdb_content,
    )

    result = amber_rename_protein_structure(
        pdb_id=pdb_id,
        protein_dir=protein_dir,
    )

    output_path = Path(result["amber_output_path"])

    out_text = output_path.read_text()
    assert "HIP" in out_text
    assert result["his_to_hip"] == 1


def test_output_files_created(tmp_path: Path) -> None:
    pdb_id = "TEST4"
    protein_dir = tmp_path / pdb_id
    components_dir = protein_dir / "components"

    pdb_content = """\
ATOM      1  N   HIS A   1      0.000   0.000   0.000
ATOM      2  CA  HIS A   1      0.000   0.000   0.000
ATOM      3  HD1 HIS A   1      0.000   0.000   0.000
END
"""
    _write_pdb(
        components_dir / f"{pdb_id}_proteinH.pdb",
        pdb_content,
    )

    result = amber_rename_protein_structure(
        pdb_id=pdb_id,
        protein_dir=protein_dir,
    )

    output_path = Path(result["amber_output_path"])
    stats_path = Path(result["amber_stats_json_path"])

    assert output_path.exists()
    assert stats_path.exists()
