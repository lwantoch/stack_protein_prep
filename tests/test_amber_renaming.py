from __future__ import annotations

import json
from pathlib import Path

import pytest
from Bio.PDB import PDBParser

from stack_protein_preparation.amber_renaming import (
    _atom_names,
    _distance,
    _iter_atoms_altloc_ok,
    amber_rename_protein_structure,
    amber_rename_selected_structure,
    amber_rename_variant_structure,
    build_standard_amber_output_path,
    build_variant_amber_output_path,
    rename_structure_by_hydrogens,
    sanitize_variant_label,
    write_amber_renamed_pdb,
)


def _make_atom_line(
    serial: int,
    atom_name: str,
    residue_name: str,
    chain_id: str,
    residue_number: int,
    x: float,
    y: float,
    z: float,
    element: str,
    *,
    altloc: str = " ",
    record_name: str = "ATOM",
) -> str:
    return (
        f"{record_name:<6}"
        f"{serial:5d} "
        f"{atom_name:>4}"
        f"{altloc:1}"
        f"{residue_name:>3} "
        f"{chain_id:1}"
        f"{residue_number:4d}    "
        f"{x:8.3f}{y:8.3f}{z:8.3f}"
        f"{1.00:6.2f}{20.00:6.2f}          "
        f"{element:>2}"
    )


def _write_pdb(path: Path, lines: list[str]) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def _parse_structure(pdb_path: Path):
    parser = PDBParser(QUIET=True)
    return parser.get_structure("prot", str(pdb_path))


def _parse_first_model(pdb_path: Path):
    structure = _parse_structure(pdb_path)
    return next(structure.get_models())


def _parse_first_chain(pdb_path: Path):
    model = _parse_first_model(pdb_path)
    return next(model.get_chains())


def _residue_by_resseq(chain, resseq: int):
    for residue in chain.get_residues():
        if residue.id[1] == resseq:
            return residue
    raise KeyError(resseq)


def _residue_names(chain) -> list[str]:
    return [residue.get_resname().strip() for residue in chain.get_residues()]


def _build_histidine_pdb(path: Path, *, his_variant: str) -> Path:
    """
    his_variant:
    - "hid" -> HD1 present
    - "hie" -> HE2 present
    - "hip" -> HD1 and HE2 present
    - "none" -> neither present
    """
    lines = [
        _make_atom_line(1, "N", "HIS", "A", 1, 0.000, 0.000, 0.000, "N"),
        _make_atom_line(2, "CA", "HIS", "A", 1, 1.460, 0.000, 0.000, "C"),
        _make_atom_line(3, "C", "HIS", "A", 1, 2.900, 0.000, 0.000, "C"),
        _make_atom_line(4, "O", "HIS", "A", 1, 3.900, 0.000, 0.000, "O"),
        _make_atom_line(5, "ND1", "HIS", "A", 1, 1.200, 1.100, 0.000, "N"),
        _make_atom_line(6, "NE2", "HIS", "A", 1, 2.200, 1.100, 0.000, "N"),
    ]

    serial = 7
    if his_variant in {"hid", "hip"}:
        lines.append(
            _make_atom_line(serial, "HD1", "HIS", "A", 1, 1.200, 2.000, 0.000, "H")
        )
        serial += 1

    if his_variant in {"hie", "hip"}:
        lines.append(
            _make_atom_line(serial, "HE2", "HIS", "A", 1, 2.200, 2.000, 0.000, "H")
        )

    lines += ["TER", "END"]
    return _write_pdb(path, lines)


def _build_asp_glu_pdb(path: Path) -> Path:
    lines = [
        _make_atom_line(1, "N", "ASP", "A", 1, 0.000, 0.000, 0.000, "N"),
        _make_atom_line(2, "CA", "ASP", "A", 1, 1.460, 0.000, 0.000, "C"),
        _make_atom_line(3, "C", "ASP", "A", 1, 2.900, 0.000, 0.000, "C"),
        _make_atom_line(4, "O", "ASP", "A", 1, 3.900, 0.000, 0.000, "O"),
        _make_atom_line(5, "OD1", "ASP", "A", 1, 1.100, 1.200, 0.000, "O"),
        _make_atom_line(6, "OD2", "ASP", "A", 1, 2.100, 1.200, 0.000, "O"),
        _make_atom_line(7, "HD2", "ASP", "A", 1, 2.100, 2.100, 0.000, "H"),
        _make_atom_line(8, "N", "GLU", "A", 2, 4.300, 0.000, 0.000, "N"),
        _make_atom_line(9, "CA", "GLU", "A", 2, 5.760, 0.000, 0.000, "C"),
        _make_atom_line(10, "C", "GLU", "A", 2, 7.200, 0.000, 0.000, "C"),
        _make_atom_line(11, "O", "GLU", "A", 2, 8.200, 0.000, 0.000, "O"),
        _make_atom_line(12, "OE1", "GLU", "A", 2, 5.500, 1.200, 0.000, "O"),
        _make_atom_line(13, "OE2", "GLU", "A", 2, 6.500, 1.200, 0.000, "O"),
        _make_atom_line(14, "HE2", "GLU", "A", 2, 6.500, 2.100, 0.000, "H"),
        "TER",
        "END",
    ]
    return _write_pdb(path, lines)


def _build_cysteine_pdb(path: Path, *, with_hg: bool) -> Path:
    lines = [
        _make_atom_line(1, "N", "CYS", "A", 1, 0.000, 0.000, 0.000, "N"),
        _make_atom_line(2, "CA", "CYS", "A", 1, 1.460, 0.000, 0.000, "C"),
        _make_atom_line(3, "C", "CYS", "A", 1, 2.900, 0.000, 0.000, "C"),
        _make_atom_line(4, "O", "CYS", "A", 1, 3.900, 0.000, 0.000, "O"),
        _make_atom_line(5, "SG", "CYS", "A", 1, 1.460, 1.800, 0.000, "S"),
    ]
    if with_hg:
        lines.append(_make_atom_line(6, "HG", "CYS", "A", 1, 1.460, 2.800, 0.000, "H"))
    lines += ["TER", "END"]
    return _write_pdb(path, lines)


def _build_disulfide_pdb(path: Path) -> Path:
    """
    Two cysteines with SG-SG distance 2.0 A -> both should become CYX.
    """
    lines = [
        _make_atom_line(1, "N", "CYS", "A", 1, 0.000, 0.000, 0.000, "N"),
        _make_atom_line(2, "CA", "CYS", "A", 1, 1.460, 0.000, 0.000, "C"),
        _make_atom_line(3, "C", "CYS", "A", 1, 2.900, 0.000, 0.000, "C"),
        _make_atom_line(4, "O", "CYS", "A", 1, 3.900, 0.000, 0.000, "O"),
        _make_atom_line(5, "SG", "CYS", "A", 1, 0.000, 2.000, 0.000, "S"),
        _make_atom_line(6, "N", "CYM", "A", 2, 4.500, 0.000, 0.000, "N"),
        _make_atom_line(7, "CA", "CYM", "A", 2, 5.960, 0.000, 0.000, "C"),
        _make_atom_line(8, "C", "CYM", "A", 2, 7.400, 0.000, 0.000, "C"),
        _make_atom_line(9, "O", "CYM", "A", 2, 8.400, 0.000, 0.000, "O"),
        _make_atom_line(10, "SG", "CYM", "A", 2, 2.000, 2.000, 0.000, "S"),
        "TER",
        "END",
    ]
    return _write_pdb(path, lines)


def _build_water_hetero_altloc_pdb(path: Path) -> Path:
    lines = [
        _make_atom_line(1, "N", "HIS", "A", 1, 0.000, 0.000, 0.000, "N"),
        _make_atom_line(2, "CA", "HIS", "A", 1, 1.460, 0.000, 0.000, "C"),
        _make_atom_line(3, "C", "HIS", "A", 1, 2.900, 0.000, 0.000, "C"),
        _make_atom_line(4, "O", "HIS", "A", 1, 3.900, 0.000, 0.000, "O"),
        _make_atom_line(5, "ND1", "HIS", "A", 1, 1.200, 1.100, 0.000, "N"),
        _make_atom_line(6, "NE2", "HIS", "A", 1, 2.200, 1.100, 0.000, "N"),
        _make_atom_line(7, "HD1", "HIS", "A", 1, 1.200, 2.000, 0.000, "H", altloc="B"),
        _make_atom_line(8, "HE2", "HIS", "A", 1, 2.200, 2.000, 0.000, "H", altloc="A"),
        _make_atom_line(
            9,
            "O",
            "HOH",
            "A",
            2,
            10.000,
            0.000,
            0.000,
            "O",
            record_name="HETATM",
        ),
        _make_atom_line(
            10,
            "C1",
            "LIG",
            "A",
            3,
            11.000,
            0.000,
            0.000,
            "C",
            record_name="HETATM",
        ),
        "TER",
        "END",
    ]
    return _write_pdb(path, lines)


def test_iter_atoms_altloc_ok_keeps_blank_and_a_only(tmp_path: Path) -> None:
    pdb_path = _build_water_hetero_altloc_pdb(tmp_path / "altloc_test.pdb")
    chain = _parse_first_chain(pdb_path)
    his_residue = _residue_by_resseq(chain, 1)

    kept_atom_names = {
        atom.get_name().strip() for atom in _iter_atoms_altloc_ok(his_residue)
    }

    assert "HE2" in kept_atom_names
    assert "HD1" not in kept_atom_names


def test_atom_names_uses_altloc_filter(tmp_path: Path) -> None:
    pdb_path = _build_water_hetero_altloc_pdb(tmp_path / "atom_names_test.pdb")
    chain = _parse_first_chain(pdb_path)
    his_residue = _residue_by_resseq(chain, 1)

    atom_name_set = _atom_names(his_residue)

    assert "HE2" in atom_name_set
    assert "HD1" not in atom_name_set


def test_distance_returns_expected_value() -> None:
    assert _distance((0.0, 0.0, 0.0), (0.0, 3.0, 4.0)) == pytest.approx(5.0)


def test_sanitize_variant_label_normalizes_symbols() -> None:
    assert sanitize_variant_label("modeller complete") == "modeller_complete"
    assert sanitize_variant_label("alpha/fold") == "alpha_fold"
    assert sanitize_variant_label("gaps") == "gaps"


def test_sanitize_variant_label_raises_for_empty_result() -> None:
    with pytest.raises(ValueError):
        sanitize_variant_label("///")


def test_build_standard_amber_output_path_returns_expected_path(
    tmp_path: Path,
) -> None:
    output_path = build_standard_amber_output_path(
        pdb_id="1ABC",
        protein_dir=tmp_path / "1ABC",
    )

    assert output_path == (
        tmp_path / "1ABC" / "components" / "1ABC_protein_as_Amber.pdb"
    )


def test_build_variant_amber_output_path_returns_expected_path(
    tmp_path: Path,
) -> None:
    output_path = build_variant_amber_output_path(
        pdb_id="1ABC",
        protein_dir=tmp_path / "1ABC",
        variant_label="modeller complete",
    )

    assert output_path == (
        tmp_path / "1ABC" / "components" / "1ABC_protein_as_Amber_modeller_complete.pdb"
    )


def test_rename_structure_by_hydrogens_renames_his_to_hid(tmp_path: Path) -> None:
    pdb_path = _build_histidine_pdb(tmp_path / "his_hid.pdb", his_variant="hid")
    structure = _parse_structure(pdb_path)

    stats = rename_structure_by_hydrogens(structure)

    residue = next(structure.get_residues())
    assert residue.get_resname().strip() == "HID"
    assert stats["HIS->HID"] == 1


def test_rename_structure_by_hydrogens_renames_his_to_hie(tmp_path: Path) -> None:
    pdb_path = _build_histidine_pdb(tmp_path / "his_hie.pdb", his_variant="hie")
    structure = _parse_structure(pdb_path)

    stats = rename_structure_by_hydrogens(structure)

    residue = next(structure.get_residues())
    assert residue.get_resname().strip() == "HIE"
    assert stats["HIS->HIE"] == 1


def test_rename_structure_by_hydrogens_renames_his_to_hip(tmp_path: Path) -> None:
    pdb_path = _build_histidine_pdb(tmp_path / "his_hip.pdb", his_variant="hip")
    structure = _parse_structure(pdb_path)

    stats = rename_structure_by_hydrogens(structure)

    residue = next(structure.get_residues())
    assert residue.get_resname().strip() == "HIP"
    assert stats["HIS->HIP"] == 1


def test_rename_structure_by_hydrogens_counts_his_without_protons_in_strict_mode(
    tmp_path: Path,
) -> None:
    pdb_path = _build_histidine_pdb(tmp_path / "his_none.pdb", his_variant="none")
    structure = _parse_structure(pdb_path)

    stats = rename_structure_by_hydrogens(structure, strict_his=True)

    residue = next(structure.get_residues())
    assert residue.get_resname().strip() == "HIS"
    assert stats["HIS_no_protons"] == 1


def test_rename_structure_by_hydrogens_renames_asp_and_glu_protonated_forms(
    tmp_path: Path,
) -> None:
    pdb_path = _build_asp_glu_pdb(tmp_path / "asp_glu.pdb")
    structure = _parse_structure(pdb_path)

    stats = rename_structure_by_hydrogens(structure)

    chain = next(structure.get_chains())
    residue_names = _residue_names(chain)

    assert residue_names == ["ASH", "GLH"]
    assert stats["ASP->ASH"] == 1
    assert stats["GLU->GLH"] == 1


def test_rename_structure_by_hydrogens_renames_cys_to_cym_when_deprotonated(
    tmp_path: Path,
) -> None:
    pdb_path = _build_cysteine_pdb(tmp_path / "cys_no_hg.pdb", with_hg=False)
    structure = _parse_structure(pdb_path)

    stats = rename_structure_by_hydrogens(structure)

    residue = next(structure.get_residues())
    assert residue.get_resname().strip() == "CYM"
    assert stats["CYS->CYM"] == 1


def test_rename_structure_by_hydrogens_keeps_cys_when_hg_present(
    tmp_path: Path,
) -> None:
    pdb_path = _build_cysteine_pdb(tmp_path / "cys_hg.pdb", with_hg=True)
    structure = _parse_structure(pdb_path)

    stats = rename_structure_by_hydrogens(structure)

    residue = next(structure.get_residues())
    assert residue.get_resname().strip() == "CYS"
    assert stats["CYS kept"] == 1


def test_rename_structure_by_hydrogens_renames_disulfide_pair_to_cyx(
    tmp_path: Path,
) -> None:
    pdb_path = _build_disulfide_pdb(tmp_path / "disulfide.pdb")
    structure = _parse_structure(pdb_path)

    stats = rename_structure_by_hydrogens(
        structure,
        disulf_min=1.8,
        disulf_max=2.2,
    )

    chain = next(structure.get_chains())
    residue_names = _residue_names(chain)

    assert residue_names == ["CYX", "CYX"]
    assert stats["CYS->CYX"] == 1
    assert stats["CYM->CYX"] == 1


def test_rename_structure_by_hydrogens_ignores_water_and_hetero_residues(
    tmp_path: Path,
) -> None:
    pdb_path = _build_water_hetero_altloc_pdb(tmp_path / "water_hetero.pdb")
    structure = _parse_structure(pdb_path)

    stats = rename_structure_by_hydrogens(structure)

    model = next(structure.get_models())
    chain = next(model.get_chains())
    residue_names = [res.get_resname().strip() for res in chain.get_residues()]

    assert residue_names[0] == "HIE"
    assert residue_names[1] == "HOH"
    assert residue_names[2] == "LIG"
    assert stats["HIS->HIE"] == 1


def test_write_amber_renamed_pdb_writes_output_stats_json_and_log(
    tmp_path: Path,
) -> None:
    input_pdb_path = _build_histidine_pdb(tmp_path / "input_his.pdb", his_variant="hip")
    output_pdb_path = tmp_path / "output_his_as_Amber.pdb"
    log_path = tmp_path / "amber_renamed.log"

    stats = write_amber_renamed_pdb(
        input_pdb_path=input_pdb_path,
        output_pdb_path=output_pdb_path,
        log_path=log_path,
    )

    assert output_pdb_path.exists()
    assert log_path.exists()
    assert stats["HIS->HIP"] == 1

    stats_path = tmp_path / "output_his_as_Amber_stats.json"
    assert stats_path.exists()

    stats_data = json.loads(stats_path.read_text(encoding="utf-8"))
    assert stats_data["HIS->HIP"] == 1

    log_text = log_path.read_text(encoding="utf-8")
    assert "input=" in log_text
    assert "output=" in log_text
    assert "HIS->HIP=1" in log_text


def test_write_amber_renamed_pdb_raises_for_missing_input(tmp_path: Path) -> None:
    with pytest.raises(FileNotFoundError):
        write_amber_renamed_pdb(
            input_pdb_path=tmp_path / "missing.pdb",
            output_pdb_path=tmp_path / "out.pdb",
        )


def test_amber_rename_selected_structure_runs_on_explicit_input_and_returns_summary(
    tmp_path: Path,
) -> None:
    input_pdb_path = _build_asp_glu_pdb(tmp_path / "input_selected.pdb")
    output_pdb_path = tmp_path / "output_selected_as_Amber.pdb"
    log_path = tmp_path / "amber_renamed.log"

    result = amber_rename_selected_structure(
        input_pdb_path=input_pdb_path,
        output_pdb_path=output_pdb_path,
        strict_his=False,
        disulf_min=1.8,
        disulf_max=2.2,
        variant_label="modeller_complete",
        log_path=log_path,
    )

    assert output_pdb_path.exists()
    assert Path(result["amber_input_path"]) == input_pdb_path.resolve()
    assert Path(result["amber_output_path"]) == output_pdb_path.resolve()
    assert Path(result["amber_log_path"]) == log_path.resolve()
    assert result["amber_renaming_success"] is True
    assert result["amber_variant"] == "modeller_complete"
    assert result["asp_to_ash"] == 1
    assert result["glu_to_glh"] == 1
    assert result["his_to_hid"] == 0
    assert result["his_to_hie"] == 0
    assert result["his_to_hip"] == 0
    assert result["cys_to_cym"] == 0
    assert result["cym_to_cys"] == 0
    assert result["cys_to_cyx"] == 0
    assert result["cym_to_cyx"] == 0

    stats_json_path = output_pdb_path.with_name(f"{output_pdb_path.stem}_stats.json")
    assert Path(result["amber_stats_json_path"]) == stats_json_path
    assert stats_json_path.exists()


def test_amber_rename_variant_structure_writes_variant_specific_output(
    tmp_path: Path,
) -> None:
    pdb_id = "1ABC"
    protein_dir = tmp_path / pdb_id
    input_pdb_path = _build_asp_glu_pdb(tmp_path / "variant_input.pdb")

    result = amber_rename_variant_structure(
        pdb_id=pdb_id,
        protein_dir=protein_dir,
        variant_label="alphafold_complete",
        input_pdb_path=input_pdb_path,
    )

    output_pdb_path = (
        protein_dir / "components" / f"{pdb_id}_protein_as_Amber_alphafold_complete.pdb"
    )
    stats_json_path = (
        protein_dir
        / "components"
        / f"{pdb_id}_protein_as_Amber_alphafold_complete_stats.json"
    )
    log_path = protein_dir / "components" / "amber_renamed.log"

    assert output_pdb_path.exists()
    assert stats_json_path.exists()
    assert log_path.exists()
    assert result["amber_renaming_success"] is True
    assert result["amber_variant"] == "alphafold_complete"
    assert Path(result["amber_output_path"]) == output_pdb_path.resolve()
    assert Path(result["amber_stats_json_path"]) == stats_json_path.resolve()
    assert Path(result["amber_log_path"]) == log_path.resolve()
    assert result["asp_to_ash"] == 1
    assert result["glu_to_glh"] == 1


def test_amber_rename_protein_structure_uses_standard_paths_and_returns_summary(
    tmp_path: Path,
) -> None:
    pdb_id = "1ABC"
    protein_dir = tmp_path / pdb_id
    components_dir = protein_dir / "components"
    input_pdb_path = components_dir / f"{pdb_id}_proteinH.pdb"

    _build_asp_glu_pdb(input_pdb_path)

    result = amber_rename_protein_structure(
        pdb_id=pdb_id,
        protein_dir=protein_dir,
    )

    output_pdb_path = components_dir / f"{pdb_id}_protein_as_Amber.pdb"
    stats_json_path = components_dir / f"{pdb_id}_protein_as_Amber_stats.json"
    log_path = components_dir / "amber_renamed.log"

    assert result["amber_renaming_success"] is True
    assert Path(result["amber_input_path"]) == input_pdb_path.resolve()
    assert Path(result["amber_output_path"]) == output_pdb_path.resolve()
    assert Path(result["amber_stats_json_path"]) == stats_json_path.resolve()
    assert Path(result["amber_log_path"]) == log_path.resolve()

    assert output_pdb_path.exists()
    assert stats_json_path.exists()
    assert log_path.exists()

    assert result["asp_to_ash"] == 1
    assert result["glu_to_glh"] == 1
    assert result["his_to_hid"] == 0
    assert result["his_to_hie"] == 0
    assert result["his_to_hip"] == 0
    assert result["cys_to_cym"] == 0
    assert result["cym_to_cys"] == 0
    assert result["cys_to_cyx"] == 0
    assert result["cym_to_cyx"] == 0
