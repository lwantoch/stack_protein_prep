from __future__ import annotations

from pathlib import Path

import pytest

from stack_protein_preparation.prepared_structure import (
    PreparedStructureSummary,
    _is_atom_or_hetatm_record,
    _read_atom_lines_from_pdb,
    _renumber_atom_serial,
    _write_merged_pdb_sections,
    build_prepared_structure,
    build_prepared_structure_for_pdb_directory,
    get_default_ligand_input_path,
    get_default_metals_input_path,
    get_default_prepared_protein_input_path,
    get_default_water_input_path,
    get_prepared_structure_output_path,
    sanitize_variant_label,
    summarize_prepared_structure,
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
    record_name: str = "ATOM",
) -> str:
    return (
        f"{record_name:<6}"
        f"{serial:5d} "
        f"{atom_name:>4}"
        f" "
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


def _build_protein_pdb(path: Path) -> Path:
    lines = [
        _make_atom_line(10, "N", "ALA", "A", 1, 0.0, 0.0, 0.0, "N"),
        _make_atom_line(11, "CA", "ALA", "A", 1, 1.4, 0.0, 0.0, "C"),
        _make_atom_line(12, "C", "ALA", "A", 1, 2.8, 0.0, 0.0, "C"),
        _make_atom_line(13, "O", "ALA", "A", 1, 3.8, 0.0, 0.0, "O"),
        "TER",
        "END",
    ]
    return _write_pdb(path, lines)


def _build_water_pdb(path: Path) -> Path:
    lines = [
        _make_atom_line(40, "O", "HOH", "W", 10, 4.0, 4.0, 4.0, "O", "HETATM"),
        "TER",
        "END",
    ]
    return _write_pdb(path, lines)


def _build_ligand_pdb(path: Path) -> Path:
    lines = [
        _make_atom_line(50, "C1", "LIG", "L", 20, 5.0, 5.0, 5.0, "C", "HETATM"),
        _make_atom_line(51, "O1", "LIG", "L", 20, 6.0, 5.0, 5.0, "O", "HETATM"),
        "TER",
        "END",
    ]
    return _write_pdb(path, lines)


def _build_metals_pdb(path: Path) -> Path:
    lines = [
        _make_atom_line(60, "ZN", "ZN", "M", 30, 7.0, 7.0, 7.0, "ZN", "HETATM"),
        "TER",
        "END",
    ]
    return _write_pdb(path, lines)


def test_sanitize_variant_label_normalizes_symbols() -> None:
    assert sanitize_variant_label("alpha fold") == "alpha_fold"
    assert sanitize_variant_label("modeller/complete") == "modeller_complete"


def test_sanitize_variant_label_raises_for_empty_result() -> None:
    with pytest.raises(ValueError):
        sanitize_variant_label("///")


def test_default_component_paths_are_correct(tmp_path: Path) -> None:
    pdb_directory = tmp_path / "1ABC"

    assert get_default_prepared_protein_input_path(pdb_directory, "1ABC") == (
        pdb_directory / "components" / "1ABC_protein_internal_capped.pdb"
    )
    assert get_default_water_input_path(pdb_directory, "1ABC") == (
        pdb_directory / "components" / "1ABC_water.pdb"
    )
    assert get_default_ligand_input_path(pdb_directory, "1ABC") == (
        pdb_directory / "components" / "1ABC_ligand.pdb"
    )
    assert get_default_metals_input_path(pdb_directory, "1ABC") == (
        pdb_directory / "components" / "1ABC_metals.pdb"
    )


def test_get_prepared_structure_output_path_without_gaps(tmp_path: Path) -> None:
    output_path = get_prepared_structure_output_path(
        pdb_directory=tmp_path / "1ABC",
        pdb_id="1ABC",
        had_gaps=False,
    )
    assert output_path == tmp_path / "1ABC" / "prepared" / "1ABC.pdb"


def test_get_prepared_structure_output_path_with_gaps_variant_gaps(
    tmp_path: Path,
) -> None:
    output_path = get_prepared_structure_output_path(
        pdb_directory=tmp_path / "1ABC",
        pdb_id="1ABC",
        had_gaps=True,
        structure_variant="gaps",
    )
    assert output_path == tmp_path / "1ABC" / "prepared" / "gaps" / "1ABC.pdb"


def test_get_prepared_structure_output_path_with_gaps_variant_complete(
    tmp_path: Path,
) -> None:
    output_path = get_prepared_structure_output_path(
        pdb_directory=tmp_path / "1ABC",
        pdb_id="1ABC",
        had_gaps=True,
        structure_variant="complete",
    )
    assert output_path == tmp_path / "1ABC" / "prepared" / "complete" / "1ABC.pdb"


def test_get_prepared_structure_output_path_raises_for_missing_variant(
    tmp_path: Path,
) -> None:
    with pytest.raises(ValueError):
        get_prepared_structure_output_path(
            pdb_directory=tmp_path / "1ABC",
            pdb_id="1ABC",
            had_gaps=True,
            structure_variant=None,
        )


def test_is_atom_or_hetatm_record_detects_valid_lines() -> None:
    assert _is_atom_or_hetatm_record("ATOM      1  N   ALA A   1") is True
    assert _is_atom_or_hetatm_record("HETATM    2  ZN   ZN  M   1") is True
    assert _is_atom_or_hetatm_record("TER") is False
    assert _is_atom_or_hetatm_record("END") is False


def test_read_atom_lines_from_pdb_reads_only_atom_and_hetatm(tmp_path: Path) -> None:
    pdb_path = tmp_path / "mixed.pdb"
    _write_pdb(
        pdb_path,
        [
            _make_atom_line(1, "N", "ALA", "A", 1, 0, 0, 0, "N"),
            "TER",
            _make_atom_line(2, "ZN", "ZN", "M", 2, 1, 1, 1, "ZN", "HETATM"),
            "END",
        ],
    )

    lines = _read_atom_lines_from_pdb(pdb_path)

    assert len(lines) == 2
    assert lines[0].startswith("ATOM")
    assert lines[1].startswith("HETATM")


def test_read_atom_lines_from_pdb_raises_for_missing_file(tmp_path: Path) -> None:
    with pytest.raises(FileNotFoundError):
        _read_atom_lines_from_pdb(tmp_path / "missing.pdb")


def test_renumber_atom_serial_replaces_serial_field() -> None:
    line = _make_atom_line(77, "N", "ALA", "A", 1, 0, 0, 0, "N")
    renumbered = _renumber_atom_serial(line, 5)

    assert renumbered.startswith("ATOM      5")
    assert renumbered[12:] == line[12:]


def test_write_merged_pdb_sections_writes_sections_in_expected_order(tmp_path: Path) -> None:
    output_pdb_path = tmp_path / "merged.pdb"

    protein_lines = [_make_atom_line(10, "N", "ALA", "A", 1, 0, 0, 0, "N")]
    water_lines = [_make_atom_line(20, "O", "HOH", "W", 2, 1, 1, 1, "O", "HETATM")]
    ligand_lines = [_make_atom_line(30, "C1", "LIG", "L", 3, 2, 2, 2, "C", "HETATM")]
    metals_lines = [_make_atom_line(40, "ZN", "ZN", "M", 4, 3, 3, 3, "ZN", "HETATM")]

    n_written = _write_merged_pdb_sections(
        output_pdb_path=output_pdb_path,
        ordered_sections=[protein_lines, water_lines, ligand_lines, metals_lines],
    )

    assert n_written == 4

    text_lines = output_pdb_path.read_text(encoding="utf-8").splitlines()

    atom_lines = [line for line in text_lines if line.startswith(("ATOM", "HETATM"))]
    assert len(atom_lines) == 4
    assert atom_lines[0].startswith("ATOM      1")
    assert atom_lines[1].startswith("HETATM    2")
    assert atom_lines[2].startswith("HETATM    3")
    assert atom_lines[3].startswith("HETATM    4")

    assert text_lines[-1] == "END"
    assert text_lines.count("TER") == 4


def test_build_prepared_structure_with_only_protein(tmp_path: Path) -> None:
    protein_input_path = _build_protein_pdb(tmp_path / "protein.pdb")
    output_pdb_path = tmp_path / "prepared" / "1ABC.pdb"

    summary = build_prepared_structure(
        output_pdb_path=output_pdb_path,
        protein_input_path=protein_input_path,
    )

    assert summary == PreparedStructureSummary(
        pdb_id="1ABC",
        output_pdb_path=output_pdb_path,
        protein_input_path=protein_input_path,
        water_input_path=None,
        ligand_input_path=None,
        metals_input_path=None,
        had_gaps=False,
        structure_variant=None,
        water_included=False,
        ligand_included=False,
        metals_included=False,
        n_atom_records_written=4,
    )

    text_lines = output_pdb_path.read_text(encoding="utf-8").splitlines()
    atom_lines = [line for line in text_lines if line.startswith(("ATOM", "HETATM"))]
    assert len(atom_lines) == 4
    assert text_lines[-2] == "TER"
    assert text_lines[-1] == "END"


def test_build_prepared_structure_includes_optional_components(tmp_path: Path) -> None:
    protein_input_path = _build_protein_pdb(tmp_path / "protein.pdb")
    water_input_path = _build_water_pdb(tmp_path / "water.pdb")
    ligand_input_path = _build_ligand_pdb(tmp_path / "ligand.pdb")
    metals_input_path = _build_metals_pdb(tmp_path / "metals.pdb")
    output_pdb_path = tmp_path / "prepared" / "1ABC.pdb"

    summary = build_prepared_structure(
        output_pdb_path=output_pdb_path,
        protein_input_path=protein_input_path,
        water_input_path=water_input_path,
        ligand_input_path=ligand_input_path,
        metals_input_path=metals_input_path,
    )

    assert summary.water_included is True
    assert summary.ligand_included is True
    assert summary.metals_included is True
    assert summary.water_input_path == water_input_path
    assert summary.ligand_input_path == ligand_input_path
    assert summary.metals_input_path == metals_input_path
    assert summary.n_atom_records_written == 8

    text_lines = output_pdb_path.read_text(encoding="utf-8").splitlines()
    atom_lines = [line for line in text_lines if line.startswith(("ATOM", "HETATM"))]

    assert len(atom_lines) == 8
    assert atom_lines[0][17:20].strip() == "ALA"
    assert atom_lines[4][17:20].strip() == "HOH"
    assert atom_lines[5][17:20].strip() == "LIG"
    assert atom_lines[7][17:20].strip() == "ZN"


def test_build_prepared_structure_ignores_missing_optional_components(
    tmp_path: Path,
) -> None:
    protein_input_path = _build_protein_pdb(tmp_path / "protein.pdb")
    output_pdb_path = tmp_path / "prepared" / "1ABC.pdb"

    summary = build_prepared_structure(
        output_pdb_path=output_pdb_path,
        protein_input_path=protein_input_path,
        water_input_path=tmp_path / "missing_water.pdb",
        ligand_input_path=tmp_path / "missing_ligand.pdb",
        metals_input_path=tmp_path / "missing_metals.pdb",
    )

    assert summary.water_included is False
    assert summary.ligand_included is False
    assert summary.metals_included is False
    assert summary.water_input_path is None
    assert summary.ligand_input_path is None
    assert summary.metals_input_path is None
    assert summary.n_atom_records_written == 4


def test_build_prepared_structure_raises_for_missing_protein_input(
    tmp_path: Path,
) -> None:
    with pytest.raises(FileNotFoundError):
        build_prepared_structure(
            output_pdb_path=tmp_path / "out.pdb",
            protein_input_path=tmp_path / "missing_protein.pdb",
        )


def test_build_prepared_structure_raises_for_empty_protein_records(
    tmp_path: Path,
) -> None:
    protein_input_path = _write_pdb(tmp_path / "empty_protein.pdb", ["TER", "END"])

    with pytest.raises(ValueError):
        build_prepared_structure(
            output_pdb_path=tmp_path / "out.pdb",
            protein_input_path=protein_input_path,
        )


def test_build_prepared_structure_for_pdb_directory_without_gaps_uses_default_paths(
    tmp_path: Path,
) -> None:
    pdb_directory = tmp_path / "1ABC"
    components_dir = pdb_directory / "components"

    _build_protein_pdb(components_dir / "1ABC_protein_internal_capped.pdb")
    _build_water_pdb(components_dir / "1ABC_water.pdb")
    _build_ligand_pdb(components_dir / "1ABC_ligand.pdb")
    _build_metals_pdb(components_dir / "1ABC_metals.pdb")

    summary = build_prepared_structure_for_pdb_directory(
        pdb_directory=pdb_directory,
        pdb_id="1ABC",
        had_gaps=False,
    )

    assert summary.output_pdb_path == pdb_directory / "prepared" / "1ABC.pdb"
    assert summary.had_gaps is False
    assert summary.structure_variant is None
    assert summary.water_included is True
    assert summary.ligand_included is True
    assert summary.metals_included is True


def test_build_prepared_structure_for_pdb_directory_with_gaps_complete_variant(
    tmp_path: Path,
) -> None:
    pdb_directory = tmp_path / "1ABC"
    protein_input_path = _build_protein_pdb(tmp_path / "complete_variant_protein.pdb")

    summary = build_prepared_structure_for_pdb_directory(
        pdb_directory=pdb_directory,
        pdb_id="1ABC",
        had_gaps=True,
        structure_variant="complete",
        protein_input_path=protein_input_path,
    )

    assert summary.output_pdb_path == (
        pdb_directory / "prepared" / "complete" / "1ABC.pdb"
    )
    assert summary.had_gaps is True
    assert summary.structure_variant == "complete"


def test_build_prepared_structure_for_pdb_directory_with_gaps_gaps_variant(
    tmp_path: Path,
) -> None:
    pdb_directory = tmp_path / "1ABC"
    protein_input_path = _build_protein_pdb(tmp_path / "gaps_variant_protein.pdb")

    summary = build_prepared_structure_for_pdb_directory(
        pdb_directory=pdb_directory,
        pdb_id="1ABC",
        had_gaps=True,
        structure_variant="gaps",
        protein_input_path=protein_input_path,
    )

    assert summary.output_pdb_path == pdb_directory / "prepared" / "gaps" / "1ABC.pdb"
    assert summary.had_gaps is True
    assert summary.structure_variant == "gaps"


def test_build_prepared_structure_for_pdb_directory_raises_for_invalid_gap_variant(
    tmp_path: Path,
) -> None:
    with pytest.raises(ValueError):
        build_prepared_structure_for_pdb_directory(
            pdb_directory=tmp_path / "1ABC",
            pdb_id="1ABC",
            had_gaps=True,
            structure_variant=None,
            protein_input_path=_build_protein_pdb(tmp_path / "protein.pdb"),
        )


def test_summarize_prepared_structure_formats_expected_text(tmp_path: Path) -> None:
    summary = PreparedStructureSummary(
        pdb_id="1ABC",
        output_pdb_path=tmp_path / "prepared" / "1ABC.pdb",
        protein_input_path=tmp_path / "components" / "1ABC_protein_internal_capped.pdb",
        water_input_path=tmp_path / "components" / "1ABC_water.pdb",
        ligand_input_path=None,
        metals_input_path=tmp_path / "components" / "1ABC_metals.pdb",
        had_gaps=True,
        structure_variant="gaps",
        water_included=True,
        ligand_included=False,
        metals_included=True,
        n_atom_records_written=6,
    )

    text = summarize_prepared_structure(summary)

    assert "pdb_id=1ABC" in text
    assert "structure_variant=gaps" in text
    assert "water_included=True" in text
    assert "ligand_included=False" in text
    assert "metals_included=True" in text
    assert "n_atom_records_written=6" in text
