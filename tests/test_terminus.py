from __future__ import annotations

from pathlib import Path

import pytest
from Bio.PDB import PDBParser

from stack_protein_preparation.terminus import (
    build_default_terminus_output_path,
    build_variant_terminus_output_path,
    convert_protein_termini_to_amber,
    convert_variant_protein_termini_to_amber,
    sanitize_variant_label,
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
) -> str:
    return (
        f"ATOM  {serial:5d} {atom_name:>4} {residue_name:>3} {chain_id:1}"
        f"{residue_number:4d}    "
        f"{x:8.3f}{y:8.3f}{z:8.3f}"
        f"{1.00:6.2f}{20.00:6.2f}          "
        f"{element:>2}"
    )


def _write_test_pdb(path: Path, lines: list[str]) -> Path:
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def _build_three_residue_chain_pdb(path: Path) -> Path:
    lines = [
        _make_atom_line(1, "N", "ALA", "A", 1, 0.000, 0.000, 0.000, "N"),
        _make_atom_line(2, "H1", "ALA", "A", 1, -0.700, 0.300, 0.000, "H"),
        _make_atom_line(3, "H2", "ALA", "A", 1, -0.700, -0.300, 0.000, "H"),
        _make_atom_line(4, "H3", "ALA", "A", 1, 0.000, 0.000, 0.900, "H"),
        _make_atom_line(5, "CA", "ALA", "A", 1, 1.460, 0.000, 0.000, "C"),
        _make_atom_line(6, "C", "ALA", "A", 1, 2.900, 0.000, 0.000, "C"),
        _make_atom_line(7, "O", "ALA", "A", 1, 3.900, 0.000, 0.000, "O"),
        _make_atom_line(8, "N", "GLY", "A", 2, 4.200, 0.000, 0.000, "N"),
        _make_atom_line(9, "CA", "GLY", "A", 2, 5.660, 0.000, 0.000, "C"),
        _make_atom_line(10, "C", "GLY", "A", 2, 7.100, 0.000, 0.000, "C"),
        _make_atom_line(11, "O", "GLY", "A", 2, 8.100, 0.000, 0.000, "O"),
        _make_atom_line(12, "N", "SER", "A", 3, 8.400, 0.000, 0.000, "N"),
        _make_atom_line(13, "CA", "SER", "A", 3, 9.860, 0.000, 0.000, "C"),
        _make_atom_line(14, "C", "SER", "A", 3, 11.300, 0.000, 0.000, "C"),
        _make_atom_line(15, "O", "SER", "A", 3, 12.300, 0.000, 0.000, "O"),
        _make_atom_line(16, "OXT", "SER", "A", 3, 11.300, 1.250, 0.000, "O"),
        "TER",
        "END",
    ]
    return _write_test_pdb(path, lines)


def _build_single_residue_chain_pdb(path: Path) -> Path:
    lines = [
        _make_atom_line(1, "N", "GLY", "A", 7, 0.000, 0.000, 0.000, "N"),
        _make_atom_line(2, "H1", "GLY", "A", 7, -0.700, 0.300, 0.000, "H"),
        _make_atom_line(3, "CA", "GLY", "A", 7, 1.460, 0.000, 0.000, "C"),
        _make_atom_line(4, "C", "GLY", "A", 7, 2.900, 0.000, 0.000, "C"),
        _make_atom_line(5, "O", "GLY", "A", 7, 3.900, 0.000, 0.000, "O"),
        "TER",
        "END",
    ]
    return _write_test_pdb(path, lines)


def _build_chain_with_internal_caps_pdb(path: Path) -> Path:
    lines = [
        _make_atom_line(1, "N", "ALA", "A", 1, 0.000, 0.000, 0.000, "N"),
        _make_atom_line(2, "CA", "ALA", "A", 1, 1.460, 0.000, 0.000, "C"),
        _make_atom_line(3, "C", "ALA", "A", 1, 2.900, 0.000, 0.000, "C"),
        _make_atom_line(4, "O", "ALA", "A", 1, 3.900, 0.000, 0.000, "O"),
        _make_atom_line(5, "CH3", "NME", "A", 2, 4.300, 0.000, 0.000, "C"),
        _make_atom_line(6, "N", "NME", "A", 2, 5.400, 0.000, 0.000, "N"),
        _make_atom_line(7, "CH3", "ACE", "A", 3, 6.700, 0.000, 0.000, "C"),
        _make_atom_line(8, "C", "ACE", "A", 3, 7.900, 0.000, 0.000, "C"),
        _make_atom_line(9, "O", "ACE", "A", 3, 8.900, 0.000, 0.000, "O"),
        _make_atom_line(10, "N", "SER", "A", 4, 9.300, 0.000, 0.000, "N"),
        _make_atom_line(11, "CA", "SER", "A", 4, 10.760, 0.000, 0.000, "C"),
        _make_atom_line(12, "C", "SER", "A", 4, 12.200, 0.000, 0.000, "C"),
        _make_atom_line(13, "O", "SER", "A", 4, 13.200, 0.000, 0.000, "O"),
        _make_atom_line(14, "OXT", "SER", "A", 4, 12.200, 1.250, 0.000, "O"),
        "TER",
        "END",
    ]
    return _write_test_pdb(path, lines)


def _parse_first_chain(pdb_path: Path):
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("test", str(pdb_path))
    model = next(structure.get_models())
    return next(model.get_chains())


def _residue_names(chain) -> list[str]:
    return [residue.get_resname().strip() for residue in chain.get_residues()]


def _residues(chain):
    return list(chain.get_residues())


def _atom_names(residue) -> set[str]:
    return {atom.get_name().strip() for atom in residue.get_atoms()}


def test_sanitize_variant_label_normalizes_symbols() -> None:
    assert sanitize_variant_label("modeller complete") == "modeller_complete"
    assert sanitize_variant_label("alpha/fold") == "alpha_fold"


def test_sanitize_variant_label_raises_for_empty_result() -> None:
    with pytest.raises(ValueError):
        sanitize_variant_label("///")


def test_build_default_terminus_output_path_returns_expected_path(
    tmp_path: Path,
) -> None:
    output_path = build_default_terminus_output_path(
        pdb_directory=tmp_path / "1ABC",
        pdb_id="1ABC",
    )
    assert output_path == (
        tmp_path / "1ABC" / "components" / "1ABC_protein_amber_termini.pdb"
    )


def test_build_variant_terminus_output_path_returns_expected_path(
    tmp_path: Path,
) -> None:
    output_path = build_variant_terminus_output_path(
        pdb_directory=tmp_path / "1ABC",
        pdb_id="1ABC",
        variant_label="modeller complete",
    )
    assert output_path == (
        tmp_path
        / "1ABC"
        / "components"
        / "1ABC_protein_amber_termini_modeller_complete.pdb"
    )


def test_convert_protein_termini_to_amber_converts_true_chain_ends(
    tmp_path: Path,
) -> None:
    input_pdb_path = _build_three_residue_chain_pdb(
        tmp_path / "input_three_residues.pdb"
    )
    output_pdb_path = tmp_path / "output_three_residues.pdb"

    result_list = convert_protein_termini_to_amber(
        input_pdb_path=input_pdb_path,
        output_pdb_path=output_pdb_path,
    )

    assert output_pdb_path.exists()
    assert len(result_list) == 1

    result = result_list[0]
    assert result.chain_id == "A"
    assert result.first_resseq == 1
    assert result.last_resseq == 3
    assert result.first_old_resname == "ALA"
    assert result.first_new_resname == "NALA"
    assert result.last_old_resname == "SER"
    assert result.last_new_resname == "CSER"
    assert result.single_residue_chain is False

    output_chain = _parse_first_chain(output_pdb_path)
    residue_list = _residues(output_chain)

    assert _residue_names(output_chain) == ["ALA", "GLY", "SER"]

    first_residue = residue_list[0]
    last_residue = residue_list[-1]

    assert {"H1", "H2", "H3"}.issubset(_atom_names(first_residue))
    assert "O" in _atom_names(last_residue)
    assert "OXT" in _atom_names(last_residue)


def test_convert_protein_termini_to_amber_handles_single_residue_chain(
    tmp_path: Path,
) -> None:
    input_pdb_path = _build_single_residue_chain_pdb(
        tmp_path / "input_single_residue.pdb"
    )
    output_pdb_path = tmp_path / "output_single_residue.pdb"

    result_list = convert_protein_termini_to_amber(
        input_pdb_path=input_pdb_path,
        output_pdb_path=output_pdb_path,
    )

    assert output_pdb_path.exists()
    assert len(result_list) == 1

    result = result_list[0]
    assert result.chain_id == "A"
    assert result.first_resseq == 7
    assert result.last_resseq == 7
    assert result.first_old_resname == "GLY"
    assert result.last_old_resname == "GLY"
    assert result.first_new_resname == "GLY"
    assert result.last_new_resname == "GLY"
    assert result.single_residue_chain is True

    output_chain = _parse_first_chain(output_pdb_path)
    residue_list = _residues(output_chain)

    assert _residue_names(output_chain) == ["GLY"]

    residue = residue_list[0]
    atom_name_set = _atom_names(residue)

    assert "H1" in atom_name_set
    assert "O" in atom_name_set
    assert "OXT" in atom_name_set


def test_convert_protein_termini_to_amber_does_not_modify_internal_ace_nme(
    tmp_path: Path,
) -> None:
    input_pdb_path = _build_chain_with_internal_caps_pdb(
        tmp_path / "input_internal_caps.pdb"
    )
    output_pdb_path = tmp_path / "output_internal_caps.pdb"

    result_list = convert_protein_termini_to_amber(
        input_pdb_path=input_pdb_path,
        output_pdb_path=output_pdb_path,
    )

    assert len(result_list) == 1
    assert result_list[0].first_new_resname == "NALA"
    assert result_list[0].last_new_resname == "CSER"

    output_chain = _parse_first_chain(output_pdb_path)
    residue_name_list = _residue_names(output_chain)

    assert residue_name_list == ["ALA", "NME", "ACE", "SER"]


def test_convert_protein_termini_to_amber_output_text_keeps_chain_and_resseq_positions_parseable(
    tmp_path: Path,
) -> None:
    input_pdb_path = _build_three_residue_chain_pdb(tmp_path / "input_text_check.pdb")
    output_pdb_path = tmp_path / "output_text_check.pdb"

    convert_protein_termini_to_amber(
        input_pdb_path=input_pdb_path,
        output_pdb_path=output_pdb_path,
    )

    text_lines = output_pdb_path.read_text(encoding="utf-8").splitlines()
    atom_lines = [line for line in text_lines if line.startswith("ATOM")]

    assert atom_lines, "No ATOM lines written"

    for line in atom_lines:
        chain_id = line[21].strip()
        residue_number_text = line[22:26].strip()

        assert chain_id == "A"
        assert residue_number_text.isdigit()

    output_chain = _parse_first_chain(output_pdb_path)
    residue_list = _residues(output_chain)

    assert residue_list[0].id[1] == 1
    assert residue_list[1].id[1] == 2
    assert residue_list[2].id[1] == 3


def test_convert_protein_termini_to_amber_raises_for_missing_input(
    tmp_path: Path,
) -> None:
    missing_input_path = tmp_path / "does_not_exist.pdb"
    output_pdb_path = tmp_path / "unused_output.pdb"

    with pytest.raises(FileNotFoundError):
        convert_protein_termini_to_amber(
            input_pdb_path=missing_input_path,
            output_pdb_path=output_pdb_path,
        )


def test_convert_variant_protein_termini_to_amber_writes_variant_specific_output(
    tmp_path: Path,
) -> None:
    pdb_id = "1ABC"
    pdb_directory = tmp_path / pdb_id
    input_pdb_path = _build_three_residue_chain_pdb(tmp_path / "variant_input.pdb")

    result = convert_variant_protein_termini_to_amber(
        pdb_id=pdb_id,
        pdb_directory=pdb_directory,
        variant_label="alphafold_complete",
        input_pdb_path=input_pdb_path,
    )

    output_pdb_path = (
        pdb_directory
        / "components"
        / f"{pdb_id}_protein_amber_termini_alphafold_complete.pdb"
    )

    assert output_pdb_path.exists()
    assert result["amber_termini_success"] is True
    assert result["amber_termini_variant"] == "alphafold_complete"
    assert Path(result["amber_termini_input_path"]) == input_pdb_path.resolve()
    assert Path(result["amber_termini_output_path"]) == output_pdb_path.resolve()
    assert result["amber_termini_chain_count"] == 1
