from __future__ import annotations

from pathlib import Path

import pytest
from Bio.PDB import PDBParser

from stack_protein_preparation.cap import (
    ChainCappingSummary,
    Fragment,
    InternalGapBoundary,
    _has_peptide_bond,
    _is_polymer_residue,
    cap_internal_gap_boundaries,
    cap_internal_gaps_for_pdb_directory,
    cap_internal_gaps_for_variant,
    cap_internal_gaps_in_structure,
    find_internal_fragments,
    find_internal_gap_boundaries,
    get_default_internal_capped_output_path,
    get_variant_internal_capped_output_path,
    sanitize_variant_label,
    summarize_capping_results,
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


def _build_contiguous_two_residue_chain(path: Path) -> Path:
    lines = [
        _make_atom_line(1, "N", "ALA", "A", 1, 0.000, 0.000, 0.000, "N"),
        _make_atom_line(2, "CA", "ALA", "A", 1, 1.460, 0.000, 0.000, "C"),
        _make_atom_line(3, "C", "ALA", "A", 1, 2.900, 0.000, 0.000, "C"),
        _make_atom_line(4, "O", "ALA", "A", 1, 3.900, 0.000, 0.000, "O"),
        _make_atom_line(5, "N", "GLY", "A", 2, 4.200, 0.000, 0.000, "N"),
        _make_atom_line(6, "H1", "GLY", "A", 2, 4.700, 0.700, 0.000, "H"),
        _make_atom_line(7, "CA", "GLY", "A", 2, 5.660, 0.000, 0.000, "C"),
        _make_atom_line(8, "C", "GLY", "A", 2, 7.100, 0.000, 0.000, "C"),
        _make_atom_line(9, "O", "GLY", "A", 2, 8.100, 0.000, 0.000, "O"),
        _make_atom_line(10, "OXT", "GLY", "A", 2, 7.100, 1.200, 0.000, "O"),
        "TER",
        "END",
    ]
    return _write_pdb(path, lines)


def _build_internal_gap_chain(path: Path) -> Path:
    """
    Chain A with two fragments:
    residue 1 -> residue 5 not peptide-bonded
    """
    lines = [
        _make_atom_line(1, "N", "ALA", "A", 1, 0.000, 0.000, 0.000, "N"),
        _make_atom_line(2, "H1", "ALA", "A", 1, -0.700, 0.300, 0.000, "H"),
        _make_atom_line(3, "CA", "ALA", "A", 1, 1.460, 0.000, 0.000, "C"),
        _make_atom_line(4, "C", "ALA", "A", 1, 2.900, 0.000, 0.000, "C"),
        _make_atom_line(5, "O", "ALA", "A", 1, 3.900, 0.000, 0.000, "O"),
        _make_atom_line(6, "OXT", "ALA", "A", 1, 2.900, 1.200, 0.000, "O"),
        _make_atom_line(7, "N", "SER", "A", 5, 15.000, 0.000, 0.000, "N"),
        _make_atom_line(8, "H1", "SER", "A", 5, 15.500, 0.700, 0.000, "H"),
        _make_atom_line(9, "CA", "SER", "A", 5, 16.460, 0.000, 0.000, "C"),
        _make_atom_line(10, "C", "SER", "A", 5, 17.900, 0.000, 0.000, "C"),
        _make_atom_line(11, "O", "SER", "A", 5, 18.900, 0.000, 0.000, "O"),
        _make_atom_line(12, "OXT", "SER", "A", 5, 17.900, 1.200, 0.000, "O"),
        "TER",
        "END",
    ]
    return _write_pdb(path, lines)


def _build_non_polymer_only(path: Path) -> Path:
    lines = [
        _make_atom_line(1, "O", "HOH", "A", 1, 0.0, 0.0, 0.0, "O", "HETATM"),
        _make_atom_line(2, "ZN", "ZN", "A", 2, 2.0, 0.0, 0.0, "ZN", "HETATM"),
        "TER",
        "END",
    ]
    return _write_pdb(path, lines)


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


def test_get_default_internal_capped_output_path_returns_expected_path(
    tmp_path: Path,
) -> None:
    output_path = get_default_internal_capped_output_path(
        pdb_directory=tmp_path / "1ABC",
        pdb_id="1ABC",
    )
    assert output_path == (
        tmp_path / "1ABC" / "components" / "1ABC_protein_internal_capped.pdb"
    )


def test_get_variant_internal_capped_output_path_returns_expected_path(
    tmp_path: Path,
) -> None:
    output_path = get_variant_internal_capped_output_path(
        pdb_directory=tmp_path / "1ABC",
        pdb_id="1ABC",
        variant_label="alphafold complete",
    )
    assert output_path == (
        tmp_path
        / "1ABC"
        / "components"
        / "1ABC_protein_internal_capped_alphafold_complete.pdb"
    )


def test_is_polymer_residue_distinguishes_standard_and_hetero(tmp_path: Path) -> None:
    pdb_path = _build_internal_gap_chain(tmp_path / "polymer_check.pdb")
    chain = _parse_first_chain(pdb_path)
    residues = _residues(chain)

    assert _is_polymer_residue(residues[0]) is True
    assert _is_polymer_residue(residues[1]) is True

    non_polymer_path = _build_non_polymer_only(tmp_path / "non_polymer_check.pdb")
    chain2 = _parse_first_chain(non_polymer_path)
    residues2 = _residues(chain2)

    assert _is_polymer_residue(residues2[0]) is False


def test_has_peptide_bond_detects_connected_and_disconnected_pairs(
    tmp_path: Path,
) -> None:
    contiguous_path = _build_contiguous_two_residue_chain(tmp_path / "contiguous.pdb")
    chain = _parse_first_chain(contiguous_path)
    residues = _residues(chain)

    assert _has_peptide_bond(residues[0], residues[1], peptide_bond_max_distance=1.8)

    gap_path = _build_internal_gap_chain(tmp_path / "gap_chain.pdb")
    chain_gap = _parse_first_chain(gap_path)
    residues_gap = _residues(chain_gap)

    assert not _has_peptide_bond(
        residues_gap[0],
        residues_gap[1],
        peptide_bond_max_distance=1.8,
    )


def test_find_internal_fragments_returns_single_fragment_for_contiguous_chain(
    tmp_path: Path,
) -> None:
    pdb_path = _build_contiguous_two_residue_chain(tmp_path / "single_fragment.pdb")
    chain = _parse_first_chain(pdb_path)

    fragments = find_internal_fragments(chain, peptide_bond_max_distance=1.8)

    assert len(fragments) == 1
    assert fragments[0] == Fragment(
        chain_id="A",
        start_index=0,
        end_index=1,
        start_resseq=1,
        end_resseq=2,
        start_resname="ALA",
        end_resname="GLY",
    )


def test_find_internal_fragments_returns_two_fragments_for_internal_gap(
    tmp_path: Path,
) -> None:
    pdb_path = _build_internal_gap_chain(tmp_path / "two_fragments.pdb")
    chain = _parse_first_chain(pdb_path)

    fragments = find_internal_fragments(chain, peptide_bond_max_distance=1.8)

    assert len(fragments) == 2
    assert fragments[0].start_resseq == 1
    assert fragments[0].end_resseq == 1
    assert fragments[1].start_resseq == 5
    assert fragments[1].end_resseq == 5


def test_find_internal_gap_boundaries_returns_expected_boundary(tmp_path: Path) -> None:
    pdb_path = _build_internal_gap_chain(tmp_path / "boundary.pdb")
    chain = _parse_first_chain(pdb_path)

    boundaries = find_internal_gap_boundaries(chain, peptide_bond_max_distance=1.8)

    assert boundaries == [
        InternalGapBoundary(
            chain_id="A",
            left_residue_index=0,
            right_residue_index=1,
            left_resseq=1,
            right_resseq=5,
            left_resname="ALA",
            right_resname="SER",
        )
    ]


def test_cap_internal_gap_boundaries_adds_ace_and_nme(tmp_path: Path) -> None:
    pdb_path = _build_internal_gap_chain(tmp_path / "cap_one_chain.pdb")
    chain = _parse_first_chain(pdb_path)

    capped_chain, summary = cap_internal_gap_boundaries(
        chain=chain,
        serial_counter=[0],
        peptide_bond_max_distance=1.8,
    )

    residue_names = _residue_names(capped_chain)
    residues = _residues(capped_chain)

    assert residue_names == ["ALA", "NME", "ACE", "SER"]
    assert summary.chain_id == "A"
    assert summary.n_fragments == 2
    assert summary.n_internal_boundaries == 1
    assert summary.n_boundaries_capped == 1

    ala_residue = residues[0]
    ser_residue = residues[-1]
    nme_residue = residues[1]
    ace_residue = residues[2]

    assert "OXT" not in _atom_names(ala_residue)
    assert "H1" not in _atom_names(ser_residue)
    assert {"N", "CH3"}.issubset(_atom_names(nme_residue))
    assert {"CH3", "C", "O"}.issubset(_atom_names(ace_residue))


def test_cap_internal_gaps_in_structure_writes_output_and_returns_summary(
    tmp_path: Path,
) -> None:
    input_pdb_path = _build_internal_gap_chain(tmp_path / "input_cap_structure.pdb")
    output_pdb_path = tmp_path / "output_cap_structure.pdb"

    summary_list = cap_internal_gaps_in_structure(
        input_pdb_path=input_pdb_path,
        output_pdb_path=output_pdb_path,
        peptide_bond_max_distance=1.8,
    )

    assert output_pdb_path.exists()
    assert len(summary_list) == 1
    assert summary_list[0].n_internal_boundaries == 1
    assert summary_list[0].n_boundaries_capped == 1

    output_chain = _parse_first_chain(output_pdb_path)
    assert _residue_names(output_chain) == ["ALA", "NME", "ACE", "SER"]


def test_cap_internal_gaps_in_structure_raises_for_missing_input(
    tmp_path: Path,
) -> None:
    with pytest.raises(FileNotFoundError):
        cap_internal_gaps_in_structure(
            input_pdb_path=tmp_path / "missing.pdb",
            output_pdb_path=tmp_path / "out.pdb",
        )


def test_cap_internal_gaps_for_pdb_directory_uses_standard_output_path(
    tmp_path: Path,
) -> None:
    pdb_id = "1ABC"
    pdb_directory = tmp_path / pdb_id
    input_pdb_path = _build_internal_gap_chain(tmp_path / "directory_input.pdb")

    output_pdb_path, summary_list = cap_internal_gaps_for_pdb_directory(
        pdb_directory=pdb_directory,
        pdb_id=pdb_id,
        input_pdb_path=input_pdb_path,
        peptide_bond_max_distance=1.8,
    )

    expected_output = (
        pdb_directory / "components" / f"{pdb_id}_protein_internal_capped.pdb"
    )

    assert output_pdb_path == expected_output
    assert output_pdb_path.exists()
    assert len(summary_list) == 1
    assert summary_list[0].n_internal_boundaries == 1


def test_cap_internal_gaps_for_variant_writes_variant_specific_output(
    tmp_path: Path,
) -> None:
    pdb_id = "1ABC"
    pdb_directory = tmp_path / pdb_id
    input_pdb_path = _build_internal_gap_chain(tmp_path / "variant_input.pdb")

    result = cap_internal_gaps_for_variant(
        pdb_directory=pdb_directory,
        pdb_id=pdb_id,
        variant_label="alphafold_complete",
        input_pdb_path=input_pdb_path,
        peptide_bond_max_distance=1.8,
    )

    expected_output = (
        pdb_directory
        / "components"
        / f"{pdb_id}_protein_internal_capped_alphafold_complete.pdb"
    )

    assert expected_output.exists()
    assert result["internal_capping_success"] is True
    assert result["internal_capping_variant"] == "alphafold_complete"
    assert Path(result["internal_capping_input_path"]) == input_pdb_path.resolve()
    assert Path(result["internal_capping_output_path"]) == expected_output.resolve()
    assert result["internal_capping_chain_count"] == 1
    assert result["internal_capping_total_boundaries"] == 1
    assert result["internal_capping_total_capped_boundaries"] == 1


def test_summarize_capping_results_formats_expected_text() -> None:
    summary_list = [
        ChainCappingSummary(
            chain_id="A",
            n_fragments=2,
            n_internal_boundaries=1,
            n_boundaries_capped=1,
        ),
        ChainCappingSummary(
            chain_id="B",
            n_fragments=1,
            n_internal_boundaries=0,
            n_boundaries_capped=0,
        ),
    ]

    text = summarize_capping_results(summary_list)

    assert "chain=A n_fragments=2 n_internal_boundaries=1 n_boundaries_capped=1" in text
    assert "chain=B n_fragments=1 n_internal_boundaries=0 n_boundaries_capped=0" in text
