"""
/home/grheco/repositorios/stack_protein_prep/tests/test_gaps.py

Tests for internal residue-numbering gap detection.

Responsibilities
----------------
- verify simple gap detection from minimal PDB inputs
- verify chain-wise grouping
- verify summary statistics
- verify that non-protein residues are ignored
- verify that protein-like nonstandard residues do not create artificial gaps

Notes
-----
- these tests write temporary minimal PDB files
- openmm.app.PDBFile is used by the production code, so the PDB text must be
  valid enough for OpenMM to parse
"""

from __future__ import annotations

from pathlib import Path

from stack_protein_preparation.gaps import (
    debug_excluded_residues,
    gaps_by_chain_for_pdb,
    gaps_for_pdb,
    summarize_gaps,
)


def _make_atom_line(
    serial: int,
    atom_name: str,
    resname: str,
    chain_id: str,
    resseq: int,
    x: float,
    y: float,
    z: float,
    element: str,
) -> str:
    """
    Build one minimally valid ATOM/HETATM-style PDB line.

    This keeps the formatting explicit so temporary test PDBs stay readable.
    """
    return (
        f"ATOM  {serial:5d} {atom_name:^4s} {resname:>3s} {chain_id:1s}"
        f"{resseq:4d}    "
        f"{x:8.3f}{y:8.3f}{z:8.3f}"
        f"{1.00:6.2f}{20.00:6.2f}          "
        f"{element:>2s}"
    )


def _make_hetatm_line(
    serial: int,
    atom_name: str,
    resname: str,
    chain_id: str,
    resseq: int,
    x: float,
    y: float,
    z: float,
    element: str,
) -> str:
    """
    Build one minimally valid HETATM line.
    """
    return (
        f"HETATM{serial:5d} {atom_name:^4s} {resname:>3s} {chain_id:1s}"
        f"{resseq:4d}    "
        f"{x:8.3f}{y:8.3f}{z:8.3f}"
        f"{1.00:6.2f}{20.00:6.2f}          "
        f"{element:>2s}"
    )


def _write_test_pdb(pdb_path: Path, lines: list[str]) -> Path:
    """
    Write a minimal PDB file with END at the end.
    """
    pdb_path.write_text("\n".join(lines + ["END", ""]), encoding="utf-8")
    return pdb_path


def test_gaps_for_pdb_returns_empty_when_numbering_is_continuous(
    tmp_path: Path,
) -> None:
    """
    No numbering discontinuity -> no gaps.
    """
    pdb_path = tmp_path / "continuous.pdb"

    lines = [
        _make_atom_line(1, "N", "ALA", "A", 1, 0.0, 0.0, 0.0, "N"),
        _make_atom_line(2, "CA", "ALA", "A", 1, 1.0, 0.0, 0.0, "C"),
        _make_atom_line(3, "N", "GLY", "A", 2, 2.0, 0.0, 0.0, "N"),
        _make_atom_line(4, "CA", "GLY", "A", 2, 3.0, 0.0, 0.0, "C"),
        _make_atom_line(5, "N", "SER", "A", 3, 4.0, 0.0, 0.0, "N"),
        _make_atom_line(6, "CA", "SER", "A", 3, 5.0, 0.0, 0.0, "C"),
    ]
    _write_test_pdb(pdb_path, lines)

    assert gaps_for_pdb(pdb_path) == {}


def test_gaps_for_pdb_detects_single_gap(tmp_path: Path) -> None:
    """
    Residues 1 and 4 -> missing 2-3.
    """
    pdb_path = tmp_path / "single_gap.pdb"

    lines = [
        _make_atom_line(1, "N", "ALA", "A", 1, 0.0, 0.0, 0.0, "N"),
        _make_atom_line(2, "CA", "ALA", "A", 1, 1.0, 0.0, 0.0, "C"),
        _make_atom_line(3, "N", "GLY", "A", 4, 2.0, 0.0, 0.0, "N"),
        _make_atom_line(4, "CA", "GLY", "A", 4, 3.0, 0.0, 0.0, "C"),
    ]
    _write_test_pdb(pdb_path, lines)

    gaps = gaps_for_pdb(pdb_path)

    assert gaps == [
        {
            "chain": "A",
            "prev_res": 1,
            "prev_resname": "ALA",
            "start": 2,
            "end": 3,
            "next_res": 4,
            "next_resname": "GLY",
            "gap_size": 2,
        }
    ]


def test_gaps_for_pdb_detects_multiple_gaps_in_one_chain(tmp_path: Path) -> None:
    """
    Residues 10, 11, 15, 20 -> gaps 12-14 and 16-19.
    """
    pdb_path = tmp_path / "multiple_gaps.pdb"

    lines = [
        _make_atom_line(1, "N", "ALA", "A", 10, 0.0, 0.0, 0.0, "N"),
        _make_atom_line(2, "CA", "ALA", "A", 10, 1.0, 0.0, 0.0, "C"),
        _make_atom_line(3, "N", "GLY", "A", 11, 2.0, 0.0, 0.0, "N"),
        _make_atom_line(4, "CA", "GLY", "A", 11, 3.0, 0.0, 0.0, "C"),
        _make_atom_line(5, "N", "SER", "A", 15, 4.0, 0.0, 0.0, "N"),
        _make_atom_line(6, "CA", "SER", "A", 15, 5.0, 0.0, 0.0, "C"),
        _make_atom_line(7, "N", "VAL", "A", 20, 6.0, 0.0, 0.0, "N"),
        _make_atom_line(8, "CA", "VAL", "A", 20, 7.0, 0.0, 0.0, "C"),
    ]
    _write_test_pdb(pdb_path, lines)

    gaps = gaps_for_pdb(pdb_path)

    assert gaps == [
        {
            "chain": "A",
            "prev_res": 11,
            "prev_resname": "GLY",
            "start": 12,
            "end": 14,
            "next_res": 15,
            "next_resname": "SER",
            "gap_size": 3,
        },
        {
            "chain": "A",
            "prev_res": 15,
            "prev_resname": "SER",
            "start": 16,
            "end": 19,
            "next_res": 20,
            "next_resname": "VAL",
            "gap_size": 4,
        },
    ]


def test_gaps_by_chain_for_pdb_separates_chains(tmp_path: Path) -> None:
    """
    Gaps must be detected independently per chain.
    """
    pdb_path = tmp_path / "two_chains.pdb"

    lines = [
        _make_atom_line(1, "N", "ALA", "A", 1, 0.0, 0.0, 0.0, "N"),
        _make_atom_line(2, "CA", "ALA", "A", 1, 1.0, 0.0, 0.0, "C"),
        _make_atom_line(3, "N", "GLY", "A", 3, 2.0, 0.0, 0.0, "N"),
        _make_atom_line(4, "CA", "GLY", "A", 3, 3.0, 0.0, 0.0, "C"),
        _make_atom_line(5, "N", "SER", "B", 10, 4.0, 0.0, 0.0, "N"),
        _make_atom_line(6, "CA", "SER", "B", 10, 5.0, 0.0, 0.0, "C"),
        _make_atom_line(7, "N", "VAL", "B", 11, 6.0, 0.0, 0.0, "N"),
        _make_atom_line(8, "CA", "VAL", "B", 11, 7.0, 0.0, 0.0, "C"),
    ]
    _write_test_pdb(pdb_path, lines)

    grouped = gaps_by_chain_for_pdb(pdb_path)

    assert grouped == {
        "A": [
            {
                "chain": "A",
                "prev_res": 1,
                "prev_resname": "ALA",
                "start": 2,
                "end": 2,
                "next_res": 3,
                "next_resname": "GLY",
                "gap_size": 1,
            }
        ]
    }


def test_nonprotein_residues_are_excluded_from_gap_detection(tmp_path: Path) -> None:
    """
    Waters / ligands should not count as protein residues and should appear in debug output.
    """
    pdb_path = tmp_path / "excluded_residues.pdb"

    lines = [
        _make_atom_line(1, "N", "ALA", "A", 1, 0.0, 0.0, 0.0, "N"),
        _make_atom_line(2, "CA", "ALA", "A", 1, 1.0, 0.0, 0.0, "C"),
        _make_hetatm_line(3, "O", "HOH", "A", 2, 2.0, 0.0, 0.0, "O"),
        _make_hetatm_line(4, "ZN", "ZN", "A", 3, 3.0, 0.0, 0.0, "Zn"),
        _make_atom_line(5, "N", "GLY", "A", 4, 4.0, 0.0, 0.0, "N"),
        _make_atom_line(6, "CA", "GLY", "A", 4, 5.0, 0.0, 0.0, "C"),
    ]
    _write_test_pdb(pdb_path, lines)

    gaps = gaps_for_pdb(pdb_path)
    excluded = debug_excluded_residues(pdb_path)

    assert gaps == [
        {
            "chain": "A",
            "prev_res": 1,
            "prev_resname": "ALA",
            "start": 2,
            "end": 3,
            "next_res": 4,
            "next_resname": "GLY",
            "gap_size": 2,
        }
    ]
    assert excluded == {"A": [("HOH", "2"), ("ZN", "3")]}


def test_protein_like_nonstandard_residue_does_not_create_artificial_gap(
    tmp_path: Path,
) -> None:
    """
    MSE should count as protein-like residue and therefore bridge numbering normally.
    """
    pdb_path = tmp_path / "mse_is_protein_like.pdb"

    lines = [
        _make_atom_line(1, "N", "ALA", "A", 1, 0.0, 0.0, 0.0, "N"),
        _make_atom_line(2, "CA", "ALA", "A", 1, 1.0, 0.0, 0.0, "C"),
        _make_atom_line(3, "N", "MSE", "A", 2, 2.0, 0.0, 0.0, "N"),
        _make_atom_line(4, "CA", "MSE", "A", 2, 3.0, 0.0, 0.0, "C"),
        _make_atom_line(5, "N", "GLY", "A", 3, 4.0, 0.0, 0.0, "N"),
        _make_atom_line(6, "CA", "GLY", "A", 3, 5.0, 0.0, 0.0, "C"),
    ]
    _write_test_pdb(pdb_path, lines)

    assert gaps_for_pdb(pdb_path) == []


def test_summarize_gaps_returns_expected_statistics(tmp_path: Path) -> None:
    """
    Summary should report counts and aggregate sizes correctly.
    """
    pdb_path = tmp_path / "summary_case.pdb"

    lines = [
        _make_atom_line(1, "N", "ALA", "A", 5, 0.0, 0.0, 0.0, "N"),
        _make_atom_line(2, "CA", "ALA", "A", 5, 1.0, 0.0, 0.0, "C"),
        _make_atom_line(3, "N", "GLY", "A", 8, 2.0, 0.0, 0.0, "N"),
        _make_atom_line(4, "CA", "GLY", "A", 8, 3.0, 0.0, 0.0, "C"),
        _make_atom_line(5, "N", "SER", "B", 20, 4.0, 0.0, 0.0, "N"),
        _make_atom_line(6, "CA", "SER", "B", 20, 5.0, 0.0, 0.0, "C"),
        _make_atom_line(7, "N", "VAL", "B", 22, 6.0, 0.0, 0.0, "N"),
        _make_atom_line(8, "CA", "VAL", "B", 22, 7.0, 0.0, 0.0, "C"),
    ]
    _write_test_pdb(pdb_path, lines)

    summary = summarize_gaps(pdb_path)

    assert summary == {
        "has_gaps": True,
        "n_gaps": 2,
        "max_gap_size": 2,
        "total_missing_residues": 3,
        "gap_sizes": [2, 1],
        "chains_with_gaps": ["A", "B"],
    }
