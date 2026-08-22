"""Tests for stack_protein_preparation._filler_loop_refine.

Direct coverage for the utility helpers.  The MODELLER-heavy
refine_loops_via_modeller entry point is exercised end-to-end in the
bench (see fruton_bench48_full.py), but its supporting geometry helper
_chirality_d_count is small enough to unit-test here.
"""
from __future__ import annotations

from pathlib import Path

import pytest

from stack_protein_preparation._filler_loop_refine import (
    LoopRefineResult,
    _chirality_d_count,
)


def _atom(serial: int, name: str, resname: str, chain: str, resnum: int,
          x: float, y: float, z: float, element: str = "C") -> str:
    return (
        f"ATOM  {serial:5d} {name:>4} {resname:>3} {chain:1}"
        f"{resnum:4d}    "
        f"{x:8.3f}{y:8.3f}{z:8.3f}"
        f"  1.00  0.00           {element:>2}"
    )


def _write(pdb_path: Path, atoms: list[str]) -> Path:
    pdb_path.parent.mkdir(parents=True, exist_ok=True)
    pdb_path.write_text("\n".join(atoms) + "\nEND\n", encoding="utf-8")
    return pdb_path


# ---------------------------------------------------------- _chirality_d_count


def test_chirality_d_count_l_alanine_returns_zero(tmp_path: Path) -> None:
    # L-Ala: signed volume ((N-CA) x (C-CA)) . (CB-CA) > 0.
    atoms = [
        _atom(1, "N",  "ALA", "A", 1, -0.500, 1.400, 0.000, "N"),
        _atom(2, "CA", "ALA", "A", 1,  0.000, 0.000, 0.000, "C"),
        _atom(3, "C",  "ALA", "A", 1,  1.500, 0.000, 0.000, "C"),
        _atom(4, "CB", "ALA", "A", 1, -0.500, -0.800, -1.200, "C"),
        _atom(5, "O",  "ALA", "A", 1,  2.000, 1.100, 0.000, "O"),
    ]
    p = _write(tmp_path / "l.pdb", atoms)
    assert _chirality_d_count(p) == 0


def test_chirality_d_count_d_alanine_returns_one(tmp_path: Path) -> None:
    # Mirror CB across CA plane -> D configuration.
    atoms = [
        _atom(1, "N",  "ALA", "A", 1, -0.500, 1.400, 0.000, "N"),
        _atom(2, "CA", "ALA", "A", 1,  0.000, 0.000, 0.000, "C"),
        _atom(3, "C",  "ALA", "A", 1,  1.500, 0.000, 0.000, "C"),
        _atom(4, "CB", "ALA", "A", 1, -0.500, -0.800, +1.200, "C"),  # flipped
        _atom(5, "O",  "ALA", "A", 1,  2.000, 1.100, 0.000, "O"),
    ]
    p = _write(tmp_path / "d.pdb", atoms)
    assert _chirality_d_count(p) == 1


def test_chirality_d_count_skips_glycine(tmp_path: Path) -> None:
    # Gly has no CB -> should be skipped.  Even with fabricated CB, GLY
    # is explicitly excluded.
    atoms = [
        _atom(1, "N",  "GLY", "A", 1, -0.500, 1.400, 0.000, "N"),
        _atom(2, "CA", "GLY", "A", 1,  0.000, 0.000, 0.000, "C"),
        _atom(3, "C",  "GLY", "A", 1,  1.500, 0.000, 0.000, "C"),
    ]
    p = _write(tmp_path / "gly.pdb", atoms)
    assert _chirality_d_count(p) == 0


def test_chirality_d_count_skips_incomplete_sidechain(tmp_path: Path) -> None:
    # No CB present -> skipped even if resname != GLY.
    atoms = [
        _atom(1, "N",  "ALA", "A", 1, -0.500, 1.400, 0.000, "N"),
        _atom(2, "CA", "ALA", "A", 1,  0.000, 0.000, 0.000, "C"),
        _atom(3, "C",  "ALA", "A", 1,  1.500, 0.000, 0.000, "C"),
    ]
    p = _write(tmp_path / "no_cb.pdb", atoms)
    assert _chirality_d_count(p) == 0


def test_chirality_d_count_missing_file_returns_minus_one(tmp_path: Path) -> None:
    assert _chirality_d_count(tmp_path / "nonexistent.pdb") == -1


def test_chirality_d_count_skips_hetatm(tmp_path: Path) -> None:
    # HETATM record shouldn't be counted even if it has full backbone+CB
    atoms = [
        f"HETATM    1  N   LIG A   1      -0.500   1.400   0.000  1.00  0.00           N",
        f"HETATM    2  CA  LIG A   1       0.000   0.000   0.000  1.00  0.00           C",
        f"HETATM    3  C   LIG A   1       1.500   0.000   0.000  1.00  0.00           C",
        f"HETATM    4  CB  LIG A   1      -0.500  -0.800   1.200  1.00  0.00           C",
    ]
    p = _write(tmp_path / "het.pdb", atoms)
    assert _chirality_d_count(p) == 0


# ---------------------------------------------------------- LoopRefineResult.to_dict


def test_loop_refine_result_to_dict():
    r = LoopRefineResult(
        output_pdb_path=Path("/tmp/x.pdb"),
        ran=True,
        n_conformers_built=3,
        n_conformers_kept=2,
        best_dope=-1234.5,
        diagnostics=["conformer 1: DOPE=-1234.5"],
    )
    d = r.to_dict()
    assert d["ran"] is True
    assert d["n_conformers_built"] == 3
    assert d["n_conformers_kept"] == 2
    assert d["best_dope"] == -1234.5
    assert d["diagnostics"] == ["conformer 1: DOPE=-1234.5"]


def test_loop_refine_result_default_fallback():
    r = LoopRefineResult(output_pdb_path=Path("/tmp/x.pdb"), ran=False)
    d = r.to_dict()
    assert d["ran"] is False
    assert d["n_conformers_built"] == 0
    assert d["n_conformers_kept"] == 0
    assert d["best_dope"] is None
    assert d["fallback_reason"] is None or isinstance(d["fallback_reason"], str)
