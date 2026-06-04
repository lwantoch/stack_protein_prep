"""Tests for insertion_codes — pure-Python PDB insertion-code helpers."""
from __future__ import annotations

from pathlib import Path

import pytest

from stack_protein_preparation.insertion_codes import (
    _infer_pdb_id_from_path,
    _python_delinsertion,
    pdb_has_insertion_codes,
    process_pdb_for_delinsertion,
)


# ---------------------------------------------------------------------------
# PDB line helpers
# ---------------------------------------------------------------------------

def _atom_line(
    serial: int,
    resname: str,
    chain: str,
    resseq: int,
    icode: str = " ",
) -> str:
    return (
        f"ATOM  {serial:5d}  CA  {resname:3s} {chain}{resseq:4d}{icode}   "
        f"   1.000   2.000   3.000  1.00  0.00           C  \n"
    )


def _write_pdb(tmp_path: Path, name: str, lines: list[str]) -> Path:
    p = tmp_path / name
    p.write_text("".join(lines) + "END\n", encoding="utf-8")
    return p


# ---------------------------------------------------------------------------
# pdb_has_insertion_codes
# ---------------------------------------------------------------------------

def test_no_insertion_codes_returns_false(tmp_path: Path) -> None:
    pdb = _write_pdb(tmp_path, "clean.pdb", [
        _atom_line(1, "ALA", "A", 1),
        _atom_line(2, "GLY", "A", 2),
    ])
    assert pdb_has_insertion_codes(pdb) is False


def test_insertion_code_detected(tmp_path: Path) -> None:
    pdb = _write_pdb(tmp_path, "ins.pdb", [
        _atom_line(1, "ALA", "A", 10),
        _atom_line(2, "GLY", "A", 10, icode="A"),
    ])
    assert pdb_has_insertion_codes(pdb) is True


def test_non_atom_lines_ignored(tmp_path: Path) -> None:
    pdb = _write_pdb(tmp_path, "remark.pdb", [
        "REMARK  10A something\n",
        _atom_line(1, "ALA", "A", 1),
    ])
    assert pdb_has_insertion_codes(pdb) is False


# ---------------------------------------------------------------------------
# _python_delinsertion
# ---------------------------------------------------------------------------

def _build_pdb_text(*lines: str) -> str:
    return "".join(lines) + "END\n"


def test_delinsertion_removes_insertion_code(tmp_path: Path) -> None:
    text = _build_pdb_text(
        _atom_line(1, "ALA", "A", 10),
        _atom_line(2, "GLY", "A", 10, icode="A"),
        _atom_line(3, "SER", "A", 11),
    )
    result = _python_delinsertion(text)
    # insertion code column (index 26) should be space for all lines
    for line in result.splitlines():
        if line.startswith("ATOM") or line.startswith("HETATM"):
            assert line[26] == " ", f"Non-space insertion code in: {line!r}"


def test_delinsertion_renumbers_sequentially(tmp_path: Path) -> None:
    text = _build_pdb_text(
        _atom_line(1, "ALA", "A", 5),
        _atom_line(2, "GLY", "A", 5, icode="A"),
        _atom_line(3, "SER", "A", 6),
    )
    result = _python_delinsertion(text)
    atom_lines = [l for l in result.splitlines() if l.startswith("ATOM")]
    resnums = [int(l[22:26]) for l in atom_lines]
    assert resnums == [5, 6, 7]


def test_delinsertion_leaves_clean_pdb_unchanged(tmp_path: Path) -> None:
    text = _build_pdb_text(
        _atom_line(1, "ALA", "A", 1),
        _atom_line(2, "GLY", "A", 2),
    )
    result = _python_delinsertion(text)
    atom_lines = [l for l in result.splitlines() if l.startswith("ATOM")]
    assert [int(l[22:26]) for l in atom_lines] == [1, 2]


# ---------------------------------------------------------------------------
# _infer_pdb_id_from_path
# ---------------------------------------------------------------------------

def test_infer_pdb_id_from_proteins_layout(tmp_path: Path) -> None:
    fake_path = tmp_path / "proteins" / "3OLL" / "3OLL.pdb"
    fake_path.parent.mkdir(parents=True)
    fake_path.touch()
    assert _infer_pdb_id_from_path(fake_path) == "3OLL"


def test_infer_pdb_id_falls_back_to_parent(tmp_path: Path) -> None:
    fake_path = tmp_path / "2AFX" / "input.pdb"
    fake_path.parent.mkdir(parents=True)
    fake_path.touch()
    assert _infer_pdb_id_from_path(fake_path) == "2AFX"


# ---------------------------------------------------------------------------
# process_pdb_for_delinsertion
# ---------------------------------------------------------------------------

def test_process_copies_clean_pdb(tmp_path: Path) -> None:
    pdb = _write_pdb(tmp_path, "1ABC.pdb", [_atom_line(1, "ALA", "A", 1)])
    out = tmp_path / "out" / "1ABC_delinsertion.pdb"
    result = process_pdb_for_delinsertion(pdb, out)
    assert result["had_insertion_codes"] is False
    assert result["status"] == "none"
    assert out.is_file()


def test_process_runs_delinsertion_when_icode_present(tmp_path: Path) -> None:
    pdb = _write_pdb(tmp_path, "1INS.pdb", [
        _atom_line(1, "ALA", "A", 10),
        _atom_line(2, "GLY", "A", 10, icode="A"),
    ])
    out = tmp_path / "out" / "1INS_delinsertion.pdb"
    result = process_pdb_for_delinsertion(pdb, out)
    assert result["had_insertion_codes"] is True
    assert out.is_file()
    content = out.read_text()
    for line in content.splitlines():
        if line.startswith("ATOM") or line.startswith("HETATM"):
            assert line[26] == " "
