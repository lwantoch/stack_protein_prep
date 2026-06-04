"""Tests for fragment_split — backbone-connectivity fragment detection."""
from __future__ import annotations

from pathlib import Path

import pytest

from stack_protein_preparation.fragment_split import (
    AtomRecord,
    BoundaryRecord,
    FragmentRecord,
    ResidueRecord,
    _infer_pdb_id_from_input_path,
    assess_peptide_connection,
    read_residue_records_from_protein_pdb,
)

DEFAULT_CUTOFF = 1.8


# ---------------------------------------------------------------------------
# _infer_pdb_id_from_input_path
# ---------------------------------------------------------------------------

def test_infer_pdb_id_strips_protein_suffix(tmp_path: Path) -> None:
    p = tmp_path / "3OLL_protein.pdb"
    assert _infer_pdb_id_from_input_path(p) == "3OLL"


def test_infer_pdb_id_strips_protein_midfix(tmp_path: Path) -> None:
    p = tmp_path / "1ABC_protein_monomer.pdb"
    assert _infer_pdb_id_from_input_path(p) == "1ABC"


def test_infer_pdb_id_plain_stem(tmp_path: Path) -> None:
    p = tmp_path / "2AFX_something.pdb"
    assert _infer_pdb_id_from_input_path(p) == "2AFX"


# ---------------------------------------------------------------------------
# ResidueRecord helpers
# ---------------------------------------------------------------------------

def _make_residue(
    chain: str,
    resseq: int,
    resname: str = "ALA",
    atoms: dict[str, tuple[float, float, float]] | None = None,
) -> ResidueRecord:
    r = ResidueRecord(chain_id=chain, resseq=resseq, icode="", resname=resname)
    for name, (x, y, z) in (atoms or {}).items():
        r.add_atom(AtomRecord(
            atom_name=name,
            altloc="",
            x=x,
            y=y,
            z=z,
            line=f"ATOM      1  {name:<4s}{resname:3s} {chain}{resseq:4d}    {x:8.3f}{y:8.3f}{z:8.3f}\n",
        ))
    return r


def test_residue_label_includes_chain_and_resseq() -> None:
    r = _make_residue("B", 42, "GLY")
    assert "B" in r.residue_label
    assert "42" in r.residue_label


def test_residue_label_blank_chain() -> None:
    r = _make_residue("", 1, "ALA")
    assert "(blank)" in r.chain_label


def test_residue_get_preferred_atom_none_when_missing() -> None:
    r = _make_residue("A", 1)
    assert r.get_preferred_atom("CA") is None


def test_residue_atom_count() -> None:
    r = _make_residue("A", 1, atoms={"CA": (0, 0, 0), "CB": (1, 0, 0)})
    assert r.atom_count() == 2


# ---------------------------------------------------------------------------
# FragmentRecord helpers
# ---------------------------------------------------------------------------

def test_fragment_residue_count() -> None:
    residues = [_make_residue("A", i) for i in range(1, 4)]
    frag = FragmentRecord(
        chain_id="A",
        fragment_index_global=1,
        fragment_index_in_chain=1,
        residues=residues,
        output_pdb_path=Path("out.pdb"),
    )
    assert frag.residue_count == 3
    assert frag.start_resseq == 1
    assert frag.end_resseq == 3


# ---------------------------------------------------------------------------
# assess_peptide_connection
# ---------------------------------------------------------------------------

def test_connected_residues_return_true() -> None:
    left = _make_residue("A", 1, atoms={"C": (0.0, 0.0, 0.0)})
    right = _make_residue("A", 2, atoms={"N": (1.3, 0.0, 0.0)})
    connected, boundary = assess_peptide_connection(left, right, DEFAULT_CUTOFF)
    assert connected is True
    assert boundary is None


def test_too_far_residues_return_boundary(tmp_path: Path) -> None:
    left = _make_residue("A", 1, atoms={"C": (0.0, 0.0, 0.0)})
    right = _make_residue("A", 5, atoms={"N": (3.0, 0.0, 0.0)})
    connected, boundary = assess_peptide_connection(left, right, DEFAULT_CUTOFF)
    assert connected is False
    assert boundary is not None
    assert boundary.cn_distance is not None
    assert boundary.cn_distance > DEFAULT_CUTOFF


def test_missing_left_c_returns_boundary() -> None:
    left = _make_residue("A", 1)   # no C atom
    right = _make_residue("A", 2, atoms={"N": (1.3, 0.0, 0.0)})
    connected, boundary = assess_peptide_connection(left, right, DEFAULT_CUTOFF)
    assert connected is False
    assert boundary is not None
    assert boundary.reason == "missing_left_backbone_C"


def test_missing_right_n_returns_boundary() -> None:
    left = _make_residue("A", 1, atoms={"C": (0.0, 0.0, 0.0)})
    right = _make_residue("A", 2)  # no N atom
    connected, boundary = assess_peptide_connection(left, right, DEFAULT_CUTOFF)
    assert connected is False
    assert boundary is not None
    assert boundary.reason == "missing_right_backbone_N"


def test_boundary_records_residue_jump() -> None:
    left = _make_residue("A", 1, atoms={"C": (0.0, 0.0, 0.0)})
    right = _make_residue("A", 10, atoms={"N": (5.0, 0.0, 0.0)})
    _, boundary = assess_peptide_connection(left, right, DEFAULT_CUTOFF)
    assert boundary is not None
    assert boundary.residue_jump == 9


# ---------------------------------------------------------------------------
# read_residue_records_from_protein_pdb
# ---------------------------------------------------------------------------

def _atom_line(serial: int, atom: str, resname: str, chain: str, resseq: int,
               x: float, y: float, z: float) -> str:
    return (
        f"ATOM  {serial:5d}  {atom:<4s}{resname:3s} {chain}{resseq:4d}    "
        f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00           C  \n"
    )


def _write_pdb(path: Path, lines: list[str]) -> Path:
    path.write_text("".join(lines) + "END\n", encoding="utf-8")
    return path


def test_reads_two_residues(tmp_path: Path) -> None:
    pdb = _write_pdb(tmp_path / "test.pdb", [
        _atom_line(1, "CA", "ALA", "A", 1, 1.0, 0.0, 0.0),
        _atom_line(2, "CA", "GLY", "A", 2, 2.0, 0.0, 0.0),
    ])
    residues = read_residue_records_from_protein_pdb(pdb)
    assert len(residues) == 2
    assert residues[0].resname == "ALA"
    assert residues[1].resname == "GLY"


def test_ignores_hetatm_lines(tmp_path: Path) -> None:
    lines = [
        _atom_line(1, "CA", "ALA", "A", 1, 1.0, 0.0, 0.0),
        f"HETATM{2:5d}  ZN  ZN  A{2:4d}       4.0     5.0     6.0   1.00  0.00          ZN  \n",
    ]
    pdb = _write_pdb(tmp_path / "test.pdb", lines)
    residues = read_residue_records_from_protein_pdb(pdb)
    assert len(residues) == 1
    assert residues[0].resname == "ALA"


def test_multiple_atoms_same_residue_grouped(tmp_path: Path) -> None:
    pdb = _write_pdb(tmp_path / "test.pdb", [
        _atom_line(1, "N", "ALA", "A", 1, 0.0, 0.0, 0.0),
        _atom_line(2, "CA", "ALA", "A", 1, 1.0, 0.0, 0.0),
        _atom_line(3, "C", "ALA", "A", 1, 2.0, 0.0, 0.0),
    ])
    residues = read_residue_records_from_protein_pdb(pdb)
    assert len(residues) == 1
    assert residues[0].atom_count() == 3
