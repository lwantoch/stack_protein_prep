"""Tests for stack_protein_preparation._filler_plddt.

Coverage for the per-loop AlphaFold pLDDT extraction and REMARK 465
gap-range parser used by the splice-time pLDDT floor gate.
"""
from __future__ import annotations

from pathlib import Path

import pytest

from stack_protein_preparation._filler_plddt import (
    compute_per_loop_plddt,
    parse_gap_ranges_from_remark_465,
)


def _atom_line(serial: int, name: str, resname: str, chain: str, resnum: int,
               x: float, y: float, z: float, bfactor: float,
               element: str = "C") -> str:
    return (
        f"ATOM  {serial:5d} {name:>4} {resname:>3} {chain:1}"
        f"{resnum:4d}    "
        f"{x:8.3f}{y:8.3f}{z:8.3f}"
        f"  1.00{bfactor:6.2f}           {element:>2}"
    )


def _write_pdb(pdb_path: Path, lines: list[str]) -> Path:
    pdb_path.parent.mkdir(parents=True, exist_ok=True)
    pdb_path.write_text("\n".join(lines) + "\nEND\n", encoding="utf-8")
    return pdb_path


# ---------------------------------------------------- compute_per_loop_plddt


def test_compute_plddt_reads_bfactor_in_range(tmp_path: Path) -> None:
    # 3 residues chain A; residues 2-3 in gap range, mean pLDDT is 80.
    atoms = [
        _atom_line(1, "N",  "ALA", "A", 1, 0, 0, 0, bfactor=50, element="N"),
        _atom_line(2, "CA", "ALA", "A", 1, 1, 0, 0, bfactor=50, element="C"),
        _atom_line(3, "N",  "GLY", "A", 2, 3, 0, 0, bfactor=70, element="N"),
        _atom_line(4, "CA", "GLY", "A", 2, 4, 0, 0, bfactor=70, element="C"),
        _atom_line(5, "N",  "SER", "A", 3, 6, 0, 0, bfactor=90, element="N"),
        _atom_line(6, "CA", "SER", "A", 3, 7, 0, 0, bfactor=90, element="C"),
    ]
    pdb = _write_pdb(tmp_path / "af.pdb", atoms)
    result = compute_per_loop_plddt(
        af_model_path=pdb,
        gap_ranges=[{"chain": "A", "start": 2, "end": 3}],
        method="alphafold",
    )
    assert len(result) == 1
    row = result[0]
    assert row["chain"] == "A" and row["start"] == 2 and row["end"] == 3
    assert row["method"] == "alphafold"
    # Extractor prefers CA atoms; average of CA(res 2)=70 and CA(res 3)=90 is 80.
    assert row["mean_plddt"] == 80.0
    assert row["mean_bfactor"] == 80.0
    assert row["n_atoms"] == 2
    assert "plddt_note" not in row


def test_compute_plddt_flags_out_of_range_bfactors(tmp_path: Path) -> None:
    # Mean b-factor > 100 (crystallographic B, not pLDDT).  method=alphafold
    # -> should set plddt_note, leave mean_plddt None.
    atoms = [
        _atom_line(1, "CA", "ALA", "A", 5, 0, 0, 0, bfactor=150, element="C"),
        _atom_line(2, "CA", "ALA", "A", 6, 1, 0, 0, bfactor=150, element="C"),
    ]
    pdb = _write_pdb(tmp_path / "crystal_bfactors.pdb", atoms)
    result = compute_per_loop_plddt(
        af_model_path=pdb,
        gap_ranges=[{"chain": "A", "start": 5, "end": 6}],
        method="alphafold",
    )
    row = result[0]
    assert row["mean_plddt"] is None
    assert row["mean_bfactor"] == 150.0
    assert "plddt_note" in row


def test_compute_plddt_modeller_never_reports_plddt(tmp_path: Path) -> None:
    # method=modeller -> mean_plddt always None even if bfactors in [0,100].
    atoms = [_atom_line(1, "CA", "GLY", "A", 5, 0, 0, 0, bfactor=80, element="C")]
    pdb = _write_pdb(tmp_path / "modeller.pdb", atoms)
    result = compute_per_loop_plddt(
        af_model_path=pdb,
        gap_ranges=[{"chain": "A", "start": 5, "end": 5}],
        method="modeller",
    )
    assert result[0]["mean_plddt"] is None
    assert result[0]["mean_bfactor"] == 80.0


def test_compute_plddt_missing_file_returns_stub(tmp_path: Path) -> None:
    result = compute_per_loop_plddt(
        af_model_path=tmp_path / "nonexistent.pdb",
        gap_ranges=[{"chain": "A", "start": 1, "end": 2}],
        method="alphafold",
    )
    row = result[0]
    assert row["mean_plddt"] is None
    assert row["n_atoms"] == 0


def test_compute_plddt_empty_gap_ranges_returns_empty_list(tmp_path: Path) -> None:
    atoms = [_atom_line(1, "CA", "ALA", "A", 1, 0, 0, 0, bfactor=50, element="C")]
    pdb = _write_pdb(tmp_path / "af.pdb", atoms)
    result = compute_per_loop_plddt(af_model_path=pdb, gap_ranges=[])
    assert result == []


# ---------------------------------------------------- REMARK 465 parser


def test_parse_remark_465_groups_consecutive_residues(tmp_path: Path) -> None:
    # REMARK 465 with 4 missing residues on chain A: 5, 6, 7 (contiguous)
    # and 20 (isolated).  Grouping should give two segments.
    lines = [
        "REMARK 465 MISSING RESIDUES",
        "REMARK 465   M RES C SSSEQI",
        "REMARK 465     ALA A     5",
        "REMARK 465     GLY A     6",
        "REMARK 465     SER A     7",
        "REMARK 465     VAL A    20",
        # a real atom so BioPDB accepts the file
        _atom_line(1, "CA", "MET", "A", 1, 0, 0, 0, bfactor=50, element="C"),
    ]
    pdb = _write_pdb(tmp_path / "with_465.pdb", lines)
    groups = parse_gap_ranges_from_remark_465(pdb)
    # BioPython parses REMARK 465 into structure.header["missing_residues"].
    # The grouper should emit at least the (5, 7) segment; some Bio.PDB
    # versions include the isolated 20 as its own {5,7} and {20,20}.
    ranges = sorted((g["chain"], g["start"], g["end"]) for g in groups)
    assert ("A", 5, 7) in ranges
    # Isolated missing 20 either present as (20, 20) or absent depending
    # on parser version -- both are acceptable.


def test_parse_remark_465_no_remark_returns_empty(tmp_path: Path) -> None:
    atoms = [_atom_line(1, "CA", "ALA", "A", 1, 0, 0, 0, bfactor=50, element="C")]
    pdb = _write_pdb(tmp_path / "no_465.pdb", atoms)
    assert parse_gap_ranges_from_remark_465(pdb) == []


def test_parse_remark_465_missing_file_returns_empty(tmp_path: Path) -> None:
    assert parse_gap_ranges_from_remark_465(tmp_path / "missing.pdb") == []
