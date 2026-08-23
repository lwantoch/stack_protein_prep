"""Tests for stack_protein_preparation._blind_test."""
from __future__ import annotations

from pathlib import Path

import pytest

from stack_protein_preparation import _blind_test as bt


def _write_crystal(path: Path) -> Path:
    """3-residue crystal PDB, chain A residues 1-3, all ALA."""
    text = (
        "ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00  0.00           N\n"
        "ATOM      2  CA  ALA A   1       1.458   0.000   0.000  1.00  0.00           C\n"
        "ATOM      3  C   ALA A   1       2.009   1.420   0.000  1.00  0.00           C\n"
        "ATOM      4  O   ALA A   1       1.251   2.390   0.000  1.00  0.00           O\n"
        "ATOM      5  N   ALA A   2       3.332   1.540   0.000  1.00  0.00           N\n"
        "ATOM      6  CA  ALA A   2       4.000   2.828   0.000  1.00  0.00           C\n"
        "ATOM      7  C   ALA A   2       5.500   2.700   0.000  1.00  0.00           C\n"
        "ATOM      8  O   ALA A   2       6.000   1.600   0.000  1.00  0.00           O\n"
        "ATOM      9  N   ALA A   3       6.200   3.700   0.000  1.00  0.00           N\n"
        "ATOM     10  CA  ALA A   3       7.500   4.000   0.000  1.00  0.00           C\n"
        "ATOM     11  C   ALA A   3       8.100   5.400   0.000  1.00  0.00           C\n"
        "ATOM     12  O   ALA A   3       7.500   6.500   0.000  1.00  0.00           O\n"
        "TER\nEND\n"
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text)
    return path


def _write_filled(
    path: Path, ca2_offset: tuple[float, float, float] = (0.0, 0.0, 0.0),
) -> Path:
    """3-residue filled PDB, with residue-2 CA shifted by ca2_offset."""
    dx, dy, dz = ca2_offset
    text = (
        "ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00  0.00           N\n"
        "ATOM      2  CA  ALA A   1       1.458   0.000   0.000  1.00  0.00           C\n"
        "ATOM      3  C   ALA A   1       2.009   1.420   0.000  1.00  0.00           C\n"
        "ATOM      4  O   ALA A   1       1.251   2.390   0.000  1.00  0.00           O\n"
        "ATOM      5  N   ALA A   2       3.332   1.540   0.000  1.00  0.00           N\n"
        f"ATOM      6  CA  ALA A   2     {4.000 + dx:7.3f} {2.828 + dy:7.3f} {0.000 + dz:7.3f}  1.00  0.00           C\n"
        "ATOM      7  C   ALA A   2       5.500   2.700   0.000  1.00  0.00           C\n"
        "ATOM      8  O   ALA A   2       6.000   1.600   0.000  1.00  0.00           O\n"
        "ATOM      9  N   ALA A   3       6.200   3.700   0.000  1.00  0.00           N\n"
        "ATOM     10  CA  ALA A   3       7.500   4.000   0.000  1.00  0.00           C\n"
        "ATOM     11  C   ALA A   3       8.100   5.400   0.000  1.00  0.00           C\n"
        "ATOM     12  O   ALA A   3       7.500   6.500   0.000  1.00  0.00           O\n"
        "TER\nEND\n"
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text)
    return path


# ---------------------------------------------------------------------------
# GapRange
# ---------------------------------------------------------------------------


def test_gaprange_contains():
    r = bt.GapRange("A", 5, 10)
    assert r.contains("A", 5)
    assert r.contains("A", 10)
    assert not r.contains("A", 4)
    assert not r.contains("A", 11)
    assert not r.contains("B", 7)


def test_gaprange_size():
    assert bt.GapRange("A", 5, 10).size() == 6
    assert bt.GapRange("A", 5, 5).size() == 1


# ---------------------------------------------------------------------------
# mask_crystal_pdb
# ---------------------------------------------------------------------------


def test_mask_removes_matching_residues(tmp_path: Path):
    crystal = _write_crystal(tmp_path / "crystal.pdb")
    masked = tmp_path / "masked.pdb"
    bt.mask_crystal_pdb(crystal, [bt.GapRange("A", 2, 2)], masked)
    text = masked.read_text()
    assert "ALA A   1" in text
    assert "ALA A   2" not in text
    assert "ALA A   3" in text


def test_mask_preserves_non_atom_lines(tmp_path: Path):
    crystal = tmp_path / "crystal.pdb"
    crystal.write_text(
        "HEADER   TEST\n"
        "REMARK   nothing to see\n"
        "ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00  0.00           C\n"
        "TER\nEND\n"
    )
    masked = tmp_path / "masked.pdb"
    bt.mask_crystal_pdb(crystal, [bt.GapRange("A", 1, 1)], masked)
    text = masked.read_text()
    assert "HEADER" in text
    assert "REMARK" in text
    assert "TER" in text
    assert "END" in text
    assert "ATOM " not in text  # only ATOM row was removed


def test_mask_empty_range_list_is_identity(tmp_path: Path):
    crystal = _write_crystal(tmp_path / "crystal.pdb")
    masked = tmp_path / "masked.pdb"
    bt.mask_crystal_pdb(crystal, [], masked)
    assert masked.read_text() == crystal.read_text()


def test_mask_range_missing_from_crystal_is_noop(tmp_path: Path):
    crystal = _write_crystal(tmp_path / "crystal.pdb")
    masked = tmp_path / "masked.pdb"
    bt.mask_crystal_pdb(crystal, [bt.GapRange("Z", 99, 100)], masked)
    assert masked.read_text() == crystal.read_text()


# ---------------------------------------------------------------------------
# score_blind_fill
# ---------------------------------------------------------------------------


def test_score_blind_fill_zero_offset_gives_zero_rmsd(tmp_path: Path):
    crystal = _write_crystal(tmp_path / "crystal.pdb")
    filled = _write_filled(tmp_path / "filled.pdb", ca2_offset=(0, 0, 0))
    r = bt.score_blind_fill(crystal, filled, [bt.GapRange("A", 2, 2)])
    assert r.ran is True
    assert r.passed is True
    assert len(r.ranges) == 1
    assert r.ranges[0].n_matched() == 1
    assert abs(r.overall_ca_rmsd_A()) < 1e-6


def test_score_blind_fill_offset_gives_expected_ca_distance(tmp_path: Path):
    crystal = _write_crystal(tmp_path / "crystal.pdb")
    # Shift CA of residue 2 by (0.5, 0, 0) → Cα distance should be 0.5
    filled = _write_filled(tmp_path / "filled.pdb", ca2_offset=(0.5, 0, 0))
    r = bt.score_blind_fill(crystal, filled, [bt.GapRange("A", 2, 2)])
    d = r.ranges[0].residues[0].ca_distance_A
    assert abs(d - 0.5) < 1e-3


def test_score_blind_fill_missing_filled_residue_marks_unmatched(tmp_path: Path):
    crystal = _write_crystal(tmp_path / "crystal.pdb")
    # Emit a filled PDB that only has residues 1 and 3 (residue 2 missing)
    filled = tmp_path / "filled.pdb"
    filled.write_text(
        "ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00  0.00           N\n"
        "ATOM      2  CA  ALA A   1       1.458   0.000   0.000  1.00  0.00           C\n"
        "ATOM      9  N   ALA A   3       6.200   3.700   0.000  1.00  0.00           N\n"
        "ATOM     10  CA  ALA A   3       7.500   4.000   0.000  1.00  0.00           C\n"
        "TER\nEND\n"
    )
    r = bt.score_blind_fill(crystal, filled, [bt.GapRange("A", 2, 2)])
    assert r.ran is True
    # Range scored, but the residue has ca_distance_A = None
    assert r.ranges[0].residues[0].ca_distance_A is None
    assert r.ranges[0].n_matched() == 0
    # overall_ca_rmsd_A returns None when nothing matched
    assert r.overall_ca_rmsd_A() is None


def test_score_blind_fill_missing_crystal_range_is_skipped(tmp_path: Path):
    """If crystal itself lacked the residue, we skip silently."""
    crystal = _write_crystal(tmp_path / "crystal.pdb")
    filled = _write_filled(tmp_path / "filled.pdb")
    r = bt.score_blind_fill(crystal, filled, [bt.GapRange("A", 99, 100)])
    # Range appears but residues list is empty (crystal lacked them)
    assert len(r.ranges) == 1
    assert r.ranges[0].residues == []


def test_score_blind_fill_fail_open_missing_files(tmp_path: Path):
    r = bt.score_blind_fill(
        tmp_path / "nope.pdb", tmp_path / "nope.pdb",
        [bt.GapRange("A", 1, 5)],
    )
    assert r.ran is False
    assert r.passed is False
    assert "not found" in r.fallback_reason


def test_score_blind_fill_backbone_rmsd_matches_ca_at_backbone_offset(tmp_path: Path):
    """When only CA moves, backbone RMSD < CA distance (N/C/O share full coords)."""
    crystal = _write_crystal(tmp_path / "crystal.pdb")
    filled = _write_filled(tmp_path / "filled.pdb", ca2_offset=(0.5, 0, 0))
    r = bt.score_blind_fill(crystal, filled, [bt.GapRange("A", 2, 2)])
    entry = r.ranges[0].residues[0]
    # backbone_rmsd = sqrt((0.5^2 + 0 + 0 + 0) / 4) = 0.25
    assert entry.backbone_rmsd_A == pytest.approx(0.25, abs=1e-3)


# ---------------------------------------------------------------------------
# Result convenience methods
# ---------------------------------------------------------------------------


def test_result_to_dict_serialisable(tmp_path: Path):
    crystal = _write_crystal(tmp_path / "crystal.pdb")
    filled = _write_filled(tmp_path / "filled.pdb", ca2_offset=(0.1, 0, 0))
    r = bt.score_blind_fill(crystal, filled, [bt.GapRange("A", 2, 3)])
    import json as _json
    d = r.to_dict()
    s = _json.dumps(d)  # must not raise
    assert "per_range" in s
    assert "overall_ca_rmsd_A" in s


def test_summarise_produces_readable_line(tmp_path: Path):
    crystal = _write_crystal(tmp_path / "crystal.pdb")
    filled = _write_filled(tmp_path / "filled.pdb", ca2_offset=(0.1, 0, 0))
    r = bt.score_blind_fill(crystal, filled, [bt.GapRange("A", 2, 2)])
    s = bt.summarise(r)
    assert "blind-test" in s
    assert "Cα-RMSD" in s or "Å" in s


def test_summarise_when_skipped(tmp_path: Path):
    r = bt.score_blind_fill(tmp_path / "nope1", tmp_path / "nope2",
                            [bt.GapRange("A", 1, 5)])
    assert "skipped" in bt.summarise(r)


def test_residue_score_resname_match_flag():
    good = bt.ResidueBlindScore(
        chain="A", resnum=1,
        crystal_resname="ALA", filled_resname="ALA",
        ca_distance_A=0.0, backbone_rmsd_A=0.0,
    )
    bad = bt.ResidueBlindScore(
        chain="A", resnum=1,
        crystal_resname="ALA", filled_resname="GLY",
        ca_distance_A=0.0, backbone_rmsd_A=0.0,
    )
    missing = bt.ResidueBlindScore(
        chain="A", resnum=1,
        crystal_resname="ALA", filled_resname=None,
        ca_distance_A=None, backbone_rmsd_A=None,
    )
    assert good.resname_matches() is True
    assert bad.resname_matches() is False
    assert missing.resname_matches() is False


def test_range_score_backbone_rmsd_aggregate_matches_manual(tmp_path: Path):
    """Two residues, backbone RMSDs 0.25 and 0.0 → aggregate ≈ 0.177."""
    crystal = _write_crystal(tmp_path / "crystal.pdb")
    filled = _write_filled(tmp_path / "filled.pdb", ca2_offset=(0.5, 0, 0))
    r = bt.score_blind_fill(crystal, filled, [bt.GapRange("A", 2, 3)])
    rng = r.ranges[0]
    # sqrt((0.25^2 + 0^2)/2) = 0.25/sqrt(2) ≈ 0.1768
    assert rng.backbone_rmsd_A() == pytest.approx(0.25 / (2 ** 0.5), abs=1e-3)
