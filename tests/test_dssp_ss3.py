"""Tests for stack_protein_preparation._dssp_ss3.

Never launches real DSSP -- shutil.which is monkey-patched to force
the fail-open path, and cross_check_idr_vs_local_ss is tested with
synthetic DSSPResult objects.
"""
from __future__ import annotations

from pathlib import Path

import pytest

from stack_protein_preparation import _dssp_ss3 as ds


def _fake_dssp_result(ss3_map: dict[tuple[str, int], str]) -> ds.DSSPResult:
    return ds.DSSPResult(
        ran=True, passed=True, ss3_by_residue=dict(ss3_map),
        diagnostics=["synthetic"],
    )


# ---------------------------------------------------------------------------
# fail-open paths
# ---------------------------------------------------------------------------


def test_missing_pdb_returns_fail_open(tmp_path: Path):
    r = ds.compute_dssp_ss3(tmp_path / "nope.pdb")
    assert r.ran is False
    assert "not found" in r.fallback_reason


def test_missing_dssp_binary_returns_fail_open(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch,
):
    # Force shutil.which to return None for every name
    import shutil as _sh
    monkeypatch.setattr(_sh, "which", lambda name: None)
    fake_pdb = tmp_path / "prot.pdb"
    fake_pdb.write_text("ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00  0.00           N\nEND\n")
    r = ds.compute_dssp_ss3(fake_pdb)
    assert r.ran is False
    assert "dssp" in r.fallback_reason.lower()


# ---------------------------------------------------------------------------
# DSSPResult convenience methods
# ---------------------------------------------------------------------------


def test_dssp_result_counts_ss3():
    r = _fake_dssp_result({
        ("A", 1): "H", ("A", 2): "H", ("A", 3): "E", ("A", 4): "C",
        ("A", 5): "C", ("A", 6): "H",
    })
    assert r.n_helix() == 3
    assert r.n_sheet() == 1
    assert r.n_coil() == 2


def test_dssp_result_to_dict_contains_metrics():
    r = _fake_dssp_result({("A", 1): "H", ("A", 2): "C"})
    d = r.to_dict()
    assert d["n_residues"] == 2
    assert d["n_helix"] == 1
    assert d["n_coil"] == 1
    assert d["ran"] is True


def test_dssp_8to3_mapping():
    """H/G/I -> H, E/B -> E, everything else -> C."""
    assert ds._DSSP8_TO_SS3["H"] == "H"
    assert ds._DSSP8_TO_SS3["G"] == "H"
    assert ds._DSSP8_TO_SS3["I"] == "H"
    assert ds._DSSP8_TO_SS3["E"] == "E"
    assert ds._DSSP8_TO_SS3["B"] == "E"
    assert ds._DSSP8_TO_SS3["T"] == "C"
    assert ds._DSSP8_TO_SS3["S"] == "C"
    assert ds._DSSP8_TO_SS3["-"] == "C"


def test_summarise_ss3_composition():
    r = _fake_dssp_result({
        ("A", i): "H" for i in range(1, 88)  # 87 H
    } | {
        ("A", i): "E" for i in range(88, 133)  # 45 E
    } | {
        ("A", i): "C" for i in range(133, 213)  # 80 C
    })
    s = ds.summarise_ss3_composition(r)
    assert "212 residues" in s
    assert "H 87" in s
    assert "E 45" in s
    assert "C 80" in s


def test_summarise_empty_dssp():
    r = ds.DSSPResult(ran=True, passed=False, ss3_by_residue={})
    assert "no residues assigned" in ds.summarise_ss3_composition(r)


# ---------------------------------------------------------------------------
# cross_check_idr_vs_local_ss
# ---------------------------------------------------------------------------


def test_cross_check_empty_when_dssp_not_ran():
    r = ds.DSSPResult(ran=False, passed=False)
    disagreements = ds.cross_check_idr_vs_local_ss(
        r,
        gap_ranges_by_chain=[("A", 50, 60)],
        idr_flag_by_gap={("A", 50, 60): True},
    )
    assert disagreements == []


def test_cross_check_override_when_flank_is_helix():
    """Flanking crystal is helix -> advise override the IDR reject."""
    r = _fake_dssp_result({
        ("A", 47): "H", ("A", 48): "H", ("A", 49): "H",  # left flank
        ("A", 61): "H", ("A", 62): "H", ("A", 63): "H",  # right flank
    })
    disagreements = ds.cross_check_idr_vs_local_ss(
        r,
        gap_ranges_by_chain=[("A", 50, 60)],
        idr_flag_by_gap={("A", 50, 60): True},
        flank_window=3,
    )
    assert len(disagreements) >= 1
    for d in disagreements:
        assert d.ss3 == "H"
        assert "OVERRIDING" in d.advice or "override" in d.advice.lower()


def test_cross_check_keep_reject_when_flank_is_coil():
    """Flanking crystal is all coil -> IDR reject stands."""
    r = _fake_dssp_result({
        ("A", 47): "C", ("A", 48): "C", ("A", 49): "C",
        ("A", 61): "C", ("A", 62): "C", ("A", 63): "C",
    })
    disagreements = ds.cross_check_idr_vs_local_ss(
        r,
        gap_ranges_by_chain=[("A", 50, 60)],
        idr_flag_by_gap={("A", 50, 60): True},
    )
    assert len(disagreements) == 1
    assert disagreements[0].ss3 == "C"
    assert "stands" in disagreements[0].advice.lower()


def test_cross_check_skips_non_idr_gaps():
    """Non-IDR gap should never appear in the disagreement list."""
    r = _fake_dssp_result({("A", 47): "H", ("A", 61): "H"})
    disagreements = ds.cross_check_idr_vs_local_ss(
        r,
        gap_ranges_by_chain=[("A", 50, 60)],
        idr_flag_by_gap={("A", 50, 60): False},  # NOT IDR
    )
    assert disagreements == []


def test_cross_check_mixed_flank_finds_structured_residues():
    """One helix + several coils in flank -> emit advisory for the H."""
    r = _fake_dssp_result({
        ("A", 47): "C", ("A", 48): "H", ("A", 49): "C",  # 48 is H
        ("A", 61): "C", ("A", 62): "C", ("A", 63): "C",
    })
    disagreements = ds.cross_check_idr_vs_local_ss(
        r,
        gap_ranges_by_chain=[("A", 50, 60)],
        idr_flag_by_gap={("A", 50, 60): True},
    )
    # Exactly one structured residue (H at 48) -> one override advisory
    assert len(disagreements) == 1
    assert disagreements[0].ss3 == "H"
    assert disagreements[0].resnum == 48


def test_cross_check_flank_window_argument():
    """Small window misses far-away H; large window catches it."""
    r = _fake_dssp_result({("A", 55): "H"})  # inside gap boundary
    # Actually this is INSIDE the gap; the flank check needs to see OUTSIDE
    # Fix: put H outside the gap
    r = _fake_dssp_result({("A", 45): "H"})  # 5 residues before gap start
    small = ds.cross_check_idr_vs_local_ss(
        r,
        gap_ranges_by_chain=[("A", 50, 60)],
        idr_flag_by_gap={("A", 50, 60): True},
        flank_window=3,
    )
    large = ds.cross_check_idr_vs_local_ss(
        r,
        gap_ranges_by_chain=[("A", 50, 60)],
        idr_flag_by_gap={("A", 50, 60): True},
        flank_window=10,
    )
    # With window=3, residue 45 is out of reach (50-3=47).
    # With window=10, residue 45 is in reach (50-10=40).
    small_has_h = any(d.ss3 == "H" for d in small)
    large_has_h = any(d.ss3 == "H" for d in large)
    assert small_has_h is False
    assert large_has_h is True
