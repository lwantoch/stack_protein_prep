"""Tests for stack_protein_preparation._omega_scan."""
from __future__ import annotations

import math
from pathlib import Path

import pytest

from stack_protein_preparation import _omega_scan as os_


def _write_two_residues_trans(path: Path) -> Path:
    """Two residues in TRANS peptide bond (ω ≈ 180°)."""
    # Standard planar trans-ALA - ALA geometry.  Coordinates crafted so that
    # CA(1)-C(1)-N(2)-CA(2) is ~180° (all atoms coplanar).
    text = (
        "ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00  0.00           N\n"
        "ATOM      2  CA  ALA A   1       1.458   0.000   0.000  1.00  0.00           C\n"
        "ATOM      3  C   ALA A   1       2.009   1.420   0.000  1.00  0.00           C\n"
        "ATOM      4  O   ALA A   1       1.251   2.390   0.000  1.00  0.00           O\n"
        "ATOM      5  N   ALA A   2       3.332   1.540   0.000  1.00  0.00           N\n"
        "ATOM      6  CA  ALA A   2       4.000   2.828   0.000  1.00  0.00           C\n"
        "ATOM      7  C   ALA A   2       5.500   2.700   0.000  1.00  0.00           C\n"
        "ATOM      8  O   ALA A   2       6.000   1.600   0.000  1.00  0.00           O\n"
        "TER\nEND\n"
    )
    path.write_text(text)
    return path


def _write_two_residues_cis_pro(path: Path) -> Path:
    """Two residues where residue 2 is PRO and ω ≈ 0° (cis)."""
    text = (
        "ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00  0.00           N\n"
        "ATOM      2  CA  ALA A   1       1.458   0.000   0.000  1.00  0.00           C\n"
        "ATOM      3  C   ALA A   1       2.009   1.420   0.000  1.00  0.00           C\n"
        "ATOM      4  O   ALA A   1       1.251   2.390   0.000  1.00  0.00           O\n"
        # For cis, place N(2) and CA(2) on the SAME side as CA(1) w.r.t. C(1)-N(2)
        "ATOM      5  N   PRO A   2       3.332   1.540   0.000  1.00  0.00           N\n"
        "ATOM      6  CA  PRO A   2       3.900   0.200   0.000  1.00  0.00           C\n"
        "ATOM      7  C   PRO A   2       5.400   0.200   0.000  1.00  0.00           C\n"
        "ATOM      8  O   PRO A   2       6.000  -0.900   0.000  1.00  0.00           O\n"
        "TER\nEND\n"
    )
    path.write_text(text)
    return path


def _write_non_consecutive_pair(path: Path) -> Path:
    """Two residues but numbered 1 and 5 -> no ω computed."""
    text = (
        "ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00  0.00           N\n"
        "ATOM      2  CA  ALA A   1       1.458   0.000   0.000  1.00  0.00           C\n"
        "ATOM      3  C   ALA A   1       2.009   1.420   0.000  1.00  0.00           C\n"
        "ATOM      4  O   ALA A   1       1.251   2.390   0.000  1.00  0.00           O\n"
        "ATOM      5  N   ALA A   5       3.332   1.540   0.000  1.00  0.00           N\n"
        "ATOM      6  CA  ALA A   5       4.000   2.828   0.000  1.00  0.00           C\n"
        "ATOM      7  C   ALA A   5       5.500   2.700   0.000  1.00  0.00           C\n"
        "ATOM      8  O   ALA A   5       6.000   1.600   0.000  1.00  0.00           O\n"
        "TER\nEND\n"
    )
    path.write_text(text)
    return path


# ---------------------------------------------------------------------------
# Classification helper
# ---------------------------------------------------------------------------


def test_classify_trans():
    assert os_._classify(180.0, "ALA") == "trans"
    assert os_._classify(-179.0, "ALA") == "trans"
    assert os_._classify(160.0, "ALA") == "trans"


def test_classify_cis_pro():
    assert os_._classify(0.0, "PRO") == "cis_pro"
    assert os_._classify(15.0, "PRO") == "cis_pro"
    assert os_._classify(-25.0, "PRO") == "cis_pro"


def test_classify_cis_nonpro():
    assert os_._classify(0.0, "ALA") == "cis_nonpro"
    assert os_._classify(20.0, "GLY") == "cis_nonpro"


def test_classify_non_planar():
    assert os_._classify(60.0, "ALA") == "non_planar"
    assert os_._classify(-90.0, "ALA") == "non_planar"
    assert os_._classify(140.0, "ALA") == "non_planar"


def test_classify_boundary_150():
    """|ω| = 150 exact -> not trans (must be > 150)."""
    assert os_._classify(150.0, "ALA") == "non_planar"
    assert os_._classify(150.001, "ALA") == "trans"


def test_classify_boundary_30():
    """|ω| = 30 exact -> not cis (must be < 30)."""
    assert os_._classify(30.0, "PRO") == "non_planar"
    assert os_._classify(29.999, "PRO") == "cis_pro"


# ---------------------------------------------------------------------------
# scan_omega_dihedrals I/O
# ---------------------------------------------------------------------------


def test_missing_pdb_fail_open(tmp_path: Path):
    r = os_.scan_omega_dihedrals(tmp_path / "nope.pdb")
    assert r.ran is False
    assert r.passed is False
    assert "not found" in r.fallback_reason


def test_bad_pdb_returns_ran_but_not_passed(tmp_path: Path):
    """Parseable but empty-of-residues PDB -> passed=False."""
    (tmp_path / "empty.pdb").write_text("HEADER\nEND\n")
    r = os_.scan_omega_dihedrals(tmp_path / "empty.pdb")
    assert r.ran is True
    assert r.passed is False
    assert r.entries == []


def test_trans_pair_produces_one_entry(tmp_path: Path):
    pdb = _write_two_residues_trans(tmp_path / "trans.pdb")
    r = os_.scan_omega_dihedrals(pdb)
    assert r.ran is True
    assert r.passed is True
    assert len(r.entries) == 1
    e = r.entries[0]
    assert e.chain == "A"
    assert e.resnum_i == 1 and e.resnum_j == 2
    assert e.resname_i == "ALA" and e.resname_j == "ALA"
    assert e.kind == "trans"
    assert abs(abs(e.omega_deg) - 180.0) < 5.0


def test_non_consecutive_pair_skipped(tmp_path: Path):
    pdb = _write_non_consecutive_pair(tmp_path / "gap.pdb")
    r = os_.scan_omega_dihedrals(pdb)
    assert r.ran is True
    # residues numbered 1 and 5 -> skipped
    assert r.entries == []


# ---------------------------------------------------------------------------
# Result convenience methods
# ---------------------------------------------------------------------------


def _fake_result_with_kinds(kinds: list[str]) -> os_.OmegaScanResult:
    entries = [
        os_.OmegaEntry(
            chain="A", resnum_i=i, resname_i="ALA",
            resnum_j=i + 1, resname_j="ALA",
            omega_deg=180.0 if k == "trans" else 0.0,
            kind=k,
        )
        for i, k in enumerate(kinds, start=1)
    ]
    return os_.OmegaScanResult(ran=True, passed=bool(entries), entries=entries)


def test_result_counts():
    r = _fake_result_with_kinds(
        ["trans"] * 100 + ["cis_pro"] * 3 + ["cis_nonpro"] * 1 + ["non_planar"] * 2
    )
    assert r.n_trans() == 100
    assert r.n_cis_pro() == 3
    assert r.n_cis_nonpro() == 1
    assert r.n_non_planar() == 2


def test_result_to_dict():
    r = _fake_result_with_kinds(["trans", "cis_pro"])
    d = r.to_dict()
    assert d["n_total"] == 2
    assert d["n_trans"] == 1
    assert d["n_cis_pro"] == 1
    assert d["n_cis_nonpro"] == 0
    assert d["n_non_planar"] == 0


def test_omega_values_extracts_series():
    r = _fake_result_with_kinds(["trans", "trans", "cis_nonpro"])
    vals = r.omega_values()
    assert len(vals) == 3
    assert all(isinstance(v, float) for v in vals)


def test_summarise_string_format():
    r = _fake_result_with_kinds(["trans"] * 199 + ["cis_pro"] * 1)
    s = os_.summarise(r)
    assert "200 pairs" in s
    assert "trans 199" in s
    assert "cis-Pro 1" in s


def test_summarise_empty_scan():
    r = os_.OmegaScanResult(ran=True, passed=False)
    assert "no residue pairs" in os_.summarise(r)


def test_aggregate_omega_values_skips_failed():
    good = _fake_result_with_kinds(["trans", "trans"])
    bad = os_.OmegaScanResult(ran=False, passed=False)
    vals = os_.aggregate_omega_values([good, bad, good])
    assert len(vals) == 4  # 2 + 0 + 2
