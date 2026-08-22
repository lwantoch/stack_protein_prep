"""Tests for stack_protein_preparation._filler_quality_check.

Exercises each independent metric on hand-constructed synthetic PDBs so a
regression in Rama-region rectangles, peptide-bond thresholds, chirality
sign, clash detection, or omega classification is caught by ci.
"""
from __future__ import annotations

from pathlib import Path

import pytest

from stack_protein_preparation._filler_quality_check import (
    QualityReport,
    _vdw_clash_threshold,
    check_model_quality,
)


def _write(pdb_path: Path, atoms: list[str]) -> Path:
    pdb_path.parent.mkdir(parents=True, exist_ok=True)
    pdb_path.write_text("\n".join(atoms) + "\nEND\n", encoding="utf-8")
    return pdb_path


def _atom(serial: int, name: str, resname: str, chain: str, resnum: int,
          x: float, y: float, z: float, element: str = "C") -> str:
    return (
        f"ATOM  {serial:5d} {name:>4} {resname:>3} {chain:1}"
        f"{resnum:4d}    "
        f"{x:8.3f}{y:8.3f}{z:8.3f}"
        f"  1.00  0.00           {element:>2}"
    )


# ----------------------------------------------------------------- vdW helper


def test_vdw_clash_threshold_known_pairs():
    # C(1.70) + N(1.55) - 0.5 = 2.75
    assert _vdw_clash_threshold("C", "N") == pytest.approx(2.75, abs=1e-6)
    # C-C = 2.90; O-O = 2.54
    assert _vdw_clash_threshold("C", "C") == pytest.approx(2.90, abs=1e-6)
    assert _vdw_clash_threshold("O", "O") == pytest.approx(2.54, abs=1e-6)
    # Unknown element falls back to 1.7 + 1.7 - 0.5 = 2.9 (safe over-estimate)
    assert _vdw_clash_threshold("X", "Y") == pytest.approx(2.90, abs=1e-6)


# ---------------------------------------------------------- empty / minimal


def test_empty_report_percentages_and_ratios():
    r = QualityReport()
    assert r.rama_favoured_pct() == 0.0
    assert r.rama_outlier_pct() == 0.0
    assert r.peptide_bond_tolerated_pct() == 0.0
    assert r.clashscore_per_1000_atoms() == 0.0


def test_check_model_quality_minimal_1_residue_pdb(tmp_path: Path) -> None:
    # Single ALA: no peptide bonds, no ω, no chirality (Gly-only-safe)
    p = _write(tmp_path / "ala1.pdb", [
        _atom(1, "N", "ALA", "A", 1, 0, 0, 0, "N"),
        _atom(2, "CA", "ALA", "A", 1, 1.46, 0, 0, "C"),
        _atom(3, "C", "ALA", "A", 1, 2.9, 0, 0, "C"),
        _atom(4, "O", "ALA", "A", 1, 3.5, 1.2, 0, "O"),
    ])
    r = check_model_quality(p)
    assert r.n_residues == 1
    assert r.n_peptide_bonds == 0
    assert r.n_ca_chirality_checked == 0  # missing CB
    assert r.n_omega_checked == 0
    assert r.n_vdw_clashes == 0


# ---------------------------------------------------------- peptide-bond classification


def test_peptide_bond_ideal_broken_and_tolerated_bins(tmp_path: Path) -> None:
    # 3-residue chain: bond 1 ideal (1.33), bond 2 broken (2.5 A).
    atoms = [
        _atom(1, "N", "ALA", "A", 1, 0, 0, 0, "N"),
        _atom(2, "CA", "ALA", "A", 1, 1.46, 0, 0, "C"),
        _atom(3, "C", "ALA", "A", 1, 2.9, 0, 0, "C"),
        _atom(4, "O", "ALA", "A", 1, 3.5, 1.2, 0, "O"),
        # res 2 with N at 1.33 A (ideal) from res 1 C
        _atom(5, "N", "GLY", "A", 2, 4.23, 0, 0, "N"),
        _atom(6, "CA", "GLY", "A", 2, 5.69, 0, 0, "C"),
        _atom(7, "C", "GLY", "A", 2, 7.13, 0, 0, "C"),
        _atom(8, "O", "GLY", "A", 2, 7.73, 1.2, 0, "O"),
        # res 3 with N at 2.5 A (broken) from res 2 C
        _atom(9, "N", "SER", "A", 3, 9.63, 0, 0, "N"),
        _atom(10, "CA", "SER", "A", 3, 11.09, 0, 0, "C"),
        _atom(11, "C", "SER", "A", 3, 12.53, 0, 0, "C"),
        _atom(12, "O", "SER", "A", 3, 13.13, 1.2, 0, "O"),
    ]
    p = _write(tmp_path / "peptide.pdb", atoms)
    r = check_model_quality(p)
    assert r.n_peptide_bonds == 2
    assert r.n_peptide_bonds_ideal == 1     # 1 -> 2 at 1.33 A
    assert r.n_peptide_bonds_broken == 1     # 2 -> 3 at 2.5 A
    # broken_peptide_bonds should have the offender labelled
    assert len(r.broken_peptide_bonds) == 1
    assert "GLY A:2→SER A:3" in r.broken_peptide_bonds[0]


def test_peptide_bond_skipped_across_resnum_gap(tmp_path: Path) -> None:
    # 2 residues numbered 1 and 5 (gap of 3) -- should NOT report a broken bond
    atoms = [
        _atom(1, "N", "ALA", "A", 1, 0, 0, 0, "N"),
        _atom(2, "CA", "ALA", "A", 1, 1.46, 0, 0, "C"),
        _atom(3, "C", "ALA", "A", 1, 2.9, 0, 0, "C"),
        _atom(4, "O", "ALA", "A", 1, 3.5, 1.2, 0, "O"),
        _atom(5, "N", "SER", "A", 5, 20, 0, 0, "N"),
        _atom(6, "CA", "SER", "A", 5, 21.46, 0, 0, "C"),
        _atom(7, "C", "SER", "A", 5, 22.9, 0, 0, "C"),
        _atom(8, "O", "SER", "A", 5, 23.5, 1.2, 0, "O"),
    ]
    p = _write(tmp_path / "gap.pdb", atoms)
    r = check_model_quality(p)
    assert r.n_peptide_bonds == 0  # numbering gap -> skipped
    assert r.n_peptide_bonds_broken == 0


# ---------------------------------------------------------- chirality


def test_chirality_l_amino_acid_passes(tmp_path: Path) -> None:
    # Standard L-Ala geometry: N, CA at origin, C along +x, CB below (-y, -z).
    # Signed volume ((N-CA) x (C-CA)) . (CB-CA) should be positive for L.
    atoms = [
        _atom(1, "N", "ALA", "A", 1, -0.500, 1.400, 0.000, "N"),
        _atom(2, "CA", "ALA", "A", 1, 0.000, 0.000, 0.000, "C"),
        _atom(3, "C", "ALA", "A", 1, 1.500, 0.000, 0.000, "C"),
        _atom(4, "CB", "ALA", "A", 1, -0.500, -0.800, -1.200, "C"),
        _atom(5, "O", "ALA", "A", 1, 2.000, 1.100, 0.000, "O"),
    ]
    p = _write(tmp_path / "l_ala.pdb", atoms)
    r = check_model_quality(p)
    assert r.n_ca_chirality_checked == 1
    assert r.n_ca_chirality_outliers == 0  # L configuration


def test_chirality_d_amino_acid_flagged(tmp_path: Path) -> None:
    # D-Ala: mirror CB across the CA plane -> sign flips.
    atoms = [
        _atom(1, "N", "ALA", "A", 1, -0.500, 1.400, 0.000, "N"),
        _atom(2, "CA", "ALA", "A", 1, 0.000, 0.000, 0.000, "C"),
        _atom(3, "C", "ALA", "A", 1, 1.500, 0.000, 0.000, "C"),
        _atom(4, "CB", "ALA", "A", 1, -0.500, -0.800, +1.200, "C"),  # flipped z
        _atom(5, "O", "ALA", "A", 1, 2.000, 1.100, 0.000, "O"),
    ]
    p = _write(tmp_path / "d_ala.pdb", atoms)
    r = check_model_quality(p)
    assert r.n_ca_chirality_checked == 1
    assert r.n_ca_chirality_outliers == 1  # D flagged
    assert "signed_vol" in r.ca_chirality_outlier_residues[0]


# ---------------------------------------------------------- omega classification


def test_omega_trans_natural_peptide(tmp_path: Path) -> None:
    # Build a trans peptide: CA1, C1, N2, CA2 all on same line and coplanar
    # -> omega ~ 180 deg (or -180).
    atoms = [
        _atom(1, "N", "ALA", "A", 1, -1.500, 0.000, 0.000, "N"),
        _atom(2, "CA", "ALA", "A", 1, 0.000, 0.000, 0.000, "C"),
        _atom(3, "C", "ALA", "A", 1, 1.500, 0.000, 0.000, "C"),
        _atom(4, "O", "ALA", "A", 1, 2.000, 1.100, 0.000, "O"),
        _atom(5, "N", "GLY", "A", 2, 2.830, 0.000, 0.000, "N"),
        # place CA2 along +x from N so omega = 180 deg
        _atom(6, "CA", "GLY", "A", 2, 4.290, 0.000, 0.000, "C"),
        _atom(7, "C", "GLY", "A", 2, 5.790, 0.000, 0.000, "C"),
        _atom(8, "O", "GLY", "A", 2, 6.390, 1.100, 0.000, "O"),
    ]
    p = _write(tmp_path / "trans.pdb", atoms)
    r = check_model_quality(p)
    assert r.n_omega_checked == 1
    assert r.n_omega_trans == 1
    assert r.n_omega_non_planar == 0
    assert r.n_omega_cis_nonpro == 0


# ---------------------------------------------------------- relative gate


def test_relative_gate_no_change_passes():
    b = QualityReport()
    b.n_residues = 100; b.n_rama_favoured = 95; b.n_rama_outlier = 1
    b.n_heavy_atoms = 800; b.n_clash_pairs = 2
    f = QualityReport()
    f.n_residues = 100; f.n_rama_favoured = 95; f.n_rama_outlier = 1
    f.n_heavy_atoms = 800; f.n_clash_pairs = 2
    passed, reasons = f.passes_relative_gate(b)
    assert passed
    assert reasons == []


def test_relative_gate_gained_broken_bond_fails():
    b = QualityReport(n_residues=100, n_peptide_bonds_broken=0)
    f = QualityReport(n_residues=100, n_peptide_bonds_broken=1)
    passed, reasons = f.passes_relative_gate(b)
    assert not passed
    assert any("Peptide bonds broken gained 1" in r for r in reasons)


def test_relative_gate_gained_chirality_d_fails():
    b = QualityReport(n_residues=100, n_ca_chirality_outliers=0)
    f = QualityReport(n_residues=100, n_ca_chirality_outliers=1)
    passed, reasons = f.passes_relative_gate(b)
    assert not passed
    assert any("Chirality D-outliers gained 1" in r for r in reasons)


def test_relative_gate_gained_cis_nonpro_fails():
    b = QualityReport(n_residues=100, n_omega_cis_nonpro=0)
    f = QualityReport(n_residues=100, n_omega_cis_nonpro=1)
    passed, reasons = f.passes_relative_gate(b)
    assert not passed
    assert any("cis-nonPro peptide bonds gained 1" in r for r in reasons)


def test_one_line_summary_contains_all_key_metrics():
    r = QualityReport()
    r.n_residues = 500; r.n_rama_favoured = 480; r.n_rama_outlier = 5
    r.n_heavy_atoms = 4000; r.n_vdw_clashes = 20
    r.n_peptide_bonds_broken = 0; r.n_omega_cis_nonpro = 0
    r.n_omega_non_planar = 0; r.n_ca_chirality_outliers = 0
    s = r.one_line_summary()
    assert "500 aa" in s
    assert "clashscore=5.0" in s
    assert "broken_bnd=0" in s
    assert "D-chir=0" in s
