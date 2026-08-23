"""Tests for stack_protein_preparation._metal_coord_scan."""
from __future__ import annotations

from pathlib import Path

import pytest

from stack_protein_preparation import _metal_coord_scan as mcs


# ---------------------------------------------------------------------------
# PDB fixtures
# ---------------------------------------------------------------------------


def _write_zn_his_pdb(path: Path, zn_distance: float = 2.10) -> Path:
    """Minimal PDB: HIS with NE2 + a Zn ion at ``zn_distance`` Å from NE2."""
    # NE2 at origin; place Zn on +x at zn_distance
    path.write_text(
        "ATOM      1  N   HIS A   1       0.000   0.000   3.000  1.00  0.00           N\n"
        "ATOM      2  CA  HIS A   1       0.000   0.000   4.000  1.00  0.00           C\n"
        "ATOM      3  C   HIS A   1       0.000   0.000   5.000  1.00  0.00           C\n"
        "ATOM      4  O   HIS A   1       1.000   0.000   5.000  1.00  0.00           O\n"
        "ATOM      5  CB  HIS A   1       0.000   1.000   4.000  1.00  0.00           C\n"
        "ATOM      6  ND1 HIS A   1       0.000   2.000   4.000  1.00  0.00           N\n"
        "ATOM      7  NE2 HIS A   1       0.000   0.000   0.000  1.00  0.00           N\n"
        f"HETATM    8 ZN    ZN B 101    {zn_distance: >8.3f}   0.000   0.000  1.00  0.00          ZN\n"
        "END\n"
    )
    return path


def _write_two_metals_pdb(path: Path) -> Path:
    """Two-metal case: Zn coordinated by CYS-SG, Mg coordinated by ASP-OD1."""
    path.write_text(
        "ATOM      1  N   CYS A   1       0.000   0.000   3.000  1.00  0.00           N\n"
        "ATOM      2  CA  CYS A   1       0.000   0.000   4.000  1.00  0.00           C\n"
        "ATOM      3  C   CYS A   1       0.000   0.000   5.000  1.00  0.00           C\n"
        "ATOM      4  O   CYS A   1       1.000   0.000   5.000  1.00  0.00           O\n"
        "ATOM      5  CB  CYS A   1       0.000   1.000   4.000  1.00  0.00           C\n"
        "ATOM      6  SG  CYS A   1       0.000   0.000   0.000  1.00  0.00           S\n"
        "ATOM      7  N   ASP A   2      10.000   0.000   3.000  1.00  0.00           N\n"
        "ATOM      8  CA  ASP A   2      10.000   0.000   4.000  1.00  0.00           C\n"
        "ATOM      9  C   ASP A   2      10.000   0.000   5.000  1.00  0.00           C\n"
        "ATOM     10  O   ASP A   2      11.000   0.000   5.000  1.00  0.00           O\n"
        "ATOM     11  CB  ASP A   2      10.000   1.000   4.000  1.00  0.00           C\n"
        "ATOM     12  CG  ASP A   2      10.000   2.000   4.000  1.00  0.00           C\n"
        "ATOM     13  OD1 ASP A   2      10.000   0.000   0.000  1.00  0.00           O\n"
        "ATOM     14  OD2 ASP A   2      10.000   3.000   5.000  1.00  0.00           O\n"
        "HETATM   15 ZN    ZN B 101       2.300   0.000   0.000  1.00  0.00          ZN\n"
        "HETATM   16 MG    MG B 102      12.000   0.000   0.000  1.00  0.00          MG\n"
        "END\n"
    )
    return path


def _write_no_metal_pdb(path: Path) -> Path:
    path.write_text(
        "ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00  0.00           N\n"
        "ATOM      2  CA  ALA A   1       1.458   0.000   0.000  1.00  0.00           C\n"
        "ATOM      3  C   ALA A   1       2.009   1.420   0.000  1.00  0.00           C\n"
        "END\n"
    )
    return path


# ---------------------------------------------------------------------------
# constants
# ---------------------------------------------------------------------------


def test_donor_atoms_include_std_histidines():
    for name in ("HIS", "HID", "HIE", "HIP"):
        assert mcs.DONOR_ATOMS_BY_RESNAME[name] == {"ND1", "NE2"}


def test_metal_elements_include_common_transition_metals():
    for m in ("ZN", "MG", "MN", "FE", "CU", "CA"):
        assert m in mcs.METAL_ELEMENTS


# ---------------------------------------------------------------------------
# scan_metal_coordination — fail-open
# ---------------------------------------------------------------------------


def test_missing_pdb_returns_fail_open(tmp_path: Path):
    r = mcs.scan_metal_coordination(tmp_path / "nope.pdb")
    assert r.ran is False
    assert "not found" in r.fallback_reason


def test_no_metal_returns_zero_contacts(tmp_path: Path):
    p = _write_no_metal_pdb(tmp_path / "no_metal.pdb")
    r = mcs.scan_metal_coordination(p)
    assert r.ran is True and r.passed is True
    assert r.n_metals == 0
    assert r.contacts == []


# ---------------------------------------------------------------------------
# scan_metal_coordination — happy path
# ---------------------------------------------------------------------------


def test_zn_his_single_contact(tmp_path: Path):
    p = _write_zn_his_pdb(tmp_path / "zn_his.pdb", zn_distance=2.10)
    r = mcs.scan_metal_coordination(p)
    assert r.ran is True
    assert r.n_metals == 1
    assert r.n_contacts() == 1
    c = r.contacts[0]
    assert c.key.metal_element == "ZN"
    assert c.key.donor_resname == "HIS"
    assert c.key.donor_atom == "NE2"
    assert abs(c.distance_A - 2.10) < 0.05


def test_zn_beyond_cutoff_no_contact(tmp_path: Path):
    p = _write_zn_his_pdb(tmp_path / "zn_far.pdb", zn_distance=5.00)
    r = mcs.scan_metal_coordination(p, coord_cutoff_A=3.0)
    assert r.ran is True
    assert r.n_metals == 1
    assert r.n_contacts() == 0


def test_two_metals_two_contacts(tmp_path: Path):
    p = _write_two_metals_pdb(tmp_path / "two.pdb")
    r = mcs.scan_metal_coordination(p, coord_cutoff_A=3.0)
    assert r.n_metals == 2
    # Zn↔SG (2.3 Å) and Mg↔OD1 (2.0 Å) both within cutoff
    assert r.n_contacts() == 2
    by_element = r.by_metal_element()
    assert by_element == {"ZN": 1, "MG": 1}


def test_result_to_dict(tmp_path: Path):
    p = _write_zn_his_pdb(tmp_path / "zn.pdb")
    r = mcs.scan_metal_coordination(p)
    d = r.to_dict()
    assert d["n_metals"] == 1
    assert d["n_contacts"] == 1
    assert d["by_metal_element"]["ZN"] == 1


def test_summarise_message_for_metals(tmp_path: Path):
    p = _write_two_metals_pdb(tmp_path / "two.pdb")
    r = mcs.scan_metal_coordination(p)
    s = mcs.summarise(r)
    assert "2 metals" in s
    assert "ZN 1" in s and "MG 1" in s


def test_summarise_no_metal(tmp_path: Path):
    p = _write_no_metal_pdb(tmp_path / "empty.pdb")
    r = mcs.scan_metal_coordination(p)
    assert "no metal" in mcs.summarise(r).lower()


# ---------------------------------------------------------------------------
# compare_metal_coordination
# ---------------------------------------------------------------------------


def _fake_contact(distance: float, donor_atom: str = "NE2") -> mcs.MetalDonorContact:
    return mcs.MetalDonorContact(
        key=mcs.DonorKey(
            metal_chain="B", metal_resnum=101, metal_element="ZN",
            donor_chain="A", donor_resnum=1, donor_resname="HIS",
            donor_atom=donor_atom,
        ),
        distance_A=distance,
    )


def _fake_result(contacts: list[mcs.MetalDonorContact]) -> mcs.MetalCoordScanResult:
    return mcs.MetalCoordScanResult(
        ran=True, passed=True, n_metals=1, contacts=contacts,
    )


def test_compare_matches_same_donor():
    crystal = _fake_result([_fake_contact(2.05)])
    fruton = _fake_result([_fake_contact(2.31)])
    deltas = mcs.compare_metal_coordination(crystal, fruton)
    assert len(deltas) == 1
    assert deltas[0].status() == "preserved"
    assert abs(deltas[0].delta_A() - 0.26) < 1e-6


def test_compare_flags_lost_donor():
    crystal = _fake_result([_fake_contact(2.05, donor_atom="NE2")])
    fruton = _fake_result([])  # donor displaced beyond cutoff
    deltas = mcs.compare_metal_coordination(crystal, fruton)
    assert deltas[0].status() == "lost"
    assert deltas[0].delta_A() is None


def test_compare_flags_gained_donor():
    crystal = _fake_result([])
    fruton = _fake_result([_fake_contact(2.20)])
    deltas = mcs.compare_metal_coordination(crystal, fruton)
    assert deltas[0].status() == "gained"


def test_compare_preserves_multiple_donors_sorted():
    crystal = _fake_result([
        _fake_contact(2.05, donor_atom="NE2"),
        _fake_contact(2.15, donor_atom="ND1"),
    ])
    fruton = _fake_result([
        _fake_contact(2.30, donor_atom="NE2"),
        _fake_contact(2.05, donor_atom="ND1"),
    ])
    deltas = mcs.compare_metal_coordination(crystal, fruton)
    assert [d.key.donor_atom for d in deltas] == ["ND1", "NE2"]  # sorted alpha


# ---------------------------------------------------------------------------
# aggregate_deltas
# ---------------------------------------------------------------------------


def test_aggregate_skips_lost_and_gained():
    deltas_pdb1 = [
        mcs.MetalDonorDelta(
            key=mcs.DonorKey("B", 101, "ZN", "A", 1, "HIS", "NE2"),
            distance_crystal_A=2.05, distance_fruton_A=2.10,
        ),
        mcs.MetalDonorDelta(
            key=mcs.DonorKey("B", 101, "ZN", "A", 2, "HIS", "ND1"),
            distance_crystal_A=2.20, distance_fruton_A=None,  # lost
        ),
    ]
    deltas_pdb2 = [
        mcs.MetalDonorDelta(
            key=mcs.DonorKey("B", 102, "MG", "A", 5, "ASP", "OD1"),
            distance_crystal_A=None, distance_fruton_A=2.00,  # gained
        ),
    ]
    values = mcs.aggregate_deltas([deltas_pdb1, deltas_pdb2])
    assert len(values) == 1  # only the one preserved donor
    assert abs(values[0] - 0.05) < 1e-6
