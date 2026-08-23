"""Tests for _metal_coord_protonation.

Exhaustive coverage per the user mandate 'genau musst du jede Möglichkeit
betrachten' — one test per amino acid protonation case + bidentate/mono/
bridging edge cases + REMARK 620 interaction + PDB rewrite.
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import pytest

from stack_protein_preparation import _metal_coord_protonation as mp
from stack_protein_preparation._component_confidence import Confidence


# --- fixture: synthetic MetalCoordScanResult ---------------------------------

@dataclass
class _FakeKey:
    metal_chain: str; metal_resnum: int; metal_element: str
    donor_chain: str; donor_resnum: int; donor_resname: str; donor_atom: str


@dataclass
class _FakeContact:
    key: _FakeKey; distance_A: float


@dataclass
class _FakeScan:
    ran: bool = True
    passed: bool = True
    contacts: list = None


def _c(donor_resname, donor_atom, dist, resnum=100, chain="A",
       metal="ZN", m_chain="B", m_resnum=500):
    return _FakeContact(
        _FakeKey(m_chain, m_resnum, metal, chain, resnum, donor_resname, donor_atom),
        dist,
    )


# ---------------------------------------------------------------------------
# HIS: REMARK 620 wins when present, else geometry
# ---------------------------------------------------------------------------


def test_his_remark620_takes_precedence():
    scan = _FakeScan(contacts=[_c("HIS", "NE2", 2.05), _c("HIS", "ND1", 3.5)])
    r620 = {("A", 100, ""): "HID"}
    ovs = mp.derive_overrides_from_geometry("/nonexistent.pdb", r620, "", scan)
    assert len(ovs) == 1
    assert ovs[0].forced_resname == "HID"
    assert ovs[0].confidence is Confidence.HIGH
    assert ovs[0].provenance == "remark620"


def test_his_geometry_fallback_ne2_closer_gives_hid():
    scan = _FakeScan(contacts=[_c("HIS", "NE2", 2.1), _c("HIS", "ND1", 3.5)])
    ovs = mp.derive_overrides_from_geometry("/nonexistent.pdb", {}, "", scan)
    assert ovs[0].forced_resname == "HID"      # NE2→metal → HID
    assert ovs[0].donor_atom == "NE2"
    assert ovs[0].confidence is Confidence.MEDIUM
    assert ovs[0].provenance == "geometry"


def test_his_geometry_fallback_nd1_closer_gives_hie():
    scan = _FakeScan(contacts=[_c("HIS", "NE2", 3.5), _c("HIS", "ND1", 2.15)])
    ovs = mp.derive_overrides_from_geometry("/nonexistent.pdb", {}, "", scan)
    assert ovs[0].forced_resname == "HIE"      # ND1→metal → HIE
    assert ovs[0].donor_atom == "ND1"


def test_his_donor_beyond_cutoff_skipped():
    scan = _FakeScan(contacts=[_c("HIS", "NE2", 3.6), _c("HIS", "ND1", 3.7)])
    ovs = mp.derive_overrides_from_geometry("/nonexistent.pdb", {}, "", scan)
    assert ovs == []


# ---------------------------------------------------------------------------
# CYS: always CYM when metal within cutoff
# ---------------------------------------------------------------------------


def test_cys_sg_close_becomes_cym():
    scan = _FakeScan(contacts=[_c("CYS", "SG", 2.3)])
    ovs = mp.derive_overrides_from_geometry("/nonexistent.pdb", {}, "", scan)
    assert ovs[0].forced_resname == "CYM"
    assert ovs[0].confidence is Confidence.HIGH


def test_cyx_metal_contact_medium_confidence_with_warning():
    """Disulfide Cys with metal contact → unusual, verify."""
    scan = _FakeScan(contacts=[_c("CYX", "SG", 2.3)])
    ovs = mp.derive_overrides_from_geometry("/nonexistent.pdb", {}, "", scan)
    assert ovs[0].forced_resname == "CYM"
    assert ovs[0].confidence is Confidence.MEDIUM
    assert "verify" in ovs[0].suggested_action.lower()


def test_cys_sg_beyond_cutoff_skipped():
    scan = _FakeScan(contacts=[_c("CYS", "SG", 3.4)])   # S cutoff is 3.2
    ovs = mp.derive_overrides_from_geometry("/nonexistent.pdb", {}, "", scan)
    assert ovs == []


# ---------------------------------------------------------------------------
# ASP / GLU: bidentate vs monodentate; paper evidence for ASH/GLH
# ---------------------------------------------------------------------------


def test_asp_bidentate_high_confidence():
    scan = _FakeScan(contacts=[_c("ASP", "OD1", 2.4), _c("ASP", "OD2", 2.5)])
    ovs = mp.derive_overrides_from_geometry("/nonexistent.pdb", {}, "", scan)
    assert ovs[0].forced_resname == "ASP"
    assert ovs[0].bidentate is True
    assert ovs[0].confidence is Confidence.HIGH


def test_asp_monodentate_no_paper_hint_default_deprot_medium():
    scan = _FakeScan(contacts=[_c("ASP", "OD1", 2.4)])
    ovs = mp.derive_overrides_from_geometry("/nonexistent.pdb", {}, "", scan)
    assert ovs[0].forced_resname == "ASP"
    assert ovs[0].bidentate is False
    assert ovs[0].confidence is Confidence.MEDIUM


def test_asp_monodentate_with_paper_ash_hint_becomes_ash_low_conf():
    scan = _FakeScan(contacts=[_c("ASP", "OD1", 2.5)])
    ovs = mp.derive_overrides_from_geometry(
        "/nonexistent.pdb", {}, "the buried Asp45 is in the protonated ASH form", scan,
    )
    assert ovs[0].forced_resname == "ASH"
    assert ovs[0].confidence is Confidence.LOW


def test_glu_bidentate_high_confidence():
    scan = _FakeScan(contacts=[_c("GLU", "OE1", 2.4), _c("GLU", "OE2", 2.4)])
    ovs = mp.derive_overrides_from_geometry("/nonexistent.pdb", {}, "", scan)
    assert ovs[0].forced_resname == "GLU"
    assert ovs[0].bidentate is True
    assert ovs[0].confidence is Confidence.HIGH


def test_glu_monodentate_paper_glh_becomes_glh():
    scan = _FakeScan(contacts=[_c("GLU", "OE1", 2.5)])
    ovs = mp.derive_overrides_from_geometry(
        "/nonexistent.pdb", {}, "the neutral carboxylate form of Glu at pH...", scan,
    )
    assert ovs[0].forced_resname == "GLH"


# ---------------------------------------------------------------------------
# TYR, LYS, SEC — the exotic-but-real donors
# ---------------------------------------------------------------------------


def test_tyr_oh_close_becomes_tym():
    scan = _FakeScan(contacts=[_c("TYR", "OH", 2.2, metal="FE")])
    ovs = mp.derive_overrides_from_geometry("/nonexistent.pdb", {}, "", scan)
    assert ovs[0].forced_resname == "TYM"
    assert ovs[0].confidence is Confidence.MEDIUM
    assert "tyrosinate" in ovs[0].reason.lower()


def test_lys_nz_close_becomes_lyn():
    scan = _FakeScan(contacts=[_c("LYS", "NZ", 2.4, metal="NI")])
    ovs = mp.derive_overrides_from_geometry("/nonexistent.pdb", {}, "", scan)
    assert ovs[0].forced_resname == "LYN"
    assert "neutral amine" in ovs[0].reason.lower()


def test_sec_selenolate_low_confidence():
    scan = _FakeScan(contacts=[_c("SEC", "SE", 2.6, metal="MO")])
    ovs = mp.derive_overrides_from_geometry("/nonexistent.pdb", {}, "", scan)
    assert ovs[0].forced_resname == "SEC"
    assert ovs[0].confidence is Confidence.LOW
    assert "selenolate" in ovs[0].reason.lower()


# ---------------------------------------------------------------------------
# Residues that DO NOT need overrides
# ---------------------------------------------------------------------------


def test_met_no_override_generated():
    """Met-Sδ coordinates but is already neutral thioether — no override."""
    scan = _FakeScan(contacts=[_c("MET", "SD", 2.3, metal="CU")])
    ovs = mp.derive_overrides_from_geometry("/nonexistent.pdb", {}, "", scan)
    assert ovs == []


def test_asn_no_override_generated():
    scan = _FakeScan(contacts=[_c("ASN", "OD1", 2.5)])
    ovs = mp.derive_overrides_from_geometry("/nonexistent.pdb", {}, "", scan)
    assert ovs == []


# ---------------------------------------------------------------------------
# Bridging metals (same donor coordinates 2 metals)
# ---------------------------------------------------------------------------


def test_asp_bridging_two_metals_marked_bridging():
    scan = _FakeScan(contacts=[
        _c("ASP", "OD1", 2.4, metal="ZN", m_resnum=500),
        _c("ASP", "OD2", 2.5, metal="ZN", m_resnum=501),
    ])
    ovs = mp.derive_overrides_from_geometry("/nonexistent.pdb", {}, "", scan)
    assert ovs[0].bridging_metals >= 1  # counted per-atom occurrences


def test_multiple_residues_each_get_override():
    scan = _FakeScan(contacts=[
        _c("HIS", "NE2", 2.1, resnum=94),
        _c("HIS", "ND1", 2.15, resnum=96),
        _c("CYS", "SG", 2.3, resnum=124),
    ])
    ovs = mp.derive_overrides_from_geometry("/nonexistent.pdb", {}, "", scan)
    assert len(ovs) == 3
    forced = {o.forced_resname for o in ovs}
    assert forced == {"HID", "HIE", "CYM"}


# ---------------------------------------------------------------------------
# ComponentConfidence bridge
# ---------------------------------------------------------------------------


def test_to_component_confidence_shape():
    scan = _FakeScan(contacts=[_c("HIS", "NE2", 2.05, resnum=94)])
    ov = mp.derive_overrides_from_geometry("/nonexistent.pdb", {}, "", scan)[0]
    cc = ov.to_component_confidence()
    assert cc.component_type == "protonation"
    assert "HIS94A" in cc.name and "HID" in cc.name
    assert cc.confidence is Confidence.MEDIUM
    assert cc.details["donor_atom"] == "NE2"
    assert cc.details["distance_A"] == 2.05


# ---------------------------------------------------------------------------
# PDB rewrite
# ---------------------------------------------------------------------------


def test_apply_overrides_rewrites_resname(tmp_path: Path):
    pdb = tmp_path / "in.pdb"
    pdb.write_text(
        "ATOM      1  N   HIS A  94       0.000   0.000   0.000  1.00  0.00           N\n"
        "ATOM      2  CA  HIS A  94       1.458   0.000   0.000  1.00  0.00           C\n"
        "ATOM      3  N   CYS A 124       5.000   0.000   0.000  1.00  0.00           N\n"
        "END\n"
    )
    scan = _FakeScan(contacts=[
        _c("HIS", "NE2", 2.05, resnum=94),
        _c("CYS", "SG", 2.3, resnum=124),
    ])
    ovs = mp.derive_overrides_from_geometry(pdb, {}, "", scan)
    out = mp.apply_overrides_to_pdb(pdb, tmp_path / "out.pdb", ovs)
    text = out.read_text()
    assert "HID A  94" in text
    assert "CYM A 124" in text
    assert "HIS A  94" not in text
    assert "CYS A 124" not in text


def test_apply_overrides_preserves_non_target_residues(tmp_path: Path):
    """Every ATOM line without a matching override must be unchanged."""
    pdb = tmp_path / "in.pdb"
    original_atom = "ATOM      1  N   ALA A  50       0.000   0.000   0.000  1.00  0.00           N"
    pdb.write_text(original_atom + "\nEND\n")
    ovs = []  # empty
    out = mp.apply_overrides_to_pdb(pdb, tmp_path / "out.pdb", ovs)
    assert original_atom in out.read_text()


# ---------------------------------------------------------------------------
# Summary
# ---------------------------------------------------------------------------


def test_summarise_empty():
    assert "no metal-coord" in mp.summarise_overrides([])


def test_summarise_counts_by_type():
    scan = _FakeScan(contacts=[
        _c("HIS", "NE2", 2.05, resnum=94),
        _c("HIS", "NE2", 2.05, resnum=96),
        _c("CYS", "SG", 2.3, resnum=124),
    ])
    ovs = mp.derive_overrides_from_geometry("/nonexistent.pdb", {}, "", scan)
    s = mp.summarise_overrides(ovs)
    assert "HID×2" in s or "HID" in s
    assert "CYM" in s
    assert "MEDIUM=" in s or "HIGH=" in s
