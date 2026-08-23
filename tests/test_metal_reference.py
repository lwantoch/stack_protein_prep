"""Tests for stack_protein_preparation._metal_reference.

Verifies the 89-row TS metal reference CSV loads cleanly, all fields
parse into the expected types, canonical PDB resnames resolve to
sensible entries, and the validation-oracle correctly flags obvious
mismatches.
"""
from __future__ import annotations

import pytest

from stack_protein_preparation import _metal_reference


def test_table_loads_with_expected_row_count():
    resnames = _metal_reference.all_known_pdb_resnames()
    # 89 data rows in the CSV
    assert len(resnames) >= 80, f"Expected ~89 rows, got {len(resnames)}"


def test_zn_lookup_returns_zinc_entry():
    ref = _metal_reference.lookup_metal("ZN")
    assert ref is not None
    assert ref.element == "Zn"
    assert ref.oxidation_state == "+2"
    assert 4 in ref.coord_numbers
    assert "tetrahedral" in {g.lower() for g in ref.geometries}


def test_fe2_and_fe3_distinct():
    fe2 = _metal_reference.lookup_metal("FE2")
    fe3 = _metal_reference.lookup_metal("FE3")
    assert fe2 is not None and fe3 is not None
    assert fe2.oxidation_state == "+2"
    assert fe3.oxidation_state == "+3"


def test_cu1_soft_donor_preference():
    cu1 = _metal_reference.lookup_metal("CU1")
    assert cu1 is not None
    assert cu1.oxidation_state == "+1"
    assert "soft" in cu1.donor_preference_hsab.lower()


def test_lookup_case_insensitive():
    assert _metal_reference.lookup_metal("zn") == _metal_reference.lookup_metal("ZN")


def test_unknown_resname_returns_none():
    assert _metal_reference.lookup_metal("XXX") is None
    assert _metal_reference.lookup_metal("") is None


def test_validate_ok_for_canonical_zinc_site():
    ok, reasons = _metal_reference.validate_coordination(
        "ZN", coord_number=4, geometry="tetrahedral",
        donor_elements=("N", "N", "S", "S"),
    )
    assert ok is True, f"Canonical Zn²⁺ tetrahedral (2His+2Cys) should validate, got: {reasons}"


def test_validate_flags_wrong_coord_number():
    # Zn²⁺ with CN=8 is not a biological Zn coord environment
    ok, reasons = _metal_reference.validate_coordination(
        "ZN", coord_number=8, geometry="square_antiprismatic",
    )
    assert ok is False
    assert any("coord number" in r for r in reasons)


def test_validate_flags_wrong_geometry():
    # Zn²⁺ with CN=4 but linear geometry (Hg-like) is anomalous
    ok, reasons = _metal_reference.validate_coordination(
        "ZN", coord_number=4, geometry="linear",
    )
    assert ok is False
    assert any("geometry" in r for r in reasons)


def test_validate_unknown_resname_returns_advisory():
    ok, reasons = _metal_reference.validate_coordination(
        "UNKNOWN_METAL", coord_number=4,
    )
    assert ok is True  # fail-open when unknown
    assert any("no validation performed" in r for r in reasons)


def test_describe_returns_readable_summary():
    summary = _metal_reference.describe("ZN")
    assert "ZN" in summary
    assert "Zn+2" in summary
    assert "tetrahedral" in summary.lower() or "octahedral" in summary.lower()


def test_all_rows_have_element_and_charge():
    for resname in _metal_reference.all_known_pdb_resnames():
        ref = _metal_reference.lookup_metal(resname)
        assert ref is not None
        assert ref.element, f"row {resname} has no element"
        assert ref.oxidation_state or ref.common_charge, f"row {resname} has no charge"


def test_iron_sulfur_clusters_present():
    """SF4 / FES / F3S iron-sulfur cluster resnames should be in the table."""
    for cluster in ("SF4", "FES", "F3S"):
        ref = _metal_reference.lookup_metal(cluster)
        assert ref is not None, f"iron-sulfur cluster {cluster!r} missing from table"


def test_lanthanides_covered():
    """At least one Ln³⁺ ion (Gd) should be present -- MRI/probe use."""
    gd = _metal_reference.lookup_metal("GD")
    assert gd is not None
    assert "3" in gd.oxidation_state or "3" in gd.common_charge
