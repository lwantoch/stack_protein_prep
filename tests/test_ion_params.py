"""Tests for stack_protein_preparation._ion_params.

Loads the 34-row Li-Merz 12-6-4 reference CSV and verifies the
public lookup API (per-resname params, coverage boolean, tleap
source-line dedup, MCPB flag, describe).
"""
from __future__ import annotations

import pytest

from stack_protein_preparation import _ion_params


def test_table_loads_with_expected_row_count():
    ions = _ion_params.all_known_ions()
    assert len(ions) >= 30, f"expected ~34 rows, got {len(ions)}"


def test_zn_lookup_returns_zinc_entry():
    ip = _ion_params.lookup_ion("ZN")
    assert ip is not None
    assert ip.element == "Zn"
    assert ip.formal_charge == "+2"
    assert "12-6-4" in ip.leaprc_route


def test_fe_default_is_fe3():
    """PDB resname 'FE' defaults to Fe(III) per common practice."""
    ip = _ion_params.lookup_ion("FE")
    assert ip is not None
    assert "+3" in ip.formal_charge


def test_fe2_and_fe3_distinct_from_fe():
    ip2 = _ion_params.lookup_ion("FE2")
    ip3 = _ion_params.lookup_ion("FE3")
    assert ip2 is not None and ip3 is not None
    assert ip2.formal_charge == "+2"
    assert ip3.formal_charge == "+3"


def test_lookup_case_insensitive():
    assert _ion_params.lookup_ion("zn") == _ion_params.lookup_ion("ZN")


def test_lookup_unknown_returns_none():
    assert _ion_params.lookup_ion("XXX") is None
    assert _ion_params.lookup_ion("") is None
    assert _ion_params.lookup_ion(None) is None  # type: ignore[arg-type]


def test_covered_by_12_6_4():
    assert _ion_params.covered_by_12_6_4("ZN") is True
    assert _ion_params.covered_by_12_6_4("GD") is True
    assert _ion_params.covered_by_12_6_4("XXX") is False


def test_tleap_source_lines_dedup():
    """Multiple ions in the same leaprc should produce one source line."""
    lines = _ion_params.tleap_source_lines(["ZN", "MG", "CA", "FE2"])
    joined = "\n".join(lines)
    # all covered by the same base 12-6-4 leaprc -> exactly one entry
    assert joined.count("source leaprc.water.tip3p_ion12-6-4") == 1


def test_tleap_source_lines_mixed_base_and_lanthanide():
    """GD requires the HFE-tuned lanthanide extension in addition."""
    lines = _ion_params.tleap_source_lines(["ZN", "GD"])
    assert any("12-6-4" in ln and "HFE" not in ln for ln in lines)
    assert any("HFE" in ln for ln in lines)


def test_tleap_source_lines_empty_on_unknown_only():
    assert _ion_params.tleap_source_lines(["XXX", "YYY"]) == []


def test_tleap_source_lines_preserves_input_order_for_new_leaprcs():
    """First unique leaprc appears first in output."""
    zn_first = _ion_params.tleap_source_lines(["ZN", "GD"])
    gd_first = _ion_params.tleap_source_lines(["GD", "ZN"])
    # Both should contain both source lines, but order differs
    assert set(zn_first) == set(gd_first)
    assert zn_first[0] != gd_first[0]


def test_flag_mcpb_required_catches_zn_type1_cu():
    """Catalytic Zn / type-1 Cu need MCPB.py bonded modeling."""
    flagged = _ion_params.flag_mcpb_required(["ZN", "CU", "MG", "NA"])
    # Zn note says "catalytic Zn sites usually need MCPB.py bonded model"
    assert "ZN" in flagged
    # Cu note says "type-1 sites need MCPB"
    assert "CU" in flagged
    # Mg and Na don't need MCPB
    assert "MG" not in flagged
    assert "NA" not in flagged


def test_describe_readable_summary():
    s = _ion_params.describe("ZN")
    assert "ZN" in s
    assert "Zn+2" in s
    assert "12-6-4" in s


def test_lanthanides_have_hfe_leaprc():
    for ln in ("CE", "GD", "EU", "YB"):
        ip = _ion_params.lookup_ion(ln)
        assert ip is not None
        assert "HFE" in ip.leaprc_route


def test_alkali_ions_present():
    """Standard water leaprc covers alkali; still worth being listed for
    the tleap-source dedup path."""
    for a in ("NA", "K", "LI"):
        ip = _ion_params.lookup_ion(a)
        assert ip is not None
        assert ip.formal_charge == "+1"


def test_needs_mcpb_property():
    zn = _ion_params.lookup_ion("ZN")
    assert zn is not None
    assert zn.needs_mcpb is True
    na = _ion_params.lookup_ion("NA")
    assert na is not None
    assert na.needs_mcpb is False
