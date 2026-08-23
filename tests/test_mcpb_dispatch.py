"""Tests for stack_protein_preparation._mcpb_dispatch."""
from __future__ import annotations

import pytest

from stack_protein_preparation import _mcpb_dispatch as dispatch
from stack_protein_preparation._component_confidence import Confidence
from stack_protein_preparation._metal_reference import MetalReference


def _fake_ref(
    resname="ZN", element="Zn", oxid=2, coord=(4,),
    route="12-6-4", notes="Zn2+ Li-Merz",
) -> MetalReference:
    return MetalReference(
        pdb_resname=resname, element=element, oxidation_state=oxid,
        common_charge=oxid, coord_numbers=coord, geometries=("tetrahedral",),
        donor_preference_hsab="borderline", typical_distances_angstrom={"N": (2.0, 2.2), "O": (2.0, 2.3)},
        spin_state_default="diamagnetic", force_field_route=route,
        notes=notes,
    )


# ---------------------------------------------------------------------------
# Env-var override
# ---------------------------------------------------------------------------


def test_env_override_gaussian_forces_dft(monkeypatch):
    monkeypatch.setenv("FRUTON_MCPB_ENGINE", "gaussian")
    d = dispatch.decide_mcpb_route("ZN", "no keywords", metal_reference=_fake_ref())
    assert d.tier == dispatch.TIER_DFT_BONDED
    assert d.confidence is Confidence.HIGH
    assert d.method_label == "dft_gaussian"


def test_env_override_orca_forces_dft(monkeypatch):
    monkeypatch.setenv("FRUTON_MCPB_ENGINE", "orca")
    d = dispatch.decide_mcpb_route("ZN", "", metal_reference=_fake_ref())
    assert d.tier == dispatch.TIER_DFT_BONDED
    assert d.method_label == "dft_orca"


def test_env_override_xtb_forces_xtb(monkeypatch):
    monkeypatch.setenv("FRUTON_MCPB_ENGINE", "xtb")
    d = dispatch.decide_mcpb_route("ZN", "", metal_reference=_fake_ref())
    assert d.tier == dispatch.TIER_XTB_BONDED
    assert d.confidence is Confidence.MEDIUM


def test_env_override_unknown_value_ignored(monkeypatch):
    monkeypatch.setenv("FRUTON_MCPB_ENGINE", "nonsense")
    d = dispatch.decide_mcpb_route("ZN", "", metal_reference=_fake_ref())
    # Falls through to default path — Li-Merz for ZN
    assert d.tier == dispatch.TIER_NONBONDED


# ---------------------------------------------------------------------------
# Nonbonded (Li-Merz) route
# ---------------------------------------------------------------------------


def test_zn_labile_no_paper_flag_goes_nonbonded(monkeypatch):
    monkeypatch.delenv("FRUTON_MCPB_ENGINE", raising=False)
    d = dispatch.decide_mcpb_route(
        "ZN", "Zn coordinates 4 waters in structural role",
        metal_reference=_fake_ref(route="12-6-4"),
    )
    assert d.tier == dispatch.TIER_NONBONDED
    assert d.confidence is Confidence.HIGH
    assert d.method_label == "LiMerz_12_6_4"


def test_zn_no_route_no_paper_flag_default_nonbonded(monkeypatch):
    monkeypatch.delenv("FRUTON_MCPB_ENGINE", raising=False)
    d = dispatch.decide_mcpb_route(
        "ZN", "", metal_reference=_fake_ref(route=""),
    )
    assert d.tier == dispatch.TIER_NONBONDED
    # No explicit route → MEDIUM confidence conservative default
    assert d.confidence is Confidence.MEDIUM
    assert d.method_label == "LiMerz_12_6_4_default"


# ---------------------------------------------------------------------------
# xtb bonded route
# ---------------------------------------------------------------------------


def test_zn_catalytic_paper_flag_upgrades_to_xtb(monkeypatch):
    monkeypatch.delenv("FRUTON_MCPB_ENGINE", raising=False)
    d = dispatch.decide_mcpb_route(
        "ZN",
        "H94 is the catalytic residue coordinating the Zn2+ ion",
        metal_reference=_fake_ref(route="12-6-4"),
    )
    assert d.tier == dispatch.TIER_XTB_BONDED
    assert d.confidence is Confidence.MEDIUM
    assert d.method_label == "xtb_gfn2"


def test_reference_requires_mcpb_route_selects_xtb(monkeypatch):
    monkeypatch.delenv("FRUTON_MCPB_ENGINE", raising=False)
    d = dispatch.decide_mcpb_route(
        "CU", "structural site",
        metal_reference=_fake_ref(resname="CU", element="Cu", route="mcpb-required"),
    )
    assert d.tier == dispatch.TIER_XTB_BONDED


# ---------------------------------------------------------------------------
# DFT-only route
# ---------------------------------------------------------------------------


def test_iron_sulfur_cluster_forces_dft(monkeypatch):
    monkeypatch.delenv("FRUTON_MCPB_ENGINE", raising=False)
    d = dispatch.decide_mcpb_route(
        "FE",
        "the [4Fe-4S] cluster mediates electron transfer",
        metal_reference=_fake_ref(resname="FE", element="Fe", route="mcpb-required"),
    )
    assert d.tier == dispatch.TIER_DFT_BONDED
    assert "spin-crossover" in d.reason.lower() or "fe-s" in d.reason.lower() or "radical" in d.reason.lower()
    assert d.confidence is Confidence.MEDIUM
    assert "orca" in d.suggested_action.lower()


def test_spin_crossover_forces_dft(monkeypatch):
    monkeypatch.delenv("FRUTON_MCPB_ENGINE", raising=False)
    d = dispatch.decide_mcpb_route(
        "FE",
        "displays high-spin/low-spin transition on ligand binding",
        metal_reference=_fake_ref(resname="FE", element="Fe"),
    )
    assert d.tier == dispatch.TIER_DFT_BONDED


def test_mixed_valence_forces_dft(monkeypatch):
    monkeypatch.delenv("FRUTON_MCPB_ENGINE", raising=False)
    d = dispatch.decide_mcpb_route(
        "CU", "the CuA site is a mixed-valence dicopper centre",
        metal_reference=_fake_ref(resname="CU", element="Cu"),
    )
    assert d.tier == dispatch.TIER_DFT_BONDED


# ---------------------------------------------------------------------------
# Manual route (unknown metal)
# ---------------------------------------------------------------------------


def test_unknown_metal_goes_manual(monkeypatch):
    monkeypatch.delenv("FRUTON_MCPB_ENGINE", raising=False)
    d = dispatch.decide_mcpb_route("XYZ", "", metal_reference=None)
    assert d.tier == dispatch.TIER_MANUAL
    assert d.confidence is Confidence.LOW
    assert "not in" in d.reason
    assert "add" in d.suggested_action.lower() or "supply" in d.suggested_action.lower()


# ---------------------------------------------------------------------------
# Bridge to ComponentConfidence
# ---------------------------------------------------------------------------


def test_to_component_confidence_populates_all_fields(monkeypatch):
    monkeypatch.delenv("FRUTON_MCPB_ENGINE", raising=False)
    d = dispatch.decide_mcpb_route(
        "ZN", "catalytic Zn2+", metal_reference=_fake_ref(),
    )
    cc = d.to_component_confidence(name="Zn A501")
    assert cc.component_type == "metal"
    assert cc.name == "Zn A501"
    assert cc.method == d.method_label
    assert cc.confidence == d.confidence
    assert cc.reason == d.reason
    assert cc.details["tier"] == d.tier
    assert cc.details["element"] == d.element


def test_decisions_always_have_reasons():
    """Every code path must produce a non-empty reason for reviewer defensibility."""
    scenarios = [
        ("ZN", "", _fake_ref(route="12-6-4")),
        ("ZN", "catalytic", _fake_ref(route="12-6-4")),
        ("FE", "[4Fe-4S]", _fake_ref(resname="FE", element="Fe", route="mcpb-required")),
        ("XYZ", "", None),
    ]
    for resname, evidence, ref in scenarios:
        d = dispatch.decide_mcpb_route(resname, evidence, metal_reference=ref)
        assert d.reason and len(d.reason) > 10, f"skinny reason for {resname}: {d.reason!r}"
