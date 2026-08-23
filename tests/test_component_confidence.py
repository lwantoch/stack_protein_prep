"""Tests for stack_protein_preparation._component_confidence."""
from __future__ import annotations

import json

import pytest

from stack_protein_preparation._component_confidence import (
    Confidence,
    ComponentConfidence,
    ProteinDeliveryReport,
    summarise_delivery,
)


# ---------------------------------------------------------------------------
# Confidence enum
# ---------------------------------------------------------------------------


def test_confidence_worst_of_returns_max():
    assert Confidence.worst_of([Confidence.HIGH, Confidence.MEDIUM]) is Confidence.MEDIUM
    assert Confidence.worst_of([Confidence.HIGH, Confidence.LOW]) is Confidence.LOW
    assert Confidence.worst_of([Confidence.HIGH, Confidence.FAILED]) is Confidence.FAILED
    assert Confidence.worst_of([Confidence.HIGH, Confidence.HIGH]) is Confidence.HIGH


def test_confidence_worst_of_empty_defaults_high():
    """No components → we haven't attempted anything, but by definition we
    haven't broken anything either, so HIGH."""
    assert Confidence.worst_of([]) is Confidence.HIGH


# ---------------------------------------------------------------------------
# ComponentConfidence validation
# ---------------------------------------------------------------------------


def test_component_requires_reason():
    with pytest.raises(ValueError, match="reason must not be empty"):
        ComponentConfidence(
            component_type="metal", name="Zn A501",
            confidence=Confidence.HIGH, reason="",
        )


def test_component_rejects_unknown_type():
    with pytest.raises(ValueError, match="unknown component_type"):
        ComponentConfidence(
            component_type="not_a_component", name="x",
            confidence=Confidence.HIGH, reason="test",
        )


def test_component_accepts_string_confidence():
    """JSON round-trip robustness: 'high' string becomes Confidence.HIGH."""
    c = ComponentConfidence(
        component_type="metal", name="Zn A501",
        confidence="high", reason="Li-Merz",
    )
    assert c.confidence is Confidence.HIGH


def test_component_to_dict_shape():
    c = ComponentConfidence(
        component_type="gap_fill", name="chain A gap 45-58",
        confidence=Confidence.MEDIUM,
        reason="mean pLDDT 65, gap length 14",
        suggested_action="visually inspect the 14-residue loop",
        method="MODELLER_slow",
        details={"n_conformers_kept": 3, "mean_plddt": 65.2},
    )
    d = c.to_dict()
    assert d["confidence"] == "medium"
    assert d["method"] == "MODELLER_slow"
    assert d["details"]["mean_plddt"] == 65.2


# ---------------------------------------------------------------------------
# ProteinDeliveryReport overall_status
# ---------------------------------------------------------------------------


def _c(conf: Confidence, comp_type: str = "gap_fill", name: str = "x") -> ComponentConfidence:
    return ComponentConfidence(
        component_type=comp_type, name=name,
        confidence=conf, reason="test-only",
        suggested_action="review" if conf is not Confidence.HIGH else "",
    )


def test_status_full_confidence_when_all_high():
    r = ProteinDeliveryReport(
        pdb="8ABC",
        components=[_c(Confidence.HIGH), _c(Confidence.HIGH)],
    )
    assert r.overall_status == "delivered_full_confidence"


def test_status_with_notes_when_any_medium_only():
    r = ProteinDeliveryReport(
        pdb="8ABC",
        components=[_c(Confidence.HIGH), _c(Confidence.MEDIUM)],
    )
    assert r.overall_status == "delivered_with_notes"


def test_status_needs_review_when_any_low():
    r = ProteinDeliveryReport(
        pdb="8ABC",
        components=[_c(Confidence.HIGH), _c(Confidence.MEDIUM), _c(Confidence.LOW)],
    )
    assert r.overall_status == "delivered_needs_review"


def test_status_needs_review_when_component_failed_but_model_written():
    """USER REFRAME: a failed component doesn't mean the whole model is
    not_delivered — as long as the pipeline still produced a partial
    model, it's delivered_needs_review."""
    r = ProteinDeliveryReport(
        pdb="8ABC", model_written=True,
        components=[_c(Confidence.HIGH), _c(Confidence.FAILED, "metal", "Fe A700")],
    )
    assert r.overall_status == "delivered_needs_review"


def test_status_not_delivered_only_when_pipeline_crashed():
    r = ProteinDeliveryReport(
        pdb="8ABC", model_written=False,
        components=[_c(Confidence.HIGH)],
    )
    assert r.overall_status == "not_delivered"


# ---------------------------------------------------------------------------
# action items + component_type_counts
# ---------------------------------------------------------------------------


def test_action_items_only_include_non_high():
    r = ProteinDeliveryReport(
        pdb="8ABC",
        components=[
            _c(Confidence.HIGH, "gap_fill", "gap 45-50"),
            _c(Confidence.LOW, "metal", "Fe A700"),
            _c(Confidence.MEDIUM, "cofactor", "HEM A600"),
        ],
    )
    items = r.action_items()
    assert len(items) == 2
    assert any("Fe A700" in it for it in items)
    assert any("HEM A600" in it for it in items)


def test_component_type_counts_populated():
    r = ProteinDeliveryReport(
        pdb="8ABC",
        components=[
            _c(Confidence.HIGH, "gap_fill", "g1"),
            _c(Confidence.HIGH, "gap_fill", "g2"),
            _c(Confidence.MEDIUM, "gap_fill", "g3"),
            _c(Confidence.HIGH, "metal", "Zn A501"),
        ],
    )
    counts = r.component_type_counts()
    assert counts["gap_fill"]["high"] == 2
    assert counts["gap_fill"]["medium"] == 1
    assert counts["metal"]["high"] == 1


# ---------------------------------------------------------------------------
# to_dict + JSON round-trip
# ---------------------------------------------------------------------------


def test_to_dict_is_json_serialisable():
    r = ProteinDeliveryReport(
        pdb="8ABC",
        components=[
            _c(Confidence.HIGH, "metal", "Zn A501"),
            _c(Confidence.LOW, "gap_fill", "g1"),
        ],
        tleap_loads=True, md_deck_written=True,
    )
    s = json.dumps(r.to_dict())
    reloaded = json.loads(s)
    assert reloaded["pdb"] == "8ABC"
    assert reloaded["overall_status"] == "delivered_needs_review"
    assert len(reloaded["components"]) == 2


# ---------------------------------------------------------------------------
# summarise_delivery (bench-wide)
# ---------------------------------------------------------------------------


def test_summarise_delivery_counts_by_status():
    reports = [
        ProteinDeliveryReport(pdb="A", components=[_c(Confidence.HIGH)]),
        ProteinDeliveryReport(pdb="B", components=[_c(Confidence.MEDIUM)]),
        ProteinDeliveryReport(pdb="C", components=[_c(Confidence.LOW)]),
        ProteinDeliveryReport(pdb="D", model_written=False),
    ]
    s = summarise_delivery(reports)
    assert s["n_proteins"] == 4
    assert s["by_overall_status"]["delivered_full_confidence"] == 1
    assert s["by_overall_status"]["delivered_with_notes"] == 1
    assert s["by_overall_status"]["delivered_needs_review"] == 1
    assert s["by_overall_status"]["not_delivered"] == 1


def test_summarise_delivery_component_totals():
    reports = [
        ProteinDeliveryReport(
            pdb="A",
            components=[
                _c(Confidence.HIGH, "metal", "Zn A501"),
                _c(Confidence.HIGH, "metal", "Zn A502"),
                _c(Confidence.MEDIUM, "cofactor", "HEM A600"),
                _c(Confidence.LOW, "gap_fill", "g1"),
            ],
        ),
    ]
    s = summarise_delivery(reports)
    assert s["component_totals"]["high"] == 2
    assert s["component_totals"]["medium"] == 1
    assert s["component_totals"]["low"] == 1
    assert s["component_totals"]["total"] == 4
