"""Tests for stack_protein_preparation._adaptive_refine_policy.

Encodes the reviewer-defensibility invariants: no hand-picked
per-benchmark thresholds, all decisions derived from either the
protein's own crystal baseline or an INDEPENDENT reference population.
"""
from __future__ import annotations

from dataclasses import dataclass

import pytest

from stack_protein_preparation import _adaptive_refine_policy as arp


@dataclass
class _FakeQC:
    """Just the fields the policy inspects."""
    n_clash_pairs: int = 0
    n_omega_non_planar: int = 0
    n_omega_cis_nonpro: int = 0
    n_peptide_bonds_broken: int = 0


# Test-only reference (small, hand-crafted, NOT the 199-crystal packaged one).
_TEST_REFERENCE = {
    "n_clash_pairs":              {"p50": 1, "p90":  7, "p95": 12, "p99": 30, "max": 33},
    "n_omega_non_planar":         {"p50": 0, "p90":  1, "p95":  1, "p99":  2, "max":  3},
    "n_omega_cis_nonpro":         {"p50": 0, "p90":  1, "p95":  2, "p99":  3, "max":  5},
    "n_peptide_bonds_broken":     {"p50": 0, "p90":  1, "p95":  2, "p99":  7, "max": 10},
}


# ---------------------------------------------------------------------------
# packaged reference data
# ---------------------------------------------------------------------------


def test_packaged_reference_loads():
    ref = arp.load_reference_percentiles()
    # sanity: has all four metrics + provenance
    for metric in arp._REFERENCE_METRICS:
        assert metric in ref
        assert "p99" in ref[metric]
    assert "_provenance" in ref


# ---------------------------------------------------------------------------
# no-regression path (fruton not worse than crystal)
# ---------------------------------------------------------------------------


def test_no_regression_no_retry():
    b = _FakeQC(n_clash_pairs=3, n_omega_non_planar=1)
    f = _FakeQC(n_clash_pairs=3, n_omega_non_planar=1)  # identical
    d = arp.decide_refine_action(b, f, reference_percentiles=_TEST_REFERENCE)
    assert d.action == "no_retry"
    assert d.regressed_metrics == {}
    assert d.exceeded_ceilings == {}


def test_fruton_improves_no_retry():
    """FRUTON output BETTER than crystal → no retry."""
    b = _FakeQC(n_clash_pairs=10, n_omega_non_planar=2)
    f = _FakeQC(n_clash_pairs=3, n_omega_non_planar=0)
    d = arp.decide_refine_action(b, f, reference_percentiles=_TEST_REFERENCE)
    assert d.action == "no_retry"


# ---------------------------------------------------------------------------
# regression → retry_slow (Principle 1)
# ---------------------------------------------------------------------------


def test_any_omega_gain_triggers_retry_slow():
    """+1 ω non-planar vs crystal → retry_slow, no magic threshold."""
    b = _FakeQC(n_omega_non_planar=0)
    f = _FakeQC(n_omega_non_planar=1)
    d = arp.decide_refine_action(b, f, reference_percentiles=_TEST_REFERENCE)
    assert d.action == "retry_slow"
    assert "n_omega_non_planar" in d.regressed_metrics
    assert d.regressed_metrics["n_omega_non_planar"] == (0, 1)


def test_any_clash_gain_triggers_retry_slow():
    b = _FakeQC(n_clash_pairs=3)
    f = _FakeQC(n_clash_pairs=5)
    d = arp.decide_refine_action(b, f, reference_percentiles=_TEST_REFERENCE)
    assert d.action == "retry_slow"


def test_broken_bond_gain_triggers_retry_slow():
    b = _FakeQC(n_peptide_bonds_broken=0)
    f = _FakeQC(n_peptide_bonds_broken=1)
    d = arp.decide_refine_action(b, f, reference_percentiles=_TEST_REFERENCE)
    assert d.action == "retry_slow"


def test_cis_nonpro_gain_triggers_retry_slow():
    b = _FakeQC(n_omega_cis_nonpro=0)
    f = _FakeQC(n_omega_cis_nonpro=1)
    d = arp.decide_refine_action(b, f, reference_percentiles=_TEST_REFERENCE)
    assert d.action == "retry_slow"


def test_multiple_regressions_all_reported():
    b = _FakeQC(n_clash_pairs=3, n_omega_non_planar=0)
    f = _FakeQC(n_clash_pairs=5, n_omega_non_planar=1)
    d = arp.decide_refine_action(b, f, reference_percentiles=_TEST_REFERENCE)
    assert d.action == "retry_slow"
    assert set(d.regressed_metrics) == {"n_clash_pairs", "n_omega_non_planar"}


# ---------------------------------------------------------------------------
# ceiling → skip_ceiling (Principle 2)
# ---------------------------------------------------------------------------


def test_clash_above_p99_reference_skips_ceiling():
    """4AT5 (clash_gain 582) — fruton clashes way beyond p99=30."""
    b = _FakeQC(n_clash_pairs=3)
    f = _FakeQC(n_clash_pairs=585)
    d = arp.decide_refine_action(b, f, reference_percentiles=_TEST_REFERENCE)
    assert d.action == "skip_ceiling"
    assert "n_clash_pairs" in d.exceeded_ceilings
    fruton_val, ref_val = d.exceeded_ceilings["n_clash_pairs"]
    assert fruton_val == 585 and ref_val == 30


def test_clash_at_reference_p99_boundary_does_not_skip():
    """Boundary check: exactly at p99 is NOT beyond it — retry, not skip."""
    b = _FakeQC(n_clash_pairs=3)
    f = _FakeQC(n_clash_pairs=30)  # exactly p99
    d = arp.decide_refine_action(b, f, reference_percentiles=_TEST_REFERENCE)
    assert d.action == "retry_slow"


def test_ceiling_takes_precedence_over_regression():
    """When both regression AND ceiling breached, ceiling wins (skip)."""
    b = _FakeQC(n_clash_pairs=1, n_omega_non_planar=0)
    f = _FakeQC(n_clash_pairs=100, n_omega_non_planar=1)
    d = arp.decide_refine_action(b, f, reference_percentiles=_TEST_REFERENCE)
    assert d.action == "skip_ceiling"


def test_configurable_ceiling_percentile_p90():
    """Users can choose a stricter ceiling (p90 instead of p99)."""
    b = _FakeQC(n_clash_pairs=1)
    f = _FakeQC(n_clash_pairs=15)  # above p90=7 but below p99=30
    d_p99 = arp.decide_refine_action(b, f, reference_percentiles=_TEST_REFERENCE, ceiling_percentile="p99")
    d_p90 = arp.decide_refine_action(b, f, reference_percentiles=_TEST_REFERENCE, ceiling_percentile="p90")
    assert d_p99.action == "retry_slow"
    assert d_p90.action == "skip_ceiling"


# ---------------------------------------------------------------------------
# reviewer-defensibility invariants
# ---------------------------------------------------------------------------


def test_thresholds_come_from_reference_data_not_hardcoded():
    """Policy body must not embed magic numbers matched to bench-observed values.

    Strips docstrings/comments before checking so historical mentions in prose
    don't false-positive.
    """
    import ast, inspect
    tree = ast.parse(inspect.getsource(arp))
    # Drop docstrings from module + every function/class
    for node in ast.walk(tree):
        if isinstance(node, (ast.Module, ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef)):
            if (node.body
                    and isinstance(node.body[0], ast.Expr)
                    and isinstance(node.body[0].value, ast.Constant)
                    and isinstance(node.body[0].value.value, str)):
                node.body.pop(0)
    code_only = ast.unparse(tree)
    # Strip inline comments
    code_only = "\n".join(
        (line.split("#", 1)[0] if not line.strip().startswith("#") else "")
        for line in code_only.splitlines()
    )
    # These are the case-specific magic thresholds we removed.
    for magic in ("> 200", ">= 3", ">= 1", "> 20 or"):
        assert magic not in code_only, (
            f"magic threshold {magic!r} snuck back into policy code"
        )


def test_regressed_scenarios_cover_all_gate_metrics():
    """Sanity: every metric in _REFERENCE_METRICS is checked by the policy."""
    for metric in arp._REFERENCE_METRICS:
        b = _FakeQC()
        f = _FakeQC(**{metric: 1})  # +1 on just this one
        d = arp.decide_refine_action(b, f, reference_percentiles=_TEST_REFERENCE)
        assert metric in d.regressed_metrics, f"policy failed to detect {metric} regression"


def test_decision_to_dict_round_trip():
    b = _FakeQC(n_clash_pairs=1)
    f = _FakeQC(n_clash_pairs=5)
    d = arp.decide_refine_action(b, f, reference_percentiles=_TEST_REFERENCE)
    payload = d.to_dict()
    assert payload["action"] == "retry_slow"
    assert payload["regressed_metrics"]["n_clash_pairs"] == [1, 5]


def test_missing_metric_field_defaults_to_zero():
    """Robustness: partial QC objects (e.g. legacy formats) don't crash."""
    class _Partial:
        n_clash_pairs = 5  # only one field
    b = _Partial()
    f = _FakeQC(n_clash_pairs=8)
    d = arp.decide_refine_action(b, f, reference_percentiles=_TEST_REFERENCE)
    assert d.action == "retry_slow"
