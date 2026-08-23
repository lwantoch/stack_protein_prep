"""Tests for stack_protein_preparation._audit_report."""
from __future__ import annotations

import csv
import json
from pathlib import Path

import pytest

from stack_protein_preparation import _audit_report as ar


_TEST_REFERENCE = {
    "n_clash_pairs":              {"p50": 1, "p90":  7, "p95": 12, "p99": 30, "max": 33},
    "n_omega_non_planar":         {"p50": 0, "p90":  1, "p95":  1, "p99":  2, "max":  3},
    "n_omega_cis_nonpro":         {"p50": 0, "p90":  1, "p95":  2, "p99":  3, "max":  5},
    "n_peptide_bonds_broken":     {"p50": 0, "p90":  1, "p95":  2, "p99":  7, "max": 10},
}


def _bench_row(pdb, **kw):
    return {
        "pdb": pdb,
        "gate_pass": kw.get("gate_pass", True),
        "n_gaps": kw.get("n_gaps", 1),
        "delta_n": kw.get("delta_n", 0),
        "rolled_gaps": kw.get("rolled_gaps", 0),
        "n_clash_pairs": kw.get("n_clash_pairs", 0),
        "n_omega_non_planar": kw.get("n_omega_non_planar", 0),
        "n_omega_cis_nonpro": kw.get("n_omega_cis_nonpro", 0),
        "n_peptide_bonds_broken": kw.get("n_peptide_bonds_broken", 0),
        "reasons": kw.get("reasons", []),
        "non_planar_omega_residues": kw.get("non_planar_omega_residues", []),
    }


# ---------------------------------------------------------------------------
# tier classification (no bench-specific magic — all bands from reference)
# ---------------------------------------------------------------------------


def test_volltreffer_when_all_metrics_in_p90_and_no_regression():
    r = _bench_row("8ABC", n_clash_pairs=3, n_omega_non_planar=0)
    row = ar.classify_protein(
        r,
        crystal_metrics={"n_clash_pairs": 3, "n_omega_non_planar": 0,
                         "n_omega_cis_nonpro": 0, "n_peptide_bonds_broken": 0},
        reference_percentiles=_TEST_REFERENCE,
    )
    assert row.tier == "volltreffer"
    assert row.action == "accept"


def test_wahrscheinlich_ok_when_metric_in_p90_p95_band():
    r = _bench_row("8ABC", n_clash_pairs=10)  # p90=7, p95=12 → p90-p95
    row = ar.classify_protein(r, reference_percentiles=_TEST_REFERENCE)
    assert row.tier == "wahrscheinlich_ok"


def test_wahrscheinlich_ok_on_small_regression_only():
    """+1 clash vs crystal but everything still in p90 → wahrscheinlich_ok."""
    r = _bench_row("8ABC", n_clash_pairs=4)
    row = ar.classify_protein(
        r, crystal_metrics={"n_clash_pairs": 3, "n_omega_non_planar": 0,
                            "n_omega_cis_nonpro": 0, "n_peptide_bonds_broken": 0},
        reference_percentiles=_TEST_REFERENCE,
    )
    assert row.tier == "wahrscheinlich_ok"


def test_grenzwertig_when_metric_in_p95_p99_band():
    r = _bench_row("8ABC", n_clash_pairs=20)  # p95=12 p99=30 → p95-p99
    row = ar.classify_protein(r, reference_percentiles=_TEST_REFERENCE)
    assert row.tier == "grenzwertig"


def test_rejected_when_above_p99():
    r = _bench_row("8ABC", n_clash_pairs=50)  # > p99=30
    row = ar.classify_protein(r, reference_percentiles=_TEST_REFERENCE)
    assert row.tier == "rejected"
    assert row.action == "reject_or_reroll"


def test_rejected_when_gate_fails_even_if_metrics_ok():
    r = _bench_row("8ABC", gate_pass=False, n_clash_pairs=2)
    row = ar.classify_protein(r, reference_percentiles=_TEST_REFERENCE)
    assert row.tier == "rejected"


def test_rejected_when_error_field_present():
    r = {"pdb": "8ABC", "error": "MODELLER failed on chain A"}
    row = ar.classify_protein(r, reference_percentiles=_TEST_REFERENCE)
    assert row.tier == "rejected"
    assert "MODELLER" in row.notes


# ---------------------------------------------------------------------------
# per-metric band assignment
# ---------------------------------------------------------------------------


def test_band_boundaries():
    r90 = {"p90": 7, "p95": 12, "p99": 30}
    assert ar._band_for(7, r90) == "in_p90"
    assert ar._band_for(8, r90) == "p90-p95"
    assert ar._band_for(12, r90) == "p90-p95"
    assert ar._band_for(13, r90) == "p95-p99"
    assert ar._band_for(30, r90) == "p95-p99"
    assert ar._band_for(31, r90) == "above_p99"


def test_metric_bands_populated_in_row():
    r = _bench_row("8ABC", n_clash_pairs=8)
    row = ar.classify_protein(r, reference_percentiles=_TEST_REFERENCE)
    assert row.metrics["n_clash_pairs"]["band"] == "p90-p95"
    assert row.metrics["n_clash_pairs"]["ref_p90"] == 7


# ---------------------------------------------------------------------------
# build_audit_report + tier_summary
# ---------------------------------------------------------------------------


def test_build_audit_report_end_to_end():
    bench = [
        _bench_row("A", n_clash_pairs=3),   # volltreffer (no crystal, so no regression signal)
        _bench_row("B", n_clash_pairs=10),  # wahrscheinlich_ok
        _bench_row("C", n_clash_pairs=20),  # grenzwertig
        _bench_row("D", n_clash_pairs=50),  # rejected
    ]
    rows = ar.build_audit_report(bench, reference_percentiles=_TEST_REFERENCE)
    tiers = [r.tier for r in rows]
    assert tiers == ["volltreffer", "wahrscheinlich_ok", "grenzwertig", "rejected"]


def test_tier_summary_percentages():
    bench = [_bench_row(f"P{i}", n_clash_pairs=3) for i in range(3)]  # 3 volltreffer
    bench.append(_bench_row("X", n_clash_pairs=50))  # 1 rejected
    rows = ar.build_audit_report(bench, reference_percentiles=_TEST_REFERENCE)
    s = ar.tier_summary(rows)
    assert s["total"] == 4
    assert s["counts"]["volltreffer"] == 3
    assert s["counts"]["rejected"] == 1
    assert s["percentages"]["volltreffer"] == pytest.approx(75.0)


# ---------------------------------------------------------------------------
# CSV export — the primary user deliverable
# ---------------------------------------------------------------------------


def test_write_audit_csv_opens_cleanly(tmp_path: Path):
    bench = [
        _bench_row("8ABC", n_clash_pairs=3),
        _bench_row("8XYZ", n_clash_pairs=50, gate_pass=False,
                   reasons=["clash gained 47"],
                   non_planar_omega_residues=["A/123 GLY→ALA: ω=-45.2°"]),
    ]
    rows = ar.build_audit_report(bench, reference_percentiles=_TEST_REFERENCE)
    out = ar.write_audit_csv(rows, tmp_path / "audit.csv")
    assert out.is_file()
    # Re-read and verify
    with out.open() as fh:
        reader = csv.DictReader(fh)
        loaded = list(reader)
    assert len(loaded) == 2
    assert loaded[0]["pdb"] == "8ABC"
    assert loaded[0]["tier"] == "volltreffer"
    assert loaded[1]["tier"] == "rejected"
    assert "A/123" in loaded[1]["suspect_omega_non_planar"]
    assert "clash gained 47" in loaded[1]["reasons"]


def test_write_audit_csv_includes_all_ref_metric_columns(tmp_path: Path):
    rows = ar.build_audit_report([_bench_row("A")], reference_percentiles=_TEST_REFERENCE)
    out = ar.write_audit_csv(rows, tmp_path / "audit.csv")
    header = out.read_text().splitlines()[0]
    for metric in ar._REFERENCE_METRICS:
        for suffix in ("_crystal", "_fruton", "_delta", "_ref_p90", "_ref_p99", "_band"):
            assert f"{metric}{suffix}" in header


def test_write_audit_csv_handles_empty_input(tmp_path: Path):
    out = ar.write_audit_csv([], tmp_path / "audit.csv")
    assert out.is_file()
    assert "pdb" in out.read_text()


# ---------------------------------------------------------------------------
# markdown summary
# ---------------------------------------------------------------------------


def test_write_audit_markdown_lists_rejected(tmp_path: Path):
    bench = [
        _bench_row("A", n_clash_pairs=3),   # volltreffer
        _bench_row("X", n_clash_pairs=50, gate_pass=False,
                   reasons=["clash gained 47"]),
    ]
    rows = ar.build_audit_report(bench, reference_percentiles=_TEST_REFERENCE)
    out = ar.write_audit_markdown_summary(rows, tmp_path / "sum.md")
    text = out.read_text()
    assert "volltreffer" in text
    assert "rejected" in text
    assert "**X**" in text
    assert "clash gained 47" in text


# ---------------------------------------------------------------------------
# reviewer-defensibility invariants
# ---------------------------------------------------------------------------


def test_no_magic_thresholds_in_classifier():
    """Tier assignment must derive purely from reference bands + regression;
    no bench-specific literals."""
    import ast, inspect
    tree = ast.parse(inspect.getsource(ar.classify_protein))
    # Remove docstring
    if (tree.body and isinstance(tree.body[0], ast.FunctionDef)
            and tree.body[0].body and isinstance(tree.body[0].body[0], ast.Expr)):
        tree.body[0].body.pop(0)
    src = ast.unparse(tree)
    for magic in ("> 200", ">= 3", ">= 1", "> 20 or"):
        assert magic not in src, f"magic threshold {magic!r} in classifier"


def test_baseline_delta_populated_when_crystal_metrics_provided():
    r = _bench_row("A", n_clash_pairs=5)
    crystal = {"n_clash_pairs": 3, "n_omega_non_planar": 0,
               "n_omega_cis_nonpro": 0, "n_peptide_bonds_broken": 0}
    row = ar.classify_protein(r, crystal_metrics=crystal,
                              reference_percentiles=_TEST_REFERENCE)
    assert row.metrics["n_clash_pairs"]["delta"] == 2
    assert row.metrics["n_clash_pairs"]["crystal"] == 3
    assert row.metrics["n_clash_pairs"]["fruton"] == 5
