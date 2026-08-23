"""FRUTON audit report — the primary user-facing deliverable.

REFRAMED 2026-08-23 (user feedback): FRUTON's value proposition is NOT
'produces perfect gap-fills' but 'produces honest, per-protein audit
documentation that flags every suspicious case for reviewer inspection.'

This module classifies each processed protein into one of four tiers
based on FRUTON's own quality-check output + the same INDEPENDENT
reference distribution (data/baseline_reference_percentiles.json) that
_adaptive_refine_policy uses.  The tier + a per-protein suspects table
becomes a CSV that opens in Excel.

Tiers (progressively lower confidence, uses the same 199-crystal
reference distribution so nothing here is bench-specific):

  volltreffer
      No quality-metric regression vs the protein's own crystal
      baseline AND all metrics inside the reference p90 band.
      "Ship it, MMBSA-ready."

  wahrscheinlich_ok
      Small regression on ≥1 metric but everything still inside the
      reference p95 band; or zero-regression with a couple metrics
      inside p95-p99.  "Likely fine, quick eyeball recommended."

  grenzwertig
      Regression or output between reference p95 and p99.
      "Borderline — user should manually inspect the flagged residues."

  rejected
      Any metric exceeds the reference p99 band OR the FRUTON gate
      already rejected the fill.  "Do not use as-is; either drop the
      fill (REMARK 465) or re-refine externally."

The classifier + CSV writer are pure Python + stdlib (uses ``csv`` from
the standard library; no pandas dependency).  License-free.
"""
from __future__ import annotations

import csv
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Iterable

from stack_protein_preparation._adaptive_refine_policy import (
    _REFERENCE_METRICS,
    load_reference_percentiles,
)


TIER_ORDER = ("volltreffer", "wahrscheinlich_ok", "grenzwertig", "rejected")

# Per-tier reviewer-facing suggested action.
TIER_ACTION = {
    "volltreffer":      "accept",
    "wahrscheinlich_ok": "accept_with_note",
    "grenzwertig":      "manual_review",
    "rejected":         "reject_or_reroll",
}


@dataclass
class ProteinAuditRow:
    """One protein's audit-report row.  Contains everything a reviewer
    needs to triage the model without opening a PDB viewer."""
    pdb: str
    tier: str
    gate_pass: bool
    action: str
    n_gaps: int
    delta_n: int
    n_rolled_gaps: int
    # Per-metric values + baseline + reference-percentile band
    metrics: dict[str, dict[str, Any]] = field(default_factory=dict)
    # Which residues each gate flagged
    suspect_residues: dict[str, list[str]] = field(default_factory=dict)
    reasons: list[str] = field(default_factory=list)
    notes: str = ""

    def to_csv_row(self) -> dict[str, Any]:
        base = {
            "pdb": self.pdb,
            "tier": self.tier,
            "action": self.action,
            "gate_pass": self.gate_pass,
            "n_gaps": self.n_gaps,
            "delta_n": self.delta_n,
            "n_rolled_gaps": self.n_rolled_gaps,
        }
        for metric in _REFERENCE_METRICS:
            m = self.metrics.get(metric, {})
            base[f"{metric}_crystal"] = m.get("crystal", "")
            base[f"{metric}_fruton"] = m.get("fruton", "")
            base[f"{metric}_delta"] = m.get("delta", "")
            base[f"{metric}_ref_p90"] = m.get("ref_p90", "")
            base[f"{metric}_ref_p99"] = m.get("ref_p99", "")
            base[f"{metric}_band"] = m.get("band", "")  # in_p90 / p90-p95 / p95-p99 / above_p99
        # Suspect-residue lists (comma-joined for CSV)
        for cat in ("omega_non_planar", "omega_cis_nonpro", "rama_outlier",
                    "chirality", "clash"):
            base[f"suspect_{cat}"] = "; ".join(self.suspect_residues.get(cat, []))
        base["reasons"] = " | ".join(self.reasons)
        base["notes"] = self.notes
        return base


def _band_for(value: int, ref_entry: dict[str, int]) -> str:
    if value <= ref_entry.get("p90", 0):
        return "in_p90"
    if value <= ref_entry.get("p95", 0):
        return "p90-p95"
    if value <= ref_entry.get("p99", 0):
        return "p95-p99"
    return "above_p99"


def classify_protein(
    result: dict[str, Any],
    crystal_metrics: dict[str, int] | None = None,
    reference_percentiles: dict[str, Any] | None = None,
) -> ProteinAuditRow:
    """Turn one per-protein bench result dict into a ProteinAuditRow.

    Args:
        result: bench-output row (as emitted by fruton_bench_mmbsa200_full).
            Must contain at least ``pdb``; other fields default when absent.
        crystal_metrics: optional dict of raw-crystal quality metrics for
            THIS protein (e.g. from baseline_quality_check_full).  When
            provided, delta = fruton − crystal is filled in per metric.
        reference_percentiles: optional override; loads packaged data
            (199-crystal reference) when None.
    """
    if reference_percentiles is None:
        reference_percentiles = load_reference_percentiles()

    pdb = str(result.get("pdb", "?"))
    gate_pass = bool(result.get("gate_pass"))

    per_metric: dict[str, dict[str, Any]] = {}
    highest_band = "in_p90"
    band_rank = {"in_p90": 0, "p90-p95": 1, "p95-p99": 2, "above_p99": 3}
    any_regression = False
    for metric in _REFERENCE_METRICS:
        fruton_val = int(result.get(metric, 0))
        crystal_val = int((crystal_metrics or {}).get(metric, 0))
        ref_entry = reference_percentiles.get(metric, {})
        band = _band_for(fruton_val, ref_entry) if ref_entry else "in_p90"
        per_metric[metric] = {
            "crystal": crystal_val if crystal_metrics else "",
            "fruton": fruton_val,
            "delta": (fruton_val - crystal_val) if crystal_metrics else "",
            "ref_p90": ref_entry.get("p90", ""),
            "ref_p99": ref_entry.get("p99", ""),
            "band": band,
        }
        if band_rank[band] > band_rank[highest_band]:
            highest_band = band
        if crystal_metrics and fruton_val > crystal_val:
            any_regression = True

    # Tier assignment — no bench-specific magic numbers, only reference bands
    # + protein-relative regression signal + upstream gate decision.
    if not gate_pass or highest_band == "above_p99":
        tier = "rejected"
    elif highest_band == "p95-p99":
        tier = "grenzwertig"
    elif highest_band == "p90-p95" or any_regression:
        tier = "wahrscheinlich_ok"
    else:
        tier = "volltreffer"

    suspect: dict[str, list[str]] = {
        "omega_non_planar": list(result.get("non_planar_omega_residues", []) or []),
        "omega_cis_nonpro": list(result.get("cis_nonpro_omega_residues", []) or []),
        "rama_outlier":     list(result.get("rama_outlier_residues", []) or []),
        "chirality":        list(result.get("chirality_flagged_residues", []) or []),
        "clash":            list(result.get("clash_examples", []) or []),
    }

    notes = ""
    if "error" in result:
        tier = "rejected"
        notes = f"pipeline error: {result['error']}"

    return ProteinAuditRow(
        pdb=pdb,
        tier=tier,
        gate_pass=gate_pass,
        action=TIER_ACTION[tier],
        n_gaps=int(result.get("n_gaps", 0)),
        delta_n=int(result.get("delta_n", 0)),
        n_rolled_gaps=int(result.get("rolled_gaps", 0)),
        metrics=per_metric,
        suspect_residues=suspect,
        reasons=list(result.get("reasons", []) or []),
        notes=notes,
    )


def build_audit_report(
    bench_results: Iterable[dict[str, Any]],
    baseline_qc_results: Iterable[dict[str, Any]] | None = None,
    reference_percentiles: dict[str, Any] | None = None,
) -> list[ProteinAuditRow]:
    """Classify a whole bench in one call.

    baseline_qc_results is the output of baseline_quality_check_full;
    when provided, per-metric delta vs crystal is populated.
    """
    baseline_by_pdb: dict[str, dict[str, int]] = {}
    for x in (baseline_qc_results or []):
        pid = x.get("pdb")
        if pid:
            baseline_by_pdb[str(pid).upper()] = {
                m: int(x.get(m, 0)) for m in _REFERENCE_METRICS
            }
    return [
        classify_protein(
            r,
            crystal_metrics=baseline_by_pdb.get(str(r.get("pdb", "")).upper()),
            reference_percentiles=reference_percentiles,
        )
        for r in bench_results
    ]


def tier_summary(rows: Iterable[ProteinAuditRow]) -> dict[str, Any]:
    counts = {t: 0 for t in TIER_ORDER}
    for row in rows:
        counts[row.tier] = counts.get(row.tier, 0) + 1
    total = sum(counts.values())
    pct = {t: (100 * counts[t] / total if total else 0.0) for t in TIER_ORDER}
    return {"total": total, "counts": counts, "percentages": pct}


def write_audit_csv(rows: Iterable[ProteinAuditRow], output_path: str | Path) -> Path:
    """Emit a CSV that opens directly in Excel; one row per protein."""
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    rows = list(rows)
    if not rows:
        # Still create the file with just a header hint
        output_path.write_text("pdb,tier,action,gate_pass,n_gaps,delta_n,notes\n")
        return output_path
    fieldnames = list(rows[0].to_csv_row().keys())
    with output_path.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=fieldnames)
        w.writeheader()
        for row in rows:
            w.writerow(row.to_csv_row())
    return output_path


def write_audit_markdown_summary(
    rows: Iterable[ProteinAuditRow],
    output_path: str | Path,
) -> Path:
    """Reviewer-facing one-page markdown: tier breakdown + rejected list."""
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    rows = list(rows)
    summary = tier_summary(rows)

    lines: list[str] = []
    lines.append("# FRUTON audit summary")
    lines.append("")
    lines.append(f"**{summary['total']} proteins** classified against the "
                 "199-crystal reference distribution")
    lines.append("")
    lines.append("| Tier | Action | n | % |")
    lines.append("|---|---|---:|---:|")
    for t in TIER_ORDER:
        lines.append(
            f"| {t} | {TIER_ACTION[t]} | {summary['counts'][t]} | "
            f"{summary['percentages'][t]:.1f} |"
        )
    lines.append("")

    for tier in ("rejected", "grenzwertig"):
        tier_rows = [r for r in rows if r.tier == tier]
        if not tier_rows:
            continue
        lines.append(f"## {tier} ({len(tier_rows)})")
        lines.append("")
        for r in tier_rows:
            reasons = "; ".join(r.reasons[:3]) if r.reasons else r.notes
            lines.append(f"- **{r.pdb}** — Δn={r.delta_n:+d} rolled={r.n_rolled_gaps} — {reasons}")
        lines.append("")

    output_path.write_text("\n".join(lines) + "\n")
    return output_path
