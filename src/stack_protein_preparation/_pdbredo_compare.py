"""PDB-REDO comparison harness (Nature R2 concern).

Reviewer perspective: "you preserved chirality / peptide bonds / metal
geometry — but how do your final models compare to the community
gold-standard re-refinement (PDB-REDO) of the same PDB entries?"

This module is analysis-only.  It does not fetch PDB-REDO (that
requires network + user-side batching); it consumes a JSON that the
user prepares from the PDB-REDO REST API or their local mirror.
Given:

- a FRUTON bench-results JSON (list of per-protein dicts with metrics),
- a PDB-REDO metrics JSON (dict of pdb_id → metrics),

it produces per-PDB, per-metric deltas and aggregate summary stats.

Canonical PDB-REDO metric names (matching their JSON output):
- clashscore
- rama_favoured  (% Ramachandran favoured)
- rama_outlier   (% Ramachandran outlier)
- rotamer_outlier
- rms_bond
- rms_angle
- r_free

Pure stdlib.  License-free.
"""
from __future__ import annotations

import json
import statistics
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable


# Metric field-name mapping: FRUTON bench-json key -> canonical PDB-REDO name.
# Extend as new gate metrics are wired into the bench output.
FRUTON_TO_CANONICAL: dict[str, str] = {
    "clash": "clashscore",
    "n_rama_outlier_pct": "rama_outlier",
    "n_rama_favoured_pct": "rama_favoured",
    "n_rotamer_outlier_pct": "rotamer_outlier",
    "rms_bond_A": "rms_bond",
    "rms_angle_deg": "rms_angle",
    "r_free": "r_free",
}

# Metrics for which "smaller is better" (used to pick sign of Δ).
LOWER_IS_BETTER: set[str] = {
    "clashscore", "rama_outlier", "rotamer_outlier",
    "rms_bond", "rms_angle", "r_free",
}

# Metrics for which "larger is better".
HIGHER_IS_BETTER: set[str] = {"rama_favoured"}


@dataclass
class PdbRedoMetrics:
    """Canonical per-PDB PDB-REDO metric bundle."""
    pdb_id: str
    clashscore: float | None = None
    rama_favoured: float | None = None
    rama_outlier: float | None = None
    rotamer_outlier: float | None = None
    rms_bond: float | None = None
    rms_angle: float | None = None
    r_free: float | None = None

    def to_dict(self) -> dict:
        return {
            "pdb_id": self.pdb_id,
            "clashscore": self.clashscore,
            "rama_favoured": self.rama_favoured,
            "rama_outlier": self.rama_outlier,
            "rotamer_outlier": self.rotamer_outlier,
            "rms_bond": self.rms_bond,
            "rms_angle": self.rms_angle,
            "r_free": self.r_free,
        }


@dataclass
class MetricDelta:
    pdb_id: str
    metric: str
    fruton_value: float | None
    pdbredo_value: float | None

    def delta(self) -> float | None:
        """FRUTON − PDB-REDO (raw arithmetic delta)."""
        if self.fruton_value is None or self.pdbredo_value is None:
            return None
        return self.fruton_value - self.pdbredo_value

    def fruton_better(self) -> bool | None:
        """True if FRUTON's value is better than PDB-REDO's for this metric.

        Returns None when either value is missing or when the metric has
        no known better-direction convention.
        """
        d = self.delta()
        if d is None:
            return None
        if self.metric in LOWER_IS_BETTER:
            return d < 0
        if self.metric in HIGHER_IS_BETTER:
            return d > 0
        return None


@dataclass
class AggregateMetricStats:
    metric: str
    n_pairs: int
    mean_delta: float
    median_delta: float
    n_fruton_better: int
    n_pdbredo_better: int
    n_tied: int

    def to_dict(self) -> dict:
        return {
            "metric": self.metric,
            "n_pairs": self.n_pairs,
            "mean_delta": round(self.mean_delta, 4),
            "median_delta": round(self.median_delta, 4),
            "n_fruton_better": self.n_fruton_better,
            "n_pdbredo_better": self.n_pdbredo_better,
            "n_tied": self.n_tied,
        }


def load_pdbredo_metrics_json(path: str | Path) -> dict[str, PdbRedoMetrics]:
    """Load a mapping ``{pdb_id: {metric: value, ...}}`` into typed records.

    The JSON must be a dict keyed by PDB ID (case-insensitive; upper-cased
    internally).  Unknown metric fields are dropped silently — we do not
    complain about PDB-REDO adding new fields we haven't wired yet.
    """
    path = Path(path)
    with path.open() as fh:
        raw = json.load(fh)
    if not isinstance(raw, dict):
        raise ValueError(f"{path} is not a dict of pdb_id -> metrics")

    known = {
        "clashscore", "rama_favoured", "rama_outlier",
        "rotamer_outlier", "rms_bond", "rms_angle", "r_free",
    }
    out: dict[str, PdbRedoMetrics] = {}
    for pdb_id, metrics in raw.items():
        if not isinstance(metrics, dict):
            continue
        kept = {k: v for k, v in metrics.items() if k in known and v is not None}
        out[pdb_id.upper()] = PdbRedoMetrics(pdb_id=pdb_id.upper(), **kept)
    return out


def _canonical_fruton_metrics(record: dict) -> dict[str, float]:
    """Translate a FRUTON per-PDB dict into canonical-metric-name values.

    Missing fields are silently omitted so that partial per-PDB records
    still contribute to whatever metrics they do have.
    """
    out: dict[str, float] = {}
    for fruton_key, canonical in FRUTON_TO_CANONICAL.items():
        val = record.get(fruton_key)
        if val is None:
            continue
        try:
            out[canonical] = float(val)
        except (TypeError, ValueError):
            continue
    return out


def compare_bench_to_pdbredo(
    fruton_bench: list[dict],
    pdbredo_metrics: dict[str, PdbRedoMetrics],
) -> list[MetricDelta]:
    """Emit one MetricDelta per (pdb_id, metric) pair where both values exist.

    Args:
        fruton_bench: list of per-protein dicts as produced by
            scripts/CESGA_SLURM/fruton_bench_mmbsa200_full.py.
        pdbredo_metrics: mapping produced by ``load_pdbredo_metrics_json``.
    """
    deltas: list[MetricDelta] = []
    for record in fruton_bench:
        pdb_id = (record.get("pdb") or "").upper()
        if not pdb_id:
            continue
        pdbredo = pdbredo_metrics.get(pdb_id)
        if pdbredo is None:
            continue
        fruton_metrics = _canonical_fruton_metrics(record)
        for metric, fruton_val in fruton_metrics.items():
            pdbredo_val = getattr(pdbredo, metric, None)
            if pdbredo_val is None:
                continue
            deltas.append(MetricDelta(
                pdb_id=pdb_id,
                metric=metric,
                fruton_value=fruton_val,
                pdbredo_value=pdbredo_val,
            ))
    return deltas


def aggregate_deltas_by_metric(
    deltas: Iterable[MetricDelta],
) -> dict[str, AggregateMetricStats]:
    """Group deltas by metric and return per-metric summary stats."""
    by_metric: dict[str, list[MetricDelta]] = {}
    for d in deltas:
        by_metric.setdefault(d.metric, []).append(d)

    out: dict[str, AggregateMetricStats] = {}
    for metric, ds in by_metric.items():
        raw_values = [d.delta() for d in ds if d.delta() is not None]
        if not raw_values:
            continue
        n_fruton_better = sum(1 for d in ds if d.fruton_better() is True)
        n_pdbredo_better = sum(1 for d in ds if d.fruton_better() is False)
        n_tied = len(ds) - n_fruton_better - n_pdbredo_better
        out[metric] = AggregateMetricStats(
            metric=metric,
            n_pairs=len(raw_values),
            mean_delta=statistics.fmean(raw_values),
            median_delta=statistics.median(raw_values),
            n_fruton_better=n_fruton_better,
            n_pdbredo_better=n_pdbredo_better,
            n_tied=n_tied,
        )
    return out


def format_comparison_table(
    stats: dict[str, AggregateMetricStats],
) -> str:
    """Reviewer-friendly markdown row per metric."""
    if not stats:
        return "(no overlapping metrics between FRUTON bench and PDB-REDO reference)\n"

    lines = [
        "| Metric | n | ⟨Δ⟩ | median Δ | FRUTON better | PDB-REDO better | tied |",
        "|---|---:|---:|---:|---:|---:|---:|",
    ]
    for metric in sorted(stats):
        s = stats[metric]
        lines.append(
            f"| {metric} | {s.n_pairs} | {s.mean_delta:+.3f} | "
            f"{s.median_delta:+.3f} | {s.n_fruton_better} | "
            f"{s.n_pdbredo_better} | {s.n_tied} |"
        )
    return "\n".join(lines) + "\n"


def summarise_delta_direction(
    stats: dict[str, AggregateMetricStats],
) -> list[str]:
    """One reviewer-facing sentence per metric describing the direction."""
    lines: list[str] = []
    for metric in sorted(stats):
        s = stats[metric]
        winner = (
            "FRUTON"
            if s.n_fruton_better > s.n_pdbredo_better
            else "PDB-REDO"
            if s.n_pdbredo_better > s.n_fruton_better
            else "tie"
        )
        lines.append(
            f"{metric}: {winner} wins on {max(s.n_fruton_better, s.n_pdbredo_better)}"
            f"/{s.n_pairs} PDBs (⟨Δ⟩ = {s.mean_delta:+.3f})"
        )
    return lines
