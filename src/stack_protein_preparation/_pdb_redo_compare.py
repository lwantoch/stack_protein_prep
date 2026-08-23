"""PDB-REDO cross-comparison harness (Nature R2 concern).

Reviewer perspective: FRUTON produces gap-filled + refined models from
crystal templates + AF splices.  A natural sanity check is:
> does the FRUTON output compare favourably to the same PDB entry's
  PDB-REDO re-refinement?

PDB-REDO (Joosten et al. 2014, https://pdb-redo.eu) is the de-facto
community benchmark for improving deposited crystal-structure quality:
it re-refines every PDB entry with modern software + geometric
restraints and publishes per-metric before/after JSON summaries.

This module is the *analysis* layer.  It does not fetch PDB-REDO --
that requires an external HTTP call the user runs on demand.  Given
two directories of per-PDB metric JSONs (one from FRUTON, one from
PDB-REDO), it:

1. Names the metric fields canonically (METRIC_NAMES).
2. Loads matched per-PDB records by shared pdb_id.
3. Computes per-metric per-PDB Δ (fruton − pdb_redo).  Negative Δ
   means FRUTON is closer to zero (usually better) for outlier /
   clash metrics; the direction depends on the metric.

Pure Python + stdlib.  License-free.
"""
from __future__ import annotations

import json
import statistics
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable


# Canonical metric catalogue.  Each entry lists the direction that is
# "better" (lower_is_better=True means smaller values are better --
# applies to outlier counts, clashscore, non-planar-ω, cis-nonPro).
METRIC_NAMES: dict[str, dict[str, object]] = {
    "n_rama_outlier":     {"units": "count",       "lower_is_better": True,
                           "reviewer_label": "Ramachandran outliers"},
    "n_clash_pairs":      {"units": "count",       "lower_is_better": True,
                           "reviewer_label": "Steric clash pairs"},
    "n_omega_non_planar": {"units": "count",       "lower_is_better": True,
                           "reviewer_label": "Non-planar ω peptide bonds"},
    "n_omega_cis_nonpro": {"units": "count",       "lower_is_better": True,
                           "reviewer_label": "cis-nonPro peptide bonds"},
    "clashscore":         {"units": "per 1000 at", "lower_is_better": True,
                           "reviewer_label": "MolProbity clashscore"},
    "rama_favoured_pct":  {"units": "%",           "lower_is_better": False,
                           "reviewer_label": "Ramachandran favoured %"},
    "rotamer_outlier_pct":{"units": "%",           "lower_is_better": True,
                           "reviewer_label": "Rotamer outlier %"},
    "bond_rms_z":         {"units": "z",           "lower_is_better": True,
                           "reviewer_label": "Bond RMS Z"},
    "angle_rms_z":        {"units": "z",           "lower_is_better": True,
                           "reviewer_label": "Angle RMS Z"},
}


@dataclass
class PerPdbMetricDelta:
    """One (pdb, metric) delta between FRUTON and PDB-REDO."""
    pdb: str
    metric: str
    fruton_value: float | None
    pdb_redo_value: float | None

    def delta(self) -> float | None:
        if self.fruton_value is None or self.pdb_redo_value is None:
            return None
        return float(self.fruton_value) - float(self.pdb_redo_value)

    def better_side(self, lower_is_better: bool) -> str:
        d = self.delta()
        if d is None:
            return "unavailable"
        if d == 0:
            return "tie"
        if lower_is_better:
            return "fruton" if d < 0 else "pdb_redo"
        # higher is better (e.g. rama_favoured_pct)
        return "fruton" if d > 0 else "pdb_redo"


@dataclass
class PerMetricSummary:
    metric: str
    reviewer_label: str
    units: str
    lower_is_better: bool
    n_matched: int
    n_fruton_better: int
    n_pdb_redo_better: int
    n_tie: int
    n_unavailable: int
    mean_delta: float | None
    median_delta: float | None

    def to_dict(self) -> dict:
        return {
            "metric": self.metric,
            "reviewer_label": self.reviewer_label,
            "units": self.units,
            "lower_is_better": self.lower_is_better,
            "n_matched": self.n_matched,
            "n_fruton_better": self.n_fruton_better,
            "n_pdb_redo_better": self.n_pdb_redo_better,
            "n_tie": self.n_tie,
            "n_unavailable": self.n_unavailable,
            "mean_delta": (None if self.mean_delta is None
                           else round(self.mean_delta, 4)),
            "median_delta": (None if self.median_delta is None
                             else round(self.median_delta, 4)),
        }


@dataclass
class PdbRedoComparison:
    per_metric: list[PerMetricSummary]
    per_pdb_deltas: dict[str, list[PerPdbMetricDelta]] = field(default_factory=dict)

    def to_dict(self) -> dict:
        return {
            "per_metric": [m.to_dict() for m in self.per_metric],
            "per_pdb_deltas": {
                pdb: [
                    {
                        "metric": d.metric,
                        "fruton_value": d.fruton_value,
                        "pdb_redo_value": d.pdb_redo_value,
                        "delta": d.delta(),
                    }
                    for d in deltas
                ]
                for pdb, deltas in self.per_pdb_deltas.items()
            },
        }


def load_metrics_dir(directory: str | Path) -> dict[str, dict]:
    """Load ``<pdb_id>.json`` files under ``directory`` into ``{pdb: metrics}``.

    Files must be JSON objects (not lists) whose top level is the metric
    map.  The ``pdb_id`` is the file stem (case-preserving; uppercased
    for comparison consistency).
    """
    d = Path(directory)
    if not d.is_dir():
        raise FileNotFoundError(f"metrics directory not found: {d}")
    out: dict[str, dict] = {}
    for path in sorted(d.glob("*.json")):
        with path.open() as fh:
            payload = json.load(fh)
        if not isinstance(payload, dict):
            raise ValueError(f"{path} is not a JSON object at top level")
        out[path.stem.upper()] = payload
    return out


def _to_float_or_none(value) -> float | None:
    if value is None:
        return None
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def build_comparison(
    fruton_dir: str | Path,
    pdb_redo_dir: str | Path,
    metrics: Iterable[str] | None = None,
) -> PdbRedoComparison:
    """Assemble a PdbRedoComparison from two directories of per-PDB JSONs.

    Args:
        fruton_dir: path holding ``<PDB>.json`` records emitted by
            FRUTON's quality-check summariser.
        pdb_redo_dir: path holding matching ``<PDB>.json`` records
            extracted from PDB-REDO's per-entry summary JSON.
        metrics: subset of ``METRIC_NAMES`` keys to compare; ``None``
            compares every metric present in either source.
    """
    fruton = load_metrics_dir(fruton_dir)
    redo = load_metrics_dir(pdb_redo_dir)

    if metrics is None:
        metric_keys = list(METRIC_NAMES.keys())
    else:
        metric_keys = list(metrics)
        for k in metric_keys:
            if k not in METRIC_NAMES:
                raise ValueError(
                    f"unknown metric {k!r}; known: {sorted(METRIC_NAMES)}"
                )

    shared_pdbs = sorted(set(fruton) & set(redo))

    per_metric: list[PerMetricSummary] = []
    per_pdb: dict[str, list[PerPdbMetricDelta]] = {p: [] for p in shared_pdbs}

    for m in metric_keys:
        spec = METRIC_NAMES[m]
        lower_is_better = bool(spec["lower_is_better"])
        deltas_for_metric: list[float] = []
        n_fruton = n_redo = n_tie = n_unavailable = 0
        for pdb in shared_pdbs:
            f_val = _to_float_or_none(fruton[pdb].get(m))
            r_val = _to_float_or_none(redo[pdb].get(m))
            entry = PerPdbMetricDelta(
                pdb=pdb, metric=m,
                fruton_value=f_val, pdb_redo_value=r_val,
            )
            per_pdb[pdb].append(entry)
            side = entry.better_side(lower_is_better)
            if side == "unavailable":
                n_unavailable += 1
            elif side == "tie":
                n_tie += 1
            elif side == "fruton":
                n_fruton += 1
            else:
                n_redo += 1
            d = entry.delta()
            if d is not None:
                deltas_for_metric.append(d)
        per_metric.append(PerMetricSummary(
            metric=m,
            reviewer_label=str(spec["reviewer_label"]),
            units=str(spec["units"]),
            lower_is_better=lower_is_better,
            n_matched=len(shared_pdbs),
            n_fruton_better=n_fruton,
            n_pdb_redo_better=n_redo,
            n_tie=n_tie,
            n_unavailable=n_unavailable,
            mean_delta=(statistics.fmean(deltas_for_metric) if deltas_for_metric else None),
            median_delta=(statistics.median(deltas_for_metric) if deltas_for_metric else None),
        ))

    return PdbRedoComparison(per_metric=per_metric, per_pdb_deltas=per_pdb)


def format_comparison_table(comp: PdbRedoComparison) -> str:
    """Reviewer-facing markdown table.  One row per metric."""
    lines: list[str] = []
    lines.append(
        "| Metric | Direction | n_matched | FRUTON better | PDB-REDO better | tie | unavail | mean Δ | median Δ |"
    )
    lines.append(
        "|---|---|---:|---:|---:|---:|---:|---:|---:|"
    )
    for m in comp.per_metric:
        direction = "↓ better" if m.lower_is_better else "↑ better"
        mean_str = "n/a" if m.mean_delta is None else f"{m.mean_delta:+.3f}"
        med_str = "n/a" if m.median_delta is None else f"{m.median_delta:+.3f}"
        lines.append(
            f"| {m.reviewer_label} | {direction} | {m.n_matched} | "
            f"{m.n_fruton_better} | {m.n_pdb_redo_better} | "
            f"{m.n_tie} | {m.n_unavailable} | {mean_str} | {med_str} |"
        )
    return "\n".join(lines) + "\n"
