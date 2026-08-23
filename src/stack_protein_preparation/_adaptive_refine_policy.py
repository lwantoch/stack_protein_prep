"""Principle-driven adaptive-refine trigger policy.

REPLACES the hand-picked thresholds that were tuned to the specific 48
AF-available proteins in the MMBSA_200 bench (2026-08-23 iter-1..5).
Those were case-specific overfitting: reviewers would rightly reject
"clash > 200" and "ω_np ≥ 1" as magic numbers matched to observed
failure counts on the same set being benchmarked.

Principle 1 (trigger slow-refine)
    Trigger a slow-refine retry whenever FRUTON's fast-refine output is
    WORSE than the SAME protein's own crystal baseline on ANY of the
    quality metrics we care about (clash pairs, non-planar ω,
    cis-nonPro, broken peptide bonds).  This is a data-independent
    protein-relative comparison — no magic thresholds.

Principle 2 (ceiling — skip slow retry)
    Skip the slow-refine retry when the fast-refine output is *worse
    than any deposited crystal in the reference population*.  The
    ceiling comes from the p99 of a 199-crystal INDEPENDENT reference
    distribution stored in ``data/baseline_reference_percentiles.json``
    (source: FRUTON-NEW/ raw crystals, disjoint from the 48 bench
    proteins in intent — raw crystals vs crystal+AF-fill+refine).
    Rationale: if the fast output is worse than 99 % of independently-
    deposited crystals, MODELLER slow-refine's 5-conformer search will
    almost certainly not recover it — spending compute is wasted and
    (empirically) blocks the SLURM queue by stalling for >1 h.

Neither principle references the 48 bench proteins.  A new benchmark
set drops in without any threshold re-tuning.
"""
from __future__ import annotations

import json
from dataclasses import dataclass
from importlib import resources
from typing import Any


_REFERENCE_METRICS = (
    "n_clash_pairs",
    "n_omega_non_planar",
    "n_omega_cis_nonpro",
    "n_peptide_bonds_broken",
)


def load_reference_percentiles() -> dict[str, Any]:
    """Load the 199-crystal reference percentile table (package data)."""
    text = (
        resources.files("stack_protein_preparation")
        .joinpath("data/baseline_reference_percentiles.json")
        .read_text()
    )
    return json.loads(text)


@dataclass
class RefinePolicyDecision:
    action: str  # "retry_slow" | "skip_ceiling" | "no_retry"
    regressed_metrics: dict[str, tuple[int, int]]  # {metric: (baseline, fruton)}
    exceeded_ceilings: dict[str, tuple[int, int]]  # {metric: (fruton_value, p99_reference)}
    reason: str

    def to_dict(self) -> dict:
        return {
            "action": self.action,
            "regressed_metrics": {k: list(v) for k, v in self.regressed_metrics.items()},
            "exceeded_ceilings": {k: list(v) for k, v in self.exceeded_ceilings.items()},
            "reason": self.reason,
        }


def decide_refine_action(
    crystal_qc: Any,
    fruton_qc: Any,
    reference_percentiles: dict[str, Any] | None = None,
    ceiling_percentile: str = "p99",
) -> RefinePolicyDecision:
    """Compare a fast-refine output vs its crystal baseline + a reference
    population; return whether to retry with slow-refine, skip due to
    ceiling, or do nothing.

    Args:
        crystal_qc: quality-check result on the raw crystal (has integer
            attributes for each metric in ``_REFERENCE_METRICS``).
        fruton_qc: quality-check result on the fast-refined FRUTON output.
        reference_percentiles: optional override (mostly for tests); when
            None, loads the packaged 199-crystal table.
        ceiling_percentile: which percentile ('p90'/'p95'/'p99'/'max') to
            treat as the ceiling.  Default 'p99'.

    Returns a RefinePolicyDecision.
    """
    if reference_percentiles is None:
        reference_percentiles = load_reference_percentiles()

    regressed: dict[str, tuple[int, int]] = {}
    exceeded: dict[str, tuple[int, int]] = {}

    for metric in _REFERENCE_METRICS:
        baseline = int(getattr(crystal_qc, metric, 0))
        fruton = int(getattr(fruton_qc, metric, 0))
        if fruton > baseline:
            regressed[metric] = (baseline, fruton)

        ref = reference_percentiles.get(metric)
        if isinstance(ref, dict) and ceiling_percentile in ref:
            ceiling_val = int(ref[ceiling_percentile])
            if fruton > ceiling_val:
                exceeded[metric] = (fruton, ceiling_val)

    if exceeded:
        return RefinePolicyDecision(
            action="skip_ceiling",
            regressed_metrics=regressed,
            exceeded_ceilings=exceeded,
            reason=(
                f"fruton output exceeds {ceiling_percentile} of 199-crystal "
                f"reference on {sorted(exceeded)}; slow-refine unlikely to "
                f"recover, deferring to rollback"
            ),
        )
    if regressed:
        return RefinePolicyDecision(
            action="retry_slow",
            regressed_metrics=regressed,
            exceeded_ceilings={},
            reason=(
                f"fruton fast-refine regressed on {sorted(regressed)} vs "
                f"crystal baseline; retrying with slow-refine (5 conformers)"
            ),
        )
    return RefinePolicyDecision(
        action="no_retry",
        regressed_metrics={},
        exceeded_ceilings={},
        reason="fruton fast-refine did not regress any quality metric",
    )
