"""Per-family stratified analysis of FRUTON bench results (Nature R1).

Reviewer concern (Nature R1): pipeline-wide aggregate hides class
imbalance.  A 90% pass rate on the MMBSA_200 bench might come from
95% on kinases (which have plentiful crystal templates) and 60% on
GPCRs (which have flexible ICLs the crystals rarely resolve).

This module is the *analysis* layer.  It accepts:

1. A bench-results JSON (list of per-protein dicts, matching the
   shape produced by scripts/CESGA_SLURM/fruton_bench_mmbsa200_full.py).
2. A ``family_by_pdb`` mapping ``{"8ABC": "kinase", ...}`` supplied
   by the user (or built up over time in the repo; see
   scripts/build_family_assignment.py — future).

Produces per-family aggregate + per-family per-PDB details for a
reviewer figure.

Pure Python + stdlib.  License-free.
"""
from __future__ import annotations

import json
import statistics
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable


# Canonical family taxonomy.  User-supplied mappings should use these
# labels so the aggregate report stays comparable across bench runs.
# ``__unassigned__`` is a synthetic bucket for PDBs missing from the map.
FAMILY_LABELS: tuple[str, ...] = (
    "kinase",
    "gpcr",
    "protease",
    "metalloenzyme",
    "phosphatase",
    "nuclear_receptor",
    "hydrolase",
    "transferase",
    "oxidoreductase",
    "cofactor_dependent",  # NAD/FAD/HEM/SAM/PLP-carrying enzymes not else classified
    "multi_domain",
    "other",
    "__unassigned__",
)


@dataclass
class FamilyAggregate:
    family: str
    n_pdbs: int
    n_passed: int
    pass_rate: float
    mean_delta_n: float
    total_delta_n: int
    mean_clash: float
    mean_brk: float
    pdb_ids: list[str] = field(default_factory=list)

    def to_dict(self) -> dict:
        return {
            "family": self.family,
            "n_pdbs": self.n_pdbs,
            "n_passed": self.n_passed,
            "pass_rate": round(self.pass_rate, 4),
            "mean_delta_n": round(self.mean_delta_n, 3),
            "total_delta_n": self.total_delta_n,
            "mean_clash": round(self.mean_clash, 3),
            "mean_brk": round(self.mean_brk, 3),
            "pdb_ids": list(self.pdb_ids),
        }


@dataclass
class FamilyStratification:
    per_family: list[FamilyAggregate]
    overall: FamilyAggregate

    def to_dict(self) -> dict:
        return {
            "per_family": [a.to_dict() for a in self.per_family],
            "overall": self.overall.to_dict(),
        }


def load_family_mapping(path: str | Path) -> dict[str, str]:
    """Load a JSON ``{pdb_id: family_label}`` mapping.

    PDB IDs are uppercased; unknown family labels are kept as-is (the
    reviewer can inspect them in the diagnostic column even if they
    fall outside FAMILY_LABELS).
    """
    p = Path(path)
    with p.open() as fh:
        raw = json.load(fh)
    if not isinstance(raw, dict):
        raise ValueError(f"{p} is not a JSON object at top level")
    out: dict[str, str] = {}
    for pdb, family in raw.items():
        if not isinstance(family, str):
            continue
        out[str(pdb).upper()] = family
    return out


def _safe_mean(values: Iterable[float]) -> float:
    values = list(values)
    return statistics.fmean(values) if values else 0.0


def _aggregate(name: str, records: list[dict]) -> FamilyAggregate:
    n = len(records)
    n_passed = sum(1 for r in records if r.get("gate_pass"))
    return FamilyAggregate(
        family=name,
        n_pdbs=n,
        n_passed=n_passed,
        pass_rate=(n_passed / n) if n else 0.0,
        mean_delta_n=_safe_mean(int(r.get("delta_n", 0)) for r in records),
        total_delta_n=sum(int(r.get("delta_n", 0)) for r in records),
        mean_clash=_safe_mean(float(r.get("clash", 0)) for r in records),
        mean_brk=_safe_mean(float(r.get("brk", 0)) for r in records),
        pdb_ids=sorted({str(r.get("pdb", "")).upper() for r in records if r.get("pdb")}),
    )


def stratify_bench_by_family(
    bench_results: list[dict],
    family_by_pdb: dict[str, str],
) -> FamilyStratification:
    """Group ``bench_results`` by ``family_by_pdb`` and emit aggregates.

    PDBs missing from ``family_by_pdb`` land in ``__unassigned__``.
    An empty family is dropped from the per-family list (so the
    reviewer table only carries populated rows).
    """
    by_family: dict[str, list[dict]] = {}
    for r in bench_results:
        pdb = str(r.get("pdb", "")).upper()
        if not pdb:
            continue
        family = family_by_pdb.get(pdb, "__unassigned__")
        by_family.setdefault(family, []).append(r)

    per_family: list[FamilyAggregate] = []
    for family, recs in by_family.items():
        agg = _aggregate(family, recs)
        if agg.n_pdbs > 0:
            per_family.append(agg)

    # Sort: canonical families first (in FAMILY_LABELS order), then
    # any unknown labels alphabetically, then __unassigned__ last.
    order = {name: i for i, name in enumerate(FAMILY_LABELS)}

    def sort_key(a: FamilyAggregate) -> tuple:
        if a.family == "__unassigned__":
            return (2, "")
        if a.family in order:
            return (0, order[a.family])
        return (1, a.family)

    per_family.sort(key=sort_key)
    overall = _aggregate("overall", bench_results)
    return FamilyStratification(per_family=per_family, overall=overall)


def format_stratification_table(strat: FamilyStratification) -> str:
    """Reviewer-facing markdown table with one row per family."""
    lines: list[str] = []
    lines.append(
        "| Family | n_pdbs | pass % | rescued | ⟨ΔN⟩ | ⟨clash⟩ | ⟨brk⟩ |"
    )
    lines.append("|---|---:|---:|---:|---:|---:|---:|")
    for a in strat.per_family + [strat.overall]:
        lines.append(
            f"| {a.family} | {a.n_pdbs} | {100 * a.pass_rate:.1f}% | "
            f"{a.total_delta_n} | {a.mean_delta_n:.2f} | "
            f"{a.mean_clash:.2f} | {a.mean_brk:.2f} |"
        )
    return "\n".join(lines) + "\n"


def imbalance_summary(strat: FamilyStratification) -> str:
    """Short reviewer-facing highlight of best / worst family + overall."""
    if not strat.per_family:
        return "no populated families"
    best = max(strat.per_family, key=lambda a: a.pass_rate)
    worst = min(strat.per_family, key=lambda a: a.pass_rate)
    return (
        f"overall pass rate {100 * strat.overall.pass_rate:.1f}% "
        f"(n={strat.overall.n_pdbs}); "
        f"best family '{best.family}' {100 * best.pass_rate:.1f}% "
        f"(n={best.n_pdbs}); "
        f"worst family '{worst.family}' {100 * worst.pass_rate:.1f}% "
        f"(n={worst.n_pdbs})"
    )
