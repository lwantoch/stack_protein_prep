"""Family-stratified bench analysis (Nature R1 concern).

Reviewer perspective: "your bench-wide pass rate is X%, but which
protein families make that up?  Kinases are trivial to fill — how do
you do on multi-domain transferases, GPCRs, or intrinsically-flexible
scaffolds?"

Given a FRUTON bench-results JSON + a family-mapping CSV
(``pdb_id,family``), this module partitions the bench into per-family
buckets and computes the same aggregate stats as _ablation.compute_variant_aggregate,
one row per family.

Pure stdlib.  License-free.
"""
from __future__ import annotations

import csv
import statistics
from dataclasses import dataclass, field
from pathlib import Path


UNKNOWN_FAMILY = "unassigned"


@dataclass
class FamilyAggregate:
    family: str
    n_pdbs: int
    n_passed: int
    pass_rate: float
    mean_clash: float
    mean_brk: float
    total_delta_n: int
    mean_delta_n: float
    member_pdbs: list[str] = field(default_factory=list)

    def to_dict(self) -> dict:
        return {
            "family": self.family,
            "n_pdbs": self.n_pdbs,
            "n_passed": self.n_passed,
            "pass_rate": round(self.pass_rate, 4),
            "mean_clash": round(self.mean_clash, 3),
            "mean_brk": round(self.mean_brk, 3),
            "total_delta_n": self.total_delta_n,
            "mean_delta_n": round(self.mean_delta_n, 3),
            "member_pdbs": list(self.member_pdbs),
        }


def load_family_mapping(path: str | Path) -> dict[str, str]:
    """Read a CSV with header ``pdb_id,family`` into a dict.

    PDB IDs are upper-cased so lookups work regardless of the bench
    JSON casing.  Rows with missing family field are skipped silently.
    """
    path = Path(path)
    mapping: dict[str, str] = {}
    with path.open(newline="") as fh:
        reader = csv.DictReader(fh)
        if reader.fieldnames is None:
            return mapping
        # Accept flexible column names: {pdb, pdb_id} × {family, class, group}.
        pdb_col = next(
            (c for c in ("pdb_id", "pdb") if c in reader.fieldnames), None,
        )
        fam_col = next(
            (c for c in ("family", "class", "group") if c in reader.fieldnames),
            None,
        )
        if pdb_col is None or fam_col is None:
            raise ValueError(
                f"family mapping {path}: expected columns pdb_id,family "
                f"(or pdb,class/group); got {reader.fieldnames}"
            )
        for row in reader:
            pdb = (row.get(pdb_col) or "").strip().upper()
            fam = (row.get(fam_col) or "").strip()
            if not pdb or not fam:
                continue
            mapping[pdb] = fam
    return mapping


def _safe_mean(values: list[float]) -> float:
    return statistics.fmean(values) if values else 0.0


def stratify_bench_by_family(
    bench_results: list[dict],
    family_mapping: dict[str, str],
    include_unassigned: bool = True,
) -> list[FamilyAggregate]:
    """Group per-protein records by family and aggregate.

    Args:
        bench_results: list of per-protein dicts as produced by
            scripts/CESGA_SLURM/fruton_bench_mmbsa200_full.py.
        family_mapping: ``{pdb_id_upper: family_name}`` from
            ``load_family_mapping``.
        include_unassigned: if True, PDBs not in the mapping are
            bucketed as ``UNKNOWN_FAMILY``; if False, they are dropped.

    Returns list of FamilyAggregate, sorted by family name (with
    ``unassigned`` last if present).
    """
    by_family: dict[str, list[dict]] = {}
    for record in bench_results:
        pdb = str(record.get("pdb") or "").upper()
        if not pdb:
            continue
        fam = family_mapping.get(pdb)
        if fam is None:
            if not include_unassigned:
                continue
            fam = UNKNOWN_FAMILY
        by_family.setdefault(fam, []).append(record)

    out: list[FamilyAggregate] = []
    for fam, records in by_family.items():
        clashes = [float(r.get("clash", 0)) for r in records]
        brks = [float(r.get("brk", 0)) for r in records]
        delta_ns = [int(r.get("delta_n", 0)) for r in records]
        n_passed = sum(1 for r in records if r.get("gate_pass"))
        pdb_ids = [str(r.get("pdb", "")).upper() for r in records]
        out.append(FamilyAggregate(
            family=fam,
            n_pdbs=len(records),
            n_passed=n_passed,
            pass_rate=(n_passed / len(records)) if records else 0.0,
            mean_clash=_safe_mean(clashes),
            mean_brk=_safe_mean(brks),
            total_delta_n=sum(delta_ns),
            mean_delta_n=_safe_mean(delta_ns),
            member_pdbs=sorted(pdb_ids),
        ))

    def _sort_key(fa: FamilyAggregate) -> tuple[int, str]:
        # Unassigned always last.
        return (1 if fa.family == UNKNOWN_FAMILY else 0, fa.family.lower())

    out.sort(key=_sort_key)
    return out


def format_family_table(aggregates: list[FamilyAggregate]) -> str:
    """Reviewer-facing markdown: one row per family."""
    if not aggregates:
        return "(no family aggregates)\n"
    lines = [
        "| Family | n_pdbs | pass % | rescued | ⟨clash⟩ | ⟨brk⟩ | ⟨ΔN⟩ |",
        "|---|---:|---:|---:|---:|---:|---:|",
    ]
    for a in aggregates:
        lines.append(
            f"| {a.family} | {a.n_pdbs} | {100 * a.pass_rate:.1f}% | "
            f"{a.total_delta_n} | {a.mean_clash:.2f} | "
            f"{a.mean_brk:.2f} | {a.mean_delta_n:.2f} |"
        )
    return "\n".join(lines) + "\n"


def hardest_families(
    aggregates: list[FamilyAggregate],
    min_family_size: int = 3,
    top_n: int = 5,
) -> list[FamilyAggregate]:
    """Return the ``top_n`` families with lowest pass_rate (min-size gate).

    A family with fewer than ``min_family_size`` members is excluded
    so single-outlier families don't dominate the "hardest" list.
    """
    eligible = [a for a in aggregates
                if a.n_pdbs >= min_family_size and a.family != UNKNOWN_FAMILY]
    return sorted(eligible, key=lambda a: a.pass_rate)[:top_n]
