"""Ablation-study harness for FRUTON pipeline (Nature R4 / JCTC R2).

Reviewer concern: FRUTON stacks multiple gates (pLDDT, IDR, clash,
chirality, ω-planarity, MolProbity-style reviewer gate, adaptive
rollback).  How much does each individual gate contribute to the
bench-wide pass rate + rescue count?  Ablation = disable one gate at a
time, rerun the bench, tabulate deltas.

This module is the *analysis* layer — it does not re-execute the
pipeline (that requires SLURM + hours of GPU time).  Instead it:

1. Names the ablate-able gates canonically (GATE_NAMES).
2. Loads a directory of per-variant bench JSONs
   (baseline.json, no_plddt.json, no_idr.json, ...) each carrying
   the same per-PDB result shape produced by
   scripts/CESGA_SLURM/fruton_bench_mmbsa200_full.py.
3. Computes per-variant aggregates + per-PDB delta-vs-baseline.
4. Emits a comparison table (list of dicts) for CSV / matplotlib.

Pure Python + stdlib.  License-free.
"""
from __future__ import annotations

import json
import statistics
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable


# Canonical list of gates that can be toggled off in an ablation run.
# Each name maps to a plain-language description a reviewer would
# recognise; the CLI runner (future) is expected to accept these as
# --disable-gate flags.
GATE_NAMES: dict[str, str] = {
    "plddt": "AF pLDDT window-mean gate (rejects fills where per-residue AF "
             "confidence is below threshold)",
    "idr":   "UniProt / MobiDB IDR gate (rejects fills whose gap window "
             "overlaps an annotated intrinsically-disordered region)",
    "clash": "post-splice clash gate (rejects fills that introduce serious "
             "non-bonded overlaps)",
    "chirality_d": "MODELLER LoopModel D-chirality reject (drops conformers "
             "whose Cα improper dihedral flips to D)",
    "chirality_improper": "improper-dihedral chirality gate (post-fill "
             "check on |χ_N-CA-C-CB| ≈ 120°)",
    "omega": "ω peptide-bond planarity gate (rejects new cis-nonPro / "
             "non-planar peptide bonds)",
    "rollback": "adaptive rollback (un-fittable gaps to REMARK 465 instead "
             "of broken bonds)",
    "loop_refine": "MODELLER LoopModel refinement (adaptive fast → slow)",
    "reviewer_gate": "MolProbity-style reviewer gate (Ramachandran outliers, "
             "clashscore, ω non-planar count)",
}


@dataclass
class VariantAggregate:
    """Per-variant statistics computed across a bench JSON."""
    variant: str
    n_pdbs: int
    n_passed: int
    pass_rate: float
    mean_clash: float
    mean_brk: float
    total_delta_n: int
    mean_delta_n: float
    mean_refine_seconds: float
    disabled_gate: str | None = None

    def to_dict(self) -> dict:
        return {
            "variant": self.variant,
            "disabled_gate": self.disabled_gate,
            "n_pdbs": self.n_pdbs,
            "n_passed": self.n_passed,
            "pass_rate": round(self.pass_rate, 4),
            "mean_clash": round(self.mean_clash, 3),
            "mean_brk": round(self.mean_brk, 3),
            "total_delta_n": self.total_delta_n,
            "mean_delta_n": round(self.mean_delta_n, 3),
            "mean_refine_seconds": round(self.mean_refine_seconds, 3),
        }


@dataclass
class PerPdbDelta:
    """One protein's before/after delta between baseline and a variant."""
    pdb: str
    baseline_gate_pass: bool
    variant_gate_pass: bool
    baseline_delta_n: int
    variant_delta_n: int
    baseline_clash: int
    variant_clash: int

    def sign(self) -> str:
        """+1 rescued (baseline FAIL, variant PASS), -1 lost, 0 unchanged."""
        if self.variant_gate_pass and not self.baseline_gate_pass:
            return "rescued"
        if self.baseline_gate_pass and not self.variant_gate_pass:
            return "lost"
        return "unchanged"


@dataclass
class AblationComparison:
    baseline: VariantAggregate
    variants: list[VariantAggregate]
    per_pdb_deltas: dict[str, list[PerPdbDelta]] = field(default_factory=dict)

    def to_dict(self) -> dict:
        return {
            "baseline": self.baseline.to_dict(),
            "variants": [v.to_dict() for v in self.variants],
            "per_pdb_deltas_by_variant": {
                vname: [
                    {
                        "pdb": d.pdb,
                        "baseline_gate_pass": d.baseline_gate_pass,
                        "variant_gate_pass": d.variant_gate_pass,
                        "baseline_delta_n": d.baseline_delta_n,
                        "variant_delta_n": d.variant_delta_n,
                        "baseline_clash": d.baseline_clash,
                        "variant_clash": d.variant_clash,
                        "sign": d.sign(),
                    }
                    for d in deltas
                ]
                for vname, deltas in self.per_pdb_deltas.items()
            },
        }


def load_bench_json(path: str | Path) -> list[dict]:
    """Load a bench-results JSON (list of per-protein dicts) or raise."""
    path = Path(path)
    with path.open() as fh:
        data = json.load(fh)
    if not isinstance(data, list):
        raise ValueError(f"{path} is not a list of per-protein records")
    return data


def _safe_mean(values: Iterable[float]) -> float:
    values = list(values)
    return statistics.fmean(values) if values else 0.0


def compute_variant_aggregate(
    variant_name: str,
    results: list[dict],
    disabled_gate: str | None = None,
) -> VariantAggregate:
    """Roll a per-protein result list into a single aggregate row."""
    n = len(results)
    n_passed = sum(1 for r in results if r.get("gate_pass"))
    return VariantAggregate(
        variant=variant_name,
        disabled_gate=disabled_gate,
        n_pdbs=n,
        n_passed=n_passed,
        pass_rate=(n_passed / n) if n else 0.0,
        mean_clash=_safe_mean(r.get("clash", 0) for r in results),
        mean_brk=_safe_mean(r.get("brk", 0) for r in results),
        total_delta_n=sum(int(r.get("delta_n", 0)) for r in results),
        mean_delta_n=_safe_mean(int(r.get("delta_n", 0)) for r in results),
        mean_refine_seconds=_safe_mean(
            float(r.get("refine_seconds", 0.0)) for r in results
        ),
    )


def diff_per_pdb(
    baseline: list[dict],
    variant: list[dict],
) -> list[PerPdbDelta]:
    """Emit one PerPdbDelta per PDB that appears in both lists."""
    by_id = {r["pdb"]: r for r in baseline if "pdb" in r}
    out: list[PerPdbDelta] = []
    for r in variant:
        pid = r.get("pdb")
        if pid is None or pid not in by_id:
            continue
        b = by_id[pid]
        out.append(PerPdbDelta(
            pdb=pid,
            baseline_gate_pass=bool(b.get("gate_pass")),
            variant_gate_pass=bool(r.get("gate_pass")),
            baseline_delta_n=int(b.get("delta_n", 0)),
            variant_delta_n=int(r.get("delta_n", 0)),
            baseline_clash=int(b.get("clash", 0)),
            variant_clash=int(r.get("clash", 0)),
        ))
    return out


def build_comparison(
    baseline_json: str | Path,
    variant_jsons: dict[str, str | Path],
    disabled_gate_by_variant: dict[str, str] | None = None,
) -> AblationComparison:
    """Assemble an AblationComparison from a baseline JSON + variant JSONs.

    Args:
        baseline_json: path to the reference bench-results JSON.
        variant_jsons: mapping ``{variant_label: json_path}``.
        disabled_gate_by_variant: optional mapping labelling which gate is
            disabled in each variant (keys of ``GATE_NAMES``).
    """
    baseline_results = load_bench_json(baseline_json)
    baseline_agg = compute_variant_aggregate("baseline", baseline_results)

    variants: list[VariantAggregate] = []
    per_pdb: dict[str, list[PerPdbDelta]] = {}
    disabled_gate_by_variant = disabled_gate_by_variant or {}
    for name, path in variant_jsons.items():
        var_results = load_bench_json(path)
        gate = disabled_gate_by_variant.get(name)
        if gate is not None and gate not in GATE_NAMES:
            raise ValueError(
                f"variant {name!r}: disabled_gate={gate!r} not in GATE_NAMES "
                f"(known: {sorted(GATE_NAMES)})"
            )
        variants.append(compute_variant_aggregate(name, var_results, disabled_gate=gate))
        per_pdb[name] = diff_per_pdb(baseline_results, var_results)

    return AblationComparison(
        baseline=baseline_agg,
        variants=variants,
        per_pdb_deltas=per_pdb,
    )


def format_comparison_table(comp: AblationComparison) -> str:
    """Reviewer-facing markdown table: one row per variant, columns for
    pass_rate, total rescued residues, mean clash, mean broken bonds.
    """
    rows = [comp.baseline, *comp.variants]
    lines: list[str] = []
    lines.append("| Variant | Disabled gate | n_pdbs | pass % | rescued | ⟨clash⟩ | ⟨brk⟩ |")
    lines.append("|---|---|---:|---:|---:|---:|---:|")
    for a in rows:
        gate_label = a.disabled_gate or "-"
        lines.append(
            f"| {a.variant} | {gate_label} | {a.n_pdbs} | "
            f"{100 * a.pass_rate:.1f}% | {a.total_delta_n} | "
            f"{a.mean_clash:.2f} | {a.mean_brk:.2f} |"
        )
    return "\n".join(lines) + "\n"


def rescued_and_lost_counts(comp: AblationComparison) -> dict[str, tuple[int, int]]:
    """Per variant: how many proteins moved FAIL→PASS ('rescued') and PASS→FAIL ('lost')."""
    out: dict[str, tuple[int, int]] = {}
    for vname, deltas in comp.per_pdb_deltas.items():
        rescued = sum(1 for d in deltas if d.sign() == "rescued")
        lost = sum(1 for d in deltas if d.sign() == "lost")
        out[vname] = (rescued, lost)
    return out
