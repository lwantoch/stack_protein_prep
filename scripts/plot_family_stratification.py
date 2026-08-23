#!/usr/bin/env python3
"""Plot FRUTON pass rate + rescued residues stratified by protein family.

Reviewer figure for Nature R1: pipeline-wide aggregate hides class
imbalance.  Break the bench down by kinase / GPCR / protease /
metalloenzyme / … and show pass rate + residues rescued per family.

Usage:
    python plot_family_stratification.py \\
        --bench-json runs/mmbsa200_baseline.json \\
        --family-map metrics/family_by_pdb.json \\
        --outdir figures/nature_r1_family/

Produces:
    family_stratification.png     grouped bar (pass rate + rescued
                                  per family, overall reference line)
    family_stratification.md      markdown table (one row per family
                                  + overall)
    family_stratification.json    full report

Pure Python + matplotlib.  License-free.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

from stack_protein_preparation._family_stratification import (
    format_stratification_table,
    imbalance_summary,
    load_family_mapping,
    stratify_bench_by_family,
)


def _parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--bench-json", type=Path, required=True,
                   help="Bench-results JSON (list of per-protein dicts).")
    p.add_argument("--family-map", type=Path, required=True,
                   help="JSON {pdb_id: family_label} mapping.")
    p.add_argument("--outdir", type=Path, required=True)
    p.add_argument(
        "--title",
        default="FRUTON bench, stratified by family — MMBSA_200",
    )
    return p.parse_args(argv)


def _plot(strat, outpath: Path, title: str) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    outpath.parent.mkdir(parents=True, exist_ok=True)

    families = [a.family for a in strat.per_family]
    n_pdbs = [a.n_pdbs for a in strat.per_family]
    pass_pct = [100 * a.pass_rate for a in strat.per_family]
    rescued = [a.total_delta_n for a in strat.per_family]
    labels = [f"{f}\nn={n}" for f, n in zip(families, n_pdbs)]

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(max(9.0, 1.0 * len(families) + 3), 4.5))

    # Left: pass rate per family + overall reference
    ax1.bar(labels, pass_pct, color="#2b6cb0", edgecolor="none")
    ax1.axhline(100 * strat.overall.pass_rate, color="black",
                linestyle="--", linewidth=0.9,
                label=f"overall {100 * strat.overall.pass_rate:.1f}%")
    ax1.set_ylabel("gate PASS rate (%)")
    ax1.set_ylim(0, 100)
    ax1.set_title("Pass rate per family")
    ax1.legend(fontsize=8, loc="lower right")
    for tick in ax1.get_xticklabels():
        tick.set_rotation(30)
        tick.set_horizontalalignment("right")

    # Right: total residues rescued per family
    ax2.bar(labels, rescued, color="#c05621", edgecolor="none")
    ax2.set_ylabel("total residues rescued (Σ ΔN)")
    ax2.set_title("Residues rescued per family")
    for tick in ax2.get_xticklabels():
        tick.set_rotation(30)
        tick.set_horizontalalignment("right")

    fig.suptitle(title)
    fig.tight_layout()
    fig.savefig(outpath, dpi=150)
    plt.close(fig)


def main(argv: list[str] | None = None) -> int:
    ns = _parse_args(argv)
    if not ns.bench_json.is_file():
        print(f"[family_strat] ERROR: bench json not found: {ns.bench_json}",
              file=sys.stderr)
        return 2
    if not ns.family_map.is_file():
        print(f"[family_strat] ERROR: family map not found: {ns.family_map}",
              file=sys.stderr)
        return 2

    try:
        with ns.bench_json.open() as fh:
            bench = json.load(fh)
        if not isinstance(bench, list):
            raise ValueError(f"{ns.bench_json} is not a list")
        family_map = load_family_mapping(ns.family_map)
    except (ValueError, json.JSONDecodeError) as exc:
        print(f"[family_strat] ERROR: {exc}", file=sys.stderr)
        return 2

    ns.outdir.mkdir(parents=True, exist_ok=True)
    strat = stratify_bench_by_family(bench, family_map)

    (ns.outdir / "family_stratification.json").write_text(
        json.dumps(strat.to_dict(), indent=2)
    )
    (ns.outdir / "family_stratification.md").write_text(
        format_stratification_table(strat)
    )
    _plot(strat, ns.outdir / "family_stratification.png", title=ns.title)

    print(f"[family_strat] {imbalance_summary(strat)}")
    print(f"[family_strat] -> {ns.outdir / 'family_stratification.png'}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
