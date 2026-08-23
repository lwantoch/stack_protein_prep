#!/usr/bin/env python3
"""Plot FRUTON vs PDB-REDO per-metric comparison (Nature R2).

Reviewer figure for Nature R2: how does FRUTON's output compare
to the same PDB entry's PDB-REDO re-refinement?

Usage:
    python plot_pdb_redo_comparison.py \\
        --fruton-dir metrics/mmbsa200_fruton/ \\
        --pdb-redo-dir metrics/mmbsa200_pdb_redo/ \\
        --outdir figures/nature_r2_pdb_redo/

Each directory holds <PDB>.json files with the metric fields listed
in _pdb_redo_compare.METRIC_NAMES.

Produces:
    pdb_redo_comparison.png     grouped bar plot: n_fruton_better vs
                                n_pdb_redo_better per metric, tie in
                                grey, unavailable in light grey.
    pdb_redo_comparison.md      reviewer-facing markdown table.
    pdb_redo_comparison.json    machine-readable full report.

Pure Python + matplotlib.  License-free.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

from stack_protein_preparation._pdb_redo_compare import (
    METRIC_NAMES,
    build_comparison,
    format_comparison_table,
)


def _parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--fruton-dir", type=Path, required=True)
    p.add_argument("--pdb-redo-dir", type=Path, required=True)
    p.add_argument("--outdir", type=Path, required=True)
    p.add_argument(
        "--metric",
        action="append",
        default=None,
        help="Restrict to specific metric keys (repeatable). "
             "Default: every metric in METRIC_NAMES.",
    )
    p.add_argument(
        "--title",
        default="FRUTON vs PDB-REDO — MMBSA_200 bench",
    )
    return p.parse_args(argv)


def _plot(comp, outpath: Path, title: str) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    outpath.parent.mkdir(parents=True, exist_ok=True)

    labels = [m.reviewer_label for m in comp.per_metric]
    fruton_wins = [m.n_fruton_better for m in comp.per_metric]
    redo_wins = [m.n_pdb_redo_better for m in comp.per_metric]
    ties = [m.n_tie for m in comp.per_metric]
    unavailable = [m.n_unavailable for m in comp.per_metric]

    if not labels:
        # Nothing to plot; produce a minimal figure so the reviewer
        # sees "no matched pdbs" instead of a missing file.
        fig, ax = plt.subplots(figsize=(6, 3))
        ax.text(0.5, 0.5, "no matched PDBs", ha="center", va="center",
                transform=ax.transAxes)
        ax.axis("off")
        fig.savefig(outpath, dpi=150)
        plt.close(fig)
        return

    fig, ax = plt.subplots(figsize=(max(9.0, 0.8 * len(labels) + 3), 5.0))
    x_positions = list(range(len(labels)))
    bar_h = 0.6

    # Stack ordering: FRUTON-better (blue), PDB-REDO-better (orange), tie (grey), unavailable (light grey).
    bottom = [0] * len(labels)
    ax.barh(x_positions, fruton_wins, height=bar_h, left=bottom,
            color="#2b6cb0", edgecolor="none", label="FRUTON better")
    bottom = [b + v for b, v in zip(bottom, fruton_wins)]
    ax.barh(x_positions, redo_wins, height=bar_h, left=bottom,
            color="#c05621", edgecolor="none", label="PDB-REDO better")
    bottom = [b + v for b, v in zip(bottom, redo_wins)]
    ax.barh(x_positions, ties, height=bar_h, left=bottom,
            color="#a0aec0", edgecolor="none", label="tie")
    bottom = [b + v for b, v in zip(bottom, ties)]
    ax.barh(x_positions, unavailable, height=bar_h, left=bottom,
            color="#e2e8f0", edgecolor="none", label="unavailable")

    ax.set_yticks(x_positions)
    ax.set_yticklabels(labels)
    ax.invert_yaxis()
    ax.set_xlabel("PDBs")
    ax.set_title(title)
    ax.legend(loc="lower right", fontsize=9)

    fig.tight_layout()
    fig.savefig(outpath, dpi=150)
    plt.close(fig)


def main(argv: list[str] | None = None) -> int:
    ns = _parse_args(argv)
    metrics = ns.metric or list(METRIC_NAMES.keys())
    ns.outdir.mkdir(parents=True, exist_ok=True)
    try:
        comp = build_comparison(ns.fruton_dir, ns.pdb_redo_dir, metrics=metrics)
    except (FileNotFoundError, ValueError) as exc:
        print(f"[pdb_redo] ERROR: {exc}", file=sys.stderr)
        return 2

    (ns.outdir / "pdb_redo_comparison.json").write_text(
        json.dumps(comp.to_dict(), indent=2)
    )
    (ns.outdir / "pdb_redo_comparison.md").write_text(
        format_comparison_table(comp)
    )
    _plot(comp, ns.outdir / "pdb_redo_comparison.png", title=ns.title)

    n_matched = comp.per_metric[0].n_matched if comp.per_metric else 0
    print(
        f"[pdb_redo] matched {n_matched} PDBs across {len(comp.per_metric)} "
        f"metrics -> {ns.outdir / 'pdb_redo_comparison.png'}"
    )
    for m in comp.per_metric:
        print(
            f"[pdb_redo]   {m.reviewer_label:<32}: "
            f"FRUTON {m.n_fruton_better:>3}  PDB-REDO {m.n_pdb_redo_better:>3}  "
            f"tie {m.n_tie:>3}  unavail {m.n_unavailable:>3}"
        )
    return 0


if __name__ == "__main__":
    sys.exit(main())
