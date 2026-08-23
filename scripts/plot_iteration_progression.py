#!/usr/bin/env python3
"""Plot FRUTON iter-over-iter gate PASS % and non-planar ω trend.

Killer publication figure for the paper introduction: at a glance a
reviewer sees the 89.6 % → 91.3 % → 97.8 % → 97.9 % gate-PASS
progression across the iter-2 → iter-3 → iter-4 → iter-5-mini
hardening loop, plus the non-planar ω count halving from 6 → 3.

Usage:
    python plot_iteration_progression.py \\
        --history iteration_history.json \\
        --outdir figures/paper_intro/

``iteration_history.json`` schema:
    [
        {"iter_label": "iter-2", "n_pdbs": 48, "gate_pass_pct": 89.6,
         "n_omega_non_planar": null, "notes": "43/48 from BASELINE_STATS.md"},
        {"iter_label": "iter-3", "n_pdbs": 46, "gate_pass_pct": 91.3,
         "n_omega_non_planar": 6, "notes": "bench_20260823_1513"},
        ...
    ]

``n_omega_non_planar`` may be ``null`` for iterations where we do not
have the ω figure data.  The bench_output/*/PROVENANCE.md files are
the ground-truth source for these numbers; this script trusts the
input JSON and does not re-parse them.

Produces:
    iteration_progression.png    two-panel matplotlib line chart:
                                 left  = gate PASS % vs iter, with
                                         annotations per point.
                                 right = non-planar ω count vs iter,
                                         log-scale y for the trend.
    iteration_progression.csv    machine-readable per-iter roll-up.

Pure Python + matplotlib.  License-free.
"""
from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path


def _parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--history", type=Path, required=True,
                   help="iteration_history.json — list of per-iter dicts.")
    p.add_argument("--outdir", type=Path, required=True)
    p.add_argument(
        "--title",
        default="FRUTON hardening progression — MMBSA_200 bench",
    )
    return p.parse_args(argv)


def load_history(path: Path) -> list[dict]:
    with path.open() as fh:
        data = json.load(fh)
    if not isinstance(data, list):
        raise ValueError(f"{path} is not a JSON list")
    for i, entry in enumerate(data):
        if "iter_label" not in entry:
            raise ValueError(f"entry {i}: missing iter_label")
        if "gate_pass_pct" not in entry:
            raise ValueError(f"entry {i}: missing gate_pass_pct")
    return data


def _write_csv(history: list[dict], outpath: Path) -> None:
    outpath.parent.mkdir(parents=True, exist_ok=True)
    fields = ["iter_label", "n_pdbs", "gate_pass_pct",
              "n_omega_non_planar", "notes"]
    with outpath.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, extrasaction="ignore")
        w.writeheader()
        for entry in history:
            w.writerow({k: entry.get(k, "") for k in fields})


def _plot(history: list[dict], outpath: Path, title: str) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    outpath.parent.mkdir(parents=True, exist_ok=True)

    labels = [str(e["iter_label"]) for e in history]
    x = list(range(len(labels)))
    pass_pcts = [float(e["gate_pass_pct"]) for e in history]

    # non-planar ω may be null for some iterations; keep parallel lists
    # and drop nulls in the right-panel plotting.
    omega_x: list[int] = []
    omega_y: list[int] = []
    for i, e in enumerate(history):
        v = e.get("n_omega_non_planar")
        if v is not None:
            omega_x.append(i)
            omega_y.append(int(v))

    fig, (ax_l, ax_r) = plt.subplots(1, 2, figsize=(11.5, 4.5))

    # Left panel: gate PASS % over iterations
    ax_l.plot(x, pass_pcts, "-o", color="#2b6cb0", linewidth=2.0, markersize=8)
    for xi, pct, lbl in zip(x, pass_pcts, labels):
        ax_l.annotate(
            f"{pct:.1f} %",
            xy=(xi, pct),
            xytext=(0, 8),
            textcoords="offset points",
            ha="center",
            fontsize=9,
            fontweight="bold",
        )
    ax_l.set_xticks(x)
    ax_l.set_xticklabels(labels)
    ax_l.set_ylabel("Gate PASS rate (%)")
    ax_l.set_ylim(min(pass_pcts) - 5, 101)
    ax_l.set_title("Gate PASS rate per iteration")
    ax_l.axhline(100, color="grey", linestyle=":", linewidth=0.8)
    ax_l.grid(True, axis="y", linestyle=":", alpha=0.5)

    # Right panel: non-planar ω count over iterations (log-y)
    if omega_x:
        ax_r.plot(
            omega_x, omega_y, "-o",
            color="#c05621", linewidth=2.0, markersize=8,
        )
        for xi, v in zip(omega_x, omega_y):
            ax_r.annotate(
                str(v), xy=(xi, max(v, 1)),
                xytext=(0, 8), textcoords="offset points",
                ha="center", fontsize=9, fontweight="bold",
            )
        ax_r.set_yscale("symlog", linthresh=1)
        ax_r.set_ylim(0, max(max(omega_y) * 2, 10))
    else:
        ax_r.text(
            0.5, 0.5, "no ω-planarity data yet",
            transform=ax_r.transAxes, ha="center", va="center",
        )

    ax_r.set_xticks(x)
    ax_r.set_xticklabels(labels)
    ax_r.set_ylabel("Non-planar ω peptide bonds (count)")
    ax_r.set_title("Non-planar ω per iteration (log-y)")
    ax_r.grid(True, axis="y", linestyle=":", alpha=0.5)

    fig.suptitle(title)
    fig.tight_layout()
    fig.savefig(outpath, dpi=150)
    plt.close(fig)


def main(argv: list[str] | None = None) -> int:
    ns = _parse_args(argv)
    if not ns.history.is_file():
        print(f"[iter_progress] ERROR: history not found: {ns.history}",
              file=sys.stderr)
        return 2

    try:
        history = load_history(ns.history)
    except (ValueError, json.JSONDecodeError) as exc:
        print(f"[iter_progress] ERROR: {exc}", file=sys.stderr)
        return 2

    if not history:
        print("[iter_progress] history is empty; nothing to plot",
              file=sys.stderr)
        return 2

    ns.outdir.mkdir(parents=True, exist_ok=True)
    _write_csv(history, ns.outdir / "iteration_progression.csv")
    _plot(history, ns.outdir / "iteration_progression.png", title=ns.title)

    first = history[0]
    last = history[-1]
    print(
        f"[iter_progress] {len(history)} iterations plotted "
        f"({first['iter_label']} {first['gate_pass_pct']:.1f}% -> "
        f"{last['iter_label']} {last['gate_pass_pct']:.1f}%) "
        f"-> {ns.outdir / 'iteration_progression.png'}"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
