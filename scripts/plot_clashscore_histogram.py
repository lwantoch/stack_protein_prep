#!/usr/bin/env python3
"""Plot bench-wide clashscore distribution (JCTC R2).

Reviewer figure: histogram of clash-metric distribution across the
MMBSA_200 bench.  Reads FRUTON's per-protein quality-report JSONs
(the ``n_clash_pairs`` and ``n_vdw_clashes`` fields emitted by
_filler_quality_check.check_model_quality) either directly from a
bench-results list or from a directory of per-PDB JSONs.

Usage (bench list):
    python plot_clashscore_histogram.py \\
        --bench-json runs/mmbsa200_baseline.json \\
        --outdir figures/jctc_r2_clash/

Usage (per-PDB directory):
    python plot_clashscore_histogram.py \\
        --pdb-metrics-dir metrics/mmbsa200_fruton/ \\
        --outdir figures/jctc_r2_clash/

Produces:
    clashscore_histogram.png   two-panel histogram (n_clash_pairs +
                               n_vdw_clashes), with vertical
                               reference lines at 4 and 20 (the
                               MolProbity 90th/99th percentile
                               benchmarks per Chen et al. 2010).
    clashscore_histogram.csv   per-PDB (pdb_id, n_clash_pairs,
                               n_vdw_clashes).
    summary.txt                one-liner + percentiles.

Pure Python + matplotlib.  License-free.
"""
from __future__ import annotations

import argparse
import csv
import json
import statistics
import sys
from pathlib import Path


def _parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    src = p.add_mutually_exclusive_group(required=True)
    src.add_argument(
        "--bench-json", type=Path,
        help="Bench-results JSON (list of per-PDB dicts).",
    )
    src.add_argument(
        "--pdb-metrics-dir", type=Path,
        help="Directory of <PDB>.json quality-report files.",
    )
    p.add_argument("--outdir", type=Path, required=True)
    p.add_argument(
        "--bin-width",
        type=int,
        default=1,
        help="Histogram bin width (default 1).",
    )
    p.add_argument(
        "--x-max",
        type=int,
        default=100,
        help="Histogram x-axis upper bound (default 100).",
    )
    p.add_argument(
        "--log-y",
        action="store_true",
        help="Log-scale y-axis (highlights tail).",
    )
    p.add_argument(
        "--title",
        default="Clashscore distribution — MMBSA_200 bench",
    )
    return p.parse_args(argv)


def _load_from_bench_json(path: Path) -> list[dict]:
    with path.open() as fh:
        data = json.load(fh)
    if not isinstance(data, list):
        raise ValueError(f"{path} is not a JSON list of per-PDB records")
    return data


def _load_from_metrics_dir(path: Path) -> list[dict]:
    if not path.is_dir():
        raise FileNotFoundError(f"metrics dir not found: {path}")
    out: list[dict] = []
    for p in sorted(path.glob("*.json")):
        with p.open() as fh:
            payload = json.load(fh)
        if not isinstance(payload, dict):
            continue
        # Backfill pdb_id from file stem if missing.
        payload.setdefault("pdb", p.stem.upper())
        out.append(payload)
    return out


def _to_int(x, default: int | None = None) -> int | None:
    if x is None:
        return default
    try:
        return int(x)
    except (TypeError, ValueError):
        return default


def _percentile(sorted_values: list[int], pct: float) -> float:
    if not sorted_values:
        return 0.0
    if pct <= 0:
        return float(sorted_values[0])
    if pct >= 100:
        return float(sorted_values[-1])
    k = (len(sorted_values) - 1) * pct / 100.0
    lo = int(k)
    hi = min(lo + 1, len(sorted_values) - 1)
    frac = k - lo
    return sorted_values[lo] * (1 - frac) + sorted_values[hi] * frac


def _write_csv(rows: list[dict], outpath: Path) -> None:
    outpath.parent.mkdir(parents=True, exist_ok=True)
    fields = ["pdb_id", "n_clash_pairs", "n_vdw_clashes"]
    with outpath.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields)
        w.writeheader()
        w.writerows(rows)


def _plot(
    n_clash_pairs: list[int],
    n_vdw_clashes: list[int],
    outpath: Path,
    bin_width: int,
    x_max: int,
    log_y: bool,
    title: str,
) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    outpath.parent.mkdir(parents=True, exist_ok=True)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11.5, 4.5))
    bins = list(range(0, x_max + bin_width, bin_width))

    for ax, values, ylabel, colour, header in (
        (ax1, n_clash_pairs, "PDBs", "steelblue", "n_clash_pairs (all)"),
        (ax2, n_vdw_clashes, "PDBs", "darkorange", "n_vdw_clashes (MolProbity-style)"),
    ):
        clipped = [min(v, x_max) for v in values]
        ax.hist(clipped, bins=bins, color=colour, edgecolor="none")
        if log_y:
            ax.set_yscale("log")
        # Reference lines: MolProbity 90th ~= 4, 99th ~= 20 (Chen et al. 2010).
        ax.axvline(4, color="grey", linewidth=0.9, linestyle=":",
                   label="MolProbity 90th ≈ 4")
        ax.axvline(20, color="grey", linewidth=0.9, linestyle=":",
                   label="MolProbity 99th ≈ 20")
        ax.set_xlim(0, x_max)
        ax.set_xlabel(header)
        ax.set_ylabel(ylabel)
        ax.set_title(header)
        if values:
            s_vals = sorted(values)
            p50 = _percentile(s_vals, 50)
            p90 = _percentile(s_vals, 90)
            p99 = _percentile(s_vals, 99)
            mn = statistics.fmean(values)
            n_within_4 = sum(1 for v in values if v <= 4)
            legend = (
                f"n = {len(values)}\n"
                f"mean = {mn:.2f}\n"
                f"p50 = {p50:.1f}\n"
                f"p90 = {p90:.1f}\n"
                f"p99 = {p99:.1f}\n"
                f"≤ 4 : {n_within_4} ({100 * n_within_4 / len(values):.1f}%)"
            )
            ax.text(
                0.97, 0.97, legend,
                transform=ax.transAxes,
                fontsize=8, family="monospace",
                verticalalignment="top", horizontalalignment="right",
                bbox={"boxstyle": "round,pad=0.4",
                      "facecolor": "white", "alpha": 0.85},
            )
        ax.legend(loc="lower right", fontsize=7)

    fig.suptitle(title)
    fig.tight_layout()
    fig.savefig(outpath, dpi=150)
    plt.close(fig)


def build_summary(
    n_clash_pairs: list[int], n_vdw_clashes: list[int]
) -> str:
    lines: list[str] = []
    for label, values in (
        ("n_clash_pairs", n_clash_pairs),
        ("n_vdw_clashes", n_vdw_clashes),
    ):
        if not values:
            lines.append(f"{label}: no data")
            continue
        s = sorted(values)
        lines.append(
            f"{label}: n={len(values)} "
            f"mean={statistics.fmean(values):.2f} "
            f"p50={_percentile(s, 50):.1f} "
            f"p90={_percentile(s, 90):.1f} "
            f"p99={_percentile(s, 99):.1f} "
            f"max={s[-1]} "
            f"(≤ 4 in {sum(1 for v in values if v <= 4)}/{len(values)})"
        )
    return "\n".join(lines) + "\n"


def main(argv: list[str] | None = None) -> int:
    ns = _parse_args(argv)
    try:
        if ns.bench_json is not None:
            if not ns.bench_json.is_file():
                print(f"[clash] ERROR: bench not found: {ns.bench_json}", file=sys.stderr)
                return 2
            records = _load_from_bench_json(ns.bench_json)
        else:
            records = _load_from_metrics_dir(ns.pdb_metrics_dir)
    except (FileNotFoundError, ValueError, json.JSONDecodeError) as exc:
        print(f"[clash] ERROR: {exc}", file=sys.stderr)
        return 2

    ns.outdir.mkdir(parents=True, exist_ok=True)

    rows: list[dict] = []
    n_clash_pairs_vals: list[int] = []
    n_vdw_clashes_vals: list[int] = []
    for r in records:
        pdb = str(r.get("pdb", "?"))
        n_cp = _to_int(r.get("n_clash_pairs", r.get("clash")))
        n_vdw = _to_int(r.get("n_vdw_clashes"))
        rows.append({
            "pdb_id": pdb,
            "n_clash_pairs": "" if n_cp is None else n_cp,
            "n_vdw_clashes": "" if n_vdw is None else n_vdw,
        })
        if n_cp is not None:
            n_clash_pairs_vals.append(n_cp)
        if n_vdw is not None:
            n_vdw_clashes_vals.append(n_vdw)

    _write_csv(rows, ns.outdir / "clashscore_histogram.csv")
    summary = build_summary(n_clash_pairs_vals, n_vdw_clashes_vals)
    (ns.outdir / "summary.txt").write_text(summary)
    _plot(
        n_clash_pairs_vals, n_vdw_clashes_vals,
        ns.outdir / "clashscore_histogram.png",
        bin_width=ns.bin_width, x_max=ns.x_max, log_y=ns.log_y,
        title=ns.title,
    )

    print(f"[clash] {len(records)} records; "
          f"clash_pairs n={len(n_clash_pairs_vals)}, "
          f"vdw_clashes n={len(n_vdw_clashes_vals)} "
          f"-> {ns.outdir / 'clashscore_histogram.png'}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
