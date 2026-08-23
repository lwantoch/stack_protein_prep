#!/usr/bin/env python3
"""Blind-test bench figure: Cα-RMSD of FRUTON fills vs held-out crystal (Nature R3).

Reviewer-facing figure for the crystal-held-out blind test:
- histogram of per-residue Cα distances for held-out residues
- side panel: per-PDB overall Cα-RMSD

Usage:
    python plot_blind_test_rmsd.py \\
        --bench-spec metrics/blind_test_bench.json \\
        --outdir figures/nature_r3_blind/

bench_spec.json format:
[
    {
        "pdb_id": "8ABC",
        "crystal_pdb": "bench/8ABC/input_crystal.pdb",
        "filled_pdb":  "bench/8ABC/final_model.pdb",
        "held_out": [{"chain": "A", "first_resnum": 50, "last_resnum": 60}]
    },
    ...
]

Produces:
    blind_test_rmsd.png      histogram (per-residue Cα distance) +
                             per-PDB Cα-RMSD box (log-y optional).
    blind_test_rmsd.csv      one row per held-out residue.
    per_pdb_summary.txt      per-PDB one-liner.

Pure Python + Bio.PDB + matplotlib.  License-free.
"""
from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path

from stack_protein_preparation._blind_test import (
    GapRange,
    score_blind_fill,
    summarise,
)


def _parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--bench-spec", type=Path, required=True,
                   help="JSON list of per-PDB blind-test entries (see docstring).")
    p.add_argument("--outdir", type=Path, required=True)
    p.add_argument(
        "--bin-width",
        type=float,
        default=0.2,
        help="Histogram bin width in Å (default 0.2).",
    )
    p.add_argument(
        "--x-max",
        type=float,
        default=8.0,
        help="Histogram x-axis upper bound in Å (default 8.0).",
    )
    p.add_argument(
        "--log-y",
        action="store_true",
        help="Log-scale the histogram y-axis (highlights the tail).",
    )
    p.add_argument(
        "--title",
        default="Blind-test Cα-RMSD — MMBSA_200",
    )
    return p.parse_args(argv)


def _load_bench_spec(path: Path) -> list[dict]:
    with path.open() as fh:
        data = json.load(fh)
    if not isinstance(data, list):
        raise ValueError(f"{path} is not a JSON list")
    return data


def _entry_to_ranges(entry: dict) -> list[GapRange]:
    ranges: list[GapRange] = []
    for r in entry.get("held_out", []):
        try:
            ranges.append(GapRange(
                chain=str(r["chain"]),
                first_resnum=int(r["first_resnum"]),
                last_resnum=int(r["last_resnum"]),
            ))
        except (KeyError, TypeError, ValueError):
            continue
    return ranges


def _plot(
    per_residue_ca: list[float],
    per_pdb_rmsd: list[float],
    outpath: Path,
    bin_width: float,
    x_max: float,
    log_y: bool,
    title: str,
) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    outpath.parent.mkdir(parents=True, exist_ok=True)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11.5, 4.5))

    # Left: per-residue Cα distance histogram
    n_bins = max(1, int(x_max / bin_width))
    bins = [i * bin_width for i in range(n_bins + 1)]
    ax1.hist(
        [v for v in per_residue_ca if v <= x_max],
        bins=bins, color="steelblue", edgecolor="none",
    )
    if log_y:
        ax1.set_yscale("log")
    ax1.axvline(1.0, color="grey", linewidth=0.9, linestyle=":", label="1.0 Å")
    ax1.axvline(2.5, color="grey", linewidth=0.9, linestyle=":", label="2.5 Å")
    ax1.set_xlim(0.0, x_max)
    ax1.set_xlabel("per-residue Cα distance vs held-out crystal (Å)")
    ax1.set_ylabel("count" + (" (log)" if log_y else ""))
    ax1.set_title("Per-residue Cα distance")
    if per_residue_ca:
        import statistics
        mean_d = statistics.fmean(per_residue_ca)
        med_d = statistics.median(per_residue_ca)
        n_within_1 = sum(1 for v in per_residue_ca if v <= 1.0)
        n_within_25 = sum(1 for v in per_residue_ca if v <= 2.5)
        legend = (
            f"n = {len(per_residue_ca)}\n"
            f"mean = {mean_d:.3f} Å\n"
            f"median = {med_d:.3f} Å\n"
            f"≤ 1.0 Å : {n_within_1} ({100 * n_within_1 / len(per_residue_ca):.1f}%)\n"
            f"≤ 2.5 Å : {n_within_25} ({100 * n_within_25 / len(per_residue_ca):.1f}%)"
        )
        ax1.text(
            0.97, 0.97, legend,
            transform=ax1.transAxes,
            fontsize=8, family="monospace",
            verticalalignment="top", horizontalalignment="right",
            bbox={"boxstyle": "round,pad=0.4",
                  "facecolor": "white", "alpha": 0.85},
        )

    # Right: per-PDB overall Cα-RMSD as sorted bar (helps reviewer see spread)
    if per_pdb_rmsd:
        sorted_rmsd = sorted(per_pdb_rmsd)
        ax2.bar(range(len(sorted_rmsd)), sorted_rmsd,
                color="darkorange", edgecolor="none")
        ax2.set_xlabel("bench PDB (sorted by RMSD)")
        ax2.set_ylabel("Cα-RMSD (Å)")
        ax2.set_title("Per-PDB overall Cα-RMSD")
        ax2.axhline(1.0, color="grey", linewidth=0.9, linestyle=":")
        ax2.axhline(2.5, color="grey", linewidth=0.9, linestyle=":")
    else:
        ax2.text(0.5, 0.5, "no PDBs scored",
                 ha="center", va="center", transform=ax2.transAxes)
        ax2.set_title("Per-PDB overall Cα-RMSD")

    fig.suptitle(title)
    fig.tight_layout()
    fig.savefig(outpath, dpi=150)
    plt.close(fig)


def main(argv: list[str] | None = None) -> int:
    ns = _parse_args(argv)
    if not ns.bench_spec.is_file():
        print(f"[blind_test] ERROR: spec not found: {ns.bench_spec}",
              file=sys.stderr)
        return 2
    try:
        spec = _load_bench_spec(ns.bench_spec)
    except (ValueError, json.JSONDecodeError) as exc:
        print(f"[blind_test] ERROR: {exc}", file=sys.stderr)
        return 2

    ns.outdir.mkdir(parents=True, exist_ok=True)

    per_residue_ca: list[float] = []
    per_pdb_rmsd: list[float] = []
    per_pdb_summary: list[str] = []
    rows: list[dict] = []

    for entry in spec:
        pid = str(entry.get("pdb_id", "?"))
        crystal = Path(entry.get("crystal_pdb", ""))
        filled = Path(entry.get("filled_pdb", ""))
        ranges = _entry_to_ranges(entry)
        if not ranges:
            per_pdb_summary.append(f"[SKIP] {pid}: no held_out ranges")
            continue

        result = score_blind_fill(crystal, filled, ranges)
        if not result.ran:
            per_pdb_summary.append(f"[SKIP] {pid}: {result.fallback_reason}")
            continue

        per_pdb_summary.append(f"[OK]   {pid}: {summarise(result)}")

        for rng in result.ranges:
            for r in rng.residues:
                if r.ca_distance_A is not None:
                    per_residue_ca.append(r.ca_distance_A)
                rows.append({
                    "pdb_id": pid,
                    "chain": r.chain,
                    "resnum": r.resnum,
                    "crystal_resname": r.crystal_resname,
                    "filled_resname": r.filled_resname or "",
                    "ca_distance_A": ("" if r.ca_distance_A is None
                                      else f"{r.ca_distance_A:.3f}"),
                    "backbone_rmsd_A": ("" if r.backbone_rmsd_A is None
                                        else f"{r.backbone_rmsd_A:.3f}"),
                    "resname_match": "1" if r.resname_matches() else "0",
                })
        overall = result.overall_ca_rmsd_A()
        if overall is not None:
            per_pdb_rmsd.append(overall)

    fields = ["pdb_id", "chain", "resnum", "crystal_resname", "filled_resname",
              "ca_distance_A", "backbone_rmsd_A", "resname_match"]
    with (ns.outdir / "blind_test_rmsd.csv").open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields)
        w.writeheader()
        w.writerows(rows)
    (ns.outdir / "per_pdb_summary.txt").write_text(
        "\n".join(per_pdb_summary) + "\n"
    )
    _plot(
        per_residue_ca, per_pdb_rmsd,
        ns.outdir / "blind_test_rmsd.png",
        bin_width=ns.bin_width, x_max=ns.x_max, log_y=ns.log_y,
        title=ns.title,
    )

    print(
        f"[blind_test] pdbs_scored={len(per_pdb_rmsd)} "
        f"residues_scored={len(per_residue_ca)} "
        f"-> {ns.outdir / 'blind_test_rmsd.png'}"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
