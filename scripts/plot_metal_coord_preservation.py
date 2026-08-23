#!/usr/bin/env python3
"""Plot metal-coordinating side-chain preservation across a bench (JACS R3).

Reviewer figure: for every protein in the bench that has a metal ion,
compute per-donor Δd = distance(FRUTON model) - distance(crystal), then
plot the bench-wide distribution.  A preservation-preserving pipeline
peaks sharply at 0.0 Å; a pipeline that pushes side chains around
smears the distribution.

Usage:
    python plot_metal_coord_preservation.py \\
        --pair 8ABC=crystals/8ABC.pdb:fruton/8ABC/final_model.pdb \\
        --pair 8XYZ=crystals/8XYZ.pdb:fruton/8XYZ/final_model.pdb \\
        --outdir figures/jacs_r3_metal_preservation/

Each --pair is ``label=crystal_pdb:fruton_pdb``.

Produces:
    metal_coord_preservation.png     Δd distribution histogram + per-donor
                                     count of lost/gained/preserved.
    metal_coord_preservation.csv     one row per matched donor.
    per_pdb_summary.txt              one line per input pair.

Pure Python + matplotlib.  License-free.
"""
from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

from stack_protein_preparation._metal_coord_scan import (
    aggregate_deltas,
    compare_metal_coordination,
    scan_metal_coordination,
    summarise,
)


def _parse_pair(spec: str) -> tuple[str, Path, Path]:
    if "=" not in spec or ":" not in spec.split("=", 1)[1]:
        raise argparse.ArgumentTypeError(
            f"--pair must be label=crystal_pdb:fruton_pdb, got {spec!r}"
        )
    label, rest = spec.split("=", 1)
    crystal, fruton = rest.split(":", 1)
    return label.strip(), Path(crystal.strip()), Path(fruton.strip())


def _parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument(
        "--pair",
        action="append",
        required=True,
        help="Format: label=crystal_pdb:fruton_pdb. Repeatable.",
    )
    p.add_argument("--outdir", type=Path, required=True)
    p.add_argument(
        "--coord-cutoff-A",
        type=float,
        default=3.0,
        help="Metal↔donor cutoff (default 3.0 per Harding 2001).",
    )
    p.add_argument(
        "--bin-width-A",
        type=float,
        default=0.05,
    )
    p.add_argument(
        "--title",
        default="Metal-coordinating side-chain preservation — Δd = FRUTON − crystal",
    )
    return p.parse_args(argv)


def _plot(
    deltas: list[float],
    lost: int,
    gained: int,
    preserved: int,
    outpath: Path,
    bin_width: float,
    title: str,
) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    outpath.parent.mkdir(parents=True, exist_ok=True)

    fig, ax = plt.subplots(figsize=(8.5, 4.5))
    if deltas:
        max_abs = max(abs(min(deltas)), abs(max(deltas)), 0.5)
        edge = max_abs + bin_width
        n_bins = int((2 * edge) / bin_width) + 1
        bins = [-edge + i * bin_width for i in range(n_bins + 1)]
        ax.hist(deltas, bins=bins, color="teal", edgecolor="none")
    ax.axvline(0.0, color="black", linewidth=0.7)
    for x in (-0.5, -0.3, 0.3, 0.5):
        ax.axvline(x, color="grey", linestyle=":", linewidth=0.7)
    ax.set_xlabel("Δd  metal–donor  [Å]  (FRUTON − crystal)")
    ax.set_ylabel("count (matched donors)")
    ax.set_title(title)

    n_total = preserved + lost + gained
    if n_total or deltas:
        legend = (
            f"matched donors    {preserved}\n"
            f"lost   (crystal→FRUTON out of cutoff)   {lost}\n"
            f"gained (FRUTON→crystal absent)          {gained}\n"
            f"Δd stats: n={len(deltas)}"
        )
        if deltas:
            import statistics
            legend += (
                f"\n  mean {statistics.fmean(deltas):+.3f} Å"
                f"\n  median {statistics.median(deltas):+.3f} Å"
                f"\n  max |Δd| {max(abs(d) for d in deltas):.3f} Å"
            )
        ax.text(
            0.02, 0.98, legend,
            transform=ax.transAxes,
            fontsize=8, family="monospace",
            verticalalignment="top",
            bbox={"boxstyle": "round,pad=0.4", "facecolor": "white", "alpha": 0.85},
        )

    fig.tight_layout()
    fig.savefig(outpath, dpi=150)
    plt.close(fig)


def _write_csv(rows: list[dict], outpath: Path) -> None:
    outpath.parent.mkdir(parents=True, exist_ok=True)
    fields = [
        "pdb_id", "metal_chain", "metal_resnum", "metal_element",
        "donor_chain", "donor_resnum", "donor_resname", "donor_atom",
        "distance_crystal_A", "distance_fruton_A", "delta_A", "status",
    ]
    with outpath.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields)
        w.writeheader()
        w.writerows(rows)


def main(argv: list[str] | None = None) -> int:
    ns = _parse_args(argv)

    ns.outdir.mkdir(parents=True, exist_ok=True)

    rows: list[dict] = []
    all_delta_lists: list[list] = []
    lost_total = gained_total = preserved_total = 0
    per_pdb_lines: list[str] = []

    for spec in ns.pair:
        label, crystal_path, fruton_path = _parse_pair(spec)
        cry = scan_metal_coordination(crystal_path, coord_cutoff_A=ns.coord_cutoff_A)
        fru = scan_metal_coordination(fruton_path, coord_cutoff_A=ns.coord_cutoff_A)
        if not cry.ran or not fru.ran:
            per_pdb_lines.append(
                f"[SKIP] {label}: crystal={summarise(cry)} | fruton={summarise(fru)}"
            )
            continue
        deltas = compare_metal_coordination(cry, fru)
        all_delta_lists.append(deltas)
        for d in deltas:
            rows.append({
                "pdb_id": label,
                "metal_chain": d.key.metal_chain,
                "metal_resnum": d.key.metal_resnum,
                "metal_element": d.key.metal_element,
                "donor_chain": d.key.donor_chain,
                "donor_resnum": d.key.donor_resnum,
                "donor_resname": d.key.donor_resname,
                "donor_atom": d.key.donor_atom,
                "distance_crystal_A":
                    "" if d.distance_crystal_A is None else f"{d.distance_crystal_A:.3f}",
                "distance_fruton_A":
                    "" if d.distance_fruton_A is None else f"{d.distance_fruton_A:.3f}",
                "delta_A": "" if d.delta_A() is None else f"{d.delta_A():+.3f}",
                "status": d.status(),
            })
            s = d.status()
            if s == "preserved":
                preserved_total += 1
            elif s == "lost":
                lost_total += 1
            elif s == "gained":
                gained_total += 1
        per_pdb_lines.append(f"[OK]   {label}: {summarise(cry)}  vs  {summarise(fru)}")

    delta_values = aggregate_deltas(all_delta_lists)
    _write_csv(rows, ns.outdir / "metal_coord_preservation.csv")
    _plot(
        delta_values,
        lost=lost_total,
        gained=gained_total,
        preserved=preserved_total,
        outpath=ns.outdir / "metal_coord_preservation.png",
        bin_width=ns.bin_width_A,
        title=ns.title,
    )
    (ns.outdir / "per_pdb_summary.txt").write_text("\n".join(per_pdb_lines) + "\n")

    print(
        f"[plot_metal] pairs={len(ns.pair)} preserved={preserved_total} "
        f"lost={lost_total} gained={gained_total} "
        f"-> {ns.outdir / 'metal_coord_preservation.png'}"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
