#!/usr/bin/env python3
"""Plot metal-coordinating side-chain preservation across a bench (JACS R3).

Reviewer figure for JACS R3: does FRUTON preserve the geometry of
metal-coordinating side chains?  A displaced His-Nε or Asp-Oδ shifts
catalytic geometry and can invalidate the model biologically.

Compares each crystal input PDB against the FRUTON output PDB:
- runs scan_metal_coordination on both
- pairs donor atoms by identity (chain, resnum, atomname, metal)
- emits Δd = distance(FRUTON) − distance(crystal) per matched donor
- histograms Δd bench-wide + per-donor-type breakdown

Usage:
    python plot_metal_donor_preservation.py \\
        --crystal-glob 'bench/*/input_crystal.pdb' \\
        --fruton-glob  'bench/*/final_model.pdb' \\
        --pair-by parent \\
        --outdir figures/jacs_r3_metal/

Produces:
    metal_donor_delta.png     matplotlib Δd histogram + per-donor-type
                              inset bar plot
    metal_donor_delta.csv     one row per matched donor: pdb_id, metal,
                              donor_resname, donor_atom, d_crystal,
                              d_fruton, delta_A, status
    per_pdb_summary.txt       per-PDB one-liner

Pure Python + Bio.PDB + matplotlib.  License-free.
"""
from __future__ import annotations

import argparse
import csv
import glob
import sys
from collections import defaultdict
from pathlib import Path

from stack_protein_preparation._metal_coord_scan import (
    compare_metal_coordination,
    scan_metal_coordination,
    summarise,
)


def _parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--crystal-glob", required=True,
                   help="Glob for crystal (input) PDBs.")
    p.add_argument("--fruton-glob", required=True,
                   help="Glob for FRUTON output PDBs.")
    p.add_argument(
        "--pair-by",
        default="parent",
        choices=("parent", "stem"),
        help="'parent' pairs by parent-directory name (default), "
             "'stem' pairs by file-stem.",
    )
    p.add_argument("--outdir", type=Path, required=True)
    p.add_argument(
        "--coord-cutoff-a",
        type=float,
        default=3.0,
        help="Metal-donor coordination cutoff in Angstroms (default 3.0).",
    )
    p.add_argument(
        "--bin-width",
        type=float,
        default=0.05,
        help="Δd histogram bin width in Å (default 0.05).",
    )
    p.add_argument(
        "--x-range",
        type=float,
        default=1.0,
        help="Symmetric x-axis limit in Å (default ±1.0).",
    )
    p.add_argument(
        "--title",
        default="Metal-donor preservation — MMBSA_200 bench",
    )
    return p.parse_args(argv)


def _key(path: Path, mode: str) -> str:
    return path.parent.name if mode == "parent" else path.stem


def _index_by_key(paths: list[Path], mode: str) -> dict[str, Path]:
    out: dict[str, Path] = {}
    for p in paths:
        k = _key(p, mode)
        if k not in out:
            out[k] = p
    return out


def _plot(
    all_deltas: list[float],
    per_donor_type_deltas: dict[str, list[float]],
    outpath: Path,
    bin_width: float,
    x_range: float,
    title: str,
) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    outpath.parent.mkdir(parents=True, exist_ok=True)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11.5, 4.5))

    # Left: bench-wide Δd histogram
    n_bins = max(1, int(2 * x_range / bin_width))
    bins = [-x_range + i * bin_width for i in range(n_bins + 1)]
    ax1.hist(all_deltas, bins=bins, color="steelblue", edgecolor="none")
    ax1.axvline(0.0, color="black", linewidth=0.8, linestyle="-")
    for x in (-0.3, 0.3):
        ax1.axvline(x, color="grey", linewidth=0.8, linestyle=":")
    ax1.set_xlim(-x_range, x_range)
    ax1.set_xlabel("Δd = d(FRUTON) − d(crystal)  [Å]")
    ax1.set_ylabel("matched donor contacts")
    ax1.set_title("Δd distribution (all donors)")

    if all_deltas:
        import statistics
        mean_d = statistics.fmean(all_deltas)
        med_d = statistics.median(all_deltas)
        n_within_03 = sum(1 for d in all_deltas if abs(d) <= 0.3)
        legend = (
            f"n = {len(all_deltas)}\n"
            f"mean Δd = {mean_d:+.3f} Å\n"
            f"median Δd = {med_d:+.3f} Å\n"
            f"|Δd| ≤ 0.3 Å : {n_within_03} ({100 * n_within_03 / len(all_deltas):.1f}%)"
        )
        ax1.text(
            0.02, 0.98, legend,
            transform=ax1.transAxes,
            fontsize=8, family="monospace",
            verticalalignment="top",
            bbox={"boxstyle": "round,pad=0.4",
                  "facecolor": "white", "alpha": 0.85},
        )

    # Right: per-donor-type mean |Δd| bar
    if per_donor_type_deltas:
        import statistics
        items = sorted(
            per_donor_type_deltas.items(),
            key=lambda it: -len(it[1]),
        )
        labels = [f"{k}\nn={len(v)}" for k, v in items]
        means_abs = [
            statistics.fmean(abs(d) for d in v) if v else 0.0
            for _, v in items
        ]
        ax2.bar(labels, means_abs, color="darkorange", edgecolor="none")
        ax2.set_ylabel("mean |Δd|  [Å]")
        ax2.set_title("Per donor-atom type")
        ax2.axhline(0.3, color="grey", linewidth=0.8, linestyle=":")
        for tick in ax2.get_xticklabels():
            tick.set_rotation(30)
            tick.set_horizontalalignment("right")
    else:
        ax2.text(
            0.5, 0.5, "no matched donors",
            transform=ax2.transAxes, ha="center", va="center",
        )
        ax2.set_title("Per donor-atom type")

    fig.suptitle(title)
    fig.tight_layout()
    fig.savefig(outpath, dpi=150)
    plt.close(fig)


def main(argv: list[str] | None = None) -> int:
    ns = _parse_args(argv)

    crystal_paths = sorted({Path(p) for p in glob.glob(ns.crystal_glob) if Path(p).is_file()})
    fruton_paths = sorted({Path(p) for p in glob.glob(ns.fruton_glob) if Path(p).is_file()})
    if not crystal_paths or not fruton_paths:
        print("[metal_preserve] no crystal or FRUTON PDBs matched", file=sys.stderr)
        return 2

    crystal_by_key = _index_by_key(crystal_paths, ns.pair_by)
    fruton_by_key = _index_by_key(fruton_paths, ns.pair_by)
    shared_keys = sorted(set(crystal_by_key) & set(fruton_by_key))
    if not shared_keys:
        print(
            f"[metal_preserve] no matching PDB pairs (pair-by={ns.pair_by})",
            file=sys.stderr,
        )
        return 2

    ns.outdir.mkdir(parents=True, exist_ok=True)

    rows: list[dict] = []
    all_deltas: list[float] = []
    per_donor_type_deltas: dict[str, list[float]] = defaultdict(list)
    per_pdb_summary: list[str] = []
    n_lost = n_gained = n_preserved = 0

    for pid in shared_keys:
        c_res = scan_metal_coordination(crystal_by_key[pid], coord_cutoff_A=ns.coord_cutoff_a)
        f_res = scan_metal_coordination(fruton_by_key[pid], coord_cutoff_A=ns.coord_cutoff_a)
        if not (c_res.ran and f_res.ran):
            per_pdb_summary.append(f"[SKIP] {pid}: {c_res.fallback_reason or f_res.fallback_reason}")
            continue
        deltas = compare_metal_coordination(c_res, f_res)
        if not deltas:
            per_pdb_summary.append(f"[OK]   {pid}: no metal contacts")
            continue

        for d in deltas:
            status = d.status()
            if status == "preserved":
                n_preserved += 1
            elif status == "lost":
                n_lost += 1
            else:
                n_gained += 1
            delta_val = d.delta_A()
            rows.append({
                "pdb_id": pid,
                "metal_element": d.key.metal_element,
                "metal_chain": d.key.metal_chain,
                "metal_resnum": d.key.metal_resnum,
                "donor_chain": d.key.donor_chain,
                "donor_resnum": d.key.donor_resnum,
                "donor_resname": d.key.donor_resname,
                "donor_atom": d.key.donor_atom,
                "d_crystal_A": ("" if d.distance_crystal_A is None
                                else f"{d.distance_crystal_A:.3f}"),
                "d_fruton_A":  ("" if d.distance_fruton_A is None
                                else f"{d.distance_fruton_A:.3f}"),
                "delta_A":     ("" if delta_val is None else f"{delta_val:+.3f}"),
                "status":      status,
            })
            if delta_val is not None:
                all_deltas.append(delta_val)
                donor_type = f"{d.key.donor_resname}-{d.key.donor_atom}"
                per_donor_type_deltas[donor_type].append(delta_val)

        per_pdb_summary.append(
            f"[OK]   {pid}: crystal→{summarise(c_res)} | fruton→{summarise(f_res)}"
        )

    outdir = ns.outdir
    fields = ["pdb_id", "metal_element", "metal_chain", "metal_resnum",
              "donor_chain", "donor_resnum", "donor_resname", "donor_atom",
              "d_crystal_A", "d_fruton_A", "delta_A", "status"]
    with (outdir / "metal_donor_delta.csv").open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields)
        w.writeheader()
        w.writerows(rows)
    (outdir / "per_pdb_summary.txt").write_text("\n".join(per_pdb_summary) + "\n")
    _plot(
        all_deltas, dict(per_donor_type_deltas),
        outdir / "metal_donor_delta.png",
        bin_width=ns.bin_width, x_range=ns.x_range, title=ns.title,
    )

    print(
        f"[metal_preserve] n_pdbs={len(shared_keys)} "
        f"donors preserved={n_preserved} lost={n_lost} gained={n_gained} "
        f"-> {outdir / 'metal_donor_delta.png'}"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
