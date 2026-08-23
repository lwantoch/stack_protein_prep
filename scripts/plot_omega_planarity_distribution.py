#!/usr/bin/env python3
"""Plot ω peptide-bond planarity distribution across a bench of PDBs.

Reviewer figure for JACS R1 concern: does FRUTON preserve peptide-bond
planarity across the MMBSA_200 benchmark set?

Produces:
    <outdir>/omega_distribution.png        matplotlib histogram, 1° bins,
                                           annotates trans / cis-Pro /
                                           cis-nonPro / non-planar counts.
    <outdir>/omega_distribution.csv        one row per residue pair:
                                           pdb_id, chain, resnum_i,
                                           resnum_j, resname_j, omega_deg,
                                           kind

Usage:
    python plot_omega_planarity_distribution.py \\
        --pdb-glob '/path/to/bench/*/final_model.pdb' \\
        --outdir figures/jacs_r1_omega/

Free / license-independent: pure Python + Bio.PDB + matplotlib.  Never
requires phenix / MolProbity / any commercial tool.
"""
from __future__ import annotations

import argparse
import csv
import glob
import sys
from pathlib import Path

from stack_protein_preparation._omega_scan import (
    scan_omega_dihedrals,
    summarise,
)


def _parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument(
        "--pdb-glob",
        action="append",
        required=True,
        help="Glob pattern for input PDBs; can be repeated.",
    )
    p.add_argument("--outdir", type=Path, required=True)
    p.add_argument(
        "--pdb-id-from",
        default="parent",
        choices=("parent", "stem"),
        help="'parent' = use parent directory name as pdb_id "
             "(good for bench/8XYZ/final_model.pdb); "
             "'stem' = use file stem (good for bench/8XYZ.pdb).",
    )
    p.add_argument(
        "--bin-width",
        type=float,
        default=1.0,
        help="Histogram bin width in degrees (default 1.0).",
    )
    p.add_argument(
        "--log-y",
        action="store_true",
        help="Use log scale on the y-axis (highlights the rare non-planar tail).",
    )
    p.add_argument(
        "--title",
        default="ω peptide-bond planarity — MMBSA_200 bench",
    )
    return p.parse_args(argv)


def _collect_pdb_paths(patterns: list[str]) -> list[Path]:
    paths: list[Path] = []
    for pat in patterns:
        for match in glob.glob(pat):
            p = Path(match)
            if p.is_file():
                paths.append(p)
    return sorted(set(paths))


def _pdb_id(path: Path, mode: str) -> str:
    return path.parent.name if mode == "parent" else path.stem


def _write_csv(rows: list[dict], outpath: Path) -> None:
    outpath.parent.mkdir(parents=True, exist_ok=True)
    fields = ["pdb_id", "chain", "resnum_i", "resnum_j",
              "resname_i", "resname_j", "omega_deg", "kind"]
    with outpath.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields)
        w.writeheader()
        w.writerows(rows)


def _plot(
    omega_values: list[float],
    counts: dict[str, int],
    outpath: Path,
    bin_width: float,
    log_y: bool,
    title: str,
) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    outpath.parent.mkdir(parents=True, exist_ok=True)

    fig, ax = plt.subplots(figsize=(9.0, 4.5))
    bins = [-180.0 + i * bin_width for i in range(int(360 / bin_width) + 1)]
    ax.hist(omega_values, bins=bins, color="steelblue", edgecolor="none")
    if log_y:
        ax.set_yscale("log")

    # Reference lines: cis boundaries ±30, trans boundaries ±150
    for x in (-150.0, -30.0, 30.0, 150.0):
        ax.axvline(x, color="grey", linestyle=":", linewidth=0.8)

    ax.set_xlim(-180.0, 180.0)
    ax.set_xlabel("ω dihedral CA(i)-C(i)-N(i+1)-CA(i+1) [°]")
    ax.set_ylabel("count" + (" (log)" if log_y else ""))
    ax.set_title(title)

    n_total = sum(counts.values())
    if n_total:
        legend = (
            f"n={n_total} pairs\n"
            f"trans        {counts.get('trans', 0):>6}  "
            f"({100 * counts.get('trans', 0) / n_total:>5.2f}%)\n"
            f"cis-Pro      {counts.get('cis_pro', 0):>6}  "
            f"({100 * counts.get('cis_pro', 0) / n_total:>5.2f}%)\n"
            f"cis-nonPro   {counts.get('cis_nonpro', 0):>6}  "
            f"({100 * counts.get('cis_nonpro', 0) / n_total:>5.2f}%)\n"
            f"non-planar   {counts.get('non_planar', 0):>6}  "
            f"({100 * counts.get('non_planar', 0) / n_total:>5.2f}%)"
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


def main(argv: list[str] | None = None) -> int:
    ns = _parse_args(argv)
    pdb_paths = _collect_pdb_paths(ns.pdb_glob)
    if not pdb_paths:
        print("[plot_omega] no PDBs matched the given glob(s)", file=sys.stderr)
        return 2

    ns.outdir.mkdir(parents=True, exist_ok=True)

    rows: list[dict] = []
    all_omegas: list[float] = []
    counts: dict[str, int] = {
        "trans": 0, "cis_pro": 0, "cis_nonpro": 0, "non_planar": 0,
    }

    per_pdb_summary: list[str] = []
    for path in pdb_paths:
        result = scan_omega_dihedrals(path)
        if not result.ran or not result.passed:
            per_pdb_summary.append(
                f"[SKIP] {path}: {result.fallback_reason or 'no entries'}"
            )
            continue
        pid = _pdb_id(path, ns.pdb_id_from)
        for entry in result.entries:
            rows.append({
                "pdb_id": pid,
                "chain": entry.chain,
                "resnum_i": entry.resnum_i,
                "resnum_j": entry.resnum_j,
                "resname_i": entry.resname_i,
                "resname_j": entry.resname_j,
                "omega_deg": f"{entry.omega_deg:.3f}",
                "kind": entry.kind,
            })
            all_omegas.append(entry.omega_deg)
            counts[entry.kind] = counts.get(entry.kind, 0) + 1
        per_pdb_summary.append(f"[OK]   {pid}: {summarise(result)}")

    _write_csv(rows, ns.outdir / "omega_distribution.csv")
    _plot(
        all_omegas, counts,
        ns.outdir / "omega_distribution.png",
        bin_width=ns.bin_width,
        log_y=ns.log_y,
        title=ns.title,
    )

    (ns.outdir / "per_pdb_summary.txt").write_text(
        "\n".join(per_pdb_summary) + "\n"
    )

    print(f"[plot_omega] n_pdbs_scanned={len(pdb_paths)} "
          f"n_pairs={sum(counts.values())} "
          f"cis_nonpro={counts.get('cis_nonpro', 0)} "
          f"non_planar={counts.get('non_planar', 0)} "
          f"-> {ns.outdir / 'omega_distribution.png'}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
