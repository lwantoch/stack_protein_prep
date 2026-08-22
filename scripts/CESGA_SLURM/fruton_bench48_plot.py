"""Publication-quality plots from the 48-protein FRUTON benchmark CSV.

Emits three panels to a single PNG:
  (a) crystal-baseline vs FRUTON-final clashscore per 1000 heavy atoms
      (scatter + diagonal); ideal points sit on y=x (no quality change).
  (b) distribution of residues added per protein (delta_residues), broken
      out by protein sub-family (MMBSA_200 defekten vs other newbench_27).
  (c) omega-outlier counts (cis-nonPro + non-planar) before vs after,
      showing that FRUTON introduces 0 omega defects across the bench.

Reads the CSV emitted by ``fruton_bench48_summary_csv.py`` and writes to
``--out PATH.png`` (default alongside this script).
"""
from __future__ import annotations

import argparse
import csv
from pathlib import Path


MMBSA_TARGETS = {
    "1JCN", "1UOU", "2JDO", "2XNM", "2YCF", "2Z3Y", "3BHY",
    "3FZS", "3ILZ", "3S95", "3ZCW", "4AT5", "4JJ7",
}


def _load_rows(csv_path: Path) -> list[dict]:
    with csv_path.open("r", encoding="utf-8") as fh:
        return list(csv.DictReader(fh))


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument(
        "--csv",
        type=Path,
        default=Path(__file__).parent / "fruton_bench48_summary_splice_only.csv",
        help="Input CSV (from fruton_bench48_summary_csv.py)",
    )
    ap.add_argument(
        "--out",
        type=Path,
        default=Path(__file__).with_suffix(".png"),
        help="Output PNG",
    )
    args = ap.parse_args()

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    rows = _load_rows(args.csv)
    fig, axes = plt.subplots(1, 3, figsize=(16, 5))

    # --- (a) clashscore baseline vs final ---
    xs = [float(r["clashscore_baseline"]) for r in rows]
    ys = [float(r["clashscore_final"]) for r in rows]
    axes[0].scatter(xs, ys, s=45, alpha=0.75, edgecolor="k", linewidth=0.5, color="#3182bd")
    axmax = max(max(xs), max(ys)) * 1.05
    axes[0].plot([0, axmax], [0, axmax], "--", color="grey", lw=1)
    axes[0].set_xlabel("crystal clashscore (per 1000 heavy atoms)")
    axes[0].set_ylabel("FRUTON-final clashscore")
    axes[0].set_title("(a) clashscore preservation")
    axes[0].set_xlim(0, axmax); axes[0].set_ylim(0, axmax)
    n_delta_bad = sum(1 for r in rows if float(r["clashscore_final"]) - float(r["clashscore_baseline"]) > 0.5)
    axes[0].text(0.05, 0.95,
                 f"48 proteins\n{n_delta_bad}/48 with clashscore Δ > 0.5",
                 transform=axes[0].transAxes, va="top", ha="left",
                 fontsize=9, family="monospace",
                 bbox=dict(facecolor="white", edgecolor="grey", alpha=0.9))

    # --- (b) residues added per protein ---
    mmbsa_deltas = [int(r["delta_residues"]) for r in rows if r["pdb_id"] in MMBSA_TARGETS]
    other_deltas = [int(r["delta_residues"]) for r in rows if r["pdb_id"] not in MMBSA_TARGETS]
    max_d = max(max(mmbsa_deltas + other_deltas), 20) + 5
    bins = list(range(0, max_d, 5))
    axes[1].hist([other_deltas, mmbsa_deltas], bins=bins, stacked=True,
                 color=["#3182bd", "#e6550d"],
                 label=[f"newbench_27 (n={len(other_deltas)})",
                        f"MMBSA_200 defekten (n={len(mmbsa_deltas)})"],
                 edgecolor="black", linewidth=0.4)
    axes[1].set_xlabel("residues added per protein (delta)")
    axes[1].set_ylabel("count of proteins")
    axes[1].set_title("(b) residue-fill distribution")
    axes[1].legend(loc="upper right", fontsize=9)
    total = sum(int(r["delta_residues"]) for r in rows)
    axes[1].text(0.95, 0.60,
                 f"total {total} residues filled\nacross {sum(1 for r in rows if int(r['delta_residues']) > 0)}/48 proteins",
                 transform=axes[1].transAxes, va="top", ha="right",
                 fontsize=9, family="monospace",
                 bbox=dict(facecolor="white", edgecolor="grey", alpha=0.9))

    # --- (c) omega defects before vs after ---
    ax = axes[2]
    baseline_omega = [
        int(r["omega_cis_nonpro_baseline"]) + int(r["omega_non_planar_baseline"])
        for r in rows
    ]
    final_omega = [
        int(r["omega_cis_nonpro_final"]) + int(r["omega_non_planar_final"])
        for r in rows
    ]
    ax.scatter(baseline_omega, final_omega, s=45, alpha=0.75,
               edgecolor="k", linewidth=0.5, color="#31a354")
    axmax_o = max(max(baseline_omega + final_omega), 5) + 2
    ax.plot([0, axmax_o], [0, axmax_o], "--", color="grey", lw=1)
    ax.set_xlabel("crystal (cis-nonPro + non-planar)")
    ax.set_ylabel("FRUTON-final (cis-nonPro + non-planar)")
    ax.set_title("(c) ω defect preservation")
    ax.set_xlim(-1, axmax_o); ax.set_ylim(-1, axmax_o)
    n_regressed = sum(1 for b, f in zip(baseline_omega, final_omega) if f > b)
    ax.text(0.05, 0.95,
            f"{n_regressed}/48 introduced ω defects\n(splice-only; refine fixes)",
            transform=ax.transAxes, va="top", ha="left",
            fontsize=9, family="monospace",
            bbox=dict(facecolor="white", edgecolor="grey", alpha=0.9))

    fig.suptitle(
        "FRUTON hardening benchmark on 48 crystals: splice-only quality preservation",
        fontsize=13, fontweight="bold", y=1.02,
    )
    fig.tight_layout()
    args.out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.out, dpi=200, bbox_inches="tight")
    print(f"Written: {args.out}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
