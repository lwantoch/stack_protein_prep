#!/usr/bin/env python3
"""Plot FRUTON gate-ablation comparison across bench variants.

Reviewer figure for Nature R4 / JCTC R2: which gate contributes how
much to the observed bench-wide pass rate + rescue count?

Usage:
    python plot_ablation_comparison.py \\
        --baseline runs/mmbsa200_baseline.json \\
        --variant no_plddt=runs/mmbsa200_no_plddt.json:plddt \\
        --variant no_idr=runs/mmbsa200_no_idr.json:idr \\
        --variant no_omega=runs/mmbsa200_no_omega.json:omega \\
        --outdir figures/nature_r4_ablation/

Each --variant argument is ``label=json_path[:disabled_gate]``.  The
``:disabled_gate`` suffix is optional and must match one of the keys in
_ablation.GATE_NAMES if given.

Produces:
    ablation_comparison.png     grouped bar plot (pass_rate + rescued
                                per variant, baseline highlighted)
    ablation_comparison.md      reviewer-facing markdown table
    ablation_comparison.json    machine-readable full report

Pure Python + matplotlib.  License-free.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

from stack_protein_preparation._ablation import (
    build_comparison,
    format_comparison_table,
    rescued_and_lost_counts,
)


def _parse_variant(spec: str) -> tuple[str, str, str | None]:
    """Parse ``label=json_path[:disabled_gate]``."""
    if "=" not in spec:
        raise argparse.ArgumentTypeError(
            f"--variant must be label=json_path[:disabled_gate], got {spec!r}"
        )
    label, rest = spec.split("=", 1)
    if ":" in rest:
        path, gate = rest.rsplit(":", 1)
    else:
        path, gate = rest, None
    return label.strip(), path.strip(), gate.strip() if gate else None


def _parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--baseline", type=Path, required=True)
    p.add_argument(
        "--variant",
        action="append",
        required=True,
        help="Format: label=json_path[:disabled_gate]. Repeatable.",
    )
    p.add_argument("--outdir", type=Path, required=True)
    p.add_argument(
        "--title",
        default="FRUTON gate ablation — MMBSA_200 bench",
    )
    return p.parse_args(argv)


def _plot(comp, outpath: Path, title: str) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    outpath.parent.mkdir(parents=True, exist_ok=True)

    labels = ["baseline"] + [v.variant for v in comp.variants]
    pass_pct = [100 * comp.baseline.pass_rate] + [
        100 * v.pass_rate for v in comp.variants
    ]
    rescued_total = [comp.baseline.total_delta_n] + [
        v.total_delta_n for v in comp.variants
    ]

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11.5, 4.5))

    colours = ["#2b6cb0"] + ["#c05621"] * len(comp.variants)
    ax1.bar(labels, pass_pct, color=colours, edgecolor="none")
    ax1.set_ylabel("gate PASS rate (%)")
    ax1.set_ylim(0, 100)
    ax1.set_title("Pass rate per variant")
    for tick in ax1.get_xticklabels():
        tick.set_rotation(30)
        tick.set_horizontalalignment("right")

    ax2.bar(labels, rescued_total, color=colours, edgecolor="none")
    ax2.set_ylabel("total residues rescued (Σ ΔN)")
    ax2.set_title("Rescued residues per variant")
    for tick in ax2.get_xticklabels():
        tick.set_rotation(30)
        tick.set_horizontalalignment("right")

    fig.suptitle(title)
    fig.tight_layout()
    fig.savefig(outpath, dpi=150)
    plt.close(fig)


def main(argv: list[str] | None = None) -> int:
    ns = _parse_args(argv)

    variant_jsons: dict[str, Path] = {}
    disabled_gate_by_variant: dict[str, str] = {}
    for spec in ns.variant:
        label, path_str, gate = _parse_variant(spec)
        variant_jsons[label] = Path(path_str)
        if gate:
            disabled_gate_by_variant[label] = gate

    ns.outdir.mkdir(parents=True, exist_ok=True)

    try:
        comp = build_comparison(
            ns.baseline, variant_jsons,
            disabled_gate_by_variant=disabled_gate_by_variant,
        )
    except (FileNotFoundError, ValueError) as exc:
        print(f"[plot_ablation] ERROR: {exc}", file=sys.stderr)
        return 2

    (ns.outdir / "ablation_comparison.json").write_text(
        json.dumps(comp.to_dict(), indent=2)
    )
    (ns.outdir / "ablation_comparison.md").write_text(
        format_comparison_table(comp)
    )
    _plot(comp, ns.outdir / "ablation_comparison.png", title=ns.title)

    counts = rescued_and_lost_counts(comp)
    print(
        f"[plot_ablation] baseline pass={comp.baseline.n_passed}/"
        f"{comp.baseline.n_pdbs} ({100 * comp.baseline.pass_rate:.1f}%)"
    )
    for name, (rescued, lost) in counts.items():
        print(f"[plot_ablation]   variant {name}: rescued={rescued} lost={lost}")
    print(f"[plot_ablation] -> {ns.outdir / 'ablation_comparison.png'}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
