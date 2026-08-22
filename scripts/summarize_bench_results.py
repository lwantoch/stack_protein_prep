#!/usr/bin/env python3
"""Summarise a FRUTON benchmark JSON as CSV + markdown for publication tables.

Reads a JSON produced by ``scripts/CESGA_SLURM/fruton_bench48_full.py``
(a list of per-protein dicts) and writes:

  <stem>.csv       -- per-protein table, machine-readable
  <stem>.md        -- reviewer-facing markdown summary with top-fills
                      block, ceiling-case list, and aggregate statistics

Usage:
    python scripts/summarize_bench_results.py <bench_results.json>

Writes to the same directory as the input JSON.
"""
from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path


CSV_COLUMNS = (
    "pdb",
    "base_n",
    "final_n",
    "delta_n",
    "n_gaps",
    "brk",
    "clash",
    "rolled_gaps",
    "refine_seconds",
    "gate_pass",
    "reasons",
)


def write_csv(results: list[dict], output_path: Path) -> None:
    """Emit a per-protein CSV; missing columns default to empty."""
    with output_path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh)
        writer.writerow(CSV_COLUMNS)
        for r in results:
            row: list = []
            for col in CSV_COLUMNS:
                val = r.get(col, "")
                if isinstance(val, list):
                    val = "; ".join(str(x) for x in val) if val else ""
                row.append(val)
            writer.writerow(row)


def summarise_markdown(results: list[dict], input_path: Path) -> str:
    """Reviewer-facing markdown table + aggregates."""
    n = len(results)
    passed = sum(1 for r in results if r.get("gate_pass"))
    filled = [r for r in results if r.get("delta_n", 0) > 0]
    total_residues = sum(r.get("delta_n", 0) for r in results)
    zero_broken = sum(1 for r in results if r.get("brk", 0) == 0)
    zero_clash = sum(1 for r in results if r.get("clash", 0) == 0)

    lines: list[str] = []
    lines.append(f"# FRUTON benchmark summary — {input_path.name}")
    lines.append("")
    lines.append("## Aggregate")
    lines.append("")
    lines.append(f"| Metric | Value |")
    lines.append(f"|---|---|")
    lines.append(f"| Proteins tested | {n} |")
    lines.append(f"| Relative gate PASS | {passed}/{n} ({100 * passed / n:.1f}%) |")
    lines.append(f"| Proteins with fills | {len(filled)}/{n} |")
    lines.append(f"| Total residues rescued | **{total_residues}** |")
    lines.append(f"| Zero broken peptide bonds | {zero_broken}/{n} |")
    lines.append(f"| Zero clash gain | {zero_clash}/{n} |")
    lines.append("")

    if filled:
        top = sorted(filled, key=lambda r: -r.get("delta_n", 0))[:15]
        lines.append("## Top fills")
        lines.append("")
        lines.append("| PDB | ΔN | broken | clash gain | gate |")
        lines.append("|---|---|---|---|---|")
        for r in top:
            gate = "PASS" if r.get("gate_pass") else "FAIL"
            lines.append(
                f"| {r['pdb']} | +{r.get('delta_n', 0)} | "
                f"{r.get('brk', 0)} | {r.get('clash', 0)} | {gate} |"
            )
        lines.append("")

    fails = [r for r in results if not r.get("gate_pass")]
    if fails:
        lines.append(f"## Failed proteins ({len(fails)})")
        lines.append("")
        for r in fails:
            reasons = r.get("reasons") or []
            lines.append(f"- **{r['pdb']}**: ΔN={r.get('delta_n', 0)}, "
                         f"brk={r.get('brk', 0)}, clash={r.get('clash', 0)}")
            for reason in reasons[:3]:
                lines.append(f"  - {reason}")
        lines.append("")

    return "\n".join(lines) + "\n"


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("json_path", type=Path, help="Path to bench results JSON")
    args = parser.parse_args()

    if not args.json_path.is_file():
        print(f"ERROR: {args.json_path} not found", file=sys.stderr)
        return 1

    with args.json_path.open() as fh:
        results = json.load(fh)
    if not isinstance(results, list):
        print(f"ERROR: expected a list at top level; got {type(results).__name__}", file=sys.stderr)
        return 1

    stem = args.json_path.stem
    csv_path = args.json_path.with_name(f"{stem}.csv")
    md_path = args.json_path.with_name(f"{stem}.md")

    write_csv(results, csv_path)
    md_path.write_text(summarise_markdown(results, args.json_path), encoding="utf-8")

    print(f"wrote {csv_path}")
    print(f"wrote {md_path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
