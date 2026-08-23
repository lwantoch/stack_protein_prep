#!/usr/bin/env python3
"""Aggregate per-PDB result JSONs (from the SLURM array bench) into a
single ``fruton_bench_mmbsa200_results.json`` + emit summary + kick
off standard analysis figures.

Usage:
    python fruton_bench_mmbsa200_aggregate.py <BENCH_OUT_DIR>

Reads every ``<PDB>.json`` in ``<BENCH_OUT_DIR>``, concatenates into a
list, writes ``fruton_bench_mmbsa200_results.json``, prints summary.
"""
from __future__ import annotations

import json
import statistics
import sys
from pathlib import Path


def main(argv: list[str]) -> int:
    if len(argv) < 2:
        print("usage: fruton_bench_mmbsa200_aggregate.py <BENCH_OUT_DIR>",
              file=sys.stderr)
        return 2
    out_dir = Path(argv[1])
    if not out_dir.is_dir():
        print(f"[aggregate] directory not found: {out_dir}", file=sys.stderr)
        return 2

    per_pdb_files = sorted(
        p for p in out_dir.glob("*.json")
        if p.name not in {
            "fruton_bench_mmbsa200_results.json",
            "fruton_bench_mmbsa200_partial.json",
        }
    )
    if not per_pdb_files:
        print(f"[aggregate] no per-PDB JSONs in {out_dir}", file=sys.stderr)
        return 2

    results: list[dict] = []
    for p in per_pdb_files:
        try:
            results.append(json.loads(p.read_text()))
        except Exception as exc:
            print(f"[aggregate] SKIP {p.name}: {exc}", file=sys.stderr)

    n_ok = sum(1 for r in results if "error" not in r)
    n_pass = sum(1 for r in results if r.get("gate_pass"))
    total_fills = sum(int(r.get("delta_n", 0)) for r in results if "delta_n" in r)
    mean_clash = (
        statistics.fmean(int(r.get("clash", 0)) for r in results if "error" not in r)
        if n_ok else 0.0
    )

    combined = out_dir / "fruton_bench_mmbsa200_results.json"
    combined.write_text(json.dumps(results, indent=2, default=str))

    print(f"[aggregate] n_pdbs={len(results)} ok={n_ok} pass={n_pass} "
          f"total_fills={total_fills} mean_clash={mean_clash:.2f}")
    print(f"[aggregate] -> {combined}")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
