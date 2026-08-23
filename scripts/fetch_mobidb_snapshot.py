#!/usr/bin/env python3
"""Refresh src/stack_protein_preparation/data/mobidb_snapshot.json for the
UniProt accessions used across FRUTON benchmarks.

Rationale (JCTC R4 reproducibility): the pipeline's IDR gate calls
MobiDB per-protein at runtime.  Caching those responses in-repo
gives us (1) deterministic reruns, (2) no network dependency during
CI, (3) a citation-able snapshot of the disorder annotations used at
publication time.

Usage:
    python scripts/fetch_mobidb_snapshot.py
        [--accessions P04637 P10636 ...]
        [--from-benchmark  /path/to/*.json]

Live-fetches each accession with use_cache=False, prints a diff
against the existing snapshot, and (unless --dry-run) writes the
updated snapshot.  Safe to re-run; idempotent per accession.
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

# Make src/ importable without an editable install
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from stack_protein_preparation import _uniprot_idr  # noqa: E402


def collect_accessions_from_benchmark(bench_json: Path) -> set[str]:
    """Best-effort extraction of UniProt accessions from a FRUTON bench
    results JSON.  Falls back to empty set on unexpected shape."""
    try:
        data = json.loads(bench_json.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return set()
    out: set[str] = set()
    for entry in (data if isinstance(data, list) else [data]):
        if not isinstance(entry, dict):
            continue
        acc = entry.get("uniprot_id") or entry.get("uniprot") or entry.get("accession")
        if isinstance(acc, str) and acc:
            out.add(acc.strip().upper())
    return out


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--accessions", nargs="*", default=[],
                        help="Explicit UniProt accessions to refresh")
    parser.add_argument("--from-benchmark", type=Path, default=None,
                        help="Read accessions from a bench_results.json")
    parser.add_argument("--dry-run", action="store_true",
                        help="Print diff without writing")
    parser.add_argument("--sleep", type=float, default=0.3,
                        help="Sleep seconds between API calls (default 0.3)")
    args = parser.parse_args()

    targets: set[str] = set(a.strip().upper() for a in args.accessions if a.strip())
    if args.from_benchmark and args.from_benchmark.is_file():
        targets |= collect_accessions_from_benchmark(args.from_benchmark)
    if not targets:
        print("No accessions supplied; nothing to refresh.", file=sys.stderr)
        return 2

    existing = _uniprot_idr._load_cache()
    changes: list[tuple[str, str]] = []
    n_ok = n_fail = 0
    for i, acc in enumerate(sorted(targets), 1):
        regions = _uniprot_idr.fetch_uniprot_disorder_regions(acc, use_cache=False)
        if regions is None:
            print(f"[{i}/{len(targets)}] {acc}: FETCH FAILED (network / unknown accession)")
            n_fail += 1
        else:
            old = existing.get(acc)
            if old != regions:
                changes.append((acc, f"{old!r} -> {regions!r}"))
                if not args.dry_run:
                    _uniprot_idr.save_cache_entry(acc, regions)
            n_ok += 1
            marker = "  (changed)" if old != regions else ""
            print(f"[{i}/{len(targets)}] {acc}: {len(regions)} region(s){marker}")
        if args.sleep and i < len(targets):
            time.sleep(args.sleep)

    print(f"\n{n_ok} ok, {n_fail} failed, {len(changes)} changed")
    if args.dry_run:
        print("(dry run -- snapshot NOT written)")
    else:
        print(f"snapshot: {_uniprot_idr._CACHE_PATH}")
    return 0 if n_fail == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
