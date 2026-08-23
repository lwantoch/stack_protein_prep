#!/usr/bin/env python3
"""Report coverage of family_by_pdb_seed.json against a bench-results JSON.

Reviewer-facing helper for Nature R1 stratified analysis: given a
bench-results JSON (list of per-PDB records), print which PDB IDs are
already covered by the shipped ``family_by_pdb_seed.json`` and which
ones need a manual family assignment before ``plot_family_stratification.py``
can produce a complete figure.

Usage:
    python scripts/build_family_assignment.py \\
        --bench-json runs/mmbsa200_baseline.json \\
        [--extra-map metrics/my_extra_families.json] \\
        [--emit merged_family_map.json]

Prints:
    covered:    <n> PDBs already have a canonical family
    unassigned: <n> PDBs missing from every supplied map

If ``--emit`` is passed, writes a merged ``{pdb_id: family_label}``
JSON to that path.  Unassigned PDBs land under the ``__unassigned__``
label so the plot script can still render them (Nature R1 wants an
honest bucket for uncategorised proteins, not a silent drop).

Pure Python + stdlib.  License-free.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

from stack_protein_preparation._family_stratification import (
    FAMILY_LABELS,
    load_family_mapping,
)


_SEED_PATH = (
    Path(__file__).resolve().parents[1]
    / "src" / "stack_protein_preparation" / "data" / "family_by_pdb_seed.json"
)


def _parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--bench-json", type=Path, required=True,
                   help="Bench-results JSON (list of {pdb, ...} records).")
    p.add_argument(
        "--seed", type=Path, default=_SEED_PATH,
        help="Family seed JSON (defaults to shipped family_by_pdb_seed.json).",
    )
    p.add_argument(
        "--extra-map", type=Path, action="append", default=None,
        help="Additional per-PDB family map JSONs (repeatable). Values "
             "in later maps overwrite earlier ones.",
    )
    p.add_argument(
        "--emit", type=Path, default=None,
        help="Write merged mapping (including __unassigned__ for missing "
             "PDBs) to this path.",
    )
    p.add_argument(
        "--fill-unassigned",
        action="store_true",
        help="Include __unassigned__ entries in the emitted merged map "
             "(otherwise emit only the covered entries).",
    )
    return p.parse_args(argv)


def _load_bench_pdb_ids(path: Path) -> list[str]:
    """Return uppercase, deduplicated bench PDB ids preserving first-seen order."""
    with path.open() as fh:
        data = json.load(fh)
    if not isinstance(data, list):
        raise ValueError(f"{path} is not a list of records")
    seen: set[str] = set()
    out: list[str] = []
    for record in data:
        pid = str(record.get("pdb", "")).strip().upper()
        if not pid or pid in seen:
            continue
        seen.add(pid)
        out.append(pid)
    return out


def _merge_maps(seed: dict[str, str], extras: list[dict[str, str]]) -> dict[str, str]:
    """Merge maps left-to-right; later entries override earlier ones."""
    merged = dict(seed)
    for extra in extras:
        for k, v in extra.items():
            merged[k.upper()] = v
    return merged


def build_coverage_report(
    bench_pdb_ids: list[str],
    merged_map: dict[str, str],
) -> dict:
    """Compute coverage stats + per-family + unassigned lists."""
    covered = [p for p in bench_pdb_ids if p in merged_map]
    unassigned = [p for p in bench_pdb_ids if p not in merged_map]
    by_family: dict[str, list[str]] = {}
    for p in covered:
        by_family.setdefault(merged_map[p], []).append(p)
    return {
        "n_bench_pdbs": len(bench_pdb_ids),
        "n_covered": len(covered),
        "n_unassigned": len(unassigned),
        "coverage_pct": (100.0 * len(covered) / len(bench_pdb_ids)) if bench_pdb_ids else 0.0,
        "by_family": {k: sorted(v) for k, v in sorted(by_family.items())},
        "unassigned": sorted(unassigned),
    }


def format_report(report: dict) -> str:
    lines: list[str] = []
    lines.append(
        f"bench_pdbs={report['n_bench_pdbs']}  "
        f"covered={report['n_covered']} ({report['coverage_pct']:.1f}%)  "
        f"unassigned={report['n_unassigned']}"
    )
    if report["by_family"]:
        lines.append("")
        lines.append("Per-family coverage:")
        for family, pdbs in report["by_family"].items():
            lines.append(f"  {family:<20} n={len(pdbs)}  {', '.join(pdbs[:8])}"
                         + ("…" if len(pdbs) > 8 else ""))
    if report["unassigned"]:
        lines.append("")
        lines.append(f"Unassigned PDBs (need manual family label): "
                     f"{len(report['unassigned'])}")
        preview = report["unassigned"][:15]
        lines.append("  " + ", ".join(preview) + ("…" if len(report["unassigned"]) > 15 else ""))
    return "\n".join(lines)


def main(argv: list[str] | None = None) -> int:
    ns = _parse_args(argv)
    if not ns.bench_json.is_file():
        print(f"[family] ERROR: bench json not found: {ns.bench_json}",
              file=sys.stderr)
        return 2
    if not ns.seed.is_file():
        print(f"[family] ERROR: seed not found: {ns.seed}", file=sys.stderr)
        return 2

    try:
        bench_pdb_ids = _load_bench_pdb_ids(ns.bench_json)
        seed_map = load_family_mapping(ns.seed)
        extras = [load_family_mapping(p) for p in (ns.extra_map or [])]
    except (ValueError, FileNotFoundError, json.JSONDecodeError) as exc:
        print(f"[family] ERROR: {exc}", file=sys.stderr)
        return 2

    merged = _merge_maps(seed_map, extras)
    # Guard: any user-supplied label that is not canonical (except
    # __unassigned__) → warn but keep the mapping.
    canonical = set(FAMILY_LABELS)
    non_canonical = {p: fam for p, fam in merged.items() if fam not in canonical}
    if non_canonical:
        print(
            f"[family] WARN: {len(non_canonical)} PDBs use non-canonical labels: "
            + ", ".join(f"{p}={f}" for p, f in list(non_canonical.items())[:6])
            + ("…" if len(non_canonical) > 6 else ""),
            file=sys.stderr,
        )

    report = build_coverage_report(bench_pdb_ids, merged)
    print(format_report(report))

    if ns.emit is not None:
        ns.emit.parent.mkdir(parents=True, exist_ok=True)
        payload: dict[str, str] = {p: merged[p] for p in bench_pdb_ids if p in merged}
        if ns.fill_unassigned:
            for p in bench_pdb_ids:
                payload.setdefault(p, "__unassigned__")
        ns.emit.write_text(
            json.dumps(payload, indent=2, sort_keys=True) + "\n"
        )
        print(f"[family] wrote merged map -> {ns.emit}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
