#!/usr/bin/env python
"""Run the Meeko gate on every FRUTON-delivered receptor and report a table.

Compares against ``newbench_27/delivered/MEEKO_STATUS.csv`` to verify the
classifier's error_class matches the historical labels. A green run means
the gate reliably reproduces the observed failure modes on this bench and
is safe to wire into the pipeline.

Usage:
    pixi run python scripts/meeko_gate_regression.py \\
        --delivered-root /mnt/netapp1/Store_othcxlwa/newbench_27/delivered \\
        --status-csv    /mnt/netapp1/Store_othcxlwa/newbench_27/delivered/MEEKO_STATUS.csv
"""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

# Make the module import work from a repo checkout without install
HERE = Path(__file__).resolve().parent
SRC = HERE.parent / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from stack_protein_preparation._meeko_gate import (  # noqa: E402
    MeekoErrorClass,
    validate_pdb_for_meeko,
)


# Map historical csv err_class -> our normalized MeekoErrorClass
_HISTORICAL_TO_OURS: dict[str, MeekoErrorClass] = {
    "valence": "valence",
    "adjacent_smarts": "adjacent_smarts",
    "no_template": "no_template",
    "other": "other",  # historical "other" often == "padding" in practice
}


def _cmp(row_historical: str, row_ours: str | None, ok_ours: bool) -> str:
    """Return AGREE / DIFF / MISS marker for the diff column."""
    if ok_ours:
        return "AGREE" if row_historical == "" else "DIFF (we say ok, csv says fail)"
    if row_historical == "":
        return "DIFF (we say fail, csv says ok)"
    expected = _HISTORICAL_TO_OURS.get(row_historical, "other")
    # padding is a stricter partition of "other" — if we say padding and csv
    # says other, that's an IMPROVEMENT not a diff
    if expected == "other" and row_ours == "padding":
        return "AGREE+ (refined)"
    return "AGREE" if row_ours == expected else f"DIFF (ours={row_ours}, csv={row_historical})"


def main() -> int:
    p = argparse.ArgumentParser()
    p.add_argument("--delivered-root", type=Path, required=True)
    p.add_argument("--status-csv", type=Path, required=True)
    p.add_argument("--pdb-filename", default="receptor.pdb",
                   help="Which file inside <delivered>/<PDB>/ to check (default receptor.pdb)")
    args = p.parse_args()

    if not args.status_csv.is_file():
        print(f"status csv not found: {args.status_csv}", file=sys.stderr)
        return 2

    historical: dict[str, tuple[str, str]] = {}   # pdb_id -> (status, err_class)
    with args.status_csv.open() as fh:
        for row in csv.DictReader(fh):
            historical[row["pdb_id"]] = (row["meeko_status"], row.get("err_class", ""))

    print(f"{'PDB':6} {'status_ours':11} {'error_class':16} {'expected(csv)':16} {'verdict'}")
    print("-" * 80)
    agree = diff = missing = 0
    for pdb_id, (hist_status, hist_err) in sorted(historical.items()):
        pdb_path = args.delivered_root / pdb_id / args.pdb_filename
        if not pdb_path.is_file():
            print(f"{pdb_id:6} {'MISSING':11} {'-':16} {hist_err or '-':16} MISS")
            missing += 1
            continue
        res = validate_pdb_for_meeko(pdb_path)
        status_ours = "OK" if res.ok else "FAIL"
        verdict = _cmp(hist_err, res.error_class, res.ok)
        print(f"{pdb_id:6} {status_ours:11} {res.error_class or '-':16} {hist_err or '-':16} {verdict}")
        if verdict.startswith("AGREE"):
            agree += 1
        else:
            diff += 1

    total = len(historical)
    print("-" * 80)
    print(f"agree: {agree}/{total}   diff: {diff}   missing: {missing}")
    return 0 if diff == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
