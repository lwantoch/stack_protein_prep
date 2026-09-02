#!/usr/bin/env python
"""End-to-end regression: apply the full protonation.py retry chain to every
delivered receptor and report the final verdict per PDB.

For each receptor.pdb we run:
  1. strict meeko gate
  2. if error_class == "valence" → openmm H-repair → re-run strict gate
  3. if still failing → tolerant gate (allow_bad_res=True) → "degraded_pass"
     or "fail"

This measures the true rescue rate of the changes wired into
``protonate_selected_structure`` without needing to re-run the whole pipeline
on every target.

NOTE: the input files here are FRUTON-delivered receptor.pdb outputs from a
previous pipeline run — that is what protonation.py itself writes out. The
retry chain does NOT re-run pdb2gmx; it only replaces H atoms and re-checks
the gate. If a receptor had a defect that pdb2gmx introduced and it needs a
full re-protonation, that will not be caught here — but the openmm-only
retry is the same one wired into protonation.py so what we measure IS what
will happen in production.

Usage:
    pixi run python scripts/meeko_retry_chain_regression.py \\
        --delivered-root /mnt/netapp1/Store_othcxlwa/newbench_27/delivered \\
        --status-csv    /mnt/netapp1/Store_othcxlwa/newbench_27/delivered/MEEKO_STATUS.csv \\
        [--work-copy]      # copy each pdb to a tmpdir before repair (recommended)
"""

from __future__ import annotations

import argparse
import csv
import shutil
import sys
import tempfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
SRC = HERE.parent / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from stack_protein_preparation._meeko_gate import validate_pdb_for_meeko  # noqa: E402
from stack_protein_preparation._openmm_h_repair import repair_hydrogens_openmm  # noqa: E402


def _classify_and_retry(pdb_path: Path, *, ph: float) -> tuple[str, str]:
    """Return (verdict, detail). verdict in
    {strict_pass, repaired_strict_pass, degraded_pass, fail, gate_error}.
    """
    pre = validate_pdb_for_meeko(pdb_path)
    if pre.ok:
        return "strict_pass", ""
    if pre.error_class == "valence":
        r = repair_hydrogens_openmm(pdb_path, ph=ph)
        if r.ok:
            post = validate_pdb_for_meeko(pdb_path)
            if post.ok:
                return "repaired_strict_pass", f"H_added={r.n_hydrogens_added}"
    deg = validate_pdb_for_meeko(pdb_path, allow_bad_res=True)
    if deg.ok:
        return "degraded_pass", f"strict_err={pre.error_class}"
    return "fail", f"{pre.error_class or 'other'}: {pre.message[:60]}"


def main() -> int:
    p = argparse.ArgumentParser()
    p.add_argument("--delivered-root", type=Path, required=True)
    p.add_argument("--status-csv", type=Path, required=True)
    p.add_argument("--pdb-filename", default="receptor.pdb")
    p.add_argument("--ph", type=float, default=7.4)
    p.add_argument("--work-copy", action="store_true",
                   help="Copy each pdb to a tmpdir before repair "
                        "(prevents in-place modification of the delivered file)")
    args = p.parse_args()

    if not args.status_csv.is_file():
        print(f"status csv not found: {args.status_csv}", file=sys.stderr)
        return 2

    historical: dict[str, tuple[str, str]] = {}
    with args.status_csv.open() as fh:
        for row in csv.DictReader(fh):
            historical[row["pdb_id"]] = (row["meeko_status"], row.get("err_class", ""))

    verdicts: dict[str, int] = {}
    print(f"{'PDB':6} {'verdict':22} {'hist_err':16} detail")
    print("-" * 100)

    with tempfile.TemporaryDirectory(prefix="meeko_retry_") as tmpdir:
        tmp_root = Path(tmpdir)
        for pdb_id, (hist_status, hist_err) in sorted(historical.items()):
            src_path = args.delivered_root / pdb_id / args.pdb_filename
            if not src_path.is_file():
                print(f"{pdb_id:6} {'MISSING':22} {hist_err or '-':16} file not found")
                verdicts["missing"] = verdicts.get("missing", 0) + 1
                continue
            if args.work_copy:
                work_path = tmp_root / f"{pdb_id}.pdb"
                shutil.copy(src_path, work_path)
            else:
                work_path = src_path
            verdict, detail = _classify_and_retry(work_path, ph=args.ph)
            verdicts[verdict] = verdicts.get(verdict, 0) + 1
            print(f"{pdb_id:6} {verdict:22} {hist_err or '-':16} {detail}")

    total = sum(verdicts.values())
    print("-" * 100)
    for v in ("strict_pass", "repaired_strict_pass", "degraded_pass",
              "fail", "gate_error", "missing"):
        n = verdicts.get(v, 0)
        pct = 100 * n / total if total else 0
        print(f"  {v:22} {n:3d} / {total} ({pct:.0f}%)")
    usable = (verdicts.get("strict_pass", 0)
              + verdicts.get("repaired_strict_pass", 0)
              + verdicts.get("degraded_pass", 0))
    print(f"  TOTAL usable         {usable:3d} / {total} ({100*usable/total if total else 0:.0f}%)")
    return 0 if verdicts.get("fail", 0) + verdicts.get("gate_error", 0) == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
