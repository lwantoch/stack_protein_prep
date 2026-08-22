"""Publication-tier CSV summary for the 48-protein FRUTON benchmark.

Emits one row per protein with the reviewer-facing metrics from
``QualityReport``: PDB id, n_residues, delta (residues added), Rama
favoured/outlier %, MolProbity-style clashscore per 1000 heavy atoms,
broken peptide bonds, cis-nonPro + non-planar omega counts, D-chirality
outliers, and the relative-gate pass/fail decision.

Reads spliced+refined+rollback final PDBs from any directory matching
``final_model_dir/{PDB}_final.pdb`` and compares them to the raw crystal
under ``/mnt/netapp1/Store_othcxlwa/FRUTON-NEW/{PDB}/{PDB}.pdb``.  Falls
back to splice-only output when no final is available.

Usage:
    python fruton_bench48_summary_csv.py [--out /path/to/summary.csv]
                                          [--models-dir /path/to/finals/]
"""
from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent / "src"))
from stack_protein_preparation._filler_af_splice import (  # noqa: E402
    _detect_missing_windows,
    _protein_residue_map,
    splice_af_gaps_into_crystal,
)
from stack_protein_preparation._filler_quality_check import (  # noqa: E402
    check_model_quality,
)
from Bio.PDB import PDBParser  # noqa: E402


CRYSTAL_ROOT = Path("/mnt/netapp1/Store_othcxlwa/FRUTON-NEW")


def _find_af_aligned(pdb_id: str) -> Path | None:
    candidates = sorted(
        (CRYSTAL_ROOT / pdb_id).glob(
            "fasta/alignments/filler/*/alphafold/alphafold_aligned_model.pdb"
        )
    )
    return candidates[0] if candidates else None


def _generate_splice_only(pdb_id: str, out_dir: Path) -> Path | None:
    crystal = CRYSTAL_ROOT / pdb_id / f"{pdb_id}.pdb"
    af = _find_af_aligned(pdb_id)
    if not crystal.is_file() or af is None:
        return None
    out = out_dir / f"{pdb_id}_splice.pdb"
    splice_af_gaps_into_crystal(crystal, af, out)
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument(
        "--out",
        type=Path,
        default=Path(__file__).with_suffix(".csv"),
        help="Output CSV path (default: alongside this script)",
    )
    ap.add_argument(
        "--models-dir",
        type=Path,
        default=None,
        help="Directory of pre-computed *_final.pdb files (from full-pipeline "
        "bench). If missing, splice-only is regenerated per protein.",
    )
    args = ap.parse_args()

    proteins = sorted(
        p.name for p in CRYSTAL_ROOT.iterdir() if p.is_dir() and (p / f"{p.name}.pdb").is_file()
    )
    # Only include proteins whose FRUTON-NEW tree has an AF alignment.
    proteins = [p for p in proteins if _find_af_aligned(p) is not None]

    import tempfile
    with tempfile.TemporaryDirectory(prefix="fruton_csv_") as _tmp:
        tmp = Path(_tmp)
        rows: list[dict] = []
        for pdb_id in proteins:
            crystal = CRYSTAL_ROOT / pdb_id / f"{pdb_id}.pdb"
            baseline = check_model_quality(crystal)
            final_path = None
            if args.models_dir is not None:
                candidate = args.models_dir / f"{pdb_id}_final.pdb"
                if candidate.is_file():
                    final_path = candidate
            if final_path is None:
                final_path = _generate_splice_only(pdb_id, tmp)
            final = check_model_quality(final_path) if final_path else baseline
            passed, reasons = final.passes_relative_gate(baseline)
            rows.append({
                "pdb_id": pdb_id,
                "n_residues_baseline": baseline.n_residues,
                "n_residues_final": final.n_residues,
                "delta_residues": final.n_residues - baseline.n_residues,
                "rama_favoured_pct_baseline": round(baseline.rama_favoured_pct(), 2),
                "rama_favoured_pct_final": round(final.rama_favoured_pct(), 2),
                "rama_outlier_pct_baseline": round(baseline.rama_outlier_pct(), 3),
                "rama_outlier_pct_final": round(final.rama_outlier_pct(), 3),
                "clashscore_baseline": round(baseline.clashscore_per_1000_atoms(), 2),
                "clashscore_final": round(final.clashscore_per_1000_atoms(), 2),
                "broken_peptide_bonds_baseline": baseline.n_peptide_bonds_broken,
                "broken_peptide_bonds_final": final.n_peptide_bonds_broken,
                "omega_cis_nonpro_baseline": baseline.n_omega_cis_nonpro,
                "omega_cis_nonpro_final": final.n_omega_cis_nonpro,
                "omega_non_planar_baseline": baseline.n_omega_non_planar,
                "omega_non_planar_final": final.n_omega_non_planar,
                "d_chirality_baseline": baseline.n_ca_chirality_outliers,
                "d_chirality_final": final.n_ca_chirality_outliers,
                "gate_pass": passed,
                "fail_reasons": "; ".join(reasons) if reasons else "",
            })
            print(f"{pdb_id}: {final.one_line_summary()} {'PASS' if passed else 'FAIL'}", flush=True)

        args.out.parent.mkdir(parents=True, exist_ok=True)
        with args.out.open("w", newline="", encoding="utf-8") as fh:
            writer = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
            writer.writeheader()
            writer.writerows(rows)
        n_pass = sum(1 for r in rows if r["gate_pass"])
        total_delta = sum(r["delta_residues"] for r in rows)
        print(f"\n=== Summary ===")
        print(f"CSV: {args.out}")
        print(f"PASS: {n_pass}/{len(rows)}")
        print(f"Total residues filled: {total_delta}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
