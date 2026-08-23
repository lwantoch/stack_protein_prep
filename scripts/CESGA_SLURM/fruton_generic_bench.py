"""Generic FRUTON bench driver — no hardcoded protein set.

USER MANDATE 2026-08-23: 'keine spezifische Lösungen, immer allgemein.'
Same pipeline runs on train (MMBSA_75), val (MMBSA_125), holdout (+30
external) via `--split` — no per-set branching, no per-set tuning.

Emits per-protein:
    <BENCH_OUT_DIR>/<PDB>.json   — full stage-by-stage record + audit
    <BENCH_OUT_DIR>/<PDB>.frcmod — (future) metal frcmod when produced
And at the end (aggregator step):
    <BENCH_OUT_DIR>/audit_report.csv   — the primary user deliverable
    <BENCH_OUT_DIR>/audit_summary.md   — tier breakdown + rejected list

Two invocation modes:

    Per-index (SLURM array):
        BENCH_OUT_DIR=... BENCH_INDEX=$SLURM_ARRAY_TASK_ID \\
            python fruton_generic_bench.py --split train

    Serial (local):
        BENCH_OUT_DIR=... python fruton_generic_bench.py --split train
"""
from __future__ import annotations

import argparse
import json
import os
import signal
import sys
import tempfile
import time
from pathlib import Path

# Package on PYTHONPATH
sys.path.insert(0, '/mnt/netapp2/Store_uni/home/otras/hcx/lwa/repos/FRUTON/stack_protein_prep/src')

from Bio.PDB import PDBParser

from stack_protein_preparation._adaptive_refine_policy import decide_refine_action
from stack_protein_preparation._bench_split import load_split, filter_to_af_ready
from stack_protein_preparation._component_confidence import (
    ComponentConfidence,
    Confidence,
    ProteinDeliveryReport,
)
from stack_protein_preparation._filler_af_splice import (
    _detect_missing_windows,
    _protein_residue_map,
    rollback_bad_gap_fills,
    splice_af_gaps_into_crystal,
)
from stack_protein_preparation._filler_loop_refine import refine_loops_via_modeller
from stack_protein_preparation._filler_quality_check import check_model_quality


FRUTON_ROOT_DEFAULT = Path("/mnt/netapp1/Store_othcxlwa/FRUTON-NEW")
PER_PROTEIN_TIMEOUT_SECONDS = 900


class TimeoutError_(Exception): ...
def _timeout_handler(signum, frame): raise TimeoutError_("per-protein timeout")


def _resolve_paths(pdb_id: str, fruton_root: Path) -> tuple[Path, Path | None]:
    """Return (crystal_pdb, af_alignment_pdb-or-None)."""
    crystal = fruton_root / pdb_id / f"{pdb_id}.pdb"
    if not crystal.is_file():
        # Try lowercase / different casings sometimes seen on disk
        for variant in (pdb_id.lower(), pdb_id.upper()):
            alt = fruton_root / variant / f"{variant}.pdb"
            if alt.is_file():
                crystal = alt
                break
    af_matches = list(fruton_root.glob(
        f"{pdb_id}/fasta/alignments/filler/*/alphafold/alphafold_aligned_model.pdb"
    ))
    return crystal, af_matches[0] if af_matches else None


def _collect_gap_confidence(
    n_gaps: int, n_rolled: int, delta_n: int,
) -> list[ComponentConfidence]:
    """Turn aggregate gap-fill counters into per-gap ComponentConfidence records.

    Per-gap identity isn't tracked in the current pipeline output — this
    coarse aggregation gives one summary record so the audit CSV isn't
    empty for the gap component.  A follow-up wiring pass in
    _filler_alphafold + _filler_loop_refine will emit one record per
    individual gap with its pLDDT + refine level + geometry check.
    """
    if n_gaps == 0:
        return []
    if n_rolled == n_gaps and delta_n == 0:
        # All gaps rolled back to REMARK 465.  Reviewer-relevant point: the
        # SHIPPED model has NO broken residues — the pipeline was honest and
        # dropped every gap it could not close cleanly.  That is a valid
        # 'wahrscheinlich_ok' outcome (Confidence.MEDIUM), NOT a failure.
        # User just knows: this model has N gap windows still as REMARK 465
        # and may need a different AF template or manual bridging to fill
        # them.  Distinct from actual quality regressions.
        return [ComponentConfidence(
            component_type="gap_fill",
            name=f"aggregate ({n_gaps} gaps)",
            confidence=Confidence.MEDIUM,
            reason=(
                f"all {n_gaps} gap window(s) rolled back to REMARK 465; "
                f"model is geometrically clean but shorter than AF-target"
            ),
            suggested_action=(
                "if the missing residues matter for downstream MD, try a "
                "different AF alignment template or manual loop-bridging; "
                "otherwise the shipped shortened model is MD-ready"
            ),
            method="rollback_all",
        )]
    if n_rolled > 0:
        return [ComponentConfidence(
            component_type="gap_fill",
            name=f"aggregate ({n_gaps} gaps, {n_rolled} rolled)",
            confidence=Confidence.MEDIUM,
            reason=f"{n_rolled}/{n_gaps} gaps rolled back; {delta_n} residues rescued",
            suggested_action=f"review the {n_rolled} rolled-back windows in the pipeline log",
            method="MODELLER_partial",
        )]
    return [ComponentConfidence(
        component_type="gap_fill",
        name=f"aggregate ({n_gaps} gaps)",
        confidence=Confidence.HIGH,
        reason=f"{n_gaps} gap(s) rescued cleanly, {delta_n} residues gained",
        method="MODELLER_full",
    )]


def _collect_quality_confidence(
    baseline_qc, final_qc,
) -> list[ComponentConfidence]:
    """One record summarising the final model's quality regression."""
    dc = final_qc.n_clash_pairs - baseline_qc.n_clash_pairs
    dω = final_qc.n_omega_non_planar - baseline_qc.n_omega_non_planar
    dbrk = final_qc.n_peptide_bonds_broken - baseline_qc.n_peptide_bonds_broken
    if dc <= 0 and dω <= 0 and dbrk <= 0:
        return [ComponentConfidence(
            component_type="gap_fill",
            name="global quality vs crystal",
            confidence=Confidence.HIGH,
            reason="fruton model did not regress any tracked quality metric",
            method="check_model_quality",
            details={"dclash": dc, "domega_np": dω, "dbroken": dbrk},
        )]
    regressed = []
    if dc > 0: regressed.append(f"clash +{dc}")
    if dω > 0: regressed.append(f"ω_np +{dω}")
    if dbrk > 0: regressed.append(f"broken +{dbrk}")
    return [ComponentConfidence(
        component_type="gap_fill",
        name="global quality vs crystal",
        confidence=Confidence.MEDIUM if (dc <= 5 and dω <= 1 and dbrk == 0) else Confidence.LOW,
        reason=f"regressed vs crystal baseline: {'; '.join(regressed)}",
        suggested_action=(
            "inspect the flagged residues in the pipeline log; if "
            "reviewer-critical, consider manual refinement of the affected "
            "loop or metal coordination sphere"
        ),
        method="check_model_quality",
        details={"dclash": dc, "domega_np": dω, "dbroken": dbrk},
    )]


def _process_one(pdb_id: str, tmp: Path, fruton_root: Path) -> dict:
    crystal, af = _resolve_paths(pdb_id, fruton_root)
    if not crystal.is_file():
        return {
            "pdb": pdb_id,
            "error": "NO_CRYSTAL",
            "delivery": ProteinDeliveryReport(pdb=pdb_id, model_written=False).to_dict(),
        }
    if af is None:
        # We can still produce a "delivered" model — just no gap fill possible.
        # For now, mark as skipped so the aggregator distinguishes.
        return {
            "pdb": pdb_id,
            "note": "NO_AF_ALIGNMENT",
            "delivery": ProteinDeliveryReport(
                pdb=pdb_id, model_written=False,
                notes="no AF alignment available; gap fill skipped",
            ).to_dict(),
        }

    components: list[ComponentConfidence] = []
    try:
        sp = tmp / f"{pdb_id}_spliced.pdb"
        ref = tmp / f"{pdb_id}_refined.pdb"
        fin = tmp / f"{pdb_id}_final.pdb"

        splice_af_gaps_into_crystal(crystal, af, sp, enable_rollback=False)

        cs = PDBParser(QUIET=True).get_structure("c", str(crystal))
        As = PDBParser(QUIET=True).get_structure("a", str(af))
        cmap = _protein_residue_map(cs); amap = _protein_residue_map(As)
        gaps = []
        for cid, ab in amap.items():
            cb = cmap.get(cid, {})
            for lo, hi in _detect_missing_windows(set(cb.keys()), sorted(ab.keys())):
                gaps.append((cid, lo, hi))

        signal.signal(signal.SIGALRM, _timeout_handler)
        signal.alarm(PER_PROTEIN_TIMEOUT_SECONDS)
        try:
            t0 = time.time()
            refine_loops_via_modeller(
                input_pdb_path=sp, output_pdb_path=ref,
                gap_ranges_by_chain=gaps,
                n_conformers=3, refine_level="fast",
                reject_new_chirality_d=True,
            )
            # Principle-driven adaptive slow retry (regression-based, no
            # magic thresholds).  Debug-print the decision so the SLURM
            # log shows whether the trigger fired for a given PDB.
            try:
                b_qc = check_model_quality(crystal)
                f_qc_fast = check_model_quality(ref)
                decision = decide_refine_action(b_qc, f_qc_fast)
                print(
                    f"    refine-policy ({pdb_id}): {decision.action} "
                    f"regressed={sorted(decision.regressed_metrics)} "
                    f"ceiling={sorted(decision.exceeded_ceilings)}",
                    flush=True,
                )
                if decision.action == "retry_slow":
                    print(f"    retry_slow ({pdb_id}): running slow-refine 5 conformers",
                          flush=True)
                    refine_loops_via_modeller(
                        input_pdb_path=sp, output_pdb_path=ref,
                        gap_ranges_by_chain=gaps,
                        n_conformers=5, refine_level="slow",
                        reject_new_chirality_d=True,
                    )
                    print(f"    retry_slow ({pdb_id}): slow-refine done", flush=True)
            except Exception as _e:
                # Do NOT swallow silently — reviewer-defensibility requires
                # the log to say what went wrong when the policy failed.
                import traceback
                print(f"    refine-policy FAILED ({pdb_id}): {type(_e).__name__}: {_e}",
                      flush=True)
                print(traceback.format_exc(), flush=True)
            refine_time = time.time() - t0
        except TimeoutError_:
            refine_time = float(PER_PROTEIN_TIMEOUT_SECONDS)
            ref.write_bytes(sp.read_bytes())
        finally:
            signal.alarm(0)

        _, rolled = rollback_bad_gap_fills(ref, fin, gaps)
        b_qc = check_model_quality(crystal)
        f_qc = check_model_quality(fin)
        gate_pass, reasons = f_qc.passes_relative_gate(b_qc)

        # Emit component confidence records
        components.extend(_collect_gap_confidence(
            n_gaps=len(gaps),
            n_rolled=len(rolled),
            delta_n=f_qc.n_residues - b_qc.n_residues,
        ))
        components.extend(_collect_quality_confidence(b_qc, f_qc))

        # Assemble delivery report
        model_written = fin.is_file()
        # A gate FAIL is NOT "not_delivered" — it's the model + component notes.
        # not_delivered is only when we couldn't produce a file at all.
        delivery = ProteinDeliveryReport(
            pdb=pdb_id,
            components=components,
            model_written=model_written,
            notes=("; ".join(reasons)[:200]) if reasons else "",
        )

        return {
            "pdb": pdb_id,
            "n_gaps": len(gaps),
            "base_n": b_qc.n_residues,
            "final_n": f_qc.n_residues,
            "delta_n": f_qc.n_residues - b_qc.n_residues,
            "brk": f_qc.n_peptide_bonds_broken,
            "clash": f_qc.n_clash_pairs,
            "n_omega_non_planar": f_qc.n_omega_non_planar,
            "n_omega_cis_nonpro": f_qc.n_omega_cis_nonpro,
            "n_rama_outlier": f_qc.n_rama_outlier,
            "n_clash_pairs": f_qc.n_clash_pairs,
            "n_peptide_bonds_broken": f_qc.n_peptide_bonds_broken,
            "rolled_gaps": len(rolled),
            "refine_seconds": refine_time,
            "gate_pass": gate_pass,
            "reasons": reasons,
            "delivery": delivery.to_dict(),
        }
    except Exception as e:
        return {
            "pdb": pdb_id,
            "error": f"{type(e).__name__}: {str(e)[:200]}",
            "delivery": ProteinDeliveryReport(
                pdb=pdb_id, model_written=False,
                notes=f"pipeline exception: {type(e).__name__}",
            ).to_dict(),
        }


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--split", required=True, choices=("train", "val", "holdout"))
    ap.add_argument("--index", type=int, default=None,
                    help="1-based index into the AF-ready subset (SLURM array).")
    ap.add_argument("--fruton-root", type=Path, default=FRUTON_ROOT_DEFAULT)
    args = ap.parse_args()

    out_env = os.environ.get("BENCH_OUT_DIR")
    if out_env:
        tmp = Path(out_env); tmp.mkdir(parents=True, exist_ok=True)
    else:
        tmp = Path(tempfile.mkdtemp(prefix=f"fruton_{args.split}_"))

    all_ids = load_split(args.split)
    ready = filter_to_af_ready(all_ids, args.fruton_root)
    print(f"[bench] split={args.split} total_in_split={len(all_ids)} af_ready={len(ready)} out={tmp}", flush=True)

    idx_env = os.environ.get("BENCH_INDEX")
    single = args.index or (int(idx_env) if idx_env else None)

    if single is not None:
        if not (1 <= single <= len(ready)):
            print(f"[bench] index {single} out of range 1..{len(ready)}", flush=True)
            return 2
        pid = ready[single - 1]
        row = _process_one(pid, tmp, args.fruton_root)
        (tmp / f"{row['pdb']}.json").write_text(json.dumps(row, indent=2, default=str))
        status = row.get("delivery", {}).get("overall_status", "?")
        print(f"[bench] {row['pdb']} → {status}", flush=True)
        return 0 if "error" not in row else 1

    rows = []
    for pid in ready:
        row = _process_one(pid, tmp, args.fruton_root)
        rows.append(row)
        status = row.get("delivery", {}).get("overall_status", "?")
        print(f"[bench] {row['pdb']} → {status}", flush=True)
        (tmp / f"{row['pdb']}.json").write_text(json.dumps(row, indent=2, default=str))
    (tmp / "combined_results.json").write_text(json.dumps(rows, indent=2, default=str))
    return 0


if __name__ == "__main__":
    sys.exit(main())
