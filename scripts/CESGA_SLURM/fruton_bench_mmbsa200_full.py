"""Full-pipeline FRUTON benchmark: splice → refine → rollback.

Runs one PDB per invocation when ``--index <1-based>`` (or ``BENCH_INDEX``
env var) is provided.  Without it, iterates over all 48 AF-available
proteins in /mnt/netapp1/Store_othcxlwa/FRUTON-NEW/.

Output paths are controlled by the ``BENCH_OUT_DIR`` env var.  When
set, per-protein artefacts and per-PDB result JSONs land there; the
serial mode also writes a combined ``fruton_bench_mmbsa200_results.json``
at the end.

Single-PDB mode (SLURM array style):

    BENCH_OUT_DIR=/lustre/.../artefacts BENCH_INDEX=$SLURM_ARRAY_TASK_ID \\
        python fruton_bench_mmbsa200_full.py

emits ``<BENCH_OUT_DIR>/<PDB>.json`` — no combined results file,
that is the aggregator's job.
"""
import argparse, os, sys, tempfile, time, json, signal
from pathlib import Path
sys.path.insert(0, '/mnt/netapp2/Store_uni/home/otras/hcx/lwa/repos/FRUTON/stack_protein_prep/src')
from stack_protein_preparation._filler_af_splice import (
    splice_af_gaps_into_crystal, rollback_bad_gap_fills,
    _detect_missing_windows, _protein_residue_map,
)
from stack_protein_preparation._filler_quality_check import check_model_quality
from stack_protein_preparation._filler_loop_refine import refine_loops_via_modeller
from Bio.PDB import PDBParser


class TimeoutError_(Exception): pass
def _timeout(signum, frame): raise TimeoutError_("per-protein timeout")


def _resolve_out_dir() -> tuple[Path, Path, Path]:
    """Return (tmpdir, results_json_path, partial_json_path)."""
    out_env = os.environ.get("BENCH_OUT_DIR")
    if out_env:
        d = Path(out_env)
        d.mkdir(parents=True, exist_ok=True)
        return d, d / "fruton_bench_mmbsa200_results.json", d / "fruton_bench_mmbsa200_partial.json"
    d = Path(tempfile.mkdtemp(prefix="fruton_bench_mmbsa200_"))
    return d, Path("/tmp/fruton_bench_mmbsa200_results.json"), Path("/tmp/fruton_bench_mmbsa200_partial.json")


def _process_one(af: Path, tmp: Path) -> dict:
    pdb = af.parts[-7]
    C = Path(f"/mnt/netapp1/Store_othcxlwa/FRUTON-NEW/{pdb}/{pdb}.pdb")
    if not C.exists():
        return {"pdb": pdb, "error": "NO_CRYS"}
    try:
        sp = tmp / f'{pdb}_spliced.pdb'
        ref = tmp / f'{pdb}_refined.pdb'
        fin = tmp / f'{pdb}_final.pdb'
        splice_af_gaps_into_crystal(C, af, sp, enable_rollback=False)
        cs = PDBParser(QUIET=True).get_structure('c', str(C))
        As = PDBParser(QUIET=True).get_structure('a', str(af))
        cmap = _protein_residue_map(cs); amap = _protein_residue_map(As)
        gaps = []
        for cid, ab in amap.items():
            cb = cmap.get(cid, {})
            for lo, hi in _detect_missing_windows(set(cb.keys()), sorted(ab.keys())):
                gaps.append((cid, lo, hi))

        signal.signal(signal.SIGALRM, _timeout); signal.alarm(900)
        try:
            t0 = time.time()
            refine_loops_via_modeller(
                input_pdb_path=sp, output_pdb_path=ref,
                gap_ranges_by_chain=gaps,
                n_conformers=3, refine_level='fast',
                reject_new_chirality_d=True,
            )
            # Adaptive fast → slow retry (mirrors _filler_alphafold.py trigger,
            # a75e077 2026-08-23).  Compare ω-non-planar + clash gain of the
            # fast-refined structure vs the crystal baseline; if either
            # threshold is exceeded, retry with slow + 5 conformers.
            #
            # 2026-08-23 iter-3 finding: 8Q68 rescued (ω_np_gain=24 → 0) at
            # threshold 3, but 4 singleton fails (4X7Q, 5U5T, 6PYR, 6Q7D)
            # all had ω_np_gain=1 and slipped below threshold.  Iter-4
            # lowers ω threshold to 1 so ANY new non-planar bond triggers
            # slow retry; 8Q68 still catches on clash_gain>20 too.
            try:
                _b_qc = check_model_quality(C)
                _f_qc = check_model_quality(ref)
                _clash_gain = _f_qc.n_clash_pairs - _b_qc.n_clash_pairs
                _omega_gain = _f_qc.n_omega_non_planar - _b_qc.n_omega_non_planar
                # Ceiling guard (iter-4 finding): 4AT5 clash_gain=582 and
                # 5HJS clash_gain=1000 both stalled MODELLER slow-refine
                # indefinitely (>1 h) with no productive convergence — the
                # splice was too broken for LoopModel to rescue.  Skip
                # slow retry in that regime and let downstream rollback
                # aggressively drop residues via REMARK 465 instead.
                if _clash_gain > 200:
                    print(f"    catastrophic clash_gain={_clash_gain} > 200 -- skipping slow retry, deferring to rollback ({pdb})", flush=True)
                elif _clash_gain > 20 or _omega_gain >= 1:
                    print(f"    adaptive-slow retry ({pdb}): clash_gain={_clash_gain} omega_gain={_omega_gain}", flush=True)
                    refine_loops_via_modeller(
                        input_pdb_path=sp, output_pdb_path=ref,
                        gap_ranges_by_chain=gaps,
                        n_conformers=5, refine_level='slow',
                        reject_new_chirality_d=True,
                    )
            except Exception as _exc:
                print(f"    adaptive-slow check failed ({pdb}): {_exc!r}", flush=True)
            refine_time = time.time() - t0
        except TimeoutError_:
            refine_time = 900.0
            ref.write_bytes(sp.read_bytes())
        finally:
            signal.alarm(0)

        _, rolled = rollback_bad_gap_fills(ref, fin, gaps)
        b = check_model_quality(C)
        s = check_model_quality(sp)
        rf = check_model_quality(ref)
        f = check_model_quality(fin)
        p, reasons = f.passes_relative_gate(b)
        return {
            'pdb': pdb, 'n_gaps': len(gaps),
            'base_n': b.n_residues, 'splice_n': s.n_residues,
            'refine_n': rf.n_residues, 'final_n': f.n_residues,
            'delta_n': f.n_residues - b.n_residues,
            'brk': f.n_peptide_bonds_broken,
            'clash': f.n_clash_pairs,
            'rolled_gaps': len(rolled),
            'refine_seconds': refine_time,
            'gate_pass': p, 'reasons': reasons,
            # extra fields for downstream figures + PDB-REDO comparison
            'n_omega_non_planar': f.n_omega_non_planar,
            'n_omega_cis_nonpro': f.n_omega_cis_nonpro,
            'n_rama_outlier': f.n_rama_outlier,
            'n_clash_pairs': f.n_clash_pairs,
        }
    except Exception as e:
        return {'pdb': pdb, 'error': f"{type(e).__name__}: {e}"}


def _header_line() -> str:
    return f"{'PDB':<5} {'gaps':>4} {'base_n':>6} {'splc_n':>6} {'refn_n':>6} {'fin_n':>6} {'Δn':>4} {'brk':>3} {'clash':>5} {'rolled':>6} {'rfn_s':>6} {'gate':>5}"


def _fmt_row(r: dict) -> str:
    if 'error' in r:
        return f"{r['pdb']:<5} ERR: {r['error'][:80]}"
    gate = 'PASS' if r.get('gate_pass') else 'FAIL'
    return (
        f"{r['pdb']:<5} {r.get('n_gaps', 0):>4} {r.get('base_n', 0):>6} "
        f"{r.get('splice_n', 0):>6} {r.get('refine_n', 0):>6} "
        f"{r.get('final_n', 0):>6} {r.get('delta_n', 0):>+4d} "
        f"{r.get('brk', 0):>3} {r.get('clash', 0):>5} "
        f"{r.get('rolled_gaps', 0):>6} {r.get('refine_seconds', 0.0):>6.0f} "
        f"{gate:>5}"
    )


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--index", type=int, default=None,
                    help="1-based index of a single PDB to process (SLURM array mode).")
    args = ap.parse_args()

    tmp, results_json, partial_json = _resolve_out_dir()
    af_files = sorted(Path('/mnt/netapp1/Store_othcxlwa/FRUTON-NEW').glob(
        '*/fasta/alignments/filler/*/alphafold/alphafold_aligned_model.pdb'
    ))
    print(f"tmpdir: {tmp}\nAF-available: {len(af_files)}\nresults: {results_json}", flush=True)

    idx_env = os.environ.get("BENCH_INDEX")
    single_index = args.index or (int(idx_env) if idx_env else None)

    if single_index is not None:
        if not (1 <= single_index <= len(af_files)):
            print(f"[bench] index {single_index} out of range 1..{len(af_files)}", flush=True)
            return 2
        af = af_files[single_index - 1]
        print(_header_line(), flush=True)
        print("-" * len(_header_line()), flush=True)
        row = _process_one(af, tmp)
        print(_fmt_row(row), flush=True)
        pid = row.get("pdb", f"idx{single_index}")
        (tmp / f"{pid}.json").write_text(json.dumps(row, indent=2, default=str))
        print(f"[bench] wrote {tmp / f'{pid}.json'}", flush=True)
        return 0 if 'error' not in row else 1

    # Serial mode: iterate over all PDBs (original behaviour).
    print(_header_line(), flush=True)
    print("-" * len(_header_line()), flush=True)
    results = []
    for i, af in enumerate(af_files, 1):
        row = _process_one(af, tmp)
        results.append(row)
        print(_fmt_row(row), flush=True)
        if i % 5 == 0:
            partial_json.write_text(json.dumps(results, indent=2, default=str))

    n_pass = sum(1 for r in results if r.get('gate_pass'))
    n_ok = len([r for r in results if 'error' not in r])
    total_fills = sum(r.get('delta_n', 0) for r in results if 'delta_n' in r)
    print(f"\n=== SUMMARY ===\nPASS: {n_pass}/{n_ok}, Total fills: {total_fills}", flush=True)
    results_json.write_text(json.dumps(results, indent=2, default=str))
    print(f"Results: {results_json}", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
