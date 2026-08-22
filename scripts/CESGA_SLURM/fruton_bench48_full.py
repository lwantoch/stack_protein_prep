"""Full-pipeline FRUTON benchmark: splice → refine → rollback on all 48 AF-
available proteins in /mnt/netapp1/Store_othcxlwa/FRUTON-NEW/.
"""
import sys, tempfile, time, json, signal
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


af_files = sorted(Path('/mnt/netapp1/Store_othcxlwa/FRUTON-NEW').glob('*/fasta/alignments/filler/*/alphafold/alphafold_aligned_model.pdb'))
tmp = Path(tempfile.mkdtemp(prefix='fruton_bench48_full_'))
print(f"tmpdir: {tmp}\nAF-available: {len(af_files)}", flush=True)

header = f"{'PDB':<5} {'gaps':>4} {'base_n':>6} {'splc_n':>6} {'refn_n':>6} {'fin_n':>6} {'Δn':>4} {'brk':>3} {'clash':>5} {'rolled':>6} {'rfn_s':>6} {'gate':>5}"
print(header, flush=True); print('-' * len(header), flush=True)

results = []
for i, af in enumerate(af_files, 1):
    pdb = af.parts[-7]
    C = Path(f'/mnt/netapp1/Store_othcxlwa/FRUTON-NEW/{pdb}/{pdb}.pdb')
    if not C.exists():
        print(f"{pdb:<5} NO_CRYS", flush=True); continue
    try:
        sp = tmp / f'{pdb}_spliced.pdb'; ref = tmp / f'{pdb}_refined.pdb'; fin = tmp / f'{pdb}_final.pdb'
        splice_af_gaps_into_crystal(C, af, sp, enable_rollback=False)
        cs = PDBParser(QUIET=True).get_structure('c', str(C))
        As = PDBParser(QUIET=True).get_structure('a', str(af))
        cmap = _protein_residue_map(cs); amap = _protein_residue_map(As)
        gaps = []
        for cid, ab in amap.items():
            cb = cmap.get(cid, {})
            for lo, hi in _detect_missing_windows(set(cb.keys()), sorted(ab.keys())):
                gaps.append((cid, lo, hi))

        # Per-protein 10-min timeout
        signal.signal(signal.SIGALRM, _timeout); signal.alarm(600)
        try:
            t0 = time.time()
            _r = refine_loops_via_modeller(input_pdb_path=sp, output_pdb_path=ref,
                                            gap_ranges_by_chain=gaps,
                                            n_conformers=3, refine_level='fast',
                                            reject_new_chirality_d=True)
            refine_time = time.time() - t0
        except TimeoutError_:
            refine_time = 600.0
            ref.write_bytes(sp.read_bytes())
            _r = None
        finally:
            signal.alarm(0)

        _, rolled = rollback_bad_gap_fills(ref, fin, gaps)
        b = check_model_quality(C); s = check_model_quality(sp); rf = check_model_quality(ref); f = check_model_quality(fin)
        p, reasons = f.passes_relative_gate(b)
        row = {
            'pdb': pdb, 'n_gaps': len(gaps),
            'base_n': b.n_residues, 'splice_n': s.n_residues, 'refine_n': rf.n_residues, 'final_n': f.n_residues,
            'delta_n': f.n_residues - b.n_residues, 'brk': f.n_peptide_bonds_broken, 'clash': f.n_clash_pairs,
            'rolled_gaps': len(rolled), 'refine_seconds': refine_time,
            'gate_pass': p, 'reasons': reasons,
        }
        results.append(row)
        print(f"{pdb:<5} {len(gaps):>4} {b.n_residues:>6} {s.n_residues:>6} {rf.n_residues:>6} {f.n_residues:>6} "
              f"{f.n_residues-b.n_residues:>+4d} {f.n_peptide_bonds_broken:>3} {f.n_clash_pairs:>5} "
              f"{len(rolled):>6} {refine_time:>6.0f} {'PASS' if p else 'FAIL':>5}", flush=True)
    except Exception as e:
        print(f"{pdb:<5} ERR: {type(e).__name__}: {str(e)[:80]}", flush=True)
        results.append({'pdb': pdb, 'error': f"{type(e).__name__}: {e}"})

    if i % 5 == 0:
        Path('/tmp/fruton_bench48_full_partial.json').write_text(json.dumps(results, indent=2, default=str))

n_pass = sum(1 for r in results if r.get('gate_pass'))
n_ok = len([r for r in results if 'error' not in r])
total_fills = sum(r.get('delta_n', 0) for r in results if 'delta_n' in r)
print(f"\n=== SUMMARY ===\nPASS: {n_pass}/{n_ok}, Total fills: {total_fills}", flush=True)

Path('/tmp/fruton_bench48_full_results.json').write_text(json.dumps(results, indent=2, default=str))
print(f"Results: /tmp/fruton_bench48_full_results.json", flush=True)
