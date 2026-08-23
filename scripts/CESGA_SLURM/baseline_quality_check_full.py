"""Baseline crystal-quality-check on ALL FRUTON-NEW crystals (208).

No MODELLER, no AF splice — just check_model_quality on the raw crystal.
Gives the pre-FRUTON metric distribution across the full set so we can
compare against FRUTON output + PDB-REDO downstream.

BENCH_OUT_DIR + BENCH_INDEX / --index behave the same as
fruton_bench_mmbsa200_full.py.  Emits ``<PDB>.json`` per protein.
"""
import argparse, os, sys, tempfile, time, json
from pathlib import Path
sys.path.insert(0, '/mnt/netapp2/Store_uni/home/otras/hcx/lwa/repos/FRUTON/stack_protein_prep/src')
from stack_protein_preparation._filler_quality_check import check_model_quality


CRYSTAL_ROOT = Path('/mnt/netapp1/Store_othcxlwa/FRUTON-NEW')


def _resolve_out_dir() -> Path:
    out_env = os.environ.get("BENCH_OUT_DIR")
    if out_env:
        d = Path(out_env); d.mkdir(parents=True, exist_ok=True); return d
    return Path(tempfile.mkdtemp(prefix="baseline_qc_"))


def _list_crystals() -> list[Path]:
    return sorted(p for p in CRYSTAL_ROOT.glob('*/*.pdb')
                  if p.stem == p.parent.name)


def _process(pdb_path: Path) -> dict:
    pdb = pdb_path.stem
    try:
        t0 = time.time()
        q = check_model_quality(pdb_path)
        dt = time.time() - t0
        return {
            'pdb': pdb,
            'n_residues': q.n_residues,
            'n_peptide_bonds_broken': q.n_peptide_bonds_broken,
            'n_clash_pairs': q.n_clash_pairs,
            'clash': q.n_clash_pairs,  # alias for downstream tools
            'n_omega_checked': q.n_omega_checked,
            'n_omega_non_planar': q.n_omega_non_planar,
            'n_omega_cis_nonpro': q.n_omega_cis_nonpro,
            'n_omega_cis_pro': q.n_omega_cis_pro,
            'n_rama_favoured': q.n_rama_favoured,
            'n_rama_allowed': q.n_rama_allowed,
            'n_rama_outlier': q.n_rama_outlier,
            'qc_seconds': dt,
        }
    except Exception as e:
        return {'pdb': pdb, 'error': f"{type(e).__name__}: {e}"}


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--index", type=int, default=None,
                    help="1-based index of a single PDB (SLURM array).")
    args = ap.parse_args()

    out_dir = _resolve_out_dir()
    crystals = _list_crystals()
    print(f"out_dir: {out_dir}\ncrystals: {len(crystals)}", flush=True)

    idx_env = os.environ.get("BENCH_INDEX")
    single = args.index or (int(idx_env) if idx_env else None)

    if single is not None:
        if not (1 <= single <= len(crystals)):
            print(f"[qc] index {single} out of range 1..{len(crystals)}", flush=True)
            return 2
        p = crystals[single - 1]
        row = _process(p)
        (out_dir / f"{row['pdb']}.json").write_text(json.dumps(row, indent=2, default=str))
        print(f"[qc] {row['pdb']}: "
              f"n={row.get('n_residues','?')} brk={row.get('n_peptide_bonds_broken','?')} "
              f"clash={row.get('n_clash_pairs','?')} ω_np={row.get('n_omega_non_planar','?')} "
              f"rama_out={row.get('n_rama_outlier','?')}", flush=True)
        return 0 if 'error' not in row else 1

    # Serial mode
    rows = []
    for p in crystals:
        row = _process(p)
        rows.append(row)
        print(f"{row['pdb']:<5} n={row.get('n_residues','?'):>4} clash={row.get('n_clash_pairs','?'):>4}", flush=True)
    (out_dir / 'baseline_quality_check_results.json').write_text(json.dumps(rows, indent=2, default=str))
    print(f"[qc] results: {out_dir / 'baseline_quality_check_results.json'}", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
