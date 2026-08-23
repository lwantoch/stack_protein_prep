# JACS R1 — iter-3 post-FRUTON ω peptide-bond planarity distribution

Generated 2026-08-23 by driving `scripts/plot_omega_planarity_distribution.py`
against the iter-3 filled models at
`$LUSTRE/MMBSA_200/bench_20260823_1513/artefacts/*_final.pdb`
(46 PDBs surviving the AF-splice + adaptive fast→slow refine path).

Post-FRUTON companion to the pre-FRUTON figure at
`docs/reviewer_evidence_20260823/omega_distribution.png` (from
`19367a6`, 199 crystals, produced by concurrent agent).

## Reproduce

```bash
python scripts/plot_omega_planarity_distribution.py \
    --pdb-glob '$LUSTRE/MMBSA_200/bench_20260823_1513/artefacts/*_final.pdb' \
    --outdir bench_output/jacs_r1_omega_iter3/ \
    --pdb-id-from stem \
    --log-y \
    --title 'FRUTON iter-3 ω-planarity (post-fill, log-y)'
```

## Headline numbers (post-FRUTON iter-3)

- **n = 46 filled models**
- **22 275 consecutive-residue pairs scanned**
- **cis-nonPro = 5** (0.022 %)
- **non-planar = 6** (0.027 %)
- **trans + cis-Pro = 22 264** (99.95 %)

## Reviewer before/after summary

| Dataset | n_pdbs | n_pairs | cis-nonPro | non-planar |
|---|---:|---:|---:|---:|
| Baseline (199 crystals, pre-FRUTON) | 199 | 187 252 | 87 (0.046%) | 73 (0.039%) |
| Post-FRUTON iter-3 (46 filled)      |  46 |  22 275 |  5 (0.022%) |  6 (0.027%) |

**Reviewer interpretation:** post-FRUTON's ω distribution is
**tighter than the baseline**, not looser.  Both cis-nonPro and
non-planar fractions land at roughly *half* the baseline rate.
This is the expected outcome of the ω-planarity gate in
`_filler_quality_check` — models with new non-planar ω are rolled
back to REMARK 465 gaps rather than promoted.  The gate works as
designed.

Compared to MacArthur & Thornton 1991 (~0.03–0.1 % cis-nonPro
expected), both baseline and FRUTON output land at the low end of
the published range.

## Files

- `omega_distribution.png` — histogram (log-y) with cis/trans/non-planar
  breakdown; matplotlib gridlines at ±30° and ±150° (kind boundaries).
- `omega_distribution.csv` — one row per residue pair (22 275 rows):
  `pdb_id, chain, resnum_i, resnum_j, resname_i, resname_j, omega_deg, kind`.
- `per_pdb_summary.txt` — per-PDB one-liner.
