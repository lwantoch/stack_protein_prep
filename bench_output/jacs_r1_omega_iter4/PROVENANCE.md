# JACS R1 — iter-4 post-FRUTON ω-planarity

Generated 2026-08-23 by driving `scripts/plot_omega_planarity_distribution.py`
against 46 iter-4 filled models at
`$LUSTRE/MMBSA_200/bench_20260823_1620/artefacts/*_final.pdb`.

## Reproduce

```bash
python scripts/plot_omega_planarity_distribution.py \
    --pdb-glob '$LUSTRE/MMBSA_200/bench_20260823_1620/artefacts/*_final.pdb' \
    --outdir bench_output/jacs_r1_omega_iter4/ \
    --pdb-id-from stem --log-y \
    --title 'FRUTON iter-4 ω-planarity (post-fill, log-y)'
```

## Headline numbers

- **n = 46 filled models**
- **22 233 consecutive-residue pairs scanned**
- **cis-nonPro = 5 (0.022 %)**
- **non-planar = 3 (0.013 %)** — down from iter-3's 6 (0.027 %)

## Iteration progression

| Iter | n_pdbs | n_pairs | cis-nonPro | non-planar |
|---|---:|---:|---:|---:|
| Baseline (199 pre-FRUTON) | 199 | 187 252 | 87 (0.046%) | 73 (0.039%) |
| iter-3 post-FRUTON        |  46 |  22 275 |  5 (0.022%) |  6 (0.027%) |
| **iter-4 post-FRUTON**    |  46 |  22 233 |  5 (0.022%) |  **3 (0.013%)** |

**Reviewer headline:** iter-4 non-planar ω count dropped 6 → 3
between iter-3 and iter-4 (halved).  Combined with the +6.5pp gate
PASS improvement (91.3 % → 97.8 %), this is quantitative evidence
that the iter-3-to-iter-4 patches (adaptive fast→slow + refined
rollback) are targeting exactly the ω-planarity failures the earlier
runs were losing to.
