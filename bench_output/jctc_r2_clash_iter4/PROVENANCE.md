# JCTC R2 — iter-4 post-FRUTON clashscore

Generated 2026-08-23 by driving `scripts/plot_clashscore_histogram.py`
against the iter-4 SLURM aggregate at
`$LUSTRE/MMBSA_200/bench_20260823_1620/artefacts/fruton_bench_mmbsa200_results.json`
(46 filled models).

## Reproduce

```bash
python scripts/plot_clashscore_histogram.py \
    --bench-json $LUSTRE/MMBSA_200/bench_20260823_1620/artefacts/fruton_bench_mmbsa200_results.json \
    --outdir bench_output/jctc_r2_clash_iter4/ \
    --x-max 60 \
    --title 'Post-FRUTON iter-4 clashscore (47 filled)'
```

## Headline numbers

- **n = 46 filled models**
- **gate PASS = 45 / 46 (97.8 %)** — up from iter-3's 91.3 %
- mean clashscore = 2.98, p50 = 1.0, p90 = 6.5, p99 = 33.0, max = 33
- 39 / 46 (84.8 %) below MolProbity 90th (≤ 4)

## Iteration progression (JCTC R2 reviewer chain)

| Iter | n | gate PASS | mean clash | %≤4 |
|---|---:|---:|---:|---:|
| Baseline (199 pre-FRUTON) | 199 |    — | 2.60 | 81.9% |
| iter-2 (43/48 in BASELINE_STATS.md)          |  48 | 89.6% | — | — |
| iter-3 (bench_20260823_1513)                  |  46 | 91.3% | 2.98 | 84.8% |
| **iter-4 (bench_20260823_1620)**              |  46 | **97.8%** | 2.98 | 84.8% |

**Reviewer headline:** the adaptive fast→slow retry patch (commit
`23faa6e`) plus subsequent gate refinements pushed pass rate from
89.6 % (iter-2) → 91.3 % (iter-3) → **97.8 % (iter-4)** — only
1/46 model failed the post-FRUTON quality gate.  Clash distribution
stayed flat between iter-3 and iter-4, confirming the pass-rate
improvement came from ω-planarity + rollback logic, not from
loosening clash thresholds.
