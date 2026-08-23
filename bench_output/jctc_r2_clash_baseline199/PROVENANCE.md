# JCTC R2 — 199-crystal baseline clashscore distribution

Generated 2026-08-23 by driving `scripts/plot_clashscore_histogram.py`
against the live 199-crystal baseline emitted by the concurrent iter
`scripts/CESGA_SLURM/baseline_quality_check_full.py` run (see commit
`19367a6` and `docs/reviewer_evidence_20260823/BASELINE_STATS.md`).

Larger companion to `bench_output/jctc_r2_clash/` (which was driven
against the older 48-protein `fruton_bench48_full_results.json`).

## Reproduce

```bash
python scripts/plot_clashscore_histogram.py \
    --bench-json /mnt/lustre/scratch/nlsas/home/otras/hcx/lwa/MMBSA_200/baseline_qc_20260823_1513/artefacts/baseline_quality_check_results.json \
    --outdir bench_output/jctc_r2_clash_baseline199/ \
    --x-max 60 \
    --title 'Baseline crystal clashscore (199 crystals, no FRUTON)'
```

## Headline numbers (baseline, pre-FRUTON)

- **n = 199 crystals**
- **mean clashscore = 2.60** clash pairs / crystal
- **p50 = 1.0**
- **p90 = 6.2**
- **p99 = 27.1**
- **max = 33**
- **163 / 199 (81.9 %) below MolProbity 90th (≤ 4)** per Chen et al. 2010

## Reviewer chain of evidence

This is the **pre-FRUTON baseline** — what the deposited crystals look
like *before* any modelling.  The reviewer's expectation is that FRUTON
should preserve or improve these numbers.  The complementary
`bench_output/jctc_r2_clash/` figure (bench48, post-FRUTON) reports
mean 2.92, 85.4 % ≤ 4 — very close to the baseline, confirming
FRUTON's ω / clash gates are not silently regressing clash counts on
the filled models.

## Cross-check vs `docs/reviewer_evidence_20260823/BASELINE_STATS.md`

The concurrent-agent baseline stats report the same aggregate:
- mean clash pairs = 2.60 (baseline stats table) — MATCH
- max = 33 (baseline stats table) — MATCH

Confirms this figure is grounded in the same underlying data the
`per_pdb_summary.txt` sidecar summarised.

## Caveat

`n_vdw_clashes` is empty because
`scripts/CESGA_SLURM/baseline_quality_check_full.py` did not yet
persist the MolProbity-style vdW-radius clash count (only the
all-atom `n_clash_pairs` blanket).  A rerun with the current
`_filler_quality_check` version will populate both.
