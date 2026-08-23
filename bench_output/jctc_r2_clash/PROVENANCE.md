# JCTC R2 — bench48 clashscore distribution

Generated 2026-08-23 by driving `scripts/plot_clashscore_histogram.py`
against the shipped `scripts/CESGA_SLURM/fruton_bench48_full_results.json`
(48-protein FRUTON benchmark, source-of-truth for the current
publication draft's baseline numbers).

## Reproduce

```bash
python scripts/plot_clashscore_histogram.py \
    --bench-json scripts/CESGA_SLURM/fruton_bench48_full_results.json \
    --outdir bench_output/jctc_r2_clash/ \
    --x-max 60
```

## Headline numbers

- n = 48 proteins
- mean clashscore = 2.92 clash pairs / structure
- p50 = 1.0
- p90 = 6.3
- p99 = 33.0
- max = 33
- 41 / 48 (85.4%) below MolProbity 90th percentile (≤ 4) per Chen et al. 2010

## Files

- `clashscore_histogram.png` — bench-wide histogram with MolProbity 90th (4) and 99th (20) reference lines.
- `clashscore_histogram.csv` — per-PDB (pdb_id, n_clash_pairs, n_vdw_clashes).
- `summary.txt` — percentile roll-up.

## Caveat

`n_vdw_clashes` is empty for this bench because the older bench48
runs pre-date the `_filler_quality_check` version that emits the
MolProbity-style vdW-radius clash count. Only `n_clash_pairs` (all-atom
< 2.0 Å blanket) is available here. A rerun on the full MMBSA_200
bench with the current pipeline will populate both.
