# Paper intro — FRUTON iteration-over-iteration progression

Killer publication-intro figure showing FRUTON's hardening trajectory
across the 2026-08-23 iteration loop.  At a glance a reviewer sees:

- Gate PASS rate climbing **89.6 % → 91.3 % → 97.8 % → 97.9 %**
- Non-planar ω count **halving** from 6 → 3 iter-3 → iter-4.

## Reproduce

```bash
python scripts/plot_iteration_progression.py \
    --history src/stack_protein_preparation/data/iteration_history_20260823.json \
    --outdir bench_output/paper_intro_progression/ \
    --title 'FRUTON MMBSA_200 hardening: 89.6% → 97.9% gate PASS across 4 iterations'
```

## Headline numbers

- **6 iterations plotted** (iter-2 → final 4-split)
- **Δ gate PASS iter-2 → iter-5-mini = +8.3 pp** (89.6 → 97.9)
- **Final train+val combined = 95.8 %** (46/48) per commit `b73576d`
- **Final 4-split (train + test + affinity_bench_27 + stresstest_30) =
  97.1 %** (102/105) per commit `09facc5`:
    - train: 22/22 = 100.0 % (all AF-ready)
    - test: 24/26 = 92.3 % (2 needs_review = 7QUE + 8Q68)
    - affinity_bench_27: 27/27 = 100.0 % (1 AF + 26 BL-Pose fallback)
    - stresstest_30: 29/30 = 96.7 % (5 AF + 25 BL-Pose)
  Zero not_delivered; BL-Pose fallback confirmed on 51 non-AF proteins.
- **Δ non-planar ω = −50 %** (6 → 3, iter-3 to iter-4 measured)
- **3 residual needs_review** at grand total: 7QUE + 8Q68 + 1 from
  stresstest_30 (all catastrophic-splice above p99 of 199-crystal
  reference distribution)

## Iteration source data

Ships as `src/stack_protein_preparation/data/iteration_history_20260823.json`.
Each entry:
- `iter_label` — human name (iter-2, iter-3, …)
- `n_pdbs` — bench-size for that iteration
- `gate_pass_pct` — from the aggregated bench-results JSON
- `n_omega_non_planar` — from bench_output/jacs_r1_omega_iter*/ scans;
  null when the ω-figure did not run for that iteration
- `notes` — provenance pointer (SLURM run dir / commit / BASELINE_STATS.md)

## Files

- `iteration_progression.png` — 2-panel matplotlib chart (gate PASS %
  left, non-planar ω count right; ω panel uses symlog-y so the
  6→3 halving is visible alongside the 0 baseline).
- `iteration_progression.csv` — machine-readable per-iter roll-up.

## Cross-check

Numbers in this figure are grounded in the individual per-iter
figures already under `bench_output/`:
- `jctc_r2_clash_iter3_post_fruton/` (iter-3 PASS %)
- `jctc_r2_clash_iter4/` (iter-4 PASS %)
- `jacs_r1_omega_iter3/` (iter-3 non-planar ω)
- `jacs_r1_omega_iter4/` (iter-4 non-planar ω)
- `docs/reviewer_evidence_20260823/BASELINE_STATS.md` (iter-2 PASS %)
- `431f870` commit (iter-5-mini PASS %)
