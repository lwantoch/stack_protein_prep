# MMBSA_200 Baseline Reviewer Statistics — 2026-08-23

Pre-FRUTON crystal-only quality-check across all 199 crystals available
in `/mnt/netapp1/Store_othcxlwa/FRUTON-NEW/`.  No modelling was applied
— these numbers describe the raw deposited coordinates.

## Aggregate statistics

| Metric             | Mean  | Median | Max   |
|--------------------|------:|-------:|------:|
| n_residues         |  474  |   342  |  1858 |
| clash pairs        |  2.60 |     1  |    33 |
| ω non-planar       |  0.18 |     0  |     3 |
| ω cis-nonPro       |  0.22 |     0  |     5 |
| Rama outliers      |  8.3  |     6  |    42 |
| broken bonds       |  0.38 |     0  |    10 |

## Reviewer-quality fraction

| Property                    |  n / 199 |    %  |
|-----------------------------|---------:|------:|
| 0 broken peptide bonds      |   169    | 84.9  |
| 0 non-planar ω              |   172    | 86.4  |
| 0 cis-nonPro ω              |   171    | 85.9  |

## ω distribution across all deposited pairs

Scanning every consecutive-residue peptide bond across 199 crystals:

- **187 252 residue pairs** scanned
- **87 cis-nonPro** (0.046 %) — MacArthur & Thornton 1991 report ~0.03–0.1 %
- **73 non-planar** (0.039 %) — expected ~0 for well-refined structures

Both figures are well below the 1 % reviewer red-flag threshold, which
means the raw crystal set is a strong starting point.  FRUTON must not
push these numbers up materially — the ω-planarity gate is what enforces
that on the filled models.

## Data provenance

- Source PDBs: `/mnt/netapp1/Store_othcxlwa/FRUTON-NEW/*/<PDB>.pdb`
- Baseline QC JSON: `$LUSTRE/MMBSA_200/baseline_qc_20260823_1513/artefacts/baseline_quality_check_results.json`
- ω figure: `omega_distribution.png` (log-y bench-wide histogram)
- Per-pair CSV: `$LUSTRE/MMBSA_200/figures_20260823/omega_distribution.csv`
- Per-PDB summary: `per_pdb_summary.txt`

Reproduce with:
```
pixi run python scripts/CESGA_SLURM/baseline_quality_check_full.py \\
    # writes $BENCH_OUT_DIR/baseline_quality_check_results.json
```

## Next iteration signals for reviewer defensibility

The FRUTON pipeline iterations of 2026-08-23 (iter-1, iter-2) both
scored 43/48 PASS = 89.6 % on the AF-available subset.  The 5 FAILs
were all ω-planarity gate rejects; iter-3 (bench_20260823_1513) tests
the adaptive fast→slow retry patch (commit 23faa6e) that should
rescue those.
