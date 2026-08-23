# JACS R3 — iter-3 metal-donor preservation

Generated 2026-08-23 by driving `scripts/plot_metal_donor_preservation.py`
against paired crystal ↔ FRUTON iter-3 filled models.  Uses
`_metal_coord_scan.compare_metal_coordination` to compute per-donor
Δd = d(FRUTON) − d(crystal) for every metal-donor pair.

## Reproduce

Set up matched stem-name symlinks (crystals live under
`FRUTON-NEW/<PDB>/<PDB>.pdb`; filled PDBs under
`bench_20260823_1513/artefacts/<PDB>_final.pdb`), then:

```bash
python scripts/plot_metal_donor_preservation.py \
    --crystal-glob '<crystal_stems>/*.pdb' \
    --fruton-glob '<filled_stems>/*.pdb' \
    --pair-by stem \
    --outdir bench_output/jacs_r3_metal_iter3/ \
    --title 'Iter-3 metal-donor preservation (Δd vs crystal, 46 filled)'
```

## Headline numbers

- **46 iter-3 filled models paired against crystals**
- **11 crystals carry metal ions** (25 %); 35 metal-free.
- **58 metal-donor contacts scored** across the 11 metallo-enzymes.
- **mean Δd = +0.000 Å**
- **max |Δd| = 0.000 Å**
- **58 / 58 (100 %) preserved** — no `lost` (crystal-only) or
  `gained` (FRUTON-only) contacts.

## Reviewer headline

FRUTON **preserved every single metal-donor contact bit-for-bit**
across all 11 metallo-enzyme benchmarks in iter-3.  Every His-Nε,
Asp-Oδ, Glu-Oε, Cys-Sγ, Thr-Oγ, and Asn-Oδ donor sits at exactly
the same distance to its metal ion in the filled model as in the
crystal.

Why: FRUTON's fill logic only introduces new residues at gap
positions and only refines new + immediately-flanking residues via
MODELLER LoopModel.  Metal ions and their coordinating side chains
almost always sit far from the gaps (the crystal resolved them, so
they're not in a REMARK 465 stretch).  The 100 % Δd = 0 result is
the direct empirical confirmation.

## Per-PDB scores (metal-carrying subset)

| PDB  | metals | donors | mean Δd | status |
|------|-------:|-------:|--------:|--------|
| 2YCF | 1 Mg   |     2  |   0.000 | preserved |
| 3S95 | 2 Na   |     2  |   0.000 | preserved |
| 3ZCW | 1 Mg   |     1  |   0.000 | preserved |
| 6VPM | 2 Mg   |     4  |   0.000 | preserved |
| 6YA6 | 2 Zn   |     8  |   0.000 | preserved |
| 7B7R | 2 Mg   |     4  |   0.000 | preserved |
| 7SXF | 1 —    |     3  |   0.000 | preserved |
| 7UL2 | 1 —    |     2  |   0.000 | preserved |
| 8C2Z | 3 —    |     5  |   0.000 | preserved |
| 8DSC | 1 —    |     0  |     n/a | metal only, no donor in cutoff |
| 9E3M | 7 —    |    27  |   0.000 | preserved |

## Files

- `metal_donor_delta.png` — Δd histogram (all zeros here) + per-donor-type mean |Δd| bar.
- `metal_donor_delta.csv` — one row per matched donor contact (58 rows), cols: pdb_id, metal_element, metal_chain, metal_resnum, donor_chain, donor_resnum, donor_resname, donor_atom, d_crystal_A, d_fruton_A, delta_A, status.
- `per_pdb_summary.txt` — per-PDB one-liner.

## Cross-check for reviewer

If a reviewer challenges "FRUTON might disturb the active site":
open `metal_donor_delta.csv`, sort by `abs(delta_A)`, and observe
that every non-empty row has `delta_A = +0.000`.  The empirical Δd
tolerance of Harding 2001 (metal-donor ± 0.3 Å) is trivially
satisfied.
