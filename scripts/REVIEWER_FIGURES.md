# FRUTON reviewer figures — script index

Each script here answers one **reviewer concern** from the target
journals (Nature Methods, JCTC, JACS).  All are pure-Python +
matplotlib, license-free, and never touch the network or require GPU
tools.  They consume artefacts the pipeline already emits (bench
results JSON, quality-report JSON, per-PDB PDBs).

Run any script with `--help` for the full argument list.

---

## Reviewer figure map

| # | Script | Reviewer concern | Deliverable |
|---|---|---|---|
| **JACS R1** | `plot_omega_planarity_distribution.py` | Does FRUTON preserve peptide-bond planarity? | ω dihedral histogram bench-wide (trans / cis-Pro / cis-nonPro / non-planar breakdown, ref lines ±30° / ±150°) |
| **JACS R2** | `demo_pdbfixer_vs_fruton_protonation.py` + `data/pdbfixer_removal_justification.md` | Chemistry-defensible reason to reject pdbfixer | 5-mode failure catalogue + side-by-side verdict CLI (naive pH=7 rule vs FRUTON paper-override extractor, no pdbfixer install needed) |
| **JACS R3** | `plot_metal_donor_preservation.py` | Are metal-coord side chains displaced? | Bench-wide Δd histogram + per-donor-type mean \|Δd\|; ref line at ±0.3 Å (Harding 2001) |
| **JACS R4** | *(module `_dssp_ss3.py`, no plot)* | Does MobiDB IDR flag ignore local crystal structure? | Bio.PDB.DSSP wrapper + `cross_check_idr_vs_local_ss` advisory ("OVERRIDE IDR reject" when flanks are H / E) |
| **JCTC R2** | `plot_clashscore_histogram.py` | Bench-wide clashscore distribution | Two-panel histogram (n_clash_pairs + n_vdw_clashes) with MolProbity 90th/99th ref lines (Chen 2010) |
| **JCTC R3** | *(gate `_filler_quality_check._chirality_check`, no plot)* | Do splices flip Cα chirality to D? | Signed-volume + improper-dihedral guard (\|χ_N-CA-C-CB\| ≈ 120°) |
| **JCTC R4** | *(module `_uniprot_idr.py`)* | MobiDB API dependency = reproducibility risk | Local `mobidb_snapshot.json` cache + `use_cache=True` default |
| **Nature R1** | `plot_family_stratification.py` | Aggregate hides class imbalance | Per-family pass rate + rescued residues bar (kinase / GPCR / protease / metalloenzyme / …); overall reference line |
| **Nature R2** | `plot_pdb_redo_comparison.py` | Community re-refinement benchmark | Per-metric stacked bar (FRUTON better / PDB-REDO better / tie / unavailable) across 9 canonical metrics |
| **Nature R3** | `plot_blind_test_rmsd.py` | Self-consistency ≠ ground truth | Hide-crystal driver (`_blind_test.mask_crystal_pdb`) + Cα-RMSD histogram vs held-out crystal + per-PDB sorted RMSD bar |
| **Nature R4** | `plot_ablation_comparison.py` | Which gate contributes what? | JSON-driven per-variant bench comparison (pass rate + rescued residues, 9 named toggle-able gates) |
| USER R1 | *(helper `_metall_params_helpers._propose_water_for_incomplete_coordination`)* | Metal geometry needs water | Auto-suggest water placement for under-coordinated metal ions |
| USER R2 | *(module `_paper_override_suggest.py`)* | Active-site protonation from literature | Regex-based extractor writing `.suggested.json` (user reviews before promotion to real override) |

---

## Example CLIs

Each script accepts `-h` for the full arg list; here are typical
invocations against the MMBSA_200 bench.

### JACS R1 — ω peptide-bond planarity

```bash
python scripts/plot_omega_planarity_distribution.py \
    --pdb-glob 'bench/*/final_model.pdb' \
    --outdir figures/jacs_r1_omega/ \
    --log-y
```

### JACS R2 — pdbfixer-vs-FRUTON protonation demo

```bash
python scripts/demo_pdbfixer_vs_fruton_protonation.py
# Or with a real per-protein evidence file:
python scripts/demo_pdbfixer_vs_fruton_protonation.py \
    --evidence bench/8CA2/paper_evidence.md
```

### JACS R3 — Metal-donor preservation

```bash
python scripts/plot_metal_donor_preservation.py \
    --crystal-glob 'bench/*/input_crystal.pdb' \
    --fruton-glob  'bench/*/final_model.pdb' \
    --pair-by parent \
    --outdir figures/jacs_r3_metal/
```

### JCTC R2 — Bench-wide clashscore histogram

```bash
python scripts/plot_clashscore_histogram.py \
    --bench-json runs/mmbsa200_baseline.json \
    --outdir figures/jctc_r2_clash/
# Or from per-PDB quality reports:
python scripts/plot_clashscore_histogram.py \
    --pdb-metrics-dir metrics/mmbsa200_fruton/ \
    --outdir figures/jctc_r2_clash/ \
    --log-y
```

### Nature R1 — Family-stratified analysis

```bash
python scripts/plot_family_stratification.py \
    --bench-json runs/mmbsa200_baseline.json \
    --family-map src/stack_protein_preparation/data/family_by_pdb_seed.json \
    --outdir figures/nature_r1_family/
```

### Nature R2 — PDB-REDO comparison

```bash
python scripts/plot_pdb_redo_comparison.py \
    --fruton-dir   metrics/mmbsa200_fruton/ \
    --pdb-redo-dir metrics/mmbsa200_pdb_redo/ \
    --outdir figures/nature_r2/
```

### Nature R3 — Blind-test crystal-held-out

```bash
python scripts/plot_blind_test_rmsd.py \
    --bench-spec metrics/blind_test_bench.json \
    --outdir figures/nature_r3/ \
    --log-y
```

The bench-spec is a list of
`{pdb_id, crystal_pdb, filled_pdb, held_out: [{chain, first_resnum, last_resnum}, …]}`.
Use `_blind_test.mask_crystal_pdb` to prepare each masked input
before running FRUTON.

### Nature R4 — Gate ablation

```bash
python scripts/plot_ablation_comparison.py \
    --baseline runs/mmbsa200_baseline.json \
    --variant no_plddt=runs/mmbsa200_no_plddt.json:plddt \
    --variant no_idr=runs/mmbsa200_no_idr.json:idr \
    --variant no_omega=runs/mmbsa200_no_omega.json:omega \
    --outdir figures/nature_r4_ablation/
```

Each `--variant` is `label=json_path[:disabled_gate]`; the optional
`:disabled_gate` label must be one of the keys in
`_ablation.GATE_NAMES`.

---

## Utility scripts (not reviewer-facing)

| Script | Purpose |
|---|---|
| `fruton.py` | Entry point — runs the full FRUTON pipeline on one protein |
| `stage_for_gbsa.py` | Prepares a FRUTON output for downstream GBSA scoring |
| `summarize_bench_results.py` | CSV + markdown roll-up of a bench JSON |
| `fetch_mobidb_snapshot.py` | Populates `data/mobidb_snapshot.json` for the IDR cache |

---

## Design invariants

- **License-free**: pure Python + Bio.PDB + matplotlib + AmberTools
  (all BSD-3 / MIT / GPL-compatible).  Never requires phenix,
  Schrödinger, pdbfixer, RFdiffusion2, or Boltz-2 to run.
- **Fail-open**: any external subprocess (DSSP, sander, antechamber,
  MODELLER, RFdiffusion2) that is missing degrades gracefully; the
  script returns a `ran=False + fallback_reason` result instead of
  raising.
- **Reviewer chain-of-evidence**: each figure ships alongside a CSV
  or JSON sidecar so the reviewer can drill from the aggregate figure
  to the per-PDB / per-residue underlying data.
