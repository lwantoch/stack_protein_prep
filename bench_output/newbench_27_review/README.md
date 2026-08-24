# newbench_27 — human review bundle

All 27 delivered models from the affinity_bench_27 (aka newbench_27) run,
consolidated for one-command PyMOL review over the juplaunch tunnel.

Source SLURM run:
`$LUSTRE/MMBSA_200/generic_affinity_bench_27_20260824_0918/artefacts/`

## What's here

- `crystals/<PDB>_crystal.pdb` → symlink to
  `/mnt/netapp1/Store_othcxlwa/newbench_27/<PDB>/<PDB>.pdb`
- `models/<PDB>_delivered.pdb` → the delivery for that PDB:
    - **1K4Y** — real FRUTON AF-splice output (`*_final.pdb`)
    - **26 others** — BL-Pose fallback: crystal shipped as-is (delivery
      == crystal; no AF gap-fill was possible)
- `quality_reports/<PDB>.json` — per-PDB confidence + action items
- `MANIFEST.csv` / `MANIFEST.md` — per-PDB roll-up (route, status,
  confidence, notes)
- `load_all_27.pml` — PyMOL script to open all 27 in one shot

## Bench-wide status

- **27 / 27 delivered_with_notes** (zero not_delivered)
- **All at medium confidence**
- **1 via AF route (1K4Y)** — real gap fill
- **26 via BL-Pose fallback** — crystal shipped as-is because no
  AF alignment was possible for those UniProt entries; downstream MD
  or docking can use the crystal directly.

## How to review from your laptop (via juplaunch tunnel)

Assuming `~/bin/juplaunch` is already set up (see memory
`[[reference_juplaunch]]`):

1. On laptop, start the reverse-tunnel PyMOL session:
   ```bash
   juplaunch     # opens Jupyter with -R 9123 reverse tunnel
   ```

2. Sync this directory to your laptop (once):
   ```bash
   rsync -avL cesga:/mnt/netapp2/Store_uni/home/otras/hcx/lwa/repos/FRUTON/stack_protein_prep/bench_output/newbench_27_review/ \
              ~/fruton_review_27/
   ```

3. From the `~/fruton_review_27/` directory on the laptop, launch PyMOL:
   ```bash
   cd ~/fruton_review_27/
   pymol load_all_27.pml
   ```

4. Objects appear in the PyMOL sidebar grouped by PDB id
   (`1K4Y` group → `1K4Y_crystal` + `1K4Y_delivered`). Toggle groups
   to inspect one entry at a time. Cyan cartoons = crystal reference,
   magenta = FRUTON delivery. Metals rendered as spheres.

5. For BL-Pose fallback cases, magenta and cyan overlap exactly (the
   delivery is literally the crystal). For 1K4Y, look for gap regions
   where magenta ≠ cyan — those are FRUTON's fills.

## Where to focus review

- **1K4Y** — the only entry with real gap-fill; inspect gap regions and
  active site preservation.
- **Cofactor-bearing entries** (5M7U carries XHA×2, others may carry
  HEM / NAD / SAM / ATP) — inspect `quality_reports/<PDB>.json` action
  items for cofactor frcmod review notes.
- **Any entry flagged "needs_review"** in `MANIFEST.md` — none currently
  (all 27 are `delivered_with_notes`).

## Regenerating this bundle

The bundle is built from symlinks; if the SLURM run drops a new
`_final.pdb` (e.g. AF becomes available for one of the 26 fallback
cases), rerun:

```bash
cd stack_protein_prep/
# re-runs the symlink-refresh + MANIFEST rebuild
bash scripts/refresh_newbench_27_review.sh    # (future — currently ad-hoc)
```

Or delete `bench_output/newbench_27_review/` and re-run the ad-hoc
commit-time script from git log.

## Caveats

- Symlinks point to `$LUSTRE` (CESGA scratch) — they only resolve on
  CESGA nodes. The `rsync -L` above dereferences them so laptop-side
  files are real copies.
- 26/27 are BL-Pose fallback = crystal-as-is; from a "did FRUTON
  disturb anything?" perspective, only 1K4Y is diagnostically
  meaningful in this bench. The BL-Pose route's value is that
  downstream MMBSA / MD can proceed for all 27 without a bad AF
  splice ever getting shipped.
