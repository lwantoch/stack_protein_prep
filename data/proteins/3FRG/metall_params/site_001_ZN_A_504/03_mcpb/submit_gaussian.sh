#!/usr/bin/env bash
# =============================================================================
# Gaussian HPC submission — 3FRG_ZN_A_504
# =============================================================================
# Submits three independent Gaussian jobs from step02_gaussian/.
# Run ONLY after bash commands.sh has completed successfully.
# Edit --ntasks, --mem, --time to match your cluster.
# After all jobs complete:  bash run_after_gaussian.sh
# =============================================================================
set -euo pipefail
cd step02_gaussian
GRP="3FRG_ZN_A_504"

# --- 1. Geometry optimisation (small model) ---------------------------------
sbatch --job-name="${GRP}_opt" --ntasks=8 --mem=16G --time=04:00:00 \
    --wrap="g16 < ${GRP}_small_opt.com > ${GRP}_small_opt.log \
         && formchk ${GRP}_small_opt.chk ${GRP}_small_opt.fchk"

# --- 2. Force-constant calculation (small model) ----------------------------
sbatch --job-name="${GRP}_fc" --ntasks=8 --mem=16G --time=04:00:00 \
    --wrap="g16 < ${GRP}_small_fc.com > ${GRP}_small_fc.log"

# --- 3. RESP charges (large model) ------------------------------------------
sbatch --job-name="${GRP}_mk" --ntasks=8 --mem=32G --time=08:00:00 \
    --wrap="g16 < ${GRP}_large_mk.com > ${GRP}_large_mk.log"

echo "Submitted 3 Gaussian jobs.  Monitor with:  squeue -u \$USER"
echo "After all complete (from 03_mcpb/):  bash run_after_gaussian.sh"
