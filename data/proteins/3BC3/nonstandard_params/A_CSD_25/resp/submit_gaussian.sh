#!/usr/bin/env bash
# =============================================================================
# Gaussian HPC submission — CSD
# =============================================================================
# Submit geometry optimisation and ESP jobs from step02_gaussian/.
# Run ONLY after bash commands.sh has completed successfully.
# Edit --ntasks, --mem, --time to match your cluster.
#
# Note: submit ESP job AFTER optimisation completes (it reads the chk file).
# After all jobs complete:  bash run_after_gaussian.sh
# =============================================================================
set -euo pipefail
cd step02_gaussian

# --- 1. Geometry optimisation (HF/6-31G*) -----------------------------------
JOB_OPT=$(sbatch --parsable \
    --job-name="CSD_opt" --ntasks=8 --mem=16G --time=08:00:00 \
    --wrap="g16 < CSD_opt.com > CSD_opt.log")
echo "Submitted opt job: ${JOB_OPT}"

# --- 2. ESP single-point (depends on opt via chk) ---------------------------
# The ESP .com reads the optimised geometry from the chk file.
# It is safe to submit immediately — Slurm will hold it until opt finishes.
sbatch --dependency=afterok:${JOB_OPT} \
    --job-name="CSD_esp" --ntasks=8 --mem=16G --time=04:00:00 \
    --wrap="g16 < CSD_esp.com > CSD_esp.log"

echo "Submitted ESP job (depends on opt)."
echo "Monitor with:  squeue -u \$USER"
echo "After both complete:  bash run_after_gaussian.sh"
