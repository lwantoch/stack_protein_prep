#!/usr/bin/env bash
# =============================================================================
# Gaussian HPC submission — 3BC3 | CSD RESP
# =============================================================================
# Submits two jobs from step02_gaussian/ to CESGA:
#   1. CSD_opt  HF/6-31G* geometry optimisation
#   2. CSD_esp  HF/6-31G* ESP single-point (submitted with --dependency=afterok
#               on CSD_opt, so it uses the optimised geometry via Geom=AllCheck)
#
# Run ONLY after commands.sh has completed successfully.
# After all jobs complete: bash run_after_gaussian.sh
# =============================================================================
set -euo pipefail

REPO_ROOT="$(git -C "$(dirname "${BASH_SOURCE[0]}")" rev-parse --show-toplevel)"
WORKER="${REPO_ROOT}/scripts/CESGA_SLURM/run_gaussian.sh"
LOG_DIR="${REPO_ROOT}/scripts/CESGA_SLURM/logs"
mkdir -p "$LOG_DIR"

STEP02="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/step02_gaussian"

# ---- 1. Geometry optimisation -----------------------------------------------
JOB_OPT=$(sbatch --parsable \
    --job-name="CSD_opt" \
    --partition=compute \
    --nodes=1 --ntasks=1 --cpus-per-task=16 \
    --mem=32G --time=24:00:00 \
    --output="${LOG_DIR}/CSD_opt_%j.out" \
    --error="${LOG_DIR}/CSD_opt_%j.err" \
    --requeue \
    --wrap="bash ${WORKER} ${STEP02}/CSD_opt.com")
echo "Submitted opt:  job ${JOB_OPT}"

# ---- 2. ESP single-point (depends on opt finishing) -------------------------
JOB_ESP=$(sbatch --parsable \
    --job-name="CSD_esp" \
    --dependency=afterok:${JOB_OPT} \
    --partition=compute \
    --nodes=1 --ntasks=1 --cpus-per-task=16 \
    --mem=32G --time=24:00:00 \
    --output="${LOG_DIR}/CSD_esp_%j.out" \
    --error="${LOG_DIR}/CSD_esp_%j.err" \
    --requeue \
    --wrap="bash ${WORKER} ${STEP02}/CSD_esp.com")
echo "Submitted esp:  job ${JOB_ESP}  (depends on ${JOB_OPT})"

echo ""
echo "Monitor with: squeue -u \$USER"
echo "After both complete: bash run_after_gaussian.sh"
