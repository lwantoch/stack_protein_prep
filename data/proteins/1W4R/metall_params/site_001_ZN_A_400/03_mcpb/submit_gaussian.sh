#!/usr/bin/env bash
# =============================================================================
# Gaussian HPC submission — 1W4R | ZN MCPB
# =============================================================================
# Submits three Gaussian jobs from step02_gaussian/ to CESGA as a job array:
#   1W4R_ZN_A_400_small_opt.com   geometry optimisation (small model)
#   1W4R_ZN_A_400_small_fc.com    force constants / Hessian (small model)
#   1W4R_ZN_A_400_large_mk.com    RESP charges (large model)
#
# All three jobs run in parallel (no inter-job dependencies).
# Run ONLY after commands.sh has completed successfully.
#
# Usage:
#   bash submit_gaussian.sh              # submit all missing jobs
#   bash submit_gaussian.sh --only-failed  # resubmit only failed jobs
#
# After all jobs complete: bash run_after_gaussian.sh
# (run_after_gaussian.sh will run formchk on small_opt.chk automatically)
# =============================================================================
set -euo pipefail

REPO_ROOT="$(git -C "$(dirname "${BASH_SOURCE[0]}")" rev-parse --show-toplevel)"
CENTRAL="${REPO_ROOT}/scripts/CESGA_SLURM/submit_gaussian.sh"

exec bash "$CENTRAL" \
    -r "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/step02_gaussian" \
    --partition compute \
    --time 24:00:00 \
    --mem 32G \
    --cpus 16 \
    "$@"
