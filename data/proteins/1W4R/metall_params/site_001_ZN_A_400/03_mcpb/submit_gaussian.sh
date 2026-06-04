#!/usr/bin/env bash
# =============================================================================
# Gaussian HPC submission — 1W4R_ZN_A_400
# =============================================================================
# Submits the three MCPB Gaussian jobs with the correct SLURM dependency:
#   small_opt  ──┐ (parallel)
#   large_mk     │
#   small_fc  ←──┘ afterok:small_opt  (needs small_opt.chk)
#
# Run ONLY after bash commands.sh has completed successfully.
# After all jobs complete:  bash run_after_gaussian.sh
# =============================================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
if REPO_ROOT="$(git -C "$SCRIPT_DIR" rev-parse --show-toplevel 2>/dev/null)"; then
    CENTRAL="${REPO_ROOT}/scripts/CESGA_SLURM/submit_gaussian.sh"
else
    CENTRAL=""
    DIR="$SCRIPT_DIR"
    for _ in 1 2 3 4 5 6 7 8; do
        if [[ -f "${DIR}/scripts/CESGA_SLURM/submit_gaussian.sh" ]]; then
            CENTRAL="${DIR}/scripts/CESGA_SLURM/submit_gaussian.sh"
            break
        fi
        DIR="$(dirname "$DIR")"
    done
fi
[[ -n "${CENTRAL:-}" && -f "$CENTRAL" ]] || { echo "ERROR: central submit_gaussian.sh not found" >&2; exit 2; }

GAUSS_DIR="$SCRIPT_DIR/step02_gaussian"
SLURM_OPTS=(--partition short --time 06:00:00 --mem 32G --cpus 16)

# 1. small_opt — geometry optimisation (produces .chk needed by small_fc)
OPT_OUT=$(bash "$CENTRAL" "$GAUSS_DIR/1W4R_ZN_A_400_small_opt.com" "${SLURM_OPTS[@]}" "$@")
echo "$OPT_OUT"
OPT_JOB=$(echo "$OPT_OUT" | grep -oP '(?<=Submitted batch job )\d+' || true)
if [[ -z "$OPT_JOB" ]]; then
    echo "WARNING: could not parse small_opt job ID — small_fc will be submitted without dependency" >&2
fi

# 2. large_mk — ESP calculation, independent of small_opt
bash "$CENTRAL" "$GAUSS_DIR/1W4R_ZN_A_400_large_mk.com" "${SLURM_OPTS[@]}" "$@"

# 3. small_fc — force constants, depends on small_opt.chk
if [[ -n "$OPT_JOB" ]]; then
    bash "$CENTRAL" "$GAUSS_DIR/1W4R_ZN_A_400_small_fc.com" "${SLURM_OPTS[@]}"         --dependency "afterok:$OPT_JOB" "$@"
else
    bash "$CENTRAL" "$GAUSS_DIR/1W4R_ZN_A_400_small_fc.com" "${SLURM_OPTS[@]}" "$@"
fi

echo ""
echo "Submitted: small_opt (job $OPT_JOB), large_mk, small_fc (afterok:$OPT_JOB)"
echo "Once all complete:  bash run_after_gaussian.sh"
