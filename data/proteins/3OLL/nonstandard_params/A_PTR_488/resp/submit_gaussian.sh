#!/usr/bin/env bash
# =============================================================================
# Gaussian HPC submission — PTR
# =============================================================================
# Delegates to the central CESGA Gaussian submission script.
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
exec bash "$CENTRAL" -r "$SCRIPT_DIR/step02_gaussian" --partition short --time 06:00:00 --mem 4000M --cpus 32 --gpus 1 "$@"
