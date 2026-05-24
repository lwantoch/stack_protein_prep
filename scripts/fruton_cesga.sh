#!/usr/bin/env bash
# =============================================================================
# FRUTON CESGA launcher
#
# Loads required CESGA modules then runs FRUTON via pixi.
# $HOME/bin/fruton should point here:
#
#   ln -sf "$STORE/stack_protein_prep/scripts/fruton_cesga.sh" "$HOME/bin/fruton"
#
# Module versions are pinned here — update and git pull to change them.
# =============================================================================
set -euo pipefail

FRUTON_ROOT="$(cd "$(dirname "$(readlink -f "${BASH_SOURCE[0]}")")/.." && pwd)"

# ---- Initialize Lmod if module function not yet available ------------------
# Needed when the script is called non-interactively (sbatch, cron, etc.)
if ! declare -f module &>/dev/null; then
    for _init in \
        /usr/share/lmod/lmod/init/bash \
        /usr/local/lmod/lmod/init/bash \
        /etc/profile.d/lmod.sh \
        /etc/profile.d/z00_lmod.sh; do
        [[ -f "$_init" ]] && source "$_init" && break
    done
fi

# ---- Load CESGA modules ----------------------------------------------------
# Load Amber and MAFFT (both cesga/2020) first, then swap to cesga/2025 for
# GROMACS. After the swap Amber and MAFFT become "Inactive" in Lmod but their
# binaries remain in PATH and are fully functional.
if declare -f module &>/dev/null; then
    module load cesga/2020 gcc/system openmpi/4.0.5_ft3_cuda amber/20.13-AmberTools-22.2
    module load mafft/7.525-with-extensions
    module load cesga/2025 gcc/system openmpi/4.1.8 gromacs/2025.3
    # Gaussian is loaded per-job by run_gaussian.sh (licensed, node-only)
else
    echo "WARNING: module system not found — GROMACS/AmberTools/MAFFT may be missing from PATH" >&2
fi

# ---- Run FRUTON via pixi ---------------------------------------------------
exec pixi run --manifest-path "$FRUTON_ROOT/pixi.toml" \
    python "$FRUTON_ROOT/scripts/fruton.py" "$@"
