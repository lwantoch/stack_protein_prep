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
# Amber and MAFFT require cesga/2020; GROMACS requires cesga/2025.  The swap
# to cesga/2025 marks amber and mafft "Inactive" in Lmod AND removes their
# PATH/LD_LIBRARY_PATH entries.  We snapshot the paths before the swap and
# restore them manually afterwards so all tools remain accessible.
if declare -f module &>/dev/null; then
    module load cesga/2020 gcc/system openmpi/4.0.5_ft3_cuda amber/20.13-AmberTools-22.2
    module load mafft/7.525-with-extensions

    # Snapshot env vars that the cesga/2025 swap will remove
    _AMBERHOME="${AMBERHOME:-}"
    _AMBER_PYTHONPATH="${PYTHONPATH:-}"      # amber module adds pymsmt etc. to PYTHONPATH
    _MAFFT_BIN="$(command -v mafft 2>/dev/null || true)"
    _MAFFT_BINARIES="${MAFFT_BINARIES:-}"   # companion binaries dir used by mafft wrapper

    # Swap to cesga/2025 for GROMACS (deactivates amber + mafft modules)
    module load cesga/2025 gcc/system openmpi/4.1.8 gromacs/2025.3

    # Re-inject Amber paths stripped by the cesga swap
    if [[ -n "$_AMBERHOME" ]]; then
        export AMBERHOME="$_AMBERHOME"
        export PATH="$_AMBERHOME/bin:$PATH"
        [[ -d "$_AMBERHOME/lib" ]] && \
            export LD_LIBRARY_PATH="$_AMBERHOME/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
    fi
    # Re-inject PYTHONPATH for AmberTools Python packages (pymsmt, etc.)
    [[ -n "$_AMBER_PYTHONPATH" ]] && \
        export PYTHONPATH="${_AMBER_PYTHONPATH}${PYTHONPATH:+:$PYTHONPATH}"
    # Re-inject MAFFT path and companion-binaries var if not already present
    if [[ -n "$_MAFFT_BIN" ]]; then
        _MAFFT_DIR="$(dirname "$_MAFFT_BIN")"
        [[ ":$PATH:" != *":$_MAFFT_DIR:"* ]] && export PATH="$_MAFFT_DIR:$PATH"
    fi
    [[ -n "$_MAFFT_BINARIES" ]] && export MAFFT_BINARIES="$_MAFFT_BINARIES"
    # Gaussian is loaded per-job by run_gaussian.sh (licensed, node-only)
else
    echo "WARNING: module system not found — GROMACS/AmberTools/MAFFT may be missing from PATH" >&2
fi

# ---- Run FRUTON via pixi ---------------------------------------------------
exec pixi run --manifest-path "$FRUTON_ROOT/pixi.toml" \
    python "$FRUTON_ROOT/scripts/fruton.py" "$@"
