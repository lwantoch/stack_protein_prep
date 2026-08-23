#!/usr/bin/env bash
#SBATCH --job-name=fruton_mmbsa200
#SBATCH --partition=medium
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=1-00:00:00
# --output/--error injected at submit time (see submit block below)
#
# Full MMBSA_200 FRUTON re-optimisation bench, artefacts on LUSTRE scratch.
#
# Submit with:
#   TS=$(date +%Y%m%d_%H%M)
#   OUT=$LUSTRE/MMBSA_200/bench_${TS}
#   mkdir -p $OUT/{logs,artefacts}
#   sbatch \
#     --output=$OUT/logs/slurm.out \
#     --error=$OUT/logs/slurm.err \
#     --export=ALL,BENCH_OUT_DIR=$OUT/artefacts \
#     scripts/CESGA_SLURM/fruton_bench_mmbsa200_sbatch.sh
set -euo pipefail

FRUTON_ROOT="/mnt/netapp2/Store_uni/home/otras/hcx/lwa/repos/FRUTON/stack_protein_prep"
BENCH_DRIVER="$FRUTON_ROOT/scripts/CESGA_SLURM/fruton_bench_mmbsa200_full.py"

: "${BENCH_OUT_DIR:?BENCH_OUT_DIR must be exported (per-protein artefacts + results JSON go here)}"
echo "[fruton bench] node=$(hostname) start=$(date -Is)"
echo "[fruton bench] BENCH_OUT_DIR=${BENCH_OUT_DIR}"
echo "[fruton bench] FRUTON_ROOT=${FRUTON_ROOT}"

# --- Load CESGA modules the pipeline expects (matches fruton_cesga.sh) -----
if ! declare -f module &>/dev/null; then
    for _init in \
        /usr/share/lmod/lmod/init/bash \
        /usr/local/lmod/lmod/init/bash \
        /etc/profile.d/lmod.sh \
        /etc/profile.d/z00_lmod.sh; do
        [[ -f "$_init" ]] && source "$_init" && break
    done
fi

if declare -f module &>/dev/null; then
    module load cesga/2020 gcc/system openmpi/4.0.5_ft3_cuda \
                amber/20.13-AmberTools-22.2 mafft/7.525-with-extensions || true
fi

# --- Sanity report ---------------------------------------------------------
echo "[fruton bench] python=$(command -v python3 || echo missing)"
echo "[fruton bench] pixi=$(command -v pixi || echo missing)"
echo "[fruton bench] modeller? $(python3 -c 'import modeller; print(modeller.info.version)' 2>/dev/null || echo n/a)"
mkdir -p "${BENCH_OUT_DIR}"

# --- Run via pixi so the stack_protein_preparation venv is picked up -------
cd "$FRUTON_ROOT"
BENCH_OUT_DIR="${BENCH_OUT_DIR}" pixi run python "$BENCH_DRIVER"
RC=$?

echo "[fruton bench] driver exit=${RC}  finish=$(date -Is)"
exit $RC
