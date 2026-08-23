#!/usr/bin/env bash
#SBATCH --job-name=fruton_mm200a
#SBATCH --partition=medium
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=0-12:00:00
# --array / --output / --error injected at submit time.
#
# Parallel per-protein SLURM array driver for the MMBSA_200 bench.
# One array task = one PDB.  Artefacts + <PDB>.json on LUSTRE scratch.
#
# Submit with:
#   TS=$(date +%Y%m%d_%H%M)
#   OUT=$LUSTRE/MMBSA_200/bench_${TS}
#   mkdir -p $OUT/{logs,artefacts}
#   # 48 AF-available PDBs today; %8 caps concurrent tasks.
#   sbatch \
#     --array=1-48%8 \
#     --output=$OUT/logs/task_%A_%a.out \
#     --error=$OUT/logs/task_%A_%a.err \
#     --export=ALL,BENCH_OUT_DIR=$OUT/artefacts \
#     scripts/CESGA_SLURM/fruton_bench_mmbsa200_array.sh
set -euo pipefail

FRUTON_ROOT="/mnt/netapp2/Store_uni/home/otras/hcx/lwa/repos/FRUTON/stack_protein_prep"
BENCH_DRIVER="$FRUTON_ROOT/scripts/CESGA_SLURM/fruton_bench_mmbsa200_full.py"

: "${BENCH_OUT_DIR:?BENCH_OUT_DIR must be exported}"
: "${SLURM_ARRAY_TASK_ID:?this script requires SLURM array}"

echo "[fruton bench] node=$(hostname) task=${SLURM_ARRAY_TASK_ID} start=$(date -Is)"
echo "[fruton bench] BENCH_OUT_DIR=${BENCH_OUT_DIR}"

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

mkdir -p "${BENCH_OUT_DIR}"

cd "$FRUTON_ROOT"
BENCH_OUT_DIR="${BENCH_OUT_DIR}" BENCH_INDEX="${SLURM_ARRAY_TASK_ID}" \
    pixi run python "$BENCH_DRIVER" --index "${SLURM_ARRAY_TASK_ID}"
RC=$?

echo "[fruton bench] task=${SLURM_ARRAY_TASK_ID} exit=${RC} finish=$(date -Is)"
exit $RC
