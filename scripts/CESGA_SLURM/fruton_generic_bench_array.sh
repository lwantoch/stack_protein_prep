#!/usr/bin/env bash
#SBATCH --job-name=fruton_gen
#SBATCH --partition=medium
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=0-12:00:00
# --array / --output / --error injected at submit time.  Split via env.
#
# Submit with:
#   TS=$(date +%Y%m%d_%H%M)
#   SPLIT=train                # or val / holdout
#   OUT=$LUSTRE/MMBSA_200/generic_${SPLIT}_${TS}
#   mkdir -p $OUT/{logs,artefacts}
#   sbatch --array=1-22%8 \
#     --output=$OUT/logs/task_%A_%a.out \
#     --error=$OUT/logs/task_%A_%a.err \
#     --export=ALL,BENCH_OUT_DIR=$OUT/artefacts,BENCH_SPLIT=$SPLIT \
#     scripts/CESGA_SLURM/fruton_generic_bench_array.sh
set -euo pipefail

FRUTON_ROOT="/mnt/netapp2/Store_uni/home/otras/hcx/lwa/repos/FRUTON/stack_protein_prep"
DRIVER="$FRUTON_ROOT/scripts/CESGA_SLURM/fruton_generic_bench.py"

: "${BENCH_OUT_DIR:?BENCH_OUT_DIR must be exported}"
: "${BENCH_SPLIT:?BENCH_SPLIT must be exported (train|val|holdout)}"
: "${SLURM_ARRAY_TASK_ID:?requires SLURM array}"

echo "[bench] node=$(hostname) task=${SLURM_ARRAY_TASK_ID} split=${BENCH_SPLIT} start=$(date -Is)"

if ! declare -f module &>/dev/null; then
    for _init in /usr/share/lmod/lmod/init/bash /etc/profile.d/lmod.sh; do
        [[ -f "$_init" ]] && source "$_init" && break
    done
fi
if declare -f module &>/dev/null; then
    module load cesga/2020 gcc/system openmpi/4.0.5_ft3_cuda \
                amber/20.13-AmberTools-22.2 mafft/7.525-with-extensions || true
    # Load xtb for MCPB Tier 2 (when the metal dispatcher picks xtb).
    module load cesga/2020 xtb/6.4.0 2>/dev/null || true
fi

mkdir -p "${BENCH_OUT_DIR}"

cd "$FRUTON_ROOT"
BENCH_OUT_DIR="${BENCH_OUT_DIR}" BENCH_INDEX="${SLURM_ARRAY_TASK_ID}" \
    pixi run python "$DRIVER" --split "${BENCH_SPLIT}" --index "${SLURM_ARRAY_TASK_ID}"
RC=$?
echo "[bench] task=${SLURM_ARRAY_TASK_ID} exit=${RC} finish=$(date -Is)"
exit $RC
