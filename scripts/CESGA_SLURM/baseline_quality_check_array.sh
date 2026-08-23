#!/usr/bin/env bash
#SBATCH --job-name=fruton_qc
#SBATCH --partition=medium
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=0-01:00:00
# Baseline crystal quality-check on all 208 FRUTON-NEW crystals.
# One array task = one crystal.  Fast (no MODELLER); expect < 5 min each.
#
# Submit with:
#   TS=$(date +%Y%m%d_%H%M)
#   OUT=$LUSTRE/MMBSA_200/baseline_qc_${TS}
#   mkdir -p $OUT/{logs,artefacts}
#   sbatch --array=1-208%16 \
#     --output=$OUT/logs/task_%A_%a.out --error=$OUT/logs/task_%A_%a.err \
#     --export=ALL,BENCH_OUT_DIR=$OUT/artefacts \
#     scripts/CESGA_SLURM/baseline_quality_check_array.sh
set -euo pipefail

FRUTON_ROOT="/mnt/netapp2/Store_uni/home/otras/hcx/lwa/repos/FRUTON/stack_protein_prep"
DRIVER="$FRUTON_ROOT/scripts/CESGA_SLURM/baseline_quality_check_full.py"

: "${BENCH_OUT_DIR:?BENCH_OUT_DIR must be exported}"
: "${SLURM_ARRAY_TASK_ID:?requires SLURM array}"

echo "[qc] node=$(hostname) task=${SLURM_ARRAY_TASK_ID} start=$(date -Is)"

if ! declare -f module &>/dev/null; then
    for _init in /usr/share/lmod/lmod/init/bash /etc/profile.d/lmod.sh; do
        [[ -f "$_init" ]] && source "$_init" && break
    done
fi
declare -f module &>/dev/null && module load cesga/2020 gcc/system || true

mkdir -p "${BENCH_OUT_DIR}"
cd "$FRUTON_ROOT"
BENCH_OUT_DIR="${BENCH_OUT_DIR}" BENCH_INDEX="${SLURM_ARRAY_TASK_ID}" \
    pixi run python "$DRIVER" --index "${SLURM_ARRAY_TASK_ID}"
RC=$?
echo "[qc] task=${SLURM_ARRAY_TASK_ID} exit=${RC} finish=$(date -Is)"
exit $RC
