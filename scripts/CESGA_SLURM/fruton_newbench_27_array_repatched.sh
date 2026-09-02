#!/bin/bash
#SBATCH -J fruton_nb27_rp
#SBATCH -p medium
#SBATCH -N 1
#SBATCH -c 4
#SBATCH --mem 16G
#SBATCH -t 8:00:00
#SBATCH -o /mnt/netapp1/Store_othcxlwa/newbench_27/logs/array_repatched/nb27rp-%A_%a.log
#SBATCH -e /mnt/netapp1/Store_othcxlwa/newbench_27/logs/array_repatched/nb27rp-%A_%a.log
#SBATCH --array=0-26
#
# NewBench 27 -- RE-RUN with all 8 fruton_fix_* patches applied.
# --------------------------------------------------------------
# All 8 patches (geometric_gap_detection, phosaa14sb, forcefield_ptm,
# per_loop_plddt, atom_rename, cofactor_params, active_site_overrides,
# af_splice) applied to feature/modeller-af-hybrid on 2026-08-20.
#
# We re-run from step 8 (gap_detection) since that is the first stage
# whose behaviour changed. Steps 1-7 outputs from the earlier run
# (job 9145049) are preserved on the checkpoints.

set -euo pipefail

DATA_DIR=/mnt/netapp1/Store_othcxlwa/newbench_27
ISO=$DATA_DIR/isolated

PDB=$(awk -v i="$SLURM_ARRAY_TASK_ID" '$1==i{print $2}' "$ISO/pdb_index.txt")
if [[ -z "$PDB" ]]; then
    echo "ERROR: no PDB for array index $SLURM_ARRAY_TASK_ID" >&2
    exit 2
fi

export PATH="$HOME/bin:$PATH"
TASK_DIR="$ISO/$PDB"

mkdir -p /mnt/netapp1/Store_othcxlwa/newbench_27/logs/array_repatched

echo "SLURM_ARRAY_JOB_ID=$SLURM_ARRAY_JOB_ID  TASK_ID=$SLURM_ARRAY_TASK_ID"
echo "PDB       : $PDB"
echo "data_dir  : $TASK_DIR"
echo "fruton    : $(which fruton || echo MISSING)"
echo "start     : $(date -Iseconds)"
echo "from-step : 8 (post-8-patch rerun)"
echo

exec fruton --data-dir "$TASK_DIR" --from-step 8
