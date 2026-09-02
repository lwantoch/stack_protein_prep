#!/bin/bash
#SBATCH -J fruton_5m7u_cat
#SBATCH -p medium
#SBATCH -N 1
#SBATCH -c 4
#SBATCH --mem 16G
#SBATCH -t 4:00:00
#SBATCH -o /mnt/netapp1/Store_othcxlwa/newbench_27/logs/5m7u_catdomain-%j.log
#SBATCH -e /mnt/netapp1/Store_othcxlwa/newbench_27/logs/5m7u_catdomain-%j.log
#
# 5M7U catalytic-domain-only prep (residues 60-395).
# ---------------------------------------------------
# Parallel companion to the full-range run (part of array 9162212, task 11).
# The 141-residue gap between resi 396 and 536 sits in the linker/helical
# bundle region -- OUTSIDE the catalytic pocket that carries the XHA ligand
# and the two catalytic aspartates D174/D175. Restricting to 60-395 skips
# that huge gap and avoids relying on AlphaFold for a 141-residue loop
# (mean plDDT 82.7, at the limit of reliable modelling).
#
# Remaining gaps in 60-395: 33-residue gap (338-370) + smaller ones -- all
# well within MODELLER's loop-modelling comfort range.
#
# Compare docking + GBSA scores between this trimmed variant and the
# full-range variant to see how much the 141-loop matters for the pocket.

set -euo pipefail

DATA_DIR=/mnt/netapp1/Store_othcxlwa/newbench_27/isolated_domain/5M7U_catdomain
export PATH="$HOME/bin:$PATH"

echo "SLURM_JOB_ID=$SLURM_JOB_ID  node=$SLURM_NODELIST"
echo "PDB       : 5M7U (catalytic domain only, range 60-395)"
echo "data_dir  : $DATA_DIR"
echo "fruton    : $(which fruton || echo MISSING)"
echo "start     : $(date -Iseconds)"
echo

exec fruton --data-dir "$DATA_DIR"
