#!/usr/bin/env bash
# =============================================================================
# MCPB.py — Steps 2-4  |  3EWJ  |  site 3EWJ_ZN_B_2
# =============================================================================
# Creates step03_amber_params/, assembles all required inputs, then runs
# MCPB.py -s 2, -s 3, -s 4, and tleap.
#
# Run ONLY after all three Gaussian jobs in step02_gaussian/ have completed.
#
# Required files in step02_gaussian/ after Gaussian:
#   3EWJ_ZN_B_2_small_opt.fchk   (from: formchk _small_opt.chk)
#   3EWJ_ZN_B_2_small_fc.log
#   3EWJ_ZN_B_2_large_mk.log
#
# Final outputs in step03_amber_params/:
#   3EWJ_ZN_B_2.frcmod            AMBER bond/angle force constants
#   3EWJ_ZN_B_2_addmed.mol2       AMBER partial charges (RESP)
#   3EWJ_ZN_B_2_tleap.in / 3EWJ_ZN_B_2_tleap.out
# =============================================================================
set -euo pipefail

# -- Create step03 directory and assemble all required inputs ----------------
mkdir -p step03_amber_params

# All step01 outputs (model PDBs, .in, mol2, MCPB intermediate files)
cp step01_gen_inputs/* step03_amber_params/

# Gaussian outputs from step02
cp step02_gaussian/3EWJ_ZN_B_2_small_opt.fchk step03_amber_params/
cp step02_gaussian/3EWJ_ZN_B_2_small_fc.log   step03_amber_params/
cp step02_gaussian/3EWJ_ZN_B_2_large_mk.log   step03_amber_params/

# -- Run MCPB.py steps 2-4 and tleap ----------------------------------------
cd step03_amber_params

# Step 2: Derive bond/angle force constants (Seminario method).
MCPB.py -i 3EWJ.in -s 2

# Step 3: Fit RESP charges.
MCPB.py -i 3EWJ.in -s 3

# Step 4: Generate LEaP topology input.
MCPB.py -i 3EWJ.in -s 4
tleap -s -f 3EWJ_ZN_B_2_tleap.in > 3EWJ_ZN_B_2_tleap.out

echo ""
echo "Parametrization complete.  Key outputs in step03_amber_params/:"
echo "  3EWJ_ZN_B_2.frcmod"
echo "  3EWJ_ZN_B_2_addmed.mol2"
echo "  3EWJ_ZN_B_2_tleap.in / 3EWJ_ZN_B_2_tleap.out"
