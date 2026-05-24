#!/usr/bin/env bash
# =============================================================================
# MCPB.py — Steps 2-4  |  1W4R  |  site 1W4R_ZN_A_400
# =============================================================================
# Creates step03_amber_params/, assembles all required inputs, then runs
# MCPB.py -s 2, -s 3, -s 4, and tleap.
#
# Run ONLY after all three Gaussian jobs in step02_gaussian/ have completed.
#
# Required files in step02_gaussian/ after Gaussian:
#   1W4R_ZN_A_400_small_opt.fchk   (from: formchk _small_opt.chk)
#   1W4R_ZN_A_400_small_fc.log
#   1W4R_ZN_A_400_large_mk.log
#
# Final outputs in step03_amber_params/:
#   1W4R_ZN_A_400.frcmod            AMBER bond/angle force constants
#   1W4R_ZN_A_400_addmed.mol2       AMBER partial charges (RESP)
#   1W4R_ZN_A_400_tleap.in / 1W4R_ZN_A_400_tleap.out
# =============================================================================
set -euo pipefail

# -- Generate formatted checkpoint (required by MCPB.py -s 2) ----------------
# formchk must be available (module load cesga/2020 g16/c1 on CESGA login node)
command -v formchk >/dev/null || { echo "ERROR: formchk not found. Load g16 module first." >&2; exit 2; }
formchk step02_gaussian/1W4R_ZN_A_400_small_opt.chk \
        step02_gaussian/1W4R_ZN_A_400_small_opt.fchk

# -- Create step03 directory and assemble all required inputs ----------------
mkdir -p step03_amber_params

# All step01 outputs (model PDBs, .in, mol2, MCPB intermediate files)
cp step01_gen_inputs/* step03_amber_params/

# Gaussian outputs from step02
cp step02_gaussian/1W4R_ZN_A_400_small_opt.fchk step03_amber_params/
cp step02_gaussian/1W4R_ZN_A_400_small_fc.log   step03_amber_params/
cp step02_gaussian/1W4R_ZN_A_400_large_mk.log   step03_amber_params/

# -- Run MCPB.py steps 2-4 and tleap ----------------------------------------
cd step03_amber_params

# Step 2: Derive bond/angle force constants (Seminario method).
MCPB.py -i 1W4R.in -s 2

# Step 3: Fit RESP charges.
MCPB.py -i 1W4R.in -s 3

# Step 4: Generate LEaP topology input.
MCPB.py -i 1W4R.in -s 4
tleap -s -f 1W4R_ZN_A_400_tleap.in > 1W4R_ZN_A_400_tleap.out

echo ""
echo "Parametrization complete.  Key outputs in step03_amber_params/:"
echo "  1W4R_ZN_A_400.frcmod"
echo "  1W4R_ZN_A_400_addmed.mol2"
echo "  1W4R_ZN_A_400_tleap.in / 1W4R_ZN_A_400_tleap.out"
