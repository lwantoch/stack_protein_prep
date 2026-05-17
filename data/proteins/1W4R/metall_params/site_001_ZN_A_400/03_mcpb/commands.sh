#!/usr/bin/env bash
# =============================================================================
# MCPB.py — Step 1  |  1W4R  |  site 1W4R_ZN_A_400
# =============================================================================
# Creates step01_gen_inputs/, copies input files, runs MCPB.py -s 1.
# Stages the three Gaussian .com files into step02_gaussian/.
#
# After this script:
#   1. Review .com files in step02_gaussian/.
#   2. Submit Gaussian jobs:  bash submit_gaussian.sh
#   3. Once all Gaussian jobs complete:  bash run_after_gaussian.sh
# =============================================================================
set -euo pipefail

# -- Create step01 directory and copy all required inputs --------------------
mkdir -p step01_gen_inputs
cp 1W4R.in 1W4R_mcpb.pdb ZN.mol2 step01_gen_inputs/

# -- Run MCPB.py step 1 (generates Gaussian input files and model PDBs) ------
cd step01_gen_inputs
MCPB.py -i 1W4R.in -s 1
cd ..

# -- Stage Gaussian input files into step02_gaussian/ ------------------------
mkdir -p step02_gaussian
cp step01_gen_inputs/1W4R_ZN_A_400_small_opt.com step02_gaussian/
cp step01_gen_inputs/1W4R_ZN_A_400_small_fc.com  step02_gaussian/
cp step01_gen_inputs/1W4R_ZN_A_400_large_mk.com  step02_gaussian/

echo ""
echo "Step 1 done.  Gaussian inputs staged in step02_gaussian/:"
echo "  1W4R_ZN_A_400_small_opt.com   geometry optimisation (small model)"
echo "  1W4R_ZN_A_400_small_fc.com    force constants / Hessian (small model)"
echo "  1W4R_ZN_A_400_large_mk.com    RESP charges (large model)"
echo ""
echo "Review each .com file, then run:  bash submit_gaussian.sh"
