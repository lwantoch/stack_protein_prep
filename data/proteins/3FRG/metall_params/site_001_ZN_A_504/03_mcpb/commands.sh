#!/usr/bin/env bash
# =============================================================================
# MCPB.py — Step 1  |  3FRG  |  site 3FRG_ZN_A_504
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
cp 3FRG.in 3FRG_mcpb.pdb ZN.mol2 step01_gen_inputs/

# -- Run MCPB.py step 1 (generates Gaussian input files and model PDBs) ------
cd step01_gen_inputs
MCPB.py -i 3FRG.in -s 1
cd ..

# -- Stage Gaussian input files into step02_gaussian/ ------------------------
mkdir -p step02_gaussian
cp step01_gen_inputs/3FRG_ZN_A_504_small_opt.com step02_gaussian/
cp step01_gen_inputs/3FRG_ZN_A_504_small_fc.com  step02_gaussian/
cp step01_gen_inputs/3FRG_ZN_A_504_large_mk.com  step02_gaussian/

echo ""
echo "Step 1 done.  Gaussian inputs staged in step02_gaussian/:"
echo "  3FRG_ZN_A_504_small_opt.com   geometry optimisation (small model)"
echo "  3FRG_ZN_A_504_small_fc.com    force constants / Hessian (small model)"
echo "  3FRG_ZN_A_504_large_mk.com    RESP charges (large model)"
echo ""
echo "Review each .com file, then run:  bash submit_gaussian.sh"
