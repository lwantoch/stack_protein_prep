#!/usr/bin/env bash
# =============================================================================
# RESP parametrization — Step 1  |  PTR
# =============================================================================
# Adds hydrogens with reduce, then uses antechamber to generate Gaussian
# input files for HF/6-31G* geometry optimisation and ESP calculation.
# Stages the .com files into step02_gaussian/ ready for HPC submission.
#
# After this script:
#   1. Review .com files in step02_gaussian/.
#   2. Submit Gaussian jobs:  bash submit_gaussian.sh
#   3. After Gaussian:  bash run_after_gaussian.sh
# =============================================================================
set -euo pipefail

mkdir -p step01_gaussian_inputs

# -- Add hydrogens -----------------------------------------------------------
reduce ../capped_model_ace_nme.pdb > step01_gaussian_inputs/PTR_capped_H.pdb

# -- Generate Gaussian input files via antechamber ---------------------------
cd step01_gaussian_inputs

# Geometry optimisation (HF/6-31G*)
antechamber \
    -i PTR_capped_H.pdb -fi pdb \
    -o PTR_opt.com -fo gcrt \
    -nc -2 -s 2 -at gaff2 \
    -gn PTR_opt \
    -gk "HF/6-31G* Opt"

# ESP single-point on optimised geometry (MK charges)
antechamber \
    -i PTR_capped_H.pdb -fi pdb \
    -o PTR_esp.com -fo gcrt \
    -nc -2 -s 2 -at gaff2 \
    -gn PTR_esp \
    -gk "HF/6-31G* Pop=MK IOp(6/33=2,6/42=6)"

cd ..

# -- Stage for Gaussian submission -------------------------------------------
mkdir -p step02_gaussian
cp step01_gaussian_inputs/PTR_opt.com step02_gaussian/
cp step01_gaussian_inputs/PTR_esp.com step02_gaussian/

echo ""
echo "Step 1 done.  Gaussian inputs staged in step02_gaussian/:"
echo "  PTR_opt.com   HF/6-31G* geometry optimisation"
echo "  PTR_esp.com   HF/6-31G* ESP (MK charges)"
echo ""
echo "Review each .com file, then run:  bash submit_gaussian.sh"
