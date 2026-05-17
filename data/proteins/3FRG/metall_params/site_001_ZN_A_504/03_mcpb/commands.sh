#!/usr/bin/env bash
set -euo pipefail

# Step 1: Generate Gaussian input files and model PDBs.
MCPB.py -i 3FRG.in -s 1

# --- Run Gaussian jobs after reviewing the generated .com files ---
# g16 < 3FRG_small_opt.com  > 3FRG_small_opt.log
# g16 < 3FRG_small_fc.com   > 3FRG_small_fc.log
# formchk 3FRG_small_opt.chk 3FRG_small_opt.fchk
# g16 < 3FRG_large_mk.com   > 3FRG_large_mk.log

# Step 2: Fit force constants from Gaussian Hessian (Seminario method).
# MCPB.py -i 3FRG.in -s 2

# Step 3: Fit RESP charges.
# MCPB.py -i 3FRG.in -s 3

# Step 4: Generate LEaP input.
# MCPB.py -i 3FRG.in -s 4
# tleap -s -f 3FRG_tleap.in > 3FRG_tleap.out
