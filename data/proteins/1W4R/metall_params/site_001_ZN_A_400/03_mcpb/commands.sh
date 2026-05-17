#!/usr/bin/env bash
set -euo pipefail

# Step 1: Generate Gaussian input files and model PDBs.
MCPB.py -i 1W4R.in -s 1

# --- Run Gaussian jobs after reviewing the generated .com files ---
# g16 < 1W4R_small_opt.com  > 1W4R_small_opt.log
# g16 < 1W4R_small_fc.com   > 1W4R_small_fc.log
# formchk 1W4R_small_opt.chk 1W4R_small_opt.fchk
# g16 < 1W4R_large_mk.com   > 1W4R_large_mk.log

# Step 2: Fit force constants from Gaussian Hessian (Seminario method).
# MCPB.py -i 1W4R.in -s 2

# Step 3: Fit RESP charges.
# MCPB.py -i 1W4R.in -s 3

# Step 4: Generate LEaP input.
# MCPB.py -i 1W4R.in -s 4
# tleap -s -f 1W4R_tleap.in > 1W4R_tleap.out
