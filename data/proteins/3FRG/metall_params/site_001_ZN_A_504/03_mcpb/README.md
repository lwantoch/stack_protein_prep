# MCPB Parametrization — 3FRG  |  site 3FRG_ZN_A_504

## Overview

MCPB.py (Metal Center Parameter Builder) derives bonded-model AMBER parameters
for transition-metal centers using quantum-chemical calculations.

Reference: doi:10.1021/acs.jcim.5b00674

## Directory layout

```
03_mcpb/                         ← run all scripts from here
├── 3FRG.in
├── 3FRG_mcpb.pdb
├── *.mol2
├── commands.sh           Step 1
├── submit_gaussian.sh    Gaussian HPC submission
├── run_after_gaussian.sh Steps 2-4
│
├── step01_gen_inputs/    created by commands.sh
│   ├── 3FRG_ZN_A_504_small_opt.com
│   ├── 3FRG_ZN_A_504_small_fc.com
│   ├── 3FRG_ZN_A_504_large_mk.com
│   └── 3FRG_ZN_A_504_small.pdb / _large.pdb  (model PDBs)
│
├── step02_gaussian/      .com files staged here; Gaussian logs written here
│   ├── 3FRG_ZN_A_504_small_opt.log / .fchk
│   ├── 3FRG_ZN_A_504_small_fc.log
│   └── 3FRG_ZN_A_504_large_mk.log
│
└── step03_amber_params/  created by run_after_gaussian.sh; final outputs here
    ├── 3FRG_ZN_A_504.frcmod
    ├── 3FRG_ZN_A_504_addmed.mol2
    └── 3FRG_ZN_A_504_tleap.in / _tleap.out
```

## Workflow (3 scripts, 4 MCPB.py invocations)

```
bash commands.sh
  → mkdir step01_gen_inputs/, run MCPB.py -s 1
  → mkdir step02_gaussian/, copy .com files

bash submit_gaussian.sh
  → sbatch all three Gaussian jobs from step02_gaussian/

bash run_after_gaussian.sh
  → mkdir step03_amber_params/, merge step01 + step02 outputs
  → MCPB.py -s 2, -s 3, -s 4, tleap
```

## Before running Step 1

Review ``3FRG.in``.  Fields auto-filled by FRUTON are marked in comments.
Fields that require manual chemistry input carry ``TODO`` placeholders:

- ``software_version`` — Gaussian version on your cluster (g09 or g16)
- ``charge_m_sm / charge_m_lm`` — total charge + spin multiplicity of the QM
  region; spin is unambiguous only for closed-shell d10 ions (Zn²⁺, Cu⁺)

## Gaussian

Gaussian is a licensed external program; FRUTON does not run it automatically.
``submit_gaussian.sh`` provides a SLURM template — adjust ``--ntasks``,
``--mem``, and ``--time`` to match your cluster's policies.

All three Gaussian jobs are independent and can be submitted simultaneously.
Steps 2-4 require all three to finish before starting.
