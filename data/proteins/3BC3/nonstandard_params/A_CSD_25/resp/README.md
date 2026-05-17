# RESP Parametrization — CSD  |  3BC3 / A_CSD_25

## Overview

This directory contains the scaffold for deriving AMBER GAFF2 bonded parameters
and HF/6-31G* RESP charges for **CSD** embedded in the protein backbone.

RESP at HF/6-31G* is required (not AM1-BCC) for consistency with ff14SB, which
derives all its backbone and side-chain charges at the same level of theory.
Reference: Maier et al. JCTC 2015; Homeyer et al. J Mol Model 2006.

## Formal charge

```
Net charge : -1   (source: hardcoded_table)
Spin       : 1 (singlet — verify for open-shell modifications)
```

## Workflow

```
bash commands.sh
  → step01_gaussian_inputs/CSD_opt.com   (HF/6-31G* Opt)
  → step01_gaussian_inputs/CSD_esp.com   (HF/6-31G* MK ESP)
  → staged to step02_gaussian/

bash submit_gaussian.sh
  → sbatch both Gaussian jobs (ESP depends on opt via SLURM dependency)

bash run_after_gaussian.sh
  → antechamber -c resp          → step03_resp_params/CSD_capped.mol2
  → parmchk2                     → step03_resp_params/CSD.frcmod
  → strip ACE/NME                → step03_resp_params/CSD.mol2
```

## Using in tleap

```tcl
source leaprc.protein.ff14SB
source leaprc.gaff2
loadamberparams step03_resp_params/CSD.frcmod
CSD = loadmol2 step03_resp_params/CSD.mol2
saveoff CSD CSD.lib
# Then in your full system:
loadoff CSD.lib
```

## References

- Homeyer et al. J Mol Model 2006 — RESP/HF/6-31G* for pSer/pThr/pTyr/pHis
- Khoury et al. JCTC 2014 — Forcefield_PTM; 32 PTMs
- Sengupta et al. JCTC 2024 — phosaa14SB/phosaa19SB
- Maier et al. JCTC 2015 — ff14SB
