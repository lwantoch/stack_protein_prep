# MCPB Parametrization — 3FRG

## What MCPB.py does

MCPB.py (Metal Center Parameter Builder) automates the bonded-model
parametrization of transition-metal centers for AMBER force fields. It:

1. Builds a "small model" (metal + first-shell donors) and a "large model"
   (small model + second-shell) from the input PDB.
2. Generates Gaussian input files for geometry optimization (small model),
   force-constant calculation (Hessian, small model), and RESP charge fitting
   (large model).
3. Uses the Seminario method to derive bond/angle force constants from the
   Gaussian Hessian.
4. Produces AMBER ``frcmod`` and ``mol2`` files for LEaP.

Reference: doi:10.1021/acs.jcim.5b00674

## Before running

Review and complete **all** ``TODO`` placeholders in ``3FRG.in``:

- ``ion_ids``         — PDB serial number of the metal atom in ``3FRG_mcpb.pdb``
- ``ion_mol2files``   — mol2 file(s) for the metal ion and non-standard ligands
- ``charge_m_sm``     — total charge and spin multiplicity for the small model
- ``charge_m_lm``     — total charge and spin multiplicity for the large model
- ``software_version``— Gaussian version available on your system (g09 or g16)
- ``cut_off``         — verify the coordination-shell cutoff is appropriate

## Gaussian

Gaussian is an external quantum-chemistry program and must be licensed and
installed separately. FRUTON does not run Gaussian automatically. After
completing step 1 (``MCPB.py -i 3FRG.in -s 1``), review the generated
``.com`` files before submitting them to Gaussian.

## Running

```bash
bash commands.sh
```
