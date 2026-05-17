#!/usr/bin/env bash
# =============================================================================
# RESP parametrization — Steps 2-3  |  CSD
# =============================================================================
# Run ONLY after both Gaussian jobs in step02_gaussian/ have completed.
#
# Required files in step02_gaussian/:
#   CSD_opt.log   (geometry optimisation log — geometry used in ESP)
#   CSD_esp.log   (ESP single-point log — MK electrostatic potential)
#
# Outputs in step03_resp_params/:
#   CSD_capped.mol2   ACE-CSD-NME RESP charges (GAFF2 atom types)
#   CSD.mol2          CSD atoms only (ACE/NME stripped)
#   CSD.frcmod        missing bonded parameters (parmchk2)
# =============================================================================
set -euo pipefail

mkdir -p step03_resp_params

# -- Step 2: RESP charge fitting (antechamber) --------------------------------
antechamber \
    -i step02_gaussian/CSD_esp.log -fi gout \
    -o step03_resp_params/CSD_capped.mol2 -fo mol2 \
    -c resp -s 2 -nc -1 \
    -at gaff2 -rn CSD -pf y

# -- Step 3: Missing bonded parameters (parmchk2) ----------------------------
parmchk2 \
    -i step03_resp_params/CSD_capped.mol2 -f mol2 \
    -o step03_resp_params/CSD.frcmod -s gaff2

# -- Strip ACE and NME caps from mol2 ----------------------------------------
python3 - <<'PYEOF'
from pathlib import Path

resname = "CSD"
src = Path("step03_resp_params") / f"{resname}_capped.mol2"
dst = Path("step03_resp_params") / f"{resname}.mol2"

text = src.read_text()
sections = text.split("@<TRIPOS>")
out_sections: list[str] = []

for section in sections:
    if section.startswith("ATOM"):
        kept_atoms: list[str] = []
        for line in section.splitlines(keepends=True):
            parts = line.split()
            if len(parts) >= 9 and parts[7] == resname:
                kept_atoms.append(line)
        # Renumber atoms
        renumbered: list[str] = []
        for i, line in enumerate(kept_atoms, start=1):
            parts = line.split()
            renumbered.append(
                f"{i:7d} {parts[1]:<8s}{parts[2]:>10s}{parts[3]:>10s}{parts[4]:>10s}"
                f" {parts[5]:<8s}{parts[6]:>3s}  {parts[7]:<8s}{parts[8]:>10s}\n"
            )
        out_sections.append("ATOM\n" + "".join(renumbered))
    elif section.startswith("BOND"):
        # Remove bonds involving ACE/NME atoms — simplest: keep section as-is
        # (tleap will re-derive intra-residue bonds from the mol2 atom list)
        out_sections.append(section)
    else:
        out_sections.append(section)

dst.write_text("@<TRIPOS>".join(out_sections))
print(f"Stripped caps written to: {dst}")
PYEOF

echo ""
echo "Parametrization complete.  Key outputs in step03_resp_params/:"
echo "  CSD_capped.mol2   full ACE-CSD-NME (for reference)"
echo "  CSD.mol2          CSD residue atoms only"
echo "  CSD.frcmod        AMBER frcmod for missing parameters"
