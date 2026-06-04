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
reduce ../capped_model_ace_nme.pdb > step01_gaussian_inputs/PTR_capped_H.pdb || true
[[ -s step01_gaussian_inputs/PTR_capped_H.pdb ]] || { echo "ERROR: reduce produced empty output" >&2; exit 1; }

# -- Adjust charge for electron parity ---------------------------------------
# Gaussian requires an even number of electrons for a closed-shell singlet.
# reduce may skip a H atom (e.g. backbone amide N-H due to clash), shifting
# the parity. Auto-correct: if (neutral_electrons + charge) is odd, add 1 to
# charge (less negative by 1) so the total is even.
CHARGE=-2
CHARGE=$(python3 - step01_gaussian_inputs/PTR_capped_H.pdb "$CHARGE" <<'PYEOF'
import sys
from pathlib import Path
_Z = {'H':1,'C':6,'N':7,'O':8,'F':9,'P':15,'S':16,'CL':17,'BR':35,'I':53,
       'ZN':30,'FE':26,'MG':12,'CA':20,'MN':25,'CU':29,'CO':27,'NI':28}
pdb_path, charge = Path(sys.argv[1]), int(sys.argv[2])
total_z = 0
for ln in pdb_path.read_text(errors='replace').splitlines():
    if not ln.startswith(('ATOM','HETATM')):
        continue
    elem = ln[76:78].strip().upper() if len(ln) > 76 else ''
    if not elem:
        elem = ''.join(c for c in ln[12:16].strip() if c.isalpha())[:2].upper()
    total_z += _Z.get(elem, 0)
total_e = total_z + (-charge)
if total_e % 2 != 0:
    charge += 1  # shift by +1 to get even count (less negative)
    import sys as _sys
    print(f"INFO: charge adjusted to {charge} for even electron count (neutral_Z={total_z})", file=_sys.stderr)
print(charge)
PYEOF
)
echo "Using charge: $CHARGE"

# -- Generate Gaussian input files via antechamber ---------------------------
cd step01_gaussian_inputs

# Geometry optimisation (HF/6-31G*)
antechamber \
    -i PTR_capped_H.pdb -fi pdb \
    -o PTR_opt.com -fo gcrt \
    -nc "$CHARGE" -s 2 -at gaff2 \
    -gn PTR_opt \
    -gk "HF/6-31G* Opt"

# ESP single-point on optimised geometry (MK charges)
antechamber \
    -i PTR_capped_H.pdb -fi pdb \
    -o PTR_esp.com -fo gcrt \
    -nc "$CHARGE" -s 2 -at gaff2 \
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
