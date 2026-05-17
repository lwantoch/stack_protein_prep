set -euo pipefail

SDF_IN="${1:?Usage: ./filter_ki_ligands.sh input.sdf [output_dir] [threshold_nM]}"
OUT_DIR="${2:-ligands_under_1000nM}"
THRESHOLD_NM="${3:-1000}"

mkdir -p "$OUT_DIR"

python - "$SDF_IN" "$OUT_DIR" "$THRESHOLD_NM" <<'PY'
from pathlib import Path
import re
import sys

from rdkit import Chem

sdf_in = Path(sys.argv[1])
out_dir = Path(sys.argv[2])
threshold_nm = float(sys.argv[3])

# Accept fields such as:
# Ki, KI, Ki_nM, experimental_Ki, exp_ki, binding_Ki, etc.
KI_FIELD_RE = re.compile(r"(?:^|[_\-\s])(exp(?:erimental)?[_\-\s]*)?ki(?:[_\-\s]|$)", re.I)
VALUE_RE = re.compile(r"([<>≤≥]?)\s*([0-9]*\.?[0-9]+(?:[eE][-+]?[0-9]+)?)\s*(pM|nM|uM|µM|mM|M)?", re.I)

UNIT_TO_NM = {
    None: 1.0,
    "pm": 1e-3,
    "nm": 1.0,
    "um": 1e3,
    "µm": 1e3,
    "mm": 1e6,
    "m": 1e9,
}

def parse_ki_nm(value: str) -> float | None:
    match = VALUE_RE.search(value)
    if not match:
        return None

    sign, number, unit = match.groups()
    number_nm = float(number) * UNIT_TO_NM.get(unit.lower() if unit else None, 1.0)

    # "<1000 nM" should pass if 1000 is at/below threshold.
    # ">1000 nM" should not pass because the true value may be larger.
    if sign in {">", "≥"}:
        return None

    return number_nm

def get_ki_nm(mol) -> float | None:
    for prop_name in mol.GetPropNames():
        if not KI_FIELD_RE.search(prop_name):
            continue
        ki_nm = parse_ki_nm(mol.GetProp(prop_name))
        if ki_nm is not None:
            return ki_nm
    return None

supplier = Chem.SDMolSupplier(str(sdf_in), removeHs=False)
kept = 0
seen = 0

for mol in supplier:
    if mol is None:
        continue

    seen += 1
    ki_nm = get_ki_nm(mol)

    if ki_nm is None or ki_nm >= threshold_nm:
        continue

    kept += 1
    out_path = out_dir / f"ligand{kept:02d}.sdf"
    writer = Chem.SDWriter(str(out_path))
    writer.write(mol)
    writer.close()

print(f"read: {seen}")
print(f"kept: {kept}")
print(f"output_dir: {out_dir}")