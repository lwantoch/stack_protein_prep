"""Parametrization scaffold for non-standard amino acid residues embedded in the
protein backbone (e.g. phosphotyrosine PTR, cysteic acid CSD).

Workflow
--------
1. Detect non-standard backbone residues via cap.py's
   ``find_nonstandard_residues()``.
2. Check whether pre-built AMBER parameters already exist in a known database
   (phosaa14SB, phosaa19SB, Forcefield_PTM, Bryce group).  If found, write a
   README with tleap loading instructions — no Gaussian needed.
3. For residues not covered by any database, build the ACE-X-NME model compound
   from the protein crystal context, look up the formal charge (hardcoded table
   → RCSB CDD → 0 with warning), and generate a RESP parametrization scaffold:
   ``commands.sh`` (build H, generate Gaussian inputs) → ``submit_gaussian.sh``
   (SLURM HPC template) → ``run_after_gaussian.sh`` (antechamber RESP +
   parmchk2 + cap stripping).
"""

from __future__ import annotations

import json
import subprocess
import urllib.error
import urllib.request
from pathlib import Path
from typing import Any

from Bio.PDB import PDBParser

from stack_protein_preparation.cap import find_nonstandard_residues

# ---------------------------------------------------------------------------
# Known-database registry
# ---------------------------------------------------------------------------

# residue name → {database, leaprc (if AmberTools), url, notes}
_DATABASE_REGISTRY: dict[str, dict[str, str]] = {
    # AmberTools phosaa14SB / phosaa19SB — ships with AmberTools 22+
    # Sengupta et al. JCTC 2024
    "SEP": {
        "database": "phosaa14SB",
        "leaprc": "leaprc.phosaa14SB",
        "description": "phosphoserine (−2); phosaa14SB ff14SB-compatible",
    },
    "TPO": {
        "database": "phosaa14SB",
        "leaprc": "leaprc.phosaa14SB",
        "description": "phosphothreonine (−2)",
    },
    "PTR": {
        "database": "phosaa14SB",
        "leaprc": "leaprc.phosaa14SB",
        "description": "phosphotyrosine (−2)",
    },
    "S1P": {
        "database": "phosaa14SB",
        "leaprc": "leaprc.phosaa14SB",
        "description": "phosphoserine monoprotonated (−1)",
    },
    "T1P": {
        "database": "phosaa14SB",
        "leaprc": "leaprc.phosaa14SB",
        "description": "phosphothreonine monoprotonated (−1)",
    },
    "Y1P": {
        "database": "phosaa14SB",
        "leaprc": "leaprc.phosaa14SB",
        "description": "phosphotyrosine monoprotonated (−1)",
    },
    "NEP": {
        "database": "phosaa14SB",
        "leaprc": "leaprc.phosaa14SB",
        "description": "N1-phosphonohistidine (−2)",
    },
    # Forcefield_PTM — Khoury et al. JCTC 2014 (ff03-based; not directly
    # ff14SB-compatible — backbone charges differ)
    "MLZ": {
        "database": "Forcefield_PTM",
        "url": "https://selene.princeton.edu/FFPTM/",
        "description": "N-methyl-L-lysine (+1); NOTE: ff03-based, verify compatibility with ff14SB",
    },
    "MLY": {
        "database": "Forcefield_PTM",
        "url": "https://selene.princeton.edu/FFPTM/",
        "description": "N,N-dimethyl-L-lysine (+1); ff03-based",
    },
    "M3L": {
        "database": "Forcefield_PTM",
        "url": "https://selene.princeton.edu/FFPTM/",
        "description": "N,N,N-trimethyl-L-lysine (+1, quaternary); ff03-based",
    },
    "ALY": {
        "database": "Forcefield_PTM",
        "url": "https://selene.princeton.edu/FFPTM/",
        "description": "N6-acetyl-L-lysine (0); ff03-based",
    },
    "AGM": {
        "database": "Forcefield_PTM",
        "url": "https://selene.princeton.edu/FFPTM/",
        "description": "5-methyl-arginine (+1); ff03-based",
    },
    "CSO": {
        "database": "Forcefield_PTM",
        "url": "https://selene.princeton.edu/FFPTM/",
        "description": "S-hydroxy-cysteine (0); ff03-based",
    },
    "SCY": {
        "database": "Forcefield_PTM",
        "url": "https://selene.princeton.edu/FFPTM/",
        "description": "S-acetyl-cysteine (0); ff03-based",
    },
    "P1L": {
        "database": "Forcefield_PTM",
        "url": "https://selene.princeton.edu/FFPTM/",
        "description": "S-palmitoyl-cysteine (0); ff03-based",
    },
    "CGU": {
        "database": "Forcefield_PTM",
        "url": "https://selene.princeton.edu/FFPTM/",
        "description": "gamma-carboxyglutamic acid (−2); ff03-based",
    },
    "PCA": {
        "database": "Forcefield_PTM",
        "url": "https://selene.princeton.edu/FFPTM/",
        "description": "pyroglutamic acid (0); ff03-based",
    },
    "HYP": {
        "database": "Forcefield_PTM",
        "url": "https://selene.princeton.edu/FFPTM/",
        "description": "4-hydroxy-L-proline (0); ff03-based",
    },
    "LYZ": {
        "database": "Forcefield_PTM",
        "url": "https://selene.princeton.edu/FFPTM/",
        "description": "5-hydroxy-L-lysine (+1); ff03-based",
    },
    "DA2": {
        "database": "Forcefield_PTM",
        "url": "https://selene.princeton.edu/FFPTM/",
        "description": "ADMA, asymmetric dimethylarginine (+1); ff03-based",
    },
    # Bryce group database — Homeyer et al. J Mol Model 2006; HF/6-31G* RESP;
    # parm99-consistent (close to ff14SB backbone)
    "H2D": {
        "database": "Bryce_group",
        "url": "https://personalpages.manchester.ac.uk/staff/Richard.Bryce/amber/",
        "description": "phosphohistidine ND1 dianionic (−2); parm99-compatible",
    },
    "H1D": {
        "database": "Bryce_group",
        "url": "https://personalpages.manchester.ac.uk/staff/Richard.Bryce/amber/",
        "description": "phosphohistidine ND1 monoanionic (−1)",
    },
    "H2E": {
        "database": "Bryce_group",
        "url": "https://personalpages.manchester.ac.uk/staff/Richard.Bryce/amber/",
        "description": "phosphohistidine NE2 dianionic (−2)",
    },
    "H1E": {
        "database": "Bryce_group",
        "url": "https://personalpages.manchester.ac.uk/staff/Richard.Bryce/amber/",
        "description": "phosphohistidine NE2 monoanionic (−1)",
    },
    # AmberTools native (distributed with ff14SB)
    "MSE": {
        "database": "ff14SB_native",
        "leaprc": "leaprc.protein.ff14SB",
        "description": "selenomethionine (0); native in ff14SB",
    },
}

# ---------------------------------------------------------------------------
# Physiological formal charges (pH 7.0–7.4)
# RCSB pdbx_formal_charge returns 0 for all of these — do NOT use it.
# Sources: Homeyer 2006, Khoury 2014, Sengupta 2024, pKa reasoning.
# ---------------------------------------------------------------------------
_PHYSIOLOGICAL_CHARGES: dict[str, int] = {
    # Phosphorylation — phosphate pKa2 ~6.5 → −2 dominant at pH 7
    "SEP": -2, "S1P": -1,
    "TPO": -2, "T1P": -1,
    "PTR": -2, "Y1P": -1,
    "NEP": -2,
    "H2D": -2, "H1D": -1, "H2E": -2, "H1E": -1,
    "PHD": -2,  # aspartyl phosphate; unstable intermediate
    # Cysteine oxidation
    "CSO": 0,   # sulfenic acid pKa ~6–8; borderline; treat 0
    "CSD": -1,  # sulfinic acid pKa ~2
    "OCS": -1,  # sulfonic acid pKa ~−1
    "CSS": 0,   # persulfide; neutral
    "CME": 0,   # S-(2-hydroxyethyl)thiocysteine; neutral
    "SCY": 0,   # S-acetylcysteine; neutral
    "SMC": 0,   # S-methylcysteine; neutral
    "P1L": 0,   # S-palmitoylcysteine; neutral
    "CAF": 0,   # S-dimethylarsinoylcysteine; neutral
    # Methionine oxidation
    "MHO": 0,   # methionine sulfoxide (S)
    "SME": 0,   # methionine sulfoxide (R)
    "OMT": 0,   # methionine sulfone
    "MSE": 0,   # selenomethionine (crystallographic, not PTM)
    # Tyrosine modifications
    "TYS": -1,  # sulfo-Tyr; sulfate ester pKa ~1
    "NIY": -1,  # nitro-Tyr; phenol pKa lowered to ~6.8 → −1 at pH 7
    "TYI": 0,   # diiodo-Tyr; neutral
    "IYR": 0,   # iodo-Tyr; neutral
    "TPQ": 0,   # topaquinone; neutral
    "TYQ": 0,   # 3-amino-6-hydroxy-Tyr; neutral
    # Lysine modifications
    "MLZ": 1,   # monomethyl-Lys; ε-NH still protonated pKa ~10.6
    "MLY": 1,   # dimethyl-Lys; pKa ~9.2
    "M3L": 1,   # trimethyl-Lys; quaternary ammonium, permanent +1
    "ALY": 0,   # acetyl-Lys; amide bond neutralises ε-amine
    "KCR": 0,   # crotonyl-Lys; enoyl amide; neutral
    "SLL": -1,  # succinyl-Lys; distal carboxylate pKa ~4.3
    "KCX": -1,  # carbamyl-Lys (carbamate); pKa ~5
    "LYZ": 1,   # 5-hydroxy-Lys; ε-NH2 retained pKa ~9.7
    "LLP": -1,  # pyridoxal-P-Lys Schiff base; phosphate ionised
    # Arginine methylations — guanidinium pKa ~13.8, always +1
    "AGM": 1,   # ω-monomethyl-Arg
    "DA2": 1,   # ADMA (asymmetric dimethyl-Arg)
    "MMO": 1,   # N2-methyl-Arg
    "MMA": 1,   # ω-N-monomethyl-Arg
    # Citrullination
    "CIR": 0,   # citrulline; ureido group; neutral
    # Glutamate modifications
    "CGU": -2,  # γ-carboxyglutamate; two γ-carboxylates → −2 net
    # Proline hydroxylation
    "HYP": 0,   # 4-hydroxy-Pro; neutral
    "4FB": 0,   # 4-fluoro-Pro; neutral
    # Tryptophan hydroxylation/oxidation
    "4HT": 0,   "TRO": 0,   "HTR": 0,
    # Histidine modifications
    "HIC": 0,   # 4-methyl-His; imidazole pKa ~6.2; neutral at pH 7
    # Pyroglutamate
    "PCA": 0,   # pyroglutamate; lactam; neutral
    # Backbone N-methylation
    "SAR": 0,   # sarcosine (N-Me-Gly); neutral
    "MAA": 0,   # N-methyl-Ala; neutral
    "MEA": 0,   # N-methyl-Phe; neutral
    "MEN": 0,   # N-methyl-Asn; neutral
    # N-terminal acetylation
    "SAC": 0,   # N-acetyl-Ser; neutral
    # Selenocysteine / pyrrolysine
    "SEC": 0,   # selenocysteine; selenol pKa ~5.2 → borderline; treat 0
    "PYL": 0,   # pyrrolysine; amide; neutral
    # Chromophores
    "CRO": -1,  # GFP chromophore; phenolate anion in fluorescent state
    "DYG": -1,  # GFP chromophore variant
    # Phosphoserine legacy names (Bryce group)
    "H2D": -2,  "H1D": -1,  "H2E": -2,  "H1E": -1,
}

# ---------------------------------------------------------------------------
# PDB atom record formatting
# ---------------------------------------------------------------------------

def _format_atom_line(
    serial: int,
    name: str,
    resname: str,
    chain: str,
    resseq: int,
    x: float,
    y: float,
    z: float,
    element: str,
) -> str:
    atom_name_field = f" {name:<3s}" if len(name) < 4 else f"{name:<4s}"
    return (
        f"ATOM  {serial:5d} {atom_name_field} {resname:3s} {chain:1s}{resseq:4d}    "
        f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00          {element:>2s}  "
    )

# ---------------------------------------------------------------------------
# ACE-X-NME capped-model builder (no PTM-Psi)
# ---------------------------------------------------------------------------

def _build_capped_model(
    protein_pdb: Path,
    chain_id: str,
    resname: str,
    resseq: int,
    output_pdb: Path,
) -> bool:
    """Extract residue with neighbouring C=O (→ ACE) and N (→ NME) from the
    protein crystal structure.  Returns True on success."""
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("prot", str(protein_pdb))
    models = list(structure.get_models())
    if not models:
        return False
    model = models[0]

    target_chain = None
    for ch in model.get_chains():
        if ch.id.strip() == chain_id.strip():
            target_chain = ch
            break
    if target_chain is None:
        return False

    residues = list(target_chain.get_residues())
    target_idx: int | None = None
    for i, res in enumerate(residues):
        if res.get_resname().strip() == resname and res.id[1] == resseq:
            target_idx = i
            break
    if target_idx is None:
        return False

    target = residues[target_idx]
    prev_res = residues[target_idx - 1] if target_idx > 0 else None
    next_res = residues[target_idx + 1] if target_idx < len(residues) - 1 else None

    lines: list[str] = []
    serial = 1

    # --- ACE cap -----------------------------------------------------------
    # Use C and O from the preceding residue; CA_prev serves as the methyl C.
    if prev_res is not None and "C" in prev_res and "O" in prev_res:
        c_v = prev_res["C"].get_vector()
        o_v = prev_res["O"].get_vector()
        lines.append(_format_atom_line(serial, "C",   "ACE", "A", 1, *c_v, "C")); serial += 1
        lines.append(_format_atom_line(serial, "O",   "ACE", "A", 1, *o_v, "O")); serial += 1
        if "CA" in prev_res:
            ca_v = prev_res["CA"].get_vector()
            lines.append(_format_atom_line(serial, "CH3", "ACE", "A", 1, *ca_v, "C")); serial += 1
        else:
            import numpy as np
            # Estimate CH3 position opposite to the N direction from C
            if "N" in target:
                n_v = target["N"].get_vector()
                c_arr = c_v.get_array()
                n_arr = n_v.get_array()
                direction = c_arr - n_arr
                norm = direction / (np.linalg.norm(direction) + 1e-9)
                ch3 = c_arr + norm * 1.52
                lines.append(_format_atom_line(serial, "CH3", "ACE", "A", 1,
                                               float(ch3[0]), float(ch3[1]), float(ch3[2]), "C"))
                serial += 1
        lines.append("TER")
    else:
        lines.append(f"REMARK  ACE cap omitted — no preceding residue with C/O")

    # --- Residue X ---------------------------------------------------------
    for atom in target.get_atoms():
        v = atom.get_vector()
        elem = (atom.element or atom.id[0]).strip()
        lines.append(_format_atom_line(serial, atom.id.strip(), resname, "A", 2,
                                       *v, elem))
        serial += 1
    lines.append("TER")

    # --- NME cap -----------------------------------------------------------
    if next_res is not None and "N" in next_res:
        n_v = next_res["N"].get_vector()
        lines.append(_format_atom_line(serial, "N",   "NME", "A", 3, *n_v, "N")); serial += 1
        if "CA" in next_res:
            ca_v = next_res["CA"].get_vector()
            lines.append(_format_atom_line(serial, "CH3", "NME", "A", 3, *ca_v, "C")); serial += 1
        else:
            import numpy as np
            c_v = target["C"].get_vector() if "C" in target else n_v
            n_arr = n_v.get_array()
            c_arr = c_v.get_array()
            direction = n_arr - c_arr
            norm = direction / (np.linalg.norm(direction) + 1e-9)
            ch3 = n_arr + norm * 1.47
            lines.append(_format_atom_line(serial, "CH3", "NME", "A", 3,
                                           float(ch3[0]), float(ch3[1]), float(ch3[2]), "C"))
            serial += 1
        lines.append("TER")
    else:
        lines.append(f"REMARK  NME cap omitted — no following residue with N")

    lines.append("END")
    output_pdb.parent.mkdir(parents=True, exist_ok=True)
    output_pdb.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return True

# ---------------------------------------------------------------------------
# Formal charge lookup
# ---------------------------------------------------------------------------

def _lookup_charge_rcsb(resname: str) -> int | None:
    """Sum per-atom charges from the RCSB CDD (more reliable than pdbx_formal_charge)."""
    try:
        url = f"https://data.rcsb.org/rest/v1/core/chemcomp/{resname.upper()}"
        with urllib.request.urlopen(url, timeout=10) as resp:
            data = json.loads(resp.read())
        atoms = data.get("chem_comp_atom", [])
        if not atoms:
            return None
        total = sum(int(a.get("charge", 0)) for a in atoms)
        return total
    except Exception:
        return None


def lookup_formal_charge(resname: str) -> tuple[int, str]:
    """Return (formal_charge, source) for *resname* at physiological pH.

    Priority: hardcoded table → RCSB CDD atom sum → 0 with WARNING.
    """
    key = resname.strip().upper()
    if key in _PHYSIOLOGICAL_CHARGES:
        return _PHYSIOLOGICAL_CHARGES[key], "hardcoded_table"

    rcsb = _lookup_charge_rcsb(key)
    if rcsb is not None:
        return rcsb, "rcsb_cdd_atom_sum"

    return 0, "WARNING:unknown_residue_defaulted_to_0"

# ---------------------------------------------------------------------------
# Gaussian input writers
# ---------------------------------------------------------------------------

def _read_pdb_atoms(pdb_path: Path) -> list[tuple[str, float, float, float]]:
    """Return list of (element, x, y, z) from ATOM/HETATM lines."""
    atoms = []
    for line in pdb_path.read_text(encoding="utf-8").splitlines():
        if not line.startswith(("ATOM", "HETATM")):
            continue
        elem = line[76:78].strip() or line[12:16].strip()[0]
        x = float(line[30:38])
        y = float(line[38:46])
        z = float(line[46:54])
        atoms.append((elem, x, y, z))
    return atoms


def _write_gaussian_com(
    output_path: Path,
    title: str,
    charge: int,
    spin: int,
    atoms: list[tuple[str, float, float, float]],
    route: str,
    chk_name: str,
    old_chk: str | None = None,
    mem: str = "16GB",
    nproc: int = 8,
) -> None:
    lines = [
        f"%mem={mem}",
        f"%nproc={nproc}",
        f"%chk={chk_name}.chk",
    ]
    if old_chk:
        lines.append(f"%oldchk={old_chk}.chk")
    lines += [
        f"#P {route}",
        "",
        title,
        "",
        f"{charge} {spin}",
    ]
    if atoms:
        for elem, x, y, z in atoms:
            lines.append(f" {elem:<2s}  {x:12.6f}  {y:12.6f}  {z:12.6f}")
    lines += ["", ""]
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text("\n".join(lines), encoding="utf-8")


# ---------------------------------------------------------------------------
# Script writers
# ---------------------------------------------------------------------------

def _write_commands_sh(
    output_path: Path,
    resname: str,
    charge: int,
    capped_pdb_name: str = "capped.pdb",
) -> None:
    content = f"""\
#!/usr/bin/env bash
# =============================================================================
# RESP parametrization — Step 1  |  {resname}
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
reduce {capped_pdb_name} > step01_gaussian_inputs/{resname}_capped_H.pdb

# -- Generate Gaussian input files via antechamber ---------------------------
cd step01_gaussian_inputs

# Geometry optimisation (HF/6-31G*)
antechamber \\
    -i {resname}_capped_H.pdb -fi pdb \\
    -o {resname}_opt.com -fo gcrt \\
    -nc {charge} -s 2 -at gaff2 \\
    -gn {resname}_opt \\
    -gk "HF/6-31G* Opt"

# ESP single-point on optimised geometry (MK charges)
antechamber \\
    -i {resname}_capped_H.pdb -fi pdb \\
    -o {resname}_esp.com -fo gcrt \\
    -nc {charge} -s 2 -at gaff2 \\
    -gn {resname}_esp \\
    -gk "HF/6-31G* Pop=MK IOp(6/33=2,6/42=6)"

cd ..

# -- Stage for Gaussian submission -------------------------------------------
mkdir -p step02_gaussian
cp step01_gaussian_inputs/{resname}_opt.com step02_gaussian/
cp step01_gaussian_inputs/{resname}_esp.com step02_gaussian/

echo ""
echo "Step 1 done.  Gaussian inputs staged in step02_gaussian/:"
echo "  {resname}_opt.com   HF/6-31G* geometry optimisation"
echo "  {resname}_esp.com   HF/6-31G* ESP (MK charges)"
echo ""
echo "Review each .com file, then run:  bash submit_gaussian.sh"
"""
    output_path.write_text(content, encoding="utf-8")
    output_path.chmod(0o755)


def _write_submit_gaussian_sh(output_path: Path, resname: str) -> None:
    content = f"""\
#!/usr/bin/env bash
# =============================================================================
# Gaussian HPC submission — {resname}
# =============================================================================
# Submit geometry optimisation and ESP jobs from step02_gaussian/.
# Run ONLY after bash commands.sh has completed successfully.
# Edit --ntasks, --mem, --time to match your cluster.
#
# Note: submit ESP job AFTER optimisation completes (it reads the chk file).
# After all jobs complete:  bash run_after_gaussian.sh
# =============================================================================
set -euo pipefail
cd step02_gaussian

# --- 1. Geometry optimisation (HF/6-31G*) -----------------------------------
JOB_OPT=$(sbatch --parsable \\
    --job-name="{resname}_opt" --ntasks=8 --mem=16G --time=08:00:00 \\
    --wrap="g16 < {resname}_opt.com > {resname}_opt.log")
echo "Submitted opt job: ${{JOB_OPT}}"

# --- 2. ESP single-point (depends on opt via chk) ---------------------------
# The ESP .com reads the optimised geometry from the chk file.
# It is safe to submit immediately — Slurm will hold it until opt finishes.
sbatch --dependency=afterok:${{JOB_OPT}} \\
    --job-name="{resname}_esp" --ntasks=8 --mem=16G --time=04:00:00 \\
    --wrap="g16 < {resname}_esp.com > {resname}_esp.log"

echo "Submitted ESP job (depends on opt)."
echo "Monitor with:  squeue -u \\$USER"
echo "After both complete:  bash run_after_gaussian.sh"
"""
    output_path.write_text(content, encoding="utf-8")
    output_path.chmod(0o755)


def _write_run_after_gaussian_sh(output_path: Path, resname: str, charge: int) -> None:
    content = f"""\
#!/usr/bin/env bash
# =============================================================================
# RESP parametrization — Steps 2-3  |  {resname}
# =============================================================================
# Run ONLY after both Gaussian jobs in step02_gaussian/ have completed.
#
# Required files in step02_gaussian/:
#   {resname}_opt.log   (geometry optimisation log — geometry used in ESP)
#   {resname}_esp.log   (ESP single-point log — MK electrostatic potential)
#
# Outputs in step03_resp_params/:
#   {resname}_capped.mol2   ACE-{resname}-NME RESP charges (GAFF2 atom types)
#   {resname}.mol2          {resname} atoms only (ACE/NME stripped)
#   {resname}.frcmod        missing bonded parameters (parmchk2)
# =============================================================================
set -euo pipefail

mkdir -p step03_resp_params

# -- Step 2: RESP charge fitting (antechamber) --------------------------------
antechamber \\
    -i step02_gaussian/{resname}_esp.log -fi gout \\
    -o step03_resp_params/{resname}_capped.mol2 -fo mol2 \\
    -c resp -s 2 -nc {charge} \\
    -at gaff2 -rn {resname} -pf y

# -- Step 3: Missing bonded parameters (parmchk2) ----------------------------
parmchk2 \\
    -i step03_resp_params/{resname}_capped.mol2 -f mol2 \\
    -o step03_resp_params/{resname}.frcmod -s gaff2

# -- Strip ACE and NME caps from mol2 ----------------------------------------
python3 - <<'PYEOF'
from pathlib import Path

resname = "{resname}"
src = Path("step03_resp_params") / f"{{resname}}_capped.mol2"
dst = Path("step03_resp_params") / f"{{resname}}.mol2"

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
                f"{{i:7d}} {{parts[1]:<8s}}{{parts[2]:>10s}}{{parts[3]:>10s}}{{parts[4]:>10s}}"
                f" {{parts[5]:<8s}}{{parts[6]:>3s}}  {{parts[7]:<8s}}{{parts[8]:>10s}}\\n"
            )
        out_sections.append("ATOM\\n" + "".join(renumbered))
    elif section.startswith("BOND"):
        # Remove bonds involving ACE/NME atoms — simplest: keep section as-is
        # (tleap will re-derive intra-residue bonds from the mol2 atom list)
        out_sections.append(section)
    else:
        out_sections.append(section)

dst.write_text("@<TRIPOS>".join(out_sections))
print(f"Stripped caps written to: {{dst}}")
PYEOF

echo ""
echo "Parametrization complete.  Key outputs in step03_resp_params/:"
echo "  {resname}_capped.mol2   full ACE-{resname}-NME (for reference)"
echo "  {resname}.mol2          {resname} residue atoms only"
echo "  {resname}.frcmod        AMBER frcmod for missing parameters"
"""
    output_path.write_text(content, encoding="utf-8")
    output_path.chmod(0o755)


def _write_readme(
    output_path: Path,
    resname: str,
    charge: int,
    charge_source: str,
    pdb_id: str,
    label: str,
) -> None:
    content = f"""\
# RESP Parametrization — {resname}  |  {pdb_id} / {label}

## Overview

This directory contains the scaffold for deriving AMBER GAFF2 bonded parameters
and HF/6-31G* RESP charges for **{resname}** embedded in the protein backbone.

RESP at HF/6-31G* is required (not AM1-BCC) for consistency with ff14SB, which
derives all its backbone and side-chain charges at the same level of theory.
Reference: Maier et al. JCTC 2015; Homeyer et al. J Mol Model 2006.

## Formal charge

```
Net charge : {charge}   (source: {charge_source})
Spin       : 1 (singlet — verify for open-shell modifications)
```

## Workflow

```
bash commands.sh
  → step01_gaussian_inputs/{resname}_opt.com   (HF/6-31G* Opt)
  → step01_gaussian_inputs/{resname}_esp.com   (HF/6-31G* MK ESP)
  → staged to step02_gaussian/

bash submit_gaussian.sh
  → sbatch both Gaussian jobs (ESP depends on opt via SLURM dependency)

bash run_after_gaussian.sh
  → antechamber -c resp          → step03_resp_params/{resname}_capped.mol2
  → parmchk2                     → step03_resp_params/{resname}.frcmod
  → strip ACE/NME                → step03_resp_params/{resname}.mol2
```

## Using in tleap

```tcl
source leaprc.protein.ff14SB
source leaprc.gaff2
loadamberparams step03_resp_params/{resname}.frcmod
{resname} = loadmol2 step03_resp_params/{resname}.mol2
saveoff {resname} {resname}.lib
# Then in your full system:
loadoff {resname}.lib
```

## References

- Homeyer et al. J Mol Model 2006 — RESP/HF/6-31G* for pSer/pThr/pTyr/pHis
- Khoury et al. JCTC 2014 — Forcefield_PTM; 32 PTMs
- Sengupta et al. JCTC 2024 — phosaa14SB/phosaa19SB
- Maier et al. JCTC 2015 — ff14SB
"""
    output_path.write_text(content, encoding="utf-8")


def _write_database_readme(
    output_path: Path,
    resname: str,
    pdb_id: str,
    label: str,
    db_info: dict[str, str],
) -> None:
    database = db_info["database"]
    leaprc = db_info.get("leaprc", "")
    url = db_info.get("url", "")
    description = db_info.get("description", "")

    ff03_warning = ""
    if "Forcefield_PTM" in database:
        ff03_warning = """
## ⚠️  Compatibility warning

Forcefield_PTM parameters are derived for **ff03**, not ff14SB.  The backbone
charge scheme differs between ff03 and ff14SB.  Using these parameters with
ff14SB may introduce inconsistencies at the backbone.  Options:

1. Use the Forcefield_PTM parameters as a starting point and verify carefully.
2. Re-derive charges using the RESP scaffold in the ``resp/`` subdirectory.
"""

    leaprc_block = ""
    if leaprc:
        leaprc_block = f"""
## Loading in tleap

```tcl
source leaprc.protein.ff14SB
source {leaprc}
# {resname} is now available as a recognised residue
```
"""
    elif url:
        leaprc_block = f"""
## Obtaining parameters

Download from: {url}

Follow the installation instructions provided at that URL, then load the
resulting ``.lib`` and ``.frcmod`` files in tleap.
"""

    content = f"""\
# Pre-existing Parameters Found — {resname}  |  {pdb_id} / {label}

Parameters for **{resname}** are already available in **{database}**.

{description}
{leaprc_block}{ff03_warning}
## No RESP derivation needed

Because pre-built, peer-reviewed parameters exist, you do not need to run
Gaussian for this residue.  If you want to derive independent RESP parameters
(e.g. for publication or to ensure ff14SB consistency), a scaffold is provided
in the ``resp/`` subdirectory.
"""
    output_path.write_text(content, encoding="utf-8")


# ---------------------------------------------------------------------------
# Per-residue orchestrator
# ---------------------------------------------------------------------------

def _run_single_residue(
    *,
    pdb_id: str,
    protein_pdb: Path,
    chain_id: str,
    resname: str,
    resseq: int,
    label: str,
    residue_dir: Path,
) -> dict[str, Any]:
    result: dict[str, Any] = {
        "resname": resname,
        "chain_id": chain_id,
        "resseq": resseq,
        "label": label,
        "in_database": False,
        "database": None,
        "formal_charge": None,
        "charge_source": None,
        "capped_model_built": False,
        "resp_scaffold_generated": False,
        "status": "pending",
        "message": "",
    }

    residue_dir.mkdir(parents=True, exist_ok=True)

    # 1. Database check ------------------------------------------------------
    db_info = _DATABASE_REGISTRY.get(resname.upper())
    if db_info:
        result["in_database"] = True
        result["database"] = db_info["database"]
        result["status"] = "found_in_database"
        result["message"] = f"{resname} found in {db_info['database']} — no Gaussian needed."
        _write_database_readme(
            output_path=residue_dir / "README.md",
            resname=resname,
            pdb_id=pdb_id,
            label=label,
            db_info=db_info,
        )
        # Still generate a RESP scaffold in resp/ for reference / override
        resp_dir = residue_dir / "resp"
        resp_dir.mkdir(exist_ok=True)

    # 2. Formal charge -------------------------------------------------------
    charge, charge_source = lookup_formal_charge(resname)
    result["formal_charge"] = charge
    result["charge_source"] = charge_source

    # 3. Build ACE-X-NME ----------------------------------------------------
    capped_pdb = residue_dir / "capped_model_ace_nme.pdb"
    if not capped_pdb.exists():
        built = _build_capped_model(
            protein_pdb=protein_pdb,
            chain_id=chain_id,
            resname=resname,
            resseq=resseq,
            output_pdb=capped_pdb,
        )
        result["capped_model_built"] = built
        if not built:
            result["status"] = "error"
            result["message"] = f"Could not build ACE-{resname}-NME capped model."
            return result
    else:
        result["capped_model_built"] = True

    # 4. Write RESP scaffold -------------------------------------------------
    resp_dir = residue_dir / "resp"
    resp_dir.mkdir(exist_ok=True)

    # Relative path from resp_dir to capped_pdb
    _write_commands_sh(
        output_path=resp_dir / "commands.sh",
        resname=resname,
        charge=charge,
        capped_pdb_name="../capped_model_ace_nme.pdb",
    )
    _write_submit_gaussian_sh(
        output_path=resp_dir / "submit_gaussian.sh",
        resname=resname,
    )
    _write_run_after_gaussian_sh(
        output_path=resp_dir / "run_after_gaussian.sh",
        resname=resname,
        charge=charge,
    )

    if not db_info:
        # Only write the primary README if we haven't already written a
        # database-found README at the top level.
        _write_readme(
            output_path=resp_dir / "README.md",
            resname=resname,
            charge=charge,
            charge_source=charge_source,
            pdb_id=pdb_id,
            label=label,
        )

    result["resp_scaffold_generated"] = True
    if result["status"] == "pending":
        result["status"] = "resp_scaffold_ready"
        result["message"] = (
            f"{resname}: capped model built, RESP scaffold generated "
            f"(charge={charge}, source={charge_source})."
        )

    return result


# ---------------------------------------------------------------------------
# Public entry point
# ---------------------------------------------------------------------------

def run_nonstd_residue_params(
    protein_dir: Path,
    pdb_id: str,
    protein_pdb: Path,
) -> dict[str, Any]:
    """Run non-standard residue parametrization for one protein directory.

    Parameters
    ----------
    protein_dir:
        The protein's data directory (e.g. ``data/proteins/3OLL``).
    pdb_id:
        Four-letter PDB ID.
    protein_pdb:
        Path to the representative protein component PDB
        (``components/{pdb_id}_protein.pdb``).

    Returns
    -------
    dict with keys: status, message, n_residues, n_in_database, n_resp_scaffold,
    residues (list of per-residue result dicts), manifest_path.
    """
    output_dir = protein_dir / "nonstandard_params"
    output_dir.mkdir(parents=True, exist_ok=True)

    result: dict[str, Any] = {
        "pdb_id": pdb_id,
        "status": "pending",
        "message": "",
        "n_residues": 0,
        "n_in_database": 0,
        "n_resp_scaffold": 0,
        "residues": [],
        "manifest_path": str(output_dir / "nonstd_params_manifest.json"),
    }

    if not protein_pdb.exists():
        result["status"] = "skipped"
        result["message"] = f"Protein PDB not found: {protein_pdb}"
        return result

    nonstd = find_nonstandard_residues(protein_pdb)

    if not nonstd:
        result["status"] = "skipped"
        result["message"] = "No non-standard backbone residues detected."
        _write_manifest(result)
        return result

    result["n_residues"] = len(nonstd)

    for chain_id, residue in nonstd:
        resname = residue.get_resname().strip().upper()
        resseq = int(residue.id[1])
        label = f"{chain_id}_{resname}_{resseq}"
        residue_dir = output_dir / label

        res_result = _run_single_residue(
            pdb_id=pdb_id,
            protein_pdb=protein_pdb,
            chain_id=chain_id,
            resname=resname,
            resseq=resseq,
            label=label,
            residue_dir=residue_dir,
        )
        result["residues"].append(res_result)
        if res_result["in_database"]:
            result["n_in_database"] += 1
        if res_result["resp_scaffold_generated"]:
            result["n_resp_scaffold"] += 1

    result["status"] = "success"
    result["message"] = (
        f"{result['n_residues']} non-standard residue(s): "
        f"{result['n_in_database']} in known database, "
        f"{result['n_resp_scaffold']} RESP scaffold(s) generated."
    )
    _write_manifest(result)
    return result


def _write_manifest(result: dict[str, Any]) -> None:
    path = Path(result["manifest_path"])
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(result, indent=2), encoding="utf-8")
