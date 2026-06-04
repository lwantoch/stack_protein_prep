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

from pathlib import Path
from typing import Any

from stack_protein_preparation.cap import find_nonstandard_residues

# Re-export everything from _nonstd_residue_params_core so existing imports keep working
from stack_protein_preparation._nonstd_residue_params_core import (  # noqa: F401 (re-exported)
    _DATABASE_REGISTRY,
    _PHYSIOLOGICAL_CHARGES,
    _format_atom_line,
    _build_capped_model,
    _lookup_charge_rcsb,
    lookup_formal_charge,
    _read_pdb_atoms,
    _write_gaussian_com,
    _write_commands_sh,
    _write_submit_gaussian_sh,
    _write_run_after_gaussian_sh,
    _write_readme,
    _write_database_readme,
    _run_single_residue,
    _write_manifest,
)


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
