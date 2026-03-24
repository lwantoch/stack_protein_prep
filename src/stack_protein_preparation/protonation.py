# /home/grheco/repositorios/stack_protein_prep/src/stack_protein_preparation/protonation.py

"""
/home/grheco/repositorios/stack_protein_prep/src/stack_protein_preparation/protonation.py

Select a prepared protein structure, protonate it with pdb2pqr/PROPKA,
and write a standardized protonated output PDB.

Responsibilities
----------------
- choose the best available structure for protonation
- prefer MODELLER output over AlphaFold fallback over original protein input
- run pdb2pqr with PROPKA-based protonation at a chosen pH
- write the protonated structure to a fixed components output path
- return a result dictionary for pipeline state tracking

Expected directory layout
-------------------------
Input:
    data/proteins/<PDB_ID>/
    ├── components/
    │   └── <PDB_ID>_protein.pdb
    ├── MODELLER/
    │   └── ...
    └── fasta/
        └── alignments/
            └── filler/
                └── ...

Possible selected input structures:
    - MODELLER final model
    - AlphaFold fallback model
    - components/<PDB_ID>_protein.pdb

Output:
    data/proteins/<PDB_ID>/components/
    └── <PDB_ID>_proteinH.pdb

Selection logic
---------------
Input priority is strictly:
1. modeller model
2. alphafold model
3. original components protein

A candidate input file is accepted only if it exists and is non-empty.

Important
---------
- the output filename is always standardized to <PDB_ID>_proteinH.pdb
- pdb2pqr output may still reflect PQR-like formatting semantics
- this module does not perform AMBER residue renaming
- atom counts are computed from ATOM/HETATM coordinate records
- this module is intended to run after filler/model-selection and before
  AMBER-specific cleanup or topology preparation
"""

from __future__ import annotations

import shutil
import subprocess
from pathlib import Path
from typing import Literal

ProtonationInputSource = Literal["protein", "modeller", "alphafold"]


def _find_pdb2pqr_executable() -> str:
    """
    Return the first available pdb2pqr executable.

    Common names differ by installation/channel.
    """
    candidates = ["pdb2pqr", "pdb2pqr30"]

    for candidate in candidates:
        if shutil.which(candidate):
            return candidate

    raise FileNotFoundError(
        "Could not find a pdb2pqr executable in PATH. "
        "Expected one of: pdb2pqr, pdb2pqr30"
    )


def count_atoms_in_structure_file(pdb_path: str | Path) -> int:
    """
    Count coordinate records in a structure-like text file.

    Notes
    -----
    This intentionally counts ATOM/HETATM lines directly instead of using a full
    structure parser, because pdb2pqr output may be PQR-like while still stored
    with a .pdb filename.
    """
    pdb_path = Path(pdb_path)

    if not pdb_path.is_file():
        raise FileNotFoundError(f"Structure file not found: {pdb_path}")

    atom_count = 0

    with pdb_path.open("r", encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("ATOM") or line.startswith("HETATM"):
                atom_count += 1

    return atom_count


def select_protonation_input(
    pdb_id: str,
    protein_dir: str | Path,
    modeller_model_path: str | Path | None = None,
    alphafold_model_path: str | Path | None = None,
) -> tuple[Path, ProtonationInputSource]:
    """
    Select the input structure for protonation.

    Preference order:
    1. modeller model, if provided and existing
    2. alphafold model, if provided and existing
    3. components/<PDBID>_protein.pdb

    Returns
    -------
    tuple[Path, ProtonationInputSource]
        (selected_input_path, source_label)
    """
    protein_dir = Path(protein_dir)
    components_dir = protein_dir / "components"
    default_protein_path = components_dir / f"{pdb_id}_protein.pdb"

    if modeller_model_path is not None:
        modeller_model_path = Path(modeller_model_path)
        if modeller_model_path.is_file() and modeller_model_path.stat().st_size > 0:
            return modeller_model_path, "modeller"

    if alphafold_model_path is not None:
        alphafold_model_path = Path(alphafold_model_path)
        if alphafold_model_path.is_file() and alphafold_model_path.stat().st_size > 0:
            return alphafold_model_path, "alphafold"

    if default_protein_path.is_file() and default_protein_path.stat().st_size > 0:
        return default_protein_path, "protein"

    raise FileNotFoundError(
        f"No valid protonation input found for {pdb_id}. "
        f"Checked modeller_model_path={modeller_model_path!s}, "
        f"alphafold_model_path={alphafold_model_path!s}, and "
        f"default protein path={default_protein_path!s}."
    )


def run_pdb2pqr_protonation(
    input_pdb: str | Path,
    output_pdb: str | Path,
    ph: float = 7.4,
    ff: str = "AMBER",
    keep_chain: bool = True,
    keep_heterogens: bool = False,
) -> subprocess.CompletedProcess[str]:
    """
    Protonate a protein structure with pdb2pqr/PROPKA and write a PDB output.

    Notes
    -----
    - pdb2pqr writes the output structure to the last positional path
    - depending on the installed version, output formatting may be closer to PQR
      semantics even if the file extension is .pdb
    """
    input_pdb = Path(input_pdb)
    output_pdb = Path(output_pdb)
    output_pdb.parent.mkdir(parents=True, exist_ok=True)

    executable = _find_pdb2pqr_executable()

    cmd = [
        executable,
        f"--ff={ff}",
        "--titration-state-method=propka",
        f"--with-ph={ph}",
    ]

    if keep_chain:
        cmd.append("--keep-chain")

    if keep_heterogens:
        cmd.append("--keep-heterogens")

    cmd.extend([str(input_pdb), str(output_pdb)])

    result = subprocess.run(
        cmd,
        check=True,
        capture_output=True,
        text=True,
    )
    return result


def protonate_protein_structure(
    pdb_id: str,
    protein_dir: str | Path,
    modeller_model_path: str | Path | None = None,
    alphafold_model_path: str | Path | None = None,
    ph: float = 7.4,
) -> dict[str, str | bool | float | int]:
    """
    Select protonation input, run pdb2pqr/propka, and save a standard output file.

    Output path is always:
        components/<PDBID>_proteinH.pdb

    Selection priority:
        modeller -> alphafold -> original components protein
    """
    protein_dir = Path(protein_dir)
    components_dir = protein_dir / "components"
    output_pdb = components_dir / f"{pdb_id}_proteinH.pdb"

    input_pdb, input_source = select_protonation_input(
        pdb_id=pdb_id,
        protein_dir=protein_dir,
        modeller_model_path=modeller_model_path,
        alphafold_model_path=alphafold_model_path,
    )

    input_atom_count = count_atoms_in_structure_file(input_pdb)

    result = run_pdb2pqr_protonation(
        input_pdb=input_pdb,
        output_pdb=output_pdb,
        ph=ph,
        ff="AMBER",
        keep_chain=True,
        keep_heterogens=False,
    )

    output_exists = output_pdb.is_file()
    output_nonempty = output_exists and output_pdb.stat().st_size > 0

    if output_nonempty:
        output_atom_count = count_atoms_in_structure_file(output_pdb)
        atom_count_increased = output_atom_count > input_atom_count
    else:
        output_atom_count = 0
        atom_count_increased = False

    protonation_success = output_nonempty

    return {
        "protonation_success": protonation_success,
        "protonation_input_source": input_source,
        "protonation_input_path": str(input_pdb),
        "protonation_output_path": str(output_pdb),
        "protonation_ph": ph,
        "input_atom_count": input_atom_count,
        "output_atom_count": output_atom_count,
        "atom_count_increased": atom_count_increased,
        "pdb2pqr_stdout": result.stdout,
        "pdb2pqr_stderr": result.stderr,
    }
