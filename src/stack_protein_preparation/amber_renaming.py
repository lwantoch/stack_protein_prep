# /home/grheco/repositorios/stack_protein_prep/src/stack_protein_preparation/amber_renaming.py

"""
/home/grheco/repositorios/stack_protein_prep/src/stack_protein_preparation/amber_renaming.py

Rename a protonated protein PDB to AMBER residue naming convention.

Responsibilities
----------------
- read a protonated protein structure from components/<PDB_ID>_proteinH.pdb
- rename selected standard residues in-place based on hydrogen presence
- detect disulfide-bonded cysteines by SG-SG distance
- write a standardized AMBER-renamed output PDB
- write renaming statistics to JSON
- optionally append a plain-text log entry

Expected directory layout
-------------------------
Input:
    data/proteins/<PDB_ID>/components/
    └── <PDB_ID>_proteinH.pdb

Output:
    data/proteins/<PDB_ID>/components/
    ├── <PDB_ID>_protein_as_Amber.pdb
    ├── <PDB_ID>_protein_as_Amber_stats.json
    └── amber_renamed.log

Renaming rules
--------------
- HIS -> HID / HIE / HIP based on hydrogens
- ASP -> ASH if protonated
- GLU -> GLH if protonated
- CYS -> CYM if deprotonated
- CYS/CYM in disulfide -> CYX

Important
---------
- no renumbering is performed
- waters are ignored
- hetero residues are ignored
- only primary altloc (' ' or 'A') atoms are considered
- this module expects a protonated input structure
- returned statistics are intended for pipeline state / XLSX summaries
"""

from __future__ import annotations

import datetime as _dt
import json
import math
from pathlib import Path
from typing import Dict, Iterable, Set, Tuple, Union

from Bio.PDB import PDBIO, PDBParser
from Bio.PDB.Atom import Atom
from Bio.PDB.Residue import Residue

ALTLOC_OK: Set[str] = {" ", "A"}  # keep primary or 'A' conformer


def _iter_atoms_altloc_ok(residue: Residue) -> Iterable[Atom]:
    """
    Yield atoms from one residue, keeping only accepted altloc states.
    """
    for atom in residue.get_atoms():
        if atom.get_altloc() in ALTLOC_OK:
            yield atom


def _atom_names(residue: Residue) -> Set[str]:
    """
    Return a set of accepted atom names for one residue.
    """
    return {atom.get_name().strip() for atom in _iter_atoms_altloc_ok(residue)}


def _distance(a, b) -> float:
    """
    Return Euclidean distance between two 3D coordinates.
    """
    dx = float(a[0] - b[0])
    dy = float(a[1] - b[1])
    dz = float(a[2] - b[2])
    return math.sqrt(dx * dx + dy * dy + dz * dz)


def rename_structure_by_hydrogens(
    structure,
    *,
    disulf_min: float = 1.8,
    disulf_max: float = 2.2,
    strict_his: bool = False,
) -> Dict[str, int]:
    """
    Rename residues in-place for AMBER-style naming.

    Rules
    -----
    - HIS -> HID / HIE / HIP based on hydrogens
    - ASP -> ASH if protonated
    - GLU -> GLH if protonated
    - CYS -> CYM if deprotonated
    - CYS/CYM in disulfide -> CYX

    Notes
    -----
    - no renumbering is performed
    - waters are ignored
    - hetero residues (ligands, ions, cofactors) are ignored
    """
    stats: Dict[str, int] = {}

    sulfur_residue_list: list[Tuple[Residue, Tuple[float, float, float]]] = []

    for model in structure:
        for chain in model:
            for residue in chain:
                resname = residue.get_resname().strip()

                if resname not in {"CYS", "CYM"}:
                    continue

                atoms = {
                    atom.get_name().strip(): atom
                    for atom in _iter_atoms_altloc_ok(residue)
                }
                sg_atom = atoms.get("SG")

                if sg_atom is None:
                    continue

                sulfur_residue_list.append((residue, tuple(sg_atom.get_coord())))

    disulfide_residue_id_set: Set[int] = set()

    for i in range(len(sulfur_residue_list)):
        residue_i, coord_i = sulfur_residue_list[i]

        for j in range(i + 1, len(sulfur_residue_list)):
            residue_j, coord_j = sulfur_residue_list[j]
            sulfur_distance = _distance(coord_i, coord_j)

            if disulf_min <= sulfur_distance <= disulf_max:
                disulfide_residue_id_set.add(id(residue_i))
                disulfide_residue_id_set.add(id(residue_j))

    for model in structure:
        for chain in model:
            for residue in chain:
                resname = residue.get_resname().strip()

                if resname in {"HOH", "WAT"}:
                    continue

                if str(residue.id[0]).startswith("H_"):
                    continue

                atom_name_set = _atom_names(residue)

                if resname in {"CYS", "CYM"}:
                    if id(residue) in disulfide_residue_id_set:
                        if resname != "CYX":
                            residue.resname = "CYX"
                            stats[f"{resname}->CYX"] = (
                                stats.get(f"{resname}->CYX", 0) + 1
                            )
                        else:
                            stats["CYX kept"] = stats.get("CYX kept", 0) + 1
                        continue

                    has_hg = ("HG" in atom_name_set) or any(
                        atom_name.endswith("HG") for atom_name in atom_name_set
                    )

                    if has_hg:
                        if resname != "CYS":
                            residue.resname = "CYS"
                            stats[f"{resname}->CYS"] = (
                                stats.get(f"{resname}->CYS", 0) + 1
                            )
                        else:
                            stats["CYS kept"] = stats.get("CYS kept", 0) + 1
                    else:
                        if resname != "CYM":
                            residue.resname = "CYM"
                            stats[f"{resname}->CYM"] = (
                                stats.get(f"{resname}->CYM", 0) + 1
                            )
                        else:
                            stats["CYM kept"] = stats.get("CYM kept", 0) + 1
                    continue

                if resname == "HIS":
                    has_hd1 = ("HD1" in atom_name_set) or any(
                        atom_name.endswith("HD1") for atom_name in atom_name_set
                    )
                    has_he2 = ("HE2" in atom_name_set) or any(
                        atom_name.endswith("HE2") for atom_name in atom_name_set
                    )

                    if has_hd1 and has_he2:
                        residue.resname = "HIP"
                        stats["HIS->HIP"] = stats.get("HIS->HIP", 0) + 1
                    elif has_hd1:
                        residue.resname = "HID"
                        stats["HIS->HID"] = stats.get("HIS->HID", 0) + 1
                    elif has_he2:
                        residue.resname = "HIE"
                        stats["HIS->HIE"] = stats.get("HIS->HIE", 0) + 1
                    else:
                        if strict_his:
                            stats["HIS_no_protons"] = stats.get("HIS_no_protons", 0) + 1
                    continue

                if resname == "ASP":
                    has_hd2 = ("HD2" in atom_name_set) or any(
                        atom_name.endswith("HD2") for atom_name in atom_name_set
                    )
                    if has_hd2:
                        residue.resname = "ASH"
                        stats["ASP->ASH"] = stats.get("ASP->ASH", 0) + 1
                    continue

                if resname == "GLU":
                    has_he2 = ("HE2" in atom_name_set) or any(
                        atom_name.endswith("HE2") for atom_name in atom_name_set
                    )
                    if has_he2:
                        residue.resname = "GLH"
                        stats["GLU->GLH"] = stats.get("GLU->GLH", 0) + 1
                    continue

    return stats


def _append_renamed_log(
    log_path: Path,
    *,
    input_pdb_path: Path,
    output_pdb_path: Path,
    stats: Dict[str, int],
    disulf_min: float,
    disulf_max: float,
    strict_his: bool,
) -> None:
    """
    Append one plain-text AMBER renaming log line.
    """
    timestamp = _dt.datetime.now().isoformat(timespec="seconds")
    parts = [f"{key}={stats[key]}" for key in sorted(stats)] if stats else ["(none)"]

    line = (
        f"[{timestamp}] input={input_pdb_path} output={output_pdb_path} "
        f"disulf_min={disulf_min} disulf_max={disulf_max} strict_his={strict_his} "
        f"stats={' '.join(parts)}\n"
    )

    log_path.parent.mkdir(parents=True, exist_ok=True)
    with log_path.open("a", encoding="utf-8") as handle:
        handle.write(line)


def write_amber_renamed_pdb(
    input_pdb_path: Union[str, Path],
    output_pdb_path: Union[str, Path],
    *,
    log_path: Union[str, Path, None] = None,
    disulf_min: float = 1.8,
    disulf_max: float = 2.2,
    strict_his: bool = False,
) -> Dict[str, int]:
    """
    Read a protonated PDB and write a new AMBER-renamed PDB.

    Typical usage
    -------------
    input : <PDBID>_proteinH.pdb
    output: <PDBID>_protein_as_Amber.pdb

    Also writes
    -----------
    <PDBID>_protein_as_Amber_stats.json
    """
    input_pdb_path = Path(input_pdb_path).resolve()
    output_pdb_path = Path(output_pdb_path).resolve()

    if not input_pdb_path.exists():
        raise FileNotFoundError(input_pdb_path)

    if log_path is None:
        log_path = output_pdb_path.parent / "amber_renamed.log"
    else:
        log_path = Path(log_path).resolve()

    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("prot", str(input_pdb_path))

    stats = rename_structure_by_hydrogens(
        structure,
        disulf_min=disulf_min,
        disulf_max=disulf_max,
        strict_his=strict_his,
    )

    output_pdb_path.parent.mkdir(parents=True, exist_ok=True)

    io = PDBIO()
    io.set_structure(structure)
    io.save(str(output_pdb_path))

    _append_renamed_log(
        log_path=log_path,
        input_pdb_path=input_pdb_path,
        output_pdb_path=output_pdb_path,
        stats=stats,
        disulf_min=disulf_min,
        disulf_max=disulf_max,
        strict_his=strict_his,
    )

    stats_path = output_pdb_path.with_name(f"{output_pdb_path.stem}_stats.json")
    stats_path.write_text(
        json.dumps(stats, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )

    return stats


def amber_rename_protein_structure(
    pdb_id: str,
    protein_dir: str | Path,
    *,
    strict_his: bool = False,
    disulf_min: float = 1.8,
    disulf_max: float = 2.2,
) -> dict[str, str | bool | int | float]:
    """
    Rename a protonated protein structure to AMBER convention.

    Expected input path:
        components/<PDBID>_proteinH.pdb

    Standardized output path:
        components/<PDBID>_protein_as_Amber.pdb
    """
    protein_dir = Path(protein_dir)
    components_dir = protein_dir / "components"

    input_pdb_path = components_dir / f"{pdb_id}_proteinH.pdb"
    output_pdb_path = components_dir / f"{pdb_id}_protein_as_Amber.pdb"
    log_path = components_dir / "amber_renamed.log"

    stats = write_amber_renamed_pdb(
        input_pdb_path=input_pdb_path,
        output_pdb_path=output_pdb_path,
        log_path=log_path,
        disulf_min=disulf_min,
        disulf_max=disulf_max,
        strict_his=strict_his,
    )

    his_to_hid = int(stats.get("HIS->HID", 0))
    his_to_hie = int(stats.get("HIS->HIE", 0))
    his_to_hip = int(stats.get("HIS->HIP", 0))

    cys_to_cym = int(stats.get("CYS->CYM", 0))
    cym_to_cys = int(stats.get("CYM->CYS", 0))
    cys_to_cyx = int(stats.get("CYS->CYX", 0))
    cym_to_cyx = int(stats.get("CYM->CYX", 0))

    asp_to_ash = int(stats.get("ASP->ASH", 0))
    glu_to_glh = int(stats.get("GLU->GLH", 0))

    amber_renaming_success = (
        output_pdb_path.is_file() and output_pdb_path.stat().st_size > 0
    )

    return {
        "amber_renaming_success": amber_renaming_success,
        "amber_input_path": str(input_pdb_path),
        "amber_output_path": str(output_pdb_path),
        "amber_log_path": str(log_path),
        "amber_strict_his": strict_his,
        "amber_disulf_min": disulf_min,
        "amber_disulf_max": disulf_max,
        "his_to_hid": his_to_hid,
        "his_to_hie": his_to_hie,
        "his_to_hip": his_to_hip,
        "cys_to_cym": cys_to_cym,
        "cym_to_cys": cym_to_cys,
        "cys_to_cyx": cys_to_cyx,
        "cym_to_cyx": cym_to_cyx,
        "asp_to_ash": asp_to_ash,
        "glu_to_glh": glu_to_glh,
        "amber_stats_json_path": str(
            output_pdb_path.with_name(f"{output_pdb_path.stem}_stats.json")
        ),
    }
