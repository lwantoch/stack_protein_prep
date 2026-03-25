"""
/home/grheco/repositorios/stack_protein_prep/src/stack_protein_preparation/terminus.py

Convert true chain termini of a protein-only PDB into AMBER-compatible termini.

Purpose
-------
- detect the first and last polymer residue in each chain
- convert the first residue to an AMBER N-terminus residue name
  (for example ALA -> NALA)
- convert the last residue to an AMBER C-terminus residue name
  (for example ALA -> CALA)
- normalize terminal atom naming for AMBER / tleap compatibility
- write the result to the standard pipeline output path

Standard output naming
----------------------
Given:
    data/proteins/<PDBID>/

this module writes by default:
    data/proteins/<PDBID>/components/<PDBID>_protein_amber_termini.pdb

Important
---------
- intended input: protein-only PDB
- recommended input in your pipeline: the protonated protein PDB
- this module is for TRUE chain termini
- internal cuts / missing internal fragments should be handled separately in cap.py
- existing AMBER-style internal residue names are preserved:
    HID, HIE, HIP, ASH, GLH, CYM, CYX, LYN, ...
- newly added atom coordinates are heuristic and should be relaxed later
- single-residue chains are chemically awkward for generic AMBER terminal naming;
  this module preserves the base residue name in that special case, but still
  normalizes the terminal atoms

Dependencies
------------
- Biopython
- numpy
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import numpy as np
from Bio.PDB import PDBIO, Atom, Chain, Model, PDBParser, Residue, Structure

SUPPORTED_POLYMER_RESIDUE_NAMES = {
    "ALA",
    "ARG",
    "ASN",
    "ASP",
    "ASH",
    "CYS",
    "CYM",
    "CYX",
    "GLN",
    "GLU",
    "GLH",
    "GLY",
    "HIS",
    "HID",
    "HIE",
    "HIP",
    "ILE",
    "LEU",
    "LYS",
    "LYN",
    "MET",
    "MSE",
    "PHE",
    "PRO",
    "SER",
    "THR",
    "TRP",
    "TYR",
    "VAL",
}

N_TERMINAL_RESNAME_MAP = {
    residue_name: f"N{residue_name}" for residue_name in SUPPORTED_POLYMER_RESIDUE_NAMES
}

C_TERMINAL_RESNAME_MAP = {
    residue_name: f"C{residue_name}" for residue_name in SUPPORTED_POLYMER_RESIDUE_NAMES
}

NTERM_HYDROGEN_ALIASES = {
    "H",
    "H1",
    "H2",
    "H3",
    "HT1",
    "HT2",
    "HT3",
    "1H",
    "2H",
    "3H",
}

CTERM_OXYGEN_ALIASES = {
    "O",
    "OXT",
    "OT1",
    "OT2",
    "O1",
    "O2",
}


@dataclass(slots=True)
class ChainTerminusResult:
    chain_id: str
    first_resseq: int | None
    last_resseq: int | None
    first_old_resname: str | None
    first_new_resname: str | None
    last_old_resname: str | None
    last_new_resname: str | None
    n_h_added: int
    oxt_added: bool
    single_residue_chain: bool


def _as_vec(values: Iterable[float]) -> np.ndarray:
    return np.asarray(list(values), dtype=float)


def _norm(vector: np.ndarray) -> float:
    return float(np.linalg.norm(vector))


def _unit(vector: np.ndarray, fallback: np.ndarray | None = None) -> np.ndarray:
    length = _norm(vector)
    if length > 1e-8:
        return vector / length

    if fallback is not None and _norm(fallback) > 1e-8:
        return fallback / _norm(fallback)

    return np.array([1.0, 0.0, 0.0], dtype=float)


def _orthogonal_unit(vector: np.ndarray) -> np.ndarray:
    vector_u = _unit(vector)
    trial = np.array([1.0, 0.0, 0.0], dtype=float)
    if abs(float(np.dot(vector_u, trial))) > 0.9:
        trial = np.array([0.0, 1.0, 0.0], dtype=float)
    return _unit(np.cross(vector_u, trial), fallback=np.array([0.0, 0.0, 1.0]))


def _safe_occupancy(atom: Atom.Atom) -> float:
    return 1.0 if atom.get_occupancy() is None else float(atom.get_occupancy())


def _safe_element(atom: Atom.Atom) -> str:
    if atom.element and atom.element.strip():
        return atom.element.strip()
    atom_name = atom.get_name().strip()
    return atom_name[0] if atom_name else "X"


def _make_atom(
    name: str,
    coord: np.ndarray,
    serial_number: int,
    element: str,
    bfactor: float = 20.0,
    occupancy: float = 1.0,
) -> Atom.Atom:
    return Atom.Atom(
        name=name,
        coord=np.asarray(coord, dtype=float),
        bfactor=bfactor,
        occupancy=occupancy,
        altloc=" ",
        fullname=f"{name:>4}",
        serial_number=serial_number,
        element=element,
    )


def _clone_atom_with_new_name(
    atom: Atom.Atom,
    serial_number: int,
    new_name: str | None = None,
) -> Atom.Atom:
    output_name = atom.get_name() if new_name is None else new_name
    return Atom.Atom(
        name=output_name,
        coord=np.array(atom.get_coord(), dtype=float),
        bfactor=float(atom.get_bfactor()),
        occupancy=_safe_occupancy(atom),
        altloc=atom.get_altloc(),
        fullname=f"{output_name:>4}",
        serial_number=serial_number,
        element=_safe_element(atom) if new_name is None else output_name[0],
    )


def _copy_residue_with_skips(
    residue: Residue.Residue,
    serial_counter: list[int],
    new_resname: str | None = None,
    skip_atom_names: set[str] | None = None,
) -> Residue.Residue:
    skip_atom_names = skip_atom_names or set()

    output_residue = Residue.Residue(
        id=residue.id,
        resname=residue.get_resname() if new_resname is None else new_resname,
        segid=residue.get_segid(),
    )

    for atom in residue.get_atoms():
        if atom.get_name().strip() in skip_atom_names:
            continue
        serial_counter[0] += 1
        output_residue.add(
            _clone_atom_with_new_name(
                atom=atom,
                serial_number=serial_counter[0],
            )
        )

    return output_residue


def _is_polymer_residue(residue: Residue.Residue) -> bool:
    if residue.id[0] != " ":
        return False
    return residue.get_resname().strip().upper() in SUPPORTED_POLYMER_RESIDUE_NAMES


def _require_atom(residue: Residue.Residue, atom_name: str) -> Atom.Atom:
    if atom_name not in residue:
        raise ValueError(
            f"Residue {residue.get_resname()} {residue.id} is missing required atom '{atom_name}'."
        )
    return residue[atom_name]


def _find_first_and_last_polymer_indices(
    chain: Chain.Chain,
) -> tuple[int, int] | None:
    residue_list = list(chain.get_residues())

    polymer_indices = [
        index
        for index, residue in enumerate(residue_list)
        if _is_polymer_residue(residue)
    ]

    if not polymer_indices:
        return None

    return polymer_indices[0], polymer_indices[-1]


def _collect_existing_nterm_hydrogens(residue: Residue.Residue) -> list[Atom.Atom]:
    hydrogen_atoms: list[Atom.Atom] = []

    for atom in residue.get_atoms():
        if atom.get_name().strip() in NTERM_HYDROGEN_ALIASES:
            hydrogen_atoms.append(atom)

    return hydrogen_atoms


def _target_nterm_h_names(base_resname: str) -> list[str]:
    if base_resname.upper() == "PRO":
        return ["H2", "H3"]
    return ["H1", "H2", "H3"]


def _build_missing_nterm_h_positions(
    residue: Residue.Residue,
    missing_count: int,
) -> list[np.ndarray]:
    if missing_count <= 0:
        return []

    n_atom = _require_atom(residue, "N")
    ca_atom = _require_atom(residue, "CA")

    n_coord = _as_vec(n_atom.coord)
    ca_coord = _as_vec(ca_atom.coord)

    bond_length = 1.01

    if residue.get_resname().strip().upper() == "PRO" and "CD" in residue:
        cd_coord = _as_vec(residue["CD"].coord)
        axis = -_unit((ca_coord - n_coord) + (cd_coord - n_coord))
    else:
        axis = -_unit(ca_coord - n_coord)

    u = _orthogonal_unit(axis)
    v = _unit(np.cross(axis, u))
    theta = np.deg2rad(55.0)

    if missing_count == 1:
        direction_list = [axis]
    elif missing_count == 2:
        direction_list = [
            _unit(np.cos(theta) * axis + np.sin(theta) * u),
            _unit(np.cos(theta) * axis - np.sin(theta) * u),
        ]
    else:
        direction_list = [
            _unit(np.cos(theta) * axis + np.sin(theta) * u),
            _unit(np.cos(theta) * axis + np.sin(theta) * (-0.5 * u + 0.8660254 * v)),
            _unit(np.cos(theta) * axis + np.sin(theta) * (-0.5 * u - 0.8660254 * v)),
        ]

    return [n_coord + bond_length * direction for direction in direction_list]


def _normalize_or_add_nterm_hydrogens(
    source_residue: Residue.Residue,
    output_residue: Residue.Residue,
    serial_counter: list[int],
) -> int:
    """
    Normalize N-terminal hydrogens to AMBER-style names.

    PRO -> H2, H3
    all others -> H1, H2, H3
    """
    target_names = _target_nterm_h_names(source_residue.get_resname().strip().upper())
    target_count = len(target_names)

    existing_h_atoms = _collect_existing_nterm_hydrogens(source_residue)[:target_count]

    names_to_remove = [
        atom.get_id()
        for atom in output_residue.get_atoms()
        if atom.get_name().strip() in NTERM_HYDROGEN_ALIASES
    ]
    for atom_name in names_to_remove:
        output_residue.detach_child(atom_name)

    for atom_index, atom in enumerate(existing_h_atoms):
        serial_counter[0] += 1
        output_residue.add(
            _clone_atom_with_new_name(
                atom=atom,
                serial_number=serial_counter[0],
                new_name=target_names[atom_index],
            )
        )

    missing_count = target_count - len(existing_h_atoms)
    if missing_count <= 0:
        return 0

    missing_coords = _build_missing_nterm_h_positions(
        residue=source_residue,
        missing_count=missing_count,
    )

    for coord_index, coord in enumerate(missing_coords):
        target_name = target_names[len(existing_h_atoms) + coord_index]
        serial_counter[0] += 1
        output_residue.add(
            _make_atom(
                name=target_name,
                coord=coord,
                serial_number=serial_counter[0],
                element="H",
            )
        )

    return missing_count


def _normalize_or_add_c_terminal_oxygens(
    source_residue: Residue.Residue,
    output_residue: Residue.Residue,
    serial_counter: list[int],
) -> bool:
    """
    Ensure AMBER-style C-terminal oxygen naming:
    - O
    - OXT
    """
    names_to_remove = [
        atom.get_id()
        for atom in output_residue.get_atoms()
        if atom.get_name().strip() in CTERM_OXYGEN_ALIASES
    ]
    for atom_name in names_to_remove:
        output_residue.detach_child(atom_name)

    oxygen_like_atoms: list[Atom.Atom] = []
    for atom in source_residue.get_atoms():
        if atom.get_name().strip() in CTERM_OXYGEN_ALIASES:
            oxygen_like_atoms.append(atom)

    if not oxygen_like_atoms:
        raise ValueError(
            f"Residue {source_residue.get_resname()} {source_residue.id} has no terminal oxygen atoms."
        )

    main_o_atom: Atom.Atom | None = None
    extra_o_atom: Atom.Atom | None = None

    for atom in oxygen_like_atoms:
        atom_name = atom.get_name().strip()
        if atom_name == "O" and main_o_atom is None:
            main_o_atom = atom
        elif atom_name != "O" and extra_o_atom is None:
            extra_o_atom = atom

    if main_o_atom is None:
        main_o_atom = oxygen_like_atoms[0]

    serial_counter[0] += 1
    output_residue.add(
        _clone_atom_with_new_name(
            atom=main_o_atom,
            serial_number=serial_counter[0],
            new_name="O",
        )
    )

    if extra_o_atom is not None:
        serial_counter[0] += 1
        output_residue.add(
            _clone_atom_with_new_name(
                atom=extra_o_atom,
                serial_number=serial_counter[0],
                new_name="OXT",
            )
        )
        return False

    c_atom = _require_atom(source_residue, "C")
    ca_atom = _require_atom(source_residue, "CA")

    c_coord = _as_vec(c_atom.coord)
    ca_coord = _as_vec(ca_atom.coord)
    o_coord = _as_vec(main_o_atom.coord)

    c_to_o = _unit(o_coord - c_coord)
    c_to_ca = _unit(ca_coord - c_coord)

    reflected_direction = _unit(
        c_to_o - 2.0 * np.dot(c_to_o, c_to_ca) * c_to_ca,
        fallback=-c_to_o,
    )
    oxt_coord = c_coord + 1.25 * reflected_direction

    serial_counter[0] += 1
    output_residue.add(
        _make_atom(
            name="OXT",
            coord=oxt_coord,
            serial_number=serial_counter[0],
            element="O",
        )
    )
    return True


def _convert_first_residue_to_nterm(
    residue: Residue.Residue,
    serial_counter: list[int],
) -> tuple[Residue.Residue, int]:
    base_resname = residue.get_resname().strip().upper()
    new_resname = N_TERMINAL_RESNAME_MAP[base_resname]

    output_residue = _copy_residue_with_skips(
        residue=residue,
        serial_counter=serial_counter,
        new_resname=new_resname,
        skip_atom_names=NTERM_HYDROGEN_ALIASES,
    )

    n_h_added = _normalize_or_add_nterm_hydrogens(
        source_residue=residue,
        output_residue=output_residue,
        serial_counter=serial_counter,
    )

    return output_residue, n_h_added


def _convert_last_residue_to_cterm(
    residue: Residue.Residue,
    serial_counter: list[int],
) -> tuple[Residue.Residue, bool]:
    base_resname = residue.get_resname().strip().upper()
    new_resname = C_TERMINAL_RESNAME_MAP[base_resname]

    output_residue = _copy_residue_with_skips(
        residue=residue,
        serial_counter=serial_counter,
        new_resname=new_resname,
        skip_atom_names=CTERM_OXYGEN_ALIASES,
    )

    oxt_added = _normalize_or_add_c_terminal_oxygens(
        source_residue=residue,
        output_residue=output_residue,
        serial_counter=serial_counter,
    )

    return output_residue, oxt_added


def _normalize_single_residue_chain(
    residue: Residue.Residue,
    serial_counter: list[int],
) -> tuple[Residue.Residue, int, bool]:
    """
    Special case: one polymer residue is simultaneously chain start and chain end.

    Generic AMBER naming with both NXXX and CXXX is not handled cleanly in tleap for
    arbitrary standard protein loading, so here we keep the base residue name but
    normalize the atom-level terminal features.
    """
    output_residue = _copy_residue_with_skips(
        residue=residue,
        serial_counter=serial_counter,
        new_resname=residue.get_resname(),
        skip_atom_names=NTERM_HYDROGEN_ALIASES | CTERM_OXYGEN_ALIASES,
    )

    n_h_added = _normalize_or_add_nterm_hydrogens(
        source_residue=residue,
        output_residue=output_residue,
        serial_counter=serial_counter,
    )
    oxt_added = _normalize_or_add_c_terminal_oxygens(
        source_residue=residue,
        output_residue=output_residue,
        serial_counter=serial_counter,
    )

    return output_residue, n_h_added, oxt_added


def convert_chain_termini_to_amber(
    chain: Chain.Chain,
    serial_counter: list[int],
) -> tuple[Chain.Chain, ChainTerminusResult]:
    residue_list = list(chain.get_residues())
    output_chain = Chain.Chain(chain.id)

    first_last_indices = _find_first_and_last_polymer_indices(chain)
    if first_last_indices is None:
        for residue in residue_list:
            output_chain.add(
                _copy_residue_with_skips(
                    residue=residue,
                    serial_counter=serial_counter,
                )
            )

        return output_chain, ChainTerminusResult(
            chain_id=chain.id,
            first_resseq=None,
            last_resseq=None,
            first_old_resname=None,
            first_new_resname=None,
            last_old_resname=None,
            last_new_resname=None,
            n_h_added=0,
            oxt_added=False,
            single_residue_chain=False,
        )

    first_index, last_index = first_last_indices
    first_residue = residue_list[first_index]
    last_residue = residue_list[last_index]

    first_old_resname = first_residue.get_resname().strip().upper()
    last_old_resname = last_residue.get_resname().strip().upper()

    n_h_added = 0
    oxt_added = False
    single_residue_chain = first_index == last_index

    for residue_index, residue in enumerate(residue_list):
        if residue_index == first_index == last_index:
            normalized_residue, n_h_added, oxt_added = _normalize_single_residue_chain(
                residue=residue,
                serial_counter=serial_counter,
            )
            output_chain.add(normalized_residue)
            continue

        if residue_index == first_index:
            nterm_residue, n_h_added = _convert_first_residue_to_nterm(
                residue=residue,
                serial_counter=serial_counter,
            )
            output_chain.add(nterm_residue)
            continue

        if residue_index == last_index:
            cterm_residue, oxt_added = _convert_last_residue_to_cterm(
                residue=residue,
                serial_counter=serial_counter,
            )
            output_chain.add(cterm_residue)
            continue

        output_chain.add(
            _copy_residue_with_skips(
                residue=residue,
                serial_counter=serial_counter,
            )
        )

    return output_chain, ChainTerminusResult(
        chain_id=chain.id,
        first_resseq=first_residue.id[1],
        last_resseq=last_residue.id[1],
        first_old_resname=first_old_resname,
        first_new_resname=(
            first_residue.get_resname().strip()
            if single_residue_chain
            else N_TERMINAL_RESNAME_MAP[first_old_resname]
        ),
        last_old_resname=last_old_resname,
        last_new_resname=(
            last_residue.get_resname().strip()
            if single_residue_chain
            else C_TERMINAL_RESNAME_MAP[last_old_resname]
        ),
        n_h_added=n_h_added,
        oxt_added=oxt_added,
        single_residue_chain=single_residue_chain,
    )


def get_default_amber_termini_output_path(
    pdb_directory: str | Path,
    pdb_id: str,
) -> Path:
    """
    Return the standard pipeline output path for the protein with AMBER termini.

    Example
    -------
    data/proteins/1ABC/ -> data/proteins/1ABC/components/1ABC_protein_amber_termini.pdb
    """
    pdb_directory = Path(pdb_directory)
    return pdb_directory / "components" / f"{pdb_id}_protein_amber_termini.pdb"


def convert_protein_termini_to_amber(
    input_pdb_path: str | Path,
    output_pdb_path: str | Path,
) -> list[ChainTerminusResult]:
    """
    Convert true chain termini of a protein-only PDB to AMBER-compatible termini.

    Parameters
    ----------
    input_pdb_path
        Input protein-only PDB path.
    output_pdb_path
        Output protein-only PDB path.

    Returns
    -------
    list[ChainTerminusResult]
        Per-chain summary.
    """
    input_pdb_path = Path(input_pdb_path)
    output_pdb_path = Path(output_pdb_path)

    if not input_pdb_path.exists():
        raise FileNotFoundError(f"Input PDB not found: {input_pdb_path}")

    parser = PDBParser(QUIET=True)
    input_structure = parser.get_structure("protein_input", str(input_pdb_path))

    if len(input_structure) == 0:
        raise ValueError(f"No model found in input PDB: {input_pdb_path}")

    input_model = next(input_structure.get_models())

    output_structure = Structure.Structure("protein_amber_termini")
    output_model = Model.Model(input_model.id)
    output_structure.add(output_model)

    serial_counter = [0]
    result_list: list[ChainTerminusResult] = []

    for chain in input_model.get_chains():
        output_chain, result = convert_chain_termini_to_amber(
            chain=chain,
            serial_counter=serial_counter,
        )
        output_model.add(output_chain)
        result_list.append(result)

    output_pdb_path.parent.mkdir(parents=True, exist_ok=True)

    io = PDBIO()
    io.set_structure(output_structure)
    io.save(str(output_pdb_path))

    return result_list


def convert_protein_termini_for_pdb_directory(
    pdb_directory: str | Path,
    pdb_id: str,
    input_pdb_path: str | Path,
) -> tuple[Path, list[ChainTerminusResult]]:
    """
    Convert protein chain termini to AMBER-compatible termini and write to the
    standard pipeline output path.
    """
    output_pdb_path = get_default_amber_termini_output_path(
        pdb_directory=pdb_directory,
        pdb_id=pdb_id,
    )

    result_list = convert_protein_termini_to_amber(
        input_pdb_path=input_pdb_path,
        output_pdb_path=output_pdb_path,
    )

    return output_pdb_path, result_list


def summarize_terminus_results(result_list: list[ChainTerminusResult]) -> str:
    summary_lines: list[str] = []

    for result in result_list:
        summary_lines.append(
            (
                f"chain={result.chain_id} "
                f"first_resseq={result.first_resseq} "
                f"last_resseq={result.last_resseq} "
                f"first_old_resname={result.first_old_resname} "
                f"first_new_resname={result.first_new_resname} "
                f"last_old_resname={result.last_old_resname} "
                f"last_new_resname={result.last_new_resname} "
                f"n_h_added={result.n_h_added} "
                f"oxt_added={result.oxt_added} "
                f"single_residue_chain={result.single_residue_chain}"
            )
        )

    return "\n".join(summary_lines)


if __name__ == "__main__":
    import argparse

    argument_parser = argparse.ArgumentParser(
        description="Convert true protein chain termini to AMBER-compatible termini."
    )
    argument_parser.add_argument(
        "pdb_directory",
        type=Path,
        help="Protein directory, for example data/proteins/1ABC",
    )
    argument_parser.add_argument(
        "pdb_id",
        type=str,
        help="PDB ID, for example 1ABC",
    )
    argument_parser.add_argument(
        "input_pdb",
        type=Path,
        help="Input protein-only PDB path.",
    )

    arguments = argument_parser.parse_args()

    output_pdb_path, result_list = convert_protein_termini_for_pdb_directory(
        pdb_directory=arguments.pdb_directory,
        pdb_id=arguments.pdb_id,
        input_pdb_path=arguments.input_pdb,
    )

    print(f"output_pdb_path={output_pdb_path}")
    print(summarize_terminus_results(result_list))
