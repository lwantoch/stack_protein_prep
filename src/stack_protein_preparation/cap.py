"""
/home/grheco/repositorios/stack_protein_prep/src/stack_protein_preparation/cap.py

Detect internal missing-segment boundaries inside protein chains and cap them
with AMBER-compatible NME / ACE residues.

Purpose
-------
- detect polymer fragments inside one chain based on actual peptide connectivity
- find internal boundaries between fragments, even if residue numbering jumps
  (for example i -> i+3 or i -> i+23)
- cap the left fragment end with NME
- cap the right fragment start with ACE
- write the capped protein to the standard pipeline output path

Standard output naming
----------------------
Given:
    data/proteins/<PDBID>/

this module writes by default:
    data/proteins/<PDBID>/components/<PDBID>_protein_internal_capped.pdb

Important
---------
- this is for internal missing segments, not true chain termini
- true chain termini should be handled separately in terminus.py
- intended input: protein-only PDB
- existing amino-acid residue names are preserved
- cap coordinates are heuristic and should be minimized later
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

N_TERMINAL_HYDROGEN_NAMES = {
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

C_TERMINAL_EXTRA_OXYGEN_NAMES = {
    "OXT",
    "OT1",
    "OT2",
}


@dataclass(slots=True)
class Fragment:
    chain_id: str
    start_index: int
    end_index: int
    start_resseq: int
    end_resseq: int
    start_resname: str
    end_resname: str


@dataclass(slots=True)
class InternalGapBoundary:
    chain_id: str
    left_residue_index: int
    right_residue_index: int
    left_resseq: int
    right_resseq: int
    left_resname: str
    right_resname: str


@dataclass(slots=True)
class ChainCappingSummary:
    chain_id: str
    n_fragments: int
    n_internal_boundaries: int
    n_boundaries_capped: int


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
    candidate = np.array([1.0, 0.0, 0.0], dtype=float)
    if abs(float(np.dot(vector_u, candidate))) > 0.9:
        candidate = np.array([0.0, 1.0, 0.0], dtype=float)
    return _unit(np.cross(vector_u, candidate), fallback=np.array([0.0, 0.0, 1.0]))


def _project_orthogonal(vector: np.ndarray, axis: np.ndarray) -> np.ndarray:
    axis_u = _unit(axis)
    return vector - np.dot(vector, axis_u) * axis_u


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


def _clone_atom(atom: Atom.Atom, serial_number: int) -> Atom.Atom:
    return Atom.Atom(
        name=atom.get_name(),
        coord=np.array(atom.get_coord(), dtype=float),
        bfactor=float(atom.get_bfactor()),
        occupancy=_safe_occupancy(atom),
        altloc=atom.get_altloc(),
        fullname=atom.get_fullname(),
        serial_number=serial_number,
        element=_safe_element(atom),
    )


def _clone_residue(
    residue: Residue.Residue,
    serial_counter: list[int],
    skip_atom_names: set[str] | None = None,
) -> Residue.Residue:
    skip_atom_names = skip_atom_names or set()

    new_residue = Residue.Residue(
        id=residue.id,
        resname=residue.get_resname(),
        segid=residue.get_segid(),
    )

    for atom in residue.get_atoms():
        if atom.get_name().strip() in skip_atom_names:
            continue
        serial_counter[0] += 1
        new_residue.add(_clone_atom(atom=atom, serial_number=serial_counter[0]))

    return new_residue


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


def _has_peptide_bond(
    left_residue: Residue.Residue,
    right_residue: Residue.Residue,
    peptide_bond_max_distance: float,
) -> bool:
    if "C" not in left_residue or "N" not in right_residue:
        return False

    left_c = _as_vec(left_residue["C"].coord)
    right_n = _as_vec(right_residue["N"].coord)
    distance = _norm(right_n - left_c)

    return distance <= peptide_bond_max_distance


def _build_ace_cap(
    right_residue: Residue.Residue,
    serial_counter: list[int],
) -> Residue.Residue:
    n_atom = _require_atom(right_residue, "N")
    ca_atom = _require_atom(right_residue, "CA")
    c_atom = _require_atom(right_residue, "C")

    n_coord = _as_vec(n_atom.coord)
    ca_coord = _as_vec(ca_atom.coord)
    c_coord = _as_vec(c_atom.coord)

    n_to_ca_unit = _unit(ca_coord - n_coord)
    ace_c_coord = n_coord - 1.33 * n_to_ca_unit

    plane_hint = _project_orthogonal(c_coord - n_coord, axis=(n_coord - ace_c_coord))
    if _norm(plane_hint) < 1e-8:
        plane_hint = _orthogonal_unit(n_to_ca_unit)

    ace_o_direction = _unit(plane_hint, fallback=_orthogonal_unit(n_to_ca_unit))
    ace_o_coord = ace_c_coord + 1.23 * ace_o_direction

    direction_to_n = _unit(n_coord - ace_c_coord)
    direction_to_o = _unit(ace_o_coord - ace_c_coord)
    ace_ch3_direction = _unit(
        -(direction_to_n + direction_to_o),
        fallback=_orthogonal_unit(direction_to_n),
    )
    ace_ch3_coord = ace_c_coord + 1.50 * ace_ch3_direction

    ace_residue = Residue.Residue(
        id=("H_ACE", right_residue.id[1], right_residue.id[2]),
        resname="ACE",
        segid=right_residue.get_segid(),
    )

    serial_counter[0] += 1
    ace_residue.add(_make_atom("CH3", ace_ch3_coord, serial_counter[0], "C"))
    serial_counter[0] += 1
    ace_residue.add(_make_atom("C", ace_c_coord, serial_counter[0], "C"))
    serial_counter[0] += 1
    ace_residue.add(_make_atom("O", ace_o_coord, serial_counter[0], "O"))

    return ace_residue


def _build_nme_cap(
    left_residue: Residue.Residue,
    serial_counter: list[int],
) -> Residue.Residue:
    ca_atom = _require_atom(left_residue, "CA")
    c_atom = _require_atom(left_residue, "C")
    o_atom = _require_atom(left_residue, "O")

    ca_coord = _as_vec(ca_atom.coord)
    c_coord = _as_vec(c_atom.coord)
    o_coord = _as_vec(o_atom.coord)

    ca_to_c_unit = _unit(c_coord - ca_coord)
    nme_n_coord = c_coord + 1.33 * ca_to_c_unit

    o_hint = _project_orthogonal(c_coord - o_coord, axis=(nme_n_coord - c_coord))
    if _norm(o_hint) < 1e-8:
        o_hint = _orthogonal_unit(ca_to_c_unit)

    nme_ch3_direction = _unit(
        (nme_n_coord - c_coord) + 0.25 * _unit(o_hint),
        fallback=(nme_n_coord - c_coord),
    )
    nme_ch3_coord = nme_n_coord + 1.46 * nme_ch3_direction

    nme_residue = Residue.Residue(
        id=("H_NME", left_residue.id[1], left_residue.id[2]),
        resname="NME",
        segid=left_residue.get_segid(),
    )

    serial_counter[0] += 1
    nme_residue.add(_make_atom("N", nme_n_coord, serial_counter[0], "N"))
    serial_counter[0] += 1
    nme_residue.add(_make_atom("CH3", nme_ch3_coord, serial_counter[0], "C"))

    return nme_residue


def get_default_internal_capped_output_path(
    pdb_directory: str | Path,
    pdb_id: str,
) -> Path:
    """
    Return the standard pipeline output path for the internally capped protein.

    Example
    -------
    data/proteins/1ABC/ -> data/proteins/1ABC/components/1ABC_protein_internal_capped.pdb
    """
    pdb_directory = Path(pdb_directory)
    return pdb_directory / "components" / f"{pdb_id}_protein_internal_capped.pdb"


def find_internal_fragments(
    chain: Chain.Chain,
    peptide_bond_max_distance: float = 1.8,
) -> list[Fragment]:
    residue_list = list(chain.get_residues())

    polymer_positions = [
        index
        for index, residue in enumerate(residue_list)
        if _is_polymer_residue(residue)
    ]

    if not polymer_positions:
        return []

    fragment_list: list[Fragment] = []

    current_start_index = polymer_positions[0]
    current_end_index = polymer_positions[0]

    for prev_index, curr_index in zip(polymer_positions[:-1], polymer_positions[1:]):
        prev_residue = residue_list[prev_index]
        curr_residue = residue_list[curr_index]

        if _has_peptide_bond(
            left_residue=prev_residue,
            right_residue=curr_residue,
            peptide_bond_max_distance=peptide_bond_max_distance,
        ):
            current_end_index = curr_index
            continue

        start_residue = residue_list[current_start_index]
        end_residue = residue_list[current_end_index]

        fragment_list.append(
            Fragment(
                chain_id=chain.id,
                start_index=current_start_index,
                end_index=current_end_index,
                start_resseq=start_residue.id[1],
                end_resseq=end_residue.id[1],
                start_resname=start_residue.get_resname().strip(),
                end_resname=end_residue.get_resname().strip(),
            )
        )

        current_start_index = curr_index
        current_end_index = curr_index

    start_residue = residue_list[current_start_index]
    end_residue = residue_list[current_end_index]

    fragment_list.append(
        Fragment(
            chain_id=chain.id,
            start_index=current_start_index,
            end_index=current_end_index,
            start_resseq=start_residue.id[1],
            end_resseq=end_residue.id[1],
            start_resname=start_residue.get_resname().strip(),
            end_resname=end_residue.get_resname().strip(),
        )
    )

    return fragment_list


def find_internal_gap_boundaries(
    chain: Chain.Chain,
    peptide_bond_max_distance: float = 1.8,
) -> list[InternalGapBoundary]:
    residue_list = list(chain.get_residues())
    fragment_list = find_internal_fragments(
        chain=chain,
        peptide_bond_max_distance=peptide_bond_max_distance,
    )

    if len(fragment_list) < 2:
        return []

    boundary_list: list[InternalGapBoundary] = []

    for left_fragment, right_fragment in zip(fragment_list[:-1], fragment_list[1:]):
        left_residue = residue_list[left_fragment.end_index]
        right_residue = residue_list[right_fragment.start_index]

        boundary_list.append(
            InternalGapBoundary(
                chain_id=chain.id,
                left_residue_index=left_fragment.end_index,
                right_residue_index=right_fragment.start_index,
                left_resseq=left_residue.id[1],
                right_resseq=right_residue.id[1],
                left_resname=left_residue.get_resname().strip(),
                right_resname=right_residue.get_resname().strip(),
            )
        )

    return boundary_list


def cap_internal_gap_boundaries(
    chain: Chain.Chain,
    serial_counter: list[int],
    peptide_bond_max_distance: float = 1.8,
) -> tuple[Chain.Chain, ChainCappingSummary]:
    residue_list = list(chain.get_residues())

    boundary_list = find_internal_gap_boundaries(
        chain=chain,
        peptide_bond_max_distance=peptide_bond_max_distance,
    )
    fragment_list = find_internal_fragments(
        chain=chain,
        peptide_bond_max_distance=peptide_bond_max_distance,
    )

    boundary_by_left_index = {x.left_residue_index: x for x in boundary_list}
    boundary_by_right_index = {x.right_residue_index: x for x in boundary_list}

    output_chain = Chain.Chain(chain.id)
    n_capped = 0

    for residue_index, residue in enumerate(residue_list):
        if residue_index in boundary_by_right_index:
            output_chain.add(
                _build_ace_cap(
                    right_residue=residue,
                    serial_counter=serial_counter,
                )
            )

        skip_atom_names: set[str] = set()

        if residue_index in boundary_by_left_index:
            skip_atom_names |= C_TERMINAL_EXTRA_OXYGEN_NAMES

        if residue_index in boundary_by_right_index:
            skip_atom_names |= N_TERMINAL_HYDROGEN_NAMES

        output_chain.add(
            _clone_residue(
                residue=residue,
                serial_counter=serial_counter,
                skip_atom_names=skip_atom_names,
            )
        )

        if residue_index in boundary_by_left_index:
            output_chain.add(
                _build_nme_cap(
                    left_residue=residue,
                    serial_counter=serial_counter,
                )
            )
            n_capped += 1

    summary = ChainCappingSummary(
        chain_id=chain.id,
        n_fragments=len(fragment_list),
        n_internal_boundaries=len(boundary_list),
        n_boundaries_capped=n_capped,
    )

    return output_chain, summary


def cap_internal_gaps_in_structure(
    input_pdb_path: str | Path,
    output_pdb_path: str | Path,
    peptide_bond_max_distance: float = 1.8,
) -> list[ChainCappingSummary]:
    input_pdb_path = Path(input_pdb_path)
    output_pdb_path = Path(output_pdb_path)

    if not input_pdb_path.exists():
        raise FileNotFoundError(f"Input PDB not found: {input_pdb_path}")

    parser = PDBParser(QUIET=True)
    input_structure = parser.get_structure("protein_input", str(input_pdb_path))

    if len(input_structure) == 0:
        raise ValueError(f"No model found in input PDB: {input_pdb_path}")

    input_model = next(input_structure.get_models())

    output_structure = Structure.Structure("protein_internal_caps")
    output_model = Model.Model(input_model.id)
    output_structure.add(output_model)

    serial_counter = [0]
    summary_list: list[ChainCappingSummary] = []

    for chain in input_model.get_chains():
        capped_chain, summary = cap_internal_gap_boundaries(
            chain=chain,
            serial_counter=serial_counter,
            peptide_bond_max_distance=peptide_bond_max_distance,
        )
        output_model.add(capped_chain)
        summary_list.append(summary)

    output_pdb_path.parent.mkdir(parents=True, exist_ok=True)

    io = PDBIO()
    io.set_structure(output_structure)
    io.save(str(output_pdb_path))

    return summary_list


def cap_internal_gaps_for_pdb_directory(
    pdb_directory: str | Path,
    pdb_id: str,
    input_pdb_path: str | Path,
    peptide_bond_max_distance: float = 1.8,
) -> tuple[Path, list[ChainCappingSummary]]:
    """
    Run internal ACE/NME capping and write to the standard pipeline output path.
    """
    output_pdb_path = get_default_internal_capped_output_path(
        pdb_directory=pdb_directory,
        pdb_id=pdb_id,
    )

    summary_list = cap_internal_gaps_in_structure(
        input_pdb_path=input_pdb_path,
        output_pdb_path=output_pdb_path,
        peptide_bond_max_distance=peptide_bond_max_distance,
    )

    return output_pdb_path, summary_list


def summarize_capping_results(summary_list: list[ChainCappingSummary]) -> str:
    return "\n".join(
        (
            f"chain={summary.chain_id} "
            f"n_fragments={summary.n_fragments} "
            f"n_internal_boundaries={summary.n_internal_boundaries} "
            f"n_boundaries_capped={summary.n_boundaries_capped}"
        )
        for summary in summary_list
    )


if __name__ == "__main__":
    import argparse

    argument_parser = argparse.ArgumentParser(
        description="Detect internal fragment boundaries and cap them with NME/ACE."
    )
    argument_parser.add_argument(
        "pdb_directory",
        type=Path,
        help="Protein directory, e.g. data/proteins/1ABC",
    )
    argument_parser.add_argument(
        "pdb_id",
        type=str,
        help="PDB ID, e.g. 1ABC",
    )
    argument_parser.add_argument(
        "input_pdb",
        type=Path,
        help="Input protein-only PDB path.",
    )
    argument_parser.add_argument(
        "--peptide-bond-max-distance",
        type=float,
        default=1.8,
        help="Maximum C-N distance treated as peptide-bonded. Default: 1.8 Å",
    )
    arguments = argument_parser.parse_args()

    output_pdb_path, summary_list = cap_internal_gaps_for_pdb_directory(
        pdb_directory=arguments.pdb_directory,
        pdb_id=arguments.pdb_id,
        input_pdb_path=arguments.input_pdb,
        peptide_bond_max_distance=arguments.peptide_bond_max_distance,
    )

    print(f"output_pdb_path={output_pdb_path}")
    print(summarize_capping_results(summary_list))
