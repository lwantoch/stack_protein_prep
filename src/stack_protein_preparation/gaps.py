# /home/grheco/repositorios/stack_protein_prep/src/stack_protein_preparation/gaps.py

"""
Detect likely internal residue gaps in protein PDB files.

Logic
-----
A gap is only reported when BOTH conditions are true:
1. there is a residue-number jump larger than 1
2. the flanking residues are not peptide-connected by C-N distance

This reduces false positives from:
- author numbering quirks
- harmless numbering skips
- modified residue naming
- some terminal / fragment-edge artifacts

Important
---------
- this module still does NOT prove a biologically missing segment
- final validation should still come from sequence alignment to UniProt
- this module is intended as a better structural pre-filter for pipeline logic

Standard use
------------
Typical pipeline entry point:
    summarize_gaps(pdb_path)

Optional range-aware use:
    summarize_gaps(pdb_path, residue_range="33-480")

Returned summary contains:
- has_gaps
- n_gaps
- max_gap_size
- total_missing_residues
- gap_sizes
- chains_with_gaps
- relevant_chain_ids

Each individual gap also contains:
- chain
- prev_res / next_res
- prev_resname / next_resname
- start / end / gap_size
- c_atom_present / n_atom_present
- peptide_cn_distance
- peptide_bond_cutoff_angstrom
- peptide_connected
- near_chain_start
- near_chain_end
- classification

Classification
--------------
- internal:
    numbering jump + not peptide-connected + not near chain edge
- terminal_like:
    numbering jump + not peptide-connected + gap near chain edge
- ambiguous:
    numbering jump + missing atoms needed for geometry check
"""

from __future__ import annotations

import math
import re
from pathlib import Path
from typing import Any

try:
    from openmm.app import PDBFile
except ImportError as exc:
    raise ImportError(
        "openmm is required for gaps.py. Install it in your environment."
    ) from exc


# ---------------------------------------------------------------------------
# Residue classification
# ---------------------------------------------------------------------------

STANDARD_AA = {
    "ALA",
    "ARG",
    "ASN",
    "ASP",
    "CYS",
    "GLN",
    "GLU",
    "GLY",
    "HIS",
    "ILE",
    "LEU",
    "LYS",
    "MET",
    "PHE",
    "PRO",
    "SER",
    "THR",
    "TRP",
    "TYR",
    "VAL",
}

AA_VARIANTS = {
    "ASH",
    "GLH",
    "HID",
    "HIE",
    "HIP",
    "CYM",
    "CYX",
    "LYN",
    "NALA",
    "NARG",
    "NASN",
    "NASP",
    "NCYS",
    "NGLN",
    "NGLU",
    "NGLY",
    "NHIS",
    "NILE",
    "NLEU",
    "NLYS",
    "NMET",
    "NPHE",
    "NPRO",
    "NSER",
    "NTHR",
    "NTRP",
    "NTYR",
    "NVAL",
    "CALA",
    "CARG",
    "CASN",
    "CASP",
    "CCYS",
    "CGLN",
    "CGLU",
    "CGLY",
    "CHIS",
    "CILE",
    "CLEU",
    "CLYS",
    "CMET",
    "CPHE",
    "CPRO",
    "CSER",
    "CTHR",
    "CTRP",
    "CTYR",
    "CVAL",
}

PROTEIN_LIKE_NONSTANDARD = {
    "MSE",
    "SEC",
    "PYL",
    "SEP",
    "TPO",
    "PTR",
    "CSO",
    "CSD",
    "OCS",
    "KCX",
    "LLP",
    "MLY",
    "HYP",
    "ACE",
    "NME",
}

PROTEIN_RESNAMES = STANDARD_AA | AA_VARIANTS | PROTEIN_LIKE_NONSTANDARD


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _strip_atom_name(atom_name: str) -> str:
    return atom_name.strip().upper()


def _distance_angstrom(atom_a: Any, atom_b: Any) -> float:
    """
    Return distance in Angstrom between two OpenMM atoms using the loaded positions.
    """
    dx = atom_a["x"] - atom_b["x"]
    dy = atom_a["y"] - atom_b["y"]
    dz = atom_a["z"] - atom_b["z"]
    return math.sqrt(dx * dx + dy * dy + dz * dz)


def _res_seq_int(residue: Any) -> int | None:
    """
    Convert residue.id to integer if possible.
    """
    try:
        return int(str(residue.id).strip())
    except (TypeError, ValueError):
        return None


def _is_chain_protein_residue(residue: Any) -> bool:
    """
    Check whether a residue belongs to the protein backbone.

    Includes common modified residues to avoid artificial gaps.
    """
    return residue.name.strip().upper() in PROTEIN_RESNAMES


def _collect_atom_positions_by_index(pdb: PDBFile) -> dict[int, dict[str, float]]:
    """
    Build atom-index -> xyz mapping in Angstrom.
    """
    positions = pdb.positions
    out: dict[int, dict[str, float]] = {}

    for atom in pdb.topology.atoms():
        pos_nm = positions[atom.index]
        out[atom.index] = {
            "x": float(pos_nm.x) * 10.0,
            "y": float(pos_nm.y) * 10.0,
            "z": float(pos_nm.z) * 10.0,
        }

    return out


def _collect_chain_protein_residues(
    chain: Any,
    atom_xyz_by_index: dict[int, dict[str, float]],
) -> list[dict[str, Any]]:
    """
    Collect normalized residue records for one chain.

    Returns
    -------
    list[dict]
        Sorted by residue number.
    """
    chain_id = str(getattr(chain, "id", None) or "?")
    out: list[dict[str, Any]] = []

    for residue in chain.residues():
        if not _is_chain_protein_residue(residue):
            continue

        resnum = _res_seq_int(residue)
        if resnum is None:
            continue

        atom_map: dict[str, dict[str, float]] = {}
        for atom in residue.atoms():
            atom_name = _strip_atom_name(atom.name)
            coords = atom_xyz_by_index.get(atom.index)
            if coords is not None:
                atom_map[atom_name] = coords

        out.append(
            {
                "chain": chain_id,
                "resnum": resnum,
                "resname": residue.name.strip().upper(),
                "atoms": atom_map,
            }
        )

    return sorted(out, key=lambda r: r["resnum"])


def _get_peptide_cn_distance(
    prev_res: dict[str, Any], next_res: dict[str, Any]
) -> float | None:
    """
    Return C(prev) -> N(next) distance in Angstrom if both atoms exist.
    """
    prev_c = prev_res["atoms"].get("C")
    next_n = next_res["atoms"].get("N")

    if prev_c is None or next_n is None:
        return None

    return _distance_angstrom(prev_c, next_n)


def _classify_gap(
    *,
    prev_idx: int,
    next_idx: int,
    n_residues_in_chain: int,
    peptide_connected: bool | None,
    edge_window: int,
) -> str:
    """
    Classify a detected numbering discontinuity.

    When backbone atoms (C, N) are absent the connectivity cannot be verified,
    so the classification falls back on the numbering jump alone and treats the
    gap as internal (ambiguous is reserved for cases where C/N atoms are
    present but bond distances are borderline).
    """
    near_chain_start = prev_idx < edge_window
    near_chain_end = (n_residues_in_chain - 1 - next_idx) < edge_window

    if peptide_connected is None:
        # Cannot verify via geometry; treat as internal gap (trust numbering).
        if near_chain_start and near_chain_end:
            return "internal"
        if near_chain_start or near_chain_end:
            return "internal"
        return "internal"

    if peptide_connected:
        return "numbering_only"

    # When both sides are near the chain edge (e.g. a 2-residue chain with a gap),
    # the gap spans the entire observable range and is treated as internal.
    if near_chain_start and near_chain_end:
        return "internal"

    if near_chain_start or near_chain_end:
        return "terminal_like"

    return "internal"


def _parse_chain_qualified_range(
    residue_range: str,
) -> tuple[str, int, str, int] | None:
    """Return (start_chain, start_resnum, end_chain, end_resnum) for 'A16-B50'.

    Returns None when the range has no chain qualifiers.
    """
    m = re.fullmatch(r"([A-Za-z])(-?\d+)-([A-Za-z])(-?\d+)", str(residue_range).strip())
    if m:
        return m.group(1).upper(), int(m.group(2)), m.group(3).upper(), int(m.group(4))
    return None


def _parse_residue_range(residue_range: str) -> tuple[int, int]:
    """
    Parse residue range like '33-480' or 'A16-B50' (chain-qualified).
    """
    stripped = str(residue_range).strip()
    # Chain-qualified format: A33-B480 — strip chain letters, keep residue numbers
    m = re.fullmatch(r"[A-Za-z](-?\d+)-[A-Za-z](-?\d+)", stripped)
    if m:
        return int(m.group(1)), int(m.group(2))
    # Legacy format: 33-480
    match = re.fullmatch(r"\s*(-?\d+)\s*-\s*(-?\d+)\s*", stripped)
    if match is None:
        raise ValueError(f"Invalid residue range: {residue_range!r}")

    start_residue = int(match.group(1))
    end_residue = int(match.group(2))

    if start_residue > end_residue:
        raise ValueError(f"Invalid residue range: start > end in {residue_range!r}")

    return start_residue, end_residue


def _choose_relevant_chain_ids_from_range(
    pdb: PDBFile,
    atom_xyz_by_index: dict[int, dict[str, float]],
    residue_range: str,
) -> tuple[str, ...]:
    """
    Choose chain(s) most relevant for the requested residue range.

    Strategy
    --------
    - When the range has explicit chain qualifiers (A16-B50), return those chains.
    - Otherwise compute overlap between requested range and observed protein residue
      numbers for each chain; keep all chains tied for best overlap.
    - If no chain overlaps, keep all chains with protein residues.
    """
    residue_numbers_by_chain: dict[str, set[int]] = {}

    for chain in pdb.topology.chains():
        residue_records = _collect_chain_protein_residues(chain, atom_xyz_by_index)
        if not residue_records:
            continue

        chain_id = residue_records[0]["chain"]
        residue_numbers_by_chain[chain_id] = {r["resnum"] for r in residue_records}

    if not residue_numbers_by_chain:
        return ()

    # When the range has explicit chain qualifiers (A16-B50), return those chains directly.
    chain_qual = _parse_chain_qualified_range(residue_range)
    if chain_qual is not None:
        start_chain, _, end_chain, _ = chain_qual
        explicit_chains = tuple(
            c for c in sorted({start_chain, end_chain})
            if c in residue_numbers_by_chain
        )
        if explicit_chains:
            return explicit_chains

    start_residue, end_residue = _parse_residue_range(residue_range)
    requested_range = set(range(start_residue, end_residue + 1))

    score_rows: list[tuple[str, int, int]] = []
    for chain_id, residue_number_set in sorted(residue_numbers_by_chain.items()):
        overlap_count = len(residue_number_set & requested_range)
        total_chain_residue_count = len(residue_number_set)
        score_rows.append((chain_id, overlap_count, total_chain_residue_count))

    best_overlap = max(row[1] for row in score_rows)
    if best_overlap == 0:
        return tuple(sorted(residue_numbers_by_chain.keys()))

    best_rows = [row for row in score_rows if row[1] == best_overlap]
    return tuple(sorted(row[0] for row in best_rows))


def _get_relevant_chain_id_set(
    pdb: PDBFile,
    atom_xyz_by_index: dict[int, dict[str, float]],
    residue_range: str,
) -> set[str]:
    if not str(residue_range).strip():
        return set()

    relevant_chain_ids = _choose_relevant_chain_ids_from_range(
        pdb=pdb,
        atom_xyz_by_index=atom_xyz_by_index,
        residue_range=residue_range,
    )
    return set(relevant_chain_ids)


# ---------------------------------------------------------------------------
# Core gap detection
# ---------------------------------------------------------------------------


def gaps_for_pdb(
    pdb_path: str | Path,
    *,
    peptide_bond_cutoff_angstrom: float = 1.8,
    edge_window: int = 1,
    include_terminal_like: bool = False,
    include_ambiguous: bool = False,
    residue_range: str = "",
) -> list[dict[str, Any]]:
    """
    Detect likely internal residue gaps.

    Parameters
    ----------
    pdb_path
        Path to the input PDB.
    peptide_bond_cutoff_angstrom
        Maximum C(prev)-N(next) distance still considered peptide-connected.
        Typical peptide bond distance is around 1.33 Å. A cutoff around 1.8 Å
        is intentionally tolerant.
    edge_window
        Number of observed residues from each chain edge considered "near edge".
    include_terminal_like
        Whether to include terminal-like discontinuities in output.
    include_ambiguous
        Whether to include ambiguous cases where geometry could not be checked.
    residue_range
        Optional requested residue range such as '33-480'. If provided, only
        chain(s) most relevant for that range are evaluated.

    Returns
    -------
    list[dict]
        Each kept gap contains:
        - chain
        - prev_res
        - prev_resname
        - start
        - end
        - next_res
        - next_resname
        - gap_size
        - c_atom_present
        - n_atom_present
        - peptide_cn_distance
        - peptide_bond_cutoff_angstrom
        - peptide_connected
        - near_chain_start
        - near_chain_end
        - classification
    """
    pdb = PDBFile(str(Path(pdb_path)))
    atom_xyz_by_index = _collect_atom_positions_by_index(pdb)
    relevant_chain_id_set = _get_relevant_chain_id_set(
        pdb=pdb,
        atom_xyz_by_index=atom_xyz_by_index,
        residue_range=residue_range,
    )

    gaps_out: list[dict[str, Any]] = []

    for chain in pdb.topology.chains():
        residues = _collect_chain_protein_residues(chain, atom_xyz_by_index)
        if not residues:
            continue

        chain_id = residues[0]["chain"]
        if relevant_chain_id_set and chain_id not in relevant_chain_id_set:
            continue

        n_residues = len(residues)
        if n_residues < 2:
            continue

        for prev_idx, (prev_res, next_res) in enumerate(zip(residues, residues[1:])):
            next_idx = prev_idx + 1

            prev_n = prev_res["resnum"]
            next_n = next_res["resnum"]

            numbering_jump = next_n - prev_n
            if numbering_jump <= 1:
                continue

            start = prev_n + 1
            end = next_n - 1

            peptide_cn_distance = _get_peptide_cn_distance(prev_res, next_res)

            if peptide_cn_distance is None:
                peptide_connected = None
            else:
                peptide_connected = peptide_cn_distance <= peptide_bond_cutoff_angstrom

            near_chain_start = prev_idx < edge_window
            near_chain_end = (n_residues - 1 - next_idx) < edge_window

            classification = _classify_gap(
                prev_idx=prev_idx,
                next_idx=next_idx,
                n_residues_in_chain=n_residues,
                peptide_connected=peptide_connected,
                edge_window=edge_window,
            )

            keep = False
            if classification == "internal":
                keep = True
            elif classification == "terminal_like" and include_terminal_like:
                keep = True
            elif classification == "ambiguous" and include_ambiguous:
                keep = True

            if not keep:
                continue

            gaps_out.append(
                {
                    "chain": prev_res["chain"],
                    "prev_res": prev_n,
                    "prev_resname": prev_res["resname"],
                    "start": start,
                    "end": end,
                    "next_res": next_n,
                    "next_resname": next_res["resname"],
                    "gap_size": end - start + 1,
                }
            )

    return gaps_out


# ---------------------------------------------------------------------------
# Grouping & summary
# ---------------------------------------------------------------------------


def gaps_by_chain_for_pdb(
    pdb_path: str | Path,
    *,
    peptide_bond_cutoff_angstrom: float = 1.8,
    edge_window: int = 1,
    include_terminal_like: bool = False,
    include_ambiguous: bool = False,
    residue_range: str = "",
) -> dict[str, list[dict[str, Any]]]:
    """
    Group kept gaps by chain.
    """
    grouped: dict[str, list[dict[str, Any]]] = {}

    for gap in gaps_for_pdb(
        pdb_path,
        peptide_bond_cutoff_angstrom=peptide_bond_cutoff_angstrom,
        edge_window=edge_window,
        include_terminal_like=include_terminal_like,
        include_ambiguous=include_ambiguous,
        residue_range=residue_range,
    ):
        grouped.setdefault(gap["chain"], []).append(gap)

    return grouped


def summarize_gaps(
    pdb_path: str | Path,
    *,
    peptide_bond_cutoff_angstrom: float = 1.8,
    edge_window: int = 1,
    include_terminal_like: bool = False,
    include_ambiguous: bool = False,
    residue_range: str = "",
) -> dict[str, Any]:
    """
    Return summary statistics for pipeline decisions.

    By default, only likely internal gaps are counted.

    If residue_range is provided, only chain(s) most relevant for that range
    are considered.
    """
    grouped = gaps_by_chain_for_pdb(
        pdb_path,
        peptide_bond_cutoff_angstrom=peptide_bond_cutoff_angstrom,
        edge_window=edge_window,
        include_terminal_like=include_terminal_like,
        include_ambiguous=include_ambiguous,
        residue_range=residue_range,
    )
    all_gaps = [g for gaps in grouped.values() for g in gaps]
    sizes = [g["gap_size"] for g in all_gaps]

    return {
        "has_gaps": bool(all_gaps),
        "n_gaps": len(all_gaps),
        "max_gap_size": max(sizes, default=0),
        "total_missing_residues": sum(sizes),
        "gap_sizes": sizes,
        "chains_with_gaps": sorted(grouped),
    }


# ---------------------------------------------------------------------------
# Debugging
# ---------------------------------------------------------------------------


def debug_excluded_residues(
    pdb_path: str | Path,
    *,
    residue_range: str = "",
) -> dict[str, list[tuple[str, str]]]:
    """
    List residues ignored during gap detection.

    Useful to debug false gaps caused by residue filtering.
    """
    pdb = PDBFile(str(Path(pdb_path)))
    atom_xyz_by_index = _collect_atom_positions_by_index(pdb)
    relevant_chain_id_set = _get_relevant_chain_id_set(
        pdb=pdb,
        atom_xyz_by_index=atom_xyz_by_index,
        residue_range=residue_range,
    )

    excluded: dict[str, list[tuple[str, str]]] = {}

    for chain in pdb.topology.chains():
        chain_id = str(getattr(chain, "id", None) or "?")
        if relevant_chain_id_set and chain_id not in relevant_chain_id_set:
            continue

        for residue in chain.residues():
            if not _is_chain_protein_residue(residue):
                excluded.setdefault(chain_id, []).append(
                    (residue.name.strip().upper(), str(residue.id))
                )

    return excluded


def debug_all_numbering_discontinuities(
    pdb_path: str | Path,
    *,
    peptide_bond_cutoff_angstrom: float = 1.8,
    edge_window: int = 1,
    residue_range: str = "",
) -> list[dict[str, Any]]:
    """
    Return ALL numbering discontinuities, including those filtered out.

    Useful to compare:
    - numbering-only cases
    - terminal-like cases
    - ambiguous cases
    - kept internal cases
    """
    pdb = PDBFile(str(Path(pdb_path)))
    atom_xyz_by_index = _collect_atom_positions_by_index(pdb)
    relevant_chain_id_set = _get_relevant_chain_id_set(
        pdb=pdb,
        atom_xyz_by_index=atom_xyz_by_index,
        residue_range=residue_range,
    )

    out: list[dict[str, Any]] = []

    for chain in pdb.topology.chains():
        residues = _collect_chain_protein_residues(chain, atom_xyz_by_index)
        if not residues:
            continue

        chain_id = residues[0]["chain"]
        if relevant_chain_id_set and chain_id not in relevant_chain_id_set:
            continue

        n_residues = len(residues)
        if n_residues < 2:
            continue

        for prev_idx, (prev_res, next_res) in enumerate(zip(residues, residues[1:])):
            next_idx = prev_idx + 1

            prev_n = prev_res["resnum"]
            next_n = next_res["resnum"]

            numbering_jump = next_n - prev_n
            if numbering_jump <= 1:
                continue

            start = prev_n + 1
            end = next_n - 1

            peptide_cn_distance = _get_peptide_cn_distance(prev_res, next_res)
            if peptide_cn_distance is None:
                peptide_connected = None
            else:
                peptide_connected = peptide_cn_distance <= peptide_bond_cutoff_angstrom

            near_chain_start = prev_idx < edge_window
            near_chain_end = (n_residues - 1 - next_idx) < edge_window

            classification = _classify_gap(
                prev_idx=prev_idx,
                next_idx=next_idx,
                n_residues_in_chain=n_residues,
                peptide_connected=peptide_connected,
                edge_window=edge_window,
            )

            out.append(
                {
                    "chain": prev_res["chain"],
                    "prev_res": prev_n,
                    "prev_resname": prev_res["resname"],
                    "start": start,
                    "end": end,
                    "next_res": next_n,
                    "next_resname": next_res["resname"],
                    "gap_size": end - start + 1,
                    "c_atom_present": "C" in prev_res["atoms"],
                    "n_atom_present": "N" in next_res["atoms"],
                    "peptide_cn_distance": peptide_cn_distance,
                    "peptide_bond_cutoff_angstrom": peptide_bond_cutoff_angstrom,
                    "peptide_connected": peptide_connected,
                    "near_chain_start": near_chain_start,
                    "near_chain_end": near_chain_end,
                    "classification": classification,
                }
            )

    return out
