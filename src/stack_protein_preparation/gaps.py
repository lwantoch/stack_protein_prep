"""
/home/grheco/repositorios/stack_protein_prep/src/stack_protein_preparation/gaps.py

Detect internal residue-numbering gaps in protein PDB files.

Responsibilities
----------------
- identify protein-like residues per chain
- detect internal numbering discontinuities (gaps)
- provide per-gap detailed context (flanking residues)
- group gaps per chain
- summarize gap statistics for pipeline decisions

Important
---------
- this module detects numbering gaps, not guaranteed biological missing residues
- author numbering in PDB files may intentionally skip residue numbers
- true missing residues must be validated later via sequence alignment (e.g. UniProt mapping)

Assumptions
-----------
- input PDB is already cleaned (e.g. insertion codes removed)
- residue numbering still reflects original author numbering (not renumbered)
- non-protein components (water, ligands, metals) are present but ignored
"""

from __future__ import annotations

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
    "ASH",  # protonated ASP
    "GLH",  # protonated GLU
    "HID",
    "HIE",
    "HIP",  # histidine variants
    "CYM",
    "CYX",  # cysteine variants
    "LYN",  # deprotonated LYS
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
}

PROTEIN_RESNAMES = STANDARD_AA | AA_VARIANTS | PROTEIN_LIKE_NONSTANDARD


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


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

    Includes modified residues to avoid artificial gaps.
    """
    return residue.name.strip().upper() in PROTEIN_RESNAMES


def _collect_chain_protein_residues(chain: Any) -> list[dict[str, Any]]:
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

        out.append(
            {
                "chain": chain_id,
                "resnum": resnum,
                "resname": residue.name.strip().upper(),
            }
        )

    return sorted(out, key=lambda r: r["resnum"])


# ---------------------------------------------------------------------------
# Core gap detection
# ---------------------------------------------------------------------------


def gaps_for_pdb(pdb_path: str | Path) -> list[dict[str, Any]]:
    """
    Detect internal residue-numbering gaps.

    Returns
    -------
    list[dict]
        Each gap contains:
        - chain
        - prev_res
        - prev_resname
        - start
        - end
        - next_res
        - next_resname
        - gap_size
    """
    pdb = PDBFile(str(Path(pdb_path)))
    gaps_out: list[dict[str, Any]] = []

    for chain in pdb.topology.chains():
        residues = _collect_chain_protein_residues(chain)

        for prev_res, curr_res in zip(residues, residues[1:]):
            prev_n = prev_res["resnum"]
            curr_n = curr_res["resnum"]

            if curr_n - prev_n <= 1:
                continue

            start = prev_n + 1
            end = curr_n - 1

            gaps_out.append(
                {
                    "chain": prev_res["chain"],
                    "prev_res": prev_n,
                    "prev_resname": prev_res["resname"],
                    "start": start,
                    "end": end,
                    "next_res": curr_n,
                    "next_resname": curr_res["resname"],
                    "gap_size": end - start + 1,
                }
            )

    return gaps_out


# ---------------------------------------------------------------------------
# Grouping & summary
# ---------------------------------------------------------------------------


def gaps_by_chain_for_pdb(pdb_path: str | Path) -> dict[str, list[dict[str, Any]]]:
    """
    Group gaps by chain.
    """
    grouped: dict[str, list[dict[str, Any]]] = {}

    for gap in gaps_for_pdb(pdb_path):
        grouped.setdefault(gap["chain"], []).append(gap)

    return grouped


def summarize_gaps(pdb_path: str | Path) -> dict[str, Any]:
    """
    Return summary statistics for pipeline decisions.
    """
    grouped = gaps_by_chain_for_pdb(pdb_path)
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


def debug_excluded_residues(pdb_path: str | Path) -> dict[str, list[tuple[str, str]]]:
    """
    List residues ignored during gap detection.

    Useful to debug false gaps caused by filtering.
    """
    pdb = PDBFile(str(Path(pdb_path)))
    excluded: dict[str, list[tuple[str, str]]] = {}

    for chain in pdb.topology.chains():
        chain_id = str(getattr(chain, "id", None) or "?")

        for residue in chain.residues():
            if not _is_chain_protein_residue(residue):
                excluded.setdefault(chain_id, []).append(
                    (residue.name.strip().upper(), str(residue.id))
                )

    return excluded
