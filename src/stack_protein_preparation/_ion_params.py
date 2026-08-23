"""Ion 12-6-4 parameter lookup for FRUTON tleap script generation.

Purpose: tell downstream tleap which ``leaprc.water.tip3p_ion12-6-4``
family covers each metal resname, and flag ions that require MCPB.py
bonded modeling instead of the LJ parameter route.

Data source: ``data/ion_12_6_4_reference.csv`` (hand-curated from
Li-Merz JCTC 2013 (divalents), 2015 (mono/tri/tetra), 2020 (Ln3+)
+ Joung-Cheatham JPCB 2008 (alkali).  Column notes on each entry
explain the caveats and biological usage per ion.

Public API:

    lookup_ion(pdb_resname)   -> IonParams | None
    covered_by_12_6_4(res)    -> bool
    tleap_source_lines()      -> list[str]  # leaprc's to source
    describe(pdb_resname)     -> str

License-free: pure Python + csv stdlib.  Loads once per process.
"""
from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path


_CSV_PATH = Path(__file__).parent / "data" / "ion_12_6_4_reference.csv"


@dataclass(frozen=True)
class IonParams:
    pdb_resname: str
    element: str
    formal_charge: str
    amber_atom_type: str
    ion_size_class: str
    li_merz_paper: str
    leaprc_route: str
    notes: str

    @property
    def needs_mcpb(self) -> bool:
        """Catalytic transition-metal sites in enzymes need MCPB.py bonded
        modeling on top of the 12-6-4 LJ params -- flagged in notes."""
        return "MCPB" in self.notes


_TABLE: dict[str, IonParams] | None = None


def _load_table() -> dict[str, IonParams]:
    if not _CSV_PATH.is_file():
        return {}
    out: dict[str, IonParams] = {}
    with _CSV_PATH.open(newline="", encoding="utf-8") as fh:
        for row in csv.DictReader(fh):
            resname = (row.get("pdb_resname") or "").strip().upper()
            if not resname:
                continue
            out[resname] = IonParams(
                pdb_resname=resname,
                element=(row.get("element") or "").strip(),
                formal_charge=(row.get("formal_charge") or "").strip(),
                amber_atom_type=(row.get("amber_atom_type") or "").strip(),
                ion_size_class=(row.get("ion_size_class") or "").strip(),
                li_merz_paper=(row.get("li_merz_paper") or "").strip(),
                leaprc_route=(row.get("leaprc_route") or "").strip(),
                notes=(row.get("notes") or "").strip(),
            )
    return out


def _get_table() -> dict[str, IonParams]:
    global _TABLE
    if _TABLE is None:
        _TABLE = _load_table()
    return _TABLE


def _reset_for_tests() -> None:
    """Force reload of the table on next lookup.  Testing only."""
    global _TABLE
    _TABLE = None


def lookup_ion(pdb_resname: str) -> IonParams | None:
    """Return IonParams for a PDB resname, or None if unknown.

    ``lookup_ion('ZN')`` -> Zn2+ params (Li-Merz 2013 12-6-4)
    ``lookup_ion('FE2')`` -> Fe2+ params
    ``lookup_ion('FE')``  -> Fe3+ params (PDB default when unlabelled)
    """
    if not pdb_resname:
        return None
    return _get_table().get(pdb_resname.strip().upper())


def covered_by_12_6_4(pdb_resname: str) -> bool:
    """True iff this resname has a 12-6-4 LJ parametrization."""
    ip = lookup_ion(pdb_resname)
    if ip is None:
        return False
    return "12-6-4" in ip.leaprc_route


def tleap_source_lines(pdb_resnames: list[str]) -> list[str]:
    """Return the deduplicated list of ``source leaprc...`` lines
    required to cover the given ion resnames.  Empty when no known
    ions are in the input.

    Example:
        >>> tleap_source_lines(["ZN", "MG", "GD"])
        ['source leaprc.water.tip3p_ion12-6-4',
         'source leaprc.water.tip3p_ion12-6-4_HFE']

    Order-preserving so tleap loads the base water first, then the
    HFE-tuned lanthanide extension (if needed).
    """
    seen: list[str] = []
    for res in pdb_resnames:
        ip = lookup_ion(res)
        if ip is None or not ip.leaprc_route:
            continue
        line = f"source {ip.leaprc_route}"
        if line not in seen:
            seen.append(line)
    return seen


def flag_mcpb_required(pdb_resnames: list[str]) -> list[str]:
    """Return the subset of input resnames whose IonParams.needs_mcpb is
    True -- these want a bonded-model MCPB.py cluster on top of LJ params.
    """
    out: list[str] = []
    for res in pdb_resnames:
        ip = lookup_ion(res)
        if ip is not None and ip.needs_mcpb:
            out.append(res)
    return out


def describe(pdb_resname: str) -> str:
    """Reviewer-friendly one-line summary."""
    ip = lookup_ion(pdb_resname)
    if ip is None:
        return f"{pdb_resname!r} not in ion 12-6-4 reference table"
    return (
        f"{ip.pdb_resname} ({ip.element}{ip.formal_charge}, "
        f"class={ip.ion_size_class}, amber={ip.amber_atom_type}, "
        f"leaprc={ip.leaprc_route}, ref={ip.li_merz_paper})"
    )


def all_known_ions() -> tuple[str, ...]:
    return tuple(sorted(_get_table().keys()))
