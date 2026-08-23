"""Metal-coordinating side-chain geometry scanner (JACS R3 concern).

Reviewer perspective: FRUTON runs energy minimisation + LoopModel loop
refinement.  Reviewers worry: does that displace side chains that
coordinate metal ions?  Even a ~0.3 Å shift on a His-Nε or Asp-Oδ
donor can distort the catalytic geometry to the point of biological
irrelevance.

This module computes, per metal:
- which side-chain donor atoms are within a coordination cutoff (default
  3.0 Å per Harding 2001 / CheckMyMetal 2008 cutoffs),
- the metal↔donor distance,
- an optional crystal-vs-FRUTON Δd for the same donor identity.

Standalone: no dependency on metalls_check / metal_hydrogen_cleanup.
Uses Bio.PDB + a local mini-catalogue of standard-AA donor atoms.
License-free (BSD-3 Bio.PDB only).
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable


# Per-residue coordinating-atom list (canonical AA donors).  Kept small
# and standalone so the module can serve as a self-contained audit tool.
DONOR_ATOMS_BY_RESNAME: dict[str, set[str]] = {
    "HIS": {"ND1", "NE2"},
    "HID": {"ND1", "NE2"},
    "HIE": {"ND1", "NE2"},
    "HIP": {"ND1", "NE2"},
    "CYS": {"SG"},
    "CYX": {"SG"},
    "CYM": {"SG"},
    "SEC": {"SE"},
    "MET": {"SD"},
    "TYR": {"OH"},
    "ASP": {"OD1", "OD2"},
    "ASH": {"OD1", "OD2"},
    "GLU": {"OE1", "OE2"},
    "GLH": {"OE1", "OE2"},
    "SER": {"OG"},
    "THR": {"OG1"},
    "ASN": {"OD1"},
    "GLN": {"OE1"},
    "LYS": {"NZ"},
    "ARG": {"NH1", "NH2"},
}

# Common structural / catalytic metals detected in the PDB.  Any atom
# whose element (upper-case) is in this set is treated as a metal.
METAL_ELEMENTS: set[str] = {
    "ZN", "MG", "CA", "MN", "FE", "CU", "NI", "CO", "CD",
    "NA", "K", "LI", "CS", "RB",
    "MO", "W", "V", "PT", "PD",
    "AU", "AG", "HG", "PB",
    "LU", "GD", "EU", "TB", "DY", "HO", "ER", "TM", "YB",
}

DEFAULT_COORD_CUTOFF_A = 3.0


@dataclass(frozen=True)
class DonorKey:
    """Identity of a donor atom (matchable across crystal + FRUTON)."""
    metal_chain: str
    metal_resnum: int
    metal_element: str
    donor_chain: str
    donor_resnum: int
    donor_resname: str
    donor_atom: str


@dataclass
class MetalDonorContact:
    key: DonorKey
    distance_A: float


@dataclass
class MetalCoordScanResult:
    ran: bool
    passed: bool
    fallback_reason: str | None = None
    contacts: list[MetalDonorContact] = field(default_factory=list)
    n_metals: int = 0

    def n_contacts(self) -> int:
        return len(self.contacts)

    def by_metal_element(self) -> dict[str, int]:
        counts: dict[str, int] = {}
        for c in self.contacts:
            counts[c.key.metal_element] = counts.get(c.key.metal_element, 0) + 1
        return counts

    def to_dict(self) -> dict:
        return {
            "ran": self.ran,
            "passed": self.passed,
            "fallback_reason": self.fallback_reason,
            "n_metals": self.n_metals,
            "n_contacts": self.n_contacts(),
            "by_metal_element": self.by_metal_element(),
        }


@dataclass
class MetalDonorDelta:
    """Crystal-vs-FRUTON Δd for a matched donor identity."""
    key: DonorKey
    distance_crystal_A: float | None
    distance_fruton_A: float | None

    def delta_A(self) -> float | None:
        if self.distance_crystal_A is None or self.distance_fruton_A is None:
            return None
        return self.distance_fruton_A - self.distance_crystal_A

    def status(self) -> str:
        if self.distance_crystal_A is None and self.distance_fruton_A is not None:
            return "gained"
        if self.distance_crystal_A is not None and self.distance_fruton_A is None:
            return "lost"
        return "preserved"


def _element_of(atom) -> str:
    e = (getattr(atom, "element", "") or "").strip().upper()
    if e:
        return e
    # Fallback: infer from atom name (rare in modern PDBs).
    name = atom.get_name().strip().upper()
    return "".join(c for c in name if c.isalpha())[:2]


def scan_metal_coordination(
    pdb_path: str | Path,
    coord_cutoff_A: float = DEFAULT_COORD_CUTOFF_A,
) -> MetalCoordScanResult:
    """Return per-metal coordination contacts under ``coord_cutoff_A``.

    Fail-open: missing PDB or missing Bio.PDB → ran=False.
    Otherwise walks the first model, indexes side-chain donor atoms
    per DONOR_ATOMS_BY_RESNAME, and for each metal emits one contact
    per donor within cutoff.
    """
    pdb_path = Path(pdb_path)
    if not pdb_path.is_file():
        return MetalCoordScanResult(
            ran=False, passed=False,
            fallback_reason=f"pdb not found: {pdb_path}",
        )

    try:
        from Bio.PDB import PDBParser, NeighborSearch
    except Exception as exc:  # noqa: BLE001
        return MetalCoordScanResult(
            ran=False, passed=False,
            fallback_reason=f"Bio.PDB import failed: {exc!r}",
        )

    try:
        struct = PDBParser(QUIET=True).get_structure("m", str(pdb_path))
        model = next(iter(struct))
    except Exception as exc:  # noqa: BLE001
        return MetalCoordScanResult(
            ran=True, passed=False,
            fallback_reason=f"PDB parse failed: {exc!r}",
        )

    donor_atoms = []          # (chain_id, resnum, resname, atom_name, Bio.PDB.Atom)
    metal_atoms = []          # (chain_id, resnum, element, Bio.PDB.Atom)

    for chain in model:
        cid = chain.id
        for res in chain:
            resname = res.resname.strip().upper()
            donor_names = DONOR_ATOMS_BY_RESNAME.get(resname, set())
            for atom in res:
                aname = atom.get_name().strip().upper()
                elem = _element_of(atom)
                if elem in METAL_ELEMENTS:
                    metal_atoms.append((cid, res.id[1], elem, atom))
                    continue
                # het_flag ' ' = standard AA; skip HOH/ligand donors on this pass
                if res.id[0] != " ":
                    continue
                if aname in donor_names:
                    donor_atoms.append((cid, res.id[1], resname, aname, atom))

    if not metal_atoms:
        return MetalCoordScanResult(
            ran=True, passed=True, n_metals=0,
        )

    all_donor_bio_atoms = [d[4] for d in donor_atoms]
    ns = NeighborSearch(all_donor_bio_atoms) if all_donor_bio_atoms else None

    contacts: list[MetalDonorContact] = []
    for m_cid, m_resnum, m_elem, m_atom in metal_atoms:
        if ns is None:
            continue
        neighbours = ns.search(m_atom.coord, coord_cutoff_A, level="A")
        for donor_atom in neighbours:
            for cid, resnum, resname, aname, ref_atom in donor_atoms:
                if ref_atom is donor_atom:
                    d = float(m_atom - donor_atom)
                    contacts.append(MetalDonorContact(
                        key=DonorKey(
                            metal_chain=m_cid, metal_resnum=m_resnum,
                            metal_element=m_elem,
                            donor_chain=cid, donor_resnum=resnum,
                            donor_resname=resname, donor_atom=aname,
                        ),
                        distance_A=d,
                    ))
                    break

    return MetalCoordScanResult(
        ran=True, passed=True,
        n_metals=len(metal_atoms),
        contacts=contacts,
    )


def compare_metal_coordination(
    crystal_result: MetalCoordScanResult,
    fruton_result: MetalCoordScanResult,
) -> list[MetalDonorDelta]:
    """Match donors by DonorKey across crystal + FRUTON, emit per-donor Δd.

    A donor present in one but not the other is emitted with the missing
    distance as None (status='lost' or 'gained').
    """
    crystal_by_key = {c.key: c.distance_A for c in crystal_result.contacts}
    fruton_by_key = {c.key: c.distance_A for c in fruton_result.contacts}
    keys = set(crystal_by_key) | set(fruton_by_key)
    return [
        MetalDonorDelta(
            key=k,
            distance_crystal_A=crystal_by_key.get(k),
            distance_fruton_A=fruton_by_key.get(k),
        )
        for k in sorted(
            keys,
            key=lambda k: (
                k.metal_chain, k.metal_resnum,
                k.donor_chain, k.donor_resnum, k.donor_atom,
            ),
        )
    ]


def summarise(result: MetalCoordScanResult) -> str:
    if not result.ran:
        return f"metal-coord scan: skipped ({result.fallback_reason})"
    if result.n_metals == 0:
        return "metal-coord scan: no metal ions in structure"
    breakdown = ", ".join(
        f"{elem} {n}" for elem, n in sorted(result.by_metal_element().items())
    )
    return (
        f"metal-coord scan: {result.n_metals} metals, "
        f"{result.n_contacts()} donor contacts ({breakdown})"
    )


def aggregate_deltas(
    deltas_per_pdb: Iterable[list[MetalDonorDelta]],
) -> list[float]:
    """Flatten preserved-donor Δd values across a bench iterable."""
    out: list[float] = []
    for deltas in deltas_per_pdb:
        for d in deltas:
            delta = d.delta_A()
            if delta is not None:
                out.append(delta)
    return out
