"""ω dihedral peptide-bond scanner (JACS R1 reviewer concern).

Standalone helper — independent of _filler_quality_check — that walks a
PDB and returns the ω dihedral (CA(i)-C(i)-N(i+1)-CA(i+1)) for every
consecutive same-chain residue pair.  Feeds the matplotlib
distribution figure at scripts/plot_omega_planarity_distribution.py.

Reviewer perspective:
- trans (|ω| > 150°) is expected for ~99% of residues
- cis (|ω| < 30°) is expected only for Pro (~5% of Pro residues in
  crystal structures per MacArthur & Thornton 1991)
- cis-nonPro > 1% or any non-planar (30° < |ω| < 150°) is a red flag
  that a modelling tool has broken peptide-bond geometry

Pure Python + Bio.PDB.  No external binaries.  License-free.
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable


TRANS_ABS_MIN_DEG = 150.0
CIS_ABS_MAX_DEG = 30.0


@dataclass
class OmegaEntry:
    chain: str
    resnum_i: int
    resname_i: str
    resnum_j: int
    resname_j: str
    omega_deg: float
    kind: str  # "trans" / "cis_pro" / "cis_nonpro" / "non_planar"


@dataclass
class OmegaScanResult:
    ran: bool
    passed: bool
    fallback_reason: str | None = None
    entries: list[OmegaEntry] = field(default_factory=list)

    def n_trans(self) -> int:
        return sum(1 for e in self.entries if e.kind == "trans")

    def n_cis_pro(self) -> int:
        return sum(1 for e in self.entries if e.kind == "cis_pro")

    def n_cis_nonpro(self) -> int:
        return sum(1 for e in self.entries if e.kind == "cis_nonpro")

    def n_non_planar(self) -> int:
        return sum(1 for e in self.entries if e.kind == "non_planar")

    def omega_values(self) -> list[float]:
        return [e.omega_deg for e in self.entries]

    def to_dict(self) -> dict:
        return {
            "ran": self.ran,
            "passed": self.passed,
            "fallback_reason": self.fallback_reason,
            "n_total": len(self.entries),
            "n_trans": self.n_trans(),
            "n_cis_pro": self.n_cis_pro(),
            "n_cis_nonpro": self.n_cis_nonpro(),
            "n_non_planar": self.n_non_planar(),
        }


def _classify(omega_deg: float, next_resname: str) -> str:
    abs_omega = abs(omega_deg)
    if abs_omega > TRANS_ABS_MIN_DEG:
        return "trans"
    if abs_omega < CIS_ABS_MAX_DEG:
        return "cis_pro" if next_resname == "PRO" else "cis_nonpro"
    return "non_planar"


def scan_omega_dihedrals(pdb_path: str | Path) -> OmegaScanResult:
    """Compute ω dihedrals for every consecutive same-chain residue pair.

    Returns fail-open result (ran=False) when the PDB is missing or
    Bio.PDB is not installed.  Otherwise walks the first model, groups
    by chain, and emits one OmegaEntry per adjacent (i, i+1) pair whose
    residue numbers differ by exactly 1.
    """
    pdb_path = Path(pdb_path)
    if not pdb_path.is_file():
        return OmegaScanResult(
            ran=False, passed=False,
            fallback_reason=f"pdb not found: {pdb_path}",
        )

    try:
        from Bio.PDB import PDBParser
        from Bio.PDB.vectors import calc_dihedral
    except Exception as exc:  # noqa: BLE001
        return OmegaScanResult(
            ran=False, passed=False,
            fallback_reason=f"Bio.PDB import failed: {exc!r}",
        )

    try:
        struct = PDBParser(QUIET=True).get_structure("m", str(pdb_path))
        model = next(iter(struct))
    except Exception as exc:  # noqa: BLE001
        return OmegaScanResult(
            ran=True, passed=False,
            fallback_reason=f"PDB parse failed: {exc!r}",
        )

    entries: list[OmegaEntry] = []
    for chain in model:
        cid = chain.id
        # Keep only standard residues (het_flag == ' ') to skip HOH / ligands.
        residues = [r for r in chain if r.id[0] == " "]
        for i in range(len(residues) - 1):
            r_i, r_j = residues[i], residues[i + 1]
            if r_j.id[1] - r_i.id[1] != 1:
                continue
            if not all(a in r_i for a in ("CA", "C")):
                continue
            if not all(a in r_j for a in ("N", "CA")):
                continue
            try:
                omega = math.degrees(calc_dihedral(
                    r_i["CA"].get_vector(),
                    r_i["C"].get_vector(),
                    r_j["N"].get_vector(),
                    r_j["CA"].get_vector(),
                ))
            except Exception:  # noqa: BLE001
                continue
            entries.append(OmegaEntry(
                chain=cid,
                resnum_i=r_i.id[1], resname_i=r_i.resname,
                resnum_j=r_j.id[1], resname_j=r_j.resname,
                omega_deg=omega,
                kind=_classify(omega, r_j.resname),
            ))

    return OmegaScanResult(
        ran=True,
        passed=len(entries) > 0,
        entries=entries,
    )


def summarise(result: OmegaScanResult) -> str:
    """Reviewer-friendly one-liner."""
    n = len(result.entries)
    if n == 0:
        return "ω scan: no residue pairs"
    return (
        f"ω scan: {n} pairs ("
        f"trans {result.n_trans()} {100 * result.n_trans() / n:.1f}%, "
        f"cis-Pro {result.n_cis_pro()}, "
        f"cis-nonPro {result.n_cis_nonpro()}, "
        f"non-planar {result.n_non_planar()})"
    )


def aggregate_omega_values(results: Iterable[OmegaScanResult]) -> list[float]:
    """Concatenate ω values across a bench-wide iterable of scans."""
    out: list[float] = []
    for r in results:
        if r.ran and r.passed:
            out.extend(r.omega_values())
    return out
