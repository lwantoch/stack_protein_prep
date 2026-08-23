"""DSSP secondary-structure computation + IDR cross-check.

Reviewer concern (JACS R4): a MobiDB-annotated intrinsically-disordered
region tells us the *sequence* is IDR-prone, but the SPECIFIC crystal
in FRUTON's input may have resolved a structured local conformation
(unusual ligand contact, stabilised by a partner chain, snapshot of a
transient fold).  The IDR gate as currently wired (via
_uniprot_idr.gap_overlaps_uniprot_idr) rejects the AF fill for such a
window even though the crystal itself carries observed structured
residues flanking the gap.

This module lets FRUTON compute the DSSP 3-state (H/E/C) secondary-
structure assignment on the crystal template and cross-check against
the MobiDB IDR annotation before deciding whether to keep the gap
fill.  Two rules a JACS reviewer would accept:

1. If crystal residues flanking the gap (± 3 residues) are H or E in
   DSSP, the gap probably folds too -- overriding the MobiDB IDR
   rejection is defensible.
2. If crystal residues flanking are all-C (coil) *and* MobiDB says
   IDR, the rejection stands.

Fail-open pattern (matches _uniprot_idr, _cofactor_params): when the
DSSP binary is not on PATH, the cross-check returns "no advisory"
so the pipeline degrades gracefully to the pure-sequence IDR gate.

Uses Bio.PDB.DSSP as the wrapper; requires the ``mkdssp`` (or
legacy ``dssp``) binary to be installed and on PATH.  Free (BSD-3)
per Zielesny/Kabsch original + Touw 2015 update.  License-free per
Davide's later refactoring pass.
"""
from __future__ import annotations

import shutil
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable


# DSSP 8-state -> 3-state mapping (Frishman & Argos 1995 canonical):
#   H, G, I  -> H (helix)
#   E, B     -> E (sheet)
#   T, S, -  -> C (coil / turn / bend / unresolved)
_DSSP8_TO_SS3 = {
    "H": "H", "G": "H", "I": "H",
    "E": "E", "B": "E",
    "T": "C", "S": "C", "-": "C", "": "C", " ": "C",
    "P": "C",  # polyproline II
}


@dataclass
class DSSPResult:
    ran: bool
    passed: bool
    fallback_reason: str | None = None
    # {(chain_id, resnum): ss3_char}
    ss3_by_residue: dict[tuple[str, int], str] = field(default_factory=dict)
    diagnostics: list[str] = field(default_factory=list)

    def n_helix(self) -> int:
        return sum(1 for v in self.ss3_by_residue.values() if v == "H")

    def n_sheet(self) -> int:
        return sum(1 for v in self.ss3_by_residue.values() if v == "E")

    def n_coil(self) -> int:
        return sum(1 for v in self.ss3_by_residue.values() if v == "C")

    def to_dict(self) -> dict:
        return {
            "ran": self.ran,
            "passed": self.passed,
            "fallback_reason": self.fallback_reason,
            "n_residues": len(self.ss3_by_residue),
            "n_helix": self.n_helix(),
            "n_sheet": self.n_sheet(),
            "n_coil": self.n_coil(),
            "diagnostics": list(self.diagnostics),
        }


@dataclass
class IdrDisagreement:
    """One residue where MobiDB IDR annotation disagrees with local SS."""
    chain: str
    resnum: int
    ss3: str                    # H, E, or C
    mobidb_annotation: str      # "idr" or "structured"
    advice: str                 # "override_idr_reject" / "keep_idr_reject" / etc.


def _resolve_dssp_binary(dssp_bin: str | None) -> str | None:
    """Try mkdssp first (modern), then dssp (legacy) unless overridden."""
    if dssp_bin:
        return shutil.which(dssp_bin)
    for name in ("mkdssp", "dssp"):
        found = shutil.which(name)
        if found:
            return found
    return None


def compute_dssp_ss3(
    pdb_path: str | Path,
    dssp_bin: str | None = None,
) -> DSSPResult:
    """Run DSSP on ``pdb_path`` and return per-residue 3-state SS.

    Returns ``ran=False`` when DSSP is not installed; the caller then
    falls back to the pure-MobiDB IDR gate.  ``passed=True`` when
    DSSP ran without exception AND at least one residue was assigned.
    """
    pdb_path = Path(pdb_path)
    if not pdb_path.is_file():
        return DSSPResult(
            ran=False, passed=False,
            fallback_reason=f"pdb not found: {pdb_path}",
        )

    resolved = _resolve_dssp_binary(dssp_bin)
    if resolved is None:
        return DSSPResult(
            ran=False, passed=False,
            fallback_reason=(
                "dssp/mkdssp not on PATH; install via 'apt install dssp' "
                "or 'conda install -c salilab dssp', or set the path via "
                "the dssp_bin argument."
            ),
        )

    try:
        from Bio.PDB import PDBParser
        from Bio.PDB.DSSP import DSSP
    except Exception as exc:  # noqa: BLE001
        return DSSPResult(
            ran=False, passed=False,
            fallback_reason=f"Bio.PDB.DSSP import failed: {exc!r}",
        )

    try:
        struct = PDBParser(QUIET=True).get_structure("m", str(pdb_path))
        model = next(iter(struct))
        dssp = DSSP(model, str(pdb_path), dssp=resolved)
    except Exception as exc:  # noqa: BLE001
        return DSSPResult(
            ran=True, passed=False,
            fallback_reason=f"DSSP execution failed: {exc!r}",
        )

    ss3: dict[tuple[str, int], str] = {}
    for key in dssp.keys():
        # Bio.PDB DSSP keys are (chain_id, ('', resnum, ''))
        try:
            chain_id = key[0]
            resnum = int(key[1][1])
        except (IndexError, ValueError):
            continue
        try:
            ss8 = dssp[key][2]
        except (IndexError, KeyError):
            continue
        ss3[(chain_id, resnum)] = _DSSP8_TO_SS3.get(ss8, "C")

    return DSSPResult(
        ran=True,
        passed=len(ss3) > 0,
        ss3_by_residue=ss3,
        diagnostics=[
            f"dssp_bin={resolved}",
            f"n_residues_assigned={len(ss3)}",
        ],
    )


def cross_check_idr_vs_local_ss(
    dssp_result: DSSPResult,
    gap_ranges_by_chain: Iterable[tuple[str, int, int]],
    idr_flag_by_gap: dict[tuple[str, int, int], bool],
    flank_window: int = 3,
) -> list[IdrDisagreement]:
    """For each gap flagged as IDR by MobiDB, inspect local SS3 on the
    flanking crystal residues (± ``flank_window``).  If any flanking
    residue is H or E, emit an advisory to override the IDR reject.

    Args:
        dssp_result: output of ``compute_dssp_ss3`` on the crystal PDB.
        gap_ranges_by_chain: iterable of ``(chain_id, first_resnum,
            last_resnum)`` tuples naming AF-inserted regions.
        idr_flag_by_gap: mapping ``{(chain, first, last): is_idr}``
            per the MobiDB gate (``True`` if the gap window overlaps
            an IDR region).
        flank_window: how many residues each side of the gap to inspect.

    Returns list of ``IdrDisagreement`` entries; empty when everything
    agrees or DSSP is unavailable.
    """
    if not dssp_result.ran or not dssp_result.passed:
        return []

    out: list[IdrDisagreement] = []
    for chain, lo, hi in gap_ranges_by_chain:
        gap_key = (chain, lo, hi)
        is_idr = idr_flag_by_gap.get(gap_key, False)
        if not is_idr:
            continue
        flanking: list[tuple[int, str]] = []
        for offset in range(1, flank_window + 1):
            for target in (lo - offset, hi + offset):
                ss = dssp_result.ss3_by_residue.get((chain, target))
                if ss:
                    flanking.append((target, ss))
        structured = [(r, s) for r, s in flanking if s in ("H", "E")]
        if structured:
            # Emit one advisory per structured flanking residue -- the
            # user then decides whether to override the IDR reject.
            for resnum, ss in structured:
                out.append(IdrDisagreement(
                    chain=chain,
                    resnum=resnum,
                    ss3=ss,
                    mobidb_annotation="idr",
                    advice=(
                        f"crystal residue {chain}:{resnum} is {ss} "
                        f"(structured) but MobiDB flags gap {chain}:"
                        f"{lo}-{hi} as IDR -- consider OVERRIDING the "
                        f"IDR reject: local structural context may be "
                        f"stabilised by ligand / partner chain."
                    ),
                ))
        else:
            out.append(IdrDisagreement(
                chain=chain,
                resnum=lo,
                ss3="C",
                mobidb_annotation="idr",
                advice=(
                    f"gap {chain}:{lo}-{hi} flanked by all-coil (or missing) "
                    f"crystal residues AND MobiDB flags as IDR -- IDR "
                    f"reject stands: no structural evidence to override."
                ),
            ))
    return out


def summarise_ss3_composition(dssp_result: DSSPResult) -> str:
    """Reviewer-friendly one-liner: '212 residues (H 87 41.0%, E 45 21.2%, C 80 37.7%)'."""
    n = len(dssp_result.ss3_by_residue)
    if n == 0:
        return "DSSP: no residues assigned"
    h = dssp_result.n_helix()
    e = dssp_result.n_sheet()
    c = dssp_result.n_coil()
    return (
        f"DSSP: {n} residues ("
        f"H {h} {100 * h / n:.1f}%, "
        f"E {e} {100 * e / n:.1f}%, "
        f"C {c} {100 * c / n:.1f}%)"
    )
