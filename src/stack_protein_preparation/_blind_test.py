"""Blind-test harness: hide crystal residues, compare FRUTON fill to
hidden ground-truth (Nature R3 concern).

Reviewer perspective: FRUTON's bench numbers ("residues rescued",
"pass rate") report *self-consistency* between the FRUTON output and
the geometry gates.  A stronger claim is:

    if we hide a well-resolved loop from the crystal and force FRUTON
    to reconstruct it from AF + splice, how close does the fill land
    to the original crystal coordinates?

Two pieces here:

1. ``mask_crystal_pdb``: writes a copy of the input crystal PDB with a
   list of ``(chain, first_resnum, last_resnum)`` ranges *removed*.
   That masked PDB feeds FRUTON as if the crystallographer never
   observed those residues.

2. ``score_blind_fill``: after the user runs FRUTON on the masked
   crystal, feeds (filled_pdb, original_crystal_pdb, held_out_ranges)
   back through this helper.  Per matched (chain, resnum) residue in
   the held-out range it computes Cα distance, backbone atom distances
   (N, CA, C, O), a residue-identity match flag, plus per-range and
   overall Cα-RMSD.

Pure Python + Bio.PDB.  License-free.
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable


@dataclass(frozen=True)
class GapRange:
    """Held-out residue range: inclusive first and last residue numbers."""
    chain: str
    first_resnum: int
    last_resnum: int

    def contains(self, chain: str, resnum: int) -> bool:
        return (
            chain == self.chain
            and self.first_resnum <= resnum <= self.last_resnum
        )

    def size(self) -> int:
        return self.last_resnum - self.first_resnum + 1


@dataclass
class ResidueBlindScore:
    chain: str
    resnum: int
    crystal_resname: str
    filled_resname: str | None
    ca_distance_A: float | None  # None when filled model lacks the residue
    backbone_rmsd_A: float | None  # over available {N,CA,C,O}

    def resname_matches(self) -> bool:
        return (
            self.filled_resname is not None
            and self.filled_resname.strip().upper()
            == self.crystal_resname.strip().upper()
        )


@dataclass
class RangeBlindScore:
    range: GapRange
    residues: list[ResidueBlindScore] = field(default_factory=list)

    def n_matched(self) -> int:
        return sum(1 for r in self.residues if r.ca_distance_A is not None)

    def ca_rmsd_A(self) -> float | None:
        vals = [r.ca_distance_A for r in self.residues if r.ca_distance_A is not None]
        if not vals:
            return None
        return math.sqrt(sum(v * v for v in vals) / len(vals))

    def backbone_rmsd_A(self) -> float | None:
        vals = [r.backbone_rmsd_A for r in self.residues if r.backbone_rmsd_A is not None]
        if not vals:
            return None
        return math.sqrt(sum(v * v for v in vals) / len(vals))

    def n_resname_mismatch(self) -> int:
        return sum(1 for r in self.residues if not r.resname_matches())


@dataclass
class BlindTestResult:
    ran: bool
    passed: bool
    fallback_reason: str | None = None
    ranges: list[RangeBlindScore] = field(default_factory=list)

    def overall_ca_rmsd_A(self) -> float | None:
        all_vals: list[float] = []
        for rng in self.ranges:
            for r in rng.residues:
                if r.ca_distance_A is not None:
                    all_vals.append(r.ca_distance_A)
        if not all_vals:
            return None
        return math.sqrt(sum(v * v for v in all_vals) / len(all_vals))

    def n_matched_total(self) -> int:
        return sum(r.n_matched() for r in self.ranges)

    def to_dict(self) -> dict:
        return {
            "ran": self.ran,
            "passed": self.passed,
            "fallback_reason": self.fallback_reason,
            "n_ranges": len(self.ranges),
            "n_matched_total": self.n_matched_total(),
            "overall_ca_rmsd_A": (
                None if self.overall_ca_rmsd_A() is None
                else round(self.overall_ca_rmsd_A(), 4)
            ),
            "per_range": [
                {
                    "chain": r.range.chain,
                    "first_resnum": r.range.first_resnum,
                    "last_resnum": r.range.last_resnum,
                    "size": r.range.size(),
                    "n_matched": r.n_matched(),
                    "ca_rmsd_A": (None if r.ca_rmsd_A() is None
                                  else round(r.ca_rmsd_A(), 4)),
                    "backbone_rmsd_A": (None if r.backbone_rmsd_A() is None
                                        else round(r.backbone_rmsd_A(), 4)),
                    "n_resname_mismatch": r.n_resname_mismatch(),
                }
                for r in self.ranges
            ],
        }


def mask_crystal_pdb(
    crystal_pdb_path: str | Path,
    ranges: Iterable[GapRange],
    output_pdb_path: str | Path,
) -> Path:
    """Write ``crystal_pdb_path`` to ``output_pdb_path`` with residues
    inside ``ranges`` removed.

    Non-ATOM lines (HETATM, TER, HEADER, REMARK, END, …) are preserved.
    ATOM lines are dropped when their (chain, resnum) matches a range.
    """
    crystal_pdb_path = Path(crystal_pdb_path)
    output_pdb_path = Path(output_pdb_path)
    output_pdb_path.parent.mkdir(parents=True, exist_ok=True)
    ranges = list(ranges)

    with crystal_pdb_path.open() as fh_in, output_pdb_path.open("w") as fh_out:
        for line in fh_in:
            if line.startswith("ATOM"):
                # PDB columns 22 (chain) and 23-26 (resnum)
                try:
                    chain = line[21]
                    resnum = int(line[22:26])
                except (IndexError, ValueError):
                    fh_out.write(line)
                    continue
                if any(r.contains(chain, resnum) for r in ranges):
                    continue
            fh_out.write(line)
    return output_pdb_path


def _load_model(pdb_path: Path):
    from Bio.PDB import PDBParser
    struct = PDBParser(QUIET=True).get_structure("m", str(pdb_path))
    return next(iter(struct))


def _index_residues(model) -> dict[tuple[str, int], object]:
    out: dict[tuple[str, int], object] = {}
    for chain in model:
        cid = chain.id
        for res in chain:
            if res.id[0] != " ":  # skip het / water
                continue
            out[(cid, res.id[1])] = res
    return out


def _atom_distance(a, b) -> float:
    return float(a - b)


def _backbone_rmsd(res_a, res_b) -> float | None:
    """RMSD over backbone atoms present in both residues."""
    values: list[float] = []
    for name in ("N", "CA", "C", "O"):
        if name in res_a and name in res_b:
            d = _atom_distance(res_a[name], res_b[name])
            values.append(d * d)
    if not values:
        return None
    return math.sqrt(sum(values) / len(values))


def score_blind_fill(
    crystal_pdb_path: str | Path,
    filled_pdb_path: str | Path,
    held_out_ranges: Iterable[GapRange],
) -> BlindTestResult:
    """Score a FRUTON fill against the held-out crystal ground truth.

    Fail-open: any missing input or Bio.PDB import failure returns
    ``ran=False`` so a bench harness can skip cleanly.
    """
    crystal_pdb_path = Path(crystal_pdb_path)
    filled_pdb_path = Path(filled_pdb_path)
    if not crystal_pdb_path.is_file():
        return BlindTestResult(
            ran=False, passed=False,
            fallback_reason=f"crystal pdb not found: {crystal_pdb_path}",
        )
    if not filled_pdb_path.is_file():
        return BlindTestResult(
            ran=False, passed=False,
            fallback_reason=f"filled pdb not found: {filled_pdb_path}",
        )

    try:
        crystal_model = _load_model(crystal_pdb_path)
        filled_model = _load_model(filled_pdb_path)
    except Exception as exc:  # noqa: BLE001
        return BlindTestResult(
            ran=False, passed=False,
            fallback_reason=f"PDB parse failed: {exc!r}",
        )

    crystal_by_key = _index_residues(crystal_model)
    filled_by_key = _index_residues(filled_model)

    ranges = list(held_out_ranges)
    range_scores: list[RangeBlindScore] = []

    for rng in ranges:
        r_scores: list[ResidueBlindScore] = []
        for resnum in range(rng.first_resnum, rng.last_resnum + 1):
            key = (rng.chain, resnum)
            crystal_res = crystal_by_key.get(key)
            if crystal_res is None:
                # Crystal itself lacked the residue — nothing to score.
                continue
            filled_res = filled_by_key.get(key)
            if filled_res is None:
                r_scores.append(ResidueBlindScore(
                    chain=rng.chain, resnum=resnum,
                    crystal_resname=crystal_res.resname,
                    filled_resname=None,
                    ca_distance_A=None,
                    backbone_rmsd_A=None,
                ))
                continue
            ca_dist = (
                _atom_distance(crystal_res["CA"], filled_res["CA"])
                if "CA" in crystal_res and "CA" in filled_res
                else None
            )
            bb_rmsd = _backbone_rmsd(crystal_res, filled_res)
            r_scores.append(ResidueBlindScore(
                chain=rng.chain, resnum=resnum,
                crystal_resname=crystal_res.resname,
                filled_resname=filled_res.resname,
                ca_distance_A=ca_dist,
                backbone_rmsd_A=bb_rmsd,
            ))
        range_scores.append(RangeBlindScore(range=rng, residues=r_scores))

    return BlindTestResult(
        ran=True,
        passed=any(r.n_matched() > 0 for r in range_scores),
        ranges=range_scores,
    )


def summarise(result: BlindTestResult) -> str:
    if not result.ran:
        return f"blind-test: skipped ({result.fallback_reason})"
    if not result.passed:
        return "blind-test: ran but zero residues matched"
    rmsd = result.overall_ca_rmsd_A()
    rmsd_str = "n/a" if rmsd is None else f"{rmsd:.3f} Å"
    return (
        f"blind-test: {result.n_matched_total()} residues scored across "
        f"{len(result.ranges)} ranges, overall Cα-RMSD = {rmsd_str}"
    )
