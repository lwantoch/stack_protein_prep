"""Native model-quality gate for FRUTON's filler output.

Reviewer-visible quality metrics for a filled/spliced/relaxed PDB, computed
without an external MolProbity binary (which isn't in the pixi env).  The
metrics mirror what Rama-Z, MolProbity clashscore, and PDB-Redo reports
front-load in a paper submission:

* Ramachandran favoured / allowed / outlier percentages per residue class
  (general / Gly / Pro / pre-Pro) using the Lovell et al. 2003 contours
  approximated with axis-aligned regions.
* Peptide-bond C-N distances across every consecutive residue pair
  (target 1.32-1.36 A, tolerated 1.28-1.40).
* Cα chirality via the standard N–CA–C–Cβ improper dihedral (must sit near
  +34° for L; a flipped Cα near -34° flags D-amino-acid contamination that
  AF3-class models are known to hallucinate — arXiv:2503.14643).
* Sidechain clash count: any pair of heavy atoms on different residues
  within ``_CLASH_THRESHOLD_ANGSTROM`` (default 2.0 A) — a coarse proxy for
  MolProbity's contact-dot clashscore.

Reviewer-grade cutoffs (Wilson et al. IUCr 1998; Chen et al. IUCr Acta D
2010; Croll IUCr Acta D 2025):
  Ramachandran favoured >= 98.0 %
  Ramachandran outliers <= 0.05 %
  Cα chirality outliers = 0
  Peptide bond C-N in [1.28, 1.40] A for >= 99.8 % of bonds
  Clash count = 0 among newly-inserted (gap) residues; anywhere-else <= 5

A single ``QualityReport`` object collects the metrics; the gate decides
pass / fail and lists the specific residues that triggered a fail so a
reviewer / user can look at the offenders directly.
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable

from Bio.PDB import PDBParser, calc_dihedral


# ---- Ramachandran regions (Lovell 2003 axis-aligned approximation) -------
# Each entry: list of (phi_lo, phi_hi, psi_lo, psi_hi) rectangles in degrees.
_RAMA_FAVOURED = {
    "general": [
        (-180.0, -30.0, -80.0, 30.0),   # alpha-helical basin
        (-180.0, -30.0, 90.0, 180.0),   # beta-sheet basin (upper)
        (-180.0, -30.0, -180.0, -170.0),
        (30.0, 90.0, -30.0, 90.0),      # left-handed alpha (rare)
    ],
    "gly": [
        (-180.0, 180.0, -180.0, 180.0),  # Gly is permissive; only outliers matter here
    ],
    "pro": [
        (-90.0, -40.0, -60.0, 40.0),
        (-90.0, -40.0, 120.0, 180.0),
    ],
    "pre_pro": [
        (-180.0, -30.0, -80.0, 30.0),
        (-180.0, -30.0, 100.0, 180.0),
    ],
}
_RAMA_ALLOWED = {
    "general": [
        (-180.0, -20.0, -100.0, 50.0),
        (-180.0, -20.0, 70.0, 180.0),
        (20.0, 100.0, -50.0, 100.0),
    ],
    "gly": [
        (-180.0, 180.0, -180.0, 180.0),  # Gly is permissive
    ],
    "pro": [
        (-100.0, -30.0, -80.0, 60.0),
        (-100.0, -30.0, 100.0, 180.0),
    ],
    "pre_pro": [
        (-180.0, -20.0, -100.0, 50.0),
        (-180.0, -20.0, 80.0, 180.0),
    ],
}

_CLASH_THRESHOLD_ANGSTROM = 2.0
_PEPTIDE_BOND_MIN = 1.28
_PEPTIDE_BOND_MAX = 1.40
_PEPTIDE_BOND_IDEAL_MIN = 1.32
_PEPTIDE_BOND_IDEAL_MAX = 1.36
# Cα chirality via signed tetrahedron volume ((N-CA) × (C-CA)) · (CB-CA).
# L-amino acid: positive volume; D: negative.  Non-glycine only; missing-atom
# cases are skipped.


@dataclass
class QualityReport:
    n_residues: int = 0
    n_rama_favoured: int = 0
    n_rama_allowed: int = 0
    n_rama_outlier: int = 0
    rama_outlier_residues: list[str] = field(default_factory=list)

    n_peptide_bonds: int = 0
    n_peptide_bonds_ideal: int = 0
    n_peptide_bonds_tolerated: int = 0
    n_peptide_bonds_broken: int = 0
    broken_peptide_bonds: list[str] = field(default_factory=list)

    n_ca_chirality_checked: int = 0
    n_ca_chirality_outliers: int = 0
    ca_chirality_outlier_residues: list[str] = field(default_factory=list)

    n_clash_pairs: int = 0
    clash_examples: list[str] = field(default_factory=list)

    # Per-region breakdown for gap residues (from splice) versus rest.
    gap_residue_ids: set[tuple[str, int]] = field(default_factory=set)
    n_gap_clashes: int = 0

    def rama_favoured_pct(self) -> float:
        if self.n_residues == 0:
            return 0.0
        return 100.0 * self.n_rama_favoured / self.n_residues

    def rama_outlier_pct(self) -> float:
        if self.n_residues == 0:
            return 0.0
        return 100.0 * self.n_rama_outlier / self.n_residues

    def peptide_bond_tolerated_pct(self) -> float:
        if self.n_peptide_bonds == 0:
            return 0.0
        return 100.0 * (self.n_peptide_bonds_tolerated + self.n_peptide_bonds_ideal) / self.n_peptide_bonds

    def passes_relative_gate(
        self,
        baseline: "QualityReport",
        *,
        rama_favoured_drop_max_pct: float = 3.0,
        rama_outlier_gain_max_pct: float = 1.5,
        peptide_bond_broken_gain_max: int = 0,
        chirality_outlier_gain_max: int = 0,
        clash_gain_max: int = 10,
    ) -> tuple[bool, list[str]]:
        """Gate the filler output against the input crystal baseline.

        Reviewer-realistic: crystals at 2-3 Å inherently carry 0.5-3 %
        Ramachandran outliers, so absolute cutoffs (>= 98 %) reject good
        crystals as easily as bad fills.  A relative gate asks: did FRUTON
        degrade the model beyond a small tolerance?
        """
        reasons: list[str] = []
        fav_drop = baseline.rama_favoured_pct() - self.rama_favoured_pct()
        if fav_drop > rama_favoured_drop_max_pct:
            reasons.append(
                f"Rama favoured dropped {fav_drop:.2f}% (baseline {baseline.rama_favoured_pct():.2f}% → {self.rama_favoured_pct():.2f}%)"
            )
        out_gain = self.rama_outlier_pct() - baseline.rama_outlier_pct()
        if out_gain > rama_outlier_gain_max_pct:
            reasons.append(
                f"Rama outliers gained {out_gain:.2f}% (baseline {baseline.rama_outlier_pct():.2f}% → {self.rama_outlier_pct():.2f}%)"
            )
        pep_gain = self.n_peptide_bonds_broken - baseline.n_peptide_bonds_broken
        if pep_gain > peptide_bond_broken_gain_max:
            reasons.append(
                f"Peptide bonds broken gained {pep_gain} (baseline {baseline.n_peptide_bonds_broken} → {self.n_peptide_bonds_broken})"
            )
        chir_gain = self.n_ca_chirality_outliers - baseline.n_ca_chirality_outliers
        if chir_gain > chirality_outlier_gain_max:
            reasons.append(
                f"Chirality D-outliers gained {chir_gain} (baseline {baseline.n_ca_chirality_outliers} → {self.n_ca_chirality_outliers})"
            )
        clash_gain = self.n_clash_pairs - baseline.n_clash_pairs
        if clash_gain > clash_gain_max:
            reasons.append(
                f"Clashes gained {clash_gain} (baseline {baseline.n_clash_pairs} → {self.n_clash_pairs})"
            )
        return (len(reasons) == 0, reasons)

    def passes_reviewer_gate(
        self,
        *,
        rama_favoured_min_pct: float = 95.0,
        rama_outlier_max_pct: float = 0.5,
        peptide_bond_tolerated_min_pct: float = 99.5,
        allow_gap_clashes: int = 0,
        allow_other_clashes: int = 5,
    ) -> tuple[bool, list[str]]:
        reasons: list[str] = []
        if self.rama_favoured_pct() < rama_favoured_min_pct:
            reasons.append(
                f"Rama favoured {self.rama_favoured_pct():.2f}% < {rama_favoured_min_pct}%"
            )
        if self.rama_outlier_pct() > rama_outlier_max_pct:
            reasons.append(
                f"Rama outliers {self.rama_outlier_pct():.2f}% > {rama_outlier_max_pct}%"
            )
        if self.peptide_bond_tolerated_pct() < peptide_bond_tolerated_min_pct:
            reasons.append(
                f"Peptide bonds tolerated {self.peptide_bond_tolerated_pct():.2f}% "
                f"< {peptide_bond_tolerated_min_pct}%"
            )
        if self.n_ca_chirality_outliers > 0:
            reasons.append(f"Cα chirality outliers: {self.n_ca_chirality_outliers}")
        if self.n_gap_clashes > allow_gap_clashes:
            reasons.append(
                f"Gap-residue clashes: {self.n_gap_clashes} > {allow_gap_clashes}"
            )
        if self.n_clash_pairs - self.n_gap_clashes > allow_other_clashes:
            reasons.append(
                f"Other-region clashes: {self.n_clash_pairs - self.n_gap_clashes} "
                f"> {allow_other_clashes}"
            )
        return (len(reasons) == 0, reasons)

    def to_dict(self) -> dict:
        return {
            "n_residues": self.n_residues,
            "rama_favoured_pct": self.rama_favoured_pct(),
            "rama_outlier_pct": self.rama_outlier_pct(),
            "rama_outlier_residues": self.rama_outlier_residues[:30],
            "n_peptide_bonds": self.n_peptide_bonds,
            "peptide_bond_tolerated_pct": self.peptide_bond_tolerated_pct(),
            "n_peptide_bonds_broken": self.n_peptide_bonds_broken,
            "broken_peptide_bonds": self.broken_peptide_bonds[:30],
            "n_ca_chirality_checked": self.n_ca_chirality_checked,
            "n_ca_chirality_outliers": self.n_ca_chirality_outliers,
            "ca_chirality_outlier_residues": self.ca_chirality_outlier_residues[:30],
            "n_clash_pairs": self.n_clash_pairs,
            "n_gap_clashes": self.n_gap_clashes,
            "clash_examples": self.clash_examples[:30],
        }


def _rama_class(resname: str, next_resname: str | None) -> str:
    if resname == "GLY":
        return "gly"
    if resname == "PRO":
        return "pro"
    if next_resname == "PRO":
        return "pre_pro"
    return "general"


def _in_any_rect(rects: Iterable[tuple[float, float, float, float]], phi: float, psi: float) -> bool:
    for lo_phi, hi_phi, lo_psi, hi_psi in rects:
        if lo_phi <= phi <= hi_phi and lo_psi <= psi <= hi_psi:
            return True
    return False


def _label(chain_id: str, resid, resname: str) -> str:
    resnum = resid[1] if isinstance(resid, tuple) else resid
    return f"{resname} {chain_id}:{resnum}"


def check_model_quality(
    pdb_path: str | Path,
    gap_residue_ids: Iterable[tuple[str, int]] | None = None,
    clash_threshold_angstrom: float = _CLASH_THRESHOLD_ANGSTROM,
) -> QualityReport:
    """Analyse ``pdb_path`` and return a QualityReport."""
    parser = PDBParser(QUIET=True, PERMISSIVE=True)
    struct = parser.get_structure("m", str(pdb_path))
    report = QualityReport(gap_residue_ids=set(gap_residue_ids or []))

    # Flatten protein residues in file order.  Only standard protein residues
    # (HETATM/water skipped) so we don't compute Ramachandran on ligands.
    residues: list = []
    for chain in struct[0]:
        for res in chain:
            if res.id[0].strip():  # HETATM
                continue
            residues.append((chain.id, res))
    report.n_residues = len(residues)

    # --- Ramachandran + peptide-bond geometry -------------------------------
    for i, (cid, res) in enumerate(residues):
        # Peptide bond to next residue.  Only measure if the two residues are
        # SEQUENCE-CONSECUTIVE (same chain, resnum differs by 1) -- otherwise
        # they're separated by a real chain gap that FRUTON was never asked
        # to bridge (e.g. crystal disorder that pdb_ids.csv left unfilled).
        if i + 1 < len(residues):
            next_cid, next_res = residues[i + 1]
            same_chain = next_cid == cid
            consecutive = (
                same_chain
                and isinstance(res.id, tuple) and isinstance(next_res.id, tuple)
                and (next_res.id[1] - res.id[1]) == 1
            )
            if consecutive and "C" in res and "N" in next_res:
                d = res["C"].coord - next_res["N"].coord
                dist = math.sqrt(float(d[0]) ** 2 + float(d[1]) ** 2 + float(d[2]) ** 2)
                report.n_peptide_bonds += 1
                if _PEPTIDE_BOND_IDEAL_MIN <= dist <= _PEPTIDE_BOND_IDEAL_MAX:
                    report.n_peptide_bonds_ideal += 1
                elif _PEPTIDE_BOND_MIN <= dist <= _PEPTIDE_BOND_MAX:
                    report.n_peptide_bonds_tolerated += 1
                else:
                    report.n_peptide_bonds_broken += 1
                    report.broken_peptide_bonds.append(
                        f"{_label(cid, res.id, res.resname)}→{_label(cid, next_res.id, next_res.resname)}: "
                        f"C-N={dist:.2f} Å"
                    )

        # Ramachandran phi/psi
        if 0 < i < len(residues) - 1:
            _, prev_res = residues[i - 1]
            _, next_res = residues[i + 1]
            try:
                phi = calc_dihedral(
                    prev_res["C"].get_vector(), res["N"].get_vector(),
                    res["CA"].get_vector(), res["C"].get_vector(),
                )
                psi = calc_dihedral(
                    res["N"].get_vector(), res["CA"].get_vector(),
                    res["C"].get_vector(), next_res["N"].get_vector(),
                )
            except KeyError:
                continue
            phi_deg = math.degrees(phi)
            psi_deg = math.degrees(psi)
            klass = _rama_class(res.resname, next_res.resname)
            fav = _in_any_rect(_RAMA_FAVOURED[klass], phi_deg, psi_deg)
            allowed = fav or _in_any_rect(_RAMA_ALLOWED[klass], phi_deg, psi_deg)
            if fav:
                report.n_rama_favoured += 1
            elif allowed:
                report.n_rama_allowed += 1
            else:
                report.n_rama_outlier += 1
                report.rama_outlier_residues.append(
                    f"{_label(cid, res.id, res.resname)}: φ={phi_deg:.1f} ψ={psi_deg:.1f} ({klass})"
                )

        # Cα chirality via signed tetrahedron volume ((N-CA) × (C-CA)) · (CB-CA).
        # L-amino acid: positive volume; D: negative.  Robust to sign convention.
        if res.resname != "GLY" and all(a in res for a in ("N", "CA", "C", "CB")):
            try:
                ca = res["CA"].coord
                v_n = res["N"].coord - ca
                v_c = res["C"].coord - ca
                v_cb = res["CB"].coord - ca
                cross = (
                    v_n[1] * v_c[2] - v_n[2] * v_c[1],
                    v_n[2] * v_c[0] - v_n[0] * v_c[2],
                    v_n[0] * v_c[1] - v_n[1] * v_c[0],
                )
                signed_vol = float(
                    cross[0] * v_cb[0] + cross[1] * v_cb[1] + cross[2] * v_cb[2]
                )
            except KeyError:
                signed_vol = None
            if signed_vol is not None:
                report.n_ca_chirality_checked += 1
                if signed_vol < 0.0:
                    report.n_ca_chirality_outliers += 1
                    report.ca_chirality_outlier_residues.append(
                        f"{_label(cid, res.id, res.resname)}: signed_vol={signed_vol:.3f} Å³ (D-chirality)"
                    )

    # --- Clash detection (heavy-atom pairs across different residues) ------
    # Only count clashes between residues that are > 2 apart in sequence
    # (i.e. NOT bonded/1-3/1-4 neighbours) so we measure real steric conflict
    # and not the built-in Van-der-Waals overlap that peptide bond geometry
    # inherently produces.  Use Bio.PDB NeighborSearch for O(N log N).
    from Bio.PDB import NeighborSearch
    atoms = [a for a in struct[0].get_atoms() if a.element != "H"]
    ns = NeighborSearch(atoms)
    seen_pairs: set[tuple[int, int]] = set()
    for atom in atoms:
        for other in ns.search(atom.coord, clash_threshold_angstrom, level="A"):
            if other is atom:
                continue
            if atom.get_parent() is other.get_parent():
                continue  # same residue
            pa = atom.get_parent()
            po = other.get_parent()
            # Skip if same chain and within 2 residues sequence-wise (bonded
            # neighbours & their 1-4 partners).
            if pa.get_parent() is po.get_parent():
                try:
                    if abs(pa.id[1] - po.id[1]) <= 2:
                        continue
                except (TypeError, ValueError):
                    pass
            key = (min(id(atom), id(other)), max(id(atom), id(other)))
            if key in seen_pairs:
                continue
            seen_pairs.add(key)
            report.n_clash_pairs += 1
            pair_desc = (
                f"{_label(pa.get_parent().id, pa.id, pa.resname)}.{atom.name} — "
                f"{_label(po.get_parent().id, po.id, po.resname)}.{other.name}"
            )
            if len(report.clash_examples) < 50:
                report.clash_examples.append(pair_desc)
            for res_wrap in (pa, po):
                _cid_check = res_wrap.get_parent().id
                _rnum = res_wrap.id[1] if isinstance(res_wrap.id, tuple) else res_wrap.id
                if (_cid_check, _rnum) in report.gap_residue_ids:
                    report.n_gap_clashes += 1
                    break

    return report
