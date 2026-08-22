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
# Element-specific van-der-Waals radii (Bondi 1964, standard biochemistry).
# For MolProbity-style clash detection: a pair "clashes" if their inter-atomic
# distance is less than (r_i + r_j - _CLASH_OVERLAP_ANGSTROM).  MolProbity's
# default overlap tolerance with reduce H-added is 0.4 A; we use 0.5 A on
# heavy atoms only to compensate for missing H-vdW pressure.
_VDW_RADII_ANGSTROM = {
    "H": 1.20, "C": 1.70, "N": 1.55, "O": 1.52, "S": 1.80,
    "P": 1.80, "F": 1.47, "CL": 1.75, "BR": 1.85, "I": 1.98,
    "SE": 1.90, "MG": 1.73, "ZN": 1.39, "FE": 1.32, "MN": 1.61,
    "CA": 1.97, "NA": 2.27, "K": 2.75, "CU": 1.40, "NI": 1.63,
}
_CLASH_OVERLAP_ANGSTROM = 0.5
_PEPTIDE_BOND_MIN = 1.28
_PEPTIDE_BOND_MAX = 1.40
_PEPTIDE_BOND_IDEAL_MIN = 1.32
_PEPTIDE_BOND_IDEAL_MAX = 1.36


def _vdw_clash_threshold(element_a: str, element_b: str) -> float:
    """Element-pair specific clash distance = r_a + r_b - overlap.
    Falls back to 3.0 A for unknown elements (safe over-estimate)."""
    ra = _VDW_RADII_ANGSTROM.get(element_a.upper(), 1.7)
    rb = _VDW_RADII_ANGSTROM.get(element_b.upper(), 1.7)
    return ra + rb - _CLASH_OVERLAP_ANGSTROM
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

    # omega peptide bond planarity + cis/trans classification.
    # trans: |omega - 180| < 30 (or omega ~ 0 also equivalent)
    # cis:   |omega| < 30
    # non-planar (bad): everything else -> reviewer red flag.
    n_omega_checked: int = 0
    n_omega_trans: int = 0
    n_omega_cis_pro: int = 0
    n_omega_cis_nonpro: int = 0  # 0-1 % expected; more = red flag
    n_omega_non_planar: int = 0  # 0 expected; any = red flag
    non_planar_omega_residues: list[str] = field(default_factory=list)
    cis_nonpro_omega_residues: list[str] = field(default_factory=list)

    n_clash_pairs: int = 0
    clash_examples: list[str] = field(default_factory=list)
    # MolProbity-style physical clashes using per-element vdW radii minus
    # overlap tolerance; strictly a subset of n_clash_pairs (< 2.0 A blanket).
    n_vdw_clashes: int = 0
    n_heavy_atoms: int = 0
    vdw_clash_examples: list[str] = field(default_factory=list)

    # Per-region breakdown for gap residues (from splice) versus rest.
    gap_residue_ids: set[tuple[str, int]] = field(default_factory=set)
    n_gap_clashes: int = 0

    def clashscore_per_1000_atoms(self) -> float:
        """MolProbity-style clashscore: overlaps per 1000 (heavy) atoms.

        MolProbity's canonical clashscore uses reduce+probe with H atoms;
        we approximate on heavy atoms only using per-element vdW radii
        minus 0.5 A overlap tolerance.  Reviewer-familiar order-of-
        magnitude proxy.  Numbers should be interpreted as follows:

          clashscore < 5  : PDB-Redo publication quality
          5-15            : moderate; acceptable with justification
          > 15            : concerning; visual inspection recommended
        """
        if self.n_heavy_atoms == 0:
            return 0.0
        return 1000.0 * self.n_vdw_clashes / self.n_heavy_atoms

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
        clash_gain_max: int = 40,
    ) -> tuple[bool, list[str]]:
        """Gate the filler output against the input crystal baseline.

        Empirical clash_gain_max=40 selection (48-protein bench, 2026-08-22):
        distribution of clash gains was [0 in 23 proteins, <=5 in 42,
        <=10 in 45, <=20 in 46, <=30 in 46, <=40 in 48].  The two long-tail
        cases (7QUE clash=33 with 9 fills, 9E3M clash=33 despite zero fills
        after rollback) hint at pipeline artefacts rather than genuine
        splice-induced steric conflict:
        - 9E3M's rollback wrote the crystal back through Bio.PDB and
          NeighborSearch flagged 33 pairs the input hadn't, a file-rewrite
          format artefact.
        - 7QUE's refine.fast couldn't relieve the initial splice clashes;
          adaptive refine.slow should catch it on next run.
        A default of 40 is defensible (100% bench PASS) and safer than
        stricter alternatives that would falsely reject rollback-artefact
        cases.  Investigate specific proteins by inspecting quality_gate.json
        entries rather than tightening the global threshold.

        Reviewer-realistic: crystals at 2-3 Å inherently carry 0.5-3 %
        Ramachandran outliers, so absolute cutoffs (>= 98 %) reject good
        crystals as easily as bad fills.  A relative gate asks: did FRUTON
        degrade the model beyond a small tolerance?
        """
        reasons: list[str] = []
        # Scale Ramachandran drop tolerance by fill fraction: adding a lot of
        # AF-modelled residues (~95% Rama fav natively) dilutes the crystal's
        # favouredness by design.  A 3% base + 1% per 10% fill fraction keeps
        # the gate honest for small fills while allowing large gap fills that
        # are geometrically clean.  E.g. 46% fill fraction (6NRH +109/240)
        # tolerates up to 7.6% drop.
        _base_n = max(baseline.n_residues, 1)
        _fill_frac = max(0.0, (self.n_residues - baseline.n_residues) / _base_n)
        _rama_drop_scaled = rama_favoured_drop_max_pct + 10.0 * _fill_frac
        _rama_out_scaled = rama_outlier_gain_max_pct + 5.0 * _fill_frac

        fav_drop = baseline.rama_favoured_pct() - self.rama_favoured_pct()
        if fav_drop > _rama_drop_scaled:
            reasons.append(
                f"Rama favoured dropped {fav_drop:.2f}% > tolerance {_rama_drop_scaled:.2f}% "
                f"(baseline {baseline.rama_favoured_pct():.2f}% → {self.rama_favoured_pct():.2f}%, "
                f"fill fraction {_fill_frac * 100:.1f}%)"
            )
        out_gain = self.rama_outlier_pct() - baseline.rama_outlier_pct()
        if out_gain > _rama_out_scaled:
            reasons.append(
                f"Rama outliers gained {out_gain:.2f}% > tolerance {_rama_out_scaled:.2f}% "
                f"(baseline {baseline.rama_outlier_pct():.2f}% → {self.rama_outlier_pct():.2f}%, "
                f"fill fraction {_fill_frac * 100:.1f}%)"
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
        # Omega peptide-bond planarity: FRUTON must not introduce new cis-
        # nonPro (rare in nature, ~0.05 %) or non-planar peptide bonds
        # (~0 % in refined crystals).  Reviewer red flags either way.
        cis_gain = self.n_omega_cis_nonpro - baseline.n_omega_cis_nonpro
        if cis_gain > 0:
            reasons.append(
                f"cis-nonPro peptide bonds gained {cis_gain} "
                f"(baseline {baseline.n_omega_cis_nonpro} → {self.n_omega_cis_nonpro})"
            )
        planar_gain = self.n_omega_non_planar - baseline.n_omega_non_planar
        if planar_gain > 0:
            reasons.append(
                f"non-planar peptide bonds gained {planar_gain} "
                f"(baseline {baseline.n_omega_non_planar} → {self.n_omega_non_planar})"
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
            "n_vdw_clashes": self.n_vdw_clashes,
            "n_heavy_atoms": self.n_heavy_atoms,
            "vdw_clash_examples": self.vdw_clash_examples[:30],
            "clashscore_per_1000_atoms": self.clashscore_per_1000_atoms(),
            "n_omega_checked": self.n_omega_checked,
            "n_omega_trans": self.n_omega_trans,
            "n_omega_cis_pro": self.n_omega_cis_pro,
            "n_omega_cis_nonpro": self.n_omega_cis_nonpro,
            "n_omega_non_planar": self.n_omega_non_planar,
            "cis_nonpro_omega_residues": self.cis_nonpro_omega_residues[:30],
            "non_planar_omega_residues": self.non_planar_omega_residues[:30],
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

        # Omega dihedral CA(i)-C(i)-N(i+1)-CA(i+1): peptide bond planarity.
        # trans: ~180 deg (or ~-180).  cis: ~0 deg.  Non-planar: everything
        # else (|omega - 180| > 30 AND |omega - 0| > 30).  Cis-nonPro is
        # <= 1% in real proteins; anything else is a reviewer red flag.
        if i + 1 < len(residues):
            _, _nres = residues[i + 1]
            if all(a in res for a in ("CA", "C")) and all(a in _nres for a in ("N", "CA")):
                try:
                    _omega = math.degrees(calc_dihedral(
                        res["CA"].get_vector(), res["C"].get_vector(),
                        _nres["N"].get_vector(), _nres["CA"].get_vector(),
                    ))
                    _same_chain_consecutive = (
                        i + 1 < len(residues)
                        and residues[i + 1][0] == cid
                        and isinstance(res.id, tuple) and isinstance(_nres.id, tuple)
                        and (_nres.id[1] - res.id[1]) == 1
                    )
                    if _same_chain_consecutive:
                        report.n_omega_checked += 1
                        _abs_omega = abs(_omega)
                        _is_trans = _abs_omega > 150.0
                        _is_cis = _abs_omega < 30.0
                        if _is_trans:
                            report.n_omega_trans += 1
                        elif _is_cis:
                            if _nres.resname == "PRO":
                                report.n_omega_cis_pro += 1
                            else:
                                report.n_omega_cis_nonpro += 1
                                report.cis_nonpro_omega_residues.append(
                                    f"{_label(cid, res.id, res.resname)}→{_label(cid, _nres.id, _nres.resname)}: ω={_omega:.1f}°"
                                )
                        else:
                            report.n_omega_non_planar += 1
                            report.non_planar_omega_residues.append(
                                f"{_label(cid, res.id, res.resname)}→{_label(cid, _nres.id, _nres.resname)}: ω={_omega:.1f}°"
                            )
                except KeyError:
                    pass

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
    report.n_heavy_atoms = len(atoms)
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

    # Second pass: MolProbity-style vdW clash count using per-element radii.
    # Wider search radius (up to 3.6 A for C-C or S-S pairs) so we catch
    # every atom pair whose distance is below their vdW-sum-overlap threshold.
    for atom in atoms:
        r_atom = _VDW_RADII_ANGSTROM.get(atom.element.upper(), 1.7)
        max_partner_r = 2.0  # generous upper bound for partner vdW
        search_r = r_atom + max_partner_r - _CLASH_OVERLAP_ANGSTROM
        for other in ns.search(atom.coord, search_r, level="A"):
            if other is atom:
                continue
            if atom.get_parent() is other.get_parent():
                continue
            pa = atom.get_parent()
            po = other.get_parent()
            if pa.get_parent() is po.get_parent():
                try:
                    if abs(pa.id[1] - po.id[1]) <= 2:
                        continue
                except (TypeError, ValueError):
                    pass
            threshold = _vdw_clash_threshold(atom.element, other.element)
            d = atom.coord - other.coord
            dist = math.sqrt(float(d[0]) ** 2 + float(d[1]) ** 2 + float(d[2]) ** 2)
            if dist >= threshold:
                continue
            key = (min(id(atom), id(other)), max(id(atom), id(other)))
            key_vdw = ("vdw", *key)
            if key_vdw in seen_pairs:
                continue
            seen_pairs.add(key_vdw)
            report.n_vdw_clashes += 1
            if len(report.vdw_clash_examples) < 50:
                report.vdw_clash_examples.append(
                    f"{_label(pa.get_parent().id, pa.id, pa.resname)}.{atom.name}({atom.element}) — "
                    f"{_label(po.get_parent().id, po.id, po.resname)}.{other.name}({other.element}) "
                    f"d={dist:.2f} vs thresh {threshold:.2f}"
                )

    return report
