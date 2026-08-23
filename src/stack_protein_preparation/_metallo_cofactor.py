"""Metallo-cofactor detection: heme systems + iron-sulfur clusters.

USER MANDATE 2026-08-24 (deep research completed):
    heme + Fe-S clusters brauchen SPEZIALBEHANDLUNG weil Fe + porphyrin
    (bzw. Fe + inorganische sulfides) EINE parametrisierungs-einheit
    sind, nicht separate metal + cofactor.  Research report at
    ``data/metallo_cofactor_frcmods.json`` documents every route.

Design (post-deep-research):

    Heme — route by AXIAL-LIGAND FINGERPRINT, cite published frcmod:
        bis-His  → Autenrieth 2004 hemall (HIGH)     — cyt b, cyt b5
        His-Met  → Autenrieth 2004 hemall (HIGH)     — cyt c soluble
        His-only → Giammona 1984 (HIGH)              — deoxy Mb/Hb
        His-O2   → hemall O2-bound (MEDIUM)          — oxy Mb/Hb
        HEC covalent (CXXCH) → hemall + auto CYM patch for 2 Cys thioethers
        Cys-thiolate → Shahrokh 2012 (MEDIUM)        — P450 family
        HEA / heme a  → Johansson 2016 (MEDIUM)      — cyt c oxidase
        no-axial      → LOW ("orphan heme?")
        Bundle URLs and citations in
        ``data/metallo_cofactor_frcmods.json``.

    Fe-S — route by CLUSTER TYPE + REDOX + Rieske detection:
        SF4/F4S ([4Fe-4S])  → Carvalho & Swart 2014 (MEDIUM)  — per redox
        F3S ([3Fe-4S])       → Carvalho & Swart 2014 (MEDIUM)  — per redox
        FES ([2Fe-2S]) Cys4 → Carvalho & Swart 2014 (MEDIUM)  — ferredoxin
        FES ([2Fe-2S]) Rieske → Molina-Molina 2014 (MEDIUM)   — 2Cys+2His
        CFN (FeMo-co)        → LOW (no reliable static frcmod, QM/MM only)
        FEO                  → Li & Merz 12-6-4 (MEDIUM)      — mononuclear

Auto-emitted protonation overrides:
    HID for axial His in heme systems
    CYM for every bridging Cys in Fe-S clusters
    CYM for the 2 CXXCH thioether Cys residues in HEC (heme c)

xtb-based de-novo RESP is NOT USED for these systems — the research
consensus (Bannwarth 2019, Carvalho 2014, Björnsson 2022) is that
GFN2-xTB is spin-restricted open-shell and cannot represent the
broken-symmetry antiferromagnetic coupling that Fe-S / heme need.
The library-lookup route is what production pipelines actually do.
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable

from stack_protein_preparation._component_confidence import (
    ComponentConfidence,
    Confidence,
)


HEME_RESNAMES = ("HEM", "HEB", "HEC", "HEA", "HEO", "HEZ", "SRM", "COH")
FES_CLUSTER_RESNAMES = {
    "SF4": ("[4Fe-4S]", 4),
    "F4S": ("[4Fe-4S]", 4),
    "F3S": ("[3Fe-4S]", 3),
    "FES": ("[2Fe-2S]", 2),
    "FE2": ("[2Fe]",    2),   # rare, di-iron non-sulfur cluster
    "FEO": ("[Fe-O]",   1),
    "CFN": ("[MoFe-7S]", None),  # nitrogenase FeMo-co (special)
    "ICS": ("[4Fe-4S]", 4),
}

# Ligand distances for classifying heme axial contacts (protein atom → Fe).
# Fe-N(His) typical 1.9-2.2 Å; Fe-S(Met/Cys) 2.2-2.5 Å; Fe-O(water/Tyr) 1.9-2.2 Å.
_AXIAL_CUTOFF_A: dict[str, float] = {
    "N": 2.6, "O": 2.6, "S": 2.9,
}


@dataclass
class HemeSystem:
    heme_resname: str        # HEM/HEB/HEC/…
    heme_chain: str
    heme_resnum: int
    fe_coord: tuple[float, float, float]
    axial_ligands: list[tuple[str, str, int, str, float]]
    # (protein_resname, protein_chain, protein_resnum, protein_atom, distance)
    is_covalent_hec: bool = False  # True if HEC + 2 Cys-thioether links found

    def axial_type(self) -> str:
        """Return a short label of the axial coordination pattern."""
        elements = sorted(
            (l[3][0] if l[3] else "?") for l in self.axial_ligands
        )
        by_res = sorted(l[0] for l in self.axial_ligands)
        if by_res == ["HIS", "HIS"]:
            return "bis-His"          # cyt b5, cyt c-like
        if by_res in (["HIS", "MET"], ["MET", "HIS"]):
            return "His-Met"          # cyt c, cyt c1
        if by_res == ["HIS"]:
            return "His-only"         # deoxy-hemoglobin / myoglobin
        if by_res == ["HIS", "SER"] or by_res == ["HIS", "TYR"]:
            return "His-hydroxy"
        if by_res in (["CYS"], ["CYS", "HIS"]):
            return "Cys-thiolate"     # P450 / chloroperoxidase
        if not by_res:
            return "no-axial"
        return "-".join(by_res)

    def to_component_confidence(self) -> ComponentConfidence:
        axial = self.axial_type()
        n_axial = len(self.axial_ligands)
        name = (
            f"{self.heme_resname} {self.heme_chain}{self.heme_resnum} "
            f"({axial}, {n_axial} axial)"
        )
        # Route + confidence per axial pattern, citing the specific published
        # frcmod that FRUTON expects the user to supply (see
        # data/metallo_cofactor_frcmods.json for URLs).
        if self.is_covalent_hec:
            conf = Confidence.MEDIUM
            method = "hemall_hec_cxxch_patch"
            reason = (
                "heme c with 2 Cys-thioether covalent links (CXXCH motif); "
                "route: hemall base + Autenrieth 2004 cyt-c patch adding "
                "CAB-Sγ + CAC-Sγ bond terms; both Cys become CYM-analog"
            )
            action = (
                "supply hemall.frcmod + apply the Cys-thioether patch; "
                "verify both bridging Cys are auto-set to CYM"
            )
        elif self.heme_resname in ("HEA", "HAS", "HEO"):
            conf = Confidence.MEDIUM
            method = "johansson_2016_heme_a"
            reason = (
                "heme a / a3 variant → cytochrome c oxidase family; "
                "route: Johansson 2016 CcO parameters (PMC4979044) "
                "with CuB coupling"
            )
            action = (
                "supply Johansson 2016 heme-a + CuA/CuB parameter bundle; "
                "verify redox state of the coupled Cu centre"
            )
        elif axial == "bis-His":
            conf = Confidence.HIGH
            method = "autenrieth_2004_hemall_bis_his"
            reason = (
                "bis-His axial → cytochrome b / b5 / neuroglobin family; "
                "Autenrieth 2004 hemall (JCC 25:1613) applies directly; "
                "pick oxidised (Fe³⁺) or reduced (Fe²⁺) charge set at build"
            )
            action = ""
        elif axial == "His-Met":
            conf = Confidence.HIGH
            method = "autenrieth_2004_hemall_his_met"
            reason = (
                "His-Met axial → cytochrome c soluble; "
                "Autenrieth 2004 hemall His-Met variant applies"
            )
            action = ""
        elif axial == "His-only":
            conf = Confidence.HIGH
            method = "giammona_1984_deoxy"
            reason = (
                "single-His axial → deoxy-myoglobin / deoxy-hemoglobin; "
                "Giammona 1984 5-coord Fe(II) HS or hemall 5-coord block"
            )
            action = ""
        elif "Cys-thiolate" in axial or (axial == "CYS" or "CYS" in axial):
            conf = Confidence.MEDIUM
            method = "shahrokh_2012_p450"
            reason = (
                "Cys-thiolate axial → cytochrome P450 family; "
                "Shahrokh 2012 (JCC 33:119, PMC3242737) 4-state bundle: "
                "pick Fe³⁺ resting / Fe²⁺ reduced / Fe²⁺-O2 oxy / Cpd I"
            )
            action = (
                "supply Shahrokh 2012 P450 frcmod (SI download); "
                "set FRUTON_HEME_STATE=[resting|reduced|oxy|cpdI] "
                "for the reactive state; axial Cys auto-set to CYM"
            )
        elif axial == "no-axial":
            conf = Confidence.LOW
            method = "orphan_heme"
            reason = (
                "no axial protein ligand within 2.6 Å of Fe — orphan heme "
                "or crystallographic-only occupancy?"
            )
            action = (
                "verify the heme is in a real binding site (not a soaked-in "
                "isolated cofactor); may need manual re-placement"
            )
        else:
            conf = Confidence.MEDIUM
            method = "hemall_generic"
            reason = (
                f"non-canonical axial coordination pattern {axial!r} — "
                f"hemall may fit but has not been validated for this topology"
            )
            action = (
                "verify axial-ligand parametrisation matches paper; "
                "consider a per-protein custom MCPB.py run"
            )

        return ComponentConfidence(
            component_type="cofactor",
            name=name,
            confidence=conf,
            reason=reason,
            suggested_action=action,
            method=method,
            details={
                "resname": self.heme_resname,
                "axial_type": axial,
                "n_axial": n_axial,
                "covalent_hec": self.is_covalent_hec,
                "axial_ligands": [
                    {
                        "resname": r, "chain": c, "resnum": rn,
                        "atom": at, "distance_A": round(d, 3),
                    }
                    for (r, c, rn, at, d) in self.axial_ligands
                ],
            },
        )

    def his_overrides(self) -> list[tuple[str, int, str, str]]:
        """Return (chain, resnum, icode, forced_resname) for axial His → HID."""
        out = []
        for r, c, rn, at, _d in self.axial_ligands:
            if r == "HIS":
                # Axial His donates its Nε (usually) to Fe → HID (Nδ-H, Nε free/to metal)
                if at == "NE2":
                    out.append((c, rn, "", "HID"))
                elif at == "ND1":
                    out.append((c, rn, "", "HIE"))
        return out


@dataclass
class FeSCluster:
    resname: str            # SF4 / F3S / FES / ...
    cluster_label: str      # "[4Fe-4S]" etc.
    chain: str
    resnum: int
    n_fe_expected: int | None
    bridging_cys_residues: list[tuple[str, int, str]] = field(default_factory=list)
    # (chain, resnum, icode) for each Cys within thiolate-bonding distance
    rieske_his_residues: list[tuple[str, int, str]] = field(default_factory=list)
    # (chain, resnum, icode) for each His within Fe-N cutoff (Rieske variant)

    def is_rieske(self) -> bool:
        """Rieske [2Fe-2S] has 2 Cys + 2 His anchors (not 4 Cys)."""
        return (
            self.resname == "FES"
            and len(self.bridging_cys_residues) == 2
            and len(self.rieske_his_residues) == 2
        )

    def to_component_confidence(self) -> ComponentConfidence:
        name = f"{self.resname} {self.chain}{self.resnum} ({self.cluster_label})"
        n_cys = len(self.bridging_cys_residues)

        # FeMo-co / FeFe-co nitrogenase: no reliable static frcmod exists.
        # Björnsson group workflow (PMC8908755) requires ChemShell QM/MM.
        if self.resname == "CFN":
            return ComponentConfidence(
                component_type="cofactor",
                name=name,
                confidence=Confidence.LOW,
                reason=(
                    "FeMo-cofactor of nitrogenase: no reliable static MM "
                    "parameter set exists.  Björnsson group ASH/ChemShell "
                    "workflow (PMC8908755) requires QM/MM, not pure MM."
                ),
                suggested_action=(
                    "route this protein to QM/MM (ASH — "
                    "https://ash.readthedocs.io/en/latest/QM-MM-protein.html); "
                    "OR exclude from the pure-MM benchmark"
                ),
                method="femoco_qmmm_only",
                details={
                    "resname": self.resname, "cluster": self.cluster_label,
                    "reason_category": "manual_review_required",
                },
            )

        if n_cys == 0 and not self.rieske_his_residues:
            return ComponentConfidence(
                component_type="cofactor",
                name=name,
                confidence=Confidence.MEDIUM,
                reason=(
                    f"{self.cluster_label} cluster detected but no bridging "
                    f"Cys/His anchors within cutoff — orphan cluster?"
                ),
                suggested_action=(
                    "verify the cluster is protein-bound in the paper; "
                    "if orphan, may need manual placement or removal"
                ),
                method="fes_orphan_cluster",
                details={"resname": self.resname},
            )

        # Rieske [2Fe-2S] (2 Cys + 2 His) → Molina-Molina 2014 custom set
        if self.is_rieske():
            return ComponentConfidence(
                component_type="cofactor",
                name=name + " [Rieske 2Cys+2His]",
                confidence=Confidence.MEDIUM,
                reason=(
                    "Rieske [2Fe-2S] variant: 2 Cys + 2 His anchors "
                    "(cf. cytochrome b6f, nitrobenzene dioxygenase).  "
                    "Route: Molina-Molina 2014 custom set "
                    "(JCTC 10:2, doi:10.1021/ct500205z)"
                ),
                suggested_action=(
                    "supply Molina-Molina 2014 Rieske frcmod; both Cys "
                    "auto-set to CYM; His axial partonation depends on "
                    "which N points at Fe (auto-detected by metal-coord "
                    "protonation module)"
                ),
                method="molina_molina_2014_rieske",
                details={
                    "resname": self.resname, "cluster": self.cluster_label,
                    "n_bridging_cys": n_cys,
                    "n_rieske_his": len(self.rieske_his_residues),
                },
            )

        # Standard [4Fe-4S] / [3Fe-4S] / [2Fe-2S]-Cys4 → Carvalho & Swart 2014.
        # Per-redox-state frcmod: FRUTON emits with default redox_state
        # '2+' unless the user supplies FRUTON_FES_REDOX env var.
        return ComponentConfidence(
            component_type="cofactor",
            name=name,
            confidence=Confidence.MEDIUM,
            reason=(
                f"{self.cluster_label} cluster with {n_cys} bridging Cys "
                f"thiolates; route: Carvalho & Swart 2014 "
                f"(JCIM 54:613, doi:10.1021/ci400718m).  Bonded model "
                f"with Seminario constants from broken-symmetry DFT; "
                f"per-redox-state frcmod (default 2+ for [4Fe-4S], "
                f"1+ for [3Fe-4S])"
            ),
            suggested_action=(
                "supply Carvalho & Swart 2014 frcmod for the appropriate "
                "cluster + redox state; verify total charge in tleap output "
                "(e.g. [4Fe-4S]²⁺ + 4×CYM(-1) = −2 net); all bridging "
                "Cys auto-set to CYM"
            ),
            method="carvalho_swart_2014",
            details={
                "resname": self.resname, "cluster": self.cluster_label,
                "n_fe_expected": self.n_fe_expected,
                "n_bridging_cys": n_cys,
                "default_redox_state": "2+" if self.resname == "SF4" else "1+",
                "bridging_cys": [
                    {"chain": c, "resnum": rn, "icode": ic}
                    for (c, rn, ic) in self.bridging_cys_residues
                ],
            },
        )
        return ComponentConfidence(
            component_type="cofactor",
            name=name,
            confidence=conf,
            reason=reason,
            suggested_action=action,
            method=method,
            details={
                "resname": self.resname,
                "cluster": self.cluster_label,
                "n_fe_expected": self.n_fe_expected,
                "n_bridging_cys": n_cys,
                "bridging_cys": [
                    {"chain": c, "resnum": rn, "icode": ic}
                    for (c, rn, ic) in self.bridging_cys_residues
                ],
            },
        )

    def cys_overrides(self) -> list[tuple[str, int, str, str]]:
        """Every bridging Cys → CYM (thiolate)."""
        return [(c, rn, ic, "CYM") for (c, rn, ic) in self.bridging_cys_residues]


def _read_hetatm_atoms(pdb_path: Path) -> list[tuple[str, str, int, str, str, tuple[float, float, float]]]:
    """Yield (resname, chain, resnum, icode, atom_name, xyz) for every HETATM.

    Standard PDB columns: resname 17-20, chain 21, resnum 22-26, icode 26.
    Atom name 12-16, coords 30-54 in 8.3f triplets.
    """
    out = []
    if not pdb_path.is_file():
        return out
    for line in pdb_path.read_text(errors="replace").splitlines():
        if not line.startswith("HETATM"):
            continue
        try:
            resname = line[17:20].strip().upper()
            chain = line[21].strip() or "?"
            resnum = int(line[22:26])
            icode = line[26].strip()
            atom = line[12:16].strip().upper()
            x = float(line[30:38]); y = float(line[38:46]); z = float(line[46:54])
        except (ValueError, IndexError):
            continue
        out.append((resname, chain, resnum, icode, atom, (x, y, z)))
    return out


def _read_atom_atoms(pdb_path: Path) -> list[tuple[str, str, int, str, str, tuple[float, float, float]]]:
    """Same as _read_hetatm_atoms but for ATOM (protein) records."""
    out = []
    if not pdb_path.is_file():
        return out
    for line in pdb_path.read_text(errors="replace").splitlines():
        if not line.startswith("ATOM"):
            continue
        try:
            resname = line[17:20].strip().upper()
            chain = line[21].strip() or "?"
            resnum = int(line[22:26])
            icode = line[26].strip()
            atom = line[12:16].strip().upper()
            x = float(line[30:38]); y = float(line[38:46]); z = float(line[46:54])
        except (ValueError, IndexError):
            continue
        out.append((resname, chain, resnum, icode, atom, (x, y, z)))
    return out


def _distance(a: tuple[float, float, float], b: tuple[float, float, float]) -> float:
    return math.sqrt(sum((a[i] - b[i]) ** 2 for i in range(3)))


def detect_heme_systems(pdb_path: str | Path) -> list[HemeSystem]:
    """Find every HEM/HEA/HEB/HEC residue + resolve its axial ligands.

    Two-step:
      1. Scan HETATM for heme residues; locate the FE atom inside each.
      2. Neighborsearch protein ATOM records within axial cutoff of Fe.
    """
    pdb_path = Path(pdb_path)
    hetatm = _read_hetatm_atoms(pdb_path)
    atom = _read_atom_atoms(pdb_path)
    if not hetatm:
        return []

    # Group heme HETATM atoms per residue
    by_heme: dict[tuple[str, str, int, str], list] = {}
    for (rn, c, num, ic, at, xyz) in hetatm:
        if rn in HEME_RESNAMES:
            by_heme.setdefault((rn, c, num, ic), []).append((at, xyz))

    hemes: list[HemeSystem] = []
    for (rn, c, num, ic), atoms in by_heme.items():
        fe_coord = None
        for at, xyz in atoms:
            if at == "FE":
                fe_coord = xyz
                break
        if fe_coord is None:
            continue           # heme without Fe? skip (unusual)

        # Axial ligands: protein atoms within cutoff
        axial: list = []
        for (a_rn, a_c, a_num, a_ic, a_at, a_xyz) in atom:
            elem = a_at[0].upper() if a_at else ""
            cutoff = _AXIAL_CUTOFF_A.get(elem)
            if not cutoff:
                continue
            # Only consider likely donor atoms
            if a_at not in ("NE2", "ND1", "SD", "SG", "OG", "OG1", "OH", "OD1", "OD2", "OE1", "OE2"):
                continue
            d = _distance(fe_coord, a_xyz)
            if d <= cutoff:
                axial.append((a_rn, a_c, a_num, a_at, d))

        # Sort by distance (closest first) — usually first 1-2 are the axials
        axial.sort(key=lambda x: x[4])
        # Truncate to the 2 closest for canonical heme (bis-axial max)
        axial = axial[:2]

        # HEC covalent detection: two Cys thioethers link to heme via CAB/CAC
        is_hec_covalent = False
        if rn == "HEC":
            # For HEC, expect 2 Cys SG atoms within 2.5 Å of heme CAB or CAC
            cab = next((xyz for at, xyz in atoms if at == "CAB"), None)
            cac = next((xyz for at, xyz in atoms if at == "CAC"), None)
            n_cys_thioether = 0
            for (a_rn, a_c, a_num, a_ic, a_at, a_xyz) in atom:
                if a_rn == "CYS" and a_at == "SG":
                    for target in (cab, cac):
                        if target and _distance(a_xyz, target) < 2.5:
                            n_cys_thioether += 1
                            break
            is_hec_covalent = (n_cys_thioether >= 2)

        hemes.append(HemeSystem(
            heme_resname=rn, heme_chain=c, heme_resnum=num,
            fe_coord=fe_coord, axial_ligands=axial,
            is_covalent_hec=is_hec_covalent,
        ))
    return hemes


def detect_fes_clusters(pdb_path: str | Path) -> list[FeSCluster]:
    """Find every [nFe-mS] cluster HETATM residue + its bridging Cys."""
    pdb_path = Path(pdb_path)
    hetatm = _read_hetatm_atoms(pdb_path)
    atom = _read_atom_atoms(pdb_path)
    if not hetatm:
        return []

    # Group cluster HETATM atoms per residue
    by_cluster: dict[tuple[str, str, int, str], list] = {}
    for (rn, c, num, ic, at, xyz) in hetatm:
        if rn in FES_CLUSTER_RESNAMES:
            by_cluster.setdefault((rn, c, num, ic), []).append((at, xyz))

    clusters: list[FeSCluster] = []
    for (rn, c, num, ic), atoms in by_cluster.items():
        label, n_fe = FES_CLUSTER_RESNAMES[rn]
        # Bridging Cys: Cys-SG within 2.6 Å of any Fe atom in the cluster
        fe_positions = [xyz for at, xyz in atoms if at.startswith("FE")]
        bridging_cys: set[tuple[str, int, str]] = set()
        # Rieske variant of [2Fe-2S]: 2 His anchors (Nε or Nδ) within 2.6 Å
        # of Fe.  Detect any His-Nring within cutoff of any Fe.
        rieske_his: set[tuple[str, int, str]] = set()
        for (a_rn, a_c, a_num, a_ic, a_at, a_xyz) in atom:
            if a_rn == "CYS" and a_at == "SG":
                for fe in fe_positions:
                    if _distance(a_xyz, fe) < 2.6:
                        bridging_cys.add((a_c, a_num, a_ic))
                        break
            elif a_rn == "HIS" and a_at in ("NE2", "ND1"):
                for fe in fe_positions:
                    if _distance(a_xyz, fe) < 2.6:
                        rieske_his.add((a_c, a_num, a_ic))
                        break
        clusters.append(FeSCluster(
            resname=rn, cluster_label=label,
            chain=c, resnum=num,
            n_fe_expected=n_fe,
            bridging_cys_residues=sorted(bridging_cys),
            rieske_his_residues=sorted(rieske_his),
        ))
    return clusters
