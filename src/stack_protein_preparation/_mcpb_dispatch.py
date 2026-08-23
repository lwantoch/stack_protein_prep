"""MCPB tier-dispatch — decides which parametrization route to take per metal.

USER MANDATE 2026-08-23: FRUTON should attempt SOMETHING for every metal
even 30-residue-gap heme-iron or Fe-S clusters, and report honest
confidence per metal so the user knows where to look.

Four tiers:

    nonbonded_1264  — Li-Merz 12-6-4 LJ params (no QM run needed).
                      Confidence HIGH.  Suitable for labile / structural
                      metals in aqueous or carboxylate coordination:
                      Na+, K+, Ca2+, Mg2+, Zn2+ with ≥4 water/COO
                      ligands, no catalytic annotation.  Reviewer-
                      accepted per Li & Merz series 2013–2020.

    xtb_bonded      — semi-empirical GFN2-xTB (Bannwarth 2019 JCTC).
                      Confidence MEDIUM.  Suitable for catalytic Zn/Cu,
                      structural Fe(II)/Fe(III), when paper indicates
                      coordination bond character.  Free (GPL), fast
                      (seconds-minutes).

    dft_bonded      — full DFT via ORCA (free) or Gaussian (commercial).
                      Confidence HIGH.  Suitable for spin-crossover Fe,
                      multi-metal clusters (Fe-S), radical mechanisms,
                      unusual oxidation states, or user-opt-in
                      (FRUTON_MCPB_ENGINE=orca/gaussian).

    manual          — FRUTON cannot automate this metal.  Confidence LOW.
                      User must supply frcmod externally.

The dispatch decision is data-driven from:
    1. The metal reference oracle (`ts_metal_reference.csv`, force_field_route field).
    2. Optional paper-evidence hints ("catalytic", "spin-crossover", "Fe-S").
    3. Optional user override (FRUTON_MCPB_ENGINE env var).

License-free (stdlib only).  Reviewer-defensible: every decision comes
with a reason string that names the reference data source.
"""
from __future__ import annotations

import os
from dataclasses import dataclass, field
from typing import Any

from stack_protein_preparation._component_confidence import (
    ComponentConfidence,
    Confidence,
)
from stack_protein_preparation._metal_reference import (
    MetalReference,
    lookup_metal,
)


TIER_NONBONDED = "nonbonded_1264"
TIER_XTB_BONDED = "xtb_bonded"
TIER_DFT_BONDED = "dft_bonded"
TIER_MANUAL = "manual"


# Keywords in paper_evidence.md that push a metal from nonbonded to bonded.
# Nothing here is bench-tuned; the terms are standard biochem vocabulary.
_BONDED_KEYWORDS = (
    "catalytic",
    "catalysis",
    "covalent bond",
    "coordination bond",
    "coord. bond",
    "cross-link",
    "cofactor",
)

# Keywords that push xtb_bonded → dft_bonded (semi-empirical inadequate).
_DFT_ONLY_KEYWORDS = (
    "spin-crossover",
    "spin crossover",
    "high-spin/low-spin",
    "iron-sulfur",
    "fe-s cluster",
    "[2fe-2s]",
    "[3fe-4s]",
    "[4fe-4s]",
    "cubane",
    "mixed-valence",
    "radical",
    "mo-nitrogenase",
    "fenitrogenase",
    "molybdopterin",
)


@dataclass
class MCPBDispatchDecision:
    """One metal's dispatch verdict."""
    metal_resname: str        # PDB residue name (ZN, MG, FE, ...)
    element: str              # element symbol
    tier: str                 # one of TIER_*
    confidence: Confidence
    method_label: str         # short label written into audit CSV
    reason: str               # why this tier
    suggested_action: str = ""
    details: dict[str, Any] = field(default_factory=dict)

    def to_component_confidence(self, name: str) -> ComponentConfidence:
        """Bridge into the pipeline-wide confidence collector."""
        return ComponentConfidence(
            component_type="metal",
            name=name,
            confidence=self.confidence,
            reason=self.reason,
            suggested_action=self.suggested_action,
            method=self.method_label,
            details={
                "tier": self.tier,
                "element": self.element,
                **self.details,
            },
        )


def _extract_engine_override() -> str | None:
    """Read $FRUTON_MCPB_ENGINE (values: 'xtb', 'orca', 'gaussian')."""
    value = os.environ.get("FRUTON_MCPB_ENGINE", "").strip().lower()
    return value if value in ("xtb", "orca", "gaussian") else None


def _mentions_any(text: str, keywords: tuple[str, ...]) -> bool:
    if not text:
        return False
    lower = text.lower()
    return any(kw in lower for kw in keywords)


def decide_mcpb_route(
    metal_resname: str,
    paper_evidence_text: str = "",
    metal_reference: MetalReference | None = None,
) -> MCPBDispatchDecision:
    """Return the parametrization dispatch for ONE metal.

    Args:
        metal_resname: PDB residue name (e.g. 'ZN', 'FE', 'MG').
        paper_evidence_text: optional free text from the protein's
            ``paper_evidence.md`` — scanned for the bonded / DFT-only
            keyword tables above.
        metal_reference: optional pre-loaded MetalReference; when None,
            looked up from the packaged 89-row oracle.
    """
    resname = metal_resname.strip().upper()
    ref = metal_reference or lookup_metal(resname)

    override = _extract_engine_override()
    if override == "gaussian" or override == "orca":
        return MCPBDispatchDecision(
            metal_resname=resname,
            element=(ref.element if ref else resname),
            tier=TIER_DFT_BONDED,
            confidence=Confidence.HIGH,
            method_label=f"dft_{override}",
            reason=f"FRUTON_MCPB_ENGINE={override} — user forced DFT",
            suggested_action="",
        )
    if override == "xtb":
        return MCPBDispatchDecision(
            metal_resname=resname,
            element=(ref.element if ref else resname),
            tier=TIER_XTB_BONDED,
            confidence=Confidence.MEDIUM,
            method_label="xtb_gfn2",
            reason="FRUTON_MCPB_ENGINE=xtb — user forced xtb bonded",
            suggested_action="verify metal-donor distances after MM minimisation",
        )

    if ref is None:
        return MCPBDispatchDecision(
            metal_resname=resname,
            element=resname,
            tier=TIER_MANUAL,
            confidence=Confidence.LOW,
            method_label="unknown_metal",
            reason=f"{resname!r} not in the 89-row metal reference oracle",
            suggested_action=(
                f"add {resname} to data/ts_metal_reference.csv with its "
                "oxidation state + coordination profile, or supply "
                "external frcmod/mol2 for this metal"
            ),
            details={"needs_reference_entry": True},
        )

    # DFT-only keywords take precedence over everything below.
    if _mentions_any(paper_evidence_text, _DFT_ONLY_KEYWORDS):
        return MCPBDispatchDecision(
            metal_resname=resname,
            element=ref.element,
            tier=TIER_DFT_BONDED,
            confidence=Confidence.MEDIUM,
            method_label="dft_required",
            reason=(
                f"paper evidence mentions DFT-only feature (spin-crossover / "
                f"Fe-S cluster / radical / mixed-valence); semi-empirical "
                f"GFN2-xTB inadequate for this chemistry"
            ),
            suggested_action=(
                "run ORCA (free-academic) or Gaussian on the extracted "
                "cluster; set FRUTON_MCPB_ENGINE=orca and re-run this stage"
            ),
            details={
                "ref_notes": ref.notes,
                "dft_reason": "paper keyword match",
            },
        )

    is_catalytic = _mentions_any(paper_evidence_text, _BONDED_KEYWORDS)
    ref_route = (ref.force_field_route or "").strip().lower()

    # Reference-oracle route wins when it explicitly says '12-6-4' AND paper
    # doesn't flag catalysis.
    if "12-6-4" in ref_route and not is_catalytic:
        return MCPBDispatchDecision(
            metal_resname=resname,
            element=ref.element,
            tier=TIER_NONBONDED,
            confidence=Confidence.HIGH,
            method_label="LiMerz_12_6_4",
            reason=(
                f"reference oracle route='{ref.force_field_route}' + no "
                f"catalytic annotation in paper evidence → Li-Merz 12-6-4 "
                f"LJ nonbonded suffices"
            ),
            suggested_action="",
            details={"reference_route": ref.force_field_route},
        )

    # If the reference oracle explicitly requires MCPB, go xtb by default;
    # elevate to DFT if paper marks it reviewer-critical.
    if "mcpb" in ref_route or is_catalytic:
        return MCPBDispatchDecision(
            metal_resname=resname,
            element=ref.element,
            tier=TIER_XTB_BONDED,
            confidence=Confidence.MEDIUM,
            method_label="xtb_gfn2",
            reason=(
                f"reference route='{ref.force_field_route}'"
                + (", paper flags catalytic role" if is_catalytic else "")
                + " → bonded model; xtb GFN2 (fast, free) as first pass"
            ),
            suggested_action=(
                "check the emitted frcmod against the crystal metal-donor "
                "distances (Harding 2001 tolerances ±0.3 Å); escalate to "
                "DFT (FRUTON_MCPB_ENGINE=orca) if geometry drifts"
            ),
            details={"reference_route": ref.force_field_route},
        )

    # Default: no explicit route + not catalytic → conservative nonbonded.
    return MCPBDispatchDecision(
        metal_resname=resname,
        element=ref.element,
        tier=TIER_NONBONDED,
        confidence=Confidence.MEDIUM,
        method_label="LiMerz_12_6_4_default",
        reason=(
            f"no explicit reference route + no catalytic paper flag; "
            f"defaulting to Li-Merz 12-6-4 LJ (conservative)"
        ),
        suggested_action=(
            "verify Li-Merz parameters exist for this metal in "
            "data/ion_12_6_4_reference.csv; if not, escalate to xtb"
        ),
        details={"reference_route": ref.force_field_route or "(unspecified)"},
    )
