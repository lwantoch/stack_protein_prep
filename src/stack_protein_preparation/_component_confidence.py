"""Per-component confidence data model — the reviewer-honest core of FRUTON.

USER REFRAME 2026-08-23: 'wir wollen das FRUTON für möglichst viele
Proteine ein fertiges MD-Model liefert. Egal ob 30-Residue-Lücke oder
Heme-Iron oder Iron-Schwefel-Cluster: der versucht das Beste daraus zu
machen. Ihr sollt nicht 100 % Erfolgsrate als Ziel sehen. Aber pro
Komponente Bericht liefern: bei dem Zn2+ sind wir uns sicher, bei
diesem +4 Metal haben wir zwar ein Model gebaut aber können nicht
100 % sagen es ist korrekt — vielleicht hast du Glück, vielleicht
musst du manuell eingreifen.'

Design:
    Each protein produces N component confidence records; each record
    names one thing FRUTON did (a specific gap, a specific metal ion,
    a specific cofactor, a specific protonation override, ...) and
    attaches a HIGH / MEDIUM / LOW confidence + a human-readable reason
    + a suggested action if not HIGH.

    The protein's overall status is derived from the min confidence
    across components — but importantly, the reviewer sees the WHOLE
    list, not just the min.  A protein with 5 HIGH gaps + 1 LOW metal
    is not "rejected" — it is "delivered_with_notes: check the Fe
    parametrisation".

Component types (from the pipeline stages FRUTON touches):
    gap_fill       — one per gap window; pLDDT + IDR + length gate it.
    metal          — one per metal ion; MCPB tier + reference oracle.
    cofactor       — one per HETATM cofactor; canonical vs antechamber.
    protonation    — one per active-site override; paper-confirmed vs
                     auto-suggested.
    nonstd_residue — one per non-standard AA; catalogue vs antechamber.
    disulfide      — one per SS bond; auto-detected vs paper-confirmed.

Nothing in this module is bench-specific.  Reviewers get the same
confidence rubric on any protein the pipeline is fed.  License-free
(stdlib only).
"""
from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum
from typing import Any, Iterable


class Confidence(str, Enum):
    """Reviewer-facing confidence bucket."""
    HIGH = "high"      # FRUTON is confident; standard reference route
    MEDIUM = "medium"  # FRUTON produced a plausible answer; user glance recommended
    LOW = "low"        # FRUTON attempted best-effort; user MUST verify
    FAILED = "failed"  # FRUTON could not produce this component at all

    @classmethod
    def worst_of(cls, values: Iterable["Confidence"]) -> "Confidence":
        """Aggregate: the min confidence across components dominates."""
        order = {cls.HIGH: 0, cls.MEDIUM: 1, cls.LOW: 2, cls.FAILED: 3}
        seq = list(values)
        if not seq:
            return cls.HIGH
        return max(seq, key=lambda c: order[c])


COMPONENT_TYPES = (
    "gap_fill",
    "metal",
    "cofactor",
    "protonation",
    "nonstd_residue",
    "disulfide",
)


@dataclass
class ComponentConfidence:
    """One thing FRUTON did, with an honest confidence + suggested action."""
    component_type: str            # one of COMPONENT_TYPES
    name: str                      # e.g. 'chain A gap 45-58', 'Zn A501', 'HEM A600'
    confidence: Confidence
    reason: str                    # WHY this confidence — must not be empty
    suggested_action: str = ""     # only meaningful when confidence != HIGH
    method: str = ""               # short label, e.g. 'LiMerz1264', 'MODELLER_slow', 'antechamber_gaff2'
    details: dict[str, Any] = field(default_factory=dict)

    def __post_init__(self):
        if not self.reason:
            raise ValueError(
                f"ComponentConfidence({self.component_type}, {self.name}): "
                "reason must not be empty — reviewers need to know WHY"
            )
        if self.component_type not in COMPONENT_TYPES:
            raise ValueError(
                f"unknown component_type {self.component_type!r}; "
                f"expected one of {COMPONENT_TYPES}"
            )
        if not isinstance(self.confidence, Confidence):
            # Accept string too, for JSON round-trips
            self.confidence = Confidence(str(self.confidence))

    def to_dict(self) -> dict[str, Any]:
        return {
            "component_type": self.component_type,
            "name": self.name,
            "confidence": self.confidence.value,
            "method": self.method,
            "reason": self.reason,
            "suggested_action": self.suggested_action,
            "details": dict(self.details),
        }


@dataclass
class ProteinDeliveryReport:
    """The reviewer-facing summary for one protein.

    `overall_status` categorises what the user gets:
        delivered_full_confidence   — every component HIGH; ship it
        delivered_with_notes        — at least one MEDIUM; user glance
        delivered_needs_review      — at least one LOW; user MUST verify
        not_delivered               — pipeline could not build a loadable
                                      model (rare — e.g. splice + all
                                      refine retries failed catastrophically)
    """
    pdb: str
    components: list[ComponentConfidence] = field(default_factory=list)
    model_written: bool = True    # False only when the pipeline crashed hard
    tleap_loads: bool | None = None  # populated by sanity-validate step
    md_deck_written: bool = False
    notes: str = ""

    @property
    def overall_confidence(self) -> Confidence:
        return Confidence.worst_of(c.confidence for c in self.components)

    @property
    def overall_status(self) -> str:
        if not self.model_written:
            return "not_delivered"
        worst = self.overall_confidence
        if worst is Confidence.FAILED:
            # A component failed but we still produced a partial model
            # → still delivered, but with a component the user must fix.
            return "delivered_needs_review"
        if worst is Confidence.LOW:
            return "delivered_needs_review"
        if worst is Confidence.MEDIUM:
            return "delivered_with_notes"
        return "delivered_full_confidence"

    def components_by_confidence(self, level: Confidence) -> list[ComponentConfidence]:
        return [c for c in self.components if c.confidence is level]

    def component_type_counts(self) -> dict[str, dict[str, int]]:
        """{component_type: {confidence: count}}."""
        out: dict[str, dict[str, int]] = {}
        for c in self.components:
            bucket = out.setdefault(c.component_type, {
                lv.value: 0 for lv in Confidence
            })
            bucket[c.confidence.value] += 1
        return out

    def action_items(self) -> list[str]:
        """Return the list of suggested user actions for non-HIGH components."""
        return [
            f"[{c.component_type} {c.name}] {c.suggested_action}"
            for c in self.components
            if c.confidence is not Confidence.HIGH and c.suggested_action
        ]

    def to_dict(self) -> dict[str, Any]:
        return {
            "pdb": self.pdb,
            "overall_status": self.overall_status,
            "overall_confidence": self.overall_confidence.value,
            "model_written": self.model_written,
            "tleap_loads": self.tleap_loads,
            "md_deck_written": self.md_deck_written,
            "notes": self.notes,
            "n_components": len(self.components),
            "confidence_breakdown": self.component_type_counts(),
            "components": [c.to_dict() for c in self.components],
            "action_items": self.action_items(),
        }


def summarise_delivery(rows: Iterable[ProteinDeliveryReport]) -> dict[str, Any]:
    """Bench-wide aggregation: how many proteins by delivery status."""
    rows = list(rows)
    n = len(rows)
    by_status: dict[str, int] = {}
    n_components_total = 0
    n_high = n_medium = n_low = n_failed = 0
    for r in rows:
        by_status[r.overall_status] = by_status.get(r.overall_status, 0) + 1
        for c in r.components:
            n_components_total += 1
            if c.confidence is Confidence.HIGH: n_high += 1
            elif c.confidence is Confidence.MEDIUM: n_medium += 1
            elif c.confidence is Confidence.LOW: n_low += 1
            else: n_failed += 1
    return {
        "n_proteins": n,
        "by_overall_status": by_status,
        "by_overall_status_pct": {
            k: (100 * v / n if n else 0.0) for k, v in by_status.items()
        },
        "component_totals": {
            "high": n_high, "medium": n_medium, "low": n_low, "failed": n_failed,
            "total": n_components_total,
        },
    }
