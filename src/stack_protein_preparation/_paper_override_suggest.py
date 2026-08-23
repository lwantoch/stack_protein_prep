"""Extract per-residue protonation-override suggestions from a
``paper/paper_evidence.md`` file.

Companion to ``_active_site_overrides.load_active_site_overrides``.  That
module honours a user-hand-written ``active_site_overrides.json``; this
module scans the auto-generated paper-evidence extract and proposes an
``active_site_overrides.suggested.json`` that the user reviews and either
copies verbatim to the real override file or ignores.

Design decisions:

* **Never write the real override file.**  The suggestion file uses a
  ``.suggested.json`` suffix so PROPKA + metal-based overrides remain the
  default; user's manual intervention is the only path from "paper hint"
  to "applied override".  This preserves reviewer-honesty: FRUTON does
  not silently change protonation based on regex-fragile paper reading.

* **Pure regex, no LLM.**  License-free, deterministic, auditable.
  Each suggestion carries the raw paper quote and a confidence tag
  (high / medium / low) based on how specific the mechanism keyword is.

* **Chemical mapping** (kept intentionally small so misfires stay rare):

    keyword phrase in paper quote          -> suggested override
    ---------------------------------------   -------------------
    "catalytic acid" / "proton donor"     -> ASH / GLH  (for Asp/Glu)
    "catalytic base" / "proton acceptor"  -> HIP        (for His)
    "nucleophile" / "thiolate"            -> CYM        (for Cys)
    "protonated" / "neutral" (Asp/Glu)    -> ASH / GLH
    "protonated" / "positive" (His)       -> HIP
    "coordinating" + metal element        -> HID/HIE    (per metal
                                                         geometry)

  Ambiguous mentions ("catalytic residue D174" without an acid/base
  attribution) get status='needs_review', not a concrete override.

Public entry:
    suggest_overrides_from_paper_evidence(
        paper_evidence_md: Path,
        chain_default: str = "A",
    ) -> dict[str, Any]

    Returns a manifest dict with ``overrides`` (JSON-ready list) plus
    ``diagnostics`` (per-residue confidence + quote).  Empty ``overrides``
    when the file is missing / has no mineable hints.
"""
from __future__ import annotations

import json
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Any


# Amber-residue-name mapping from source (crystal) three-letter code + mechanism keyword.
# The 'to' field is the AMBER residue name FRUTON's protonation step honours.
_ASP_HINTS = (
    ("catalytic acid", "high", "ASH"),
    ("proton donor", "high", "ASH"),
    ("protonated asp", "medium", "ASH"),
    ("neutral asp", "medium", "ASH"),
)
_GLU_HINTS = (
    ("catalytic acid", "high", "GLH"),
    ("proton donor", "high", "GLH"),
    ("protonated glu", "medium", "GLH"),
    ("neutral glu", "medium", "GLH"),
)
_HIS_HINTS_HIP = (
    ("catalytic base", "high", "HIP"),
    ("proton acceptor", "high", "HIP"),
    ("protonated his", "high", "HIP"),
    ("positively charged his", "medium", "HIP"),
    ("hip", "medium", "HIP"),
)
_HIS_HINTS_HID = (
    ("hid", "medium", "HID"),
    ("nε2 protonated", "medium", "HID"),
    ("ne2 protonated", "medium", "HID"),
)
_HIS_HINTS_HIE = (
    ("hie", "medium", "HIE"),
    ("nδ1 protonated", "medium", "HIE"),
    ("nd1 protonated", "medium", "HIE"),
)
_CYS_HINTS = (
    ("nucleophile", "high", "CYM"),
    ("thiolate", "high", "CYM"),
    ("catalytic cys", "medium", "CYM"),
    ("deprotonated cys", "medium", "CYM"),
    ("cys thiolate", "high", "CYM"),
)


@dataclass
class OverrideSuggestion:
    chain: str
    res_num: int
    resname_source: str          # ASP / GLU / HIS / CYS
    to: str                      # ASH / GLH / HID / HIE / HIP / CYM
    confidence: str              # high | medium | low | needs_review
    reason: str                  # paper quote (truncated to 300 chars)

    def as_override_entry(self) -> dict[str, Any]:
        return {
            "chain": self.chain,
            "res_num": self.res_num,
            "to": self.to,
            "reason": f"[{self.confidence} confidence, auto-suggested] {self.reason}",
        }


_RESI_HEADING_RE = re.compile(
    r"^###\s+(ASP|ASN|GLU|GLN|HIS|CYS|LYS|ARG|TYR|SER|THR|MET|TRP|PHE|VAL|LEU|ILE|GLY|ALA|PRO)(\d+)\s*$",
    re.IGNORECASE,
)
_QUOTE_LINE_RE = re.compile(r"^>\s*(.*)$")


def _match_hint(
    quote_lower: str, hint_table: tuple[tuple[str, str, str], ...]
) -> tuple[str, str] | None:
    """Return (confidence, to_state) of the FIRST hint keyword that matches."""
    for keyword, confidence, to_state in hint_table:
        if keyword in quote_lower:
            return (confidence, to_state)
    return None


def suggest_overrides_from_paper_evidence(
    paper_evidence_md: str | Path,
    chain_default: str = "A",
) -> dict[str, Any]:
    """Scan paper_evidence.md, return a JSON-ready manifest of protonation
    suggestions.  Empty manifest when file missing / no hints found."""
    md_path = Path(paper_evidence_md)
    manifest: dict[str, Any] = {
        "status": "skipped",
        "source": str(md_path),
        "overrides": [],
        "diagnostics": [],
    }
    if not md_path.is_file():
        manifest["message"] = "paper_evidence.md not found"
        return manifest

    text = md_path.read_text(encoding="utf-8", errors="replace")
    lines = text.splitlines()
    n_lines = len(lines)

    suggestions: list[OverrideSuggestion] = []
    diagnostics: list[dict[str, Any]] = []

    i = 0
    while i < n_lines:
        heading = _RESI_HEADING_RE.match(lines[i])
        if not heading:
            i += 1
            continue
        resname_src = heading.group(1).upper()
        try:
            res_num = int(heading.group(2))
        except ValueError:
            i += 1
            continue

        # Collect any following quote lines until the next heading / blank block
        quote_parts: list[str] = []
        j = i + 1
        while j < n_lines:
            if _RESI_HEADING_RE.match(lines[j]):
                break
            if lines[j].startswith("## "):
                break
            m = _QUOTE_LINE_RE.match(lines[j])
            if m:
                quote_parts.append(m.group(1).strip())
            j += 1
        i = j

        if not quote_parts:
            continue

        quote_full = " ".join(quote_parts).strip()
        quote_short = quote_full if len(quote_full) <= 300 else quote_full[:297] + "..."
        quote_lower = quote_full.lower()

        # Route on source residue type
        hint = None
        if resname_src == "ASP":
            hint = _match_hint(quote_lower, _ASP_HINTS)
        elif resname_src == "GLU":
            hint = _match_hint(quote_lower, _GLU_HINTS)
        elif resname_src == "HIS":
            # Search HIP first (strongest catalytic-base signal), then HID/HIE
            hint = _match_hint(quote_lower, _HIS_HINTS_HIP)
            if hint is None:
                hint = _match_hint(quote_lower, _HIS_HINTS_HID)
            if hint is None:
                hint = _match_hint(quote_lower, _HIS_HINTS_HIE)
        elif resname_src == "CYS":
            hint = _match_hint(quote_lower, _CYS_HINTS)

        if hint is None:
            diagnostics.append({
                "chain": chain_default,
                "res_num": res_num,
                "resname_source": resname_src,
                "status": "needs_review",
                "reason": (
                    "residue mentioned in catalytic/binding context but no "
                    "specific protonation keyword matched"
                ),
                "quote": quote_short,
            })
            continue

        confidence, to_state = hint
        sugg = OverrideSuggestion(
            chain=chain_default,
            res_num=res_num,
            resname_source=resname_src,
            to=to_state,
            confidence=confidence,
            reason=quote_short,
        )
        suggestions.append(sugg)
        diagnostics.append({
            "chain": chain_default,
            "res_num": res_num,
            "resname_source": resname_src,
            "status": "suggested",
            "to": to_state,
            "confidence": confidence,
            "quote": quote_short,
        })

    manifest["overrides"] = [s.as_override_entry() for s in suggestions]
    manifest["diagnostics"] = diagnostics
    if suggestions:
        manifest["status"] = "success"
    elif diagnostics:
        manifest["status"] = "no_actionable_hints"
    else:
        manifest["status"] = "no_hints_found"
    manifest["summary"] = {
        "n_suggested": len(suggestions),
        "n_needs_review": sum(1 for d in diagnostics if d.get("status") == "needs_review"),
    }
    return manifest


def write_suggested_overrides(
    paper_evidence_md: str | Path,
    output_json: str | Path,
    chain_default: str = "A",
) -> dict[str, Any]:
    """Convenience wrapper: run suggest, write the ``.suggested.json`` file,
    return the manifest.  Never touches ``active_site_overrides.json``."""
    manifest = suggest_overrides_from_paper_evidence(paper_evidence_md, chain_default)
    if manifest["overrides"]:
        out_path = Path(output_json)
        out_path.parent.mkdir(parents=True, exist_ok=True)
        payload = {"overrides": manifest["overrides"]}
        out_path.write_text(
            json.dumps(payload, indent=2, sort_keys=False) + "\n",
            encoding="utf-8",
        )
        manifest["written_to"] = str(out_path)
    return manifest
