#!/usr/bin/env python3
"""Demo: naive pH=7 protonation rule vs FRUTON's paper-evidence override.

Reviewer-facing demonstration for the pdbfixer-removal justification
(see src/stack_protein_preparation/data/pdbfixer_removal_justification.md).

We emulate a naive default-pH histidine-tautomer rule (which is what
pdbfixer.addMissingHydrogens uses when pH=7) and contrast it with
FRUTON's paper-evidence override flow.  Pure Python, no pdbfixer
dependency required.  Prints a side-by-side table for three
publication-scale scenarios.

Run:
    python scripts/demo_pdbfixer_vs_fruton_protonation.py \\
      [--evidence <path/to/paper_evidence.md>]

The default evidence corpus is baked into this script so it runs
standalone; supply --evidence to try a real per-protein file.
"""
from __future__ import annotations

import argparse
import re
import sys
import tempfile
from dataclasses import dataclass
from pathlib import Path

from stack_protein_preparation._paper_override_suggest import (
    suggest_overrides_from_paper_evidence,
)


NAIVE_PH = 7.0

# Default set of active-site HIS scenarios paraphrased from published
# literature (see JCTC / JBC references in the .md).  Coordinates are
# not needed for a protonation-rule demo.
DEFAULT_SCENARIOS: list[tuple[str, str, str, str]] = [
    (
        "carbonic_anhydrase_II_H94_Zn",
        "HIS94",
        "H94 coordinates the catalytic Zn2+ ion; the imidazole Nδ1 "
        "protonated tautomer is required for the Zn2+ binding geometry.",
        "HIE",
    ),
    (
        "trypsin_H57_catalytic_triad",
        "HIS57",
        "H57 is the catalytic base in the Ser-His-Asp triad and "
        "transiently accepts a proton from Ser195 to form the "
        "acyl-enzyme intermediate.",
        "HIP",
    ),
    (
        "papain_C25_thiolate_nucleophile",
        "CYS25",
        "C25 is the catalytic nucleophile in papain; the thiolate "
        "form (deprotonated Cys) attacks the substrate carbonyl.",
        "CYM",
    ),
    (
        "carbonic_anhydrase_II_D162_needs_review",
        "ASP162",
        "D162 sits at the entrance to the active-site channel and "
        "participates in the hydrogen-bonding network but no specific "
        "acid/base attribution has been published.",
        "ASP",
    ),
]


@dataclass
class Verdict:
    scenario: str
    residue: str
    naive_default_ph7: str
    fruton_override: str | None
    literature_recommended: str

    @property
    def naive_wrong(self) -> bool:
        return self.naive_default_ph7 != self.literature_recommended

    @property
    def fruton_correct(self) -> bool:
        if self.fruton_override is None:
            return False
        return self.fruton_override == self.literature_recommended


def naive_ph_rule(resname_source: str, pH: float = NAIVE_PH) -> str:
    """Emulate the pdbfixer / naive-default protonation rule at pH.

    HIS: below pH 6.0 → HIP; otherwise → HIE (pdbfixer default at pH 7).
    ASP: always ASP (deprotonated -COO−) above pKa ~3.7.
    GLU: always GLU above pKa ~4.3.
    CYS: always CYS (protonated -SH) below pKa ~8.3; the thiolate CYM
         form requires explicit override.
    """
    source = resname_source.upper()
    if source == "HIS":
        return "HIP" if pH < 6.0 else "HIE"
    if source == "ASP":
        return "ASP"
    if source == "GLU":
        return "GLU"
    if source == "CYS":
        return "CYS"
    return source


def _cli_table(rows: list[Verdict]) -> str:
    """Format the verdicts as a plain-text table."""
    out: list[str] = []
    out.append(
        f"{'Scenario':<38} {'Residue':<8} "
        f"{'naive pH=7':>10}  {'FRUTON':>7}  {'lit.':>6}  {'winner':>10}"
    )
    out.append("-" * 90)
    for v in rows:
        if v.fruton_correct and v.naive_wrong:
            winner = "FRUTON"
        elif v.fruton_correct and not v.naive_wrong:
            winner = "tie"
        elif not v.fruton_correct and not v.naive_wrong:
            winner = "tie(naive)"
        else:
            winner = "review!"
        out.append(
            f"{v.scenario:<38} {v.residue:<8} "
            f"{v.naive_default_ph7:>10}  "
            f"{(v.fruton_override or '(none)'):>7}  "
            f"{v.literature_recommended:>6}  "
            f"{winner:>10}"
        )
    return "\n".join(out) + "\n"


_RES_HEADING_NUM_RE = re.compile(r"^([A-Z]+)(\d+)$")


def _scenarios_as_paper_evidence_md(
    scenarios: list[tuple[str, str, str, str]],
) -> str:
    """Emit a paper_evidence.md that suggest_overrides_from_paper_evidence
    can parse: '### HIS94' heading + '> quote' body.
    """
    out_lines: list[str] = ["# Paper evidence — demo\n\n"]
    for (_scenario, residue, quote, _rec) in scenarios:
        m = _RES_HEADING_NUM_RE.match(residue)
        if not m:
            continue
        prefix, num = m.group(1), m.group(2)
        out_lines.append(f"### {prefix}{num}\n")
        out_lines.append(f"> {quote}\n\n")
    return "".join(out_lines)


def _resname_from_label(label: str) -> str:
    m = _RES_HEADING_NUM_RE.match(label)
    return m.group(1) if m else ""


def _load_evidence_path(path: Path | None, scenarios) -> Path:
    """Return a real path (real file or a tmp one built from scenarios)."""
    if path is not None:
        return path
    tmp = Path(tempfile.mkdtemp(prefix="demo_evidence_")) / "paper_evidence.md"
    tmp.write_text(_scenarios_as_paper_evidence_md(scenarios), encoding="utf-8")
    return tmp


def _residue_number(res_label: str) -> int | None:
    m = _RES_HEADING_NUM_RE.match(res_label)
    if not m:
        return None
    try:
        return int(m.group(2))
    except ValueError:
        return None


def build_verdicts(
    scenarios: list[tuple[str, str, str, str]],
    evidence_path: Path,
) -> list[Verdict]:
    manifest = suggest_overrides_from_paper_evidence(evidence_path)
    by_resnum: dict[int, str] = {}
    for entry in manifest.get("overrides", []):
        try:
            by_resnum[int(entry["res_num"])] = entry["to"]
        except (KeyError, TypeError, ValueError):
            continue

    out: list[Verdict] = []
    for (scenario, residue, _text, literature_recommended) in scenarios:
        source_resname = _resname_from_label(residue)
        naive = naive_ph_rule(source_resname, NAIVE_PH)
        resnum = _residue_number(residue)
        fruton = by_resnum.get(resnum) if resnum is not None else None
        out.append(Verdict(
            scenario=scenario,
            residue=residue,
            naive_default_ph7=naive,
            fruton_override=fruton,
            literature_recommended=literature_recommended,
        ))
    return out


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--evidence", type=Path, default=None,
                   help="Optional path to a per-protein paper_evidence.md")
    ns = p.parse_args(argv)

    evidence_path = _load_evidence_path(ns.evidence, DEFAULT_SCENARIOS)
    verdicts = build_verdicts(DEFAULT_SCENARIOS, evidence_path)
    print(_cli_table(verdicts))

    n_naive_wrong = sum(1 for v in verdicts if v.naive_wrong)
    n_fruton_correct = sum(1 for v in verdicts if v.fruton_correct)
    print(
        f"summary: naive-pH-rule wrong on {n_naive_wrong}/{len(verdicts)}; "
        f"FRUTON paper-override correct on {n_fruton_correct}/{len(verdicts)}"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
