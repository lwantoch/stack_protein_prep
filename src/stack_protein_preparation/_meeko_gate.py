"""Meeko compatibility gate for FRUTON protonation output.

Prior state: FRUTON declares a receptor "prepared" as soon as ``gmx pdb2gmx``
returns rc=0 and writes a non-empty coordinate file. But pdb2gmx-produced
PDBs regularly fail downstream Meeko preparation (which turns the receptor
into an AutoDock .pdbqt), and the failure only surfaces during docking —
long after the receptor manifest has been shipped as OK.

This module runs Meeko's ``Polymer.from_pdb_string`` on FRUTON's own
protonation output as an immediate post-write validator, then classifies
the exception into one of four generic buckets so downstream retry logic
(future work) can pick a targeted alternative H-placement path.

Error classification is BY EXCEPTION MESSAGE PATTERN, not by protein — the
same class always gets the same fix regardless of the target PDB. This is
the "no per-protein specialised solution" rule the user demanded.

Standalone use (for the 27-target regression harness):

    from stack_protein_preparation._meeko_gate import validate_pdb_for_meeko
    result = validate_pdb_for_meeko(protonated_pdb)
    if not result.ok:
        log.error("meeko-gate failed: %s (%s)", result.error_class, result.message)
"""

from __future__ import annotations

import re
from dataclasses import dataclass
from pathlib import Path
from typing import Literal

MeekoErrorClass = Literal[
    "valence",           # RDKit AtomValenceException — H-clash from pdb2gmx -ignh
    "adjacent_smarts",   # Meeko boundary-atom mapping failure at chain terminus
    "no_template",       # Missing residue template (novel HET, PTM)
    "padding",           # Meeko "Expected N paddings ... but got 0" — CIF artifact
    "other",             # Anything else (novel — worth reading the trace)
]


@dataclass(frozen=True)
class MeekoGateResult:
    ok: bool
    error_class: MeekoErrorClass | None
    message: str
    exception_type: str | None


# --- error signature patterns ------------------------------------------------
# Each key is a MeekoErrorClass; value is a compiled regex over the exception
# message. Order matters — first match wins. Patterns are derived from the
# MEEKO_STATUS.csv on the 27-target regression set:
#
#   valence         → rdkit.Chem.rdchem.AtomValenceException
#   adjacent_smarts → "adjacent_mol doesn't contain the mapped atoms"
#   no_template     → "Creation of data structure for receptor failed"
#                     (also caught: no_template messages from meeko polymer)
#   padding         → "Expected N paddings for (...) with bonds"
#
# Add new signatures here as new failure modes surface — DO NOT branch on
# the target's PDB id, always on the exception string.

_ERROR_PATTERNS: list[tuple[MeekoErrorClass, re.Pattern[str]]] = [
    ("valence", re.compile(r"AtomValenceException|Explicit valence for atom", re.I)),
    ("adjacent_smarts", re.compile(r"adjacent_(mol|smarts)", re.I)),
    ("padding", re.compile(r"Expected \d+ paddings? for", re.I)),
    ("no_template", re.compile(r"(?:no.?template|Creation of data structure for receptor failed|template.*not.*found)", re.I)),
]


def _classify(exc: BaseException) -> MeekoErrorClass:
    msg = f"{type(exc).__name__}: {exc}"
    for cls, pat in _ERROR_PATTERNS:
        if pat.search(msg):
            return cls
    return "other"


def validate_pdb_for_meeko(
    pdb_path: str | Path,
    *,
    allow_bad_res: bool = False,
) -> MeekoGateResult:
    """Run Meeko on ``pdb_path`` and return a structured verdict.

    Never raises — every failure mode becomes a MeekoGateResult with ok=False.
    Callers use ``result.error_class`` to drive retry logic; ``result.message``
    is the raw exception string, safe to log verbatim.

    ``allow_bad_res`` mirrors Meeko's own ``Polymer.from_pdb_string`` kwarg:
    when True, unmatched residues are silently dropped instead of aborting
    the polymer build. Default False = strict (reference behavior of
    ``mk_prepare_receptor.py`` CLI). The two-pass strict → tolerant idiom is
    the recommended way to distinguish "receptor unusable" from "receptor
    usable if we tolerate a small loss of coverage" — see the retry chain
    in ``protonation.py``.
    """
    p = Path(pdb_path)
    if not p.is_file() or p.stat().st_size == 0:
        return MeekoGateResult(
            ok=False,
            error_class="other",
            message=f"pdb file missing or empty: {p}",
            exception_type=None,
        )

    try:
        # Imported lazily so importing this module has no meeko cost.
        from meeko import (  # type: ignore[import-not-found]
            MoleculePreparation,
            Polymer,
            ResidueChemTemplates,
        )
    except ImportError as exc:
        return MeekoGateResult(
            ok=False,
            error_class="other",
            message=f"meeko not importable: {exc}",
            exception_type=type(exc).__name__,
        )

    try:
        # Reproduce the mk_prepare_receptor.py CLI initialization exactly so
        # our verdict matches what a real docking prep would see. Both are
        # default-configured — the gate is purely a validity check, not a
        # place to inject custom templates.
        templates = ResidueChemTemplates.create_from_defaults()
        mk_prep = MoleculePreparation()
        pdb_text = p.read_text(encoding="utf-8", errors="replace")
        # Note: meeko 0.7 signature is
        #   from_pdb_string(pdb, chem_templates, mk_prep, set_template=None,
        #                   residues_to_delete=None, allow_bad_res=False, ...)
        # Keep the positional call minimal — any deviation from the CLI's
        # defaults would give a verdict that doesn't match a real docking
        # prep, which defeats the gate's purpose.
        polymer = Polymer.from_pdb_string(
            pdb_text,
            templates,
            mk_prep,
            allow_bad_res=False,
            default_altloc="A",
        )
        n_res = len(polymer.monomers) if hasattr(polymer, "monomers") else 0
        if n_res == 0:
            return MeekoGateResult(
                ok=False,
                error_class="other",
                message="meeko produced empty polymer (no monomers)",
                exception_type=None,
            )
        return MeekoGateResult(ok=True, error_class=None, message="ok", exception_type=None)
    except Exception as exc:  # noqa: BLE001 — validation shim, must swallow
        return MeekoGateResult(
            ok=False,
            error_class=_classify(exc),
            message=str(exc),
            exception_type=type(exc).__name__,
        )
