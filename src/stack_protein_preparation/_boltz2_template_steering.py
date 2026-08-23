"""Boltz-2 template-steering fallback scaffold (Iteration 5).

Reviewer perspective: FRUTON's default AF-splice + MODELLER LoopModel
plateaus on borderline cases where the AF fill has low pLDDT AND the
crystal template doesn't provide enough context to close the loop
cleanly.  Boltz-2 (Ingraham et al. 2024) is a diffusion-based
structure-prediction model that accepts a *template PDB* as soft prior
— which is exactly the missing "steer AF toward the crystal template"
step FRUTON needs for these cases.

This module is a **scaffold** in the same fail-open style as
_rfdiffusion2_gap.py.  Public API a future integration will use:

    spec = BoltzTemplateSteeringSpec(
        chain="A",
        first_resnum=50, last_resnum=60,
        sequence="AGKGQLTVNK",          # 10 residues from UniProt
        template_pdb=crystal_pdb,       # the resolved-crystal context
    )
    result = attempt_template_steering(context_pdb, spec, output_dir)
    if result.accepted:
        # promote result.candidate_pdbs to the pipeline
        ...
    else:
        # fall back to MODELLER LoopModel + rollback
        ...

Design goals (identical to _rfdiffusion2_gap):
- Never require Boltz-2 to be installed for tests / normal runs.
- Never call the network (do not download weights).
- Fail-open cascade: bad spec → no binary → subprocess timeout →
  subprocess non-zero → no candidate PDBs → success.
- Standardised manifest.json so a downstream refactor can swap
  Boltz-2 for Boltz-3 / AlphaFold-Multimer / ESMFold without
  changing the call site.

Pure Python + stdlib.  License-free scaffolding.
"""
from __future__ import annotations

import json
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable

from stack_protein_preparation._external_generator_runner import (
    execute_generator_subprocess,
    resolve_binary_from_candidates,
)


DEFAULT_BINARY_CANDIDATES: tuple[str, ...] = (
    "boltz",
    "boltz2",
    "boltz-2",
)
DEFAULT_TIMEOUT_S = 900  # 15 min per steering pass; heavy diffusion


@dataclass(frozen=True)
class BoltzTemplateSteeringSpec:
    """One gap to steer.

    ``template_pdb`` is a resolved-crystal PDB used as the soft prior;
    Boltz-2's template channel will bias the fill toward this
    geometry when the AF fill deviates.  ``pLDDT_threshold`` gates
    whether the caller wanted to fall back to steering (informational
    only — the scaffold does not enforce a threshold here).
    """
    chain: str
    first_resnum: int
    last_resnum: int
    sequence: str
    template_pdb: Path
    pLDDT_threshold: float = 50.0

    def n_residues(self) -> int:
        return self.last_resnum - self.first_resnum + 1

    def sanity_check(self) -> list[str]:
        problems: list[str] = []
        if self.n_residues() < 1:
            problems.append("empty gap (first_resnum > last_resnum)")
        if len(self.sequence) != self.n_residues():
            problems.append(
                f"sequence length {len(self.sequence)} does not match "
                f"gap length {self.n_residues()}"
            )
        if not Path(self.template_pdb).is_file():
            problems.append(f"template pdb not found: {self.template_pdb}")
        return problems


@dataclass
class BoltzSteeringAttempt:
    ran: bool
    accepted: bool
    fallback_reason: str | None = None
    candidate_pdbs: list[Path] = field(default_factory=list)
    diagnostics: list[str] = field(default_factory=list)
    used_binary: str | None = None

    def to_dict(self) -> dict:
        return {
            "ran": self.ran,
            "accepted": self.accepted,
            "fallback_reason": self.fallback_reason,
            "used_binary": self.used_binary,
            "n_candidates": len(self.candidate_pdbs),
            "candidate_pdbs": [str(p) for p in self.candidate_pdbs],
            "diagnostics": list(self.diagnostics),
        }


def _resolve_binary(candidates: Iterable[str] | None = None) -> str | None:
    """First-on-PATH from ``candidates`` (default DEFAULT_BINARY_CANDIDATES).

    Thin wrapper preserved for test-suite import stability; new code
    should call resolve_binary_from_candidates from
    _external_generator_runner directly.
    """
    return resolve_binary_from_candidates(
        candidates or DEFAULT_BINARY_CANDIDATES,
    )


def _write_steering_manifest(
    output_dir: Path,
    spec: BoltzTemplateSteeringSpec,
    context_pdb: Path,
) -> Path:
    """Write a scaffold manifest.json the CLI would consume.

    Generator-agnostic layout so a downstream swap to Boltz-3 or an
    AlphaFold-Multimer path can extend it without breaking callers.
    """
    manifest = {
        "generator": "boltz2_template_steering",
        "context_pdb": str(context_pdb.resolve()),
        "template_pdb": str(Path(spec.template_pdb).resolve()),
        "gap": {
            "chain": spec.chain,
            "first_resnum": spec.first_resnum,
            "last_resnum": spec.last_resnum,
            "sequence": spec.sequence,
        },
        "steering": {
            "mode": "template_bias",
            "pLDDT_threshold_below_which_steering_was_triggered":
                spec.pLDDT_threshold,
        },
    }
    output_dir.mkdir(parents=True, exist_ok=True)
    path = output_dir / "boltz_manifest.json"
    path.write_text(json.dumps(manifest, indent=2))
    return path


def attempt_template_steering(
    context_pdb: str | Path,
    spec: BoltzTemplateSteeringSpec,
    output_dir: str | Path,
    n_candidates: int = 4,
    timeout_seconds: int = DEFAULT_TIMEOUT_S,
    binary_candidates: Iterable[str] | None = None,
) -> BoltzSteeringAttempt:
    """Try to generate template-steered loop candidates via Boltz-2.

    Fail-open cascade (mirrors _rfdiffusion2_gap.attempt_gap_fill):
    1. Missing context PDB → ran=False.
    2. Invalid spec (sanity_check problems) → ran=False.
    3. No binary on PATH → ran=False + fallback_reason (Boltz-2 not
       installed; pipeline degrades to MODELLER LoopModel + rollback).
    4. subprocess non-zero / timeout / OSError → ran=True + accepted=False.
    5. Subprocess succeeds but writes no candidate_*.pdb → ran=True +
       accepted=False.
    6. Success → BoltzSteeringAttempt(candidate_pdbs=[...], used_binary=...).
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    context_pdb = Path(context_pdb)

    if not context_pdb.is_file():
        return BoltzSteeringAttempt(
            ran=False, accepted=False,
            fallback_reason=f"context pdb not found: {context_pdb}",
        )

    problems = spec.sanity_check()
    if problems:
        return BoltzSteeringAttempt(
            ran=False, accepted=False,
            fallback_reason="spec invalid: " + "; ".join(problems),
            diagnostics=list(problems),
        )

    binary = _resolve_binary(binary_candidates)
    if binary is None:
        return BoltzSteeringAttempt(
            ran=False, accepted=False,
            fallback_reason=(
                "no boltz binary on PATH; tried "
                + ",".join(binary_candidates or DEFAULT_BINARY_CANDIDATES)
                + ". Falls back to MODELLER LoopModel + rollback."
            ),
        )

    manifest = _write_steering_manifest(output_dir, spec, context_pdb)

    run = execute_generator_subprocess(
        binary=binary,
        manifest_path=manifest,
        output_dir=output_dir,
        n_candidates=n_candidates,
        timeout_seconds=timeout_seconds,
        tool_name_for_messages="boltz",
    )
    return BoltzSteeringAttempt(
        ran=run.ran, accepted=run.accepted,
        fallback_reason=run.fallback_reason,
        candidate_pdbs=list(run.candidate_pdbs),
        diagnostics=list(run.diagnostics),
        used_binary=run.used_binary,
    )


def summarise(attempt: BoltzSteeringAttempt) -> str:
    if not attempt.ran:
        return f"boltz2: skipped ({attempt.fallback_reason})"
    if not attempt.accepted:
        return f"boltz2: ran but rejected ({attempt.fallback_reason})"
    return (
        f"boltz2: {len(attempt.candidate_pdbs)} candidates via "
        f"{attempt.used_binary}"
    )
