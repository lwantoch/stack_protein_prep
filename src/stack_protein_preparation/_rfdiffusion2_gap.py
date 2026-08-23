"""Constrained-Diffusion gap-filler scaffold (RFdiffusion2 preview).

Reviewer perspective: FRUTON's current gap-fill sequence is
splice-from-AF → MODELLER LoopModel → adaptive rollback.  For long
gaps flanking flexible regions the LoopModel step can plateau at
un-fittable conformers.  RFdiffusion2 (Ingraham et al. 2023 successor,
see arXiv:2510.14989 for the Constrained-Diffusion follow-up) is a
denoising diffusion generator with contact / secondary-structure /
distance constraints, well-suited to loop filling when the flanking
Cα atoms are fixed.

This module is a **scaffold**: it defines the public API a future
integration will use, and a fail-open subprocess wrapper for the
external ``rfdiffusion2`` CLI.  When the binary is not on PATH (the
default for most CESGA and laptop environments) the fallback path is
taken and the pipeline degrades to MODELLER LoopModel + rollback.

Design goals:
- Never require RFdiffusion2 to be installed for tests / normal runs.
- Never call the network.
- Emit a machine-readable manifest so the higher-level filler
  (``_filler_alphafold.py``) can decide whether to consume the
  candidates or fall back to MODELLER.
- Standardise the constraint spec so a downstream refactor
  (Davide) can swap RFdiffusion2 for RFdiffusion3 / Boltz-2 /
  Constrained-Diffusion without changing the call site.

Pure Python + stdlib.  License-free.  BSD-3-ish scaffolding.
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
    "rfdiffusion2",
    "rfdiffusion",   # legacy name
    "constrained_diffusion",  # follow-up per arXiv:2510.14989
)
DEFAULT_TIMEOUT_S = 600  # 10 min per gap; a real integration may want longer


@dataclass(frozen=True)
class GapConstraintSpec:
    """One gap to fill.

    ``anchor_ca_by_resnum`` maps residue numbers on either side of the
    gap to their Cα XYZ (from the crystal), which the diffusion model
    should honour as hard constraints.  ``sequence`` is the primary
    sequence FRUTON wants for the gap (from UniProt).
    """
    chain: str
    first_resnum: int
    last_resnum: int
    sequence: str
    anchor_ca_by_resnum: dict[int, tuple[float, float, float]] = field(default_factory=dict)

    def n_residues(self) -> int:
        return self.last_resnum - self.first_resnum + 1

    def sanity_check(self) -> list[str]:
        """Return list of human-readable warnings / hard-fail reasons."""
        problems: list[str] = []
        if self.n_residues() < 1:
            problems.append("empty gap (first_resnum > last_resnum)")
        if len(self.sequence) != self.n_residues():
            problems.append(
                f"sequence length {len(self.sequence)} does not match "
                f"gap length {self.n_residues()}"
            )
        if not self.anchor_ca_by_resnum:
            problems.append("no anchor Cα constraints supplied")
        return problems


@dataclass
class RFDiffusionAttempt:
    ran: bool                         # True if a subprocess actually launched
    accepted: bool                    # True if we got at least one candidate PDB back
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


def _write_constraint_manifest(
    output_dir: Path,
    spec: GapConstraintSpec,
    context_pdb: Path,
) -> Path:
    """Write a scaffold manifest.json the CLI would consume.

    Kept deliberately simple + generator-agnostic so a future
    integration (Boltz-2, RFdiffusion3, …) can extend it.
    """
    manifest = {
        "generator": "constrained_diffusion",
        "context_pdb": str(context_pdb.resolve()),
        "gap": {
            "chain": spec.chain,
            "first_resnum": spec.first_resnum,
            "last_resnum": spec.last_resnum,
            "sequence": spec.sequence,
        },
        "anchor_constraints": {
            "ca_xyz_by_resnum": {
                str(k): list(v) for k, v in spec.anchor_ca_by_resnum.items()
            },
        },
    }
    output_dir.mkdir(parents=True, exist_ok=True)
    path = output_dir / "gap_manifest.json"
    path.write_text(json.dumps(manifest, indent=2))
    return path


def attempt_gap_fill(
    context_pdb: str | Path,
    gap: GapConstraintSpec,
    output_dir: str | Path,
    n_candidates: int = 4,
    timeout_seconds: int = DEFAULT_TIMEOUT_S,
    binary_candidates: Iterable[str] | None = None,
) -> RFDiffusionAttempt:
    """Try to generate loop candidates via a Constrained-Diffusion CLI.

    Fail-open cascade:
    1. Reject on sanity-check problems (bad gap spec) → ran=False.
    2. If no binary on PATH → ran=False + fallback_reason.
    3. Otherwise write manifest.json, launch subprocess.
    4. If subprocess non-zero / times out → ran=True + accepted=False.
    5. If subprocess succeeds → collect ``candidate_*.pdb`` files in
       output_dir; accepted=True when count > 0.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    context_pdb = Path(context_pdb)

    if not context_pdb.is_file():
        return RFDiffusionAttempt(
            ran=False, accepted=False,
            fallback_reason=f"context pdb not found: {context_pdb}",
        )

    problems = gap.sanity_check()
    if problems:
        return RFDiffusionAttempt(
            ran=False, accepted=False,
            fallback_reason="gap spec invalid: " + "; ".join(problems),
            diagnostics=list(problems),
        )

    binary = _resolve_binary(binary_candidates)
    if binary is None:
        return RFDiffusionAttempt(
            ran=False, accepted=False,
            fallback_reason=(
                "no constrained-diffusion binary on PATH; tried "
                + ",".join(binary_candidates or DEFAULT_BINARY_CANDIDATES)
                + ". Falls back to MODELLER LoopModel + rollback."
            ),
        )

    manifest = _write_constraint_manifest(output_dir, gap, context_pdb)

    run = execute_generator_subprocess(
        binary=binary,
        manifest_path=manifest,
        output_dir=output_dir,
        n_candidates=n_candidates,
        timeout_seconds=timeout_seconds,
        tool_name_for_messages="diffusion",
    )
    return RFDiffusionAttempt(
        ran=run.ran, accepted=run.accepted,
        fallback_reason=run.fallback_reason,
        candidate_pdbs=list(run.candidate_pdbs),
        diagnostics=list(run.diagnostics),
        used_binary=run.used_binary,
    )


def summarise(attempt: RFDiffusionAttempt) -> str:
    if not attempt.ran:
        return f"rfdiffusion2: skipped ({attempt.fallback_reason})"
    if not attempt.accepted:
        return f"rfdiffusion2: ran but rejected ({attempt.fallback_reason})"
    return (
        f"rfdiffusion2: {len(attempt.candidate_pdbs)} candidates via "
        f"{attempt.used_binary}"
    )
