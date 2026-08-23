"""ESMFold single-sequence structure-prediction scaffold.

Reviewer perspective: FRUTON's default AF-splice path assumes the
UniProt→AF cache has a good model.  For borderline cases where AF's
MSA is ambiguous (multiple isoforms, novel splice variant, edge-case
species), a single-sequence predictor like ESMFold (Rives / Lin 2023)
gives a template-free fallback the caller can drop into the splice
step without waiting for a fresh AF run.

This module is the third scaffold layered on top of
_external_generator_runner (after _rfdiffusion2_gap and
_boltz2_template_steering).  It exists specifically to demonstrate
that the shared runner absorbs new generator integrations in ~50
lines instead of the ~150 the pre-refactor duplicated pattern needed.

Design (identical fail-open cascade via the shared runner):
- Never require esmfold to be installed for tests / normal runs.
- Never call the network (do not download weights).
- Fail-open when binary missing → pipeline degrades to AF-splice.
- Standardised manifest.json so ESMFold-2 / OmegaFold / MSA-Free
  fallbacks can extend without changing the call site.

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
    "esmfold",
    "esm-fold",
    "esmfold_inference",
)
DEFAULT_TIMEOUT_S = 600  # 10 min per short chain; long chains need longer


@dataclass(frozen=True)
class EsmFoldSequenceSpec:
    """One sequence to predict.

    ``chain_id`` is used as the output PDB's chain label so the
    prediction can be merged into a multi-chain assembly.
    """
    chain_id: str
    sequence: str
    max_recycles: int = 3

    def sanity_check(self) -> list[str]:
        problems: list[str] = []
        if not self.sequence:
            problems.append("empty sequence")
        if any(c not in "ACDEFGHIKLMNPQRSTVWYX" for c in self.sequence.upper()):
            bad = sorted({c for c in self.sequence.upper()
                          if c not in "ACDEFGHIKLMNPQRSTVWYX"})
            problems.append(f"non-standard amino-acid letters: {bad}")
        if self.max_recycles < 1 or self.max_recycles > 12:
            problems.append(
                f"max_recycles={self.max_recycles} outside [1, 12] "
                "(ESMFold hyperparameter range)"
            )
        return problems


@dataclass
class EsmFoldAttempt:
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
    return resolve_binary_from_candidates(
        candidates or DEFAULT_BINARY_CANDIDATES,
    )


def _write_esmfold_manifest(
    output_dir: Path,
    spec: EsmFoldSequenceSpec,
) -> Path:
    manifest = {
        "generator": "esmfold_single_sequence",
        "sequence": {
            "chain_id": spec.chain_id,
            "sequence": spec.sequence,
            "length": len(spec.sequence),
        },
        "inference": {
            "max_recycles": spec.max_recycles,
        },
    }
    output_dir.mkdir(parents=True, exist_ok=True)
    path = output_dir / "esmfold_manifest.json"
    path.write_text(json.dumps(manifest, indent=2))
    return path


def attempt_single_sequence_fold(
    spec: EsmFoldSequenceSpec,
    output_dir: str | Path,
    n_candidates: int = 1,
    timeout_seconds: int = DEFAULT_TIMEOUT_S,
    binary_candidates: Iterable[str] | None = None,
) -> EsmFoldAttempt:
    """Try to run ESMFold on a single sequence and return the candidate PDBs.

    ESMFold is deterministic for a given sequence + recycle count, so
    ``n_candidates`` is normally 1 — the CLI is expected to write
    exactly one ``candidate_1.pdb`` (or similar).  When ``n_candidates``
    > 1, the CLI is assumed to vary the recycle-count / temperature to
    produce diverse conformers.

    Fail-open cascade delegated to _external_generator_runner:
    1. Invalid spec → ran=False.
    2. No binary on PATH → ran=False with "AF-splice" fallback advisory.
    3. Subprocess timeout / non-zero / no output → ran=True + accepted=False.
    4. Success → EsmFoldAttempt with candidate PDBs.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    problems = spec.sanity_check()
    if problems:
        return EsmFoldAttempt(
            ran=False, accepted=False,
            fallback_reason="sequence spec invalid: " + "; ".join(problems),
            diagnostics=list(problems),
        )

    binary = _resolve_binary(binary_candidates)
    if binary is None:
        return EsmFoldAttempt(
            ran=False, accepted=False,
            fallback_reason=(
                "no esmfold binary on PATH; tried "
                + ",".join(binary_candidates or DEFAULT_BINARY_CANDIDATES)
                + ". Falls back to AF-splice via the UniProt→AF cache."
            ),
        )

    manifest = _write_esmfold_manifest(output_dir, spec)

    run = execute_generator_subprocess(
        binary=binary,
        manifest_path=manifest,
        output_dir=output_dir,
        n_candidates=n_candidates,
        timeout_seconds=timeout_seconds,
        tool_name_for_messages="esmfold",
    )
    return EsmFoldAttempt(
        ran=run.ran, accepted=run.accepted,
        fallback_reason=run.fallback_reason,
        candidate_pdbs=list(run.candidate_pdbs),
        diagnostics=list(run.diagnostics),
        used_binary=run.used_binary,
    )


def summarise(attempt: EsmFoldAttempt) -> str:
    if not attempt.ran:
        return f"esmfold: skipped ({attempt.fallback_reason})"
    if not attempt.accepted:
        return f"esmfold: ran but rejected ({attempt.fallback_reason})"
    return (
        f"esmfold: {len(attempt.candidate_pdbs)} candidates via "
        f"{attempt.used_binary}"
    )
