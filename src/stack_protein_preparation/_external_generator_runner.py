"""Shared subprocess-cascade helper for external generator scaffolds.

Both _rfdiffusion2_gap.attempt_gap_fill and
_boltz2_template_steering.attempt_template_steering follow the exact
same 5-layer fail-open cascade once the input has been validated:

    1. Resolve a binary from a candidate name list (or fail-open).
    2. Launch the binary as a subprocess with a JSON manifest arg +
       an output-dir arg.
    3. Handle TimeoutExpired → ran=True, accepted=False.
    4. Handle OSError (spawn failure) → ran=False, accepted=False.
    5. Handle non-zero exit → ran=True, accepted=False, stderr in msg.
    6. Handle zero exit but no candidate_*.pdb files written →
       ran=True, accepted=False.
    7. Success → collect candidate_*.pdb files.

This module extracts (1) + (2)-(7) into shared helpers so the two
scaffolds carry only their specialised Spec + manifest writer.  When
another generator (Boltz-3, RFdiffusion3, AlphaFold-Multimer,
ESMFold) needs to be added, it can wire straight into this helper.

Pure Python + stdlib.  License-free.
"""
from __future__ import annotations

import shutil
import subprocess
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable


@dataclass
class GeneratorRunResult:
    """Uniform result the runner returns to the scaffold caller."""
    ran: bool
    accepted: bool
    fallback_reason: str | None = None
    candidate_pdbs: list[Path] = field(default_factory=list)
    diagnostics: list[str] = field(default_factory=list)
    used_binary: str | None = None


def resolve_binary_from_candidates(candidates: Iterable[str]) -> str | None:
    """First-on-PATH from ``candidates``; None if all missing."""
    for name in candidates:
        found = shutil.which(name)
        if found:
            return found
    return None


def execute_generator_subprocess(
    binary: str,
    manifest_path: Path,
    output_dir: Path,
    n_candidates: int,
    timeout_seconds: int,
    tool_name_for_messages: str,
    candidate_glob: str = "candidate_*.pdb",
) -> GeneratorRunResult:
    """Launch the generator binary with ``--manifest``, ``--output-dir``,
    ``--n-candidates`` args and cascade over the standard failure modes.

    Args:
        binary: absolute path to the generator binary (from
            resolve_binary_from_candidates).
        manifest_path: JSON manifest the generator will consume.
        output_dir: where candidate PDBs are expected to appear.
        n_candidates: how many samples to request.
        timeout_seconds: hard wall-clock cap.
        tool_name_for_messages: short label (e.g. "rfdiffusion",
            "boltz") used in fallback_reason strings.
        candidate_glob: pattern to enumerate produced candidate PDBs.

    Returns a GeneratorRunResult with used_binary always set to
    ``binary`` (even on failure) so the caller can log which binary
    was invoked.
    """
    cmd = [
        binary,
        "--manifest", str(manifest_path),
        "--output-dir", str(output_dir),
        "--n-candidates", str(n_candidates),
    ]
    try:
        proc = subprocess.run(
            cmd, capture_output=True, text=True,
            timeout=timeout_seconds,
        )
    except subprocess.TimeoutExpired:
        return GeneratorRunResult(
            ran=True, accepted=False,
            fallback_reason=(
                f"{tool_name_for_messages} timed out after {timeout_seconds}s"
            ),
            used_binary=binary,
        )
    except OSError as exc:
        return GeneratorRunResult(
            ran=False, accepted=False,
            fallback_reason=f"failed to spawn {binary}: {exc!r}",
            used_binary=binary,
        )

    if proc.returncode != 0:
        return GeneratorRunResult(
            ran=True, accepted=False,
            fallback_reason=(
                f"{tool_name_for_messages} exit={proc.returncode}; "
                f"stderr={proc.stderr.strip()[:300]}"
            ),
            used_binary=binary,
            diagnostics=[proc.stdout[-500:], proc.stderr[-500:]],
        )

    candidates = sorted(output_dir.glob(candidate_glob))
    if not candidates:
        return GeneratorRunResult(
            ran=True, accepted=False,
            fallback_reason=(
                f"{tool_name_for_messages} ran but produced no "
                f"{candidate_glob}"
            ),
            used_binary=binary,
            diagnostics=[proc.stdout[-500:]],
        )

    return GeneratorRunResult(
        ran=True, accepted=True,
        candidate_pdbs=candidates,
        used_binary=binary,
        diagnostics=[
            f"n_candidates={len(candidates)}",
            f"binary={binary}",
        ],
    )
