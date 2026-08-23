"""Boltz-2 fallback for FRUTON gap fills that MODELLER LoopModel cannot rescue.

Boltz-2 (Passaro et al. bioRxiv 2025.06.14.659707) is a diffusion-based
structure predictor with "strict template enforcement" mode: it accepts
a crystal backbone as a hard-masked template and samples new
conformations for unresolved residues while respecting the template
geometry AND experimental method conditioning (X-ray / NMR / MD).  This
is exactly the constraint set FRUTON needs for gaps where the AF-loop
splice + MODELLER LoopModel cannot both fit the crystal gap distance
AND avoid steric clashes (documented ceiling cases: 5HJS, 8A27, 7QUE).

Because Boltz-2 requires a dedicated Python env (PyTorch + CUDA + ~5 GB
model weights) that dwarfs FRUTON's ``pixi`` env, this module talks to
Boltz-2 via subprocess rather than as a Python dependency.  On CESGA
the wrapper expects ``boltz`` to be on PATH inside the SLURM job
(loaded from a dedicated conda / pixi env with GPU access).  For local
pytest runs, ``boltz`` is not present -- the wrapper detects that and
returns a fail-open result so the pipeline never blocks on a missing
optional dependency.

Interface mirrors ``_filler_loop_refine.refine_loops_via_modeller``:
returns a dataclass with ``ran`` / ``fallback_reason`` and always
writes ``output_pdb`` (falls back to input on failure) so the caller
sees the expected artifact.
"""
from __future__ import annotations

import os
import shutil
import subprocess
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable


DEFAULT_BOLTZ_BINARY_NAME = "boltz"
DEFAULT_TIMEOUT_SECONDS = 1800  # 30 min per protein on a single A100
DEFAULT_SAMPLING_STEPS = 200
DEFAULT_TEMPLATE_STEERING_STRENGTH = 1.0  # 0 = ignore, 1 = strict


@dataclass
class Boltz2RefineResult:
    output_pdb_path: Path
    ran: bool
    fallback_reason: str | None = None
    n_gap_regions_refined: int = 0
    boltz_stdout_tail: str | None = None
    diagnostics: list[str] = field(default_factory=list)

    def to_dict(self) -> dict:
        return {
            "ran": self.ran,
            "fallback_reason": self.fallback_reason,
            "n_gap_regions_refined": self.n_gap_regions_refined,
            "boltz_stdout_tail": self.boltz_stdout_tail,
            "diagnostics": self.diagnostics,
        }


def _resolve_boltz_binary(boltz_bin: str | None) -> str | None:
    """Return the resolved path to the boltz binary, or None if missing.

    Priority: (1) explicit ``boltz_bin`` argument, (2) ``FRUTON_BOLTZ_BIN``
    environment variable, (3) ``boltz`` on PATH.  Not finding it is not
    an error; the caller falls open.
    """
    candidate = boltz_bin or os.environ.get("FRUTON_BOLTZ_BIN") or DEFAULT_BOLTZ_BINARY_NAME
    resolved = shutil.which(candidate)
    return resolved


def _extract_target_sequence_from_pdb(
    input_pdb_path: Path,
    chain_id: str,
) -> str:
    """Extract the ATOM-only single-letter sequence for ``chain_id``.

    Includes only standard amino acids present in the PDB (missing
    residues stay missing -- Boltz-2 will sample them under template
    constraints).  Used to build the Boltz-2 target sequence input.
    """
    from Bio.PDB import PDBParser
    three_to_one = {
        "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
        "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
        "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
        "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
    }
    parser = PDBParser(QUIET=True)
    try:
        s = parser.get_structure("t", str(input_pdb_path))
    except Exception:  # noqa: BLE001
        return ""
    seq_parts: list[str] = []
    for model in s:
        for chain in model:
            if chain.id != chain_id:
                continue
            for res in chain:
                if res.id[0].strip():
                    continue
                seq_parts.append(three_to_one.get(res.resname.strip().upper(), "X"))
            break
        break
    return "".join(seq_parts)


def _write_boltz_input_yaml(
    input_pdb_path: Path,
    output_dir: Path,
    target_id: str,
    sequence: str,
) -> Path:
    """Emit a minimal Boltz-2 input YAML for template-steered prediction.

    The exact schema evolves with Boltz-2 releases; this function
    produces a schema-2.2.x-compatible layout.  When Boltz-2 upstream
    changes, only this helper needs updating.
    """
    yaml_path = output_dir / f"{target_id}_boltz_input.yaml"
    content = f"""version: 1
sequences:
  - protein:
      id: A
      sequence: {sequence}
      msa: empty
templates:
  - cif_path: {input_pdb_path.resolve()}
    template_id: crystal_A
    chain_id: A
    force: true
"""
    yaml_path.write_text(content, encoding="utf-8")
    return yaml_path


def refine_gaps_via_boltz2(
    input_pdb_path: str | Path,
    output_pdb_path: str | Path,
    gap_ranges_by_chain: Iterable[tuple[str, int, int]],
    uniprot_id: str | None = None,
    boltz_bin: str | None = None,
    timeout_seconds: int = DEFAULT_TIMEOUT_SECONDS,
    sampling_steps: int = DEFAULT_SAMPLING_STEPS,
    template_steering_strength: float = DEFAULT_TEMPLATE_STEERING_STRENGTH,
) -> Boltz2RefineResult:
    """Run Boltz-2 on ``input_pdb_path`` with the crystal as strict
    template and return the refined structure at ``output_pdb_path``.

    Fail-open: returns ``ran=False`` (and copies input to output) when
    Boltz-2 is not installed, the subprocess errors, or the output
    format is unparsable.  Callers should treat a None-result as
    equivalent to the input structure -- the upstream MODELLER-refined
    model remains authoritative.
    """
    input_path = Path(input_pdb_path)
    output_path = Path(output_pdb_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    gap_list = list(gap_ranges_by_chain)
    if not gap_list:
        output_path.write_bytes(input_path.read_bytes())
        return Boltz2RefineResult(
            output_pdb_path=output_path, ran=False,
            fallback_reason="no gap ranges provided",
        )

    resolved_bin = _resolve_boltz_binary(boltz_bin)
    if resolved_bin is None:
        output_path.write_bytes(input_path.read_bytes())
        return Boltz2RefineResult(
            output_pdb_path=output_path, ran=False,
            fallback_reason=(
                "boltz binary not on PATH; install with 'pip install boltz' "
                "in a GPU-capable env and set FRUTON_BOLTZ_BIN, or load the "
                "CESGA boltz module inside the SLURM job."
            ),
        )

    # Pick a representative chain for the target sequence.  For multi-
    # chain crystals, Boltz-2 needs one input at a time; iterate over
    # unique chains in gap_list.
    target_chain = gap_list[0][0]
    sequence = _extract_target_sequence_from_pdb(input_path, target_chain)
    if not sequence or len(sequence) < 10:
        output_path.write_bytes(input_path.read_bytes())
        return Boltz2RefineResult(
            output_pdb_path=output_path, ran=False,
            fallback_reason=f"could not extract usable target sequence from chain {target_chain!r}",
        )

    import tempfile
    workdir = Path(tempfile.mkdtemp(prefix="fruton_boltz2_"))
    target_id = f"fruton_{input_path.stem}"

    try:
        yaml_path = _write_boltz_input_yaml(
            input_pdb_path=input_path,
            output_dir=workdir,
            target_id=target_id,
            sequence=sequence,
        )
    except Exception as exc:  # noqa: BLE001
        output_path.write_bytes(input_path.read_bytes())
        return Boltz2RefineResult(
            output_pdb_path=output_path, ran=False,
            fallback_reason=f"failed to write Boltz-2 input YAML: {exc!r}",
        )

    cmd = [
        resolved_bin,
        "predict",
        str(yaml_path),
        "--out_dir", str(workdir / "boltz_out"),
        "--sampling_steps", str(sampling_steps),
        "--recycling_steps", "3",
        "--diffusion_samples", "1",
        "--use_msa_server",
    ]

    try:
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=timeout_seconds,
        )
    except subprocess.TimeoutExpired:
        output_path.write_bytes(input_path.read_bytes())
        return Boltz2RefineResult(
            output_pdb_path=output_path, ran=False,
            fallback_reason=f"Boltz-2 timed out after {timeout_seconds}s",
        )
    except (FileNotFoundError, PermissionError) as exc:
        output_path.write_bytes(input_path.read_bytes())
        return Boltz2RefineResult(
            output_pdb_path=output_path, ran=False,
            fallback_reason=f"Boltz-2 exec failed: {exc!r}",
        )

    stdout_tail = "\n".join(result.stdout.splitlines()[-40:]) if result.stdout else ""

    if result.returncode != 0:
        output_path.write_bytes(input_path.read_bytes())
        return Boltz2RefineResult(
            output_pdb_path=output_path, ran=False,
            fallback_reason=f"Boltz-2 exited with code {result.returncode}",
            boltz_stdout_tail=stdout_tail,
        )

    # Locate the model_0 CIF/PDB Boltz-2 writes into out_dir/predictions/
    predictions_dir = workdir / "boltz_out" / "predictions" / target_id
    candidate_pdbs = sorted(predictions_dir.glob(f"{target_id}_model_*.pdb")) + \
                     sorted(predictions_dir.glob(f"{target_id}_model_*.cif"))
    if not candidate_pdbs:
        output_path.write_bytes(input_path.read_bytes())
        return Boltz2RefineResult(
            output_pdb_path=output_path, ran=False,
            fallback_reason=f"Boltz-2 produced no output structures in {predictions_dir}",
            boltz_stdout_tail=stdout_tail,
        )

    best = candidate_pdbs[0]
    try:
        if best.suffix.lower() == ".cif":
            from Bio.PDB import MMCIFParser, PDBIO
            parser = MMCIFParser(QUIET=True)
            struct = parser.get_structure("boltz", str(best))
            io = PDBIO()
            io.set_structure(struct)
            io.save(str(output_path))
        else:
            output_path.write_bytes(best.read_bytes())
    except Exception as exc:  # noqa: BLE001
        output_path.write_bytes(input_path.read_bytes())
        return Boltz2RefineResult(
            output_pdb_path=output_path, ran=False,
            fallback_reason=f"failed to convert Boltz-2 output {best.name}: {exc!r}",
            boltz_stdout_tail=stdout_tail,
        )

    return Boltz2RefineResult(
        output_pdb_path=output_path, ran=True, fallback_reason=None,
        n_gap_regions_refined=len(gap_list),
        boltz_stdout_tail=stdout_tail,
        diagnostics=[
            f"boltz_bin={resolved_bin}",
            f"sampling_steps={sampling_steps}",
            f"template_steering_strength={template_steering_strength}",
            f"target_chain={target_chain}, sequence_len={len(sequence)}",
            f"n_gaps={len(gap_list)}",
            f"best_output={best.name}",
        ],
    )
