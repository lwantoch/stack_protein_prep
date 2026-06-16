"""MODELLER-specific mechanics for the filler pipeline."""
from __future__ import annotations

import importlib.util
import os
import subprocess
import sys
from pathlib import Path
from typing import Any

from ._filler_shared import (
    DEFAULT_ALIGNMENT_FILENAME,
    DEFAULT_FINAL_MODEL_FILENAME,
    DEFAULT_MODELLER_SCRIPT_FILENAME,
    DEFAULT_SCORE_FILENAME,
    DEFAULT_STDERR_LOG_FILENAME,
    DEFAULT_STDOUT_LOG_FILENAME,
    _debug,
    cleanup_model_pdb,
)
from ._filler_analysis import (
    _extract_chain_id_from_alignment_header,
    _normalize_chain_id,
    build_modeller_template_alignment_sequence,
    extract_sequence_from_template_pdb,
    split_template_and_target_alignment_records,
    validate_template_sequence_consistency,
)

PYTHON_BIN = Path(sys.executable)


def _detect_modeller_pythonpath() -> str:
    override = os.environ.get("MODELLER_PYTHONPATH", "").strip()
    if override:
        return override
    if importlib.util.find_spec("modeller") is not None:
        return ""
    return "/usr/lib/modeller10.8/modlib:/usr/lib/modeller10.8/lib/x86_64-intel8"


MODELLER_PYTHONPATH = _detect_modeller_pythonpath()


def _write_skip_logs(output_dir: Path, message: str) -> tuple[Path, Path]:
    stdout_log = output_dir / DEFAULT_STDOUT_LOG_FILENAME
    stderr_log = output_dir / DEFAULT_STDERR_LOG_FILENAME
    stdout_log.write_text(message + "\n", encoding="utf-8")
    stderr_log.write_text("", encoding="utf-8")
    return stdout_log, stderr_log


def _strip_terminal_template_gaps(
    template_seq: str,
    target_seq: str,
) -> tuple[str, str]:
    """Remove alignment columns where the template has only terminal (leading or
    trailing) gaps.  Internal gaps are left intact so MODELLER can model them.
    This prevents MODELLER from receiving target residues that have no structural
    template at all (e.g., a signal peptide preceding the crystallised fragment)
    which would otherwise be built de-novo and corrupt residue numbering."""
    first = 0
    while first < len(template_seq) and template_seq[first] == "-":
        first += 1
    last = len(template_seq) - 1
    while last >= first and template_seq[last] == "-":
        last -= 1
    if first > last:
        raise ValueError("Template alignment sequence is entirely gaps; cannot build MODELLER alignment")
    return template_seq[first : last + 1], target_seq[first : last + 1]


def write_modeller_alignment_from_existing_alignment(
    alignment_fasta_path: Path,
    template_pdb_path: Path,
    output_dir: Path,
    template_id: str,
    target_id: str,
    chain_id: str,
) -> Path:
    (
        template_header,
        template_alignment_skeleton,
        target_header,
        target_aligned_sequence,
    ) = split_template_and_target_alignment_records(alignment_fasta_path)

    header_chain_id = _extract_chain_id_from_alignment_header(template_header)
    normalized_chain_id = _normalize_chain_id(chain_id)
    if header_chain_id is not None and header_chain_id != normalized_chain_id:
        raise ValueError(
            f"Alignment FASTA template header chain mismatch: "
            f"expected chain {normalized_chain_id!r}, "
            f"but header {template_header!r} suggests chain {header_chain_id!r}."
        )

    template_alignment_skeleton, target_aligned_sequence = _strip_terminal_template_gaps(
        template_alignment_skeleton, target_aligned_sequence
    )

    actual_template_sequence = extract_sequence_from_template_pdb(
        template_pdb_path=template_pdb_path,
    )
    validate_template_sequence_consistency(
        template_alignment_skeleton=template_alignment_skeleton,
        actual_template_sequence=actual_template_sequence,
        output_dir=output_dir,
        template_header=template_header,
        template_pdb_path=template_pdb_path,
    )
    rebuilt_template_aligned_sequence = build_modeller_template_alignment_sequence(
        aligned_template_skeleton=template_alignment_skeleton,
        actual_template_sequence=actual_template_sequence,
    )

    output_dir.mkdir(parents=True, exist_ok=True)
    alignment_file = output_dir / DEFAULT_ALIGNMENT_FILENAME
    content = (
        f">P1;{template_id}\n"
        f"structureX:{template_id}:FIRST:@:LAST:@::::\n"
        f"{rebuilt_template_aligned_sequence}*\n"
        f">P1;{target_id}\n"
        f"sequence:{target_id}:FIRST:@:LAST:@::::\n"
        f"{target_aligned_sequence}*\n"
    )
    alignment_file.write_text(content, encoding="utf-8")
    _debug(f"Wrote MODELLER alignment file to: {alignment_file}")
    return alignment_file


def _render_restraints_block(contact_atoms: list[Any]) -> str:
    """Return the Python snippet that adds position restraints to a MODELLER model."""
    if not contact_atoms:
        return ""
    rows = ", ".join(
        f"({at.atom_name!r}, {at.resseq}, {at.x}, {at.y}, {at.z}, {at.stdev})"
        for at in contact_atoms
    )
    return f"""
_CONTACT_RESTRAINTS = [{rows}]


class _RestrainedModel(automodel):
    def special_restraints(self, aln):
        rsr = self.restraints
        lookup = {{
            (name, resseq): (x, y, z, stdev)
            for name, resseq, x, y, z, stdev in _CONTACT_RESTRAINTS
        }}
        for at in self.atoms:
            key = (at.name.strip(), int(at.residue.num))
            if key not in lookup:
                continue
            x0, y0, z0, sd = lookup[key]
            rsr.add(forms.Gaussian(group=physical.xy_distance,
                feature=features.x_coordinate(at), mean=x0, stdev=sd))
            rsr.add(forms.Gaussian(group=physical.xy_distance,
                feature=features.y_coordinate(at), mean=y0, stdev=sd))
            rsr.add(forms.Gaussian(group=physical.xy_distance,
                feature=features.z_coordinate(at), mean=z0, stdev=sd))

"""


def write_modeller_script(
    output_dir: Path,
    alignment_file: Path,
    template_id: str,
    target_id: str,
    starting_model: int = 1,
    ending_model: int = 1,
    contact_atoms: list[Any] | None = None,
) -> Path:
    restraints_block = _render_restraints_block(contact_atoms or [])
    model_class = "_RestrainedModel" if restraints_block else "automodel"

    script_path = output_dir / DEFAULT_MODELLER_SCRIPT_FILENAME
    content = f"""from modeller import *
from modeller.automodel import *

log.verbose()

env = environ()
env.io.atom_files_directory = ['.']
env.io.hetatm = True
{restraints_block}
a = {model_class}(
    env,
    alnfile='{alignment_file.name}',
    knowns='{template_id}',
    sequence='{target_id}',
    assess_methods=(assess.DOPE, assess.GA341),
)

a.starting_model = {starting_model}
a.ending_model = {ending_model}

a.make()

score_file = '{DEFAULT_SCORE_FILENAME}'
with open(score_file, 'w', encoding='utf-8') as handle:
    handle.write("model_name\\tdope_score\\tga341_score\\n")
    for model in a.outputs:
        if model.get('failure') is None:
            model_name = model.get('name')
            dope_score = model.get('DOPE score')
            ga341_score = model.get('GA341 score')
            handle.write(f"{{model_name}}\\t{{dope_score}}\\t{{ga341_score}}\\n")
"""
    script_path.write_text(content, encoding="utf-8")
    return script_path


def validate_modeller_inputs(
    output_dir: Path,
    alignment_file: Path,
    template_id: str,
) -> None:
    expected_template_pdb = output_dir / f"{template_id}.pdb"
    if not alignment_file.exists():
        raise FileNotFoundError(f"MODELLER alignment file not found: {alignment_file}")
    if not expected_template_pdb.exists():
        raise FileNotFoundError(
            f"MODELLER template PDB file is missing. Expected: {expected_template_pdb}"
        )


def run_modeller_binary(
    script_path: Path,
    working_dir: Path,
) -> tuple[Path, Path]:
    if not PYTHON_BIN.exists():
        raise FileNotFoundError(f"Python binary not found: {PYTHON_BIN}")
    if not script_path.exists():
        raise FileNotFoundError(f"MODELLER script not found: {script_path}")
    command = [str(PYTHON_BIN), script_path.name]
    env = dict(os.environ)
    if MODELLER_PYTHONPATH:
        env["PYTHONPATH"] = MODELLER_PYTHONPATH
    result = subprocess.run(
        command, cwd=working_dir, capture_output=True, text=True, env=env
    )
    stdout_log = working_dir / DEFAULT_STDOUT_LOG_FILENAME
    stderr_log = working_dir / DEFAULT_STDERR_LOG_FILENAME
    stdout_log.write_text(result.stdout, encoding="utf-8")
    stderr_log.write_text(result.stderr, encoding="utf-8")
    if result.returncode != 0:
        raise RuntimeError(
            f"MODELLER failed for script {script_path.name}. "
            f"See {stdout_log} and {stderr_log}."
        )
    return stdout_log, stderr_log


def find_raw_models(output_dir: Path, target_id: str) -> tuple[Path, ...]:
    preferred_models = tuple(sorted(output_dir.glob(f"{target_id}.B*.pdb")))
    if preferred_models:
        return preferred_models
    return tuple(
        sorted(
            path
            for path in output_dir.glob(f"{target_id}*.pdb")
            if path.name != DEFAULT_FINAL_MODEL_FILENAME
        )
    )


def _parse_model_scores(score_file: Path) -> list[tuple[str, float]]:
    rows: list[tuple[str, float]] = []
    lines = score_file.read_text(encoding="utf-8").splitlines()
    if len(lines) < 2:
        return rows
    for line in lines[1:]:
        stripped = line.strip()
        if not stripped:
            continue
        parts = stripped.split("\t")
        if len(parts) != 3:
            continue
        model_name, dope_score_text, _ga341_score_text = parts
        if dope_score_text in {"None", "", "nan"}:
            continue
        try:
            dope_score = float(dope_score_text)
        except ValueError:
            continue
        rows.append((model_name, dope_score))
    return rows


def select_best_model(
    output_dir: Path,
    target_id: str,
    raw_model_paths: tuple[Path, ...],
) -> Path:
    score_file = output_dir / DEFAULT_SCORE_FILENAME
    if score_file.exists():
        score_rows = _parse_model_scores(score_file)
        if score_rows:
            best_model_name, best_dope_score = min(score_rows, key=lambda row: row[1])
            best_model_path = output_dir / best_model_name
            if best_model_path.exists():
                _debug(
                    f"Best model selected by lowest DOPE score: "
                    f"{best_model_path} (DOPE={best_dope_score})"
                )
                return best_model_path

    if len(raw_model_paths) == 1:
        return raw_model_paths[0]
    if raw_model_paths:
        return raw_model_paths[0]

    raise FileNotFoundError(
        f"No raw MODELLER models found for target {target_id!r} in {output_dir}"
    )


def select_best_model_from_scores(output_dir: Path) -> Path:
    """Return the path to the best model based solely on the score file (lowest DOPE)."""
    score_file = output_dir / DEFAULT_SCORE_FILENAME
    if not score_file.exists():
        raise FileNotFoundError(f"Score file not found: {score_file}")
    score_rows = _parse_model_scores(score_file)
    if not score_rows:
        raise ValueError(f"No valid score rows found in {score_file}")
    best_model_name, _best_dope = min(score_rows, key=lambda row: row[1])
    best_model_path = output_dir / best_model_name
    if not best_model_path.exists():
        raise FileNotFoundError(f"Best model file not found: {best_model_path}")
    return best_model_path


def standardize_model_name(
    output_dir: Path,
    target_id: str,
    raw_model_paths: tuple[Path, ...],
    final_name: str,
) -> tuple[Path, Path]:
    best_model = select_best_model(
        output_dir=output_dir,
        target_id=target_id,
        raw_model_paths=raw_model_paths,
    )
    final_model = output_dir / final_name
    cleaned_final_model = cleanup_model_pdb(
        input_model_path=best_model,
        output_model_path=final_model,
    )
    return best_model, cleaned_final_model
