from __future__ import annotations

import importlib.util
import json
import os
import re
import subprocess
import sys
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any
from urllib.error import HTTPError, URLError
from urllib.request import urlopen, urlretrieve

from Bio.Data.PDBData import protein_letters_3to1_extended
from Bio.PDB import PDBIO, PDBParser, Superimposer
from Bio.PDB.PDBIO import Select
from Bio.PDB.Polypeptide import is_aa

from stack_protein_preparation.pipeline_logging import (
    append_exception_traceback,
    append_key_value_block,
    append_log_text,
    build_module_log_paths,
    print_double_header_box,
    print_light_table_box,
    print_mixed_box,
    write_log_header,
)

PYTHON_BIN = Path(sys.executable)


def _detect_modeller_pythonpath() -> str:
    """Return PYTHONPATH entries needed to import modeller.

    Priority:
    1. MODELLER_PYTHONPATH env var — explicit override (useful for custom installs).
    2. If modeller is already importable in this Python (conda/pixi install) — no extra path needed.
    3. Legacy local fallback (/usr/lib/modeller10.8) for dev machines with a system install.
    """
    override = os.environ.get("MODELLER_PYTHONPATH", "").strip()
    if override:
        return override
    if importlib.util.find_spec("modeller") is not None:
        return ""  # already on sys.path via conda/pixi
    return "/usr/lib/modeller10.8/modlib:/usr/lib/modeller10.8/lib/x86_64-intel8"


MODELLER_PYTHONPATH = _detect_modeller_pythonpath()

MODULE_NAME = "filler"

DEFAULT_FINAL_MODEL_FILENAME = "final_filled_model.pdb"
DEFAULT_ALIGNMENT_FILENAME = "chain_alignment.ali"
DEFAULT_MODELLER_SCRIPT_FILENAME = "run_modeller.py"
DEFAULT_SCORE_FILENAME = "model_scores.tsv"
DEFAULT_STDOUT_LOG_FILENAME = "modeller_stdout.log"
DEFAULT_STDERR_LOG_FILENAME = "modeller_stderr.log"
DEFAULT_MANIFEST_FILENAME = "filler_manifest.json"
DEFAULT_ALPHAFOLD_DIRNAME = "alphafold"
DEFAULT_ALPHAFOLD_RAW_FILENAME = "alphafold_downloaded_model.pdb"
DEFAULT_ALPHAFOLD_CROPPED_FILENAME = "alphafold_cropped_model.pdb"
DEFAULT_ALPHAFOLD_ALIGNED_FILENAME = "alphafold_aligned_model.pdb"

FORCE_FIELD_RESIDUE_ALIASES = {
    "ASH": "ASP",
    "GLH": "GLU",
    "HID": "HIS",
    "HIE": "HIS",
    "HIP": "HIS",
    "CYM": "CYS",
    "CYX": "CYS",
    "LYN": "LYS",
}


@dataclass(frozen=True)
class GapRegion:
    alignment_start: int
    alignment_end: int
    gap_length: int
    is_terminal: bool
    classification: str


@dataclass(frozen=True)
class FillDecision:
    should_run_modeller: bool
    overall_classification: str
    gap_regions: tuple[GapRegion, ...]
    skip_reason: str | None
    alphafold_candidate: bool


@dataclass(frozen=True)
class FillerRunResult:
    chain_id: str
    output_dir: Path
    alignment_file: Path
    template_pdb: Path
    script_file: Path
    modeller_model_path: Path | None
    alphafold_model_path: Path | None
    final_model_path: Path | None
    raw_model_paths: tuple[Path, ...]
    stdout_log: Path
    stderr_log: Path
    skipped: bool
    skip_reason: str | None
    fill_decision: FillDecision


@dataclass(frozen=True)
class _ParsedFirstModel:
    """Store a BioPython structure and its first model.

    Filler operations are restricted to the first model because MODELLER,
    AlphaFold alignment, and the downstream GROMACS route expect one coordinate
    model. Keeping the structure and first model together avoids reparsing the
    same PDB path in related helpers. This is an internal container because
    callers should not depend on BioPython object details from this module.
    """

    structure: Any
    model: Any


class _FirstModelChainProteinSelect(Select):
    """Select one protein-like chain from the first model.

    This selector is used for chain-specific MODELLER template writing. It keeps
    only the requested chain and only protein-like residues recognized by the
    module's BioPython-backed residue logic. Ligands, waters, metals, and other
    HETATM records are deliberately excluded from the MODELLER template.
    """

    def __init__(self, *, model_id: Any, chain_id: str) -> None:
        self._model_id = model_id
        self._chain_id = _normalize_chain_id(chain_id)

    def accept_model(self, model: Any) -> bool:
        return model.id == self._model_id

    def accept_chain(self, chain: Any) -> bool:
        return _normalize_chain_id(str(chain.id)) == self._chain_id

    def accept_residue(self, residue: Any) -> bool:
        return _is_protein_atom_residue(residue)


class _FirstModelProteinRangeSelect(Select):
    """Select protein-like residues inside a residue-number interval.

    AlphaFold models are cropped to the requested residue range before they are
    aligned to the experimental template. This selector keeps the cropping logic
    in BioPython instead of copying PDB text columns manually. It assumes that
    the caller-provided residue range uses the same numbering convention as the
    downloaded AlphaFold PDB.
    """

    def __init__(
        self,
        *,
        model_id: Any,
        start_residue: int,
        end_residue: int,
    ) -> None:
        self._model_id = model_id
        self._start_residue = start_residue
        self._end_residue = end_residue

    def accept_model(self, model: Any) -> bool:
        return model.id == self._model_id

    def accept_residue(self, residue: Any) -> bool:
        if not _is_protein_atom_residue(residue):
            return False

        residue_number = int(residue.id[1])
        return self._start_residue <= residue_number <= self._end_residue


class _FirstModelProteinOnlySelect(Select):
    """Select protein-like residues from the first model.

    This selector is used to clean MODELLER or AlphaFold-derived final models
    before the model path is returned to the runner. It removes HETATM records
    and keeps only supported protein-like ATOM residues. This preserves the
    previous cleanup intent while delegating PDB parsing and writing to
    BioPython.
    """

    def __init__(self, *, model_id: Any) -> None:
        self._model_id = model_id

    def accept_model(self, model: Any) -> bool:
        return model.id == self._model_id

    def accept_residue(self, residue: Any) -> bool:
        return _is_protein_atom_residue(residue)


def _debug(message: str) -> None:
    print(f"[filler] {message}")


def _normalize_chain_id(chain_id: str) -> str:
    return chain_id.strip() or "_"


def _normalize_residue_name_for_sequence(residue_name: str) -> str:
    """Normalize workflow-specific residue aliases before sequence conversion.

    BioPython provides the standard and extended PDB residue mapping, so this
    module should not maintain a full local three-letter to one-letter table.
    The only local aliases kept here are force-field/protonation names that are
    common in prepared structures but should be interpreted as their parent
    amino-acid residues for sequence and filler decisions. For example, ``GLH``
    is treated as glutamate and maps through ``GLU``.
    """

    normalized_name = residue_name.strip().upper()
    return FORCE_FIELD_RESIDUE_ALIASES.get(normalized_name, normalized_name)


def _residue_name_to_one_letter(residue_name: str) -> str:
    """Convert a residue name to one-letter code through BioPython.

    The conversion first applies the small force-field alias map and then uses
    BioPython's extended PDB residue table. Unknown or unsupported residue names
    become ``X`` so diagnostic sequence extraction can continue. Filtering of
    whether a residue is protein-like is handled separately.
    """

    normalized_name = _normalize_residue_name_for_sequence(residue_name)
    return protein_letters_3to1_extended.get(normalized_name, "X")


def _is_supported_protein_residue_name(residue_name: str) -> bool:
    """Return whether a residue name is protein-like for filler logic."""

    normalized_name = _normalize_residue_name_for_sequence(residue_name)
    return normalized_name in protein_letters_3to1_extended


def _is_protein_atom_residue(residue: Any) -> bool:
    """Return whether a BioPython residue should be treated as protein ATOM data.

    The filler template sequence is based on ATOM protein residues only, not
    ligands or generic HETATM records. BioPython's ``is_aa(..., standard=False)``
    covers standard and many modified amino acids. The small alias fallback keeps
    force-field protonation names usable without reintroducing a full local
    amino-acid table.
    """

    if residue.id[0] != " ":
        return False

    residue_name = residue.get_resname().strip().upper()
    return is_aa(residue, standard=False) or _is_supported_protein_residue_name(
        residue_name
    )


def _residue_has_peptide_backbone(residue: Any) -> bool:
    atom_name_set = {str(atom.id).strip().upper() for atom in residue}
    return {"N", "CA", "C"}.issubset(atom_name_set)


def _residue_to_one_letter(residue: Any) -> str:
    return _residue_name_to_one_letter(residue.get_resname())


def _parse_first_model(input_pdb_path: Path) -> _ParsedFirstModel:
    if not input_pdb_path.exists():
        raise FileNotFoundError(f"PDB file not found: {input_pdb_path}")

    parser = PDBParser(QUIET=True)
    structure = parser.get_structure(input_pdb_path.stem, str(input_pdb_path))
    models = list(structure.get_models())

    if not models:
        raise ValueError(f"No coordinate models were found in {input_pdb_path}")

    return _ParsedFirstModel(structure=structure, model=models[0])


def _count_atoms_in_first_model(input_pdb_path: Path) -> int:
    if not input_pdb_path.exists():
        return 0

    try:
        parsed = _parse_first_model(input_pdb_path)
    except Exception:
        return 0

    atom_count = 0
    for chain in parsed.model:
        for residue in chain:
            atom_count += len(list(residue.get_atoms()))

    return atom_count


def _validate_written_pdb_has_atoms(output_pdb_path: Path, description: str) -> None:
    if not output_pdb_path.exists() or output_pdb_path.stat().st_size == 0:
        raise ValueError(f"{description} was not written: {output_pdb_path}")

    if _count_atoms_in_first_model(output_pdb_path) == 0:
        output_pdb_path.unlink(missing_ok=True)
        raise ValueError(f"{description} contains no atoms: {output_pdb_path}")


def _classify_gap_length(gap_length: int) -> str:
    if 1 <= gap_length <= 5:
        return "green"
    if 6 <= gap_length < 8:
        return "yellow"
    return "alphafold_candidate"


def _safe_str_path(path: Path | None) -> str | None:
    if path is None:
        return None
    return str(path)


def _infer_pdb_id_from_target_id(target_id: str, template_pdb_path: Path) -> str:
    target_id = str(target_id).strip()
    if "_" in target_id:
        return target_id.split("_", maxsplit=1)[0].strip().upper()

    stem = template_pdb_path.stem.strip().upper()
    if stem:
        return stem.split("_", maxsplit=1)[0]

    raise ValueError(
        f"Could not infer pdb_id from target_id={target_id!r} "
        f"and template_pdb_path={template_pdb_path}"
    )


def _infer_protein_data_dir(template_pdb_path: Path) -> Path:
    """Infer the ``data/proteins`` directory from a template PDB path.

    The template path usually points either to a file inside
    ``data/proteins/<PDBID>/components`` or directly inside
    ``data/proteins/<PDBID>``. The filler module needs this parent directory to
    place module logs in the shared logging tree. This helper keeps that path
    inference local and does not expose it as part of the public API.
    """

    template_pdb_path = template_pdb_path.resolve()
    protein_dir = template_pdb_path.parent

    if protein_dir.name == "components":
        protein_dir = protein_dir.parent

    return protein_dir.parent


def _build_module_log_path(
    *,
    pdb_id: str,
    template_pdb_path: Path,
) -> Path:
    protein_data_dir = _infer_protein_data_dir(template_pdb_path)
    log_paths = build_module_log_paths(
        protein_data_dir=protein_data_dir,
        pdb_id=pdb_id,
        module_name=MODULE_NAME,
        variant_label=None,
    )
    return log_paths.module_log_path


def _log_debug(log_path: Path | None, message: str) -> None:
    _debug(message)
    if log_path is not None:
        append_log_text(log_path, f"[DEBUG] {message}")


def _print_start_box(
    *,
    pdb_id: str,
    chain_id: str | None,
    template_id: str,
    target_id: str,
    alignment_directory: Path,
    template_pdb_path: Path,
    output_dir: Path,
    residue_range: str,
    uniprot_id: str | None,
    log_path: Path,
) -> None:
    print_mixed_box(
        title=f"{MODULE_NAME} · START",
        line_list=[
            f"PDB ID         : {pdb_id}",
            f"Chain input    : {chain_id if chain_id is not None else '(auto)'}",
            f"Template ID    : {template_id}",
            f"Target ID      : {target_id}",
            f"Residue range  : {residue_range or '(infer full observed range)'}",
            f"UniProt ID     : {uniprot_id or '(none)'}",
            f"Alignment dir  : {alignment_directory}",
            f"Template PDB   : {template_pdb_path}",
            f"Output dir     : {output_dir}",
            f"Log file       : {log_path}",
        ],
    )


def _print_decision_table(
    *,
    pdb_id: str,
    resolved_chain_id: str,
    fill_decision: FillDecision,
) -> None:
    row_list = [
        ("PDB ID", pdb_id),
        ("Resolved chain", resolved_chain_id),
        ("should_run_modeller", str(fill_decision.should_run_modeller)),
        ("overall_classification", fill_decision.overall_classification),
        ("alphafold_candidate", str(fill_decision.alphafold_candidate)),
        ("n_gap_regions", str(len(fill_decision.gap_regions))),
        ("skip_reason", fill_decision.skip_reason or "(none)"),
    ]

    print_light_table_box(
        title=f"{MODULE_NAME} · decision",
        row_list=row_list,
        left_header="Field",
        right_header="Value",
    )


def _print_gap_regions_table(fill_decision: FillDecision) -> None:
    if not fill_decision.gap_regions:
        print_light_table_box(
            title=f"{MODULE_NAME} · gap regions",
            row_list=[("result", "no internal gap regions")],
            left_header="Field",
            right_header="Value",
        )
        return

    row_list = []
    for index, gap_region in enumerate(fill_decision.gap_regions, start=1):
        row_list.append(
            (
                f"gap_{index}",
                (
                    f"{gap_region.alignment_start}-{gap_region.alignment_end} | "
                    f"len={gap_region.gap_length} | "
                    f"terminal={gap_region.is_terminal} | "
                    f"class={gap_region.classification}"
                ),
            )
        )

    print_light_table_box(
        title=f"{MODULE_NAME} · gap regions",
        row_list=row_list,
        left_header="Region",
        right_header="Details",
    )


def _print_end_box(
    *,
    pdb_id: str,
    resolved_chain_id: str,
    status: str,
    final_model_path: Path | None,
    modeller_model_path: Path | None,
    alphafold_model_path: Path | None,
    stdout_log: Path,
    stderr_log: Path,
    skip_reason: str | None,
) -> None:
    print_double_header_box(
        title=f"{MODULE_NAME} · END",
        line_list=[
            f"PDB ID           : {pdb_id}",
            f"Resolved chain   : {resolved_chain_id}",
            f"Status           : {status}",
            f"Final model      : {final_model_path if final_model_path else '(none)'}",
            f"MODELLER model   : {modeller_model_path if modeller_model_path else '(none)'}",
            f"AlphaFold model  : {alphafold_model_path if alphafold_model_path else '(none)'}",
            f"stdout log       : {stdout_log}",
            f"stderr log       : {stderr_log}",
            f"skip_reason      : {skip_reason or '(none)'}",
        ],
    )


def _result_to_manifest_payload(result: FillerRunResult) -> dict[str, Any]:
    return {
        "chain_id": result.chain_id,
        "output_dir": str(result.output_dir),
        "alignment_file": str(result.alignment_file),
        "template_pdb": str(result.template_pdb),
        "script_file": str(result.script_file),
        "modeller_model_path": _safe_str_path(result.modeller_model_path),
        "alphafold_model_path": _safe_str_path(result.alphafold_model_path),
        "final_model_path": _safe_str_path(result.final_model_path),
        "raw_model_paths": [str(path) for path in result.raw_model_paths],
        "stdout_log": str(result.stdout_log),
        "stderr_log": str(result.stderr_log),
        "skipped": result.skipped,
        "skip_reason": result.skip_reason,
        "fill_decision": {
            "should_run_modeller": result.fill_decision.should_run_modeller,
            "overall_classification": result.fill_decision.overall_classification,
            "gap_regions": [
                asdict(region) for region in result.fill_decision.gap_regions
            ],
            "skip_reason": result.fill_decision.skip_reason,
            "alphafold_candidate": result.fill_decision.alphafold_candidate,
        },
    }


def _write_manifest(output_dir: Path, result: FillerRunResult) -> Path:
    manifest_path = output_dir / DEFAULT_MANIFEST_FILENAME
    payload = _result_to_manifest_payload(result)
    manifest_path.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    _debug(f"Wrote filler manifest: {manifest_path}")
    return manifest_path


def _build_result(
    *,
    chain_id: str,
    output_dir: Path,
    alignment_file: Path,
    template_pdb: Path,
    script_file: Path,
    modeller_model_path: Path | None,
    alphafold_model_path: Path | None,
    final_model_path: Path | None,
    raw_model_paths: tuple[Path, ...],
    stdout_log: Path,
    stderr_log: Path,
    skipped: bool,
    skip_reason: str | None,
    fill_decision: FillDecision,
) -> FillerRunResult:
    result = FillerRunResult(
        chain_id=chain_id,
        output_dir=output_dir,
        alignment_file=alignment_file,
        template_pdb=template_pdb,
        script_file=script_file,
        modeller_model_path=modeller_model_path,
        alphafold_model_path=alphafold_model_path,
        final_model_path=final_model_path,
        raw_model_paths=raw_model_paths,
        stdout_log=stdout_log,
        stderr_log=stderr_log,
        skipped=skipped,
        skip_reason=skip_reason,
        fill_decision=fill_decision,
    )
    _write_manifest(output_dir=output_dir, result=result)
    return result


def _strip_gaps(sequence: str) -> str:
    return sequence.replace("-", "")


def _normalize_alignment_template_for_validation(
    template_alignment_skeleton: str,
) -> str:
    """Normalize the template row before comparing it to the PDB sequence.

    The current alignment semantics are that ``-`` means an alignment gap or a
    missing region, while ``X`` means an unknown residue identity. Older
    alignments can still contain terminal placeholder ``X`` characters outside
    the actual template sequence. Those terminal placeholders are tolerated for
    validation only so they do not block AlphaFold or MODELLER routing.
    """

    ungapped = _strip_gaps(template_alignment_skeleton)
    ungapped = re.sub(r"^X+", "", ungapped)
    ungapped = re.sub(r"X+$", "", ungapped)
    return ungapped


def _parse_residue_range(residue_range: str) -> tuple[int, int]:
    """Parse a residue range string; returns (start, end) integers.

    Accepts both the legacy '25-200' format and the chain-aware 'A25-B200'
    format. For chain-aware ranges only the residue numbers are returned.
    """
    stripped = str(residue_range).strip()
    # Chain-aware format: A25-B200
    m = re.fullmatch(r"[A-Za-z](-?\d+)-[A-Za-z](-?\d+)", stripped)
    if m:
        return int(m.group(1)), int(m.group(2))
    # Legacy format: 25-200
    m = re.fullmatch(r"\s*(-?\d+)\s*-\s*(-?\d+)\s*", stripped)
    if m is None:
        raise ValueError(f"Invalid residue range: {residue_range!r}")

    start_residue = int(m.group(1))
    end_residue = int(m.group(2))

    if start_residue > end_residue:
        raise ValueError(f"Invalid residue range: start > end in {residue_range!r}")

    return start_residue, end_residue


def _infer_full_residue_range_from_protein_pdb(protein_pdb_path: Path) -> str:
    """Infer the full observed residue range from a representative protein PDB.

    This helper is used when the input CSV does not provide an explicit residue
    range. The inferred range is based only on observed protein-like coordinate
    residues in the first model. It does not infer missing terminal residues from
    SEQRES or UniProt because the filler module should not silently invent a
    biological range that is not present in the current template.
    """

    parsed = _parse_first_model(protein_pdb_path)
    residue_numbers: list[int] = []

    for chain in parsed.model:
        for residue in chain:
            if not _is_protein_atom_residue(residue):
                continue
            if "CA" not in residue:
                continue

            residue_numbers.append(int(residue.id[1]))

    if not residue_numbers:
        raise ValueError(
            f"No protein CA residues found for range inference in {protein_pdb_path}"
        )

    return f"{min(residue_numbers)}-{max(residue_numbers)}"


def _resolve_residue_range_for_filler(
    *,
    residue_range: str,
    template_pdb_path: Path,
) -> str:
    """Return caller-provided range or infer the full observed template range."""

    cleaned_range = str(residue_range).strip()
    if cleaned_range:
        _parse_residue_range(cleaned_range)
        return cleaned_range

    inferred_range = _infer_full_residue_range_from_protein_pdb(template_pdb_path)
    _debug(
        "No residue range provided; inferred full observed template range: "
        f"{inferred_range}"
    )
    return inferred_range


def _derive_per_chain_effective_range(
    *,
    residue_range: str,
    chain_id: str,
    template_pdb_path: Path,
) -> str:
    """Convert a possibly chain-aware range to a plain 'start-end' for one chain.

    For a cross-chain range like 'A25-B15':
    - start chain (A): lower bound is 25, upper bound is inferred from template
    - end   chain (B): lower bound is inferred from template, upper bound is 15
    - other chains: full template range (inferred)

    For a same-chain range 'A25-A200' where chain matches: return '25-200'.
    For a legacy plain range '25-200': return it unchanged.
    """
    stripped = residue_range.strip()
    if not stripped:
        return _infer_full_residue_range_from_protein_pdb(template_pdb_path)

    m = re.fullmatch(r"([A-Za-z])(-?\d+)-([A-Za-z])(-?\d+)", stripped)
    if m is None:
        # Legacy format — return as-is
        return stripped

    start_chain = m.group(1).upper()
    start_resnum = int(m.group(2))
    end_chain = m.group(3).upper()
    end_resnum = int(m.group(4))
    norm_chain = chain_id.strip().upper()

    template_range = _infer_full_residue_range_from_protein_pdb(template_pdb_path)
    template_start, template_end = _parse_residue_range(template_range)

    if start_chain == end_chain:
        if norm_chain == start_chain:
            return f"{start_resnum}-{end_resnum}"
        return f"{template_start}-{template_end}"

    # Cross-chain
    if norm_chain == start_chain:
        return f"{start_resnum}-{template_end}"
    if norm_chain == end_chain:
        return f"{template_start}-{end_resnum}"
    return f"{template_start}-{template_end}"


def _trim_final_model_in_place(model_path: Path, effective_residue_range: str) -> None:
    """Trim ATOM records in a final model to the effective residue range, in-place.

    effective_residue_range must be a plain 'start-end' integer range (chain-aware
    ranges must already have been converted to per-chain bounds before calling this).
    No-op if effective_residue_range is empty.
    """
    if not effective_residue_range.strip():
        return

    start_resnum, end_resnum = _parse_residue_range(effective_residue_range)

    kept_lines: list[str] = []
    with model_path.open("r", encoding="utf-8") as fh:
        for line in fh:
            if line.startswith("ATOM"):
                raw_resseq = line[22:26].strip()
                try:
                    resnum = int(raw_resseq)
                except ValueError:
                    kept_lines.append(line)
                    continue
                if start_resnum <= resnum <= end_resnum:
                    kept_lines.append(line)
            else:
                kept_lines.append(line)

    model_path.write_text("".join(kept_lines), encoding="utf-8")


def _find_gap_regions_in_sequence(aligned_sequence: str) -> tuple[GapRegion, ...]:
    """Detect gap regions from ``-`` characters only.

    This is intentionally consistent with the upstream alignment semantics:
    ``-`` represents an alignment gap or missing region, while ``X`` represents
    unknown residue identity. Treating ``X`` as a gap would incorrectly turn
    unresolved or unusual residues into reconstruction targets. Terminal gap
    regions are reported but filtered out before filler decisions are made.
    """

    gap_regions: list[GapRegion] = []
    sequence_length = len(aligned_sequence)

    in_gap = False
    gap_start = -1

    for index, character in enumerate(aligned_sequence):
        if character == "-" and not in_gap:
            in_gap = True
            gap_start = index
            continue

        if character != "-" and in_gap:
            gap_end = index - 1
            gap_length = gap_end - gap_start + 1
            is_terminal = gap_start == 0 or gap_end == sequence_length - 1
            gap_regions.append(
                GapRegion(
                    alignment_start=gap_start,
                    alignment_end=gap_end,
                    gap_length=gap_length,
                    is_terminal=is_terminal,
                    classification=_classify_gap_length(gap_length),
                )
            )
            in_gap = False

    if in_gap:
        gap_end = sequence_length - 1
        gap_length = gap_end - gap_start + 1
        is_terminal = gap_start == 0 or gap_end == sequence_length - 1
        gap_regions.append(
            GapRegion(
                alignment_start=gap_start,
                alignment_end=gap_end,
                gap_length=gap_length,
                is_terminal=is_terminal,
                classification=_classify_gap_length(gap_length),
            )
        )

    return tuple(gap_regions)


def _extract_sequence_block_from_ali(
    alignment_file: Path,
    record_id: str,
) -> str:
    text = alignment_file.read_text(encoding="utf-8")
    header = f">P1;{record_id}"

    if header not in text:
        raise ValueError(
            f"Record ID '{record_id}' not found in alignment file {alignment_file}"
        )

    block = text.split(header, maxsplit=1)[1]
    lines = block.strip().splitlines()

    if len(lines) < 2:
        raise ValueError(
            f"Could not parse alignment block for '{record_id}' in {alignment_file}"
        )

    sequence_lines: list[str] = []

    for line in lines[1:]:
        stripped = line.strip()

        if stripped.startswith(">P1;"):
            break

        if "*" in stripped:
            sequence_lines.append(stripped.split("*", maxsplit=1)[0])
            break

        if stripped:
            sequence_lines.append(stripped)

    sequence = "".join(sequence_lines).replace(" ", "").replace("\n", "")
    _debug(f"Extracted aligned sequence length for {record_id}: {len(sequence)}")
    return sequence


def analyze_fill_decision_from_template_alignment(
    template_alignment_skeleton: str,
) -> FillDecision:
    all_gap_regions = _find_gap_regions_in_sequence(template_alignment_skeleton)
    internal_gap_regions = tuple(
        gap_region for gap_region in all_gap_regions if not gap_region.is_terminal
    )

    _debug(f"Total gap regions found: {len(all_gap_regions)}")
    for gap_region in all_gap_regions:
        _debug(
            "Gap region: "
            f"start={gap_region.alignment_start}, "
            f"end={gap_region.alignment_end}, "
            f"length={gap_region.gap_length}, "
            f"terminal={gap_region.is_terminal}, "
            f"classification={gap_region.classification}"
        )

    _debug(f"Internal gap regions: {len(internal_gap_regions)}")

    if not internal_gap_regions:
        return FillDecision(
            should_run_modeller=False,
            overall_classification="no_internal_gaps",
            gap_regions=(),
            skip_reason="Template alignment contains no internal gaps.",
            alphafold_candidate=False,
        )

    has_alphafold_gap = any(
        gap_region.gap_length >= 8 for gap_region in internal_gap_regions
    )
    if has_alphafold_gap:
        return FillDecision(
            should_run_modeller=False,
            overall_classification="alphafold_candidate",
            gap_regions=internal_gap_regions,
            skip_reason=(
                "At least one internal gap is 8 residues or longer. "
                "Using AlphaFold fallback if possible."
            ),
            alphafold_candidate=True,
        )

    has_yellow_gap = any(
        6 <= gap_region.gap_length < 8 for gap_region in internal_gap_regions
    )
    if has_yellow_gap:
        return FillDecision(
            should_run_modeller=True,
            overall_classification="yellow",
            gap_regions=internal_gap_regions,
            skip_reason=None,
            alphafold_candidate=False,
        )

    return FillDecision(
        should_run_modeller=True,
        overall_classification="green",
        gap_regions=internal_gap_regions,
        skip_reason=None,
        alphafold_candidate=False,
    )


def analyze_fill_decision(
    alignment_file: Path,
    template_id: str,
) -> FillDecision:
    template_sequence = _extract_sequence_block_from_ali(
        alignment_file=alignment_file,
        record_id=template_id,
    )
    return analyze_fill_decision_from_template_alignment(template_sequence)


def list_available_atom_alignment_chain_ids(
    alignment_directory: Path,
) -> tuple[str, ...]:
    chain_ids: set[str] = set()

    exact_pattern = re.compile(r"^ATOM_chain_([^_/\\]+)_vs_UniProt\.aln\.fasta$")
    alt_pattern = re.compile(r"^ATOM_vs_UniProt_chain_([^_/\\]+)\.aln\.fasta$")

    for path in alignment_directory.glob("*.aln.fasta"):
        exact_match = exact_pattern.match(path.name)
        if exact_match is not None:
            chain_ids.add(_normalize_chain_id(exact_match.group(1)))
            continue

        alt_match = alt_pattern.match(path.name)
        if alt_match is not None:
            chain_ids.add(_normalize_chain_id(alt_match.group(1)))

    return tuple(sorted(chain_ids))


def _alignment_metrics_for_chain(
    alignment_directory: Path,
    chain_id: str,
) -> dict[str, int]:
    """Compute alignment-based metrics for chain selection.

    Metrics are based on ``-`` only, consistent with the current alignment
    semantics. They are used as tie-breakers after requested-range overlap. A
    chain with more actual internal gap evidence is preferred because it is more
    likely to be the chain requiring reconstruction.
    """

    alignment_fasta_path = find_alignment_fasta_for_filler(
        alignment_directory=alignment_directory,
        chain_id=chain_id,
    )
    (_, template_alignment_skeleton), _ = read_two_sequence_fasta(alignment_fasta_path)

    all_gap_regions = _find_gap_regions_in_sequence(template_alignment_skeleton)
    internal_gap_regions = [g for g in all_gap_regions if not g.is_terminal]

    return {
        "n_internal_gaps": len(internal_gap_regions),
        "max_internal_gap": max(
            (g.gap_length for g in internal_gap_regions), default=0
        ),
        "total_internal_gap": sum(g.gap_length for g in internal_gap_regions),
    }


def _collect_protein_residue_numbers_by_chain(
    template_pdb_path: Path,
    available_chain_ids: set[str],
) -> dict[str, set[int]]:
    parsed = _parse_first_model(template_pdb_path)
    residue_numbers_by_chain: dict[str, set[int]] = {}

    for chain in parsed.model:
        chain_id = _normalize_chain_id(str(chain.id))
        if chain_id not in available_chain_ids:
            continue

        for residue in chain:
            if not _is_protein_atom_residue(residue):
                continue

            residue_number = int(residue.id[1])
            residue_numbers_by_chain.setdefault(chain_id, set()).add(residue_number)

    return residue_numbers_by_chain


def choose_filler_chain_id_from_range_and_alignments(
    template_pdb_path: Path,
    alignment_directory: Path,
    residue_range: str,
) -> str:
    """Choose the modelling chain using BioPython PDB inspection.

    The selector combines residue-number overlap with the requested or inferred
    range, alignment evidence for actual internal ``-`` gaps, and total chain
    size as a deterministic tie-breaker. If ``residue_range`` is empty, the full
    observed residue range is inferred from the current template PDB. PDB parsing
    is delegated to BioPython rather than manual column slicing.
    """

    if not template_pdb_path.exists():
        raise FileNotFoundError(f"Template PDB not found: {template_pdb_path}")

    if not alignment_directory.exists():
        raise FileNotFoundError(f"Alignment directory not found: {alignment_directory}")

    available_chain_ids = set(
        list_available_atom_alignment_chain_ids(alignment_directory)
    )
    if not available_chain_ids:
        raise ValueError(
            f"No chain-specific ATOM alignment files found in {alignment_directory}"
        )

    effective_residue_range = _resolve_residue_range_for_filler(
        residue_range=residue_range,
        template_pdb_path=template_pdb_path,
    )
    range_start, range_end = _parse_residue_range(effective_residue_range)
    requested_range = set(range(range_start, range_end + 1))

    residue_numbers_by_chain = _collect_protein_residue_numbers_by_chain(
        template_pdb_path=template_pdb_path,
        available_chain_ids=available_chain_ids,
    )

    if not residue_numbers_by_chain:
        raise ValueError(
            "No protein-like residues found for chains that also have ATOM alignments. "
            f"template_pdb={template_pdb_path}, alignment_directory={alignment_directory}"
        )

    score_rows: list[tuple[str, int, int, int, int, int]] = []
    for chain_id, residue_number_set in sorted(residue_numbers_by_chain.items()):
        overlap_count = len(residue_number_set & requested_range)
        total_chain_residue_count = len(residue_number_set)

        metrics = _alignment_metrics_for_chain(alignment_directory, chain_id)

        score_rows.append(
            (
                chain_id,
                overlap_count,
                metrics["n_internal_gaps"],
                metrics["max_internal_gap"],
                metrics["total_internal_gap"],
                total_chain_residue_count,
            )
        )

        _debug(
            "Chain scoring for filler selection: "
            f"chain={chain_id}, "
            f"range={effective_residue_range}, "
            f"range_overlap={overlap_count}, "
            f"n_internal_gaps={metrics['n_internal_gaps']}, "
            f"max_internal_gap={metrics['max_internal_gap']}, "
            f"total_internal_gap={metrics['total_internal_gap']}, "
            f"total_residues={total_chain_residue_count}"
        )

    best_overlap = max(row[1] for row in score_rows)
    if best_overlap == 0:
        raise ValueError(
            "Could not determine filler chain from residue range among aligned chains: "
            f"requested range={effective_residue_range!r}"
        )

    best_rows = [row for row in score_rows if row[1] == best_overlap]

    best_internal_gap_count = max(row[2] for row in best_rows)
    best_rows = [row for row in best_rows if row[2] == best_internal_gap_count]

    best_max_gap = max(row[3] for row in best_rows)
    best_rows = [row for row in best_rows if row[3] == best_max_gap]

    best_total_gap = max(row[4] for row in best_rows)
    best_rows = [row for row in best_rows if row[4] == best_total_gap]

    best_total_residues = max(row[5] for row in best_rows)
    best_rows = [row for row in best_rows if row[5] == best_total_residues]

    if len(best_rows) == 1:
        selected_chain = best_rows[0][0]
        _debug(
            f"Selected filler chain from residue range {effective_residue_range}: "
            f"{selected_chain}"
        )
        return selected_chain

    tied_chain_ids = ", ".join(row[0] for row in best_rows)
    raise ValueError(
        "Could not determine a unique filler chain from residue range among aligned "
        f"chains. tied chains={tied_chain_ids}, range={effective_residue_range!r}"
    )


def resolve_filler_chain_id(
    *,
    chain_id: str | None,
    alignment_directory: Path,
    template_pdb_path: Path,
    residue_range: str,
) -> str:
    if chain_id is not None and str(chain_id).strip():
        normalized_chain_id = _normalize_chain_id(chain_id)
        _debug(f"Using caller-provided filler chain ID: {normalized_chain_id}")
        return normalized_chain_id

    resolved_chain_id = choose_filler_chain_id_from_range_and_alignments(
        template_pdb_path=template_pdb_path,
        alignment_directory=alignment_directory,
        residue_range=residue_range,
    )
    _debug(f"Resolved filler chain ID from residue range/alignment: {resolved_chain_id}")
    return resolved_chain_id


def find_alignment_fasta_for_filler(
    alignment_directory: Path,
    chain_id: str,
) -> Path:
    if not alignment_directory.exists():
        raise FileNotFoundError(f"Alignment directory not found: {alignment_directory}")

    normalized_chain_id = _normalize_chain_id(chain_id)

    exact_candidates = [
        alignment_directory / f"ATOM_chain_{normalized_chain_id}_vs_UniProt.aln.fasta",
        alignment_directory / f"ATOM_vs_UniProt_chain_{normalized_chain_id}.aln.fasta",
        alignment_directory / f"chain_{normalized_chain_id}_ATOM_vs_UniProt.aln.fasta",
    ]
    for candidate in exact_candidates:
        if candidate.exists():
            _debug(f"Selected chain-specific ATOM alignment FASTA: {candidate}")
            return candidate

    candidate_list = tuple(
        sorted(
            alignment_directory.glob(
                f"*chain*{normalized_chain_id}*ATOM*vs*UniProt*.aln.fasta"
            )
        )
    )
    if len(candidate_list) == 1:
        _debug(
            f"Selected fuzzy chain-specific ATOM alignment FASTA: {candidate_list[0]}"
        )
        return candidate_list[0]

    candidate_list = tuple(
        sorted(alignment_directory.glob("*ATOM*vs*UniProt*.aln.fasta"))
    )
    if len(candidate_list) == 1:
        _debug(f"Selected fallback ATOM alignment FASTA: {candidate_list[0]}")
        return candidate_list[0]

    if len(candidate_list) > 1:
        raise ValueError(
            f"Multiple ATOM-vs-UniProt alignment FASTA files found in "
            f"{alignment_directory}, but no exact chain-specific file for "
            f"chain {normalized_chain_id!r}."
        )

    raise FileNotFoundError(
        f"No ATOM-vs-UniProt alignment FASTA found for chain "
        f"{normalized_chain_id!r} in {alignment_directory}"
    )


def read_two_sequence_fasta(
    fasta_path: Path,
) -> tuple[tuple[str, str], tuple[str, str]]:
    text = fasta_path.read_text(encoding="utf-8").strip()

    if not text:
        raise ValueError(f"FASTA file is empty: {fasta_path}")

    record_list: list[tuple[str, str]] = []
    current_header: str | None = None
    current_sequence_lines: list[str] = []

    for line in text.splitlines():
        stripped = line.strip()
        if not stripped:
            continue

        if stripped.startswith(">"):
            if current_header is not None:
                record_list.append((current_header, "".join(current_sequence_lines)))
            current_header = stripped[1:].strip()
            current_sequence_lines = []
        else:
            current_sequence_lines.append(stripped)

    if current_header is not None:
        record_list.append((current_header, "".join(current_sequence_lines)))

    if len(record_list) != 2:
        raise ValueError(
            f"Expected exactly 2 FASTA records in {fasta_path}, "
            f"but found {len(record_list)}"
        )

    _debug(
        f"Read FASTA alignment with 2 records from {fasta_path}: "
        f"{record_list[0][0]!r}, {record_list[1][0]!r}"
    )
    return record_list[0], record_list[1]


def _is_pdb_atom_header(header: str) -> bool:
    header_upper = header.upper()
    return header_upper.startswith("PDB|") or "ATOM" in header_upper


def _is_uniprot_header(header: str) -> bool:
    header_upper = header.upper()
    return (
        header_upper.startswith("SP|")
        or header_upper.startswith("TR|")
        or "UNIPROT" in header_upper
    )


def _extract_chain_id_from_alignment_header(header: str) -> str | None:
    patterns = [
        r"(?:^|\|)chain_([A-Za-z0-9_])(?:\||$)",
        r"(?:^|\|)chain-([A-Za-z0-9_])(?:\||$)",
        r"(?:^|\|)chain:([A-Za-z0-9_])(?:\||$)",
        r"(?:^|\|)chain ([A-Za-z0-9_])(?:\||$)",
    ]

    for pattern in patterns:
        match = re.search(pattern, header, flags=re.IGNORECASE)
        if match is not None:
            return _normalize_chain_id(match.group(1))

    return None


def split_template_and_target_alignment_records(
    alignment_fasta_path: Path,
) -> tuple[str, str, str, str]:
    (header_1, seq_1), (header_2, seq_2) = read_two_sequence_fasta(alignment_fasta_path)

    if _is_pdb_atom_header(header_1) and _is_uniprot_header(header_2):
        return header_1, seq_1, header_2, seq_2

    if _is_pdb_atom_header(header_2) and _is_uniprot_header(header_1):
        return header_2, seq_2, header_1, seq_1

    raise ValueError(
        "Could not determine template/target rows from FASTA headers: "
        f"{header_1!r}, {header_2!r}"
    )


def _chain_has_protein_atoms(model: Any, chain_id: str) -> bool:
    normalized_chain_id = _normalize_chain_id(chain_id)

    for chain in model:
        if _normalize_chain_id(str(chain.id)) != normalized_chain_id:
            continue

        for residue in chain:
            if _is_protein_atom_residue(residue):
                return True

    return False


def write_chain_specific_template_pdb(
    template_pdb_path: Path,
    output_dir: Path,
    template_id: str,
    chain_id: str,
) -> Path:
    """Write a chain-specific MODELLER template with BioPython.

    The function keeps only protein-like ATOM residues from the requested chain
    and writes them through ``PDBIO``. It no longer copies records by slicing PDB
    text columns. The selected template is written as ``<template_id>.pdb`` in
    the filler output directory so MODELLER can find it next to the generated
    alignment file and script.
    """

    if not template_pdb_path.exists():
        raise FileNotFoundError(f"Template PDB not found: {template_pdb_path}")

    output_dir.mkdir(parents=True, exist_ok=True)

    normalized_chain_id = _normalize_chain_id(chain_id)
    destination = output_dir / f"{template_id}.pdb"

    parsed = _parse_first_model(template_pdb_path)
    if not _chain_has_protein_atoms(parsed.model, normalized_chain_id):
        raise ValueError(
            f"No protein-like ATOM residues found for chain {normalized_chain_id!r} "
            f"in {template_pdb_path}"
        )

    io = PDBIO()
    io.set_structure(parsed.structure)
    io.save(
        str(destination),
        select=_FirstModelChainProteinSelect(
            model_id=parsed.model.id,
            chain_id=normalized_chain_id,
        ),
        preserve_atom_numbering=True,
    )

    _validate_written_pdb_has_atoms(
        output_pdb_path=destination,
        description=f"Chain-specific template PDB for chain {normalized_chain_id}",
    )

    _debug(
        f"Wrote chain-specific template PDB for chain {normalized_chain_id} "
        f"({template_id}) to: {destination}"
    )
    return destination


def extract_sequence_from_template_pdb(template_pdb_path: Path) -> str:
    """Extract a filtered peptide sequence from a template PDB with BioPython.

    Residues are read from the first model and filtered to protein-like ATOM
    residues with the backbone atoms ``N``, ``CA``, and ``C`` present. The residue
    conversion uses BioPython's extended PDB residue table plus the small
    force-field alias map in this module. Residues that are protein-like but
    cannot be converted are represented as ``X``.
    """

    if not template_pdb_path.exists():
        raise FileNotFoundError(f"Template PDB not found: {template_pdb_path}")

    parsed = _parse_first_model(template_pdb_path)

    residue_sequence_list: list[str] = []
    skipped_residue_debug_list: list[str] = []

    for chain in parsed.model:
        chain_id = _normalize_chain_id(str(chain.id))

        for residue in chain:
            if not _is_protein_atom_residue(residue):
                continue

            residue_number = residue.id[1]
            insertion_code = str(residue.id[2]).strip()
            residue_name = residue.get_resname().strip().upper()

            if not _residue_has_peptide_backbone(residue):
                skipped_residue_debug_list.append(
                    f"{chain_id}:{residue_number}{insertion_code}:{residue_name}"
                )
                continue

            residue_sequence_list.append(_residue_to_one_letter(residue))

    if not residue_sequence_list:
        raise ValueError(
            f"No peptide-like ATOM residues found in template PDB: {template_pdb_path}"
        )

    sequence = "".join(residue_sequence_list)

    _debug(
        f"Extracted filtered template PDB sequence length from {template_pdb_path}: "
        f"{len(sequence)}"
    )
    _debug(f"Template sequence tail: {sequence[-80:]}")

    if skipped_residue_debug_list:
        preview = ", ".join(skipped_residue_debug_list[:15])
        if len(skipped_residue_debug_list) > 15:
            preview += ", ..."
        _debug(
            "Skipped non-peptide ATOM residues while extracting template sequence: "
            f"{preview}"
        )

    return sequence


def _first_sequence_mismatch(
    sequence_a: str,
    sequence_b: str,
) -> tuple[int, str, str] | None:
    shortest_length = min(len(sequence_a), len(sequence_b))
    for index in range(shortest_length):
        if sequence_a[index] != sequence_b[index]:
            return index, sequence_a[index], sequence_b[index]

    if len(sequence_a) != len(sequence_b):
        a_char = (
            sequence_a[shortest_length]
            if len(sequence_a) > shortest_length
            else "<END>"
        )
        b_char = (
            sequence_b[shortest_length]
            if len(sequence_b) > shortest_length
            else "<END>"
        )
        return shortest_length, a_char, b_char

    return None


def _write_sequence_debug_files(
    output_dir: Path,
    template_alignment_skeleton: str,
    actual_template_sequence: str,
) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    (output_dir / "debug_template_alignment_skeleton.txt").write_text(
        template_alignment_skeleton + "\n",
        encoding="utf-8",
    )
    (output_dir / "debug_template_alignment_ungapped.txt").write_text(
        _strip_gaps(template_alignment_skeleton) + "\n",
        encoding="utf-8",
    )
    (output_dir / "debug_template_pdb_sequence.txt").write_text(
        actual_template_sequence + "\n",
        encoding="utf-8",
    )


def validate_template_sequence_consistency(
    template_alignment_skeleton: str,
    actual_template_sequence: str,
    output_dir: Path,
    template_header: str,
    template_pdb_path: Path,
    *,
    max_allowed_length_difference: int = 25,
) -> None:
    normalized_alignment_template = _normalize_alignment_template_for_validation(
        template_alignment_skeleton
    )

    if normalized_alignment_template == actual_template_sequence:
        return

    _write_sequence_debug_files(
        output_dir=output_dir,
        template_alignment_skeleton=template_alignment_skeleton,
        actual_template_sequence=actual_template_sequence,
    )

    mismatch = _first_sequence_mismatch(
        normalized_alignment_template,
        actual_template_sequence,
    )

    mismatch_text = "unknown"
    if mismatch is not None:
        mismatch_index, alignment_char, actual_char = mismatch
        mismatch_text = (
            f"first mismatch at ungapped index {mismatch_index}: "
            f"alignment={alignment_char!r}, pdb_sequence={actual_char!r}"
        )

    alignment_length = len(normalized_alignment_template)
    actual_length = len(actual_template_sequence)
    length_difference = abs(alignment_length - actual_length)

    message = (
        "Template sequence mismatch between alignment FASTA and chain-specific "
        "template PDB. "
        f"Alignment header={template_header!r}. "
        f"Template PDB={template_pdb_path}. "
        f"Normalized alignment ungapped length={alignment_length}. "
        f"PDB-derived length={actual_length}. "
        f"{mismatch_text}. "
        f"Debug files written in {output_dir}."
    )

    if length_difference > max_allowed_length_difference:
        raise ValueError(message)

    _debug(f"WARNING: {message}")


def build_modeller_template_alignment_sequence(
    aligned_template_skeleton: str,
    actual_template_sequence: str,
) -> str:
    """Rebuild the MODELLER template row from the real template sequence.

    Legacy alignments can contain trailing placeholder ``X`` characters after
    the true template sequence is exhausted. Those are converted to gaps so they
    do not become fake residues in the MODELLER alignment. Internal non-gap
    characters still consume real template residues, preserving the existing
    alignment skeleton.
    """

    rebuilt_character_list: list[str] = []
    template_index = 0
    actual_template_length = len(actual_template_sequence)

    for skeleton_index, character in enumerate(aligned_template_skeleton):
        if character == "-":
            rebuilt_character_list.append("-")
            continue

        if template_index < actual_template_length:
            rebuilt_character_list.append(actual_template_sequence[template_index])
            template_index += 1
            continue

        if character == "X":
            rebuilt_character_list.append("-")
            continue

        raise ValueError(
            "Alignment template skeleton expects more real template residues than "
            f"available in the template PDB sequence. "
            f"First invalid extra non-gap character: {character!r} "
            f"at skeleton index {skeleton_index}. "
            f"Actual template length: {actual_template_length}."
        )

    if template_index != actual_template_length:
        raise ValueError(
            "Alignment template skeleton consumed fewer residues than present in "
            f"the real template PDB sequence: consumed={template_index}, "
            f"actual={actual_template_length}"
        )

    rebuilt_sequence = "".join(rebuilt_character_list)

    _debug(
        f"Rebuilt template aligned sequence length: {len(rebuilt_sequence)}; "
        f"ungapped length: {len(rebuilt_sequence.replace('-', ''))}"
    )
    return rebuilt_sequence


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

    _debug(f"Template FASTA header selected: {template_header}")
    _debug(f"Target FASTA header selected: {target_header}")

    header_chain_id = _extract_chain_id_from_alignment_header(template_header)
    normalized_chain_id = _normalize_chain_id(chain_id)
    if header_chain_id is not None and header_chain_id != normalized_chain_id:
        raise ValueError(
            f"Alignment FASTA template header chain mismatch: "
            f"expected chain {normalized_chain_id!r}, "
            f"but header {template_header!r} suggests chain {header_chain_id!r}."
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


def write_modeller_script(
    output_dir: Path,
    alignment_file: Path,
    template_id: str,
    target_id: str,
    starting_model: int = 1,
    ending_model: int = 1,
) -> Path:
    script_path = output_dir / DEFAULT_MODELLER_SCRIPT_FILENAME

    content = f"""from modeller import *
from modeller.automodel import *

log.verbose()

env = environ()
env.io.atom_files_directory = ['.']
env.io.hetatm = True

a = automodel(
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
    _debug(f"Wrote MODELLER script to: {script_path}")
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

    _debug(f"Validated MODELLER template file: {expected_template_pdb}")
    _debug(f"Validated MODELLER alignment file: {alignment_file}")


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

    _debug(f"Running command: {' '.join(command)}")
    _debug(f"Using PYTHONPATH: {env['PYTHONPATH']}")

    result = subprocess.run(
        command,
        cwd=working_dir,
        capture_output=True,
        text=True,
        env=env,
    )

    stdout_log = working_dir / DEFAULT_STDOUT_LOG_FILENAME
    stderr_log = working_dir / DEFAULT_STDERR_LOG_FILENAME

    stdout_log.write_text(result.stdout, encoding="utf-8")
    stderr_log.write_text(result.stderr, encoding="utf-8")

    _debug(f"MODELLER return code: {result.returncode}")

    if result.returncode != 0:
        raise RuntimeError(
            f"MODELLER failed for script {script_path.name}. "
            f"See {stdout_log} and {stderr_log}."
        )

    return stdout_log, stderr_log


def _write_skip_logs(
    output_dir: Path,
    message: str,
) -> tuple[Path, Path]:
    stdout_log = output_dir / DEFAULT_STDOUT_LOG_FILENAME
    stderr_log = output_dir / DEFAULT_STDERR_LOG_FILENAME
    stdout_log.write_text(message + "\n", encoding="utf-8")
    stderr_log.write_text("", encoding="utf-8")
    return stdout_log, stderr_log


def find_raw_models(
    output_dir: Path,
    target_id: str,
) -> tuple[Path, ...]:
    preferred_models = tuple(sorted(output_dir.glob(f"{target_id}.B*.pdb")))
    if preferred_models:
        _debug(
            f"Found {len(preferred_models)} raw MODELLER model(s) for target {target_id}"
        )
        return preferred_models

    fallback_models = tuple(
        sorted(
            path
            for path in output_dir.glob(f"{target_id}*.pdb")
            if path.name != DEFAULT_FINAL_MODEL_FILENAME
        )
    )
    _debug(
        f"Found {len(fallback_models)} fallback raw model candidate(s) "
        f"for target {target_id}"
    )
    return fallback_models


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
            _debug(
                f"WARNING: model_scores.tsv selected {best_model_name}, "
                f"but file does not exist in {output_dir}"
            )
        else:
            _debug(f"WARNING: No valid score rows found in {score_file}")

    if len(raw_model_paths) == 1:
        _debug(f"Falling back to sole raw model: {raw_model_paths[0]}")
        return raw_model_paths[0]

    if raw_model_paths:
        _debug(
            "WARNING: Falling back to first raw model because score-based "
            f"selection failed. Selected: {raw_model_paths[0]}"
        )
        return raw_model_paths[0]

    raise FileNotFoundError(
        f"No raw MODELLER models found for target {target_id!r} in {output_dir}"
    )


def cleanup_model_pdb(
    input_model_path: Path,
    output_model_path: Path,
) -> Path:
    """Write a protein-only cleaned model with BioPython.

    The cleanup keeps the first model and removes non-protein residues through
    a ``PDBIO`` selector. This avoids direct text filtering such as
    ``line.startswith("ATOM")`` while preserving the old behavior of returning a
    protein-only PDB for downstream preparation. The output is validated so an
    empty cleaned model fails immediately.
    """

    if not input_model_path.exists():
        raise FileNotFoundError(f"Model PDB not found: {input_model_path}")

    parsed = _parse_first_model(input_model_path)

    output_model_path.parent.mkdir(parents=True, exist_ok=True)

    io = PDBIO()
    io.set_structure(parsed.structure)
    io.save(
        str(output_model_path),
        select=_FirstModelProteinOnlySelect(model_id=parsed.model.id),
        preserve_atom_numbering=True,
    )

    _validate_written_pdb_has_atoms(
        output_pdb_path=output_model_path,
        description="Cleaned final model",
    )

    _debug(f"Cleaned final model written to: {output_model_path}")
    return output_model_path


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

    _debug(f"Best raw model path: {best_model}")
    _debug(f"Final cleaned model path: {cleaned_final_model}")
    return best_model, cleaned_final_model


def download_alphafold_structure(
    uniprot_id: str,
    output_dir: Path,
) -> Path | None:
    output_dir.mkdir(parents=True, exist_ok=True)

    api_url = f"https://alphafold.ebi.ac.uk/api/prediction/{uniprot_id}"

    try:
        with urlopen(api_url) as response:
            payload = json.load(response)
    except HTTPError as exc:
        if exc.code == 404:
            return None
        raise
    except URLError:
        raise

    if not payload:
        return None

    record = payload[0]
    pdb_url = record.get("pdbUrl")

    if not pdb_url:
        return None

    target_path = output_dir / DEFAULT_ALPHAFOLD_RAW_FILENAME
    urlretrieve(pdb_url, target_path)

    if not target_path.is_file() or target_path.stat().st_size == 0:
        return None

    return target_path


def _collect_protein_ca_atoms_by_resseq(structure: Any) -> dict[int, object]:
    atom_by_resseq: dict[int, object] = {}

    for model in structure:
        for chain in model:
            for residue in chain:
                if not _is_protein_atom_residue(residue):
                    continue
                if "CA" not in residue:
                    continue

                residue_number = int(residue.id[1])
                atom_by_resseq.setdefault(residue_number, residue["CA"])

    return atom_by_resseq


def _collect_protein_ca_atoms_in_order(structure: Any) -> list[object]:
    atom_list: list[object] = []

    for model in structure:
        for chain in model:
            for residue in chain:
                if not _is_protein_atom_residue(residue):
                    continue
                if "CA" in residue:
                    atom_list.append(residue["CA"])

    return atom_list


def align_protonated_alphafold_model_to_start_pdb(
    reference_pdb_path: Path,
    mobile_pdb_path: Path,
    output_pdb_path: Path,
) -> dict[str, str | float | bool | int]:
    if not reference_pdb_path.exists():
        raise FileNotFoundError(reference_pdb_path)

    if not mobile_pdb_path.exists():
        raise FileNotFoundError(mobile_pdb_path)

    parser = PDBParser(QUIET=True)

    reference_structure = parser.get_structure("ref", str(reference_pdb_path))
    mobile_structure = parser.get_structure("mob", str(mobile_pdb_path))

    reference_by_resseq = _collect_protein_ca_atoms_by_resseq(reference_structure)
    mobile_by_resseq = _collect_protein_ca_atoms_by_resseq(mobile_structure)

    common_resseqs = sorted(set(reference_by_resseq) & set(mobile_by_resseq))

    fixed_atom_list = [reference_by_resseq[resseq] for resseq in common_resseqs]
    moving_atom_list = [mobile_by_resseq[resseq] for resseq in common_resseqs]

    pairing_mode = "residue_number"

    if len(fixed_atom_list) < 3:
        reference_ordered = _collect_protein_ca_atoms_in_order(reference_structure)
        mobile_ordered = _collect_protein_ca_atoms_in_order(mobile_structure)

        if not reference_ordered:
            raise ValueError(
                f"No protein CA atoms found in reference: {reference_pdb_path}"
            )

        if not mobile_ordered:
            raise ValueError(f"No protein CA atoms found in mobile: {mobile_pdb_path}")

        n_matched_atoms = min(len(reference_ordered), len(mobile_ordered))
        if n_matched_atoms < 3:
            raise ValueError(
                "Not enough CA atoms for reliable alignment: "
                f"reference={len(reference_ordered)}, mobile={len(mobile_ordered)}"
            )

        fixed_atom_list = reference_ordered[:n_matched_atoms]
        moving_atom_list = mobile_ordered[:n_matched_atoms]
        pairing_mode = "positional_fallback"

    superimposer = Superimposer()
    superimposer.set_atoms(fixed_atom_list, moving_atom_list)
    superimposer.apply(mobile_structure.get_atoms())

    output_pdb_path.parent.mkdir(parents=True, exist_ok=True)

    io = PDBIO()
    io.set_structure(mobile_structure)
    io.save(str(output_pdb_path), preserve_atom_numbering=True)

    return {
        "alignment_success": output_pdb_path.exists()
        and output_pdb_path.stat().st_size > 0,
        "alignment_rmsd": float(superimposer.rms),
        "alignment_output_path": str(output_pdb_path),
        "alignment_pairing_mode": pairing_mode,
        "alignment_n_pairs": len(fixed_atom_list),
    }


def _count_protein_residues_in_range(
    input_pdb_path: Path,
    start_residue: int,
    end_residue: int,
) -> int:
    parsed = _parse_first_model(input_pdb_path)
    count = 0

    for chain in parsed.model:
        for residue in chain:
            if not _is_protein_atom_residue(residue):
                continue

            residue_number = int(residue.id[1])
            if start_residue <= residue_number <= end_residue:
                count += 1

    return count


def crop_pdb_to_range(
    input_pdb_path: Path,
    output_pdb_path: Path,
    residue_range: str,
) -> Path:
    """Crop a PDB to a residue-number range using BioPython selectors."""

    if input_pdb_path is None:
        raise FileNotFoundError("AlphaFold download returned no PDB path")

    if not input_pdb_path.exists():
        raise FileNotFoundError(f"Input PDB not found: {input_pdb_path}")

    start_residue, end_residue = _parse_residue_range(residue_range)

    kept_residue_count = _count_protein_residues_in_range(
        input_pdb_path=input_pdb_path,
        start_residue=start_residue,
        end_residue=end_residue,
    )
    if kept_residue_count == 0:
        raise ValueError(
            f"No protein residues remained after cropping {input_pdb_path} "
            f"to range {residue_range!r}"
        )

    parsed = _parse_first_model(input_pdb_path)
    output_pdb_path.parent.mkdir(parents=True, exist_ok=True)

    io = PDBIO()
    io.set_structure(parsed.structure)
    io.save(
        str(output_pdb_path),
        select=_FirstModelProteinRangeSelect(
            model_id=parsed.model.id,
            start_residue=start_residue,
            end_residue=end_residue,
        ),
        preserve_atom_numbering=True,
    )

    _validate_written_pdb_has_atoms(
        output_pdb_path=output_pdb_path,
        description=f"Cropped AlphaFold PDB for range {residue_range!r}",
    )

    _debug(f"Cropped AlphaFold PDB to range {residue_range}: {output_pdb_path}")
    return output_pdb_path


def run_alphafold_fallback_for_chain(
    output_dir: Path,
    template_pdb_path: Path,
    uniprot_id: str,
    residue_range: str,
    final_model_name: str,
    model_version: int = 4,
) -> Path:
    _ = model_version

    effective_residue_range = _resolve_residue_range_for_filler(
        residue_range=residue_range,
        template_pdb_path=template_pdb_path,
    )

    alphafold_dir = output_dir / DEFAULT_ALPHAFOLD_DIRNAME

    downloaded_pdb = download_alphafold_structure(
        uniprot_id=uniprot_id,
        output_dir=alphafold_dir,
    )

    if downloaded_pdb is None:
        raise FileNotFoundError(
            f"No AlphaFold PDB available for UniProt {uniprot_id!r}"
        )

    cropped_model_path = crop_pdb_to_range(
        input_pdb_path=downloaded_pdb,
        output_pdb_path=alphafold_dir / DEFAULT_ALPHAFOLD_CROPPED_FILENAME,
        residue_range=effective_residue_range,
    )

    aligned_model_path = alphafold_dir / DEFAULT_ALPHAFOLD_ALIGNED_FILENAME
    alignment_result = align_protonated_alphafold_model_to_start_pdb(
        reference_pdb_path=template_pdb_path,
        mobile_pdb_path=cropped_model_path,
        output_pdb_path=aligned_model_path,
    )

    _debug(
        "AlphaFold alignment result: "
        f"success={alignment_result['alignment_success']}, "
        f"rmsd={alignment_result['alignment_rmsd']}, "
        f"pairs={alignment_result['alignment_n_pairs']}, "
        f"pairing_mode={alignment_result['alignment_pairing_mode']}, "
        f"output={alignment_result['alignment_output_path']}"
    )

    final_model_path = output_dir / final_model_name
    cleaned_model_path = cleanup_model_pdb(
        input_model_path=aligned_model_path,
        output_model_path=final_model_path,
    )

    return cleaned_model_path


def run_filler_for_chain(
    alignment_directory: Path,
    template_pdb_path: Path,
    output_dir: Path,
    template_id: str,
    target_id: str,
    chain_id: str | None,
    final_model_name: str | None = None,
    starting_model: int = 1,
    ending_model: int = 1,
    skip_if_no_internal_gaps: bool = True,
    uniprot_id: str | None = None,
    residue_range: str = "",
) -> FillerRunResult:
    pdb_id = _infer_pdb_id_from_target_id(
        target_id=target_id, template_pdb_path=template_pdb_path
    )
    module_log_path = _build_module_log_path(
        pdb_id=pdb_id,
        template_pdb_path=template_pdb_path,
    )

    write_log_header(
        log_path=module_log_path,
        module_name=MODULE_NAME,
        pdb_id=pdb_id,
        extra_lines=[
            f"template_id={template_id}",
            f"target_id={target_id}",
            f"input_chain_id={chain_id if chain_id is not None else '(auto)'}",
            f"alignment_directory={alignment_directory}",
            f"template_pdb_path={template_pdb_path}",
            f"output_dir={output_dir}",
            f"starting_model={starting_model}",
            f"ending_model={ending_model}",
            f"skip_if_no_internal_gaps={skip_if_no_internal_gaps}",
            f"uniprot_id={uniprot_id or '(none)'}",
            f"residue_range={residue_range or '(infer full observed range)'}",
        ],
    )

    _print_start_box(
        pdb_id=pdb_id,
        chain_id=chain_id,
        template_id=template_id,
        target_id=target_id,
        alignment_directory=alignment_directory,
        template_pdb_path=template_pdb_path,
        output_dir=output_dir,
        residue_range=residue_range,
        uniprot_id=uniprot_id,
        log_path=module_log_path,
    )

    _log_debug(module_log_path, "=== run_filler_for_chain START ===")
    _log_debug(module_log_path, f"alignment_directory: {alignment_directory}")
    _log_debug(module_log_path, f"template_pdb_path: {template_pdb_path}")
    _log_debug(module_log_path, f"output_dir: {output_dir}")
    _log_debug(module_log_path, f"template_id: {template_id}")
    _log_debug(module_log_path, f"target_id: {target_id}")
    _log_debug(module_log_path, f"chain_id (input): {chain_id}")
    _log_debug(module_log_path, f"starting_model: {starting_model}")
    _log_debug(module_log_path, f"ending_model: {ending_model}")
    _log_debug(module_log_path, f"uniprot_id: {uniprot_id}")
    _log_debug(module_log_path, f"residue_range: {residue_range}")

    if not alignment_directory.exists():
        error = FileNotFoundError(
            f"Alignment directory not found: {alignment_directory}"
        )
        append_exception_traceback(module_log_path, error)
        raise error

    if not template_pdb_path.exists():
        error = FileNotFoundError(f"Template PDB not found: {template_pdb_path}")
        append_exception_traceback(module_log_path, error)
        raise error

    output_dir.mkdir(parents=True, exist_ok=True)

    if final_model_name is None:
        final_model_name = DEFAULT_FINAL_MODEL_FILENAME

    effective_residue_range = _resolve_residue_range_for_filler(
        residue_range=residue_range,
        template_pdb_path=template_pdb_path,
    )
    _log_debug(module_log_path, f"effective_residue_range: {effective_residue_range}")

    resolved_chain_id = resolve_filler_chain_id(
        chain_id=chain_id,
        alignment_directory=alignment_directory,
        template_pdb_path=template_pdb_path,
        residue_range=effective_residue_range,
    )
    _log_debug(module_log_path, f"Resolved chain_id: {resolved_chain_id}")

    # For chain-aware ranges derive the per-chain plain interval now that we
    # know which chain we are modelling.
    effective_residue_range = _derive_per_chain_effective_range(
        residue_range=effective_residue_range,
        chain_id=resolved_chain_id,
        template_pdb_path=template_pdb_path,
    )
    _log_debug(
        module_log_path,
        f"effective_residue_range (per-chain): {effective_residue_range}",
    )

    script_path = output_dir / DEFAULT_MODELLER_SCRIPT_FILENAME

    alignment_fasta_path = find_alignment_fasta_for_filler(
        alignment_directory=alignment_directory,
        chain_id=resolved_chain_id,
    )
    _log_debug(
        module_log_path, f"Alignment FASTA selected for filler: {alignment_fasta_path}"
    )

    copied_template_pdb = write_chain_specific_template_pdb(
        template_pdb_path=template_pdb_path,
        output_dir=output_dir,
        template_id=template_id,
        chain_id=resolved_chain_id,
    )

    (
        template_header,
        template_alignment_skeleton,
        target_header,
        _target_aligned_sequence,
    ) = split_template_and_target_alignment_records(alignment_fasta_path)

    _log_debug(module_log_path, f"Template FASTA header selected: {template_header}")
    _log_debug(module_log_path, f"Target FASTA header selected: {target_header}")

    header_chain_id = _extract_chain_id_from_alignment_header(template_header)
    normalized_chain_id = _normalize_chain_id(resolved_chain_id)
    if header_chain_id is not None and header_chain_id != normalized_chain_id:
        error = ValueError(
            f"Alignment FASTA template header chain mismatch: "
            f"expected chain {normalized_chain_id!r}, "
            f"but header {template_header!r} suggests chain {header_chain_id!r}."
        )
        append_exception_traceback(module_log_path, error)
        raise error

    fill_decision = analyze_fill_decision_from_template_alignment(
        template_alignment_skeleton=template_alignment_skeleton,
    )

    append_key_value_block(
        log_path=module_log_path,
        title="FILL DECISION",
        key_value_pairs=[
            ("resolved_chain_id", resolved_chain_id),
            ("effective_residue_range", effective_residue_range),
            ("should_run_modeller", str(fill_decision.should_run_modeller)),
            ("overall_classification", fill_decision.overall_classification),
            ("alphafold_candidate", str(fill_decision.alphafold_candidate)),
            ("n_gap_regions", str(len(fill_decision.gap_regions))),
            ("skip_reason", fill_decision.skip_reason or "(none)"),
        ],
    )

    for index, gap_region in enumerate(fill_decision.gap_regions, start=1):
        append_key_value_block(
            log_path=module_log_path,
            title=f"GAP REGION {index}",
            key_value_pairs=[
                ("alignment_start", str(gap_region.alignment_start)),
                ("alignment_end", str(gap_region.alignment_end)),
                ("gap_length", str(gap_region.gap_length)),
                ("is_terminal", str(gap_region.is_terminal)),
                ("classification", gap_region.classification),
            ],
        )

    _print_decision_table(
        pdb_id=pdb_id,
        resolved_chain_id=resolved_chain_id,
        fill_decision=fill_decision,
    )
    _print_gap_regions_table(fill_decision)

    _log_debug(
        module_log_path,
        "Fill decision: "
        f"should_run_modeller={fill_decision.should_run_modeller}, "
        f"overall_classification={fill_decision.overall_classification}, "
        f"alphafold_candidate={fill_decision.alphafold_candidate}",
    )

    modeller_model_path: Path | None = None
    alphafold_model_path: Path | None = None
    raw_model_paths: tuple[Path, ...] = ()
    alignment_file = output_dir / DEFAULT_ALIGNMENT_FILENAME

    try:
        if (
            skip_if_no_internal_gaps
            and fill_decision.overall_classification == "no_internal_gaps"
        ):
            stdout_log, stderr_log = _write_skip_logs(
                output_dir=output_dir,
                message=(
                    "Filler skipped: no internal gaps in template alignment. "
                    f"Reason: {fill_decision.skip_reason}"
                ),
            )

            append_log_text(module_log_path, "[RESULT] skipped: no_internal_gaps")
            _log_debug(module_log_path, "No internal gaps found -> skipping filler run")
            _log_debug(module_log_path, "=== run_filler_for_chain END (SKIPPED) ===")

            result = _build_result(
                chain_id=resolved_chain_id,
                output_dir=output_dir,
                alignment_file=alignment_file,
                template_pdb=copied_template_pdb,
                script_file=script_path,
                modeller_model_path=None,
                alphafold_model_path=None,
                final_model_path=None,
                raw_model_paths=(),
                stdout_log=stdout_log,
                stderr_log=stderr_log,
                skipped=True,
                skip_reason=fill_decision.skip_reason,
                fill_decision=fill_decision,
            )

            _print_end_box(
                pdb_id=pdb_id,
                resolved_chain_id=resolved_chain_id,
                status="skipped",
                final_model_path=None,
                modeller_model_path=None,
                alphafold_model_path=None,
                stdout_log=stdout_log,
                stderr_log=stderr_log,
                skip_reason=fill_decision.skip_reason,
            )
            return result

        if fill_decision.alphafold_candidate:
            if not uniprot_id:
                stdout_log, stderr_log = _write_skip_logs(
                    output_dir=output_dir,
                    message=(
                        "Filler skipped: AlphaFold fallback required, but UniProt ID "
                        "is missing."
                    ),
                )
                append_log_text(
                    module_log_path, "[RESULT] skipped: missing_uniprot_for_alphafold"
                )
                _log_debug(
                    module_log_path,
                    "AlphaFold fallback needed but UniProt ID is missing",
                )
                _log_debug(
                    module_log_path, "=== run_filler_for_chain END (SKIPPED) ==="
                )

                result = _build_result(
                    chain_id=resolved_chain_id,
                    output_dir=output_dir,
                    alignment_file=alignment_file,
                    template_pdb=copied_template_pdb,
                    script_file=script_path,
                    modeller_model_path=None,
                    alphafold_model_path=None,
                    final_model_path=None,
                    raw_model_paths=(),
                    stdout_log=stdout_log,
                    stderr_log=stderr_log,
                    skipped=True,
                    skip_reason="AlphaFold fallback required, but UniProt ID is missing.",
                    fill_decision=fill_decision,
                )

                _print_end_box(
                    pdb_id=pdb_id,
                    resolved_chain_id=resolved_chain_id,
                    status="skipped",
                    final_model_path=None,
                    modeller_model_path=None,
                    alphafold_model_path=None,
                    stdout_log=stdout_log,
                    stderr_log=stderr_log,
                    skip_reason="AlphaFold fallback required, but UniProt ID is missing.",
                )
                return result

            stdout_log, stderr_log = _write_skip_logs(
                output_dir=output_dir,
                message=(
                    "MODELLER skipped: AlphaFold fallback selected because at least "
                    "one internal gap is 8 residues or longer."
                ),
            )

            alphafold_model_path = run_alphafold_fallback_for_chain(
                output_dir=output_dir,
                template_pdb_path=copied_template_pdb,
                uniprot_id=uniprot_id,
                residue_range=effective_residue_range,
                final_model_name=final_model_name,
                model_version=4,
            )

            append_key_value_block(
                log_path=module_log_path,
                title="ALPHAFOLD OUTPUT",
                key_value_pairs=[
                    ("effective_residue_range", effective_residue_range),
                    ("alphafold_model_path", str(alphafold_model_path)),
                    ("stdout_log", str(stdout_log)),
                    ("stderr_log", str(stderr_log)),
                ],
            )

            _trim_final_model_in_place(alphafold_model_path, effective_residue_range)
            _log_debug(
                module_log_path,
                f"Trimmed AlphaFold model to range {effective_residue_range!r}: {alphafold_model_path}",
            )
            _log_debug(
                module_log_path, "=== run_filler_for_chain END (ALPHAFOLD FALLBACK) ==="
            )

            result = _build_result(
                chain_id=resolved_chain_id,
                output_dir=output_dir,
                alignment_file=alignment_file,
                template_pdb=copied_template_pdb,
                script_file=script_path,
                modeller_model_path=None,
                alphafold_model_path=alphafold_model_path,
                final_model_path=alphafold_model_path,
                raw_model_paths=(),
                stdout_log=stdout_log,
                stderr_log=stderr_log,
                skipped=False,
                skip_reason=None,
                fill_decision=fill_decision,
            )

            _print_end_box(
                pdb_id=pdb_id,
                resolved_chain_id=resolved_chain_id,
                status="success_alphafold",
                final_model_path=alphafold_model_path,
                modeller_model_path=None,
                alphafold_model_path=alphafold_model_path,
                stdout_log=stdout_log,
                stderr_log=stderr_log,
                skip_reason=None,
            )
            return result

        if not fill_decision.should_run_modeller:
            stdout_log, stderr_log = _write_skip_logs(
                output_dir=output_dir,
                message=(
                    "MODELLER not run for this chain. "
                    f"Reason: {fill_decision.skip_reason or 'unspecified'}"
                ),
            )

            append_log_text(module_log_path, "[RESULT] skipped: modeller_not_run")
            _log_debug(module_log_path, "MODELLER not run for this chain")
            _log_debug(module_log_path, "=== run_filler_for_chain END (SKIPPED) ===")

            result = _build_result(
                chain_id=resolved_chain_id,
                output_dir=output_dir,
                alignment_file=alignment_file,
                template_pdb=copied_template_pdb,
                script_file=script_path,
                modeller_model_path=None,
                alphafold_model_path=None,
                final_model_path=None,
                raw_model_paths=(),
                stdout_log=stdout_log,
                stderr_log=stderr_log,
                skipped=True,
                skip_reason=fill_decision.skip_reason,
                fill_decision=fill_decision,
            )

            _print_end_box(
                pdb_id=pdb_id,
                resolved_chain_id=resolved_chain_id,
                status="skipped",
                final_model_path=None,
                modeller_model_path=None,
                alphafold_model_path=None,
                stdout_log=stdout_log,
                stderr_log=stderr_log,
                skip_reason=fill_decision.skip_reason,
            )
            return result

        alignment_file = write_modeller_alignment_from_existing_alignment(
            alignment_fasta_path=alignment_fasta_path,
            template_pdb_path=copied_template_pdb,
            output_dir=output_dir,
            template_id=template_id,
            target_id=target_id,
            chain_id=resolved_chain_id,
        )

        validate_modeller_inputs(
            output_dir=output_dir,
            alignment_file=alignment_file,
            template_id=template_id,
        )

        script_path = write_modeller_script(
            output_dir=output_dir,
            alignment_file=alignment_file,
            template_id=template_id,
            target_id=target_id,
            starting_model=starting_model,
            ending_model=ending_model,
        )

        stdout_log, stderr_log = run_modeller_binary(
            script_path=script_path,
            working_dir=output_dir,
        )

        raw_model_paths = find_raw_models(
            output_dir=output_dir,
            target_id=target_id,
        )

        _best_raw_model_path, modeller_model_path = standardize_model_name(
            output_dir=output_dir,
            target_id=target_id,
            raw_model_paths=raw_model_paths,
            final_name=final_model_name,
        )

        append_key_value_block(
            log_path=module_log_path,
            title="MODELLER OUTPUT",
            key_value_pairs=[
                ("effective_residue_range", effective_residue_range),
                ("alignment_file", str(alignment_file)),
                ("script_path", str(script_path)),
                ("stdout_log", str(stdout_log)),
                ("stderr_log", str(stderr_log)),
                ("n_raw_models", str(len(raw_model_paths))),
                ("final_modeller_model_path", str(modeller_model_path)),
            ],
        )

        _trim_final_model_in_place(modeller_model_path, effective_residue_range)
        _log_debug(
            module_log_path,
            f"Trimmed MODELLER model to range {effective_residue_range!r}: {modeller_model_path}",
        )
        _log_debug(module_log_path, f"Final MODELLER model path: {modeller_model_path}")
        _log_debug(module_log_path, "=== run_filler_for_chain END (MODELLER) ===")

        result = _build_result(
            chain_id=resolved_chain_id,
            output_dir=output_dir,
            alignment_file=alignment_file,
            template_pdb=copied_template_pdb,
            script_file=script_path,
            modeller_model_path=modeller_model_path,
            alphafold_model_path=None,
            final_model_path=modeller_model_path,
            raw_model_paths=raw_model_paths,
            stdout_log=stdout_log,
            stderr_log=stderr_log,
            skipped=False,
            skip_reason=None,
            fill_decision=fill_decision,
        )

        _print_end_box(
            pdb_id=pdb_id,
            resolved_chain_id=resolved_chain_id,
            status="success_modeller",
            final_model_path=modeller_model_path,
            modeller_model_path=modeller_model_path,
            alphafold_model_path=None,
            stdout_log=stdout_log,
            stderr_log=stderr_log,
            skip_reason=None,
        )
        return result

    except Exception as error:
        append_log_text(module_log_path, "[RESULT] FAILED")
        append_exception_traceback(module_log_path, error)
        print_double_header_box(
            title=f"{MODULE_NAME} · END",
            line_list=[
                f"PDB ID           : {pdb_id}",
                f"Resolved chain   : {resolved_chain_id if 'resolved_chain_id' in locals() else '(unresolved)'}",
                "Status           : failed",
                f"Error type       : {type(error).__name__}",
                f"Error            : {error}",
                f"Module log       : {module_log_path}",
            ],
        )
        raise