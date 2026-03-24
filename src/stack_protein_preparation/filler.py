# /home/grheco/repositorios/stack_protein_prep/src/stack_protein_preparation/filler.py

from __future__ import annotations

import json
import os
import re
import subprocess
from dataclasses import dataclass
from pathlib import Path
from urllib.error import HTTPError, URLError
from urllib.request import urlopen, urlretrieve

PYTHON_BIN = Path("/usr/bin/python3")
MODELLER_PYTHONPATH = (
    "/usr/lib/modeller10.8/modlib:/usr/lib/modeller10.8/lib/x86_64-intel8"
)

THREE_TO_ONE_RESIDUE_CODE_MAP = {
    "ALA": "A",
    "ARG": "R",
    "ASN": "N",
    "ASP": "D",
    "ASH": "D",
    "CYS": "C",
    "CYM": "C",
    "CYX": "C",
    "GLN": "Q",
    "GLU": "E",
    "GLH": "E",
    "GLY": "G",
    "HIS": "H",
    "HID": "H",
    "HIE": "H",
    "HIP": "H",
    "ILE": "I",
    "LEU": "L",
    "LYS": "K",
    "LYN": "K",
    "MET": "M",
    "MSE": "M",
    "PHE": "F",
    "PRO": "P",
    "PYL": "O",
    "SEC": "U",
    "SER": "S",
    "THR": "T",
    "TRP": "W",
    "TYR": "Y",
    "VAL": "V",
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
    final_model_path: Path | None
    raw_model_paths: tuple[Path, ...]
    stdout_log: Path
    stderr_log: Path
    skipped: bool
    skip_reason: str | None
    fill_decision: FillDecision


def _debug(message: str) -> None:
    print(f"[filler] {message}")


def _classify_gap_length(gap_length: int) -> str:
    if 1 <= gap_length <= 5:
        return "green"
    if 6 <= gap_length <= 8:
        return "yellow"
    return "alphafold_candidate"


def _find_gap_regions_in_sequence(
    aligned_sequence: str,
) -> tuple[GapRegion, ...]:
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


def _get_internal_gap_regions(
    alignment_file: Path,
    template_id: str,
) -> tuple[GapRegion, ...]:
    template_sequence = _extract_sequence_block_from_ali(
        alignment_file=alignment_file,
        record_id=template_id,
    )

    all_gap_regions = _find_gap_regions_in_sequence(template_sequence)

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
    return internal_gap_regions


def analyze_fill_decision(
    alignment_file: Path,
    template_id: str,
) -> FillDecision:
    gap_regions = _get_internal_gap_regions(
        alignment_file=alignment_file,
        template_id=template_id,
    )

    if not gap_regions:
        return FillDecision(
            should_run_modeller=False,
            overall_classification="no_internal_gaps",
            gap_regions=(),
            skip_reason="Template alignment contains no internal gaps.",
            alphafold_candidate=False,
        )

    has_large_gap = any(gap_region.gap_length > 8 for gap_region in gap_regions)
    if has_large_gap:
        return FillDecision(
            should_run_modeller=False,
            overall_classification="alphafold_candidate",
            gap_regions=gap_regions,
            skip_reason=(
                "At least one internal gap is longer than 8 residues. "
                "Using AlphaFold fallback if possible."
            ),
            alphafold_candidate=True,
        )

    has_yellow_gap = any(6 <= gap_region.gap_length <= 8 for gap_region in gap_regions)
    if has_yellow_gap:
        return FillDecision(
            should_run_modeller=True,
            overall_classification="yellow",
            gap_regions=gap_regions,
            skip_reason=None,
            alphafold_candidate=False,
        )

    return FillDecision(
        should_run_modeller=True,
        overall_classification="green",
        gap_regions=gap_regions,
        skip_reason=None,
        alphafold_candidate=False,
    )


def find_alignment_fasta_for_filler(
    alignment_directory: Path,
) -> Path:
    """
    Use only ATOM-vs-UniProt alignment for filling.
    """
    if not alignment_directory.exists():
        raise FileNotFoundError(f"Alignment directory not found: {alignment_directory}")

    preferred_pattern_list = [
        "*ATOM*chain*vs*UniProt*.aln.fasta",
        "*ATOM*vs*UniProt*.aln.fasta",
        "*ATOM*.aln.fasta",
    ]

    for pattern in preferred_pattern_list:
        candidate_list = tuple(sorted(alignment_directory.glob(pattern)))
        if candidate_list:
            _debug(f"Selected ATOM alignment FASTA: {candidate_list[0]}")
            return candidate_list[0]

    raise FileNotFoundError(
        f"No ATOM-vs-UniProt alignment FASTA found in {alignment_directory}"
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


def split_template_and_target_alignment_records(
    alignment_fasta_path: Path,
) -> tuple[str, str, str, str]:
    """
    Return:
    -------
    template_header, template_alignment_skeleton, target_header, target_aligned_sequence
    """
    (header_1, seq_1), (header_2, seq_2) = read_two_sequence_fasta(alignment_fasta_path)

    if _is_pdb_atom_header(header_1) and _is_uniprot_header(header_2):
        return header_1, seq_1, header_2, seq_2

    if _is_pdb_atom_header(header_2) and _is_uniprot_header(header_1):
        return header_2, seq_2, header_1, seq_1

    raise ValueError(
        "Could not determine template/target rows from FASTA headers: "
        f"{header_1!r}, {header_2!r}"
    )


def write_chain_specific_template_pdb(
    template_pdb_path: Path,
    output_dir: Path,
    template_id: str,
    chain_id: str,
) -> Path:
    """
    Write a chain-specific template PDB for MODELLER.

    Only ATOM records from the selected chain are kept.
    """
    if not template_pdb_path.exists():
        raise FileNotFoundError(f"Template PDB not found: {template_pdb_path}")

    output_dir.mkdir(parents=True, exist_ok=True)
    destination = output_dir / f"{template_id}.pdb"

    selected_line_list: list[str] = []

    for line in template_pdb_path.read_text(encoding="utf-8").splitlines():
        if not line.startswith("ATOM"):
            continue

        pdb_chain_id = line[21].strip() or "_"
        if pdb_chain_id != chain_id:
            continue

        selected_line_list.append(line)

    if not selected_line_list:
        raise ValueError(
            f"No ATOM records found for chain {chain_id!r} in {template_pdb_path}"
        )

    selected_line_list.append("TER")
    selected_line_list.append("END")

    destination.write_text("\n".join(selected_line_list) + "\n", encoding="utf-8")

    _debug(f"Wrote chain-specific template PDB for chain {chain_id} to: {destination}")
    return destination


def extract_sequence_from_template_pdb(
    template_pdb_path: Path,
) -> str:
    """
    Extract one-letter protein sequence directly from a chain-specific template PDB.

    Filtering rule
    --------------
    A residue is only treated as part of the protein sequence if it has a
    peptide-like backbone:
    - N
    - CA
    - C

    Mapping rule
    ------------
    - known amino-acid-like residue names are mapped via THREE_TO_ONE_RESIDUE_CODE_MAP
    - unknown but peptide-like residues are mapped to X
    - non-peptide residues are skipped completely
    """
    if not template_pdb_path.exists():
        raise FileNotFoundError(f"Template PDB not found: {template_pdb_path}")

    residue_atom_name_set_by_key: dict[tuple[str, str, str], set[str]] = {}
    residue_name_by_key: dict[tuple[str, str, str], str] = {}

    for line in template_pdb_path.read_text(encoding="utf-8").splitlines():
        if not line.startswith("ATOM"):
            continue

        residue_name = line[17:20].strip().upper()
        pdb_chain_id = line[21].strip() or "_"
        residue_number = line[22:26].strip()
        insertion_code = line[26].strip()
        atom_name = line[12:16].strip().upper()

        residue_key = (pdb_chain_id, residue_number, insertion_code)

        if residue_key not in residue_atom_name_set_by_key:
            residue_atom_name_set_by_key[residue_key] = set()

        residue_atom_name_set_by_key[residue_key].add(atom_name)
        residue_name_by_key[residue_key] = residue_name

    if not residue_atom_name_set_by_key:
        raise ValueError(f"No ATOM residues found in template PDB: {template_pdb_path}")

    sorted_residue_key_list = sorted(
        residue_atom_name_set_by_key.keys(),
        key=lambda key: (
            key[0],
            int(key[1]) if key[1].lstrip("-").isdigit() else 10**9,
            key[2],
        ),
    )

    residue_sequence_list: list[str] = []
    skipped_residue_debug_list: list[str] = []

    for residue_key in sorted_residue_key_list:
        atom_name_set = residue_atom_name_set_by_key[residue_key]
        residue_name = residue_name_by_key[residue_key]

        has_peptide_backbone = {"N", "CA", "C"}.issubset(atom_name_set)

        if not has_peptide_backbone:
            skipped_residue_debug_list.append(
                f"{residue_key[0]}:{residue_key[1]}{residue_key[2]}:{residue_name}"
            )
            continue

        residue_one_letter = THREE_TO_ONE_RESIDUE_CODE_MAP.get(residue_name, "X")
        residue_sequence_list.append(residue_one_letter)

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


def build_modeller_template_alignment_sequence(
    aligned_template_skeleton: str,
    actual_template_sequence: str,
) -> str:
    """
    Rebuild the template alignment line using the real PDB sequence.

    Rule
    ----
    - keep '-' positions from the alignment skeleton
    - for non-gap positions, consume residues from the real template PDB sequence
    - once the real template sequence is exhausted:
        - trailing 'X' placeholder positions are converted to '-'
        - any other non-gap character is treated as an error
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
) -> Path:
    """
    Build a MODELLER .ali file using:
    - template alignment skeleton from the existing ATOM-vs-UniProt alignment
    - real template sequence extracted directly from the chain-specific template PDB
    - target aligned sequence from the existing alignment
    """
    (
        template_header,
        template_alignment_skeleton,
        target_header,
        target_aligned_sequence,
    ) = split_template_and_target_alignment_records(alignment_fasta_path)

    _debug(f"Template FASTA header selected: {template_header}")
    _debug(f"Target FASTA header selected: {target_header}")

    actual_template_sequence = extract_sequence_from_template_pdb(
        template_pdb_path=template_pdb_path,
    )

    rebuilt_template_aligned_sequence = build_modeller_template_alignment_sequence(
        aligned_template_skeleton=template_alignment_skeleton,
        actual_template_sequence=actual_template_sequence,
    )

    output_dir.mkdir(parents=True, exist_ok=True)
    alignment_file = output_dir / f"{target_id}.ali"

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
    _debug(
        f"Template ungapped length in .ali: "
        f"{len(rebuilt_template_aligned_sequence.replace('-', ''))}"
    )
    _debug(
        f"Target ungapped length in .ali: "
        f"{len(target_aligned_sequence.replace('-', ''))}"
    )

    return alignment_file


def write_modeller_script(
    output_dir: Path,
    alignment_file: Path,
    template_id: str,
    target_id: str,
    starting_model: int = 1,
    ending_model: int = 1,
) -> Path:
    script_path = output_dir / "model.py"

    content = f"""from modeller import *
from modeller.automodel import *

log.none()

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

score_file = 'model_scores.tsv'
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

    stdout_log = working_dir / "modeller_stdout.log"
    stderr_log = working_dir / "modeller_stderr.log"

    stdout_log.write_text(result.stdout, encoding="utf-8")
    stderr_log.write_text(result.stderr, encoding="utf-8")

    _debug(f"MODELLER return code: {result.returncode}")
    _debug("===== MODELLER STDOUT BEGIN =====")
    print(result.stdout)
    _debug("===== MODELLER STDOUT END =====")
    _debug("===== MODELLER STDERR BEGIN =====")
    print(result.stderr)
    _debug("===== MODELLER STDERR END =====")

    if result.returncode != 0:
        raise RuntimeError(
            f"MODELLER failed for script {script_path.name}. "
            f"See {stdout_log} and {stderr_log}."
        )

    return stdout_log, stderr_log


def find_raw_models(
    output_dir: Path,
    target_id: str,
) -> tuple[Path, ...]:
    models = tuple(sorted(output_dir.glob(f"{target_id}*.pdb")))
    _debug(f"Found {len(models)} raw model(s) for target {target_id}")
    return models


def select_best_model_from_scores(
    output_dir: Path,
) -> Path:
    score_file = output_dir / "model_scores.tsv"

    if not score_file.exists():
        raise FileNotFoundError(f"Score file not found: {score_file}")

    best_model_name: str | None = None
    best_dope_score: float | None = None

    lines = score_file.read_text(encoding="utf-8").splitlines()
    if len(lines) < 2:
        raise ValueError(f"Score file is empty or invalid: {score_file}")

    for line in lines[1:]:
        stripped = line.strip()
        if not stripped:
            continue

        parts = stripped.split("\t")
        if len(parts) != 3:
            raise ValueError(f"Invalid score line in {score_file}: {line}")

        model_name, dope_score_text, ga341_score_text = parts

        if dope_score_text in {"None", "", "nan"}:
            continue

        dope_score = float(dope_score_text)

        _debug(
            f"Candidate model: {model_name}, "
            f"DOPE={dope_score_text}, GA341={ga341_score_text}"
        )

        if best_dope_score is None or dope_score < best_dope_score:
            best_dope_score = dope_score
            best_model_name = model_name

    if best_model_name is None:
        raise ValueError(f"No valid DOPE scores found in {score_file}")

    best_model_path = output_dir / best_model_name
    if not best_model_path.exists():
        raise FileNotFoundError(f"Best-scoring model file not found: {best_model_path}")

    _debug(
        f"Best model selected by lowest DOPE score: "
        f"{best_model_path} (DOPE={best_dope_score})"
    )
    return best_model_path


def cleanup_model_pdb(
    input_model_path: Path,
    output_model_path: Path,
) -> Path:
    """
    Write a cleaned protein-only final PDB.
    """
    if not input_model_path.exists():
        raise FileNotFoundError(f"Model PDB not found: {input_model_path}")

    cleaned_line_list: list[str] = []

    for line in input_model_path.read_text(encoding="utf-8").splitlines():
        if line.startswith("ATOM"):
            cleaned_line_list.append(line)

    if not cleaned_line_list:
        raise ValueError(f"No ATOM records found in model PDB: {input_model_path}")

    cleaned_line_list.append("TER")
    cleaned_line_list.append("END")

    output_model_path.write_text(
        "\n".join(cleaned_line_list) + "\n",
        encoding="utf-8",
    )

    _debug(f"Cleaned final model written to: {output_model_path}")
    return output_model_path


def standardize_model_name(
    output_dir: Path,
    final_name: str,
) -> Path:
    best_model = select_best_model_from_scores(output_dir=output_dir)
    final_model = output_dir / final_name

    cleaned_final_model = cleanup_model_pdb(
        input_model_path=best_model,
        output_model_path=final_model,
    )

    _debug(f"Final cleaned model path: {cleaned_final_model}")
    return cleaned_final_model


def download_alphafold_structure(
    uniprot_id: str,
    output_dir: Path,
) -> Path | None:
    """
    Download an AlphaFold model for one UniProt accession via the official AFDB API.

    Returns
    -------
    Path | None
        Local path to the downloaded PDB file if available, else None.
    """
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

    target_path = output_dir / f"AF-{uniprot_id}-F1-model.pdb"
    urlretrieve(pdb_url, target_path)

    if not target_path.is_file() or target_path.stat().st_size == 0:
        return None

    return target_path


def align_protonated_alphafold_model_to_start_pdb(
    reference_pdb_path: Path,
    mobile_pdb_path: Path,
    output_pdb_path: Path,
) -> dict[str, str | float | bool]:
    """
    Align an AlphaFold model to a reference PDB using CA atoms.

    Parameters
    ----------
    reference_pdb_path
        Original structure used as the spatial reference.
    mobile_pdb_path
        AlphaFold model PDB to be aligned before protonation.
    output_pdb_path
        Output aligned structure. May be equal to mobile_pdb_path.

    Returns
    -------
    dict
        alignment metadata
    """
    from Bio.PDB import PDBIO, PDBParser, Superimposer

    if not reference_pdb_path.exists():
        raise FileNotFoundError(reference_pdb_path)

    if not mobile_pdb_path.exists():
        raise FileNotFoundError(mobile_pdb_path)

    parser = PDBParser(QUIET=True)

    reference_structure = parser.get_structure("ref", str(reference_pdb_path))
    mobile_structure = parser.get_structure("mob", str(mobile_pdb_path))

    def _collect_protein_ca_atoms(structure) -> list:
        ca_atom_list: list = []

        for model in structure:
            for chain in model:
                for residue in chain:
                    if residue.id[0] != " ":
                        continue
                    if "CA" in residue:
                        ca_atom_list.append(residue["CA"])

        return ca_atom_list

    reference_ca_atom_list = _collect_protein_ca_atoms(reference_structure)
    mobile_ca_atom_list = _collect_protein_ca_atoms(mobile_structure)

    if not reference_ca_atom_list:
        raise ValueError(
            f"No protein CA atoms found in reference: {reference_pdb_path}"
        )

    if not mobile_ca_atom_list:
        raise ValueError(f"No protein CA atoms found in mobile: {mobile_pdb_path}")

    n_matched_atoms = min(len(reference_ca_atom_list), len(mobile_ca_atom_list))

    if n_matched_atoms < 3:
        raise ValueError(
            "Not enough CA atoms for reliable alignment: "
            f"reference={len(reference_ca_atom_list)}, mobile={len(mobile_ca_atom_list)}"
        )

    fixed_atom_list = reference_ca_atom_list[:n_matched_atoms]
    moving_atom_list = mobile_ca_atom_list[:n_matched_atoms]

    superimposer = Superimposer()
    superimposer.set_atoms(fixed_atom_list, moving_atom_list)
    superimposer.apply(mobile_structure.get_atoms())

    output_pdb_path.parent.mkdir(parents=True, exist_ok=True)

    io = PDBIO()
    io.set_structure(mobile_structure)
    io.save(str(output_pdb_path))

    return {
        "alignment_success": output_pdb_path.exists()
        and output_pdb_path.stat().st_size > 0,
        "alignment_rmsd": float(superimposer.rms),
        "alignment_output_path": str(output_pdb_path),
    }


def crop_pdb_to_range(
    input_pdb_path: Path,
    output_pdb_path: Path,
    residue_range: str,
) -> Path:
    """
    Crop a PDB file to the inclusive residue range 'start-end'.

    Notes
    -----
    - keeps ATOM records only
    - keeps residues whose residue number is within the requested range
    - writes TER and END
    """
    if input_pdb_path is None:
        raise FileNotFoundError("AlphaFold download returned no PDB path")

    if not input_pdb_path.exists():
        raise FileNotFoundError(f"Input PDB not found: {input_pdb_path}")

    residue_range_text = str(residue_range).strip()
    match = re.fullmatch(r"\s*(\d+)\s*-\s*(\d+)\s*", residue_range_text)
    if match is None:
        raise ValueError(f"Invalid residue range: {residue_range!r}")

    start_residue = int(match.group(1))
    end_residue = int(match.group(2))

    if start_residue > end_residue:
        raise ValueError(f"Invalid residue range: start > end in {residue_range!r}")

    kept_line_list: list[str] = []

    for line in input_pdb_path.read_text(encoding="utf-8").splitlines():
        if not line.startswith("ATOM"):
            continue

        residue_number_text = line[22:26].strip()
        if not residue_number_text.lstrip("-").isdigit():
            continue

        residue_number = int(residue_number_text)

        if start_residue <= residue_number <= end_residue:
            kept_line_list.append(line)

    if not kept_line_list:
        raise ValueError(
            f"No ATOM records remained after cropping {input_pdb_path} "
            f"to range {residue_range!r}"
        )

    kept_line_list.append("TER")
    kept_line_list.append("END")

    output_pdb_path.write_text(
        "\n".join(kept_line_list) + "\n",
        encoding="utf-8",
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
    """
    Download AlphaFold PDB, crop it to the requested CSV range,
    align it to the starting template PDB, and write a cleaned final PDB.
    """
    _ = model_version  # intentionally unused for AFDB API path

    alphafold_dir = output_dir / "alphafold"

    downloaded_pdb = download_alphafold_structure(
        uniprot_id=uniprot_id,
        output_dir=alphafold_dir,
    )

    if downloaded_pdb is None:
        raise FileNotFoundError(
            f"No AlphaFold PDB available for UniProt {uniprot_id!r}"
        )

    cropped_path = output_dir / f"cropped_{final_model_name}"
    cropped_model_path = crop_pdb_to_range(
        input_pdb_path=downloaded_pdb,
        output_pdb_path=cropped_path,
        residue_range=residue_range,
    )

    alignment_result = align_protonated_alphafold_model_to_start_pdb(
        reference_pdb_path=template_pdb_path,
        mobile_pdb_path=cropped_model_path,
        output_pdb_path=cropped_model_path,
    )

    _debug(
        "AlphaFold alignment result: "
        f"success={alignment_result['alignment_success']}, "
        f"rmsd={alignment_result['alignment_rmsd']}, "
        f"output={alignment_result['alignment_output_path']}"
    )

    final_model_path = output_dir / final_model_name
    cleaned_model_path = cleanup_model_pdb(
        input_model_path=cropped_model_path,
        output_model_path=final_model_path,
    )

    return cleaned_model_path


def run_filler_for_chain(
    alignment_directory: Path,
    template_pdb_path: Path,
    output_dir: Path,
    template_id: str,
    target_id: str,
    chain_id: str,
    final_model_name: str | None = None,
    starting_model: int = 1,
    ending_model: int = 1,
    skip_if_no_internal_gaps: bool = True,
    uniprot_id: str | None = None,
    residue_range: str = "",
) -> FillerRunResult:
    _debug("=== run_filler_for_chain START ===")
    _debug(f"alignment_directory: {alignment_directory}")
    _debug(f"template_pdb_path: {template_pdb_path}")
    _debug(f"output_dir: {output_dir}")
    _debug(f"template_id: {template_id}")
    _debug(f"target_id: {target_id}")
    _debug(f"chain_id: {chain_id}")
    _debug(f"starting_model: {starting_model}")
    _debug(f"ending_model: {ending_model}")

    if not alignment_directory.exists():
        raise FileNotFoundError(f"Alignment directory not found: {alignment_directory}")

    if not template_pdb_path.exists():
        raise FileNotFoundError(f"Template PDB not found: {template_pdb_path}")

    output_dir.mkdir(parents=True, exist_ok=True)

    if final_model_name is None:
        final_model_name = f"{target_id}_protein_mod.pdb"

    alignment_fasta_path = find_alignment_fasta_for_filler(alignment_directory)
    _debug(f"Alignment FASTA selected for filler: {alignment_fasta_path}")

    copied_template_pdb = write_chain_specific_template_pdb(
        template_pdb_path=template_pdb_path,
        output_dir=output_dir,
        template_id=template_id,
        chain_id=chain_id,
    )

    alignment_file = write_modeller_alignment_from_existing_alignment(
        alignment_fasta_path=alignment_fasta_path,
        template_pdb_path=copied_template_pdb,
        output_dir=output_dir,
        template_id=template_id,
        target_id=target_id,
    )

    fill_decision = analyze_fill_decision(
        alignment_file=alignment_file,
        template_id=template_id,
    )

    _debug(
        "Fill decision: "
        f"should_run_modeller={fill_decision.should_run_modeller}, "
        f"overall_classification={fill_decision.overall_classification}, "
        f"alphafold_candidate={fill_decision.alphafold_candidate}"
    )

    if (
        skip_if_no_internal_gaps
        and fill_decision.overall_classification == "no_internal_gaps"
    ):
        stdout_log = output_dir / "modeller_stdout.log"
        stderr_log = output_dir / "modeller_stderr.log"
        stdout_log.write_text("", encoding="utf-8")
        stderr_log.write_text("", encoding="utf-8")

        _debug("No internal gaps found -> skipping filler run")
        _debug("=== run_filler_for_chain END (SKIPPED) ===")

        return FillerRunResult(
            chain_id=chain_id,
            output_dir=output_dir,
            alignment_file=alignment_file,
            template_pdb=copied_template_pdb,
            script_file=output_dir / "model.py",
            final_model_path=None,
            raw_model_paths=(),
            stdout_log=stdout_log,
            stderr_log=stderr_log,
            skipped=True,
            skip_reason=fill_decision.skip_reason,
            fill_decision=fill_decision,
        )

    if not fill_decision.should_run_modeller:
        stdout_log = output_dir / "modeller_stdout.log"
        stderr_log = output_dir / "modeller_stderr.log"
        stdout_log.write_text("", encoding="utf-8")
        stderr_log.write_text("", encoding="utf-8")

        if fill_decision.alphafold_candidate:
            if not uniprot_id:
                _debug("AlphaFold fallback needed but UniProt ID is missing")
                _debug("=== run_filler_for_chain END (SKIPPED) ===")
                return FillerRunResult(
                    chain_id=chain_id,
                    output_dir=output_dir,
                    alignment_file=alignment_file,
                    template_pdb=copied_template_pdb,
                    script_file=output_dir / "model.py",
                    final_model_path=None,
                    raw_model_paths=(),
                    stdout_log=stdout_log,
                    stderr_log=stderr_log,
                    skipped=True,
                    skip_reason=(
                        "AlphaFold fallback required, but UniProt ID is missing."
                    ),
                    fill_decision=fill_decision,
                )

            alphafold_final_model_path = run_alphafold_fallback_for_chain(
                output_dir=output_dir,
                template_pdb_path=template_pdb_path,
                uniprot_id=uniprot_id,
                residue_range=residue_range,
                final_model_name=final_model_name,
                model_version=4,
            )

            _debug("=== run_filler_for_chain END (ALPHAFOLD FALLBACK) ===")
            return FillerRunResult(
                chain_id=chain_id,
                output_dir=output_dir,
                alignment_file=alignment_file,
                template_pdb=copied_template_pdb,
                script_file=output_dir / "model.py",
                final_model_path=alphafold_final_model_path,
                raw_model_paths=(),
                stdout_log=stdout_log,
                stderr_log=stderr_log,
                skipped=False,
                skip_reason=None,
                fill_decision=fill_decision,
            )

        _debug("MODELLER not run for this chain")
        _debug("=== run_filler_for_chain END (SKIPPED) ===")

        return FillerRunResult(
            chain_id=chain_id,
            output_dir=output_dir,
            alignment_file=alignment_file,
            template_pdb=copied_template_pdb,
            script_file=output_dir / "model.py",
            final_model_path=None,
            raw_model_paths=(),
            stdout_log=stdout_log,
            stderr_log=stderr_log,
            skipped=True,
            skip_reason=fill_decision.skip_reason,
            fill_decision=fill_decision,
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

    final_model_path = standardize_model_name(
        output_dir=output_dir,
        final_name=final_model_name,
    )

    _debug(f"Final model path: {final_model_path}")
    _debug("=== run_filler_for_chain END ===")

    return FillerRunResult(
        chain_id=chain_id,
        output_dir=output_dir,
        alignment_file=alignment_file,
        template_pdb=copied_template_pdb,
        script_file=script_path,
        final_model_path=final_model_path,
        raw_model_paths=raw_model_paths,
        stdout_log=stdout_log,
        stderr_log=stderr_log,
        skipped=False,
        skip_reason=None,
        fill_decision=fill_decision,
    )
