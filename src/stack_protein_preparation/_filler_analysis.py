"""Gap analysis, alignment parsing, and chain-selection helpers for filler."""
from __future__ import annotations

import re
from pathlib import Path
from typing import Any

from Bio.PDB import PDBIO

from ._filler_shared import (
    FillDecision,
    GapRegion,
    _FirstModelChainProteinSelect,
    _classify_gap_length,
    _debug,
    _is_protein_atom_residue,
    _is_supported_protein_residue_name,
    _normalize_chain_id,
    _parse_first_model,
    _parse_residue_range,
    _residue_has_peptide_backbone,
    _residue_name_to_one_letter,
    _residue_to_one_letter,
    _resolve_residue_range_for_filler,
    _strip_gaps,
    _validate_written_pdb_has_atoms,
)


def _find_gap_regions_in_sequence(aligned_sequence: str) -> tuple[GapRegion, ...]:
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


def _extract_sequence_block_from_ali(alignment_file: Path, record_id: str) -> str:
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
    _debug(f"Internal gap regions: {len(internal_gap_regions)}")

    if not internal_gap_regions:
        return FillDecision(
            should_run_modeller=False,
            overall_classification="no_internal_gaps",
            gap_regions=(),
            skip_reason="Template alignment contains no internal gaps.",
            alphafold_candidate=False,
        )

    has_alphafold_gap = any(gap_region.gap_length >= 8 for gap_region in internal_gap_regions)

    # Domain-only-vs-full-length guard.
    # When a crystal that resolves only one domain is aligned against the
    # full-length UniProt target, short spurious matches outside the
    # crystallised domain scatter the template residues across the whole target.
    # The stretches between those spurious fragments then look like very large
    # "internal" gaps even though they are simply un-crystallised regions of the
    # full-length protein, not missing loops to rebuild. Filling them (especially
    # via AlphaFold, which crops to the spurious span) corrupts the model — e.g.
    # 3OLL produced a 1072-residue model cropped to 370-1441 instead of the true
    # 262-498 domain. Detect this when the missing (gap) content exceeds the
    # residues the crystal actually provides, and skip filling: the downstream
    # range-trim then keeps the genuine crystal domain.
    template_residue_count = sum(1 for c in template_alignment_skeleton if c != "-")
    total_internal_gap = sum(g.gap_length for g in internal_gap_regions)
    if has_alphafold_gap and total_internal_gap > template_residue_count:
        return FillDecision(
            should_run_modeller=False,
            overall_classification="domain_only_spread",
            gap_regions=internal_gap_regions,
            skip_reason=(
                "Template resolves only part of the full-length target "
                f"(crystallised residues={template_residue_count}, internal gap "
                f"content={total_internal_gap}). The large gaps are un-crystallised "
                "regions, not missing loops, so filling is skipped to avoid grafting "
                "over residues outside the resolved domain."
            ),
            alphafold_candidate=False,
        )
    if has_alphafold_gap:
        has_small_gap = any(gap_region.gap_length < 8 for gap_region in internal_gap_regions)
        return FillDecision(
            should_run_modeller=has_small_gap,
            overall_classification="alphafold_candidate",
            gap_regions=internal_gap_regions,
            skip_reason=(
                "At least one internal gap is 8 residues or longer — AlphaFold candidate. "
                + ("Modeller will also fill the smaller gap(s)." if has_small_gap else "")
            ),
            alphafold_candidate=True,
        )

    has_yellow_gap = any(6 <= gap_region.gap_length < 8 for gap_region in internal_gap_regions)
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


def analyze_fill_decision(alignment_file: Path, template_id: str) -> FillDecision:
    template_sequence = _extract_sequence_block_from_ali(
        alignment_file=alignment_file, record_id=template_id
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
    alignment_fasta_path = find_alignment_fasta_for_filler(
        alignment_directory=alignment_directory,
        chain_id=chain_id,
    )
    (_, template_alignment_skeleton), _ = read_two_sequence_fasta(alignment_fasta_path)
    all_gap_regions = _find_gap_regions_in_sequence(template_alignment_skeleton)
    internal_gap_regions = [g for g in all_gap_regions if not g.is_terminal]
    return {
        "n_internal_gaps": len(internal_gap_regions),
        "max_internal_gap": max((g.gap_length for g in internal_gap_regions), default=0),
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
            residue_numbers_by_chain.setdefault(chain_id, set()).add(int(residue.id[1]))
    return residue_numbers_by_chain


def choose_filler_chain_id_from_range_and_alignments(
    template_pdb_path: Path,
    alignment_directory: Path,
    residue_range: str,
) -> str:
    if not template_pdb_path.exists():
        raise FileNotFoundError(f"Template PDB not found: {template_pdb_path}")
    if not alignment_directory.exists():
        raise FileNotFoundError(f"Alignment directory not found: {alignment_directory}")

    available_chain_ids = set(list_available_atom_alignment_chain_ids(alignment_directory))
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
        return best_rows[0][0]

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
        return _normalize_chain_id(chain_id)
    return choose_filler_chain_id_from_range_and_alignments(
        template_pdb_path=template_pdb_path,
        alignment_directory=alignment_directory,
        residue_range=residue_range,
    )


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
            return candidate

    candidate_list = tuple(
        sorted(alignment_directory.glob(f"*chain*{normalized_chain_id}*ATOM*vs*UniProt*.aln.fasta"))
    )
    if len(candidate_list) == 1:
        return candidate_list[0]

    candidate_list = tuple(sorted(alignment_directory.glob("*ATOM*vs*UniProt*.aln.fasta")))
    if len(candidate_list) == 1:
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
    return destination


def extract_sequence_from_template_pdb(template_pdb_path: Path) -> str:
    if not template_pdb_path.exists():
        raise FileNotFoundError(f"Template PDB not found: {template_pdb_path}")
    parsed = _parse_first_model(template_pdb_path)
    residue_sequence_list: list[str] = []
    for chain in parsed.model:
        for residue in chain:
            if not _is_protein_atom_residue(residue):
                continue
            if not _residue_has_peptide_backbone(residue):
                continue
            residue_sequence_list.append(_residue_to_one_letter(residue))
    if not residue_sequence_list:
        raise ValueError(
            f"No peptide-like ATOM residues found in template PDB: {template_pdb_path}"
        )
    return "".join(residue_sequence_list)


def _first_sequence_mismatch(
    sequence_a: str,
    sequence_b: str,
) -> tuple[int, str, str] | None:
    shortest_length = min(len(sequence_a), len(sequence_b))
    for index in range(shortest_length):
        if sequence_a[index] != sequence_b[index]:
            return index, sequence_a[index], sequence_b[index]
    if len(sequence_a) != len(sequence_b):
        a_char = sequence_a[shortest_length] if len(sequence_a) > shortest_length else "<END>"
        b_char = sequence_b[shortest_length] if len(sequence_b) > shortest_length else "<END>"
        return shortest_length, a_char, b_char
    return None


def _write_sequence_debug_files(
    output_dir: Path,
    template_alignment_skeleton: str,
    actual_template_sequence: str,
) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    (output_dir / "debug_template_alignment_skeleton.txt").write_text(
        template_alignment_skeleton + "\n", encoding="utf-8"
    )
    (output_dir / "debug_template_alignment_ungapped.txt").write_text(
        _strip_gaps(template_alignment_skeleton) + "\n", encoding="utf-8"
    )
    (output_dir / "debug_template_pdb_sequence.txt").write_text(
        actual_template_sequence + "\n", encoding="utf-8"
    )


def _normalize_alignment_template_for_validation(template_alignment_skeleton: str) -> str:
    ungapped = _strip_gaps(template_alignment_skeleton)
    ungapped = re.sub(r"^X+", "", ungapped)
    ungapped = re.sub(r"X+$", "", ungapped)
    return ungapped


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
    mismatch = _first_sequence_mismatch(normalized_alignment_template, actual_template_sequence)
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

    return "".join(rebuilt_character_list)
