# /home/grheco/repositorios/stack_protein_prep/tests/test_filler.py

from __future__ import annotations

from pathlib import Path

import pytest

from stack_protein_preparation.filler import (
    GapRegion,
    _extract_chain_id_from_alignment_header,
    _find_gap_regions_in_sequence,
    _normalize_chain_id,
    _write_skip_logs,
    analyze_fill_decision,
    build_modeller_template_alignment_sequence,
    cleanup_model_pdb,
    extract_sequence_from_template_pdb,
    find_alignment_fasta_for_filler,
    find_raw_models,
    read_two_sequence_fasta,
    select_best_model_from_scores,
    split_template_and_target_alignment_records,
    validate_template_sequence_consistency,
    write_chain_specific_template_pdb,
    write_modeller_alignment_from_existing_alignment,
)


def _write_text(path: Path, text: str) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")
    return path


def _make_atom_line(
    serial: int,
    atom_name: str,
    residue_name: str,
    chain_id: str,
    residue_number: int,
    x: float,
    y: float,
    z: float,
    element: str,
) -> str:
    return (
        f"ATOM  {serial:5d} {atom_name:>4} {residue_name:>3} {chain_id:1}"
        f"{residue_number:4d}    "
        f"{x:8.3f}{y:8.3f}{z:8.3f}"
        f"{1.00:6.2f}{20.00:6.2f}          "
        f"{element:>2}"
    )


def _build_template_pdb_with_two_chains(path: Path) -> Path:
    lines = [
        _make_atom_line(1, "N", "ALA", "A", 1, 0.0, 0.0, 0.0, "N"),
        _make_atom_line(2, "CA", "ALA", "A", 1, 1.5, 0.0, 0.0, "C"),
        _make_atom_line(3, "C", "ALA", "A", 1, 2.7, 0.0, 0.0, "C"),
        _make_atom_line(4, "N", "GLY", "A", 2, 3.8, 0.0, 0.0, "N"),
        _make_atom_line(5, "CA", "GLY", "A", 2, 5.2, 0.0, 0.0, "C"),
        _make_atom_line(6, "C", "GLY", "A", 2, 6.4, 0.0, 0.0, "C"),
        _make_atom_line(7, "N", "SER", "B", 1, 0.0, 5.0, 0.0, "N"),
        _make_atom_line(8, "CA", "SER", "B", 1, 1.5, 5.0, 0.0, "C"),
        _make_atom_line(9, "C", "SER", "B", 1, 2.7, 5.0, 0.0, "C"),
        "TER",
        "END",
    ]
    return _write_text(path, "\n".join(lines) + "\n")


def _build_chain_a_template_pdb(path: Path) -> Path:
    lines = [
        _make_atom_line(1, "N", "ALA", "A", 1, 0.0, 0.0, 0.0, "N"),
        _make_atom_line(2, "CA", "ALA", "A", 1, 1.5, 0.0, 0.0, "C"),
        _make_atom_line(3, "C", "ALA", "A", 1, 2.7, 0.0, 0.0, "C"),
        _make_atom_line(4, "N", "GLY", "A", 2, 3.8, 0.0, 0.0, "N"),
        _make_atom_line(5, "CA", "GLY", "A", 2, 5.2, 0.0, 0.0, "C"),
        _make_atom_line(6, "C", "GLY", "A", 2, 6.4, 0.0, 0.0, "C"),
        "TER",
        "END",
    ]
    return _write_text(path, "\n".join(lines) + "\n")


def test_normalize_chain_id_returns_underscore_for_blank() -> None:
    assert _normalize_chain_id("") == "_"
    assert _normalize_chain_id("   ") == "_"
    assert _normalize_chain_id("A") == "A"


def test_find_gap_regions_in_sequence_detects_terminal_and_internal_gaps() -> None:
    gap_regions = _find_gap_regions_in_sequence("--AA---A-")

    assert gap_regions == (
        GapRegion(
            alignment_start=0,
            alignment_end=1,
            gap_length=2,
            is_terminal=True,
            classification="green",
        ),
        GapRegion(
            alignment_start=4,
            alignment_end=6,
            gap_length=3,
            is_terminal=False,
            classification="green",
        ),
        GapRegion(
            alignment_start=8,
            alignment_end=8,
            gap_length=1,
            is_terminal=True,
            classification="green",
        ),
    )


def test_analyze_fill_decision_returns_no_internal_gaps_for_terminal_only_gap(
    tmp_path: Path,
) -> None:
    ali_path = tmp_path / "test.ali"
    ali_path.write_text(
        ">P1;template\n"
        "structureX:template:FIRST:@:LAST:@::::\n"
        "AAA---*\n"
        ">P1;target\n"
        "sequence:target:FIRST:@:LAST:@::::\n"
        "AAAQQQ*\n",
        encoding="utf-8",
    )

    result = analyze_fill_decision(ali_path, "template")

    assert result.should_run_modeller is False
    assert result.overall_classification == "no_internal_gaps"
    assert result.alphafold_candidate is False
    assert result.gap_regions == ()


def test_analyze_fill_decision_returns_green_for_small_internal_gap(
    tmp_path: Path,
) -> None:
    ali_path = tmp_path / "test.ali"
    ali_path.write_text(
        ">P1;template\n"
        "structureX:template:FIRST:@:LAST:@::::\n"
        "AA---AA*\n"
        ">P1;target\n"
        "sequence:target:FIRST:@:LAST:@::::\n"
        "AAQQQAA*\n",
        encoding="utf-8",
    )

    result = analyze_fill_decision(ali_path, "template")

    assert result.should_run_modeller is True
    assert result.overall_classification == "green"
    assert result.alphafold_candidate is False
    assert len(result.gap_regions) == 1
    assert result.gap_regions[0].gap_length == 3


def test_analyze_fill_decision_returns_yellow_for_medium_internal_gap(
    tmp_path: Path,
) -> None:
    ali_path = tmp_path / "test.ali"
    ali_path.write_text(
        ">P1;template\n"
        "structureX:template:FIRST:@:LAST:@::::\n"
        "AA------AA*\n"
        ">P1;target\n"
        "sequence:target:FIRST:@:LAST:@::::\n"
        "AAQQQQQQAA*\n",
        encoding="utf-8",
    )

    result = analyze_fill_decision(ali_path, "template")

    assert result.should_run_modeller is True
    assert result.overall_classification == "yellow"
    assert result.alphafold_candidate is False
    assert len(result.gap_regions) == 1
    assert result.gap_regions[0].gap_length == 6


def test_analyze_fill_decision_returns_alphafold_candidate_for_large_internal_gap(
    tmp_path: Path,
) -> None:
    # Template covers 20 + 20 = 40 residues, internal gap is 9 — this is a
    # genuine missing-loop inside a resolved domain (gap < template content),
    # so the alphafold-candidate route is the right call.  The domain-only-spread
    # guard only fires when the gap eclipses the resolved template content.
    template_block = "A" * 20
    gap_block = "-" * 9
    target_gap_residues = "Q" * 9
    ali_path = tmp_path / "test.ali"
    ali_path.write_text(
        ">P1;template\n"
        "structureX:template:FIRST:@:LAST:@::::\n"
        f"{template_block}{gap_block}{template_block}*\n"
        ">P1;target\n"
        "sequence:target:FIRST:@:LAST:@::::\n"
        f"{template_block}{target_gap_residues}{template_block}*\n",
        encoding="utf-8",
    )

    result = analyze_fill_decision(ali_path, "template")

    assert result.should_run_modeller is False
    assert result.overall_classification == "alphafold_candidate"
    assert result.alphafold_candidate is True
    assert len(result.gap_regions) == 1
    assert result.gap_regions[0].gap_length == 9


def test_analyze_fill_decision_returns_domain_only_spread_when_gap_exceeds_template(
    tmp_path: Path,
) -> None:
    # Reproduces the 3OLL pathology: template resolves only a tiny part of the
    # full-length target (4 residues) and the alignment spreads it across a
    # 9-residue internal "gap" that is really un-crystallised regions outside
    # the resolved domain.  FRUTON must refuse to fill in this case.
    ali_path = tmp_path / "test.ali"
    ali_path.write_text(
        ">P1;template\n"
        "structureX:template:FIRST:@:LAST:@::::\n"
        "AA---------AA*\n"
        ">P1;target\n"
        "sequence:target:FIRST:@:LAST:@::::\n"
        "AAQQQQQQQQQAA*\n",
        encoding="utf-8",
    )

    result = analyze_fill_decision(ali_path, "template")

    assert result.should_run_modeller is False
    assert result.overall_classification == "domain_only_spread"
    assert result.alphafold_candidate is False
    assert "Template resolves only part" in (result.skip_reason or "")


def test_find_alignment_fasta_for_filler_prefers_exact_chain_match(
    tmp_path: Path,
) -> None:
    alignment_dir = tmp_path / "alignments"
    alignment_dir.mkdir()

    exact = alignment_dir / "ATOM_chain_A_vs_UniProt.aln.fasta"
    other = alignment_dir / "ATOM_chain_B_vs_UniProt.aln.fasta"
    exact.write_text(">x\nAA\n>y\nAA\n", encoding="utf-8")
    other.write_text(">x\nBB\n>y\nBB\n", encoding="utf-8")

    selected = find_alignment_fasta_for_filler(alignment_dir, "A")
    assert selected == exact


def test_find_alignment_fasta_for_filler_uses_single_fallback_when_only_one_exists(
    tmp_path: Path,
) -> None:
    alignment_dir = tmp_path / "alignments"
    alignment_dir.mkdir()

    fallback = alignment_dir / "ATOM_vs_UniProt.aln.fasta"
    fallback.write_text(">x\nAA\n>y\nAA\n", encoding="utf-8")

    selected = find_alignment_fasta_for_filler(alignment_dir, "A")
    assert selected == fallback


def test_find_alignment_fasta_for_filler_raises_for_multiple_fallback_candidates(
    tmp_path: Path,
) -> None:
    alignment_dir = tmp_path / "alignments"
    alignment_dir.mkdir()

    (alignment_dir / "ATOM_chain_A_vs_UniProt.aln.fasta").write_text(
        ">x\nAA\n>y\nAA\n",
        encoding="utf-8",
    )
    (alignment_dir / "ATOM_chain_B_vs_UniProt.aln.fasta").write_text(
        ">x\nBB\n>y\nBB\n",
        encoding="utf-8",
    )

    with pytest.raises(ValueError):
        find_alignment_fasta_for_filler(alignment_dir, "C")


def test_read_two_sequence_fasta_reads_exactly_two_records(tmp_path: Path) -> None:
    fasta_path = tmp_path / "pair.fasta"
    fasta_path.write_text(
        ">first\nAAAA\n>second\nBBBB\n",
        encoding="utf-8",
    )

    first, second = read_two_sequence_fasta(fasta_path)

    assert first == ("first", "AAAA")
    assert second == ("second", "BBBB")


def test_read_two_sequence_fasta_raises_for_wrong_record_count(tmp_path: Path) -> None:
    fasta_path = tmp_path / "bad.fasta"
    fasta_path.write_text(">only\nAAAA\n", encoding="utf-8")

    with pytest.raises(ValueError):
        read_two_sequence_fasta(fasta_path)


def test_extract_chain_id_from_alignment_header_parses_supported_formats() -> None:
    assert _extract_chain_id_from_alignment_header("PDB|1ABC|ATOM|chain_A") == "A"
    assert _extract_chain_id_from_alignment_header("PDB|1ABC|ATOM|chain-A") == "A"
    assert _extract_chain_id_from_alignment_header("PDB|1ABC|ATOM|chain:A") == "A"
    assert _extract_chain_id_from_alignment_header("PDB|1ABC|ATOM|chain A") == "A"


def test_extract_chain_id_from_alignment_header_returns_none_for_unknown_format() -> (
    None
):
    header = "some_header_without_chain_information"
    assert _extract_chain_id_from_alignment_header(header) is None


def test_split_template_and_target_alignment_records_identifies_pdb_and_uniprot(
    tmp_path: Path,
) -> None:
    fasta_path = tmp_path / "alignment.fasta"
    fasta_path.write_text(
        ">PDB|1ABC|ATOM|chain_A\nAA-A\n>sp|Q12345|Protein|UniProt\nAACA\n",
        encoding="utf-8",
    )

    template_header, template_seq, target_header, target_seq = (
        split_template_and_target_alignment_records(fasta_path)
    )

    assert template_header == "PDB|1ABC|ATOM|chain_A"
    assert template_seq == "AA-A"
    assert target_header == "sp|Q12345|Protein|UniProt"
    assert target_seq == "AACA"


def test_write_chain_specific_template_pdb_keeps_only_requested_chain(
    tmp_path: Path,
) -> None:
    template_pdb = _build_template_pdb_with_two_chains(tmp_path / "template.pdb")
    output_dir = tmp_path / "out"

    written = write_chain_specific_template_pdb(
        template_pdb_path=template_pdb,
        output_dir=output_dir,
        template_id="template_A",
        chain_id="A",
    )

    text = written.read_text(encoding="utf-8")
    assert " A   1" in text
    assert " A   2" in text
    assert " B   1" not in text


def test_write_chain_specific_template_pdb_raises_when_chain_missing(
    tmp_path: Path,
) -> None:
    template_pdb = _build_template_pdb_with_two_chains(tmp_path / "template.pdb")

    with pytest.raises(ValueError):
        write_chain_specific_template_pdb(
            template_pdb_path=template_pdb,
            output_dir=tmp_path / "out",
            template_id="template_Z",
            chain_id="Z",
        )


def test_extract_sequence_from_template_pdb_returns_peptide_sequence(
    tmp_path: Path,
) -> None:
    template_pdb = _build_chain_a_template_pdb(tmp_path / "template_a.pdb")

    sequence = extract_sequence_from_template_pdb(template_pdb)

    assert sequence == "AG"


def test_build_modeller_template_alignment_sequence_rebuilds_using_real_template_sequence() -> (
    None
):
    rebuilt = build_modeller_template_alignment_sequence(
        aligned_template_skeleton="A-GX",
        actual_template_sequence="AGG",
    )

    assert rebuilt == "A-GG"


def test_validate_template_sequence_consistency_writes_debug_files_but_does_not_raise_for_small_difference(
    tmp_path: Path,
) -> None:
    validate_template_sequence_consistency(
        template_alignment_skeleton="AAGG",
        actual_template_sequence="AAGX",
        output_dir=tmp_path,
        template_header="PDB|1ABC|ATOM|chain_A",
        template_pdb_path=tmp_path / "template.pdb",
        max_allowed_length_difference=25,
    )

    assert (tmp_path / "debug_template_alignment_skeleton.txt").exists()
    assert (tmp_path / "debug_template_alignment_ungapped.txt").exists()
    assert (tmp_path / "debug_template_pdb_sequence.txt").exists()


def test_validate_template_sequence_consistency_raises_for_large_length_difference(
    tmp_path: Path,
) -> None:
    with pytest.raises(ValueError):
        validate_template_sequence_consistency(
            template_alignment_skeleton="A" * 50,
            actual_template_sequence="A" * 10,
            output_dir=tmp_path,
            template_header="PDB|1ABC|ATOM|chain_A",
            template_pdb_path=tmp_path / "template.pdb",
            max_allowed_length_difference=5,
        )


def test_write_modeller_alignment_from_existing_alignment_writes_ali_file(
    tmp_path: Path,
) -> None:
    alignment_fasta = tmp_path / "ATOM_chain_A_vs_UniProt.aln.fasta"
    alignment_fasta.write_text(
        ">PDB|1ABC|ATOM|chain_A\nAG-\n>sp|Q12345|Protein|UniProt\nAGQ\n",
        encoding="utf-8",
    )

    template_pdb = _build_chain_a_template_pdb(tmp_path / "template_a.pdb")
    output_dir = tmp_path / "out"

    ali_path = write_modeller_alignment_from_existing_alignment(
        alignment_fasta_path=alignment_fasta,
        template_pdb_path=template_pdb,
        output_dir=output_dir,
        template_id="template_A",
        target_id="1ABC_A",
        chain_id="A",
    )

    text = ali_path.read_text(encoding="utf-8")
    assert ">P1;template_A" in text
    assert "structureX:template_A:FIRST:@:LAST:@::::" in text
    assert ">P1;1ABC_A" in text
    assert "sequence:1ABC_A:FIRST:@:LAST:@::::" in text
    assert "AG-*".replace("*", "")[:2] in text


def test_write_modeller_alignment_from_existing_alignment_raises_for_chain_mismatch(
    tmp_path: Path,
) -> None:
    alignment_fasta = tmp_path / "ATOM_chain_B_vs_UniProt.aln.fasta"
    alignment_fasta.write_text(
        ">PDB|1ABC|ATOM|chain_B\nAG-\n>sp|Q12345|Protein|UniProt\nAGQ\n",
        encoding="utf-8",
    )

    template_pdb = _build_chain_a_template_pdb(tmp_path / "template_a.pdb")

    with pytest.raises(ValueError):
        write_modeller_alignment_from_existing_alignment(
            alignment_fasta_path=alignment_fasta,
            template_pdb_path=template_pdb,
            output_dir=tmp_path / "out",
            template_id="template_A",
            target_id="1ABC_A",
            chain_id="A",
        )


def test_write_skip_logs_creates_stdout_and_stderr_files(tmp_path: Path) -> None:
    stdout_log, stderr_log = _write_skip_logs(tmp_path, "skip message")

    assert stdout_log.exists()
    assert stderr_log.exists()
    assert stdout_log.read_text(encoding="utf-8") == "skip message\n"
    assert stderr_log.read_text(encoding="utf-8") == ""


def test_find_raw_models_prefers_modeller_b_models(tmp_path: Path) -> None:
    output_dir = tmp_path / "models"
    output_dir.mkdir()

    preferred_1 = output_dir / "target.B99990001.pdb"
    preferred_2 = output_dir / "target.B99990002.pdb"
    fallback = output_dir / "target_protein_mod.pdb"

    preferred_1.write_text("ATOM\n", encoding="utf-8")
    preferred_2.write_text("ATOM\n", encoding="utf-8")
    fallback.write_text("ATOM\n", encoding="utf-8")

    result = find_raw_models(output_dir, "target")
    assert result == (preferred_1, preferred_2)


def test_select_best_model_from_scores_returns_lowest_dope_model(
    tmp_path: Path,
) -> None:
    output_dir = tmp_path / "models"
    output_dir.mkdir()

    (output_dir / "model_scores.tsv").write_text(
        "model_name\tdope_score\tga341_score\n"
        "model1.pdb\t-10.0\t1.0\n"
        "model2.pdb\t-25.5\t1.0\n",
        encoding="utf-8",
    )
    (output_dir / "model1.pdb").write_text("ATOM\n", encoding="utf-8")
    (output_dir / "model2.pdb").write_text("ATOM\n", encoding="utf-8")

    best = select_best_model_from_scores(output_dir)
    assert best == output_dir / "model2.pdb"


def test_cleanup_model_pdb_keeps_only_atom_records_and_adds_terminators(
    tmp_path: Path,
) -> None:
    input_model = tmp_path / "input_model.pdb"
    input_model.write_text(
        "ATOM      1  N   ALA A   1      0.000   0.000   0.000  1.00 20.00           N\n"
        "HETATM    2  ZN  ZN  A 100      1.000   1.000   1.000  1.00 20.00          ZN\n"
        "ATOM      3  CA  ALA A   1      1.500   0.000   0.000  1.00 20.00           C\n"
        "END\n",
        encoding="utf-8",
    )

    output_model = tmp_path / "cleaned.pdb"
    cleanup_model_pdb(input_model, output_model)

    lines = output_model.read_text(encoding="utf-8").splitlines()
    assert lines[0].startswith("ATOM")
    assert lines[1].startswith("ATOM")
    assert "HETATM" not in "\n".join(lines)
    assert lines[-2] == "TER"
    assert lines[-1] == "END"
