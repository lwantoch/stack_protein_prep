"""
/home/grheco/repositorios/stack_protein_prep/tests/test_sequence_alignment.py
"""

from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace

import pytest

from stack_protein_preparation.sequence_alignment import (
    AlignmentJob,
    build_alignment_job,
    combine_fasta_files,
    ensure_primary_uniprot_fasta_path,
    fetch_uniprot_fasta,
    get_primary_uniprot_fasta_path,
    resolve_uniprot_accession_from_rcsb,
    run_alignments_for_pdb_directory,
    run_mafft_alignment,
)


def test_build_alignment_job_creates_expected_paths(tmp_path: Path) -> None:
    alignment_directory = tmp_path / "alignments"
    input_paths = [tmp_path / "a.fasta", tmp_path / "b.fasta"]

    job = build_alignment_job(
        alignment_name="SEQRES_vs_UniProt",
        input_fasta_paths=input_paths,
        alignment_directory=alignment_directory,
    )

    assert isinstance(job, AlignmentJob)
    assert job.alignment_name == "SEQRES_vs_UniProt"
    assert job.input_fasta_paths == input_paths
    assert (
        job.combined_input_fasta_path
        == alignment_directory / "SEQRES_vs_UniProt.input.fasta"
    )
    assert (
        job.output_alignment_fasta_path
        == alignment_directory / "SEQRES_vs_UniProt.aln.fasta"
    )


def test_get_primary_uniprot_fasta_path_returns_none_when_missing(
    tmp_path: Path,
) -> None:
    assert get_primary_uniprot_fasta_path(tmp_path) is None


def test_get_primary_uniprot_fasta_path_returns_first_sorted_match(
    tmp_path: Path,
) -> None:
    second = tmp_path / "UniProt_Q99999.fasta"
    first = tmp_path / "UniProt_A11111.fasta"
    second.write_text(">b\nBBBB\n", encoding="utf-8")
    first.write_text(">a\nAAAA\n", encoding="utf-8")

    result = get_primary_uniprot_fasta_path(tmp_path)

    assert result == first


def test_combine_fasta_files_raises_for_empty_input_list(tmp_path: Path) -> None:
    output_path = tmp_path / "combined.fasta"

    with pytest.raises(ValueError, match="empty input_fasta_paths list"):
        combine_fasta_files([], output_path)


def test_combine_fasta_files_raises_for_missing_input_file(tmp_path: Path) -> None:
    missing = tmp_path / "missing.fasta"
    output_path = tmp_path / "combined.fasta"

    with pytest.raises(FileNotFoundError, match="Input FASTA file not found"):
        combine_fasta_files([missing], output_path)


def test_combine_fasta_files_raises_for_empty_input_file(tmp_path: Path) -> None:
    empty_file = tmp_path / "empty.fasta"
    empty_file.write_text("", encoding="utf-8")
    output_path = tmp_path / "combined.fasta"

    with pytest.raises(ValueError, match="Input FASTA file is empty"):
        combine_fasta_files([empty_file], output_path)


def test_combine_fasta_files_writes_combined_output(tmp_path: Path) -> None:
    fasta_1 = tmp_path / "one.fasta"
    fasta_2 = tmp_path / "two.fasta"
    output_path = tmp_path / "nested" / "combined.fasta"

    fasta_1.write_text(">seq1\nAAAA\n", encoding="utf-8")
    fasta_2.write_text(">seq2\nBBBB\n", encoding="utf-8")

    combine_fasta_files([fasta_1, fasta_2], output_path)

    assert output_path.exists()
    assert output_path.read_text(encoding="utf-8") == ">seq1\nAAAA\n>seq2\nBBBB\n"


def test_fetch_uniprot_fasta_writes_downloaded_fasta(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    def fake_http_get_text(url: str) -> str:
        assert "P12345.fasta" in url
        return ">sp|P12345|TEST\nMSEQUENCE\n"

    monkeypatch.setattr(
        "stack_protein_preparation.sequence_alignment._http_get_text",
        fake_http_get_text,
    )

    result = fetch_uniprot_fasta("P12345", tmp_path)

    assert result == tmp_path / "UniProt_P12345.fasta"
    assert result.read_text(encoding="utf-8") == ">sp|P12345|TEST\nMSEQUENCE\n"


def test_fetch_uniprot_fasta_raises_for_empty_response(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.setattr(
        "stack_protein_preparation.sequence_alignment._http_get_text",
        lambda url: "   ",
    )

    with pytest.raises(ValueError, match="Received empty FASTA response"):
        fetch_uniprot_fasta("P12345", tmp_path)


def test_fetch_uniprot_fasta_raises_for_non_fasta_response(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.setattr(
        "stack_protein_preparation.sequence_alignment._http_get_text",
        lambda url: "not a fasta file",
    )

    with pytest.raises(ValueError, match="does not look like FASTA"):
        fetch_uniprot_fasta("P12345", tmp_path)


def test_resolve_uniprot_accession_from_rcsb_returns_first_accession(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    def fake_http_get_json(url: str) -> dict:
        if url.endswith("/entry/1ABC"):
            return {
                "rcsb_entry_container_identifiers": {"polymer_entity_ids": ["1", "2"]}
            }

        if url.endswith("/polymer_entity/1ABC/1"):
            return {
                "rcsb_polymer_entity_container_identifiers": {"uniprot_ids": ["Q12345"]}
            }

        raise AssertionError(f"Unexpected URL: {url}")

    monkeypatch.setattr(
        "stack_protein_preparation.sequence_alignment._http_get_json",
        fake_http_get_json,
    )

    result = resolve_uniprot_accession_from_rcsb("1abc")

    assert result == "Q12345"


def test_resolve_uniprot_accession_from_rcsb_returns_none_when_no_mapping(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    def fake_http_get_json(url: str) -> dict:
        if url.endswith("/entry/1ABC"):
            return {
                "rcsb_entry_container_identifiers": {"polymer_entity_ids": ["1", "2"]}
            }

        if url.endswith("/polymer_entity/1ABC/1"):
            return {}

        if url.endswith("/polymer_entity/1ABC/2"):
            return {}

        raise AssertionError(f"Unexpected URL: {url}")

    monkeypatch.setattr(
        "stack_protein_preparation.sequence_alignment._http_get_json",
        fake_http_get_json,
    )

    result = resolve_uniprot_accession_from_rcsb("1abc")

    assert result is None


def test_ensure_primary_uniprot_fasta_path_returns_existing_local_file(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    local_file = tmp_path / "UniProt_P12345.fasta"
    local_file.write_text(">x\nAAAA\n", encoding="utf-8")

    def fail(*args, **kwargs):
        raise AssertionError(
            "Remote code path should not be called when local file exists."
        )

    monkeypatch.setattr(
        "stack_protein_preparation.sequence_alignment.resolve_uniprot_accession_from_rcsb",
        fail,
    )
    monkeypatch.setattr(
        "stack_protein_preparation.sequence_alignment.fetch_uniprot_fasta",
        fail,
    )

    result = ensure_primary_uniprot_fasta_path(tmp_path, "1ABC")

    assert result == local_file


def test_ensure_primary_uniprot_fasta_path_fetches_remote_when_missing(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    expected = tmp_path / "UniProt_Q99999.fasta"

    monkeypatch.setattr(
        "stack_protein_preparation.sequence_alignment.resolve_uniprot_accession_from_rcsb",
        lambda pdb_id: "Q99999",
    )
    monkeypatch.setattr(
        "stack_protein_preparation.sequence_alignment.fetch_uniprot_fasta",
        lambda uniprot_accession, fasta_directory: expected,
    )

    result = ensure_primary_uniprot_fasta_path(tmp_path, "1ABC")

    assert result == expected


def test_ensure_primary_uniprot_fasta_path_returns_none_when_resolution_fails(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    monkeypatch.setattr(
        "stack_protein_preparation.sequence_alignment.resolve_uniprot_accession_from_rcsb",
        lambda pdb_id: None,
    )

    result = ensure_primary_uniprot_fasta_path(tmp_path, "1ABC")

    assert result is None


def test_ensure_primary_uniprot_fasta_path_returns_none_when_resolution_raises(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    def boom(pdb_id: str) -> str:
        raise RuntimeError("RCSB down")

    monkeypatch.setattr(
        "stack_protein_preparation.sequence_alignment.resolve_uniprot_accession_from_rcsb",
        boom,
    )

    result = ensure_primary_uniprot_fasta_path(tmp_path, "1ABC")

    assert result is None


def test_ensure_primary_uniprot_fasta_path_returns_none_when_fetch_raises(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    monkeypatch.setattr(
        "stack_protein_preparation.sequence_alignment.resolve_uniprot_accession_from_rcsb",
        lambda pdb_id: "Q99999",
    )

    def boom(uniprot_accession: str, fasta_directory: Path) -> Path:
        raise RuntimeError("UniProt down")

    monkeypatch.setattr(
        "stack_protein_preparation.sequence_alignment.fetch_uniprot_fasta",
        boom,
    )

    result = ensure_primary_uniprot_fasta_path(tmp_path, "1ABC")

    assert result is None


def test_run_mafft_alignment_raises_when_mafft_missing(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    input_path = tmp_path / "input.fasta"
    output_path = tmp_path / "output.aln.fasta"
    input_path.write_text(">a\nAAAA\n", encoding="utf-8")

    monkeypatch.setattr(
        "stack_protein_preparation.sequence_alignment.shutil.which",
        lambda name: None,
    )

    with pytest.raises(FileNotFoundError, match="MAFFT executable not found"):
        run_mafft_alignment(input_path, output_path)


def test_run_mafft_alignment_writes_output_when_mafft_succeeds(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    input_path = tmp_path / "input.fasta"
    output_path = tmp_path / "out" / "output.aln.fasta"
    input_path.write_text(">a\nAAAA\n>b\nBBBB\n", encoding="utf-8")

    monkeypatch.setattr(
        "stack_protein_preparation.sequence_alignment.shutil.which",
        lambda name: "/usr/bin/mafft",
    )

    def fake_run(*args, **kwargs):
        return SimpleNamespace(
            returncode=0,
            stdout=">a\nAAAA\n>b\nBBBB\n",
            stderr="",
        )

    monkeypatch.setattr(
        "stack_protein_preparation.sequence_alignment.subprocess.run",
        fake_run,
    )

    run_mafft_alignment(input_path, output_path)

    assert output_path.exists()
    assert output_path.read_text(encoding="utf-8") == ">a\nAAAA\n>b\nBBBB\n"


def test_run_mafft_alignment_raises_when_mafft_fails(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    input_path = tmp_path / "input.fasta"
    output_path = tmp_path / "output.aln.fasta"
    input_path.write_text(">a\nAAAA\n", encoding="utf-8")

    monkeypatch.setattr(
        "stack_protein_preparation.sequence_alignment.shutil.which",
        lambda name: "/usr/bin/mafft",
    )

    def fake_run(*args, **kwargs):
        return SimpleNamespace(
            returncode=1,
            stdout="",
            stderr="something went wrong",
        )

    monkeypatch.setattr(
        "stack_protein_preparation.sequence_alignment.subprocess.run",
        fake_run,
    )

    with pytest.raises(RuntimeError, match="MAFFT failed"):
        run_mafft_alignment(input_path, output_path)


def test_run_alignments_for_pdb_directory_creates_seqres_and_atom_jobs(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    pdb_directory = tmp_path / "1abc"
    fasta_directory = pdb_directory / "fasta"
    fasta_directory.mkdir(parents=True)

    seqres_path = fasta_directory / "PDB-1ABC-SEQRES.fasta"
    atom_path = fasta_directory / "PDB-1ABC-ATOM.fasta"
    uniprot_path = fasta_directory / "UniProt_Q12345.fasta"

    seqres_path.write_text(">seqres\nAAAA\n", encoding="utf-8")
    atom_path.write_text(">atom\nBBBB\n", encoding="utf-8")
    uniprot_path.write_text(">uni\nCCCC\n", encoding="utf-8")

    monkeypatch.setattr(
        "stack_protein_preparation.sequence_alignment.ensure_primary_uniprot_fasta_path",
        lambda fasta_directory, pdb_id: uniprot_path,
    )

    recorded_jobs: list[AlignmentJob] = []

    def fake_run_alignment_job(job: AlignmentJob) -> None:
        recorded_jobs.append(job)
        job.output_alignment_fasta_path.parent.mkdir(parents=True, exist_ok=True)
        job.output_alignment_fasta_path.write_text(">aln\nAAAA\n", encoding="utf-8")

    image_calls: list[tuple[Path, Path]] = []

    def fake_alignment_to_image(alignment_path: Path, output_png_path: Path) -> None:
        image_calls.append((alignment_path, output_png_path))

    monkeypatch.setattr(
        "stack_protein_preparation.sequence_alignment.run_alignment_job",
        fake_run_alignment_job,
    )
    monkeypatch.setattr(
        "stack_protein_preparation.sequence_alignment.alignment_to_image",
        fake_alignment_to_image,
    )

    run_alignments_for_pdb_directory(pdb_directory)

    assert len(recorded_jobs) == 2
    assert recorded_jobs[0].alignment_name == "SEQRES_vs_UniProt"
    assert recorded_jobs[1].alignment_name == "ATOM_vs_UniProt"
    assert recorded_jobs[0].input_fasta_paths == [seqres_path, uniprot_path]
    assert recorded_jobs[1].input_fasta_paths == [atom_path, uniprot_path]
    assert len(image_calls) == 2
    assert image_calls[0][0].name == "SEQRES_vs_UniProt.aln.fasta"
    assert image_calls[0][1].name == "SEQRES_vs_UniProt.aln.png"
    assert image_calls[1][0].name == "ATOM_vs_UniProt.aln.fasta"
    assert image_calls[1][1].name == "ATOM_vs_UniProt.aln.png"


def test_run_alignments_for_pdb_directory_skips_when_no_uniprot_available(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    pdb_directory = tmp_path / "1abc"
    fasta_directory = pdb_directory / "fasta"
    fasta_directory.mkdir(parents=True)

    (fasta_directory / "PDB-1ABC-SEQRES.fasta").write_text(
        ">seqres\nAAAA\n", encoding="utf-8"
    )

    monkeypatch.setattr(
        "stack_protein_preparation.sequence_alignment.ensure_primary_uniprot_fasta_path",
        lambda fasta_directory, pdb_id: None,
    )

    called = {"run_alignment_job": False}

    def fake_run_alignment_job(job: AlignmentJob) -> None:
        called["run_alignment_job"] = True

    monkeypatch.setattr(
        "stack_protein_preparation.sequence_alignment.run_alignment_job",
        fake_run_alignment_job,
    )

    run_alignments_for_pdb_directory(pdb_directory)

    assert called["run_alignment_job"] is False


def test_run_alignments_for_pdb_directory_returns_when_fasta_directory_missing(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    pdb_directory = tmp_path / "1abc"

    called = {"ensure": False}

    def fake_ensure(*args, **kwargs):
        called["ensure"] = True
        return None

    monkeypatch.setattr(
        "stack_protein_preparation.sequence_alignment.ensure_primary_uniprot_fasta_path",
        fake_ensure,
    )

    run_alignments_for_pdb_directory(pdb_directory)

    assert called["ensure"] is False
