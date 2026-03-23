"""
/home/grheco/repositorios/stack_protein_prep/src/stack_protein_preparation/sequence_alignment.py

Run pairwise sequence alignments with MAFFT for local protein FASTA files.

Responsibilities
----------------
- collect local FASTA files for one PDB directory
- ensure a UniProt FASTA exists locally
- split local multi-chain PDB FASTA files into chain-specific pairwise jobs
- combine selected FASTA records into temporary two-sequence input files
- run MAFFT via subprocess
- write alignment output files to a dedicated alignment directory
- optionally render alignment PNG files
- create residue-position mapping files from PDB FASTA positions to UniProt positions

Expected directory layout
-------------------------
Input:
    data/proteins/<PDB_ID>/fasta/
    ├── PDB-<PDB_ID>-SEQRES.fasta
    ├── PDB-<PDB_ID>-ATOM.fasta
    ├── UniProt_<UNIPROT_ID>.fasta
    └── ...

Output:
    data/proteins/<PDB_ID>/fasta/alignments/
    ├── SEQRES_chain_A_vs_UniProt.input.fasta
    ├── SEQRES_chain_A_vs_UniProt.aln.fasta
    ├── SEQRES_chain_A_vs_UniProt.mapping.tsv
    ├── SEQRES_chain_A_vs_UniProt.png
    ├── ATOM_chain_A_vs_UniProt.input.fasta
    ├── ATOM_chain_A_vs_UniProt.aln.fasta
    ├── ATOM_chain_A_vs_UniProt.mapping.tsv
    ├── ATOM_chain_A_vs_UniProt.png
    └── ...

Design choices
--------------
- This module performs pairwise chain-specific comparisons:
    1. each SEQRES chain vs UniProt
    2. each ATOM chain vs UniProt
- If multiple UniProt FASTA files exist, the first one in sorted order is used.
- If no UniProt FASTA file exists locally, this module tries to:
    1. resolve a UniProt accession from the PDB entry via RCSB
    2. download the UniProt FASTA
- Temporary combined input FASTA files are kept on disk intentionally so that the
  exact MAFFT input remains inspectable for debugging and reproducibility.
- Mapping files are position-based. They map residue positions in the local PDB FASTA
  sequence to residue positions in the UniProt sequence based on the alignment.
- This does NOT automatically reconstruct original author residue IDs from the PDB file.
  It maps sequence indices unless another upstream step provides explicit residue numbering.
- Alignment rendering is imported lazily so that non-visual tests do not require
  matplotlib at import time.

What this module does NOT do
----------------------------
- It does not analyze the biological meaning of the alignment.
- It does not calculate percent identity yet.
- It does not decide which UniProt entry is biologically correct if multiple
  plausible mappings exist.
- It does not parse structural residue numbering from PDB/mmCIF files.

Testing note
------------
- MAFFT execution should usually be mocked in unit tests.
- HTTP calls should usually be mocked in unit tests.
- alignment PNG rendering should remain optional in tests and at runtime.
"""

from __future__ import annotations

import json
import re
import shutil
import subprocess
import urllib.error
import urllib.request
from dataclasses import dataclass
from pathlib import Path
from typing import Any

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

ALIGNMENTS_SUBDIRECTORY_NAME = "alignments"
UNIPROT_FASTA_GLOB_PATTERN = "UniProt_*.fasta"

RCSB_ENTRY_URL_TEMPLATE = "https://data.rcsb.org/rest/v1/core/entry/{pdb_id}"
RCSB_POLYMER_ENTITY_URL_TEMPLATE = (
    "https://data.rcsb.org/rest/v1/core/polymer_entity/{pdb_id}/{entity_id}"
)
UNIPROT_FASTA_URL_TEMPLATE = "https://rest.uniprot.org/uniprotkb/{accession}.fasta"

HTTP_TIMEOUT_SECONDS = 30


# ---------------------------------------------------------------------------
# Small data containers
# ---------------------------------------------------------------------------


@dataclass
class AlignmentJob:
    """
    Definition of one MAFFT alignment job.

    Attributes
    ----------
    alignment_name
        Short descriptive job label, used for output filenames.
    input_fasta_paths
        FASTA files that will be concatenated and passed to MAFFT.
    combined_input_fasta_path
        Combined temporary FASTA file written before MAFFT execution.
    output_alignment_fasta_path
        MAFFT output alignment path in FASTA format.
    """

    alignment_name: str
    input_fasta_paths: list[Path]
    combined_input_fasta_path: Path
    output_alignment_fasta_path: Path


# ---------------------------------------------------------------------------
# Optional visualization import
# ---------------------------------------------------------------------------


def render_alignment_image(
    alignment_fasta_path: Path,
    output_png_path: Path,
) -> None:
    """
    Render one alignment PNG using the optional visualization module.

    This import is intentionally local so that importing sequence_alignment.py
    does not require matplotlib unless image rendering is actually requested.
    """
    from stack_protein_preparation.alignment_visualization import alignment_to_image

    alignment_to_image(alignment_fasta_path, output_png_path)


# ---------------------------------------------------------------------------
# Public high-level function
# ---------------------------------------------------------------------------


def run_alignments_for_pdb_directory(
    pdb_directory: Path,
    render_images: bool = True,
) -> None:
    """
    Run standard MAFFT alignments for one local PDB directory.

    Current standard jobs
    ---------------------
    - each chain in PDB-<PDB_ID>-SEQRES.fasta vs first UniProt FASTA
    - each chain in PDB-<PDB_ID>-ATOM.fasta   vs first UniProt FASTA

    Parameters
    ----------
    pdb_directory
        Directory such as:
            data/proteins/1W4R
    render_images
        If True, render alignment PNG files after successful MAFFT execution.
        If False, skip visualization entirely.
    """
    pdb_id = pdb_directory.name.upper()
    fasta_directory = pdb_directory / "fasta"
    alignment_directory = fasta_directory / ALIGNMENTS_SUBDIRECTORY_NAME

    print(f"[INFO] Preparing alignment jobs for PDB {pdb_id}")

    if not fasta_directory.exists():
        print(f"[WARNING] FASTA directory not found for {pdb_id}: {fasta_directory}")
        return

    alignment_directory.mkdir(parents=True, exist_ok=True)

    seqres_fasta_path = fasta_directory / f"PDB-{pdb_id}-SEQRES.fasta"
    atom_fasta_path = fasta_directory / f"PDB-{pdb_id}-ATOM.fasta"

    uniprot_fasta_path = ensure_primary_uniprot_fasta_path(
        fasta_directory=fasta_directory,
        pdb_id=pdb_id,
    )

    if uniprot_fasta_path is None:
        print(
            f"[WARNING] No UniProt FASTA available for {pdb_id}. "
            "Skipping alignment jobs."
        )
        return

    print(f"[INFO] Using UniProt FASTA: {uniprot_fasta_path.name}")

    alignment_job_list: list[AlignmentJob] = []

    if seqres_fasta_path.exists():
        alignment_job_list.extend(
            build_chain_specific_alignment_jobs(
                pdb_fasta_path=seqres_fasta_path,
                uniprot_fasta_path=uniprot_fasta_path,
                alignment_prefix="SEQRES",
                alignment_directory=alignment_directory,
            )
        )
    else:
        print(f"[WARNING] Missing SEQRES FASTA for {pdb_id}: {seqres_fasta_path}")

    if atom_fasta_path.exists():
        alignment_job_list.extend(
            build_chain_specific_alignment_jobs(
                pdb_fasta_path=atom_fasta_path,
                uniprot_fasta_path=uniprot_fasta_path,
                alignment_prefix="ATOM",
                alignment_directory=alignment_directory,
            )
        )
    else:
        print(f"[WARNING] Missing ATOM FASTA for {pdb_id}: {atom_fasta_path}")

    if not alignment_job_list:
        print(f"[WARNING] No alignment jobs could be created for {pdb_id}")
        return

    for alignment_job in alignment_job_list:
        run_alignment_job(alignment_job)

        alignment_path = alignment_job.output_alignment_fasta_path
        output_png_path = alignment_path.with_suffix(".png")
        mapping_path = alignment_path.with_suffix(".mapping.tsv")

        if alignment_path.exists() and alignment_path.stat().st_size > 0:
            write_alignment_mapping_file(
                alignment_fasta_path=alignment_path,
                output_mapping_tsv_path=mapping_path,
            )

            if render_images:
                render_alignment_image(
                    alignment_fasta_path=alignment_path,
                    output_png_path=output_png_path,
                )
        else:
            print(f"[WARNING] Alignment file missing or empty: {alignment_path}")

    print(f"[INFO] Finished all alignment jobs for {pdb_id}")


# ---------------------------------------------------------------------------
# UniProt FASTA handling
# ---------------------------------------------------------------------------


def ensure_primary_uniprot_fasta_path(
    fasta_directory: Path,
    pdb_id: str,
) -> Path | None:
    """
    Return a local UniProt FASTA path.

    Behavior
    --------
    1. If a local UniProt FASTA already exists, return it.
    2. Otherwise try to resolve a UniProt accession from RCSB for this PDB ID.
    3. If resolution succeeds, fetch the FASTA from UniProt and save it locally.
    4. Return the local FASTA path, or None if all attempts fail.
    """
    local_uniprot_fasta_path = get_primary_uniprot_fasta_path(fasta_directory)
    if local_uniprot_fasta_path is not None:
        return local_uniprot_fasta_path

    print(f"[INFO] No local UniProt FASTA found for {pdb_id}. Trying remote fetch.")

    try:
        uniprot_accession = resolve_uniprot_accession_from_rcsb(pdb_id)
    except Exception as exc:
        print(f"[WARNING] Failed to resolve UniProt accession for PDB {pdb_id}: {exc}")
        return None

    if uniprot_accession is None:
        print(f"[WARNING] Could not resolve UniProt accession for PDB {pdb_id}")
        return None

    print(f"[INFO] Resolved UniProt accession for {pdb_id}: {uniprot_accession}")

    try:
        return fetch_uniprot_fasta(
            uniprot_accession=uniprot_accession,
            fasta_directory=fasta_directory,
        )
    except Exception as exc:
        print(
            f"[WARNING] Failed to fetch UniProt FASTA for {pdb_id} "
            f"({uniprot_accession}): {exc}"
        )
        return None


def resolve_uniprot_accession_from_rcsb(pdb_id: str) -> str | None:
    """
    Resolve a UniProt accession for a PDB entry using RCSB Data API.

    Strategy
    --------
    - query the entry endpoint to get polymer entity IDs
    - query each polymer entity endpoint
    - look for UniProt accession(s) in common identifier locations

    Returns
    -------
    str | None
        First resolved UniProt accession, or None if nothing could be found.
    """
    pdb_id = pdb_id.upper().strip()

    entry_url = RCSB_ENTRY_URL_TEMPLATE.format(pdb_id=pdb_id)
    entry_data = _http_get_json(entry_url)

    polymer_entity_ids = (
        entry_data.get("rcsb_entry_container_identifiers", {}).get(
            "polymer_entity_ids", []
        )
        or []
    )

    for entity_id in polymer_entity_ids:
        entity_url = RCSB_POLYMER_ENTITY_URL_TEMPLATE.format(
            pdb_id=pdb_id,
            entity_id=entity_id,
        )
        entity_data = _http_get_json(entity_url)

        accession = _extract_uniprot_accession_from_polymer_entity(entity_data)
        if accession:
            return accession

    return None


def _extract_uniprot_accession_from_polymer_entity(
    entity_data: dict[str, Any],
) -> str | None:
    """
    Try several known / plausible locations for UniProt accessions in RCSB
    polymer-entity JSON.
    """
    identifier_block = entity_data.get("rcsb_polymer_entity_container_identifiers", {})

    uniprot_ids = identifier_block.get("uniprot_ids", [])
    if isinstance(uniprot_ids, list) and uniprot_ids:
        first_id = str(uniprot_ids[0]).strip()
        if first_id:
            return first_id

    candidate_strings: list[str] = []
    _collect_uniprot_like_values(entity_data, candidate_strings)

    for candidate in candidate_strings:
        cleaned = candidate.strip()
        if cleaned:
            return cleaned

    return None


def _collect_uniprot_like_values(obj: Any, output: list[str]) -> None:
    """
    Recursive fallback collector for UniProt-like fields.
    Kept intentionally conservative.
    """
    if isinstance(obj, dict):
        for key, value in obj.items():
            lower_key = str(key).lower()
            if lower_key in {
                "uniprot_id",
                "uniprot_ids",
                "accession",
                "accessions",
                "primary_accession",
            }:
                if isinstance(value, str):
                    output.append(value)
                elif isinstance(value, list):
                    for item in value:
                        if isinstance(item, str):
                            output.append(item)
                elif isinstance(value, dict):
                    _collect_uniprot_like_values(value, output)
            else:
                _collect_uniprot_like_values(value, output)

    elif isinstance(obj, list):
        for item in obj:
            _collect_uniprot_like_values(item, output)


def fetch_uniprot_fasta(
    uniprot_accession: str,
    fasta_directory: Path,
) -> Path:
    """
    Download one UniProt FASTA file and save it locally as:
        UniProt_<ACCESSION>.fasta
    """
    uniprot_accession = uniprot_accession.strip()
    fasta_directory.mkdir(parents=True, exist_ok=True)

    output_fasta_path = fasta_directory / f"UniProt_{uniprot_accession}.fasta"
    fasta_url = UNIPROT_FASTA_URL_TEMPLATE.format(accession=uniprot_accession)

    fasta_text = _http_get_text(fasta_url)

    if not fasta_text.strip():
        raise ValueError(
            f"Received empty FASTA response for UniProt accession {uniprot_accession}"
        )

    if not fasta_text.lstrip().startswith(">"):
        raise ValueError(
            f"Response for UniProt accession {uniprot_accession} "
            "does not look like FASTA."
        )

    output_fasta_path.write_text(fasta_text, encoding="utf-8")
    print(f"[INFO] Downloaded UniProt FASTA: {output_fasta_path}")

    return output_fasta_path


def _http_get_json(url: str) -> dict[str, Any]:
    """
    GET one JSON document.
    """
    request = urllib.request.Request(
        url,
        headers={"User-Agent": "stack-protein-preparation/1.0"},
    )

    try:
        with urllib.request.urlopen(request, timeout=HTTP_TIMEOUT_SECONDS) as response:
            raw_text = response.read().decode("utf-8")
    except urllib.error.HTTPError as exc:
        raise RuntimeError(f"HTTP error for {url}: {exc.code}") from exc
    except urllib.error.URLError as exc:
        raise RuntimeError(f"URL error for {url}: {exc.reason}") from exc

    try:
        data = json.loads(raw_text)
    except json.JSONDecodeError as exc:
        raise RuntimeError(f"Invalid JSON returned by {url}") from exc

    if not isinstance(data, dict):
        raise RuntimeError(f"Unexpected non-dict JSON returned by {url}")

    return data


def _http_get_text(url: str) -> str:
    """
    GET one text response.
    """
    request = urllib.request.Request(
        url,
        headers={"User-Agent": "stack-protein-preparation/1.0"},
    )

    try:
        with urllib.request.urlopen(request, timeout=HTTP_TIMEOUT_SECONDS) as response:
            return response.read().decode("utf-8")
    except urllib.error.HTTPError as exc:
        raise RuntimeError(f"HTTP error for {url}: {exc.code}") from exc
    except urllib.error.URLError as exc:
        raise RuntimeError(f"URL error for {url}: {exc.reason}") from exc


# ---------------------------------------------------------------------------
# Alignment job construction
# ---------------------------------------------------------------------------


def build_alignment_job(
    alignment_name: str,
    input_fasta_paths: list[Path],
    alignment_directory: Path,
) -> AlignmentJob:
    """
    Create one AlignmentJob object with standardized filenames.
    """
    combined_input_fasta_path = alignment_directory / f"{alignment_name}.input.fasta"
    output_alignment_fasta_path = alignment_directory / f"{alignment_name}.aln.fasta"

    return AlignmentJob(
        alignment_name=alignment_name,
        input_fasta_paths=input_fasta_paths,
        combined_input_fasta_path=combined_input_fasta_path,
        output_alignment_fasta_path=output_alignment_fasta_path,
    )


def build_chain_specific_alignment_jobs(
    pdb_fasta_path: Path,
    uniprot_fasta_path: Path,
    alignment_prefix: str,
    alignment_directory: Path,
) -> list[AlignmentJob]:
    """
    Build one pairwise alignment job per chain-specific PDB FASTA record.

    Parameters
    ----------
    pdb_fasta_path
        Multi-entry local PDB-derived FASTA file.
    uniprot_fasta_path
        UniProt FASTA file expected to contain exactly one record.
    alignment_prefix
        Prefix such as "SEQRES" or "ATOM".
    alignment_directory
        Output directory for temporary and final alignment files.

    Returns
    -------
    list[AlignmentJob]
        One job per local PDB FASTA record.
    """
    pdb_records = read_fasta_records(pdb_fasta_path)
    if not pdb_records:
        print(f"[WARNING] No FASTA records found in {pdb_fasta_path}")
        return []

    uniprot_records = read_fasta_records(uniprot_fasta_path)
    if len(uniprot_records) != 1:
        raise ValueError(
            f"Expected exactly 1 UniProt FASTA record in {uniprot_fasta_path}, "
            f"but found {len(uniprot_records)}."
        )

    uniprot_header, uniprot_sequence = uniprot_records[0]

    alignment_job_list: list[AlignmentJob] = []

    for pdb_header, pdb_sequence in pdb_records:
        chain_label = extract_chain_label_from_header(pdb_header)
        alignment_name = f"{alignment_prefix}_{chain_label}_vs_UniProt"

        combined_input_fasta_path = (
            alignment_directory / f"{alignment_name}.input.fasta"
        )
        output_alignment_fasta_path = (
            alignment_directory / f"{alignment_name}.aln.fasta"
        )

        write_two_record_fasta(
            first_header=pdb_header,
            first_sequence=pdb_sequence,
            second_header=uniprot_header,
            second_sequence=uniprot_sequence,
            output_fasta_path=combined_input_fasta_path,
        )

        alignment_job_list.append(
            AlignmentJob(
                alignment_name=alignment_name,
                input_fasta_paths=[],
                combined_input_fasta_path=combined_input_fasta_path,
                output_alignment_fasta_path=output_alignment_fasta_path,
            )
        )

    return alignment_job_list


def extract_chain_label_from_header(fasta_header: str) -> str:
    """
    Extract a filesystem-safe chain label from a FASTA header.

    Examples
    --------
    PDB|2W8Y|ATOM|chain_A -> chain_A
    PDB|2W8Y|SEQRES|chain_B -> chain_B
    """
    header_parts = fasta_header.split("|")
    if header_parts:
        last_part = header_parts[-1].strip()
        if last_part:
            safe_last_part = re.sub(r"[^A-Za-z0-9_.-]+", "_", last_part)
            if safe_last_part:
                return safe_last_part

    safe_header = re.sub(r"[^A-Za-z0-9_.-]+", "_", fasta_header).strip("_")
    return safe_header or "chain_unknown"


def write_two_record_fasta(
    first_header: str,
    first_sequence: str,
    second_header: str,
    second_sequence: str,
    output_fasta_path: Path,
) -> None:
    """
    Write exactly two FASTA records to one file.
    """
    output_fasta_path.parent.mkdir(parents=True, exist_ok=True)

    output_text = format_fasta_record(
        first_header, first_sequence
    ) + format_fasta_record(second_header, second_sequence)
    output_fasta_path.write_text(output_text, encoding="utf-8")


def format_fasta_record(header: str, sequence: str, line_width: int = 80) -> str:
    """
    Format one FASTA record with fixed sequence line width.
    """
    wrapped_sequence_lines = [
        sequence[i : i + line_width] for i in range(0, len(sequence), line_width)
    ]
    return f">{header}\n" + "\n".join(wrapped_sequence_lines) + "\n"


# ---------------------------------------------------------------------------
# Alignment job execution
# ---------------------------------------------------------------------------


def run_alignment_job(alignment_job: AlignmentJob) -> None:
    """
    Execute one alignment job:
    1. use an already prepared two-record FASTA input
    2. run MAFFT
    """
    print(f"[INFO] Running alignment job: {alignment_job.alignment_name}")

    run_mafft_alignment(
        input_fasta_path=alignment_job.combined_input_fasta_path,
        output_alignment_path=alignment_job.output_alignment_fasta_path,
    )

    print(
        f"[INFO] Alignment job finished: {alignment_job.alignment_name} -> "
        f"{alignment_job.output_alignment_fasta_path}"
    )


def run_mafft_alignment(
    input_fasta_path: Path,
    output_alignment_path: Path,
) -> None:
    """
    Run MAFFT on a combined FASTA input file and write aligned FASTA output.
    """
    print(f"[INFO] Running MAFFT on {input_fasta_path}")

    if shutil.which("mafft") is None:
        raise FileNotFoundError(
            "MAFFT executable not found in PATH. Install mafft in the pixi environment."
        )

    completed_process = subprocess.run(
        ["mafft", "--auto", str(input_fasta_path)],
        capture_output=True,
        text=True,
        check=False,
    )

    if completed_process.returncode != 0:
        raise RuntimeError(
            "MAFFT failed.\n"
            f"Input FASTA: {input_fasta_path}\n"
            f"Return code: {completed_process.returncode}\n"
            f"STDERR:\n{completed_process.stderr}"
        )

    output_alignment_path.parent.mkdir(parents=True, exist_ok=True)
    output_alignment_path.write_text(completed_process.stdout, encoding="utf-8")

    if completed_process.stderr.strip():
        print("[DEBUG] MAFFT stderr output:")
        print(completed_process.stderr.strip())

    print(f"[INFO] Alignment written to {output_alignment_path}")


# ---------------------------------------------------------------------------
# FASTA collection helpers
# ---------------------------------------------------------------------------


def get_primary_uniprot_fasta_path(fasta_directory: Path) -> Path | None:
    """
    Select the primary UniProt FASTA file from a FASTA directory.

    Current policy
    --------------
    - collect all files matching 'UniProt_*.fasta'
    - sort them lexicographically
    - return the first file
    """
    uniprot_fasta_path_list = sorted(fasta_directory.glob(UNIPROT_FASTA_GLOB_PATTERN))

    if not uniprot_fasta_path_list:
        return None

    return uniprot_fasta_path_list[0]


def combine_fasta_files(
    input_fasta_paths: list[Path],
    output_fasta_path: Path,
) -> None:
    """
    Combine multiple FASTA files into one multi-entry FASTA file.
    Kept for compatibility, although chain-specific jobs currently write
    prepared two-record inputs directly.
    """
    if not input_fasta_paths:
        raise ValueError(
            "combine_fasta_files() received an empty input_fasta_paths list."
        )

    output_fasta_path.parent.mkdir(parents=True, exist_ok=True)

    combined_fasta_chunks: list[str] = []

    for input_fasta_path in input_fasta_paths:
        if not input_fasta_path.exists():
            raise FileNotFoundError(f"Input FASTA file not found: {input_fasta_path}")

        fasta_text = input_fasta_path.read_text(encoding="utf-8").strip()

        if not fasta_text:
            raise ValueError(f"Input FASTA file is empty: {input_fasta_path}")

        combined_fasta_chunks.append(fasta_text)

    combined_fasta_text = "\n".join(combined_fasta_chunks) + "\n"
    output_fasta_path.write_text(combined_fasta_text, encoding="utf-8")

    print(f"[INFO] Combined FASTA written to {output_fasta_path}")


# ---------------------------------------------------------------------------
# Alignment parsing and mapping creation
# ---------------------------------------------------------------------------


def write_alignment_mapping_file(
    alignment_fasta_path: Path,
    output_mapping_tsv_path: Path,
) -> None:
    """
    Create a residue-position mapping TSV from a pairwise alignment FASTA file.

    Assumptions
    -----------
    - alignment contains exactly two sequences:
        1. local PDB-derived sequence
        2. UniProt sequence
    - mapping is position-based, not original PDB author-residue-ID based
    """
    records = read_fasta_records(alignment_fasta_path)

    if len(records) != 2:
        raise ValueError(
            f"Expected exactly 2 aligned FASTA records in {alignment_fasta_path}, "
            f"but found {len(records)}."
        )

    pdb_header, pdb_aligned_sequence = records[0]
    uniprot_header, uniprot_aligned_sequence = records[1]

    mapping_rows = build_pairwise_alignment_mapping_rows(
        pdb_aligned_sequence=pdb_aligned_sequence,
        uniprot_aligned_sequence=uniprot_aligned_sequence,
    )

    output_mapping_tsv_path.parent.mkdir(parents=True, exist_ok=True)

    header_lines = [
        f"# alignment_fasta={alignment_fasta_path}",
        f"# pdb_sequence_header={pdb_header}",
        f"# uniprot_sequence_header={uniprot_header}",
        (
            "# columns=alignment_column\tpdb_residue_number\tuniprot_residue_number\t"
            "pdb_residue\tuniprot_residue\trelation"
        ),
    ]

    data_lines = [
        "\t".join(
            [
                str(row["alignment_column"]),
                _format_optional_int(row["pdb_residue_number"]),
                _format_optional_int(row["uniprot_residue_number"]),
                row["pdb_residue"],
                row["uniprot_residue"],
                row["relation"],
            ]
        )
        for row in mapping_rows
    ]

    output_text = "\n".join(header_lines + data_lines) + "\n"
    output_mapping_tsv_path.write_text(output_text, encoding="utf-8")

    print(f"[INFO] Mapping file written to {output_mapping_tsv_path}")


def build_pairwise_alignment_mapping_rows(
    pdb_aligned_sequence: str,
    uniprot_aligned_sequence: str,
) -> list[dict[str, Any]]:
    """
    Build one row per alignment column.

    Row semantics
    -------------
    - pdb_residue_number:
        1-based running residue number in the non-gap PDB sequence
    - uniprot_residue_number:
        1-based running residue number in the non-gap UniProt sequence
    - relation:
        one of:
            match
            mismatch
            insertion_in_pdb
            deletion_in_pdb

    Notes
    -----
    - insertion_in_pdb:
        residue present in PDB sequence, gap in UniProt
    - deletion_in_pdb:
        gap in PDB sequence, residue present in UniProt
    """
    if len(pdb_aligned_sequence) != len(uniprot_aligned_sequence):
        raise ValueError("Aligned sequences must have identical length.")

    pdb_residue_number = 0
    uniprot_residue_number = 0
    mapping_rows: list[dict[str, Any]] = []

    for index, (pdb_residue, uniprot_residue) in enumerate(
        zip(pdb_aligned_sequence, uniprot_aligned_sequence, strict=True),
        start=1,
    ):
        pdb_position_for_row: int | None = None
        uniprot_position_for_row: int | None = None

        if pdb_residue != "-":
            pdb_residue_number += 1
            pdb_position_for_row = pdb_residue_number

        if uniprot_residue != "-":
            uniprot_residue_number += 1
            uniprot_position_for_row = uniprot_residue_number

        relation = classify_alignment_column(
            pdb_residue=pdb_residue,
            uniprot_residue=uniprot_residue,
        )

        mapping_rows.append(
            {
                "alignment_column": index,
                "pdb_residue_number": pdb_position_for_row,
                "uniprot_residue_number": uniprot_position_for_row,
                "pdb_residue": pdb_residue,
                "uniprot_residue": uniprot_residue,
                "relation": relation,
            }
        )

    return mapping_rows


def classify_alignment_column(
    pdb_residue: str,
    uniprot_residue: str,
) -> str:
    """
    Classify one alignment column.
    """
    if pdb_residue == "-" and uniprot_residue == "-":
        return "double_gap"

    if pdb_residue == "-" and uniprot_residue != "-":
        return "deletion_in_pdb"

    if pdb_residue != "-" and uniprot_residue == "-":
        return "insertion_in_pdb"

    if pdb_residue == uniprot_residue:
        return "match"

    return "mismatch"


def read_fasta_records(fasta_path: Path) -> list[tuple[str, str]]:
    """
    Read a FASTA file into a list of (header, sequence) tuples.
    """
    if not fasta_path.exists():
        raise FileNotFoundError(f"FASTA file not found: {fasta_path}")

    lines = fasta_path.read_text(encoding="utf-8").splitlines()

    records: list[tuple[str, str]] = []
    current_header: str | None = None
    current_sequence_chunks: list[str] = []

    for raw_line in lines:
        line = raw_line.strip()
        if not line:
            continue

        if line.startswith(">"):
            if current_header is not None:
                records.append((current_header, "".join(current_sequence_chunks)))

            current_header = line[1:].strip()
            current_sequence_chunks = []
        else:
            current_sequence_chunks.append(line)

    if current_header is not None:
        records.append((current_header, "".join(current_sequence_chunks)))

    return records


def _format_optional_int(value: int | None) -> str:
    """
    Format an optional integer for TSV output.
    """
    if value is None:
        return ""
    return str(value)


# ---------------------------------------------------------------------------
# Optional batch helper
# ---------------------------------------------------------------------------


def run_alignments_for_all_pdb_directories(
    protein_root_directory: Path,
    render_images: bool = True,
) -> None:
    """
    Run standard alignment jobs for all PDB subdirectories in a root protein folder.
    """
    if not protein_root_directory.exists():
        raise FileNotFoundError(
            f"Protein root directory not found: {protein_root_directory}"
        )

    pdb_directory_list = sorted(
        path for path in protein_root_directory.iterdir() if path.is_dir()
    )

    if not pdb_directory_list:
        print(f"[WARNING] No PDB subdirectories found in {protein_root_directory}")
        return

    for pdb_directory in pdb_directory_list:
        run_alignments_for_pdb_directory(
            pdb_directory=pdb_directory,
            render_images=render_images,
        )
