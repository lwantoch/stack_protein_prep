# /home/grheco/repositorios/stack_protein_prep/src/stack_protein_preparation/fasta_files.py

"""
Create FASTA files for each local PDB entry directory.

Responsibilities
----------------
- read a downloaded PDB file from disk
- extract a declared sequence from SEQRES records
- extract an observed/resolved sequence from ATOM/HETATM records
- write local PDB-derived FASTA files
- query UniProt entries linked to the PDB ID
- write one FASTA file per UniProt accession

Important design choice
-----------------------
This module writes:
    data/proteins/<PDB_ID>/fasta/PDB-<PDB_ID>-SEQRES.fasta
    data/proteins/<PDB_ID>/fasta/PDB-<PDB_ID>-ATOM.fasta
and:
    data/proteins/<PDB_ID>/fasta/UniProt_<UNIPROT_ID>.fasta

Why two local PDB FASTA files?
------------------------------
- SEQRES reflects the declared full polymer sequence in the PDB entry
- ATOM/HETATM reflects only residues actually present in the solved model

This distinction is important when later checking:
- missing residues
- unresolved loops or termini
- mutations
- engineered constructs
- sequence mismatches versus UniProt

Important residue handling rule
-------------------------------
Unknown or non-standard residues are encoded as 'X'.
They are NOT skipped.

This preserves sequence length better and makes later comparisons more
informative.

What this module does NOT do
----------------------------
- It does not align PDB vs UniProt sequences.
- It does not decide which UniProt entry is the "best" one.
- It does not validate biological correctness of the mapping.
- It does not parse mmCIF yet.

Testing note
------------
Network calls should be mocked in unit tests.
"""

from __future__ import annotations

import re
import urllib.parse
import urllib.request
from dataclasses import dataclass
from pathlib import Path

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

UNIPROT_SEARCH_BASE_URL = "https://rest.uniprot.org/uniprotkb/search"
PDB_TO_UNIPROT_QUERY_TEMPLATE = "xref:pdb-{pdb_id}"

THREE_TO_ONE_AA = {
    "ALA": "A",
    "ARG": "R",
    "ASN": "N",
    "ASP": "D",
    "CYS": "C",
    "GLN": "Q",
    "GLU": "E",
    "GLY": "G",
    "HIS": "H",
    "ILE": "I",
    "LEU": "L",
    "LYS": "K",
    "MET": "M",
    "PHE": "F",
    "PRO": "P",
    "SER": "S",
    "THR": "T",
    "TRP": "W",
    "TYR": "Y",
    "VAL": "V",
    "HID": "H",
    "HIE": "H",
    "HIP": "H",
    "CYM": "C",
    "CYX": "C",
}


# ---------------------------------------------------------------------------
# Small data containers
# ---------------------------------------------------------------------------


@dataclass
class UniProtFastaEntry:
    """
    One UniProt FASTA entry as returned by the UniProt REST API search endpoint.
    """

    accession: str
    header: str
    sequence: str


# ---------------------------------------------------------------------------
# Residue conversion helper
# ---------------------------------------------------------------------------


def convert_residue_name_to_one_letter(residue_name: str) -> str:
    """
    Convert a PDB residue name to a one-letter amino-acid code.

    Unknown or non-standard residues are converted to 'X' instead of being
    skipped. This preserves sequence information better for later comparison.
    """
    normalized_residue_name = residue_name.upper()
    return THREE_TO_ONE_AA.get(normalized_residue_name, "X")


# ---------------------------------------------------------------------------
# Public high-level function
# ---------------------------------------------------------------------------


def create_fasta_files_for_pdb_directory(pdb_directory: Path) -> None:
    """
    Create local PDB-derived FASTA and UniProt FASTA files for one PDB directory.

    Expected layout before running:
        data/proteins/<PDB_ID>/<PDB_ID>.pdb

    Expected layout after running:
        data/proteins/<PDB_ID>/fasta/PDB-<PDB_ID>-SEQRES.fasta
        data/proteins/<PDB_ID>/fasta/PDB-<PDB_ID>-ATOM.fasta
        data/proteins/<PDB_ID>/fasta/UniProt_<UNIPROT_ID>.fasta
        ...
    """
    pdb_id = pdb_directory.name.upper()
    pdb_file_path = pdb_directory / f"{pdb_id}.pdb"
    fasta_directory = pdb_directory / "fasta"

    if not pdb_file_path.exists():
        print(f"[WARNING] PDB file not found for {pdb_id}: {pdb_file_path}")
        return

    fasta_directory.mkdir(parents=True, exist_ok=True)

    print(f"[INFO] Creating FASTA files for {pdb_id}")

    seqres_chain_sequence_dict = extract_seqres_sequences_from_pdb(pdb_file_path)

    if seqres_chain_sequence_dict:
        seqres_fasta_output_path = fasta_directory / f"PDB-{pdb_id}-SEQRES.fasta"
        write_pdb_chain_fasta(
            fasta_source_label=f"PDB|{pdb_id}|SEQRES",
            chain_sequence_dict=seqres_chain_sequence_dict,
            output_fasta_path=seqres_fasta_output_path,
        )
        print(f"[INFO] Wrote SEQRES FASTA: {seqres_fasta_output_path}")
    else:
        print(f"[WARNING] No SEQRES sequence found in {pdb_file_path}")

    atom_chain_sequence_dict = extract_observed_atom_sequences_from_pdb(pdb_file_path)

    if atom_chain_sequence_dict:
        atom_fasta_output_path = fasta_directory / f"PDB-{pdb_id}-ATOM.fasta"
        write_pdb_chain_fasta(
            fasta_source_label=f"PDB|{pdb_id}|ATOM",
            chain_sequence_dict=atom_chain_sequence_dict,
            output_fasta_path=atom_fasta_output_path,
        )
        print(f"[INFO] Wrote ATOM FASTA: {atom_fasta_output_path}")
    else:
        print(f"[WARNING] No ATOM-derived sequence found in {pdb_file_path}")

    uniprot_fasta_entries = fetch_uniprot_fasta_entries_for_pdb_id(pdb_id)

    if not uniprot_fasta_entries:
        print(f"[WARNING] No UniProt FASTA entries found for PDB {pdb_id}")
        return

    for uniprot_entry in uniprot_fasta_entries:
        uniprot_output_path = (
            fasta_directory / f"UniProt_{uniprot_entry.accession}.fasta"
        )
        write_single_fasta_entry(
            header=uniprot_entry.header,
            sequence=uniprot_entry.sequence,
            output_fasta_path=uniprot_output_path,
        )
        print(f"[INFO] Wrote UniProt FASTA: {uniprot_output_path}")


# ---------------------------------------------------------------------------
# PDB -> FASTA (SEQRES)
# ---------------------------------------------------------------------------


def extract_seqres_sequences_from_pdb(pdb_file_path: Path) -> dict[str, str]:
    """
    Extract declared polymer sequences from SEQRES records.

    Important
    ---------
    - This reflects the full declared sequence in the PDB entry.
    - Residues that are unknown or non-standard are encoded as 'X'.
    - This function does NOT tell us which residues are actually resolved
      in the 3D structure.
    """
    chain_residue_name_dict: dict[str, list[str]] = {}

    with pdb_file_path.open("r", encoding="utf-8") as pdb_handle:
        for raw_line in pdb_handle:
            if not raw_line.startswith("SEQRES"):
                continue

            chain_id = raw_line[11].strip()
            if not chain_id:
                chain_id = "_"

            residue_name_block = raw_line[19:70]
            residue_name_list = residue_name_block.split()

            chain_residue_name_dict.setdefault(chain_id, []).extend(residue_name_list)

    chain_sequence_dict: dict[str, str] = {}

    for chain_id, residue_name_list in chain_residue_name_dict.items():
        one_letter_sequence = "".join(
            convert_residue_name_to_one_letter(residue_name)
            for residue_name in residue_name_list
        )

        if one_letter_sequence:
            chain_sequence_dict[chain_id] = one_letter_sequence

    return chain_sequence_dict


# ---------------------------------------------------------------------------
# PDB -> FASTA (ATOM/HETATM)
# ---------------------------------------------------------------------------


def extract_observed_atom_sequences_from_pdb(pdb_file_path: Path) -> dict[str, str]:
    """
    Extract observed/resolved sequences from ATOM/HETATM records.

    This sequence reflects only residues that are actually present in the model.
    Missing residues will therefore be absent here.

    Notes
    -----
    - Residue identity is inferred from ATOM/HETATM records.
    - Duplicate atoms of the same residue are collapsed to one residue entry.
    - Unknown/non-standard residues are encoded as 'X'.
    """
    chain_residue_entries: dict[str, list[tuple[str, str, str]]] = {}
    seen_residue_keys: set[tuple[str, str, str, str]] = set()

    with pdb_file_path.open("r", encoding="utf-8") as pdb_handle:
        for raw_line in pdb_handle:
            if not raw_line.startswith(("ATOM", "HETATM")):
                continue

            residue_name = raw_line[17:20].strip().upper()
            chain_id = raw_line[21].strip() or "_"
            residue_number = raw_line[22:26].strip()
            insertion_code = raw_line[26].strip()

            residue_key = (chain_id, residue_number, insertion_code, residue_name)

            if residue_key in seen_residue_keys:
                continue

            seen_residue_keys.add(residue_key)
            chain_residue_entries.setdefault(chain_id, []).append(
                (residue_number, insertion_code, residue_name)
            )

    chain_sequence_dict: dict[str, str] = {}

    for chain_id, residue_entry_list in chain_residue_entries.items():
        one_letter_sequence = "".join(
            convert_residue_name_to_one_letter(residue_name)
            for _, _, residue_name in residue_entry_list
        )

        if one_letter_sequence:
            chain_sequence_dict[chain_id] = one_letter_sequence

    return chain_sequence_dict


# ---------------------------------------------------------------------------
# FASTA writing for local PDB-derived sequences
# ---------------------------------------------------------------------------


def write_pdb_chain_fasta(
    fasta_source_label: str,
    chain_sequence_dict: dict[str, str],
    output_fasta_path: Path,
) -> None:
    """
    Write one FASTA file containing one entry per chain from the local PDB.

    Example headers:
        >PDB|1W4R|SEQRES|chain_A
        >PDB|1W4R|ATOM|chain_B
    """
    with output_fasta_path.open("w", encoding="utf-8") as fasta_handle:
        for chain_id in sorted(chain_sequence_dict):
            sequence = chain_sequence_dict[chain_id]
            header = f"{fasta_source_label}|chain_{chain_id}"
            fasta_handle.write(format_fasta_record(header, sequence))


# ---------------------------------------------------------------------------
# UniProt -> FASTA
# ---------------------------------------------------------------------------


def fetch_uniprot_fasta_entries_for_pdb_id(pdb_id: str) -> list[UniProtFastaEntry]:
    """
    Query UniProt for entries cross-referenced to a PDB ID and return FASTA entries.
    """
    normalized_pdb_id = pdb_id.upper()
    uniprot_query = PDB_TO_UNIPROT_QUERY_TEMPLATE.format(pdb_id=normalized_pdb_id)

    query_parameters = {
        "query": uniprot_query,
        "format": "fasta",
    }

    request_url = (
        f"{UNIPROT_SEARCH_BASE_URL}?{urllib.parse.urlencode(query_parameters)}"
    )

    print(f"[INFO] Querying UniProt for PDB {normalized_pdb_id}")
    print(f"[DEBUG] UniProt query URL: {request_url}")

    with urllib.request.urlopen(request_url, timeout=30) as response:
        fasta_text = response.read().decode("utf-8")

    return parse_uniprot_fasta_text(fasta_text)


def parse_uniprot_fasta_text(fasta_text: str) -> list[UniProtFastaEntry]:
    """
    Parse FASTA text returned by UniProt into structured entries.
    """
    fasta_text = fasta_text.strip()

    if not fasta_text:
        return []

    raw_entries = [entry for entry in fasta_text.split(">") if entry.strip()]
    parsed_entries: list[UniProtFastaEntry] = []

    for raw_entry in raw_entries:
        entry_lines = raw_entry.strip().splitlines()
        raw_header = entry_lines[0].strip()
        sequence = "".join(line.strip() for line in entry_lines[1:])

        accession = extract_uniprot_accession_from_header(raw_header)

        parsed_entries.append(
            UniProtFastaEntry(
                accession=accession,
                header=raw_header,
                sequence=sequence,
            )
        )

    return parsed_entries


def extract_uniprot_accession_from_header(fasta_header: str) -> str:
    """
    Extract UniProt accession from a FASTA header.

    Typical examples:
        sp|P69905|HBA_HUMAN Hemoglobin subunit alpha ...
        tr|A0A...|...
    """
    header_fields = fasta_header.split("|")

    if len(header_fields) >= 3:
        return header_fields[1].strip()

    safe_token = re.sub(r"[^A-Za-z0-9_.-]+", "_", fasta_header)[:40]
    return safe_token or "UNKNOWN_UNIPROT_ID"


# ---------------------------------------------------------------------------
# Generic FASTA writing helpers
# ---------------------------------------------------------------------------


def write_single_fasta_entry(
    header: str,
    sequence: str,
    output_fasta_path: Path,
) -> None:
    """
    Write exactly one FASTA entry to disk.
    """
    with output_fasta_path.open("w", encoding="utf-8") as fasta_handle:
        fasta_handle.write(format_fasta_record(header, sequence))


def format_fasta_record(header: str, sequence: str, line_width: int = 80) -> str:
    """
    Format one FASTA record with fixed sequence line width.
    """
    wrapped_sequence_lines = [
        sequence[i : i + line_width] for i in range(0, len(sequence), line_width)
    ]

    return f">{header}\n" + "\n".join(wrapped_sequence_lines) + "\n"
