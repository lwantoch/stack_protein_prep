# /home/grheco/repositorios/stack_protein_prep/src/stack_protein_preparation/fasta_files.py

"""
Create FASTA files for each local PDB entry directory.

Responsibilities
----------------
- read a downloaded PDB file from disk
- extract a declared sequence from SEQRES records
- extract an observed/resolved sequence from ATOM/HETATM records
- optionally restrict output to chains relevant for a requested residue range
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
Unknown or non-standard residues that are actually present as residues in the
PDB are encoded as 'X'. They are NOT skipped.

However:
- missing regions are NOT encoded as 'X' here
- sequence gaps must emerge later from sequence alignment as '-'

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

from stack_protein_preparation.pipeline_logging import (
    append_key_value_block,
    append_log_text,
    build_module_log_paths,
    print_double_header_box,
    print_light_table_box,
    print_mixed_box,
    write_log_header,
)

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

MODULE_NAME = "fasta_files"

UNIPROT_SEARCH_BASE_URL = "https://rest.uniprot.org/uniprotkb/search"
PDB_TO_UNIPROT_QUERY_TEMPLATE = "xref:pdb-{pdb_id}"

THREE_TO_ONE_AA = {
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
    "SER": "S",
    "THR": "T",
    "TRP": "W",
    "TYR": "Y",
    "VAL": "V",
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
# Logging / terminal helpers
# ---------------------------------------------------------------------------


def _build_log_path_for_pdb_directory(pdb_directory: Path) -> Path:
    """
    Build the central per-protein log path for this module.

    Expected:
        data/proteins/<PDB_ID>/
    Produces:
        data/proteins/logs/<PDB_ID>/fasta_files.log
    """
    pdb_directory = pdb_directory.resolve()
    protein_data_dir = pdb_directory.parent
    pdb_id = pdb_directory.name.upper()

    log_paths = build_module_log_paths(
        protein_data_dir=protein_data_dir,
        pdb_id=pdb_id,
        module_name=MODULE_NAME,
        variant_label=None,
    )
    return log_paths.module_log_path


def _print_start_box(
    *,
    pdb_id: str,
    pdb_file_path: Path,
    fasta_directory: Path,
    residue_range: str,
    log_path: Path,
) -> None:
    print_mixed_box(
        title=f"{MODULE_NAME} · START",
        line_list=[
            f"PDB ID      : {pdb_id}",
            f"PDB file    : {pdb_file_path}",
            f"FASTA dir   : {fasta_directory}",
            f"Range       : {residue_range or '(none)'}",
            f"Log file    : {log_path}",
        ],
    )


def _print_chain_selection_box(
    *,
    pdb_id: str,
    relevant_chain_ids: tuple[str, ...],
) -> None:
    print_light_table_box(
        title=f"{MODULE_NAME} · chain selection · {pdb_id}",
        row_list=[
            (
                "Selected chains",
                ", ".join(relevant_chain_ids) if relevant_chain_ids else "(all)",
            )
        ],
        left_header="Field",
        right_header="Value",
    )


def _print_outputs_table(
    *,
    pdb_id: str,
    seqres_fasta_output_path: Path | None,
    atom_fasta_output_path: Path | None,
    uniprot_output_path_list: list[Path],
) -> None:
    row_list: list[tuple[str, str]] = [
        (
            "SEQRES FASTA",
            str(seqres_fasta_output_path)
            if seqres_fasta_output_path
            else "(not written)",
        ),
        (
            "ATOM FASTA",
            str(atom_fasta_output_path) if atom_fasta_output_path else "(not written)",
        ),
        ("UniProt FASTA count", str(len(uniprot_output_path_list))),
    ]

    if uniprot_output_path_list:
        for index, uniprot_output_path in enumerate(uniprot_output_path_list, start=1):
            row_list.append((f"UniProt file {index}", str(uniprot_output_path)))

    print_light_table_box(
        title=f"{MODULE_NAME} · outputs · {pdb_id}",
        row_list=row_list,
        left_header="Output",
        right_header="Path / Value",
    )


def _print_end_box(
    *,
    pdb_id: str,
    status: str,
    seqres_chain_count: int,
    atom_chain_count: int,
    uniprot_entry_count: int,
) -> None:
    print_double_header_box(
        title=f"{MODULE_NAME} · END",
        line_list=[
            f"PDB ID          : {pdb_id}",
            f"Status          : {status}",
            f"SEQRES chains   : {seqres_chain_count}",
            f"ATOM chains     : {atom_chain_count}",
            f"UniProt entries : {uniprot_entry_count}",
        ],
    )


# ---------------------------------------------------------------------------
# Residue conversion helper
# ---------------------------------------------------------------------------


def convert_residue_name_to_one_letter(residue_name: str) -> str:
    """
    Convert a PDB residue name to a one-letter amino-acid code.

    Rules
    -----
    - standard and supported modified amino acids map to a real one-letter code
    - unknown or unsupported present residues map to 'X'
    - missing residues are NOT represented here and must not become 'X'
    """
    normalized_residue_name = residue_name.upper()
    return THREE_TO_ONE_AA.get(normalized_residue_name, "X")


# ---------------------------------------------------------------------------
# Range / chain helpers
# ---------------------------------------------------------------------------


def _parse_chain_qualified_range(
    residue_range: str,
) -> tuple[str, int, str, int] | None:
    """Return (start_chain, start_resnum, end_chain, end_resnum) for 'A16-B50'.

    Returns None when the range has no chain qualifiers.
    """
    m = re.fullmatch(r"([A-Za-z])(-?\d+)-([A-Za-z])(-?\d+)", str(residue_range).strip())
    if m:
        return m.group(1).upper(), int(m.group(2)), m.group(3).upper(), int(m.group(4))
    return None


def parse_residue_range(residue_range: str) -> tuple[int, int]:
    """
    Parse a residue range like '33-480' or 'A16-B50' (chain-qualified).
    """
    stripped = str(residue_range).strip()
    # Chain-qualified format: A33-B480 — strip chain letters, keep residue numbers
    m = re.fullmatch(r"[A-Za-z](-?\d+)-[A-Za-z](-?\d+)", stripped)
    if m:
        return int(m.group(1)), int(m.group(2))
    # Legacy format: 33-480
    match = re.fullmatch(r"\s*(-?\d+)\s*-\s*(-?\d+)\s*", stripped)
    if match is None:
        raise ValueError(f"Invalid residue range: {residue_range!r}")

    start_residue = int(match.group(1))
    end_residue = int(match.group(2))

    if start_residue > end_residue:
        raise ValueError(f"Invalid residue range: start > end in {residue_range!r}")

    return start_residue, end_residue


def collect_atom_residue_numbers_by_chain(
    pdb_file_path: Path,
) -> dict[str, set[int]]:
    """
    Collect observed residue numbers per chain from ATOM records.

    Only protein-like residues are considered.
    """
    residue_numbers_by_chain: dict[str, set[int]] = {}

    with pdb_file_path.open("r", encoding="utf-8") as pdb_handle:
        for raw_line in pdb_handle:
            if not raw_line.startswith("ATOM"):
                continue

            residue_name = raw_line[17:20].strip().upper()
            if residue_name not in THREE_TO_ONE_AA and residue_name not in {"UNK"}:
                continue

            chain_id = raw_line[21].strip() or "_"
            residue_number_text = raw_line[22:26].strip()

            if not residue_number_text.lstrip("-").isdigit():
                continue

            residue_number = int(residue_number_text)
            residue_numbers_by_chain.setdefault(chain_id, set()).add(residue_number)

    return residue_numbers_by_chain


def choose_relevant_chain_ids_from_range(
    pdb_file_path: Path,
    residue_range: str,
) -> tuple[str, ...]:
    """
    Choose the chain(s) most relevant for the requested residue range.

    Strategy
    --------
    - compute overlap between requested range and observed ATOM residue numbers
      for each chain
    - keep all chains tied for best overlap
    - if no chain overlaps, return all chains present in the PDB-derived map

    Why return multiple chains on ties?
    -----------------------------------
    At FASTA-generation stage we prefer to preserve plausible candidates rather
    than prematurely discarding one biologically relevant chain.
    Later stages can refine to a single modelling chain if needed.
    """
    residue_numbers_by_chain = collect_atom_residue_numbers_by_chain(pdb_file_path)

    if not residue_numbers_by_chain:
        return ()

    # When the range has explicit chain qualifiers (A16-B50), return those chains directly.
    chain_qual = _parse_chain_qualified_range(residue_range)
    if chain_qual is not None:
        start_chain, _, end_chain, _ = chain_qual
        explicit_chains = tuple(
            c for c in sorted({start_chain, end_chain})
            if c in residue_numbers_by_chain
        )
        if explicit_chains:
            return explicit_chains

    range_start, range_end = parse_residue_range(residue_range)
    requested_range = set(range(range_start, range_end + 1))

    score_rows: list[tuple[str, int, int]] = []
    for chain_id, residue_number_set in sorted(residue_numbers_by_chain.items()):
        overlap_count = len(residue_number_set & requested_range)
        total_chain_residue_count = len(residue_number_set)
        score_rows.append((chain_id, overlap_count, total_chain_residue_count))

    best_overlap = max(row[1] for row in score_rows)
    if best_overlap == 0:
        return tuple(sorted(residue_numbers_by_chain.keys()))

    best_rows = [row for row in score_rows if row[1] == best_overlap]
    return tuple(sorted(row[0] for row in best_rows))


def filter_chain_sequence_dict(
    chain_sequence_dict: dict[str, str],
    allowed_chain_ids: tuple[str, ...] | None,
) -> dict[str, str]:
    """
    Restrict a chain->sequence dictionary to selected chains.

    If allowed_chain_ids is empty or None, keep all chains.
    """
    if not allowed_chain_ids:
        return dict(chain_sequence_dict)

    allowed_set = set(allowed_chain_ids)
    return {
        chain_id: sequence
        for chain_id, sequence in chain_sequence_dict.items()
        if chain_id in allowed_set
    }


# ---------------------------------------------------------------------------
# Public high-level function
# ---------------------------------------------------------------------------


def create_fasta_files_for_pdb_directory(
    pdb_directory: Path,
    residue_range: str = "",
) -> None:
    """
    Create local PDB-derived FASTA and UniProt FASTA files for one PDB directory.

    Expected layout before running:
        data/proteins/<PDB_ID>/<PDB_ID>.pdb

    Expected layout after running:
        data/proteins/<PDB_ID>/fasta/PDB-<PDB_ID>-SEQRES.fasta
        data/proteins/<PDB_ID>/fasta/PDB-<PDB_ID>-ATOM.fasta
        data/proteins/<PDB_ID>/fasta/UniProt_<UNIPROT_ID>.fasta
        ...

    Parameters
    ----------
    pdb_directory:
        Protein directory containing <PDB_ID>.pdb
    residue_range:
        Optional requested residue range such as '33-480'. If provided, only
        relevant chain(s) for that range are written to the local PDB FASTA
        files.
    """
    pdb_id = pdb_directory.name.upper()
    pdb_file_path = pdb_directory / f"{pdb_id}.pdb"
    fasta_directory = pdb_directory / "fasta"
    log_path = _build_log_path_for_pdb_directory(pdb_directory)

    write_log_header(
        log_path=log_path,
        module_name=MODULE_NAME,
        pdb_id=pdb_id,
        extra_lines=[
            f"pdb_directory={pdb_directory}",
            f"pdb_file_path={pdb_file_path}",
            f"fasta_directory={fasta_directory}",
            f"residue_range={residue_range or '(none)'}",
        ],
    )

    _print_start_box(
        pdb_id=pdb_id,
        pdb_file_path=pdb_file_path,
        fasta_directory=fasta_directory,
        residue_range=residue_range,
        log_path=log_path,
    )

    if not pdb_file_path.exists():
        append_log_text(log_path, f"[WARNING] PDB file not found: {pdb_file_path}")
        print_double_header_box(
            title=f"{MODULE_NAME} · END",
            line_list=[
                f"PDB ID      : {pdb_id}",
                "Status      : missing_input",
                f"PDB file    : {pdb_file_path}",
            ],
        )
        return

    fasta_directory.mkdir(parents=True, exist_ok=True)

    relevant_chain_ids: tuple[str, ...] = ()
    if str(residue_range).strip():
        try:
            relevant_chain_ids = choose_relevant_chain_ids_from_range(
                pdb_file_path=pdb_file_path,
                residue_range=residue_range,
            )
            append_log_text(
                log_path,
                "[CHAIN SELECTION]",
            )
            append_log_text(
                log_path,
                f"relevant_chain_ids={','.join(relevant_chain_ids) if relevant_chain_ids else '(all)'}",
            )
            _print_chain_selection_box(
                pdb_id=pdb_id,
                relevant_chain_ids=relevant_chain_ids,
            )
        except Exception as error:
            append_log_text(
                log_path,
                (
                    "[WARNING] Could not determine relevant chains from range "
                    f"{residue_range!r}: {error!r}. Writing all chains instead."
                ),
            )
            relevant_chain_ids = ()

    seqres_chain_sequence_dict = extract_seqres_sequences_from_pdb(pdb_file_path)
    seqres_chain_sequence_dict = filter_chain_sequence_dict(
        chain_sequence_dict=seqres_chain_sequence_dict,
        allowed_chain_ids=relevant_chain_ids,
    )

    seqres_fasta_output_path: Path | None = None
    if seqres_chain_sequence_dict:
        seqres_fasta_output_path = fasta_directory / f"PDB-{pdb_id}-SEQRES.fasta"
        write_pdb_chain_fasta(
            fasta_source_label=f"PDB|{pdb_id}|SEQRES",
            chain_sequence_dict=seqres_chain_sequence_dict,
            output_fasta_path=seqres_fasta_output_path,
        )
        append_key_value_block(
            log_path=log_path,
            title="SEQRES FASTA",
            key_value_pairs=[
                ("chain_count", str(len(seqres_chain_sequence_dict))),
                ("chains", ", ".join(sorted(seqres_chain_sequence_dict)) or "(none)"),
                ("output_path", str(seqres_fasta_output_path)),
            ],
        )
    else:
        append_log_text(
            log_path, f"[WARNING] No SEQRES sequence found in {pdb_file_path}"
        )

    atom_chain_sequence_dict = extract_observed_atom_sequences_from_pdb(pdb_file_path)
    atom_chain_sequence_dict = filter_chain_sequence_dict(
        chain_sequence_dict=atom_chain_sequence_dict,
        allowed_chain_ids=relevant_chain_ids,
    )

    atom_fasta_output_path: Path | None = None
    if atom_chain_sequence_dict:
        atom_fasta_output_path = fasta_directory / f"PDB-{pdb_id}-ATOM.fasta"
        write_pdb_chain_fasta(
            fasta_source_label=f"PDB|{pdb_id}|ATOM",
            chain_sequence_dict=atom_chain_sequence_dict,
            output_fasta_path=atom_fasta_output_path,
        )
        append_key_value_block(
            log_path=log_path,
            title="ATOM FASTA",
            key_value_pairs=[
                ("chain_count", str(len(atom_chain_sequence_dict))),
                ("chains", ", ".join(sorted(atom_chain_sequence_dict)) or "(none)"),
                ("output_path", str(atom_fasta_output_path)),
            ],
        )
    else:
        append_log_text(
            log_path, f"[WARNING] No ATOM-derived sequence found in {pdb_file_path}"
        )

    uniprot_fasta_entries = fetch_uniprot_fasta_entries_for_pdb_id(
        pdb_id=pdb_id,
        log_path=log_path,
    )

    uniprot_output_path_list: list[Path] = []

    if not uniprot_fasta_entries:
        append_log_text(
            log_path, f"[WARNING] No UniProt FASTA entries found for PDB {pdb_id}"
        )
        _print_outputs_table(
            pdb_id=pdb_id,
            seqres_fasta_output_path=seqres_fasta_output_path,
            atom_fasta_output_path=atom_fasta_output_path,
            uniprot_output_path_list=uniprot_output_path_list,
        )
        _print_end_box(
            pdb_id=pdb_id,
            status="partial_success",
            seqres_chain_count=len(seqres_chain_sequence_dict),
            atom_chain_count=len(atom_chain_sequence_dict),
            uniprot_entry_count=0,
        )
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
        uniprot_output_path_list.append(uniprot_output_path)

    append_key_value_block(
        log_path=log_path,
        title="UNIPROT FASTA",
        key_value_pairs=[
            ("entry_count", str(len(uniprot_fasta_entries))),
            (
                "accessions",
                ", ".join(entry.accession for entry in uniprot_fasta_entries)
                or "(none)",
            ),
        ],
    )

    _print_outputs_table(
        pdb_id=pdb_id,
        seqres_fasta_output_path=seqres_fasta_output_path,
        atom_fasta_output_path=atom_fasta_output_path,
        uniprot_output_path_list=uniprot_output_path_list,
    )

    _print_end_box(
        pdb_id=pdb_id,
        status="success",
        seqres_chain_count=len(seqres_chain_sequence_dict),
        atom_chain_count=len(atom_chain_sequence_dict),
        uniprot_entry_count=len(uniprot_fasta_entries),
    )


# ---------------------------------------------------------------------------
# PDB -> FASTA (SEQRES)
# ---------------------------------------------------------------------------


def extract_seqres_sequences_from_pdb(pdb_file_path: Path) -> dict[str, str]:
    """
    Extract declared polymer sequences from SEQRES records.

    Important
    ---------
    - This reflects the full declared sequence in the PDB entry.
    - Residues that are unknown or unsupported but explicitly present in SEQRES
      are encoded as 'X'.
    - This function does NOT tell us which residues are actually resolved in the
      3D structure.
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

    This sequence reflects only residues actually present in the model.

    Important semantics
    -------------------
    - one residue entry is emitted only if the residue is present in the PDB
    - unknown or unsupported present residues become 'X'
    - missing regions are absent here and must later appear as '-' only after
      sequence alignment
    - duplicate atoms for one residue are collapsed to one residue entry
    """
    chain_residue_entries: dict[str, list[tuple[str, str, str]]] = {}
    seen_residue_site_keys: set[tuple[str, str, str]] = set()

    with pdb_file_path.open("r", encoding="utf-8") as pdb_handle:
        for raw_line in pdb_handle:
            if not raw_line.startswith("ATOM"):
                continue

            residue_name = raw_line[17:20].strip().upper()
            chain_id = raw_line[21].strip() or "_"
            residue_number = raw_line[22:26].strip()
            insertion_code = raw_line[26].strip()

            residue_site_key = (chain_id, residue_number, insertion_code)

            if residue_site_key in seen_residue_site_keys:
                continue

            seen_residue_site_keys.add(residue_site_key)
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


def fetch_uniprot_fasta_entries_for_pdb_id(
    pdb_id: str,
    log_path: Path | None = None,
) -> list[UniProtFastaEntry]:
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

    if log_path is not None:
        append_key_value_block(
            log_path=log_path,
            title="UNIPROT QUERY",
            key_value_pairs=[
                ("pdb_id", normalized_pdb_id),
                ("query", uniprot_query),
                ("request_url", request_url),
            ],
        )

    with urllib.request.urlopen(request_url, timeout=30) as response:
        fasta_text = response.read().decode("utf-8")

    parsed_entries = parse_uniprot_fasta_text(fasta_text)

    if log_path is not None:
        append_log_text(log_path, f"[UNIPROT RESULT] entry_count={len(parsed_entries)}")
        for entry in parsed_entries:
            append_log_text(
                log_path,
                f"accession={entry.accession} header={entry.header}",
            )
        append_log_text(log_path, "")

    return parsed_entries


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
