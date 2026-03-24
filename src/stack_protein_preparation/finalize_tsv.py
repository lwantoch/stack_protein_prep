"""
/home/grheco/repositorios/stack_protein_prep/src/stack_protein_preparation/finalize_tsv.py

Build a dedicated finalize-numbering TSV for the final protein model against UniProt.

Purpose
-------
- read the final model PDB
- read the existing ATOM-vs-UniProt alignment mapping file
- derive one numbering row per residue in the final model
- assign final residue numbering using UniProt positions

Current numbering policy
------------------------
- final_resseq = uniprot_resseq
- final_icode = blank
- source_category is inferred from alignment relation if possible

Important
---------
- this builder assumes that the residue order in the final model follows the
  UniProt-based target order
- it is intended as the explicit numbering source for finalize_protein.py
"""

from __future__ import annotations

import csv
import re
from pathlib import Path

from Bio.PDB import PDBParser


def _normalize_chain_id(value: str) -> str:
    chain_id = str(value).strip()
    return chain_id if chain_id else "_"


def _normalize_icode(value: str) -> str:
    icode = str(value).strip()
    return icode if icode else " "


def _parse_alignment_mapping_rows(
    alignment_mapping_tsv_path: Path,
) -> list[dict[str, str]]:
    """
    Parse the alignment mapping export written by sequence_alignment.py.

    Expected format
    ---------------
    - comment lines starting with '#'
    - one comment line of the form:
        # columns=alignment_column\tpdb_residue_number\t...
    - following data rows separated by tabs
    """
    if not alignment_mapping_tsv_path.exists():
        raise FileNotFoundError(
            f"Alignment mapping TSV not found: {alignment_mapping_tsv_path}"
        )

    with alignment_mapping_tsv_path.open("r", encoding="utf-8") as handle:
        raw_lines = [line.rstrip("\n") for line in handle if line.strip()]

    header: list[str] | None = None
    data_lines: list[str] = []

    for line in raw_lines:
        stripped = line.strip()

        if stripped.startswith("# columns="):
            header_text = stripped[len("# columns=") :]
            header = header_text.split("\t")
            continue

        if stripped.startswith("#"):
            continue

        data_lines.append(line)

    if header is None:
        raise ValueError(
            f"Missing '# columns=' header line in alignment mapping TSV: {alignment_mapping_tsv_path}"
        )

    if not data_lines:
        raise ValueError(
            f"No data rows found in alignment mapping TSV: {alignment_mapping_tsv_path}"
        )

    row_dict_list: list[dict[str, str]] = []

    for line in data_lines:
        fields = line.split("\t")

        if len(fields) != len(header):
            raise ValueError(
                f"Column count mismatch in {alignment_mapping_tsv_path}. "
                f"Header has {len(header)} columns, row has {len(fields)} fields: {line!r}"
            )

        row_dict_list.append(dict(zip(header, fields)))

    return row_dict_list


def _get_chain_id_from_mapping_filename(mapping_path: Path) -> str:
    match = re.search(r"chain_([A-Za-z0-9])", mapping_path.name)
    if match:
        return _normalize_chain_id(match.group(1))
    return "_"


def _extract_model_residue_count_for_chain(
    final_model_pdb_path: Path,
    chain_id: str,
) -> int:
    """
    Count standard protein residues in the selected chain of the final model.
    """
    if not final_model_pdb_path.exists():
        raise FileNotFoundError(f"Final model PDB not found: {final_model_pdb_path}")

    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("final_model", str(final_model_pdb_path))

    residue_count = 0

    for model in structure:
        for chain in model:
            current_chain_id = _normalize_chain_id(chain.id)
            if current_chain_id != chain_id:
                continue

            for residue in chain:
                if residue.id[0] == " ":
                    residue_count += 1

    if residue_count == 0:
        raise ValueError(
            f"No standard protein residues found in chain {chain_id!r} of {final_model_pdb_path}"
        )

    return residue_count


def build_finalize_tsv_against_uniprot(
    final_model_pdb_path: Path,
    alignment_mapping_tsv_path: Path,
    output_finalize_tsv_path: Path,
    chain_id: str | None = None,
) -> Path:
    """
    Build a finalize-numbering TSV using UniProt positions as final residue numbers.

    Output columns
    --------------
    - chain_id
    - model_seq_position
    - final_resseq
    - final_icode
    - uniprot_resseq
    - source_category

    Policy
    ------
    - only rows with a numeric UniProt residue position are usable
    - final_resseq is set equal to uniprot_resseq
    - rows are taken in alignment order
    - the first N usable rows are assigned to the N residues in the final model
    """
    row_dict_list = _parse_alignment_mapping_rows(alignment_mapping_tsv_path)

    if chain_id is None:
        chain_id = _get_chain_id_from_mapping_filename(alignment_mapping_tsv_path)

    chain_id = _normalize_chain_id(chain_id)

    final_model_residue_count = _extract_model_residue_count_for_chain(
        final_model_pdb_path=final_model_pdb_path,
        chain_id=chain_id,
    )

    usable_row_list: list[dict[str, str]] = []

    for row in row_dict_list:
        uniprot_resseq_text = str(row.get("uniprot_residue_number", "")).strip()
        relation = str(row.get("relation", "")).strip()

        if not uniprot_resseq_text:
            continue
        if uniprot_resseq_text in {"-", "X"}:
            continue

        try:
            int(uniprot_resseq_text)
        except ValueError:
            continue

        source_category = "aligned"
        if relation == "match":
            source_category = "original_pdb"
        elif relation == "deletion_in_pdb":
            source_category = "modeled_gap"
        elif relation == "insertion_in_pdb":
            source_category = "insertion_in_pdb"

        usable_row_list.append(
            {
                "uniprot_resseq": uniprot_resseq_text,
                "source_category": source_category,
            }
        )

    if len(usable_row_list) < final_model_residue_count:
        raise ValueError(
            "Not enough usable UniProt mapping rows to cover the final model. "
            f"usable_rows={len(usable_row_list)}, "
            f"final_model_residues={final_model_residue_count}, "
            f"mapping_file={alignment_mapping_tsv_path}"
        )

    selected_row_list = usable_row_list[:final_model_residue_count]

    output_finalize_tsv_path.parent.mkdir(parents=True, exist_ok=True)

    fieldnames = [
        "chain_id",
        "model_seq_position",
        "final_resseq",
        "final_icode",
        "uniprot_resseq",
        "source_category",
    ]

    with output_finalize_tsv_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()

        for index, row in enumerate(selected_row_list, start=1):
            uniprot_resseq = int(row["uniprot_resseq"])

            writer.writerow(
                {
                    "chain_id": chain_id,
                    "model_seq_position": index,
                    "final_resseq": uniprot_resseq,
                    "final_icode": "",
                    "uniprot_resseq": uniprot_resseq,
                    "source_category": row["source_category"],
                }
            )

    return output_finalize_tsv_path
