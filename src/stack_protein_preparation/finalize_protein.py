"""
/home/grheco/repositorios/stack_protein_prep/src/stack_protein_preparation/finalize_protein.py

Finalize the prepared protein structure by renumbering residues according to a
dedicated finalize-numbering TSV.

Design
------
This module no longer tries to infer final numbering from the older alignment
mapping TSV directly. Instead, it expects a purpose-built finalize TSV that
contains one row per residue in the final model.

Expected finalize TSV columns
-----------------------------
Required:
- chain_id
- model_seq_position
- final_resseq

Optional:
- final_icode
- source_category
- original_resseq
- uniprot_resseq
- final_resname

Behavior
--------
- Reads components/<PDBID>_protein_as_Amber.pdb
- Renumbers only standard protein residues (residue.id[0] == " ")
- Matches residues by:
    (chain_id, running residue index within chain, 1-based)
- Writes:
    components/<PDBID>_protein_final.pdb

Important
---------
- The finalize TSV must already represent the desired final numbering of the
  final model, including modeled residues if present.
- This module does not try to reconstruct numbering policy from alignment logic.
- This module does not alter coordinates or chemistry.
"""

from __future__ import annotations

import csv
from pathlib import Path

from Bio.PDB import PDBIO, PDBParser

REQUIRED_FINALIZE_TSV_COLUMNS = {
    "chain_id",
    "model_seq_position",
    "final_resseq",
}


def _normalize_chain_id(value: str) -> str:
    chain_id = str(value).strip()
    return chain_id if chain_id else "_"


def _normalize_icode(value: str) -> str:
    icode = str(value).strip()
    return icode if icode else " "


def load_finalize_numbering_map_from_tsv(
    finalize_tsv_path: Path,
) -> dict[tuple[str, int], tuple[int, str]]:
    """
    Load the dedicated finalize-numbering TSV.

    Returns
    -------
    dict
        Maps:
            (chain_id, model_seq_position) -> (final_resseq, final_icode)
    """
    if not finalize_tsv_path.exists():
        raise FileNotFoundError(f"Finalize TSV not found: {finalize_tsv_path}")

    numbering_map: dict[tuple[str, int], tuple[int, str]] = {}

    with finalize_tsv_path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")

        if reader.fieldnames is None:
            raise ValueError(f"Finalize TSV has no header: {finalize_tsv_path}")

        missing_columns = REQUIRED_FINALIZE_TSV_COLUMNS.difference(reader.fieldnames)
        if missing_columns:
            raise KeyError(
                f"Missing required finalize TSV columns in {finalize_tsv_path}: "
                f"{sorted(missing_columns)}"
            )

        for row in reader:
            chain_id = _normalize_chain_id(row["chain_id"])
            model_seq_position_text = str(row["model_seq_position"]).strip()
            final_resseq_text = str(row["final_resseq"]).strip()
            final_icode = _normalize_icode(row.get("final_icode", ""))

            try:
                model_seq_position = int(model_seq_position_text)
            except ValueError as error:
                raise ValueError(
                    f"Invalid model_seq_position in {finalize_tsv_path}: "
                    f"{model_seq_position_text!r}"
                ) from error

            try:
                final_resseq = int(final_resseq_text)
            except ValueError as error:
                raise ValueError(
                    f"Invalid final_resseq in {finalize_tsv_path}: "
                    f"{final_resseq_text!r}"
                ) from error

            key = (chain_id, model_seq_position)

            if key in numbering_map:
                raise ValueError(
                    f"Duplicate finalize numbering entry for {key} in {finalize_tsv_path}"
                )

            numbering_map[key] = (final_resseq, final_icode)

    if not numbering_map:
        raise ValueError(
            f"No numbering rows loaded from finalize TSV: {finalize_tsv_path}"
        )

    return numbering_map


def renumber_structure_with_finalize_map(
    input_pdb_path: Path,
    output_pdb_path: Path,
    numbering_map: dict[tuple[str, int], tuple[int, str]],
) -> dict[str, int | Path]:
    """
    Renumber a structure using the dedicated finalize-numbering map.

    Matching rule
    -------------
    Residues are matched by:
    - chain_id
    - running protein residue index within that chain (1-based)

    Only standard protein residues are renumbered:
    - residue.id[0] == " "
    """
    if not input_pdb_path.exists():
        raise FileNotFoundError(f"Input PDB not found: {input_pdb_path}")

    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("final_model", str(input_pdb_path))

    renumbered_residues = 0
    renumbered_atoms = 0
    used_key_set: set[tuple[str, int]] = set()

    for model in structure:
        for chain in model:
            chain_id = _normalize_chain_id(chain.id)
            residue_counter = 0

            for residue in chain:
                if residue.id[0] != " ":
                    continue

                residue_counter += 1
                key = (chain_id, residue_counter)

                if key not in numbering_map:
                    raise KeyError(
                        "Missing finalize numbering entry for "
                        f"chain={chain_id!r}, model_seq_position={residue_counter}"
                    )

                final_resseq, final_icode = numbering_map[key]
                residue.id = (" ", final_resseq, final_icode)

                used_key_set.add(key)
                renumbered_residues += 1
                renumbered_atoms += sum(1 for _ in residue.get_atoms())

    unused_key_list = sorted(set(numbering_map.keys()) - used_key_set)
    if unused_key_list:
        preview = ", ".join(str(key) for key in unused_key_list[:10])
        raise ValueError(
            "Finalize TSV contains entries not used by the structure. "
            f"First unused keys: {preview}"
        )

    output_pdb_path.parent.mkdir(parents=True, exist_ok=True)

    io = PDBIO()
    io.set_structure(structure)
    io.save(str(output_pdb_path))

    return {
        "output_pdb_path": output_pdb_path,
        "renumbered_residues": renumbered_residues,
        "renumbered_atoms": renumbered_atoms,
    }


def build_finalize_tsv_from_alignment_mapping(
    alignment_mapping_tsv_path: Path,
    output_finalize_tsv_path: Path,
    chain_id: str | None = None,
) -> Path:
    """
    Build a dedicated finalize TSV from the older alignment mapping TSV.

    Current policy
    --------------
    - keep only rows with relation == "match"
    - model_seq_position is assigned as a running 1-based index over kept rows
    - final_resseq is taken from pdb_residue_number
    - final_icode is set to blank

    This builder is suitable only for cases where the final model contains
    exactly the original matched PDB residues in order.

    It is NOT sufficient for modeled-gap cases unless a richer numbering policy
    has been encoded upstream.
    """
    if not alignment_mapping_tsv_path.exists():
        raise FileNotFoundError(
            f"Alignment mapping TSV not found: {alignment_mapping_tsv_path}"
        )

    if chain_id is None:
        name = alignment_mapping_tsv_path.name
        marker = "chain_"
        if marker in name:
            after = name.split(marker, 1)[1]
            chain_id = after[0]
        else:
            chain_id = "_"

    chain_id = _normalize_chain_id(chain_id)

    with alignment_mapping_tsv_path.open("r", encoding="utf-8") as handle:
        non_comment_lines = [
            line.rstrip("\n")
            for line in handle
            if line.strip() and not line.lstrip().startswith("#")
        ]

    if not non_comment_lines:
        raise ValueError(
            f"No tabular content found in alignment mapping TSV: {alignment_mapping_tsv_path}"
        )

    header = non_comment_lines[0].split("\t")
    if len(header) == 1:
        header = non_comment_lines[0].split()

    required_columns = {"pdb_residue_number", "relation"}
    missing_columns = required_columns.difference(header)
    if missing_columns:
        raise KeyError(
            f"Missing required columns in alignment mapping TSV {alignment_mapping_tsv_path}: "
            f"{sorted(missing_columns)}"
        )

    row_dict_list: list[dict[str, str]] = []
    for raw_line in non_comment_lines[1:]:
        fields = raw_line.split("\t")
        if len(fields) == 1:
            fields = raw_line.split()

        if len(fields) != len(header):
            raise ValueError(
                f"Column count mismatch in {alignment_mapping_tsv_path}. "
                f"Header has {len(header)} columns, row has {len(fields)} fields: {raw_line!r}"
            )

        row_dict_list.append(dict(zip(header, fields)))

    output_finalize_tsv_path.parent.mkdir(parents=True, exist_ok=True)

    fieldnames = [
        "chain_id",
        "model_seq_position",
        "final_resseq",
        "final_icode",
        "source_category",
        "original_resseq",
        "uniprot_resseq",
        "final_resname",
    ]

    model_seq_position = 0

    with output_finalize_tsv_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()

        for row in row_dict_list:
            relation = str(row.get("relation", "")).strip()
            if relation != "match":
                continue

            pdb_residue_number_text = str(row.get("pdb_residue_number", "")).strip()
            if not pdb_residue_number_text:
                continue

            try:
                final_resseq = int(pdb_residue_number_text)
            except ValueError as error:
                raise ValueError(
                    f"Invalid pdb_residue_number in {alignment_mapping_tsv_path}: "
                    f"{pdb_residue_number_text!r}"
                ) from error

            model_seq_position += 1

            writer.writerow(
                {
                    "chain_id": chain_id,
                    "model_seq_position": model_seq_position,
                    "final_resseq": final_resseq,
                    "final_icode": "",
                    "source_category": "original_pdb_match",
                    "original_resseq": final_resseq,
                    "uniprot_resseq": str(
                        row.get("uniprot_residue_number", "")
                    ).strip(),
                    "final_resname": str(row.get("pdb_residue", "")).strip(),
                }
            )

    return output_finalize_tsv_path


def finalize_protein_structure(
    pdb_id: str,
    protein_dir: str | Path,
    finalize_tsv_path: str | Path,
) -> dict[str, str | bool | int]:
    """
    Finalize the prepared protein structure using a dedicated finalize TSV.

    Standardized input:
        components/<PDBID>_protein_as_Amber.pdb

    Standardized output:
        components/<PDBID>_protein_final.pdb
    """
    protein_dir = Path(protein_dir)
    components_dir = protein_dir / "components"

    input_pdb_path = components_dir / f"{pdb_id}_protein_as_Amber.pdb"
    output_pdb_path = components_dir / f"{pdb_id}_protein_final.pdb"
    finalize_tsv_path = Path(finalize_tsv_path)

    numbering_map = load_finalize_numbering_map_from_tsv(finalize_tsv_path)

    result = renumber_structure_with_finalize_map(
        input_pdb_path=input_pdb_path,
        output_pdb_path=output_pdb_path,
        numbering_map=numbering_map,
    )

    return {
        "finalize_success": output_pdb_path.exists()
        and output_pdb_path.stat().st_size > 0,
        "finalize_input_path": str(input_pdb_path),
        "finalize_output_path": str(output_pdb_path),
        "finalize_tsv_path": str(finalize_tsv_path),
        "renumbered_residues": int(result["renumbered_residues"]),
        "renumbered_atoms": int(result["renumbered_atoms"]),
    }
