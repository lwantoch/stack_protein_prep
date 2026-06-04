"""Finalize the prepared protein structure by renumbering residues according to a
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
- Standard workflow input:
    components/<PDBID>_protein_as_Amber.pdb
- Standard workflow output:
    components/<PDBID>_protein_final.pdb

- Variant workflow input:
    components/<PDBID>_protein_as_Amber_<variant>.pdb
- Variant workflow output:
    components/<PDBID>_protein_final_<variant>.pdb

Matching rule
-------------
Residues are matched by:
    (chain_id, running residue index within chain, 1-based)

Only protein-like polymer residues are renumbered.

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
from typing import Iterable

from Bio.PDB import PDBIO, PDBParser

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

MODULE_NAME = "finalize_protein"

REQUIRED_FINALIZE_TSV_COLUMNS = {
    "chain_id",
    "model_seq_position",
    "final_resseq",
}

STANDARD_AMINO_ACID_RESNAMES = {
    "ALA",
    "ARG",
    "ASN",
    "ASP",
    "ASH",
    "CYS",
    "CYM",
    "CYX",
    "GLN",
    "GLU",
    "GLH",
    "GLY",
    "HIS",
    "HID",
    "HIE",
    "HIP",
    "ILE",
    "LEU",
    "LYS",
    "LYN",
    "MET",
    "PHE",
    "PRO",
    "SER",
    "THR",
    "TRP",
    "TYR",
    "VAL",
}

PROTEIN_LIKE_NONSTANDARD_RESNAMES = {
    "MSE",
    "SEC",
    "PYL",
    "SEP",
    "TPO",
    "PTR",
    "CSO",
    "CSD",
    "OCS",
    "KCX",
    "LLP",
    "MLY",
    "HYP",
}

AMBER_TERMINUS_RESNAMES = {
    "NALA",
    "NARG",
    "NASN",
    "NASP",
    "NASH",
    "NCYS",
    "NCYM",
    "NCYX",
    "NGLN",
    "NGLU",
    "NGLH",
    "NGLY",
    "NHIS",
    "NHID",
    "NHIE",
    "NHIP",
    "NILE",
    "NLEU",
    "NLYS",
    "NLYN",
    "NMET",
    "NPHE",
    "NPRO",
    "NSER",
    "NTHR",
    "NTRP",
    "NTYR",
    "NVAL",
    "CALA",
    "CARG",
    "CASN",
    "CASP",
    "CASH",
    "CCYS",
    "CCYM",
    "CCYX",
    "CGLN",
    "CGLU",
    "CGLH",
    "CGLY",
    "CHIS",
    "CHID",
    "CHIE",
    "CHIP",
    "CILE",
    "CLEU",
    "CLYS",
    "CLYN",
    "CMET",
    "CPHE",
    "CPRO",
    "CSER",
    "CTHR",
    "CTRP",
    "CTYR",
    "CVAL",
}

INTERNAL_CAP_RESNAMES = {"ACE", "NME"}

PROTEIN_POLYMER_RESNAMES = (
    STANDARD_AMINO_ACID_RESNAMES
    | PROTEIN_LIKE_NONSTANDARD_RESNAMES
    | AMBER_TERMINUS_RESNAMES
)


def _debug(message: str) -> None:
    print(f"[{MODULE_NAME}] {message}")


def sanitize_variant_label(variant_label: str) -> str:
    cleaned = "".join(
        character if character.isalnum() or character in {"_", "-"} else "_"
        for character in variant_label.strip()
    )
    cleaned = cleaned.strip("_")
    if not cleaned:
        raise ValueError(f"Invalid empty variant label derived from {variant_label!r}")
    return cleaned


def _normalize_chain_id(value: str) -> str:
    chain_id = str(value).strip()
    return chain_id if chain_id else "_"


def _normalize_icode(value: str) -> str:
    icode = str(value).strip()
    return icode if icode else " "


def _normalize_resname(value: str) -> str:
    return str(value).strip().upper()


def _infer_protein_data_dir_from_protein_dir(protein_dir: str | Path) -> Path:
    protein_dir = Path(protein_dir).resolve()
    return protein_dir.parent


def _build_module_log_path_from_protein_dir(
    *,
    protein_dir: str | Path,
    pdb_id: str,
    variant_label: str | None = None,
) -> Path:
    protein_data_dir = _infer_protein_data_dir_from_protein_dir(protein_dir)
    log_paths = build_module_log_paths(
        protein_data_dir=protein_data_dir,
        pdb_id=str(pdb_id).strip().upper(),
        module_name=MODULE_NAME,
        variant_label=variant_label,
    )
    return log_paths.module_log_path


def _build_module_log_path_from_input_path(
    *,
    input_pdb_path: str | Path,
    pdb_id: str,
    variant_label: str | None = None,
) -> Path:
    input_pdb_path = Path(input_pdb_path).resolve()
    protein_dir = input_pdb_path.parent.parent
    return _build_module_log_path_from_protein_dir(
        protein_dir=protein_dir,
        pdb_id=pdb_id,
        variant_label=variant_label,
    )


def _print_start_box(
    *,
    pdb_id: str,
    variant_label: str | None,
    input_pdb_path: Path,
    output_pdb_path: Path,
    finalize_tsv_path: Path,
    log_path: Path,
) -> None:
    print_mixed_box(
        title=f"{MODULE_NAME} · START",
        line_list=[
            f"PDB ID        : {pdb_id}",
            f"Variant       : {variant_label or '(standard)'}",
            f"Input PDB     : {input_pdb_path}",
            f"Output PDB    : {output_pdb_path}",
            f"Finalize TSV  : {finalize_tsv_path}",
            f"Log file      : {log_path}",
        ],
    )


def _print_summary_table(
    *,
    pdb_id: str,
    variant_label: str | None,
    renumbered_residues: int,
    renumbered_atoms: int,
    finalize_success: bool,
) -> None:
    print_light_table_box(
        title=f"{MODULE_NAME} · summary",
        row_list=[
            ("PDB ID", pdb_id),
            ("Variant", variant_label or "(standard)"),
            ("renumbered_residues", str(renumbered_residues)),
            ("renumbered_atoms", str(renumbered_atoms)),
            ("finalize_success", str(finalize_success)),
        ],
        left_header="Field",
        right_header="Value",
    )


def _print_end_box(
    *,
    pdb_id: str,
    variant_label: str | None,
    status: str,
    output_pdb_path: Path,
    finalize_tsv_path: Path,
    renumbered_residues: int,
    renumbered_atoms: int,
) -> None:
    print_double_header_box(
        title=f"{MODULE_NAME} · END",
        line_list=[
            f"PDB ID              : {pdb_id}",
            f"Variant             : {variant_label or '(standard)'}",
            f"Status              : {status}",
            f"Output PDB          : {output_pdb_path}",
            f"Finalize TSV        : {finalize_tsv_path}",
            f"Renumbered residues : {renumbered_residues}",
            f"Renumbered atoms    : {renumbered_atoms}",
        ],
    )


def build_standard_finalize_input_path(
    pdb_id: str,
    protein_dir: str | Path,
) -> Path:
    protein_dir = Path(protein_dir)
    return protein_dir / "components" / f"{pdb_id}_protein_as_Amber.pdb"


def build_standard_finalize_output_path(
    pdb_id: str,
    protein_dir: str | Path,
) -> Path:
    protein_dir = Path(protein_dir)
    return protein_dir / "components" / f"{pdb_id}_protein_final.pdb"


def build_variant_finalize_input_path(
    pdb_id: str,
    protein_dir: str | Path,
    variant_label: str,
) -> Path:
    protein_dir = Path(protein_dir)
    safe_variant_label = sanitize_variant_label(variant_label)
    return (
        protein_dir
        / "components"
        / f"{pdb_id}_protein_as_Amber_{safe_variant_label}.pdb"
    )


def build_variant_finalize_output_path(
    pdb_id: str,
    protein_dir: str | Path,
    variant_label: str,
) -> Path:
    protein_dir = Path(protein_dir)
    safe_variant_label = sanitize_variant_label(variant_label)
    return (
        protein_dir / "components" / f"{pdb_id}_protein_final_{safe_variant_label}.pdb"
    )


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


def _is_protein_like_polymer_residue(residue) -> bool:
    """
    Keep only true protein-like polymer residues for finalize renumbering.

    Excludes:
    - hetero residues
    - internal caps (ACE/NME)
    """
    if residue.id[0] != " ":
        return False

    resname = _normalize_resname(residue.get_resname())

    if resname in INTERNAL_CAP_RESNAMES:
        return False

    return resname in PROTEIN_POLYMER_RESNAMES


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

    Only protein-like polymer residues are renumbered.
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
                if not _is_protein_like_polymer_residue(residue):
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


def finalize_selected_structure(
    *,
    input_pdb_path: str | Path,
    output_pdb_path: str | Path,
    finalize_tsv_path: str | Path,
    variant_label: str | None = None,
) -> dict[str, str | bool | int]:
    """
    Finalize one explicitly selected structure using a dedicated finalize TSV.
    """
    input_pdb_path = Path(input_pdb_path).resolve()
    output_pdb_path = Path(output_pdb_path).resolve()
    finalize_tsv_path = Path(finalize_tsv_path).resolve()

    pdb_id = input_pdb_path.stem.split("_", maxsplit=1)[0].upper()
    module_log_path = _build_module_log_path_from_input_path(
        input_pdb_path=input_pdb_path,
        pdb_id=pdb_id,
        variant_label=variant_label,
    )

    write_log_header(
        log_path=module_log_path,
        module_name=MODULE_NAME,
        pdb_id=pdb_id,
        variant_label=variant_label,
        extra_lines=[
            f"input_pdb_path={input_pdb_path}",
            f"output_pdb_path={output_pdb_path}",
            f"finalize_tsv_path={finalize_tsv_path}",
        ],
    )

    _print_start_box(
        pdb_id=pdb_id,
        variant_label=variant_label,
        input_pdb_path=input_pdb_path,
        output_pdb_path=output_pdb_path,
        finalize_tsv_path=finalize_tsv_path,
        log_path=module_log_path,
    )

    _debug(f"Loading finalize TSV: {finalize_tsv_path}")
    append_log_text(
        module_log_path, f"[DEBUG] Loading finalize TSV: {finalize_tsv_path}"
    )

    try:
        numbering_map = load_finalize_numbering_map_from_tsv(finalize_tsv_path)

        append_key_value_block(
            log_path=module_log_path,
            title="INPUT SUMMARY",
            key_value_pairs=[
                ("pdb_id", pdb_id),
                ("variant", variant_label or "(standard)"),
                ("input_pdb_path", str(input_pdb_path)),
                ("output_pdb_path", str(output_pdb_path)),
                ("finalize_tsv_path", str(finalize_tsv_path)),
                ("numbering_map_entries", str(len(numbering_map))),
            ],
        )

        result = renumber_structure_with_finalize_map(
            input_pdb_path=input_pdb_path,
            output_pdb_path=output_pdb_path,
            numbering_map=numbering_map,
        )

        renumbered_residues = int(result["renumbered_residues"])
        renumbered_atoms = int(result["renumbered_atoms"])

        finalize_success = (
            output_pdb_path.exists() and output_pdb_path.stat().st_size > 0
        )

        output_dict: dict[str, str | bool | int] = {
            "finalize_success": finalize_success,
            "finalize_input_path": str(input_pdb_path),
            "finalize_output_path": str(output_pdb_path),
            "finalize_tsv_path": str(finalize_tsv_path),
            "renumbered_residues": renumbered_residues,
            "renumbered_atoms": renumbered_atoms,
        }

        if variant_label is not None:
            output_dict["finalize_variant"] = variant_label

        append_key_value_block(
            log_path=module_log_path,
            title="RESULT",
            key_value_pairs=[
                ("finalize_success", str(finalize_success)),
                ("renumbered_residues", str(renumbered_residues)),
                ("renumbered_atoms", str(renumbered_atoms)),
                ("output_pdb_path", str(output_pdb_path)),
            ],
        )

        _print_summary_table(
            pdb_id=pdb_id,
            variant_label=variant_label,
            renumbered_residues=renumbered_residues,
            renumbered_atoms=renumbered_atoms,
            finalize_success=finalize_success,
        )

        _print_end_box(
            pdb_id=pdb_id,
            variant_label=variant_label,
            status="success" if finalize_success else "failed",
            output_pdb_path=output_pdb_path,
            finalize_tsv_path=finalize_tsv_path,
            renumbered_residues=renumbered_residues,
            renumbered_atoms=renumbered_atoms,
        )

        return output_dict

    except Exception as exc:
        append_log_text(module_log_path, "[RESULT] FAILED")
        append_exception_traceback(module_log_path, exc)

        print_double_header_box(
            title=f"{MODULE_NAME} · END",
            line_list=[
                f"PDB ID        : {pdb_id}",
                f"Variant       : {variant_label or '(standard)'}",
                "Status        : failed",
                f"Error type    : {type(exc).__name__}",
                f"Error         : {exc}",
                f"Module log    : {module_log_path}",
            ],
        )
        raise


def finalize_variant_structure(
    pdb_id: str,
    protein_dir: str | Path,
    *,
    variant_label: str,
    finalize_tsv_path: str | Path,
    input_pdb_path: str | Path | None = None,
) -> dict[str, str | bool | int]:
    """
    Finalize one named variant using a variant-specific input/output path.

    Default input:
        components/<PDBID>_protein_as_Amber_<variant>.pdb

    Output:
        components/<PDBID>_protein_final_<variant>.pdb
    """
    protein_dir = Path(protein_dir)

    if input_pdb_path is None:
        input_pdb_path = build_variant_finalize_input_path(
            pdb_id=pdb_id,
            protein_dir=protein_dir,
            variant_label=variant_label,
        )

    output_pdb_path = build_variant_finalize_output_path(
        pdb_id=pdb_id,
        protein_dir=protein_dir,
        variant_label=variant_label,
    )

    return finalize_selected_structure(
        input_pdb_path=input_pdb_path,
        output_pdb_path=output_pdb_path,
        finalize_tsv_path=finalize_tsv_path,
        variant_label=variant_label,
    )


def finalize_protein_structure(
    pdb_id: str,
    protein_dir: str | Path,
    finalize_tsv_path: str | Path,
) -> dict[str, str | bool | int]:
    """
    Finalize the standard protein structure using a dedicated finalize TSV.

    Standardized input:
        components/<PDBID>_protein_as_Amber.pdb

    Standardized output:
        components/<PDBID>_protein_final.pdb
    """
    protein_dir = Path(protein_dir)

    input_pdb_path = build_standard_finalize_input_path(
        pdb_id=pdb_id,
        protein_dir=protein_dir,
    )
    output_pdb_path = build_standard_finalize_output_path(
        pdb_id=pdb_id,
        protein_dir=protein_dir,
    )

    return finalize_selected_structure(
        input_pdb_path=input_pdb_path,
        output_pdb_path=output_pdb_path,
        finalize_tsv_path=finalize_tsv_path,
        variant_label=None,
    )
