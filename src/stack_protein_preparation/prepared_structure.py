# /home/grheco/repositorios/stack_protein_prep/src/stack_protein_preparation/prepared_structure.py

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Literal

StructureVariant = Literal[
    "single",
    "gaps",
    "complete",
    "best_complete",
    "large_gap_complete",
]


@dataclass(slots=True)
class PreparedStructureSummary:
    pdb_id: str
    output_pdb_path: Path
    protein_input_path: Path
    protein_input_paths: tuple[Path, ...]
    water_input_path: Path | None
    ligand_input_path: Path | None
    metals_input_path: Path | None
    structure_variant: str
    protein_inputs_were_fragments: bool
    n_protein_inputs_merged: int
    water_included: bool
    ligand_included: bool
    metals_included: bool
    n_atom_records_written: int


def sanitize_variant_label(variant_label: str) -> str:
    """
    Convert a variant label into a filesystem-safe token.
    """
    cleaned = "".join(
        character if character.isalnum() or character in {"_", "-"} else "_"
        for character in variant_label.strip()
    )
    cleaned = cleaned.strip("_")
    if not cleaned:
        raise ValueError(f"Invalid empty variant label derived from {variant_label!r}")
    return cleaned


def get_default_prepared_protein_input_path(
    pdb_directory: str | Path,
    pdb_id: str,
) -> Path:
    """
    Backward-compatible default protein input path.

    Notes
    -----
    In the new pipeline, prepared_structure should preferably receive an explicit
    protein_input_path or protein_input_paths from the previous step instead of
    relying on this default.
    """
    pdb_directory = Path(pdb_directory)
    return pdb_directory / "components" / f"{pdb_id}_protein_internal_capped.pdb"


def get_default_water_input_path(
    pdb_directory: str | Path,
    pdb_id: str,
) -> Path:
    pdb_directory = Path(pdb_directory)
    return pdb_directory / "components" / f"{pdb_id}_water.pdb"


def get_default_ligand_input_path(
    pdb_directory: str | Path,
    pdb_id: str,
) -> Path:
    pdb_directory = Path(pdb_directory)
    return pdb_directory / "components" / f"{pdb_id}_ligand.pdb"


def get_default_metals_input_path(
    pdb_directory: str | Path,
    pdb_id: str,
) -> Path:
    pdb_directory = Path(pdb_directory)
    return pdb_directory / "components" / f"{pdb_id}_metals.pdb"


def get_prepared_structure_output_path(
    pdb_directory: str | Path,
    pdb_id: str,
    structure_variant: str = "single",
) -> Path:
    """
    Return the prepared structure output path for one variant.

    Convention
    ----------
    prepared/<variant>/<PDBID>.pdb
    """
    pdb_directory = Path(pdb_directory)
    safe_variant_label = sanitize_variant_label(structure_variant)
    return pdb_directory / "prepared" / safe_variant_label / f"{pdb_id}.pdb"


def _is_atom_or_hetatm_record(line: str) -> bool:
    return line.startswith("ATOM  ") or line.startswith("HETATM")


def _read_backbone_nonstd_hetatm(
    protein_component_path: Path,
    exclude_resnames: set[str],
) -> list[str]:
    """Return HETATM lines from the original protein component that are
    backbone non-standard residues (e.g. PTR, SEP) — i.e. HETATM records
    that are NOT water, metals, or explicit ligand residues.

    pdb2gmx silently strips these during protonation; this function recovers
    them so they appear in the final prepared structure.
    """
    _WATER_NAMES = {"HOH", "WAT", "SOL"}
    if not protein_component_path.is_file():
        return []
    result: list[str] = []
    for raw_line in protein_component_path.open("r", encoding="utf-8", errors="replace"):
        line = raw_line.rstrip("\n")
        if not line.startswith("HETATM"):
            continue
        resname = line[17:20].strip().upper() if len(line) >= 20 else ""
        if resname in _WATER_NAMES or resname in exclude_resnames:
            continue
        result.append(line)
    return result


def _read_atom_lines_from_pdb(pdb_path: str | Path) -> list[str]:
    pdb_path = Path(pdb_path)

    if not pdb_path.exists():
        raise FileNotFoundError(f"PDB file not found: {pdb_path}")

    atom_lines: list[str] = []

    with pdb_path.open("r", encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.rstrip("\n")
            if _is_atom_or_hetatm_record(line):
                atom_lines.append(line)

    return atom_lines


_CHAIN_ID_ALPHABET = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789"


def _renumber_atom_serial(line: str, atom_serial: int) -> str:
    return f"{line[:6]}{atom_serial:>5}{line[11:]}"


def _rewrite_chain_id(line: str, chain_id: str) -> str:
    """Replace the chain ID (column 22, 0-indexed) in an ATOM/HETATM line."""
    return f"{line[:21]}{chain_id[:1]}{line[22:]}"


def _normalize_protein_input_paths(
    protein_input_path: str | Path | None,
    protein_input_paths: list[str | Path] | tuple[str | Path, ...] | None,
) -> list[Path]:
    """
    Normalize the protein input specification.

    Exactly one of the following modes is allowed:
    - single protein_input_path
    - multiple protein_input_paths

    Returns
    -------
    list[Path]
        Ordered list of protein PDB paths to merge.
    """
    has_single = protein_input_path is not None
    has_multiple = protein_input_paths is not None and len(protein_input_paths) > 0

    if has_single and has_multiple:
        raise ValueError(
            "Provide either protein_input_path or protein_input_paths, not both."
        )

    if not has_single and not has_multiple:
        raise ValueError(
            "No protein input specified. Provide protein_input_path or protein_input_paths."
        )

    if has_single:
        return [Path(protein_input_path)]

    assert protein_input_paths is not None

    normalized_paths: list[Path] = []
    for path_like in protein_input_paths:
        normalized_path = Path(path_like)
        normalized_paths.append(normalized_path)

    if not normalized_paths:
        raise ValueError("protein_input_paths was provided but is empty.")

    return normalized_paths


def _write_merged_pdb_sections(
    output_pdb_path: str | Path,
    ordered_sections: list[list[str]],
    section_chain_ids: list[str | None] | None = None,
) -> int:
    """
    Write a merged PDB from already prepared sections.

    Each non-empty section is written as one block followed by TER.
    When *section_chain_ids* is provided it must have the same length as
    *ordered_sections*; each non-None entry rewrites the chain-ID column of
    every ATOM/HETATM line in the corresponding section.  This is required
    when multiple protein fragments share the same original chain ID: pdb2gmx
    (and BioPython) need distinct chain IDs to recognise each segment as an
    independent chain with its own N- and C-termini.
    """
    output_pdb_path = Path(output_pdb_path)
    output_pdb_path.parent.mkdir(parents=True, exist_ok=True)

    if section_chain_ids is None:
        section_chain_ids = [None] * len(ordered_sections)

    atom_serial = 0
    n_written = 0

    with output_pdb_path.open("w", encoding="utf-8") as handle:
        for section_lines, chain_id in zip(ordered_sections, section_chain_ids):
            if not section_lines:
                continue

            for line in section_lines:
                atom_serial += 1
                if chain_id is not None and _is_atom_or_hetatm_record(line):
                    line = _rewrite_chain_id(line, chain_id)
                handle.write(_renumber_atom_serial(line, atom_serial) + "\n")
                n_written += 1

            handle.write("TER\n")

        handle.write("END\n")

    return n_written


def build_prepared_structure(
    output_pdb_path: str | Path,
    protein_input_path: str | Path | None = None,
    protein_input_paths: list[str | Path] | tuple[str | Path, ...] | None = None,
    water_input_path: str | Path | None = None,
    ligand_input_path: str | Path | None = None,
    metals_input_path: str | Path | None = None,
    backbone_nonstd_input_path: str | Path | None = None,
    *,
    structure_variant: str = "single",
) -> PreparedStructureSummary:
    """
    Build one final prepared structure.

    Supported protein input modes
    -----------------------------
    1. Single protein input:
        protein_input_path=/path/to/protein.pdb

    2. Multiple fragment inputs:
        protein_input_paths=[
            /path/to/fragment_01H.pdb,
            /path/to/fragment_02H.pdb,
            ...
        ]

    In fragment mode, the fragments are merged in the exact order provided.
    Each fragment is written as a separate section followed by TER.
    """
    output_pdb_path = Path(output_pdb_path)

    resolved_protein_input_path_list = _normalize_protein_input_paths(
        protein_input_path=protein_input_path,
        protein_input_paths=protein_input_paths,
    )

    protein_sections: list[list[str]] = []
    for resolved_protein_input_path in resolved_protein_input_path_list:
        protein_atom_lines = _read_atom_lines_from_pdb(resolved_protein_input_path)

        if not protein_atom_lines:
            raise ValueError(
                f"No ATOM/HETATM records found in protein file: {resolved_protein_input_path}"
            )

        protein_sections.append(protein_atom_lines)

    water_input = Path(water_input_path) if water_input_path is not None else None
    ligand_input = Path(ligand_input_path) if ligand_input_path is not None else None
    metals_input = Path(metals_input_path) if metals_input_path is not None else None

    water_atom_lines = (
        _read_atom_lines_from_pdb(water_input)
        if water_input is not None and water_input.exists()
        else None
    )
    ligand_atom_lines = (
        _read_atom_lines_from_pdb(ligand_input)
        if ligand_input is not None and ligand_input.exists()
        else None
    )
    metals_atom_lines = (
        _read_atom_lines_from_pdb(metals_input)
        if metals_input is not None and metals_input.exists()
        else None
    )

    backbone_nonstd_lines: list[str] = []
    if backbone_nonstd_input_path is not None:
        _bnp = Path(backbone_nonstd_input_path)
        _exclude = set()
        if ligand_atom_lines:
            for _l in ligand_atom_lines:
                if len(_l) >= 20:
                    _exclude.add(_l[17:20].strip().upper())
        if metals_atom_lines:
            for _l in metals_atom_lines:
                if len(_l) >= 20:
                    _exclude.add(_l[17:20].strip().upper())
        backbone_nonstd_lines = _read_backbone_nonstd_hetatm(_bnp, _exclude)

    # When there are multiple protein fragments, assign each a unique chain ID so
    # that pdb2gmx (and BioPython) recognise each segment as an independent chain
    # with its own N- and C-termini.  A single-fragment input keeps its original
    # chain ID (None → no rewrite).
    if len(protein_sections) > 1:
        protein_chain_ids: list[str | None] = [
            _CHAIN_ID_ALPHABET[i % len(_CHAIN_ID_ALPHABET)]
            for i in range(len(protein_sections))
        ]
    else:
        protein_chain_ids = [None]

    ordered_sections: list[list[str]] = []
    if protein_sections and backbone_nonstd_lines:
        # Insert backbone non-standard residues (e.g. PTR) into the first
        # protein section sorted by (chain, resnum) so they appear at the
        # correct sequence position.
        def _sort_key(l: str) -> tuple[str, int]:
            try:
                return (l[21], int(l[22:26]))
            except (IndexError, ValueError):
                return ("", 0)
        merged = list(protein_sections[0]) + backbone_nonstd_lines
        merged.sort(key=_sort_key)
        ordered_sections.append(merged)
        ordered_sections.extend(protein_sections[1:])
    else:
        ordered_sections.extend(protein_sections)
    section_chain_ids: list[str | None] = list(protein_chain_ids)

    if water_atom_lines:
        ordered_sections.append(water_atom_lines)
        section_chain_ids.append(None)
    if ligand_atom_lines:
        ordered_sections.append(ligand_atom_lines)
        section_chain_ids.append(None)
    if metals_atom_lines:
        ordered_sections.append(metals_atom_lines)
        section_chain_ids.append(None)

    n_atom_records_written = _write_merged_pdb_sections(
        output_pdb_path=output_pdb_path,
        ordered_sections=ordered_sections,
        section_chain_ids=section_chain_ids,
    )

    return PreparedStructureSummary(
        pdb_id=output_pdb_path.stem,
        output_pdb_path=output_pdb_path,
        protein_input_path=resolved_protein_input_path_list[0],
        protein_input_paths=tuple(resolved_protein_input_path_list),
        water_input_path=water_input if water_atom_lines else None,
        ligand_input_path=ligand_input if ligand_atom_lines else None,
        metals_input_path=metals_input if metals_atom_lines else None,
        structure_variant=structure_variant,
        protein_inputs_were_fragments=len(resolved_protein_input_path_list) > 1,
        n_protein_inputs_merged=len(resolved_protein_input_path_list),
        water_included=bool(water_atom_lines),
        ligand_included=bool(ligand_atom_lines),
        metals_included=bool(metals_atom_lines),
        n_atom_records_written=n_atom_records_written,
    )


def build_prepared_structure_for_variant(
    pdb_directory: str | Path,
    pdb_id: str,
    *,
    structure_variant: str,
    protein_input_path: str | Path | None = None,
    protein_input_paths: list[str | Path] | tuple[str | Path, ...] | None = None,
    water_input_path: str | Path | None = None,
    ligand_input_path: str | Path | None = None,
    metals_input_path: str | Path | None = None,
    backbone_nonstd_input_path: str | Path | None = None,
) -> PreparedStructureSummary:
    pdb_directory = Path(pdb_directory)

    output_pdb_path = get_prepared_structure_output_path(
        pdb_directory=pdb_directory,
        pdb_id=pdb_id,
        structure_variant=structure_variant,
    )

    resolved_water_input_path = (
        Path(water_input_path)
        if water_input_path is not None
        else get_default_water_input_path(
            pdb_directory=pdb_directory,
            pdb_id=pdb_id,
        )
    )

    resolved_ligand_input_path = (
        Path(ligand_input_path)
        if ligand_input_path is not None
        else get_default_ligand_input_path(
            pdb_directory=pdb_directory,
            pdb_id=pdb_id,
        )
    )

    resolved_metals_input_path = (
        Path(metals_input_path)
        if metals_input_path is not None
        else get_default_metals_input_path(
            pdb_directory=pdb_directory,
            pdb_id=pdb_id,
        )
    )

    if protein_input_path is None and protein_input_paths is None:
        raise ValueError(
            "build_prepared_structure_for_variant requires either "
            "protein_input_path or protein_input_paths."
        )

    return build_prepared_structure(
        output_pdb_path=output_pdb_path,
        protein_input_path=protein_input_path,
        protein_input_paths=protein_input_paths,
        water_input_path=resolved_water_input_path,
        ligand_input_path=resolved_ligand_input_path,
        metals_input_path=resolved_metals_input_path,
        backbone_nonstd_input_path=backbone_nonstd_input_path,
        structure_variant=structure_variant,
    )


def build_prepared_structure_for_pdb_directory(
    pdb_directory: str | Path,
    pdb_id: str,
    had_gaps: bool | None = None,
    structure_variant: str | None = None,
    protein_input_path: str | Path | None = None,
    protein_input_paths: list[str | Path] | tuple[str | Path, ...] | None = None,
    water_input_path: str | Path | None = None,
    ligand_input_path: str | Path | None = None,
    metals_input_path: str | Path | None = None,
) -> PreparedStructureSummary:
    """
    Backward-compatible wrapper.

    Legacy mapping
    --------------
    - had_gaps=False -> variant='single'
    - had_gaps=True, structure_variant='gaps' -> variant='gaps'
    - had_gaps=True, structure_variant='complete' -> variant='complete'

    New code should prefer:
        build_prepared_structure_for_variant(...)
    """
    pdb_directory = Path(pdb_directory)

    if protein_input_path is None and protein_input_paths is None:
        resolved_protein_input_path = get_default_prepared_protein_input_path(
            pdb_directory=pdb_directory,
            pdb_id=pdb_id,
        )
        resolved_protein_input_paths = None
    else:
        resolved_protein_input_path = (
            Path(protein_input_path) if protein_input_path is not None else None
        )
        resolved_protein_input_paths = (
            [Path(path) for path in protein_input_paths]
            if protein_input_paths is not None
            else None
        )

    if structure_variant is None:
        if had_gaps is None or had_gaps is False:
            resolved_variant = "single"
        else:
            raise ValueError(
                "When had_gaps=True in backward-compatible mode, "
                "structure_variant must be provided."
            )
    else:
        resolved_variant = structure_variant

    return build_prepared_structure_for_variant(
        pdb_directory=pdb_directory,
        pdb_id=pdb_id,
        structure_variant=resolved_variant,
        protein_input_path=resolved_protein_input_path,
        protein_input_paths=resolved_protein_input_paths,
        water_input_path=water_input_path,
        ligand_input_path=ligand_input_path,
        metals_input_path=metals_input_path,
    )


def prepared_structure_summary_to_dict(
    summary: PreparedStructureSummary,
) -> dict[str, str | bool | int]:
    return {
        "prepared_structure_success": summary.output_pdb_path.is_file()
        and summary.output_pdb_path.stat().st_size > 0,
        "prepared_structure_output_path": str(summary.output_pdb_path),
        "prepared_structure_protein_input_path": str(summary.protein_input_path),
        "prepared_structure_protein_input_paths": "|".join(
            str(path) for path in summary.protein_input_paths
        ),
        "prepared_structure_n_protein_inputs_merged": summary.n_protein_inputs_merged,
        "prepared_structure_protein_inputs_were_fragments": (
            summary.protein_inputs_were_fragments
        ),
        "prepared_structure_water_input_path": (
            str(summary.water_input_path)
            if summary.water_input_path is not None
            else ""
        ),
        "prepared_structure_ligand_input_path": (
            str(summary.ligand_input_path)
            if summary.ligand_input_path is not None
            else ""
        ),
        "prepared_structure_metals_input_path": (
            str(summary.metals_input_path)
            if summary.metals_input_path is not None
            else ""
        ),
        "prepared_structure_variant": summary.structure_variant,
        "prepared_structure_water_included": summary.water_included,
        "prepared_structure_ligand_included": summary.ligand_included,
        "prepared_structure_metals_included": summary.metals_included,
        "prepared_structure_n_atom_records_written": summary.n_atom_records_written,
    }


def summarize_prepared_structure(summary: PreparedStructureSummary) -> str:
    return (
        f"pdb_id={summary.pdb_id}\n"
        f"output_pdb_path={summary.output_pdb_path}\n"
        f"protein_input_path={summary.protein_input_path}\n"
        f"protein_input_paths={[str(path) for path in summary.protein_input_paths]}\n"
        f"water_input_path={summary.water_input_path}\n"
        f"ligand_input_path={summary.ligand_input_path}\n"
        f"metals_input_path={summary.metals_input_path}\n"
        f"structure_variant={summary.structure_variant}\n"
        f"protein_inputs_were_fragments={summary.protein_inputs_were_fragments}\n"
        f"n_protein_inputs_merged={summary.n_protein_inputs_merged}\n"
        f"water_included={summary.water_included}\n"
        f"ligand_included={summary.ligand_included}\n"
        f"metals_included={summary.metals_included}\n"
        f"n_atom_records_written={summary.n_atom_records_written}"
    )


if __name__ == "__main__":
    import argparse

    argument_parser = argparse.ArgumentParser(
        description="Build the final prepared structure for one PDB directory."
    )
    argument_parser.add_argument(
        "pdb_directory",
        type=Path,
        help="Protein directory, for example data/proteins/1ABC",
    )
    argument_parser.add_argument(
        "pdb_id",
        type=str,
        help="PDB ID, for example 1ABC",
    )
    argument_parser.add_argument(
        "--structure-variant",
        type=str,
        default="single",
        help="Variant label, e.g. single, gaps, best_complete, large_gap_complete",
    )
    argument_parser.add_argument(
        "--protein-input",
        type=Path,
        default=None,
        help="Optional explicit single protein input path.",
    )
    argument_parser.add_argument(
        "--protein-inputs",
        type=Path,
        nargs="+",
        default=None,
        help=(
            "Optional explicit ordered list of protein fragment inputs. "
            "Use this for fragment-wise protonated gaps structures."
        ),
    )
    argument_parser.add_argument(
        "--water-input",
        type=Path,
        default=None,
        help="Optional explicit water input path.",
    )
    argument_parser.add_argument(
        "--ligand-input",
        type=Path,
        default=None,
        help="Optional explicit ligand input path.",
    )
    argument_parser.add_argument(
        "--metals-input",
        type=Path,
        default=None,
        help="Optional explicit metals input path.",
    )

    arguments = argument_parser.parse_args()

    if arguments.protein_input is None and arguments.protein_inputs is None:
        protein_input_path = get_default_prepared_protein_input_path(
            pdb_directory=arguments.pdb_directory,
            pdb_id=arguments.pdb_id,
        )
        protein_input_paths = None
    else:
        protein_input_path = arguments.protein_input
        protein_input_paths = arguments.protein_inputs

    summary = build_prepared_structure_for_variant(
        pdb_directory=arguments.pdb_directory,
        pdb_id=arguments.pdb_id,
        structure_variant=arguments.structure_variant,
        protein_input_path=protein_input_path,
        protein_input_paths=protein_input_paths,
        water_input_path=arguments.water_input,
        ligand_input_path=arguments.ligand_input,
        metals_input_path=arguments.metals_input,
    )

    print(summarize_prepared_structure(summary))
