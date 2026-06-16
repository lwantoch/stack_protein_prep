"""Data types and summary dataclass for prepared protein structure variants."""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Literal


@dataclass(slots=True)
class PreparedStructureSummary:
    pdb_id: str
    output_pdb_path: Path
    protein_input_path: Path
    water_input_path: Path | None
    ligand_input_path: Path | None
    metals_input_path: Path | None
    had_gaps: bool
    structure_variant: str | None
    water_included: bool
    ligand_included: bool
    metals_included: bool
    n_atom_records_written: int


def sanitize_variant_label(variant_label: str) -> str:
    """Convert a variant label into a filesystem-safe token."""
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
    """Return the default protein input path for a PDB directory."""
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
    *,
    had_gaps: bool,
    structure_variant: str | None = None,
) -> Path:
    """Return the prepared structure output path.

    Convention
    ----------
    - had_gaps=False:
        prepared/<PDBID>.pdb
    - had_gaps=True, structure_variant='gaps':
        prepared/gaps/<PDBID>.pdb
    - had_gaps=True, structure_variant='complete':
        prepared/complete/<PDBID>.pdb
    - had_gaps=True, structure_variant=None:
        raises ValueError
    """
    pdb_directory = Path(pdb_directory)

    if not had_gaps:
        return pdb_directory / "prepared" / f"{pdb_id}.pdb"

    if structure_variant is None:
        raise ValueError(
            "structure_variant must be provided when had_gaps=True."
        )

    safe_variant = sanitize_variant_label(structure_variant)
    return pdb_directory / "prepared" / safe_variant / f"{pdb_id}.pdb"


def _is_atom_or_hetatm_record(line: str) -> bool:
    return line.startswith("ATOM  ") or line.startswith("HETATM")


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


def _renumber_atom_serial(line: str, atom_serial: int) -> str:
    return f"{line[:6]}{atom_serial:>5}{line[11:]}"


def _rewrite_chain_id(line: str, chain_id: str) -> str:
    """Replace the chain ID (column 22, 0-indexed) in an ATOM/HETATM line."""
    return f"{line[:21]}{chain_id[:1]}{line[22:]}"


def _normalize_intermediate_fragment_termini(atom_lines: list[str]) -> list[str]:
    """Rename OC1→O and strip OC2 from a protonated non-terminal fragment.

    pdb2gmx adds OC1/OC2 (carboxylate oxygens) at the C-terminus of each
    independently protonated fragment.  When fragments are assembled into one
    chain, these atoms appear at internal positions where downstream pdb2gmx
    rejects them.  Renaming OC1 to the standard backbone carbonyl O and
    dropping OC2 makes each internal residue look like a normal residue.
    """
    result = []
    for line in atom_lines:
        if _is_atom_or_hetatm_record(line):
            atom_name = line[12:16].strip()
            if atom_name == "OC2":
                continue
            if atom_name == "OC1":
                line = line[:12] + " O  " + line[16:]
        result.append(line)
    return result


def _write_merged_pdb_sections(
    output_pdb_path: str | Path,
    ordered_sections: list[list[str]],
    section_chain_ids: list[str | None] | None = None,
) -> int:
    """Write a merged PDB from already prepared sections."""
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
    protein_input_paths: list[str | Path] | None = None,
    water_input_path: str | Path | None = None,
    ligand_input_path: str | Path | None = None,
    metals_input_path: str | Path | None = None,
    *,
    had_gaps: bool = False,
    structure_variant: str | None = None,
) -> PreparedStructureSummary:
    """Build one final prepared structure."""
    output_pdb_path = Path(output_pdb_path)

    if protein_input_paths is not None:
        resolved_protein = Path(protein_input_paths[0])
        protein_sections: list[list[str]] = [
            _read_atom_lines_from_pdb(Path(p)) for p in protein_input_paths
        ]
        if len(protein_sections) > 1:
            protein_sections = [
                _normalize_intermediate_fragment_termini(s) if i < len(protein_sections) - 1 else s
                for i, s in enumerate(protein_sections)
            ]
    elif protein_input_path is not None:
        resolved_protein = Path(protein_input_path)
        protein_sections = [_read_atom_lines_from_pdb(resolved_protein)]
    else:
        raise ValueError("Either protein_input_path or protein_input_paths is required.")

    if not any(protein_sections):
        raise ValueError(
            f"No ATOM/HETATM records found in protein file: {resolved_protein}"
        )

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

    ordered_sections: list[list[str]] = list(protein_sections)
    if water_atom_lines:
        ordered_sections.append(water_atom_lines)
    if ligand_atom_lines:
        ordered_sections.append(ligand_atom_lines)
    if metals_atom_lines:
        ordered_sections.append(metals_atom_lines)

    n_atom_records_written = _write_merged_pdb_sections(
        output_pdb_path=output_pdb_path,
        ordered_sections=ordered_sections,
    )

    return PreparedStructureSummary(
        pdb_id=output_pdb_path.stem,
        output_pdb_path=output_pdb_path,
        protein_input_path=resolved_protein,
        water_input_path=water_input if water_atom_lines else None,
        ligand_input_path=ligand_input if ligand_atom_lines else None,
        metals_input_path=metals_input if metals_atom_lines else None,
        had_gaps=had_gaps,
        structure_variant=structure_variant,
        water_included=bool(water_atom_lines),
        ligand_included=bool(ligand_atom_lines),
        metals_included=bool(metals_atom_lines),
        n_atom_records_written=n_atom_records_written,
    )


def build_prepared_structure_for_pdb_directory(
    pdb_directory: str | Path,
    pdb_id: str,
    had_gaps: bool = False,
    structure_variant: str | None = None,
    protein_input_path: str | Path | None = None,
    water_input_path: str | Path | None = None,
    ligand_input_path: str | Path | None = None,
    metals_input_path: str | Path | None = None,
) -> PreparedStructureSummary:
    """Build the final prepared structure for one PDB directory.

    Parameters
    ----------
    pdb_directory
        Root PDB directory, e.g. data/proteins/1ABC.
    pdb_id
        PDB ID string.
    had_gaps
        Whether the structure had gaps that required filling.
    structure_variant
        Variant label ('gaps', 'complete', etc.). Required when had_gaps=True.
    protein_input_path
        Explicit protein input path. If None, the default capped protein path is used.
    water_input_path, ligand_input_path, metals_input_path
        Explicit component paths. If None, default component paths are used.
    """
    pdb_directory = Path(pdb_directory)

    output_pdb_path = get_prepared_structure_output_path(
        pdb_directory=pdb_directory,
        pdb_id=pdb_id,
        had_gaps=had_gaps,
        structure_variant=structure_variant,
    )

    resolved_protein = (
        Path(protein_input_path)
        if protein_input_path is not None
        else get_default_prepared_protein_input_path(pdb_directory, pdb_id)
    )

    resolved_water = (
        Path(water_input_path)
        if water_input_path is not None
        else get_default_water_input_path(pdb_directory, pdb_id)
    )

    resolved_ligand = (
        Path(ligand_input_path)
        if ligand_input_path is not None
        else get_default_ligand_input_path(pdb_directory, pdb_id)
    )

    resolved_metals = (
        Path(metals_input_path)
        if metals_input_path is not None
        else get_default_metals_input_path(pdb_directory, pdb_id)
    )

    return build_prepared_structure(
        output_pdb_path=output_pdb_path,
        protein_input_path=resolved_protein,
        water_input_path=resolved_water,
        ligand_input_path=resolved_ligand,
        metals_input_path=resolved_metals,
        had_gaps=had_gaps,
        structure_variant=structure_variant,
    )


def build_prepared_structure_for_variant(
    pdb_directory: str | Path,
    pdb_id: str,
    structure_variant: str,
    protein_input_path: str | Path | None = None,
    protein_input_paths: list[str | Path] | None = None,
    water_input_path: str | Path | None = None,
    ligand_input_path: str | Path | None = None,
    metals_input_path: str | Path | None = None,
    backbone_nonstd_input_path: str | Path | None = None,
) -> PreparedStructureSummary:
    """Build the final prepared structure for one variant.

    Always writes to prepared/<variant>/<pdb_id>.pdb.
    Accepts either a single protein_input_path or a list (protein_input_paths);
    backbone_nonstd_input_path atoms are appended to the protein section.
    """
    pdb_directory = Path(pdb_directory)

    output_pdb_path = get_prepared_structure_output_path(
        pdb_directory=pdb_directory,
        pdb_id=pdb_id,
        had_gaps=True,
        structure_variant=structure_variant,
    )

    if protein_input_paths is not None:
        primary_protein_path = Path(protein_input_paths[0])
        protein_sections: list[list[str]] = [
            _read_atom_lines_from_pdb(Path(p)) for p in protein_input_paths
        ]
        if len(protein_sections) > 1:
            protein_sections = [
                _normalize_intermediate_fragment_termini(s) if i < len(protein_sections) - 1 else s
                for i, s in enumerate(protein_sections)
            ]
    elif protein_input_path is not None:
        primary_protein_path = Path(protein_input_path)
        protein_sections = [_read_atom_lines_from_pdb(primary_protein_path)]
    else:
        raise ValueError("Either protein_input_path or protein_input_paths must be provided.")

    if not any(protein_sections):
        raise ValueError(f"No ATOM/HETATM records found in protein input(s) for {pdb_id}")

    if backbone_nonstd_input_path is not None:
        backbone_path = Path(backbone_nonstd_input_path)
        if backbone_path.exists():
            protein_sections.append(_read_atom_lines_from_pdb(backbone_path))

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

    ordered_sections: list[list[str]] = list(protein_sections)
    if water_atom_lines:
        ordered_sections.append(water_atom_lines)
    if ligand_atom_lines:
        ordered_sections.append(ligand_atom_lines)
    if metals_atom_lines:
        ordered_sections.append(metals_atom_lines)

    n_atom_records_written = _write_merged_pdb_sections(
        output_pdb_path=output_pdb_path,
        ordered_sections=ordered_sections,
    )

    return PreparedStructureSummary(
        pdb_id=pdb_id,
        output_pdb_path=output_pdb_path,
        protein_input_path=primary_protein_path,
        water_input_path=water_input if water_atom_lines else None,
        ligand_input_path=ligand_input if ligand_atom_lines else None,
        metals_input_path=metals_input if metals_atom_lines else None,
        had_gaps=True,
        structure_variant=structure_variant,
        water_included=bool(water_atom_lines),
        ligand_included=bool(ligand_atom_lines),
        metals_included=bool(metals_atom_lines),
        n_atom_records_written=n_atom_records_written,
    )


def prepared_structure_summary_to_dict(
    summary: PreparedStructureSummary,
) -> dict[str, str | bool | int]:
    return {
        "prepared_structure_success": summary.output_pdb_path.is_file()
        and summary.output_pdb_path.stat().st_size > 0,
        "prepared_structure_output_path": str(summary.output_pdb_path),
        "prepared_structure_protein_input_path": str(summary.protein_input_path),
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
        "prepared_structure_had_gaps": summary.had_gaps,
        "prepared_structure_variant": summary.structure_variant or "",
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
        f"water_input_path={summary.water_input_path}\n"
        f"ligand_input_path={summary.ligand_input_path}\n"
        f"metals_input_path={summary.metals_input_path}\n"
        f"had_gaps={summary.had_gaps}\n"
        f"structure_variant={summary.structure_variant}\n"
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
        "--had-gaps",
        action="store_true",
        help="Indicate that gaps were present.",
    )
    argument_parser.add_argument(
        "--structure-variant",
        type=str,
        default=None,
        help="Variant label, e.g. gaps, complete.",
    )
    argument_parser.add_argument(
        "--protein-input",
        type=Path,
        default=None,
        help="Optional explicit single protein input path.",
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

    summary = build_prepared_structure_for_pdb_directory(
        pdb_directory=arguments.pdb_directory,
        pdb_id=arguments.pdb_id,
        had_gaps=arguments.had_gaps,
        structure_variant=arguments.structure_variant,
        protein_input_path=arguments.protein_input,
        water_input_path=arguments.water_input,
        ligand_input_path=arguments.ligand_input,
        metals_input_path=arguments.metals_input,
    )

    print(summarize_prepared_structure(summary))
