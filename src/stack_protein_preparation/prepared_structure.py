"""
/home/grheco/repositorios/stack_protein_prep/src/stack_protein_preparation/prepared_structure.py

Build the final prepared structure for one PDB directory.

Purpose
-------
- create the prepared output directory
- decide the correct final output path:
    - prepared/<PDBID>.pdb
    - prepared/gaps/<PDBID>.pdb
    - prepared/complete/<PDBID>.pdb
- merge the already prepared protein with optional crystal waters, ligands, and metals
- write the final combined PDB named exactly <PDBID>.pdb

Expected protein input
----------------------
The protein input used here should already be chemically prepared upstream, i.e.:
- protonated
- true chain termini converted to AMBER-compatible terminal residues
- internal cut points ACE/NME-capped if needed

This module does NOT modify protein chemistry.
It only assembles the final prepared structure.

Default component paths
-----------------------
Given:
    data/proteins/<PDBID>/

default inputs are:
    components/<PDBID>_protein_internal_capped.pdb
    components/<PDBID>_water.pdb
    components/<PDBID>_ligand.pdb
    components/<PDBID>_metals.pdb

Directory logic
---------------
Case 1: no gaps
    prepared/<PDBID>.pdb

Case 2: gaps present, gapped structure variant
    prepared/gaps/<PDBID>.pdb

Case 3: gaps present, completed structure variant
    prepared/complete/<PDBID>.pdb

Notes
-----
- Only ATOM/HETATM records are merged.
- Atom serials are renumbered sequentially.
- Original chain IDs, residue numbers, residue names, coordinates, occupancies,
  B-factors, elements, and charges are preserved from input records.
- Existing TER/END records from component files are ignored and rebuilt cleanly.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Literal

StructureVariant = Literal["gaps", "complete"]


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


def get_default_prepared_protein_input_path(
    pdb_directory: str | Path,
    pdb_id: str,
) -> Path:
    """
    Return the default protein input path for prepared structure assembly.

    Current default:
        components/<PDBID>_protein_internal_capped.pdb
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
    had_gaps: bool,
    structure_variant: StructureVariant | None = None,
) -> Path:
    """
    Return the final prepared PDB output path.

    Rules
    -----
    - had_gaps=False:
        prepared/<PDBID>.pdb
    - had_gaps=True and structure_variant="gaps":
        prepared/gaps/<PDBID>.pdb
    - had_gaps=True and structure_variant="complete":
        prepared/complete/<PDBID>.pdb
    """
    pdb_directory = Path(pdb_directory)

    if not had_gaps:
        return pdb_directory / "prepared" / f"{pdb_id}.pdb"

    if structure_variant not in {"gaps", "complete"}:
        raise ValueError(
            "When had_gaps=True, structure_variant must be 'gaps' or 'complete'."
        )

    return pdb_directory / "prepared" / structure_variant / f"{pdb_id}.pdb"


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
    """
    Replace atom serial in columns 7-11 of a PDB ATOM/HETATM record.
    """
    return f"{line[:6]}{atom_serial:>5}{line[11:]}"


def _write_merged_pdb(
    output_pdb_path: str | Path,
    protein_atom_lines: list[str],
    water_atom_lines: list[str] | None = None,
    ligand_atom_lines: list[str] | None = None,
    metals_atom_lines: list[str] | None = None,
) -> int:
    output_pdb_path = Path(output_pdb_path)
    output_pdb_path.parent.mkdir(parents=True, exist_ok=True)

    ordered_sections = [
        protein_atom_lines,
        water_atom_lines or [],
        ligand_atom_lines or [],
        metals_atom_lines or [],
    ]

    atom_serial = 0
    n_written = 0

    with output_pdb_path.open("w", encoding="utf-8") as handle:
        for section_lines in ordered_sections:
            if not section_lines:
                continue

            for line in section_lines:
                atom_serial += 1
                handle.write(_renumber_atom_serial(line, atom_serial) + "\n")
                n_written += 1

            handle.write("TER\n")

        handle.write("END\n")

    return n_written


def build_prepared_structure(
    output_pdb_path: str | Path,
    protein_input_path: str | Path,
    water_input_path: str | Path | None = None,
    ligand_input_path: str | Path | None = None,
    metals_input_path: str | Path | None = None,
) -> PreparedStructureSummary:
    """
    Merge prepared protein with optional waters, ligands, and metals.

    Parameters
    ----------
    output_pdb_path
        Final combined output PDB path.
    protein_input_path
        Prepared protein PDB path. This should already be protonated and
        terminally prepared upstream.
    water_input_path
        Optional crystal water PDB path.
    ligand_input_path
        Optional ligand PDB path.
    metals_input_path
        Optional metals PDB path.

    Returns
    -------
    PreparedStructureSummary
        Assembly summary.
    """
    output_pdb_path = Path(output_pdb_path)
    protein_input_path = Path(protein_input_path)

    protein_atom_lines = _read_atom_lines_from_pdb(protein_input_path)

    if not protein_atom_lines:
        raise ValueError(
            f"No ATOM/HETATM records found in protein file: {protein_input_path}"
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

    n_atom_records_written = _write_merged_pdb(
        output_pdb_path=output_pdb_path,
        protein_atom_lines=protein_atom_lines,
        water_atom_lines=water_atom_lines,
        ligand_atom_lines=ligand_atom_lines,
        metals_atom_lines=metals_atom_lines,
    )

    return PreparedStructureSummary(
        pdb_id=output_pdb_path.stem,
        output_pdb_path=output_pdb_path,
        protein_input_path=protein_input_path,
        water_input_path=water_input if water_atom_lines else None,
        ligand_input_path=ligand_input if ligand_atom_lines else None,
        metals_input_path=metals_input if metals_atom_lines else None,
        had_gaps=False,
        structure_variant=None,
        water_included=bool(water_atom_lines),
        ligand_included=bool(ligand_atom_lines),
        metals_included=bool(metals_atom_lines),
        n_atom_records_written=n_atom_records_written,
    )


def build_prepared_structure_for_pdb_directory(
    pdb_directory: str | Path,
    pdb_id: str,
    had_gaps: bool,
    structure_variant: StructureVariant | None = None,
    protein_input_path: str | Path | None = None,
    water_input_path: str | Path | None = None,
    ligand_input_path: str | Path | None = None,
    metals_input_path: str | Path | None = None,
) -> PreparedStructureSummary:
    """
    Build the final prepared structure for one PDB directory.

    Parameters
    ----------
    pdb_directory
        Protein directory, for example data/proteins/1ABC
    pdb_id
        PDB ID, for example 1ABC
    had_gaps
        Whether the protein has gaps overall.
    structure_variant
        Required only when had_gaps=True:
        - "gaps" for the gapped structure
        - "complete" for the completed structure
    protein_input_path
        Optional explicit protein input path.
        If omitted, defaults to:
            components/<PDBID>_protein_internal_capped.pdb
    water_input_path
        Optional explicit water path.
        If omitted, defaults to:
            components/<PDBID>_water.pdb
    ligand_input_path
        Optional explicit ligand path.
        If omitted, defaults to:
            components/<PDBID>_ligand.pdb
    metals_input_path
        Optional explicit metals path.
        If omitted, defaults to:
            components/<PDBID>_metals.pdb

    Returns
    -------
    PreparedStructureSummary
        Assembly summary.
    """
    pdb_directory = Path(pdb_directory)

    output_pdb_path = get_prepared_structure_output_path(
        pdb_directory=pdb_directory,
        pdb_id=pdb_id,
        had_gaps=had_gaps,
        structure_variant=structure_variant,
    )

    resolved_protein_input_path = (
        Path(protein_input_path)
        if protein_input_path is not None
        else get_default_prepared_protein_input_path(
            pdb_directory=pdb_directory,
            pdb_id=pdb_id,
        )
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

    summary = build_prepared_structure(
        output_pdb_path=output_pdb_path,
        protein_input_path=resolved_protein_input_path,
        water_input_path=resolved_water_input_path,
        ligand_input_path=resolved_ligand_input_path,
        metals_input_path=resolved_metals_input_path,
    )

    summary.had_gaps = had_gaps
    summary.structure_variant = structure_variant

    return summary


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
        help="Use gap-aware output logic.",
    )
    argument_parser.add_argument(
        "--structure-variant",
        choices=["gaps", "complete"],
        default=None,
        help="Required when --had-gaps is set.",
    )
    argument_parser.add_argument(
        "--protein-input",
        type=Path,
        default=None,
        help="Optional explicit protein input path.",
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
