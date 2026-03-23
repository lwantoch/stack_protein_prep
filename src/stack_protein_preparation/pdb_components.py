"""
/home/grheco/repositorios/stack_protein_prep/src/stack_protein_preparation/pdb_components.py

Analyze and split PDB components into protein, water, metals, ligands, and
non-ligand crystallization artifacts.

Responsibilities
----------------
- detect metals by residue name or element symbol
- detect waters
- detect protein-like nonstandard residues
- detect ligand-like residues
- exclude common crystallization / buffer / salt artifacts from ligand counts
- split one input PDB into component-specific PDB files

Important
---------
- common trapped ions and crystallization additives such as SO4 are NOT treated
  as ligands
- protein-like modified residues such as MSE are NOT treated as ligands
- protein-like modified residues are written into the protein output PDB
- common artifacts are written into a separate artifact PDB for traceability

Notes
-----
- this module uses residue-based counting, not atom-based counting
- residue identity is defined by:
    (resname, chain_id, residue_number, insertion_code)
"""

from __future__ import annotations

from collections import Counter
from pathlib import Path
from typing import Any

METALS = {
    "LI",
    "NA",
    "K",
    "RB",
    "CS",
    "MG",
    "CA",
    "SR",
    "BA",
    "ZN",
    "FE",
    "CU",
    "MN",
    "CO",
    "NI",
    "CR",
    "AL",
    "AG",
    "AU",
    "HG",
    "CD",
}

WATER_NAMES = {
    "HOH",
    "WAT",
    "H2O",
    "TIP",
    "TIP3",
    "TIP3P",
    "SOL",
}

STANDARD_AMINO_ACIDS = {
    "ALA",
    "ARG",
    "ASN",
    "ASP",
    "CYS",
    "GLN",
    "GLU",
    "GLY",
    "HIS",
    "ILE",
    "LEU",
    "LYS",
    "MET",
    "PHE",
    "PRO",
    "SER",
    "THR",
    "TRP",
    "TYR",
    "VAL",
}

STANDARD_NUCLEIC_ACIDS = {
    "A",
    "C",
    "G",
    "U",
    "T",
    "DA",
    "DC",
    "DG",
    "DT",
    "DU",
}

STANDARD_POLYMER_RESIDUES = STANDARD_AMINO_ACIDS | STANDARD_NUCLEIC_ACIDS

PROTEIN_LIKE_NONSTANDARD = {
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
    "CYM",
    "CYX",
    "ASH",
    "GLH",
    "HID",
    "HIE",
    "HIP",
    "LYN",
}

# Common crystallization additives, salts, counterions, and buffer molecules
# that should not be treated as ligands.
CRYSTALLIZATION_ARTIFACTS = {
    "SO4",
    "PO4",
    "NO3",
    "CO3",
    "SCN",
    "CL",
    "BR",
    "IOD",
    "F",
    "NH4",
    "ACT",
    "ACE",
    "FMT",
    "CIT",
    "TAR",
    "MLI",
    "MES",
    "HEP",
    "TRS",
    "BME",
    "DTT",
    "EDO",
    "GOL",
    "PEG",
    "PG4",
    "PGE",
    "MPD",
    "IPA",
    "EOH",
    "DMS",
    "DMF",
    "ACY",
}

PROTEIN_COMPONENT_RESIDUES = STANDARD_POLYMER_RESIDUES | PROTEIN_LIKE_NONSTANDARD


def _read_pdb_lines(pdb_path: Path) -> list[str]:
    with pdb_path.open("r", encoding="utf-8") as handle:
        return handle.readlines()


def _is_atom_or_hetatm(line: str) -> bool:
    return line.startswith("ATOM") or line.startswith("HETATM")


def _resname(line: str) -> str:
    return line[17:20].strip().upper()


def _chain_id(line: str) -> str:
    return line[21:22].strip() or "?"


def _resseq(line: str) -> str:
    return line[22:26].strip()


def _icode(line: str) -> str:
    return line[26:27].strip()


def _element(line: str) -> str:
    return line[76:78].strip().upper()


def _residue_identifier(line: str) -> tuple[str, str, str, str]:
    return (
        _resname(line),
        _chain_id(line),
        _resseq(line),
        _icode(line),
    )


def _is_water_line(line: str) -> bool:
    return line.startswith("HETATM") and _resname(line) in WATER_NAMES


def _is_metal_line(line: str) -> bool:
    if not line.startswith("HETATM"):
        return False

    residue_name = _resname(line)
    element_symbol = _element(line)

    return residue_name in METALS or element_symbol in METALS


def _is_protein_like_residue_line(line: str) -> bool:
    """
    Return True for residues that belong to the polymer component, even if they
    are encoded as HETATM (for example MSE).
    """
    residue_name = _resname(line)
    return residue_name in PROTEIN_COMPONENT_RESIDUES


def _is_artifact_line(line: str) -> bool:
    """
    Return True for common non-ligand crystallization / buffer / salt artifacts.
    """
    if not line.startswith("HETATM"):
        return False

    residue_name = _resname(line)
    if residue_name in CRYSTALLIZATION_ARTIFACTS:
        return True

    return False


def _count_distinct_residues(
    pdb_lines: list[str],
    predicate,
) -> dict[str, int]:
    counts: Counter[str] = Counter()
    seen: set[tuple[str, str, str, str]] = set()

    for line in pdb_lines:
        if not predicate(line):
            continue

        residue_identifier = _residue_identifier(line)
        if residue_identifier in seen:
            continue

        counts[_resname(line)] += 1
        seen.add(residue_identifier)

    return dict(counts)


def _search_metals(pdb_lines: list[str]) -> dict[str, int]:
    """
    Return:
        { metal_symbol_or_resname : number_of_distinct_metal_residues }

    Only HETATM records are considered.
    Counting is residue-based, not atom-based.
    """
    metal_counts: Counter[str] = Counter()
    seen: set[tuple[str, str, str, str]] = set()

    for line in pdb_lines:
        if not line.startswith("HETATM"):
            continue

        residue_name = _resname(line)
        element_symbol = _element(line)

        detected_metal_symbol: str | None = None

        if element_symbol in METALS:
            detected_metal_symbol = element_symbol
        elif residue_name in METALS:
            detected_metal_symbol = residue_name

        if detected_metal_symbol is None:
            continue

        residue_identifier = (
            detected_metal_symbol,
            _chain_id(line),
            _resseq(line),
            _icode(line),
        )

        if residue_identifier in seen:
            continue

        metal_counts[detected_metal_symbol] += 1
        seen.add(residue_identifier)

    return dict(metal_counts)


def _search_waters(pdb_lines: list[str]) -> dict[str, int]:
    return _count_distinct_residues(pdb_lines, _is_water_line)


def _search_artifacts(pdb_lines: list[str]) -> dict[str, int]:
    return _count_distinct_residues(pdb_lines, _is_artifact_line)


def _search_nonstandard_residues(pdb_lines: list[str]) -> dict[str, int]:
    """
    Find non-standard polymer-like residues in ATOM/HETATM records.

    This catches residues such as:
    - MSE
    - SEP
    - TPO
    - HYP
    etc.

    It does NOT include:
    - waters
    - metals
    - crystallization artifacts
    - standard amino acids / nucleic acids
    """
    counts: Counter[str] = Counter()
    seen: set[tuple[str, str, str, str]] = set()

    for line in pdb_lines:
        if not _is_atom_or_hetatm(line):
            continue

        residue_name = _resname(line)

        if residue_name not in PROTEIN_LIKE_NONSTANDARD:
            continue

        residue_identifier = _residue_identifier(line)
        if residue_identifier in seen:
            continue

        counts[residue_name] += 1
        seen.add(residue_identifier)

    return dict(counts)


def _collect_ligand_like_residues(pdb_lines: list[str]) -> dict[str, int]:
    """
    Collect true ligand-like HETATM residues.

    Excludes:
    - waters
    - metals
    - common crystallization artifacts (e.g. SO4)
    - protein-like nonstandard residues (e.g. MSE)
    """
    counts: Counter[str] = Counter()
    seen: set[tuple[str, str, str, str]] = set()

    for line in pdb_lines:
        if not line.startswith("HETATM"):
            continue

        if _is_water_line(line):
            continue
        if _is_metal_line(line):
            continue
        if _is_artifact_line(line):
            continue
        if _is_protein_like_residue_line(line):
            continue

        residue_identifier = _residue_identifier(line)
        if residue_identifier in seen:
            continue

        counts[_resname(line)] += 1
        seen.add(residue_identifier)

    return dict(counts)


def analyze_pdb_components(pdb_path: Path) -> dict[str, Any]:
    """
    Analyze a PDB file and summarize:
    - metals
    - waters
    - crystallization / buffer artifacts
    - ligand candidates
    - non-standard polymer residues
    """
    pdb_lines = _read_pdb_lines(pdb_path)

    metals = _search_metals(pdb_lines)
    waters = _search_waters(pdb_lines)
    artifacts = _search_artifacts(pdb_lines)
    ligands = _collect_ligand_like_residues(pdb_lines)
    nonstandard_residues = _search_nonstandard_residues(pdb_lines)

    return {
        "input_pdb": str(pdb_path),
        "metals": metals,
        "waters": waters,
        "artifacts": artifacts,
        "ligands": ligands,
        "nonstandard_residues": nonstandard_residues,
        "has_metals": len(metals) > 0,
        "has_waters": len(waters) > 0,
        "has_artifacts": len(artifacts) > 0,
        "has_ligands": len(ligands) > 0,
        "has_nonstandard_residues": len(nonstandard_residues) > 0,
    }


def split_pdb_components(
    pdb_path: Path,
    output_dir: Path,
    protein_stem: str | None = None,
) -> dict[str, Any]:
    """
    Split one PDB into:
    - <PDBID>_protein.pdb
    - <PDBID>_water.pdb
    - <PDBID>_ligand.pdb
    - <PDBID>_metal.pdb
    - <PDBID>_artifact.pdb

    Classification
    --------------
    - protein:
        ATOM lines + protein-like nonstandard HETATM residues
    - water:
        HETATM water residues
    - metal:
        HETATM metal residues
    - artifact:
        HETATM crystallization / buffer / salt artifacts such as SO4
    - ligand:
        remaining HETATM residues after exclusions
    """
    pdb_lines = _read_pdb_lines(pdb_path)

    if protein_stem is None:
        protein_stem = pdb_path.stem

    output_dir.mkdir(parents=True, exist_ok=True)

    protein_path = output_dir / f"{protein_stem}_protein.pdb"
    water_path = output_dir / f"{protein_stem}_water.pdb"
    ligand_path = output_dir / f"{protein_stem}_ligand.pdb"
    metal_path = output_dir / f"{protein_stem}_metal.pdb"
    artifact_path = output_dir / f"{protein_stem}_artifact.pdb"

    protein_lines: list[str] = []
    water_lines: list[str] = []
    ligand_lines: list[str] = []
    metal_lines: list[str] = []
    artifact_lines: list[str] = []

    for line in pdb_lines:
        if line.startswith("ATOM"):
            protein_lines.append(line)
            continue

        if not line.startswith("HETATM"):
            continue

        if _is_water_line(line):
            water_lines.append(line)
        elif _is_metal_line(line):
            metal_lines.append(line)
        elif _is_protein_like_residue_line(line):
            protein_lines.append(line)
        elif _is_artifact_line(line):
            artifact_lines.append(line)
        else:
            ligand_lines.append(line)

    protein_path.write_text("".join(protein_lines), encoding="utf-8")
    water_path.write_text("".join(water_lines), encoding="utf-8")
    ligand_path.write_text("".join(ligand_lines), encoding="utf-8")
    metal_path.write_text("".join(metal_lines), encoding="utf-8")
    artifact_path.write_text("".join(artifact_lines), encoding="utf-8")

    summary = analyze_pdb_components(pdb_path)
    summary.update(
        {
            "output_dir": str(output_dir),
            "protein_pdb": str(protein_path),
            "water_pdb": str(water_path),
            "ligand_pdb": str(ligand_path),
            "metal_pdb": str(metal_path),
            "artifact_pdb": str(artifact_path),
            "n_protein_lines": len(protein_lines),
            "n_water_lines": len(water_lines),
            "n_ligand_lines": len(ligand_lines),
            "n_metal_lines": len(metal_lines),
            "n_artifact_lines": len(artifact_lines),
        }
    )
    return summary
