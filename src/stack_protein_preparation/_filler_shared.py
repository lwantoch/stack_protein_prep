"""Primitive types, constants, and low-level PDB utilities for filler submodules."""
from __future__ import annotations

import math
import os
import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

from Bio.Data.PDBData import protein_letters_3to1_extended
from Bio.PDB import PDBIO, PDBParser
from Bio.PDB.PDBIO import Select
from Bio.PDB.Polypeptide import is_aa

MODULE_NAME = "filler"

DEFAULT_FINAL_MODEL_FILENAME = "final_filled_model.pdb"
DEFAULT_ALIGNMENT_FILENAME = "chain_alignment.ali"
DEFAULT_MODELLER_SCRIPT_FILENAME = "run_modeller.py"
DEFAULT_SCORE_FILENAME = "model_scores.tsv"
DEFAULT_STDOUT_LOG_FILENAME = "modeller_stdout.log"
DEFAULT_STDERR_LOG_FILENAME = "modeller_stderr.log"
DEFAULT_MANIFEST_FILENAME = "filler_manifest.json"
DEFAULT_ALPHAFOLD_DIRNAME = "alphafold"
DEFAULT_ALPHAFOLD_RAW_FILENAME = "alphafold_downloaded_model.pdb"
DEFAULT_ALPHAFOLD_CROPPED_FILENAME = "alphafold_cropped_model.pdb"
DEFAULT_ALPHAFOLD_ALIGNED_FILENAME = "alphafold_aligned_model.pdb"

FORCE_FIELD_RESIDUE_ALIASES = {
    "ASH": "ASP",
    "GLH": "GLU",
    "HID": "HIS",
    "HIE": "HIS",
    "HIP": "HIS",
    "CYM": "CYS",
    "CYX": "CYS",
    "LYN": "LYS",
}


@dataclass(frozen=True)
class GapRegion:
    alignment_start: int
    alignment_end: int
    gap_length: int
    is_terminal: bool
    classification: str


@dataclass(frozen=True)
class FillDecision:
    should_run_modeller: bool
    overall_classification: str
    gap_regions: tuple[GapRegion, ...]
    skip_reason: str | None
    alphafold_candidate: bool


@dataclass(frozen=True)
class ContactAtom:
    """One protein atom to be held near its crystal position during MODELLER runs."""
    atom_name: str   # PDB atom name, e.g. "OD2"
    resseq: int      # PDB residue sequence number
    x: float         # crystal structure coordinate
    y: float
    z: float
    stdev: float = field(default=0.2)  # Gaussian restraint stdev in Å


@dataclass(frozen=True)
class FillerRunResult:
    chain_id: str
    output_dir: Path
    alignment_file: Path
    template_pdb: Path
    script_file: Path
    modeller_model_path: Path | None
    alphafold_model_path: Path | None
    final_model_path: Path | None
    raw_model_paths: tuple[Path, ...]
    stdout_log: Path
    stderr_log: Path
    skipped: bool
    skip_reason: str | None
    fill_decision: FillDecision


@dataclass(frozen=True)
class _ParsedFirstModel:
    structure: Any
    model: Any


class _FirstModelChainProteinSelect(Select):
    def __init__(self, *, model_id: Any, chain_id: str) -> None:
        self._model_id = model_id
        self._chain_id = _normalize_chain_id(chain_id)

    def accept_model(self, model: Any) -> bool:
        return model.id == self._model_id

    def accept_chain(self, chain: Any) -> bool:
        return _normalize_chain_id(str(chain.id)) == self._chain_id

    def accept_residue(self, residue: Any) -> bool:
        return _is_protein_atom_residue(residue)


class _FirstModelProteinRangeSelect(Select):
    def __init__(self, *, model_id: Any, start_residue: int, end_residue: int) -> None:
        self._model_id = model_id
        self._start_residue = start_residue
        self._end_residue = end_residue

    def accept_model(self, model: Any) -> bool:
        return model.id == self._model_id

    def accept_residue(self, residue: Any) -> bool:
        if not _is_protein_atom_residue(residue):
            return False
        return self._start_residue <= int(residue.id[1]) <= self._end_residue


class _FirstModelProteinOnlySelect(Select):
    def __init__(self, *, model_id: Any) -> None:
        self._model_id = model_id

    def accept_model(self, model: Any) -> bool:
        return model.id == self._model_id

    def accept_residue(self, residue: Any) -> bool:
        return _is_protein_atom_residue(residue)


def _debug(message: str) -> None:
    if os.environ.get("FRUTON_DEBUG"):
        print(f"[filler] {message}")


def _normalize_chain_id(chain_id: str) -> str:
    return chain_id.strip() or "_"


def _normalize_residue_name_for_sequence(residue_name: str) -> str:
    normalized_name = residue_name.strip().upper()
    return FORCE_FIELD_RESIDUE_ALIASES.get(normalized_name, normalized_name)


def _residue_name_to_one_letter(residue_name: str) -> str:
    normalized_name = _normalize_residue_name_for_sequence(residue_name)
    return protein_letters_3to1_extended.get(normalized_name, "X")


def _is_supported_protein_residue_name(residue_name: str) -> bool:
    normalized_name = _normalize_residue_name_for_sequence(residue_name)
    return normalized_name in protein_letters_3to1_extended


def _is_protein_atom_residue(residue: Any) -> bool:
    if residue.id[0] != " ":
        return False
    residue_name = residue.get_resname().strip().upper()
    return is_aa(residue, standard=False) or _is_supported_protein_residue_name(residue_name)


def _residue_has_peptide_backbone(residue: Any) -> bool:
    atom_name_set = {str(atom.id).strip().upper() for atom in residue}
    return {"N", "CA", "C"}.issubset(atom_name_set)


def _residue_to_one_letter(residue: Any) -> str:
    return _residue_name_to_one_letter(residue.get_resname())


def _parse_first_model(input_pdb_path: Path) -> _ParsedFirstModel:
    if not input_pdb_path.exists():
        raise FileNotFoundError(f"PDB file not found: {input_pdb_path}")
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure(input_pdb_path.stem, str(input_pdb_path))
    models = list(structure.get_models())
    if not models:
        raise ValueError(f"No coordinate models were found in {input_pdb_path}")
    return _ParsedFirstModel(structure=structure, model=models[0])


def _count_atoms_in_first_model(input_pdb_path: Path) -> int:
    if not input_pdb_path.exists():
        return 0
    try:
        parsed = _parse_first_model(input_pdb_path)
    except Exception:
        return 0
    return sum(len(list(r.get_atoms())) for c in parsed.model for r in c)


def _validate_written_pdb_has_atoms(output_pdb_path: Path, description: str) -> None:
    if not output_pdb_path.exists() or output_pdb_path.stat().st_size == 0:
        raise ValueError(f"{description} was not written: {output_pdb_path}")
    if _count_atoms_in_first_model(output_pdb_path) == 0:
        output_pdb_path.unlink(missing_ok=True)
        raise ValueError(f"{description} contains no atoms: {output_pdb_path}")


def _classify_gap_length(gap_length: int) -> str:
    if 1 <= gap_length <= 5:
        return "green"
    if 6 <= gap_length < 8:
        return "yellow"
    return "alphafold_candidate"


def _safe_str_path(path: Path | None) -> str | None:
    return None if path is None else str(path)


def _strip_gaps(sequence: str) -> str:
    return sequence.replace("-", "")


def _parse_residue_range(residue_range: str) -> tuple[int, int]:
    stripped = str(residue_range).strip()
    m = re.fullmatch(r"[A-Za-z](-?\d+)-[A-Za-z](-?\d+)", stripped)
    if m:
        return int(m.group(1)), int(m.group(2))
    m = re.fullmatch(r"\s*(-?\d+)\s*-\s*(-?\d+)\s*", stripped)
    if m is None:
        raise ValueError(f"Invalid residue range: {residue_range!r}")
    start_residue = int(m.group(1))
    end_residue = int(m.group(2))
    if start_residue > end_residue:
        raise ValueError(f"Invalid residue range: start > end in {residue_range!r}")
    return start_residue, end_residue


def _infer_full_residue_range_from_protein_pdb(protein_pdb_path: Path) -> str:
    parsed = _parse_first_model(protein_pdb_path)
    residue_numbers: list[int] = []
    for chain in parsed.model:
        for residue in chain:
            if not _is_protein_atom_residue(residue):
                continue
            if "CA" not in residue:
                continue
            residue_numbers.append(int(residue.id[1]))
    if not residue_numbers:
        raise ValueError(
            f"No protein CA residues found for range inference in {protein_pdb_path}"
        )
    return f"{min(residue_numbers)}-{max(residue_numbers)}"


def _get_uniprot_range_from_mapping(mapping_path: Path) -> str | None:
    """Return 'min_uniprot-max_uniprot' for positions where the PDB has a residue."""
    if not mapping_path.is_file():
        return None
    uniprot_positions: list[int] = []
    with mapping_path.open(encoding="utf-8") as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 6 or parts[5].strip() == "deletion_in_pdb":
                continue
            try:
                uniprot_positions.append(int(parts[2]))
            except (ValueError, IndexError):
                continue
    if not uniprot_positions:
        return None
    return f"{min(uniprot_positions)}-{max(uniprot_positions)}"


def _build_pdb_to_uniprot_residue_map(
    template_pdb_path: Path,
    mapping_path: Path,
) -> dict[int, int]:
    """Map actual PDB ATOM residue numbers → UniProt residue positions.

    The mapping TSV uses a sequential PDB index (1, 2, 3 …) that does not
    equal the residue sequence numbers stored in the ATOM records.  We
    reconstruct the correspondence by reading both sources in order.
    """
    if not mapping_path.is_file():
        return {}

    parsed = _parse_first_model(template_pdb_path)
    pdb_resnums: list[int] = []
    seen: set[int] = set()
    for chain in parsed.model:
        for residue in chain:
            if not _is_protein_atom_residue(residue):
                continue
            resnum = int(residue.id[1])
            if resnum not in seen:
                seen.add(resnum)
                pdb_resnums.append(resnum)

    seq_to_uniprot: list[tuple[int, int]] = []
    with mapping_path.open(encoding="utf-8") as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 6 or parts[5].strip() == "deletion_in_pdb":
                continue
            try:
                pdb_seq_idx = int(parts[1])
                uniprot_pos = int(parts[2])
                seq_to_uniprot.append((pdb_seq_idx, uniprot_pos))
            except (ValueError, IndexError):
                continue

    result: dict[int, int] = {}
    for pdb_seq_idx, uniprot_pos in seq_to_uniprot:
        arr_idx = pdb_seq_idx - 1
        if 0 <= arr_idx < len(pdb_resnums):
            result[pdb_resnums[arr_idx]] = uniprot_pos
    return result


def _read_hetatm_xyz(pdb_path: Path) -> list[tuple[float, float, float]]:
    """Return (x, y, z) for every HETATM line in *pdb_path*."""
    result: list[tuple[float, float, float]] = []
    if not pdb_path.is_file():
        return result
    for line in pdb_path.read_text(encoding="utf-8").splitlines():
        if not line.startswith("HETATM"):
            continue
        try:
            result.append((float(line[30:38]), float(line[38:46]), float(line[46:54])))
        except (ValueError, IndexError):
            continue
    return result


def find_metal_and_ligand_contact_atoms(
    protein_pdb_path: Path,
    metals_pdb_path: Path,
    ligand_pdb_path: Path,
    cofactor_pdb_path: Path,
    metal_cutoff: float = 3.0,
    ligand_cutoff: float = 4.0,
) -> list[ContactAtom]:
    """Return protein atoms within cutoff of any metal ion or ligand/cofactor atom.

    Uses *metal_cutoff* for atoms read from *metals_pdb_path* (tighter, for
    direct coordination bonds) and *ligand_cutoff* for atoms from the ligand
    and cofactor PDB files (looser, for non-covalent contacts).  Only ATOM-
    record residues in *protein_pdb_path* are considered.  All chains present
    in that file are searched (pass a chain-specific template to restrict).
    """
    # Each entry: (x, y, z, distance_cutoff, restraint_stdev)
    hetero_points: list[tuple[float, float, float, float, float]] = []

    for path, cutoff, stdev in [
        (metals_pdb_path, metal_cutoff, 0.1),
        (ligand_pdb_path, ligand_cutoff, 0.3),
        (cofactor_pdb_path, ligand_cutoff, 0.3),
    ]:
        for x, y, z in _read_hetatm_xyz(path):
            hetero_points.append((x, y, z, cutoff, stdev))

    if not hetero_points:
        return []

    parsed = _parse_first_model(protein_pdb_path)
    seen: set[tuple[str, int]] = set()
    contacts: list[ContactAtom] = []

    for chain in parsed.model:
        for residue in chain:
            if not _is_protein_atom_residue(residue):
                continue
            resseq = int(residue.id[1])
            for atom in residue:
                atom_name = str(atom.id).strip()
                ax, ay, az = float(atom.coord[0]), float(atom.coord[1]), float(atom.coord[2])
                for hx, hy, hz, cutoff, stdev in hetero_points:
                    if math.sqrt((ax - hx) ** 2 + (ay - hy) ** 2 + (az - hz) ** 2) <= cutoff:
                        key = (atom_name, resseq)
                        if key not in seen:
                            seen.add(key)
                            contacts.append(ContactAtom(
                                atom_name=atom_name,
                                resseq=resseq,
                                x=round(ax, 3),
                                y=round(ay, 3),
                                z=round(az, 3),
                                stdev=stdev,
                            ))
                        break

    return contacts


def _resolve_residue_range_for_filler(
    *,
    residue_range: str,
    template_pdb_path: Path,
) -> str:
    cleaned_range = str(residue_range).strip()
    if cleaned_range:
        _parse_residue_range(cleaned_range)
        return cleaned_range
    inferred_range = _infer_full_residue_range_from_protein_pdb(template_pdb_path)
    _debug(
        "No residue range provided; inferred full observed template range: "
        f"{inferred_range}"
    )
    return inferred_range


def _derive_per_chain_effective_range(
    *,
    residue_range: str,
    chain_id: str,
    template_pdb_path: Path,
) -> str:
    import re as _re
    stripped = residue_range.strip()
    if not stripped:
        return _infer_full_residue_range_from_protein_pdb(template_pdb_path)
    m = _re.fullmatch(r"([A-Za-z])(-?\d+)-([A-Za-z])(-?\d+)", stripped)
    if m is None:
        return stripped
    start_chain = m.group(1).upper()
    start_resnum = int(m.group(2))
    end_chain = m.group(3).upper()
    end_resnum = int(m.group(4))
    norm_chain = chain_id.strip().upper()
    template_range = _infer_full_residue_range_from_protein_pdb(template_pdb_path)
    template_start, template_end = _parse_residue_range(template_range)
    if start_chain == end_chain:
        if norm_chain == start_chain:
            return f"{start_resnum}-{end_resnum}"
        return f"{template_start}-{template_end}"
    if norm_chain == start_chain:
        return f"{start_resnum}-{template_end}"
    if norm_chain == end_chain:
        return f"{template_start}-{end_resnum}"
    return f"{template_start}-{template_end}"


def _trim_final_model_in_place(model_path: Path, effective_residue_range: str) -> None:
    if not effective_residue_range.strip():
        return
    start_resnum, end_resnum = _parse_residue_range(effective_residue_range)
    kept_lines: list[str] = []
    with model_path.open("r", encoding="utf-8") as fh:
        for line in fh:
            if line.startswith("ATOM"):
                raw_resseq = line[22:26].strip()
                try:
                    resnum = int(raw_resseq)
                except ValueError:
                    kept_lines.append(line)
                    continue
                if start_resnum <= resnum <= end_resnum:
                    kept_lines.append(line)
            else:
                kept_lines.append(line)
    model_path.write_text("".join(kept_lines), encoding="utf-8")


def _extract_backbone_residue_numbers(pdb_path: Path) -> list[int]:
    """Return ordered residue numbers for protein residues that have a peptide backbone.

    The order matches ``extract_sequence_from_template_pdb``, so position i in
    the returned list corresponds to the i-th non-gap character in the
    MODELLER template alignment sequence.
    """
    parsed = _parse_first_model(pdb_path)
    result: list[int] = []
    seen: set[int] = set()
    for chain in parsed.model:
        for residue in chain:
            if not _is_protein_atom_residue(residue):
                continue
            if not _residue_has_peptide_backbone(residue):
                continue
            resnum = int(residue.id[1])
            if resnum not in seen:
                seen.add(resnum)
                result.append(resnum)
    return result


def cleanup_model_pdb(input_model_path: Path, output_model_path: Path) -> Path:
    """Write a protein-only cleaned model keeping only ATOM records."""
    if not input_model_path.exists():
        raise FileNotFoundError(f"Model PDB not found: {input_model_path}")
    output_model_path.parent.mkdir(parents=True, exist_ok=True)
    atom_lines: list[str] = []
    with input_model_path.open("r", encoding="utf-8") as fh:
        for raw_line in fh:
            line = raw_line.rstrip("\n")
            if line.startswith("ATOM  "):
                atom_lines.append(line)
    if not atom_lines:
        raise ValueError(f"No ATOM records found in model PDB: {input_model_path}")
    with output_model_path.open("w", encoding="utf-8") as fh:
        for line in atom_lines:
            fh.write(line + "\n")
        fh.write("TER\n")
        fh.write("END\n")
    _debug(f"Cleaned final model written to: {output_model_path}")
    return output_model_path
