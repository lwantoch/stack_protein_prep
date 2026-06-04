"""AlphaFold-specific fallback mechanics for the filler pipeline."""
from __future__ import annotations

import json
from pathlib import Path
from typing import Any
from urllib.error import HTTPError, URLError
from urllib.request import urlopen, urlretrieve

from Bio.PDB import PDBIO, PDBParser, Superimposer

from ._filler_shared import (
    DEFAULT_ALPHAFOLD_ALIGNED_FILENAME,
    DEFAULT_ALPHAFOLD_CROPPED_FILENAME,
    DEFAULT_ALPHAFOLD_DIRNAME,
    DEFAULT_ALPHAFOLD_RAW_FILENAME,
    _FirstModelProteinRangeSelect,
    _debug,
    _is_protein_atom_residue,
    _parse_first_model,
    _parse_residue_range,
    _resolve_residue_range_for_filler,
    _validate_written_pdb_has_atoms,
    cleanup_model_pdb,
)


def download_alphafold_structure(uniprot_id: str, output_dir: Path) -> Path | None:
    output_dir.mkdir(parents=True, exist_ok=True)
    api_url = f"https://alphafold.ebi.ac.uk/api/prediction/{uniprot_id}"
    try:
        with urlopen(api_url) as response:
            payload = json.load(response)
    except HTTPError as exc:
        if exc.code == 404:
            return None
        raise
    except URLError:
        raise

    if not payload:
        return None
    record = payload[0]
    pdb_url = record.get("pdbUrl")
    if not pdb_url:
        return None

    target_path = output_dir / DEFAULT_ALPHAFOLD_RAW_FILENAME
    urlretrieve(pdb_url, target_path)
    if not target_path.is_file() or target_path.stat().st_size == 0:
        return None
    return target_path


def _collect_protein_ca_atoms_by_resseq(structure: Any) -> dict[int, object]:
    atom_by_resseq: dict[int, object] = {}
    for model in structure:
        for chain in model:
            for residue in chain:
                if not _is_protein_atom_residue(residue):
                    continue
                if "CA" not in residue:
                    continue
                residue_number = int(residue.id[1])
                atom_by_resseq.setdefault(residue_number, residue["CA"])
    return atom_by_resseq


def _collect_protein_ca_atoms_in_order(structure: Any) -> list[object]:
    atom_list: list[object] = []
    for model in structure:
        for chain in model:
            for residue in chain:
                if not _is_protein_atom_residue(residue):
                    continue
                if "CA" in residue:
                    atom_list.append(residue["CA"])
    return atom_list


def align_protonated_alphafold_model_to_start_pdb(
    reference_pdb_path: Path,
    mobile_pdb_path: Path,
    output_pdb_path: Path,
) -> dict[str, str | float | bool | int]:
    if not reference_pdb_path.exists():
        raise FileNotFoundError(reference_pdb_path)
    if not mobile_pdb_path.exists():
        raise FileNotFoundError(mobile_pdb_path)

    parser = PDBParser(QUIET=True)
    reference_structure = parser.get_structure("ref", str(reference_pdb_path))
    mobile_structure = parser.get_structure("mob", str(mobile_pdb_path))

    reference_by_resseq = _collect_protein_ca_atoms_by_resseq(reference_structure)
    mobile_by_resseq = _collect_protein_ca_atoms_by_resseq(mobile_structure)
    common_resseqs = sorted(set(reference_by_resseq) & set(mobile_by_resseq))
    fixed_atom_list = [reference_by_resseq[resseq] for resseq in common_resseqs]
    moving_atom_list = [mobile_by_resseq[resseq] for resseq in common_resseqs]
    pairing_mode = "residue_number"

    if len(fixed_atom_list) < 3:
        reference_ordered = _collect_protein_ca_atoms_in_order(reference_structure)
        mobile_ordered = _collect_protein_ca_atoms_in_order(mobile_structure)
        if not reference_ordered:
            raise ValueError(f"No protein CA atoms found in reference: {reference_pdb_path}")
        if not mobile_ordered:
            raise ValueError(f"No protein CA atoms found in mobile: {mobile_pdb_path}")
        n_matched_atoms = min(len(reference_ordered), len(mobile_ordered))
        if n_matched_atoms < 3:
            raise ValueError(
                "Not enough CA atoms for reliable alignment: "
                f"reference={len(reference_ordered)}, mobile={len(mobile_ordered)}"
            )
        fixed_atom_list = reference_ordered[:n_matched_atoms]
        moving_atom_list = mobile_ordered[:n_matched_atoms]
        pairing_mode = "positional_fallback"

    superimposer = Superimposer()
    superimposer.set_atoms(fixed_atom_list, moving_atom_list)
    superimposer.apply(mobile_structure.get_atoms())

    output_pdb_path.parent.mkdir(parents=True, exist_ok=True)
    io = PDBIO()
    io.set_structure(mobile_structure)
    io.save(str(output_pdb_path), preserve_atom_numbering=True)

    return {
        "alignment_success": output_pdb_path.exists() and output_pdb_path.stat().st_size > 0,
        "alignment_rmsd": float(superimposer.rms),
        "alignment_output_path": str(output_pdb_path),
        "alignment_pairing_mode": pairing_mode,
        "alignment_n_pairs": len(fixed_atom_list),
    }


def _count_protein_residues_in_range(
    input_pdb_path: Path,
    start_residue: int,
    end_residue: int,
) -> int:
    parsed = _parse_first_model(input_pdb_path)
    count = 0
    for chain in parsed.model:
        for residue in chain:
            if not _is_protein_atom_residue(residue):
                continue
            if start_residue <= int(residue.id[1]) <= end_residue:
                count += 1
    return count


def crop_pdb_to_range(
    input_pdb_path: Path,
    output_pdb_path: Path,
    residue_range: str,
) -> Path:
    """Crop a PDB to a residue-number range using BioPython selectors."""
    if input_pdb_path is None:
        raise FileNotFoundError("AlphaFold download returned no PDB path")
    if not input_pdb_path.exists():
        raise FileNotFoundError(f"Input PDB not found: {input_pdb_path}")

    start_residue, end_residue = _parse_residue_range(residue_range)
    kept_residue_count = _count_protein_residues_in_range(
        input_pdb_path=input_pdb_path,
        start_residue=start_residue,
        end_residue=end_residue,
    )
    if kept_residue_count == 0:
        raise ValueError(
            f"No protein residues remained after cropping {input_pdb_path} "
            f"to range {residue_range!r}"
        )

    parsed = _parse_first_model(input_pdb_path)
    output_pdb_path.parent.mkdir(parents=True, exist_ok=True)
    io = PDBIO()
    io.set_structure(parsed.structure)
    io.save(
        str(output_pdb_path),
        select=_FirstModelProteinRangeSelect(
            model_id=parsed.model.id,
            start_residue=start_residue,
            end_residue=end_residue,
        ),
        preserve_atom_numbering=True,
    )
    _validate_written_pdb_has_atoms(
        output_pdb_path=output_pdb_path,
        description=f"Cropped AlphaFold PDB for range {residue_range!r}",
    )
    return output_pdb_path


def run_alphafold_fallback_for_chain(
    output_dir: Path,
    template_pdb_path: Path,
    uniprot_id: str,
    residue_range: str,
    final_model_name: str,
    model_version: int = 4,
) -> Path:
    _ = model_version
    effective_residue_range = _resolve_residue_range_for_filler(
        residue_range=residue_range,
        template_pdb_path=template_pdb_path,
    )
    alphafold_dir = output_dir / DEFAULT_ALPHAFOLD_DIRNAME
    downloaded_pdb = download_alphafold_structure(
        uniprot_id=uniprot_id, output_dir=alphafold_dir
    )
    if downloaded_pdb is None:
        raise FileNotFoundError(
            f"No AlphaFold PDB available for UniProt {uniprot_id!r}"
        )
    cropped_model_path = crop_pdb_to_range(
        input_pdb_path=downloaded_pdb,
        output_pdb_path=alphafold_dir / DEFAULT_ALPHAFOLD_CROPPED_FILENAME,
        residue_range=effective_residue_range,
    )
    aligned_model_path = alphafold_dir / DEFAULT_ALPHAFOLD_ALIGNED_FILENAME
    alignment_result = align_protonated_alphafold_model_to_start_pdb(
        reference_pdb_path=template_pdb_path,
        mobile_pdb_path=cropped_model_path,
        output_pdb_path=aligned_model_path,
    )
    _debug(
        "AlphaFold alignment result: "
        f"success={alignment_result['alignment_success']}, "
        f"rmsd={alignment_result['alignment_rmsd']}, "
        f"pairs={alignment_result['alignment_n_pairs']}, "
        f"pairing_mode={alignment_result['alignment_pairing_mode']}, "
        f"output={alignment_result['alignment_output_path']}"
    )
    final_model_path = output_dir / final_model_name
    return cleanup_model_pdb(
        input_model_path=aligned_model_path,
        output_model_path=final_model_path,
    )
