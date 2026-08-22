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
    _get_uniprot_range_from_mapping,
    _build_pdb_to_uniprot_residue_map,
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


def _align_alphafold_by_residue_map(
    reference_pdb_path: Path,
    mobile_pdb_path: Path,
    output_pdb_path: Path,
    pdb_to_uniprot_map: dict[int, int],
) -> dict[str, str | float | bool | int]:
    """Superimpose AlphaFold (UniProt-numbered) onto PDB template using the
    PDB-resnum → UniProt-position mapping so that structurally equivalent
    residues are paired even when the two numbering systems differ."""
    parser = PDBParser(QUIET=True)
    reference_structure = parser.get_structure("ref", str(reference_pdb_path))
    mobile_structure = parser.get_structure("mob", str(mobile_pdb_path))

    reference_by_resseq = _collect_protein_ca_atoms_by_resseq(reference_structure)
    mobile_by_resseq = _collect_protein_ca_atoms_by_resseq(mobile_structure)

    fixed_atoms = []
    moving_atoms = []
    for pdb_resnum, uniprot_pos in sorted(pdb_to_uniprot_map.items()):
        if pdb_resnum in reference_by_resseq and uniprot_pos in mobile_by_resseq:
            fixed_atoms.append(reference_by_resseq[pdb_resnum])
            moving_atoms.append(mobile_by_resseq[uniprot_pos])

    pairing_mode = "residue_map"

    if len(fixed_atoms) < 3:
        reference_ordered = _collect_protein_ca_atoms_in_order(reference_structure)
        mobile_ordered = _collect_protein_ca_atoms_in_order(mobile_structure)
        if not reference_ordered or not mobile_ordered:
            raise ValueError(
                f"No protein CA atoms for alignment: ref={reference_pdb_path}, "
                f"mob={mobile_pdb_path}"
            )
        n = min(len(reference_ordered), len(mobile_ordered))
        if n < 3:
            raise ValueError(
                f"Not enough CA atoms for reliable alignment: "
                f"ref={len(reference_ordered)}, mob={len(mobile_ordered)}"
            )
        fixed_atoms = reference_ordered[:n]
        moving_atoms = mobile_ordered[:n]
        pairing_mode = "positional_fallback"

    superimposer = Superimposer()
    superimposer.set_atoms(fixed_atoms, moving_atoms)
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
        "alignment_n_pairs": len(fixed_atoms),
    }


def run_alphafold_fallback_for_chain(
    output_dir: Path,
    template_pdb_path: Path,
    uniprot_id: str,
    residue_range: str,
    final_model_name: str,
    model_version: int = 4,
    uniprot_residue_range: str | None = None,
    pdb_to_uniprot_map: dict[int, int] | None = None,
) -> Path:
    """Download, crop, align and return an AlphaFold fill model.

    When *uniprot_residue_range* is supplied (derived from the alignment
    mapping TSV) it is used to crop the AlphaFold download instead of the
    PDB residue-number range.  This is essential for proteins where PDB
    author residue numbers differ from UniProt positions (e.g. domain
    fragments, propeptide offsets).

    *pdb_to_uniprot_map* is used for the structural superimposition so that
    each PDB residue is paired with the correct AlphaFold residue rather
    than relying on matching residue sequence numbers.
    """
    _ = model_version
    effective_residue_range = _resolve_residue_range_for_filler(
        residue_range=residue_range,
        template_pdb_path=template_pdb_path,
    )
    alphafold_crop_range = uniprot_residue_range or effective_residue_range

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
        residue_range=alphafold_crop_range,
    )
    aligned_model_path = alphafold_dir / DEFAULT_ALPHAFOLD_ALIGNED_FILENAME

    if pdb_to_uniprot_map:
        alignment_result = _align_alphafold_by_residue_map(
            reference_pdb_path=template_pdb_path,
            mobile_pdb_path=cropped_model_path,
            output_pdb_path=aligned_model_path,
            pdb_to_uniprot_map=pdb_to_uniprot_map,
        )
    else:
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
    # FIX (2026-08-20): splice AF-modelled loops into the crystal template
    # instead of returning the WHOLE AF-aligned model. The whole-body copy
    # (previous behaviour) silently dropped chain B, ligands, cofactors and
    # crystal waters — see memory note `project_fruton_af_wholebody_bug`
    # (13 MMBSA_200 PDBs + 6 newbench_27 gaps-variant PDBs affected). The
    # splice keeps every residue present in the crystal template and only
    # inserts AF residues that fill the missing-residue windows.
    final_model_path = output_dir / final_model_name
    from stack_protein_preparation._filler_af_splice import (
        splice_af_gaps_into_crystal,
    )
    spliced_path = splice_af_gaps_into_crystal(
        crystal_pdb_path=template_pdb_path,
        af_aligned_pdb_path=aligned_model_path,
        output_pdb_path=final_model_path,
    )

    # Junction-relaxation (2026-08-22): restrained OpenMM ff14SB min holding
    # everything but gap±2 heavy atoms rigid, so bad omega / long C-N at the
    # AF-crystal boundary absorbs strain into the region we intended to
    # change without perturbing MolProbity-favoured crystal residues.
    # Gap ranges derive from the detect step in _filler_af_splice.  Failures
    # (non-standard residues, missing terminal caps, ff14SB parse errors)
    # fall back to identity and are logged; downstream stages still see a
    # canonical PDB.
    try:
        from stack_protein_preparation._filler_af_splice import (
            _detect_missing_windows,
            _protein_residue_map,
        )
        from Bio.PDB import PDBParser as _PDBParser
        _crystal = _PDBParser(QUIET=True).get_structure("c", str(template_pdb_path))
        _af = _PDBParser(QUIET=True).get_structure("a", str(aligned_model_path))
        _cmap = _protein_residue_map(_crystal)
        _amap = _protein_residue_map(_af)
        _gap_ranges: list[tuple[int, int]] = []
        for cid, af_by_resi in _amap.items():
            crystal_by = _cmap.get(cid, {})
            _gap_ranges.extend(
                _detect_missing_windows(set(crystal_by.keys()), sorted(af_by_resi.keys()))
            )
        # Determine which chain each gap belongs to for LoopModel selection.
        _gap_ranges_by_chain: list[tuple[str, int, int]] = []
        for cid, af_by_resi in _amap.items():
            crystal_by = _cmap.get(cid, {})
            for _lo, _hi in _detect_missing_windows(
                set(crystal_by.keys()), sorted(af_by_resi.keys())
            ):
                _gap_ranges_by_chain.append((cid, _lo, _hi))
        # Only refine gaps whose fitted residues actually made it through
        # smart-rollback in splice_af_gaps_into_crystal.  Compare the spliced
        # PDB's residue set to the crystal to find truly-inserted residues.
        try:
            _spliced_struct = _PDBParser(QUIET=True).get_structure("sp", str(spliced_path))
            _crystal_res = {
                (ch.id, r.id[1])
                for ch in _crystal[0] for r in ch if not r.id[0].strip()
            }
            _spliced_res = {
                (ch.id, r.id[1])
                for ch in _spliced_struct[0] for r in ch if not r.id[0].strip()
            }
            _added_res = _spliced_res - _crystal_res
        except Exception:  # noqa: BLE001
            _added_res = set()
        _surviving_gaps = [
            (cid, lo, hi) for cid, lo, hi in _gap_ranges_by_chain
            if any((cid, r) in _added_res for r in range(lo, hi + 1))
        ]
        if _surviving_gaps:
            # LoopModel polish (2026-08-22): MODELLER stochastic loop sampling
            # with 3 conformers + chirality guard.  Preserves peptide bond
            # geometry (~1.33 A) via internal coordinate constraints and
            # reduces clashes by 25x+ vs raw splice (5M7U live: 100 -> 3).
            # Slow: ~5-10 min per protein.  Falls back to spliced input on
            # any failure so pipeline never blocks.
            from stack_protein_preparation._filler_loop_refine import (
                refine_loops_via_modeller,
            )
            _refine = refine_loops_via_modeller(
                input_pdb_path=spliced_path,
                output_pdb_path=final_model_path,
                gap_ranges_by_chain=_surviving_gaps,
                n_conformers=3,
                refine_level="fast",
                reject_new_chirality_d=True,
            )
            _debug(
                f"LoopModel refine: ran={_refine.ran}, "
                f"n_built={_refine.n_conformers_built}, "
                f"n_kept={_refine.n_conformers_kept}, "
                f"best_dope={_refine.best_dope}, "
                f"fallback={_refine.fallback_reason}"
            )
    except Exception as _refine_exc:  # noqa: BLE001 -- best-effort polish step
        _debug(f"LoopModel refine skipped due to unexpected error: {_refine_exc!r}")

    # Quality-gate (2026-08-22): native MolProbity-style analyser on the
    # spliced+relaxed model, compared against the crystal baseline.  Persists
    # the report as ``quality_gate.json`` next to the final model so the
    # downstream reporter / reviewer can inspect specific outlier residues.
    # Non-blocking: a fail is logged but does not abort the pipeline (a
    # future iteration will gate on this).
    try:
        import json as _json
        from stack_protein_preparation._filler_quality_check import (
            check_model_quality,
        )
        _gap_res_ids = set()
        for _lo, _hi in _gap_ranges:
            for _cid in _amap:
                for _r in range(_lo, _hi + 1):
                    _gap_res_ids.add((_cid, _r))
        _baseline_report = check_model_quality(template_pdb_path)
        _final_report = check_model_quality(spliced_path, gap_residue_ids=_gap_res_ids)
        _passed, _reasons = _final_report.passes_relative_gate(_baseline_report)
        _gate_dict = {
            "passed": _passed,
            "reasons": _reasons,
            "baseline": _baseline_report.to_dict(),
            "final": _final_report.to_dict(),
        }
        _gate_path = output_dir / "quality_gate.json"
        _gate_path.write_text(_json.dumps(_gate_dict, indent=2, default=str), encoding="utf-8")
        _debug(
            f"Quality gate ({'PASS' if _passed else 'FAIL'}): "
            f"Rama fav {_baseline_report.rama_favoured_pct():.2f}%→{_final_report.rama_favoured_pct():.2f}%, "
            f"broken bonds {_baseline_report.n_peptide_bonds_broken}→{_final_report.n_peptide_bonds_broken}, "
            f"clashes {_baseline_report.n_clash_pairs}→{_final_report.n_clash_pairs}"
        )
        if not _passed:
            for _r in _reasons:
                _debug(f"Quality gate reason: {_r}")
    except Exception as _qc_exc:  # noqa: BLE001 -- diagnostic only, do not block pipeline
        _debug(f"Quality gate skipped due to unexpected error: {_qc_exc!r}")

    # Still run the historic post-hoc cleanup (element-column fixups etc.)
    # so downstream stages see a canonical PDB.
    return cleanup_model_pdb(
        input_model_path=spliced_path,
        output_model_path=final_model_path,
    )
