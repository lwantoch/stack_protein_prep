"""Fill missing residues in PDB structures using MODELLER or AlphaFold2 grafting.

This module is the public entry point. Implementation is split across:
  _filler_shared   — dataclasses, constants, PDB primitives
  _filler_analysis — gap analysis, alignment, chain selection
  _filler_modeller — MODELLER mechanics
  _filler_alphafold — AlphaFold download and alignment
"""
from __future__ import annotations

from dataclasses import asdict
from pathlib import Path
from typing import Any

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

from ._filler_shared import (
    DEFAULT_ALIGNMENT_FILENAME,
    DEFAULT_FINAL_MODEL_FILENAME,
    DEFAULT_MANIFEST_FILENAME,
    DEFAULT_MODELLER_SCRIPT_FILENAME,
    FORCE_FIELD_RESIDUE_ALIASES,
    MODULE_NAME,
    ContactAtom,
    FillDecision,
    FillerRunResult,
    GapRegion,
    _ParsedFirstModel,
    _classify_gap_length,
    _count_atoms_in_first_model,
    _derive_per_chain_effective_range,
    _infer_full_residue_range_from_protein_pdb,
    _is_protein_atom_residue,
    _is_supported_protein_residue_name,
    _normalize_residue_name_for_sequence,
    _parse_first_model,
    _parse_residue_range,
    _residue_has_peptide_backbone,
    _residue_name_to_one_letter,
    _residue_to_one_letter,
    _resolve_residue_range_for_filler,
    _safe_str_path,
    _strip_gaps,
    _trim_final_model_in_place,
    _validate_written_pdb_has_atoms,
    cleanup_model_pdb,
    find_metal_and_ligand_contact_atoms,
    _get_uniprot_range_from_mapping,
    _build_pdb_to_uniprot_residue_map,
)
from ._filler_analysis import (
    _alignment_metrics_for_chain,
    _chain_has_protein_atoms,
    _collect_protein_residue_numbers_by_chain,
    _extract_chain_id_from_alignment_header,
    _extract_sequence_block_from_ali,
    _find_gap_regions_in_sequence,
    _first_sequence_mismatch,
    _is_pdb_atom_header,
    _is_uniprot_header,
    _normalize_alignment_template_for_validation,
    _normalize_chain_id,
    _write_sequence_debug_files,
    analyze_fill_decision,
    analyze_fill_decision_from_template_alignment,
    build_modeller_template_alignment_sequence,
    choose_filler_chain_id_from_range_and_alignments,
    extract_sequence_from_template_pdb,
    find_alignment_fasta_for_filler,
    list_available_atom_alignment_chain_ids,
    read_two_sequence_fasta,
    resolve_filler_chain_id,
    split_template_and_target_alignment_records,
    validate_template_sequence_consistency,
    write_chain_specific_template_pdb,
)
from ._filler_modeller import (
    MODELLER_PYTHONPATH,
    PYTHON_BIN,
    _parse_model_scores,
    _write_skip_logs,
    find_raw_models,
    run_modeller_binary,
    select_best_model,
    select_best_model_from_scores,
    standardize_model_name,
    validate_modeller_inputs,
    write_modeller_alignment_from_existing_alignment,
    write_modeller_script,
)
from ._filler_alphafold import (
    _collect_protein_ca_atoms_by_resseq,
    _collect_protein_ca_atoms_in_order,
    _count_protein_residues_in_range,
    align_protonated_alphafold_model_to_start_pdb,
    crop_pdb_to_range,
    download_alphafold_structure,
    run_alphafold_fallback_for_chain,
)

# Re-export everything so existing callers remain unaffected
__all__ = [
    "GapRegion",
    "FillDecision",
    "FillerRunResult",
    "MODULE_NAME",
    "FORCE_FIELD_RESIDUE_ALIASES",
    "PYTHON_BIN",
    "MODELLER_PYTHONPATH",
    "DEFAULT_FINAL_MODEL_FILENAME",
    "DEFAULT_ALIGNMENT_FILENAME",
    "DEFAULT_MANIFEST_FILENAME",
    "cleanup_model_pdb",
    "ContactAtom",
    "find_metal_and_ligand_contact_atoms",
    "analyze_fill_decision",
    "analyze_fill_decision_from_template_alignment",
    "choose_filler_chain_id_from_range_and_alignments",
    "resolve_filler_chain_id",
    "find_alignment_fasta_for_filler",
    "list_available_atom_alignment_chain_ids",
    "read_two_sequence_fasta",
    "split_template_and_target_alignment_records",
    "write_chain_specific_template_pdb",
    "extract_sequence_from_template_pdb",
    "validate_template_sequence_consistency",
    "build_modeller_template_alignment_sequence",
    "write_modeller_alignment_from_existing_alignment",
    "write_modeller_script",
    "validate_modeller_inputs",
    "run_modeller_binary",
    "find_raw_models",
    "select_best_model",
    "select_best_model_from_scores",
    "standardize_model_name",
    "download_alphafold_structure",
    "align_protonated_alphafold_model_to_start_pdb",
    "crop_pdb_to_range",
    "run_alphafold_fallback_for_chain",
    "run_filler_for_chain",
    "_normalize_chain_id",
    "_extract_chain_id_from_alignment_header",
    "_find_gap_regions_in_sequence",
    "_write_skip_logs",
]


def _debug(message: str) -> None:
    import os
    if os.environ.get("FRUTON_DEBUG"):
        print(f"[filler] {message}")


def _infer_pdb_id_from_target_id(target_id: str, template_pdb_path: Path) -> str:
    target_id = str(target_id).strip()
    if "_" in target_id:
        return target_id.split("_", maxsplit=1)[0].strip().upper()
    stem = template_pdb_path.stem.strip().upper()
    if stem:
        return stem.split("_", maxsplit=1)[0]
    raise ValueError(
        f"Could not infer pdb_id from target_id={target_id!r} "
        f"and template_pdb_path={template_pdb_path}"
    )


def _infer_protein_data_dir(template_pdb_path: Path) -> Path:
    template_pdb_path = template_pdb_path.resolve()
    protein_dir = template_pdb_path.parent
    if protein_dir.name == "components":
        protein_dir = protein_dir.parent
    return protein_dir.parent


def _build_module_log_path(*, pdb_id: str, template_pdb_path: Path) -> Path:
    protein_data_dir = _infer_protein_data_dir(template_pdb_path)
    log_paths = build_module_log_paths(
        protein_data_dir=protein_data_dir,
        pdb_id=pdb_id,
        module_name=MODULE_NAME,
        variant_label=None,
    )
    return log_paths.module_log_path


def _log_debug(log_path: Path | None, message: str) -> None:
    _debug(message)
    if log_path is not None:
        append_log_text(log_path, f"[DEBUG] {message}")


def _print_start_box(
    *,
    pdb_id: str,
    chain_id: str | None,
    template_id: str,
    target_id: str,
    alignment_directory: Path,
    template_pdb_path: Path,
    output_dir: Path,
    residue_range: str,
    uniprot_id: str | None,
    log_path: Path,
) -> None:
    print_mixed_box(
        title=f"{MODULE_NAME} · START",
        line_list=[
            f"PDB ID         : {pdb_id}",
            f"Chain input    : {chain_id if chain_id is not None else '(auto)'}",
            f"Template ID    : {template_id}",
            f"Target ID      : {target_id}",
            f"Residue range  : {residue_range or '(infer full observed range)'}",
            f"UniProt ID     : {uniprot_id or '(none)'}",
            f"Alignment dir  : {alignment_directory}",
            f"Template PDB   : {template_pdb_path}",
            f"Output dir     : {output_dir}",
            f"Log file       : {log_path}",
        ],
    )


def _print_decision_table(
    *,
    pdb_id: str,
    resolved_chain_id: str,
    fill_decision: FillDecision,
) -> None:
    row_list = [
        ("PDB ID", pdb_id),
        ("Resolved chain", resolved_chain_id),
        ("should_run_modeller", str(fill_decision.should_run_modeller)),
        ("overall_classification", fill_decision.overall_classification),
        ("alphafold_candidate", str(fill_decision.alphafold_candidate)),
        ("n_gap_regions", str(len(fill_decision.gap_regions))),
        ("skip_reason", fill_decision.skip_reason or "(none)"),
    ]
    print_light_table_box(
        title=f"{MODULE_NAME} · decision",
        row_list=row_list,
        left_header="Field",
        right_header="Value",
    )


def _print_gap_regions_table(fill_decision: FillDecision) -> None:
    if not fill_decision.gap_regions:
        print_light_table_box(
            title=f"{MODULE_NAME} · gap regions",
            row_list=[("result", "no internal gap regions")],
            left_header="Field",
            right_header="Value",
        )
        return
    row_list = []
    for index, gap_region in enumerate(fill_decision.gap_regions, start=1):
        row_list.append(
            (
                f"gap_{index}",
                (
                    f"{gap_region.alignment_start}-{gap_region.alignment_end} | "
                    f"len={gap_region.gap_length} | "
                    f"terminal={gap_region.is_terminal} | "
                    f"class={gap_region.classification}"
                ),
            )
        )
    print_light_table_box(
        title=f"{MODULE_NAME} · gap regions",
        row_list=row_list,
        left_header="Region",
        right_header="Details",
    )


def _print_end_box(
    *,
    pdb_id: str,
    resolved_chain_id: str,
    status: str,
    final_model_path: Path | None,
    modeller_model_path: Path | None,
    alphafold_model_path: Path | None,
    stdout_log: Path,
    stderr_log: Path,
    skip_reason: str | None,
) -> None:
    print_double_header_box(
        title=f"{MODULE_NAME} · END",
        line_list=[
            f"PDB ID           : {pdb_id}",
            f"Resolved chain   : {resolved_chain_id}",
            f"Status           : {status}",
            f"Final model      : {final_model_path if final_model_path else '(none)'}",
            f"MODELLER model   : {modeller_model_path if modeller_model_path else '(none)'}",
            f"AlphaFold model  : {alphafold_model_path if alphafold_model_path else '(none)'}",
            f"stdout log       : {stdout_log}",
            f"stderr log       : {stderr_log}",
            f"skip_reason      : {skip_reason or '(none)'}",
        ],
    )


def _result_to_manifest_payload(result: FillerRunResult) -> dict[str, Any]:
    return {
        "chain_id": result.chain_id,
        "output_dir": str(result.output_dir),
        "alignment_file": str(result.alignment_file),
        "template_pdb": str(result.template_pdb),
        "script_file": str(result.script_file),
        "modeller_model_path": _safe_str_path(result.modeller_model_path),
        "alphafold_model_path": _safe_str_path(result.alphafold_model_path),
        "final_model_path": _safe_str_path(result.final_model_path),
        "raw_model_paths": [str(path) for path in result.raw_model_paths],
        "stdout_log": str(result.stdout_log),
        "stderr_log": str(result.stderr_log),
        "skipped": result.skipped,
        "skip_reason": result.skip_reason,
        "fill_decision": {
            "should_run_modeller": result.fill_decision.should_run_modeller,
            "overall_classification": result.fill_decision.overall_classification,
            "gap_regions": [asdict(region) for region in result.fill_decision.gap_regions],
            "skip_reason": result.fill_decision.skip_reason,
            "alphafold_candidate": result.fill_decision.alphafold_candidate,
        },
    }


def _write_manifest(output_dir: Path, result: FillerRunResult) -> Path:
    import json
    manifest_path = output_dir / DEFAULT_MANIFEST_FILENAME
    payload = _result_to_manifest_payload(result)
    manifest_path.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    return manifest_path


def _build_result(
    *,
    chain_id: str,
    output_dir: Path,
    alignment_file: Path,
    template_pdb: Path,
    script_file: Path,
    modeller_model_path: Path | None,
    alphafold_model_path: Path | None,
    final_model_path: Path | None,
    raw_model_paths: tuple[Path, ...],
    stdout_log: Path,
    stderr_log: Path,
    skipped: bool,
    skip_reason: str | None,
    fill_decision: FillDecision,
) -> FillerRunResult:
    result = FillerRunResult(
        chain_id=chain_id,
        output_dir=output_dir,
        alignment_file=alignment_file,
        template_pdb=template_pdb,
        script_file=script_file,
        modeller_model_path=modeller_model_path,
        alphafold_model_path=alphafold_model_path,
        final_model_path=final_model_path,
        raw_model_paths=raw_model_paths,
        stdout_log=stdout_log,
        stderr_log=stderr_log,
        skipped=skipped,
        skip_reason=skip_reason,
        fill_decision=fill_decision,
    )
    _write_manifest(output_dir=output_dir, result=result)
    return result


def run_filler_for_chain(
    alignment_directory: Path,
    template_pdb_path: Path,
    output_dir: Path,
    template_id: str,
    target_id: str,
    chain_id: str | None,
    final_model_name: str | None = None,
    starting_model: int = 1,
    ending_model: int = 20,
    skip_if_no_internal_gaps: bool = True,
    uniprot_id: str | None = None,
    residue_range: str = "",
) -> FillerRunResult:
    pdb_id = _infer_pdb_id_from_target_id(
        target_id=target_id, template_pdb_path=template_pdb_path
    )
    module_log_path = _build_module_log_path(
        pdb_id=pdb_id, template_pdb_path=template_pdb_path
    )

    write_log_header(
        log_path=module_log_path,
        module_name=MODULE_NAME,
        pdb_id=pdb_id,
        extra_lines=[
            f"template_id={template_id}",
            f"target_id={target_id}",
            f"input_chain_id={chain_id if chain_id is not None else '(auto)'}",
            f"alignment_directory={alignment_directory}",
            f"template_pdb_path={template_pdb_path}",
            f"output_dir={output_dir}",
            f"starting_model={starting_model}",
            f"ending_model={ending_model}",
            f"skip_if_no_internal_gaps={skip_if_no_internal_gaps}",
            f"uniprot_id={uniprot_id or '(none)'}",
            f"residue_range={residue_range or '(infer full observed range)'}",
        ],
    )

    _print_start_box(
        pdb_id=pdb_id,
        chain_id=chain_id,
        template_id=template_id,
        target_id=target_id,
        alignment_directory=alignment_directory,
        template_pdb_path=template_pdb_path,
        output_dir=output_dir,
        residue_range=residue_range,
        uniprot_id=uniprot_id,
        log_path=module_log_path,
    )

    _log_debug(module_log_path, "=== run_filler_for_chain START ===")

    if not alignment_directory.exists():
        error = FileNotFoundError(f"Alignment directory not found: {alignment_directory}")
        append_exception_traceback(module_log_path, error)
        raise error

    if not template_pdb_path.exists():
        error = FileNotFoundError(f"Template PDB not found: {template_pdb_path}")
        append_exception_traceback(module_log_path, error)
        raise error

    output_dir.mkdir(parents=True, exist_ok=True)

    if final_model_name is None:
        final_model_name = DEFAULT_FINAL_MODEL_FILENAME

    effective_residue_range = _resolve_residue_range_for_filler(
        residue_range=residue_range,
        template_pdb_path=template_pdb_path,
    )
    _log_debug(module_log_path, f"effective_residue_range: {effective_residue_range}")

    resolved_chain_id = resolve_filler_chain_id(
        chain_id=chain_id,
        alignment_directory=alignment_directory,
        template_pdb_path=template_pdb_path,
        residue_range=effective_residue_range,
    )
    _log_debug(module_log_path, f"Resolved chain_id: {resolved_chain_id}")

    effective_residue_range = _derive_per_chain_effective_range(
        residue_range=effective_residue_range,
        chain_id=resolved_chain_id,
        template_pdb_path=template_pdb_path,
    )

    script_path = output_dir / DEFAULT_MODELLER_SCRIPT_FILENAME

    alignment_fasta_path = find_alignment_fasta_for_filler(
        alignment_directory=alignment_directory,
        chain_id=resolved_chain_id,
    )

    copied_template_pdb = write_chain_specific_template_pdb(
        template_pdb_path=template_pdb_path,
        output_dir=output_dir,
        template_id=template_id,
        chain_id=resolved_chain_id,
    )

    (
        template_header,
        template_alignment_skeleton,
        target_header,
        _target_aligned_sequence,
    ) = split_template_and_target_alignment_records(alignment_fasta_path)

    header_chain_id = _extract_chain_id_from_alignment_header(template_header)
    normalized_chain_id = _normalize_chain_id(resolved_chain_id)
    if header_chain_id is not None and header_chain_id != normalized_chain_id:
        error = ValueError(
            f"Alignment FASTA template header chain mismatch: "
            f"expected chain {normalized_chain_id!r}, "
            f"but header {template_header!r} suggests chain {header_chain_id!r}."
        )
        append_exception_traceback(module_log_path, error)
        raise error

    fill_decision = analyze_fill_decision_from_template_alignment(
        template_alignment_skeleton=template_alignment_skeleton,
    )

    append_key_value_block(
        log_path=module_log_path,
        title="FILL DECISION",
        key_value_pairs=[
            ("resolved_chain_id", resolved_chain_id),
            ("effective_residue_range", effective_residue_range),
            ("should_run_modeller", str(fill_decision.should_run_modeller)),
            ("overall_classification", fill_decision.overall_classification),
            ("alphafold_candidate", str(fill_decision.alphafold_candidate)),
            ("n_gap_regions", str(len(fill_decision.gap_regions))),
            ("skip_reason", fill_decision.skip_reason or "(none)"),
        ],
    )

    _print_decision_table(
        pdb_id=pdb_id,
        resolved_chain_id=resolved_chain_id,
        fill_decision=fill_decision,
    )
    _print_gap_regions_table(fill_decision)

    modeller_model_path: Path | None = None
    alphafold_model_path: Path | None = None
    raw_model_paths: tuple[Path, ...] = ()
    alignment_file = output_dir / DEFAULT_ALIGNMENT_FILENAME

    try:
        if (
            skip_if_no_internal_gaps
            and fill_decision.overall_classification == "no_internal_gaps"
        ):
            stdout_log, stderr_log = _write_skip_logs(
                output_dir=output_dir,
                message=(
                    "Filler skipped: no internal gaps in template alignment. "
                    f"Reason: {fill_decision.skip_reason}"
                ),
            )
            append_log_text(module_log_path, "[RESULT] skipped: no_internal_gaps")
            result = _build_result(
                chain_id=resolved_chain_id,
                output_dir=output_dir,
                alignment_file=alignment_file,
                template_pdb=copied_template_pdb,
                script_file=script_path,
                modeller_model_path=None,
                alphafold_model_path=None,
                final_model_path=None,
                raw_model_paths=(),
                stdout_log=stdout_log,
                stderr_log=stderr_log,
                skipped=True,
                skip_reason=fill_decision.skip_reason,
                fill_decision=fill_decision,
            )
            _print_end_box(
                pdb_id=pdb_id,
                resolved_chain_id=resolved_chain_id,
                status="skipped",
                final_model_path=None,
                modeller_model_path=None,
                alphafold_model_path=None,
                stdout_log=stdout_log,
                stderr_log=stderr_log,
                skip_reason=fill_decision.skip_reason,
            )
            return result

        if fill_decision.alphafold_candidate:
            if not uniprot_id:
                stdout_log, stderr_log = _write_skip_logs(
                    output_dir=output_dir,
                    message=(
                        "Filler skipped: AlphaFold fallback required, but UniProt ID is missing."
                    ),
                )
                append_log_text(
                    module_log_path, "[RESULT] skipped: missing_uniprot_for_alphafold"
                )
                result = _build_result(
                    chain_id=resolved_chain_id,
                    output_dir=output_dir,
                    alignment_file=alignment_file,
                    template_pdb=copied_template_pdb,
                    script_file=script_path,
                    modeller_model_path=None,
                    alphafold_model_path=None,
                    final_model_path=None,
                    raw_model_paths=(),
                    stdout_log=stdout_log,
                    stderr_log=stderr_log,
                    skipped=True,
                    skip_reason="AlphaFold fallback required, but UniProt ID is missing.",
                    fill_decision=fill_decision,
                )
                _print_end_box(
                    pdb_id=pdb_id,
                    resolved_chain_id=resolved_chain_id,
                    status="skipped",
                    final_model_path=None,
                    modeller_model_path=None,
                    alphafold_model_path=None,
                    stdout_log=stdout_log,
                    stderr_log=stderr_log,
                    skip_reason="AlphaFold fallback required, but UniProt ID is missing.",
                )
                return result

            stdout_log, stderr_log = _write_skip_logs(
                output_dir=output_dir,
                message=(
                    "MODELLER skipped: AlphaFold fallback selected because at least "
                    "one internal gap is 8 residues or longer."
                ),
            )
            mapping_path = (
                alignment_directory
                / f"ATOM_chain_{resolved_chain_id}_vs_UniProt.aln.mapping.tsv"
            )
            uniprot_residue_range = _get_uniprot_range_from_mapping(mapping_path)
            pdb_to_uniprot_map = (
                _build_pdb_to_uniprot_residue_map(copied_template_pdb, mapping_path)
                if mapping_path.is_file()
                else {}
            )
            if uniprot_residue_range:
                _log_debug(
                    module_log_path,
                    f"AlphaFold crop: using UniProt range {uniprot_residue_range!r} "
                    f"(PDB range was {effective_residue_range!r})",
                )
            alphafold_model_path = run_alphafold_fallback_for_chain(
                output_dir=output_dir,
                template_pdb_path=copied_template_pdb,
                uniprot_id=uniprot_id,
                residue_range=effective_residue_range,
                final_model_name=final_model_name,
                model_version=4,
                uniprot_residue_range=uniprot_residue_range,
                pdb_to_uniprot_map=pdb_to_uniprot_map or None,
            )
            trim_range = uniprot_residue_range or effective_residue_range
            _trim_final_model_in_place(alphafold_model_path, trim_range)
            result = _build_result(
                chain_id=resolved_chain_id,
                output_dir=output_dir,
                alignment_file=alignment_file,
                template_pdb=copied_template_pdb,
                script_file=script_path,
                modeller_model_path=None,
                alphafold_model_path=alphafold_model_path,
                final_model_path=alphafold_model_path,
                raw_model_paths=(),
                stdout_log=stdout_log,
                stderr_log=stderr_log,
                skipped=False,
                skip_reason=None,
                fill_decision=fill_decision,
            )
            _print_end_box(
                pdb_id=pdb_id,
                resolved_chain_id=resolved_chain_id,
                status="success_alphafold",
                final_model_path=alphafold_model_path,
                modeller_model_path=None,
                alphafold_model_path=alphafold_model_path,
                stdout_log=stdout_log,
                stderr_log=stderr_log,
                skip_reason=None,
            )
            return result

        if not fill_decision.should_run_modeller:
            stdout_log, stderr_log = _write_skip_logs(
                output_dir=output_dir,
                message=(
                    "MODELLER not run for this chain. "
                    f"Reason: {fill_decision.skip_reason or 'unspecified'}"
                ),
            )
            append_log_text(module_log_path, "[RESULT] skipped: modeller_not_run")
            result = _build_result(
                chain_id=resolved_chain_id,
                output_dir=output_dir,
                alignment_file=alignment_file,
                template_pdb=copied_template_pdb,
                script_file=script_path,
                modeller_model_path=None,
                alphafold_model_path=None,
                final_model_path=None,
                raw_model_paths=(),
                stdout_log=stdout_log,
                stderr_log=stderr_log,
                skipped=True,
                skip_reason=fill_decision.skip_reason,
                fill_decision=fill_decision,
            )
            _print_end_box(
                pdb_id=pdb_id,
                resolved_chain_id=resolved_chain_id,
                status="skipped",
                final_model_path=None,
                modeller_model_path=None,
                alphafold_model_path=None,
                stdout_log=stdout_log,
                stderr_log=stderr_log,
                skip_reason=fill_decision.skip_reason,
            )
            return result

        alignment_file = write_modeller_alignment_from_existing_alignment(
            alignment_fasta_path=alignment_fasta_path,
            template_pdb_path=copied_template_pdb,
            output_dir=output_dir,
            template_id=template_id,
            target_id=target_id,
            chain_id=resolved_chain_id,
        )
        validate_modeller_inputs(
            output_dir=output_dir,
            alignment_file=alignment_file,
            template_id=template_id,
        )
        components_dir = template_pdb_path.parent
        contact_atoms = find_metal_and_ligand_contact_atoms(
            protein_pdb_path=copied_template_pdb,
            metals_pdb_path=components_dir / f"{pdb_id}_metals.pdb",
            ligand_pdb_path=components_dir / f"{pdb_id}_ligand.pdb",
            cofactor_pdb_path=components_dir / f"{pdb_id}_cofactor.pdb",
        )
        if contact_atoms:
            _log_debug(
                module_log_path,
                f"Adding position restraints for {len(contact_atoms)} "
                f"metal/ligand contact atoms",
            )
        script_path = write_modeller_script(
            output_dir=output_dir,
            alignment_file=alignment_file,
            template_id=template_id,
            target_id=target_id,
            starting_model=starting_model,
            ending_model=ending_model,
            contact_atoms=contact_atoms or None,
        )
        stdout_log, stderr_log = run_modeller_binary(
            script_path=script_path, working_dir=output_dir
        )
        raw_model_paths = find_raw_models(output_dir=output_dir, target_id=target_id)
        _best_raw_model_path, modeller_model_path = standardize_model_name(
            output_dir=output_dir,
            target_id=target_id,
            raw_model_paths=raw_model_paths,
            final_name=final_model_name,
        )
        _trim_final_model_in_place(modeller_model_path, effective_residue_range)
        append_log_text(module_log_path, "[RESULT] success: modeller")

        result = _build_result(
            chain_id=resolved_chain_id,
            output_dir=output_dir,
            alignment_file=alignment_file,
            template_pdb=copied_template_pdb,
            script_file=script_path,
            modeller_model_path=modeller_model_path,
            alphafold_model_path=None,
            final_model_path=modeller_model_path,
            raw_model_paths=raw_model_paths,
            stdout_log=stdout_log,
            stderr_log=stderr_log,
            skipped=False,
            skip_reason=None,
            fill_decision=fill_decision,
        )
        _print_end_box(
            pdb_id=pdb_id,
            resolved_chain_id=resolved_chain_id,
            status="success_modeller",
            final_model_path=modeller_model_path,
            modeller_model_path=modeller_model_path,
            alphafold_model_path=None,
            stdout_log=stdout_log,
            stderr_log=stderr_log,
            skip_reason=None,
        )
        return result

    except Exception as error:
        append_log_text(module_log_path, "[RESULT] FAILED")
        append_exception_traceback(module_log_path, error)
        print_double_header_box(
            title=f"{MODULE_NAME} · END",
            line_list=[
                f"PDB ID           : {pdb_id}",
                f"Resolved chain   : {resolved_chain_id if 'resolved_chain_id' in locals() else '(unresolved)'}",
                "Status           : failed",
                f"Error type       : {type(error).__name__}",
                f"Error            : {error}",
                f"Module log       : {module_log_path}",
            ],
        )
        raise
