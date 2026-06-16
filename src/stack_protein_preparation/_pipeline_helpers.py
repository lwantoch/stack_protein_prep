"""Pure helper functions for the pipeline runner — no _log dependency."""
from __future__ import annotations

import re
from pathlib import Path
from typing import Any

from stack_protein_preparation.parameter_audit import audit_protein_residue_parameters
from stack_protein_preparation.pipeline_paths import (
    _build_monomer_protein_path,
    _build_raw_protein_path,
    _build_variant_parameter_audit_log_path,
    _build_variant_parameter_audit_manifest_path,
)
from stack_protein_preparation.pipeline_state import (
    METALS_CHECK_LOG_PATH_COLUMN_NAME,
    METALS_CHECK_MANIFEST_PATH_COLUMN_NAME,
    METALS_CHECK_STATUS_COLUMN_NAME,
    METALS_CLASS_COLUMN_NAME,
    METALS_GEOMETRY_FOUND_COLUMN_NAME,
    METALS_GEOMETRY_MATCH_COLUMN_NAME,
    METALS_GEOMETRY_PROBABLE_COLUMN_NAME,
    METALS_ION_TYPE_COLUMN_NAME,
    METALS_MODEL_READY_COLUMN_NAME,
    METALS_PARAMETER_REFERENCE_COLUMN_NAME,
    PREPARED_ALPHAFOLD_OUTPUT_PATH_COLUMN_NAME,
    PREPARED_GAPS_OUTPUT_PATH_COLUMN_NAME,
    PREPARED_MODELLER_OUTPUT_PATH_COLUMN_NAME,
    PREPARED_STRUCTURE_ACTUAL_RANGE_COLUMN_NAME,
    PREPARED_STRUCTURE_OUTPUT_PATH_COLUMN_NAME,
    PREPARED_STRUCTURE_PROTEIN_INPUT_PATH_COLUMN_NAME,
    PREPARED_STRUCTURE_STATUS_COLUMN_NAME,
    PREPARED_STRUCTURE_VARIANT_COLUMN_NAME,
    STATUS_FAILED,
    STATUS_SUCCESS,
)
from stack_protein_preparation.pipeline_variants import (
    _choose_representative_successful_variant,
    _select_default_protein_input_path,
    _update_available_models_field,
)
from stack_protein_preparation.prepared_structure import build_prepared_structure_for_variant
from stack_protein_preparation import _pipeline_config as _cfg


def _path_list_to_text(path_list: list[Path]) -> str:
    return ";".join(str(path) for path in path_list)


def _extract_fragment_paths_from_split_result(split_result: Any) -> list[Path]:
    if isinstance(split_result, dict):
        raw_fragment_paths = split_result.get("fragment_paths", [])
    else:
        raw_fragment_paths = getattr(split_result, "fragment_paths", [])
    fragment_path_list: list[Path] = []
    for raw_path in raw_fragment_paths:
        fragment_path = Path(raw_path)
        if fragment_path.exists():
            fragment_path_list.append(fragment_path)
    return fragment_path_list


def _clear_component_output_files(pdb_id: str, components_dir: Path) -> None:
    for suffix in (
        "protein",
        "protein_monomer",
        "ligand",
        "cofactor",
        "water",
        "metals",
        "nonstandard_residues",
    ):
        (components_dir / f"{pdb_id}_{suffix}.pdb").unlink(missing_ok=True)


def _extract_monomer_output_path(result: Any, fallback_path: Path) -> Path:
    if isinstance(result, str | Path):
        return Path(result)
    output_pdb = getattr(result, "output_pdb", None)
    if output_pdb is not None:
        return Path(output_pdb)
    return fallback_path


def _copy_representative_protein_component_to_monomer_path(
    pdb_id: str,
    pdb_dir: Path,
) -> Path:
    split_protein_path = _build_raw_protein_path(pdb_id=pdb_id, pdb_dir=pdb_dir)
    monomer_protein_path = _build_monomer_protein_path(pdb_id=pdb_id, pdb_dir=pdb_dir)

    if not split_protein_path.exists() or split_protein_path.stat().st_size == 0:
        raise RuntimeError(
            f"Representative protein component was not written: {split_protein_path}"
        )

    monomer_protein_path.write_text(
        split_protein_path.read_text(encoding="utf-8"),
        encoding="utf-8",
    )

    if not monomer_protein_path.exists() or monomer_protein_path.stat().st_size == 0:
        raise RuntimeError(
            f"Representative monomer protein output was not written: {monomer_protein_path}"
        )

    return monomer_protein_path


def _find_uniprot_id_for_protein(pdb_dir: Path) -> str:
    fasta_dir = pdb_dir / "fasta"
    if not fasta_dir.exists():
        return ""
    for fasta_path in sorted(fasta_dir.glob("UniProt_*.fasta")):
        match = re.match(r"^UniProt_([A-Z0-9]+)\.fasta$", fasta_path.name)
        if match:
            return match.group(1)
    return ""


def _find_template_pdb_for_filler(
    pdb_id: str,
    pdb_dir: Path,
    components_dir: Path,
) -> Path | None:
    _ = components_dir
    try:
        return _select_default_protein_input_path(pdb_id=pdb_id, pdb_dir=pdb_dir)
    except FileNotFoundError:
        return None


def _get_sanitize_result_output_path(
    sanitize_result: Any,
    fallback_path: Path,
) -> Path:
    raw_output_path = str(getattr(sanitize_result, "output_pdb_path", "")).strip()
    if raw_output_path:
        output_path = Path(raw_output_path)
        if output_path.exists() and output_path.stat().st_size > 0:
            return output_path
    return fallback_path


def _summarize_sanitize_issues(sanitize_result: Any) -> str:
    issue_list = list(getattr(sanitize_result, "issues", ()) or ())
    if not issue_list:
        return "<none>"
    seen: dict[str, int] = {}
    for issue in issue_list:
        key = f"{getattr(issue, 'severity', '')}:{getattr(issue, 'code', '')}"
        seen[key] = seen.get(key, 0) + 1
    parts = []
    for key, count in seen.items():
        parts.append(key if count == 1 else f"{key} ×{count}")
    return " | ".join(parts)


def _run_parameter_audit_for_protein(
    *,
    pdb_id: str,
    pdb_dir: Path,
    input_pdb_path: Path,
) -> Any:
    audit_manifest_path = pdb_dir / "parameter_audit.json"
    audit_log_path = pdb_dir.parent / "logs" / pdb_id / "parameter_audit.log"
    audit_log_path.parent.mkdir(parents=True, exist_ok=True)
    return audit_protein_residue_parameters(
        input_pdb_path=input_pdb_path,
        pdb_id=pdb_id,
        force_field=_cfg._force_field,
        water_model=_cfg._water_model,
        manifest_path=audit_manifest_path,
        log_path=audit_log_path,
    )


def _run_parameter_audit_for_variant(
    *,
    pdb_id: str,
    pdb_dir: Path,
    variant_label: str,
    input_pdb_path: Path,
) -> Any:
    audit_manifest_path = _build_variant_parameter_audit_manifest_path(
        pdb_id=pdb_id,
        pdb_dir=pdb_dir,
        variant_label=variant_label,
    )
    audit_log_path = _build_variant_parameter_audit_log_path(
        pdb_id=pdb_id,
        pdb_dir=pdb_dir,
        variant_label=variant_label,
    )
    audit_log_path.parent.mkdir(parents=True, exist_ok=True)
    return audit_protein_residue_parameters(
        input_pdb_path=input_pdb_path,
        pdb_id=pdb_id,
        force_field=_cfg._force_field,
        water_model=_cfg._water_model,
        manifest_path=audit_manifest_path,
        log_path=audit_log_path,
    )


def _store_metal_check_result_in_record(
    *,
    pipeline_record: dict[str, str],
    metal_check_result: Any,
) -> None:
    pipeline_record[METALS_CHECK_STATUS_COLUMN_NAME] = str(getattr(metal_check_result, "status", ""))
    pipeline_record[METALS_ION_TYPE_COLUMN_NAME] = str(getattr(metal_check_result, "ion_type_text", ""))
    pipeline_record[METALS_CLASS_COLUMN_NAME] = str(getattr(metal_check_result, "metal_class", ""))
    pipeline_record[METALS_PARAMETER_REFERENCE_COLUMN_NAME] = str(getattr(metal_check_result, "parameter_reference_text", ""))
    pipeline_record[METALS_GEOMETRY_FOUND_COLUMN_NAME] = str(getattr(metal_check_result, "found_geometry_text", ""))
    pipeline_record[METALS_GEOMETRY_PROBABLE_COLUMN_NAME] = str(getattr(metal_check_result, "probable_geometry_text", ""))
    pipeline_record[METALS_GEOMETRY_MATCH_COLUMN_NAME] = str(getattr(metal_check_result, "geometry_match", ""))
    pipeline_record[METALS_MODEL_READY_COLUMN_NAME] = "yes" if bool(getattr(metal_check_result, "model_metal_ready", False)) else "no"
    pipeline_record[METALS_CHECK_LOG_PATH_COLUMN_NAME] = str(getattr(metal_check_result, "log_path", ""))
    pipeline_record[METALS_CHECK_MANIFEST_PATH_COLUMN_NAME] = str(getattr(metal_check_result, "manifest_path", ""))


def _variant_metal_check_accepts_end_model(variant_result: dict[str, str]) -> bool:
    value = str(variant_result.get("metals.model_ready", "")).strip().lower()
    return value in {"", "yes", "true", "1"}


def _variant_audit_accepts_end_model(audit_result: Any) -> bool:
    if not bool(getattr(audit_result, "can_use_current_force_field", False)):
        return False
    if bool(getattr(audit_result, "requires_repair", False)):
        return False
    if bool(getattr(audit_result, "requires_external_parameters", False)):
        return False
    if bool(getattr(audit_result, "requires_qm_parameters", False)):
        return False
    if bool(getattr(audit_result, "requires_metal_parameters", False)):
        return False
    return True


def _compute_actual_protein_range(pdb_path: Path) -> str:
    resnums: list[int] = []
    try:
        for line in pdb_path.read_text(encoding="utf-8", errors="replace").splitlines():
            if not line.startswith("ATOM  "):
                continue
            try:
                resnums.append(int(line[22:26]))
            except (ValueError, IndexError):
                pass
    except OSError:
        return ""
    if not resnums:
        return ""
    return f"{min(resnums)}-{max(resnums)}"


def _write_dockingbox_csv(
    pdb_id: str,
    pdb_dir: Path,
    variant_label: str,
) -> Path | None:
    from stack_protein_preparation.prepared_structure import sanitize_variant_label

    ligand_path = pdb_dir / "components" / f"{pdb_id}_ligand.pdb"
    if not ligand_path.exists() or ligand_path.stat().st_size == 0:
        return None

    coords: list[tuple[float, float, float]] = []
    for line in ligand_path.read_text(encoding="utf-8").splitlines():
        if not (line.startswith("ATOM  ") or line.startswith("HETATM")):
            continue
        try:
            coords.append((float(line[30:38]), float(line[38:46]), float(line[46:54])))
        except (ValueError, IndexError):
            continue

    if not coords:
        return None

    n = len(coords)
    cx = sum(c[0] for c in coords) / n
    cy = sum(c[1] for c in coords) / n
    cz = sum(c[2] for c in coords) / n

    max_span = 0.0
    for i in range(n):
        for j in range(i + 1, n):
            dx = coords[i][0] - coords[j][0]
            dy = coords[i][1] - coords[j][1]
            dz = coords[i][2] - coords[j][2]
            d = (dx * dx + dy * dy + dz * dz) ** 0.5
            if d > max_span:
                max_span = d

    safe_variant = sanitize_variant_label(variant_label)
    out_dir = pdb_dir / "prepared" / safe_variant
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / "dockingbox.csv"
    out_path.write_text(
        f"center_x,center_y,center_z,max_span\n{cx:.3f},{cy:.3f},{cz:.3f},{max_span:.3f}\n",
        encoding="utf-8",
    )
    return out_path


def _build_prepared_structure_for_variant(
    pdb_id: str,
    pdb_dir: Path,
    variant_label: str,
    final_protein_input_path: Path | None = None,
    final_protein_input_paths: list[Path] | None = None,
) -> Path:
    metals_input_path = pdb_dir / "components" / f"{pdb_id}_metals.pdb"
    kwargs: dict[str, object] = {
        "pdb_directory": pdb_dir,
        "pdb_id": pdb_id,
        "structure_variant": variant_label,
        "water_input_path": None,
        "ligand_input_path": None,
        "metals_input_path": metals_input_path if metals_input_path.exists() else None,
    }
    if final_protein_input_paths:
        kwargs["protein_input_paths"] = final_protein_input_paths
    else:
        if final_protein_input_path is None:
            raise ValueError(
                "Either final_protein_input_path or final_protein_input_paths must be provided."
            )
        kwargs["protein_input_path"] = final_protein_input_path
    summary = build_prepared_structure_for_variant(**kwargs)
    _write_dockingbox_csv(pdb_id, pdb_dir, variant_label)
    return summary.output_pdb_path


def _build_prepared_wat_structure_for_variant(
    pdb_id: str,
    pdb_dir: Path,
    variant_label: str,
    final_protein_input_path: Path | None = None,
    final_protein_input_paths: list[Path] | None = None,
) -> Path | None:
    from stack_protein_preparation.prepared_structure import (
        build_prepared_structure,
        sanitize_variant_label,
    )

    water_path = pdb_dir / "components" / f"{pdb_id}_water.pdb"
    if not water_path.exists() or water_path.stat().st_size == 0:
        return None

    metals_input_path = pdb_dir / "components" / f"{pdb_id}_metals.pdb"
    protein_component_path = pdb_dir / "components" / f"{pdb_id}_protein.pdb"

    safe_variant = sanitize_variant_label(variant_label)
    wat_output_path = pdb_dir / "prepared" / safe_variant / f"{pdb_id}_WAT.pdb"

    kwargs: dict[str, object] = {
        "output_pdb_path": wat_output_path,
        "water_input_path": water_path,
        "ligand_input_path": None,
        "metals_input_path": metals_input_path if metals_input_path.exists() else None,
        "structure_variant": variant_label,
    }
    if final_protein_input_paths:
        kwargs["protein_input_paths"] = final_protein_input_paths
    else:
        if final_protein_input_path is None:
            raise ValueError(
                "Either final_protein_input_path or final_protein_input_paths must be provided."
            )
        kwargs["protein_input_path"] = final_protein_input_path
    summary = build_prepared_structure(**kwargs)
    return summary.output_pdb_path


def _build_prepared_cofactor_structure_for_variant(
    pdb_id: str,
    pdb_dir: Path,
    variant_label: str,
    final_protein_input_path: Path | None = None,
    final_protein_input_paths: list[Path] | None = None,
) -> Path | None:
    from stack_protein_preparation.prepared_structure import (
        build_prepared_structure,
        sanitize_variant_label,
    )

    cofactor_path = pdb_dir / "components" / f"{pdb_id}_cofactor.pdb"
    if not cofactor_path.exists() or cofactor_path.stat().st_size == 0:
        return None

    protein_component_path = pdb_dir / "components" / f"{pdb_id}_protein.pdb"
    safe_variant = sanitize_variant_label(variant_label)
    cof_output_path = pdb_dir / "prepared" / safe_variant / f"{pdb_id}_cofactor.pdb"

    kwargs: dict[str, object] = {
        "output_pdb_path": cof_output_path,
        "water_input_path": None,
        "ligand_input_path": cofactor_path,
        "metals_input_path": None,
        "structure_variant": variant_label,
    }
    if final_protein_input_paths:
        kwargs["protein_input_paths"] = final_protein_input_paths
    else:
        if final_protein_input_path is None:
            raise ValueError(
                "Either final_protein_input_path or final_protein_input_paths must be provided."
            )
        kwargs["protein_input_path"] = final_protein_input_path
    summary = build_prepared_structure(**kwargs)
    return summary.output_pdb_path


def _apply_accepted_variant_audit_decision(
    *,
    pipeline_record: dict[str, str],
    accepted_variant_result_list: list[dict[str, str]],
) -> None:
    from stack_protein_preparation.pipeline_records import _clear_prepared_variant_output_paths

    _clear_prepared_variant_output_paths(pipeline_record)

    if not accepted_variant_result_list:
        pipeline_record[PREPARED_STRUCTURE_STATUS_COLUMN_NAME] = STATUS_FAILED
        return

    representative_variant_result = _choose_representative_successful_variant(
        accepted_variant_result_list
    )

    accepted_label_list = [
        result["variant_label"] for result in accepted_variant_result_list
    ]
    pipeline_record[PREPARED_STRUCTURE_STATUS_COLUMN_NAME] = STATUS_SUCCESS
    pipeline_record[PREPARED_STRUCTURE_VARIANT_COLUMN_NAME] = ",".join(accepted_label_list)
    pipeline_record[PREPARED_STRUCTURE_PROTEIN_INPUT_PATH_COLUMN_NAME] = (
        representative_variant_result.get("final_protein_input_paths", "")
    )
    _representative_output_path = representative_variant_result.get("prepared_output_path", "")
    pipeline_record[PREPARED_STRUCTURE_OUTPUT_PATH_COLUMN_NAME] = _representative_output_path

    if _representative_output_path:
        _actual_range = _compute_actual_protein_range(Path(_representative_output_path))
        if _actual_range:
            pipeline_record[PREPARED_STRUCTURE_ACTUAL_RANGE_COLUMN_NAME] = _actual_range

    for result in accepted_variant_result_list:
        variant_label = result.get("variant_label", "")
        model_source = result.get("model_source", "")
        prepared_output_path = result.get("prepared_output_path", "")

        if variant_label == "gaps":
            pipeline_record[PREPARED_GAPS_OUTPUT_PATH_COLUMN_NAME] = prepared_output_path
        elif variant_label in {"best_complete", "large_gap_complete"}:
            if model_source == "modeller":
                pipeline_record[PREPARED_MODELLER_OUTPUT_PATH_COLUMN_NAME] = prepared_output_path
            elif model_source == "alphafold":
                pipeline_record[PREPARED_ALPHAFOLD_OUTPUT_PATH_COLUMN_NAME] = prepared_output_path

    for metal_column_name in [
        METALS_CHECK_STATUS_COLUMN_NAME,
        METALS_ION_TYPE_COLUMN_NAME,
        METALS_CLASS_COLUMN_NAME,
        METALS_PARAMETER_REFERENCE_COLUMN_NAME,
        METALS_GEOMETRY_FOUND_COLUMN_NAME,
        METALS_GEOMETRY_PROBABLE_COLUMN_NAME,
        METALS_GEOMETRY_MATCH_COLUMN_NAME,
        METALS_MODEL_READY_COLUMN_NAME,
        METALS_CHECK_LOG_PATH_COLUMN_NAME,
        METALS_CHECK_MANIFEST_PATH_COLUMN_NAME,
    ]:
        if representative_variant_result.get(metal_column_name, ""):
            pipeline_record[metal_column_name] = representative_variant_result.get(metal_column_name, "")

    _update_available_models_field(pipeline_record)
