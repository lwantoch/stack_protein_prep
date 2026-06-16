"""Stage-execution functions for the pipeline runner — all use _cfg._log."""
from __future__ import annotations

import traceback
from pathlib import Path
from typing import Any

from stack_protein_preparation.metalls_check import (
    run_metals_check_for_model,
    run_metals_inventory_for_structure,
)
from stack_protein_preparation.monomer import write_representative_monomer_unit
from stack_protein_preparation.pipeline_paths import (
    _build_representative_unit_path,
    _build_representative_sanitize_log_path,
    _build_sanitize_log_path,
    _build_sanitized_monomer_protein_path,
    _build_sanitized_prepared_variant_path,
    _build_sanitized_protonation_path,
    _build_standard_residue_repair_log_path,
    _build_standard_residue_repair_path,
    _build_prepared_variant_sanitize_log_path,
)
from stack_protein_preparation.pipeline_state import (
    STATUS_FAILED,
    STATUS_SKIPPED,
    STATUS_SUCCESS,
)
from stack_protein_preparation.pipeline_variants import _choose_variant_structure_input_path
from stack_protein_preparation.protonation import protonate_variant_structure
from stack_protein_preparation.sanitize import sanitize_pdb_for_gromacs
from stack_protein_preparation.standard_residue_repair import (
    repair_standard_residue_heavy_atoms_for_gromacs,
)
from stack_protein_preparation.fragment_split import split_gap_variant_structure
from stack_protein_preparation import _pipeline_config as _cfg
from stack_protein_preparation._pipeline_helpers import (
    _extract_fragment_paths_from_split_result,
    _extract_monomer_output_path,
    _path_list_to_text,
    _summarize_sanitize_issues,
)


def _write_representative_monomer_for_protein(
    pdb_id: str,
    pdb_dir: Path,
    input_pdb_path: Path,
) -> Path:
    components_dir = Path(pdb_dir) / "components"
    representative_unit_path = _build_representative_unit_path(
        pdb_id=pdb_id,
        pdb_dir=pdb_dir,
    )

    if not input_pdb_path.exists():
        raise FileNotFoundError(
            f"Representative-unit input file not found: {input_pdb_path}"
        )

    components_dir.mkdir(parents=True, exist_ok=True)

    result = write_representative_monomer_unit(
        input_pdb_path,
        representative_unit_path,
        keep_non_protein_hetero=True,
        keep_waters=True,
        overwrite=True,
    )

    result_path = _extract_monomer_output_path(
        result=result,
        fallback_path=representative_unit_path,
    )

    if not result_path.exists() or result_path.stat().st_size == 0:
        raise RuntimeError(
            f"Representative unit output was not written or is empty: {result_path}"
        )

    representative_chain_ids = getattr(result, "representative_chain_ids", ())
    represented_chain_groups = getattr(result, "represented_chain_groups", ())

    _cfg._log(
        "monomer_selection:success",
        [
            f"pdb_id                     : {pdb_id}",
            f"representative_unit_input  : {input_pdb_path}",
            f"representative_unit_output : {result_path}",
            f"representative_chain_ids   : {'/'.join(representative_chain_ids)}",
            f"represented_chain_groups   : {represented_chain_groups}",
            "ligand_policy              : keep non-water HETATM only for selected representative chains",
        ],
    )

    return result_path


def _append_parameter_audit_summary(
    *,
    pdb_id: str,
    audit_result: Any,
) -> None:
    _cfg._log(
        f"parameter_audit:{pdb_id}",
        [
            f"  {pdb_id}",
            f"  status      :  {getattr(audit_result, 'status', '')}  "
            f"│  gromacs_probe: {getattr(audit_result, 'gromacs_probe_success', '')}  "
            f"can_use_ff: {getattr(audit_result, 'can_use_current_force_field', '')}",
            f"  flags       :  repair: {getattr(audit_result, 'requires_repair', '')}  "
            f"ext_params: {getattr(audit_result, 'requires_external_parameters', '')}  "
            f"qm: {getattr(audit_result, 'requires_qm_parameters', '')}  "
            f"metal: {getattr(audit_result, 'requires_metal_parameters', '')}",
            f"  supported   :  {', '.join(getattr(audit_result, 'supported_nonstandard_residue_names', ())) or '(none)'}",
            f"  unsupported :  {', '.join(getattr(audit_result, 'unsupported_residue_names', ())) or '(none)'}",
        ],
    )


def _run_standard_residue_repair_for_protonation_input(
    *,
    pdb_id: str,
    pdb_dir: Path,
    variant_label: str,
    input_pdb: Path,
) -> Any:
    repair_result = repair_standard_residue_heavy_atoms_for_gromacs(
        input_pdb_path=input_pdb,
        output_pdb_path=_build_standard_residue_repair_path(
            pdb_id=pdb_id,
            pdb_dir=pdb_dir,
            variant_label=variant_label,
        ),
        log_path=_build_standard_residue_repair_log_path(
            pdb_id=pdb_id,
            pdb_dir=pdb_dir,
            variant_label=variant_label,
        ),
        force_field=_cfg._force_field,
        allow_nonstandard_residues=False,
    )
    _cfg._log(
        "standard_residue_repair",
        [
            f"  {pdb_id}  /  {variant_label}",
            f"  status    :  {getattr(repair_result, 'status', '')}  "
            f"│  attempted: {getattr(repair_result, 'repair_attempted', '')}  "
            f"success: {getattr(repair_result, 'repair_success', '')}",
            f"  missing   :  before: {getattr(repair_result, 'missing_heavy_atom_count_before', '')}  "
            f"after: {getattr(repair_result, 'missing_heavy_atom_count_after', '')}",
            f"  message   :  {getattr(repair_result, 'message', '')}",
            f"  output    :  {getattr(repair_result, 'used_output_path', '')}",
        ],
    )
    return repair_result


def _run_sanitize_for_protonation_input(
    *,
    pdb_id: str,
    pdb_dir: Path,
    variant_label: str,
    input_pdb: Path,
) -> Any:
    repair_result = _run_standard_residue_repair_for_protonation_input(
        pdb_id=pdb_id,
        pdb_dir=pdb_dir,
        variant_label=variant_label,
        input_pdb=input_pdb,
    )
    repair_input_path = Path(str(getattr(repair_result, "used_output_path", input_pdb)))

    output_pdb_path = _build_sanitized_protonation_path(
        pdb_id=pdb_id,
        pdb_dir=pdb_dir,
        variant_label=variant_label,
    )
    log_path = _build_sanitize_log_path(
        pdb_id=pdb_id,
        pdb_dir=pdb_dir,
        variant_label=variant_label,
    )

    sanitize_result = sanitize_pdb_for_gromacs(
        input_pdb_path=repair_input_path,
        output_pdb_path=output_pdb_path,
        pdb_id=pdb_id,
        variant_label=variant_label,
        force_field=_cfg._force_field,
        log_path=log_path,
        fail_on_missing_heavy_atoms=False,
        fail_on_nonstandard_residues=False,
    )

    _cfg._log(
        "sanitize:protonation_input",
        [
            f"  {pdb_id}  /  {variant_label}",
            f"  status       :  {getattr(sanitize_result, 'status', '')}  "
            f"(gromacs_ok: {getattr(sanitize_result, 'can_run_gromacs', '')})",
            f"  atoms        :  {getattr(sanitize_result, 'atom_count_in', '')} → {getattr(sanitize_result, 'atom_count_out', '')}"
            f"  │  altloc: {getattr(sanitize_result, 'selected_altloc_count', '')}  occ_norm: {getattr(sanitize_result, 'normalized_occupancy_count', '')}",
            f"  heavy atoms  :  missing: {getattr(sanitize_result, 'missing_heavy_atom_count', '')}  repaired: {getattr(sanitize_result, 'repaired_heavy_atom_count', 0)}  blocked: {getattr(sanitize_result, 'blocked_repair_count', 0)}",
            f"  nonstandard  :  {', '.join(getattr(sanitize_result, 'nonstandard_residue_names', ())) or '(none)'}",
            f"  issues       :  {_summarize_sanitize_issues(sanitize_result)}",
            f"  input        :  {input_pdb}",
            f"  output       :  {getattr(sanitize_result, 'output_pdb_path', '')}",
        ],
    )

    return sanitize_result


def _run_sanitize_for_representative_protein(
    *,
    pdb_id: str,
    pdb_dir: Path,
    input_pdb: Path,
) -> Any:
    output_pdb_path = _build_sanitized_monomer_protein_path(
        pdb_id=pdb_id,
        pdb_dir=pdb_dir,
    )
    log_path = _build_representative_sanitize_log_path(
        pdb_id=pdb_id,
        pdb_dir=pdb_dir,
    )

    sanitize_result = sanitize_pdb_for_gromacs(
        input_pdb_path=input_pdb,
        output_pdb_path=output_pdb_path,
        pdb_id=pdb_id,
        variant_label="representative_protein",
        force_field=_cfg._force_field,
        log_path=log_path,
        fail_on_missing_heavy_atoms=False,
        fail_on_nonstandard_residues=False,
    )

    _cfg._log(
        "sanitize:representative_protein",
        [
            f"  {pdb_id}  /  representative_protein",
            f"  status       :  {getattr(sanitize_result, 'status', '')}  "
            f"(gromacs_ok: {getattr(sanitize_result, 'can_run_gromacs', '')})",
            f"  atoms        :  {getattr(sanitize_result, 'atom_count_in', '')} → {getattr(sanitize_result, 'atom_count_out', '')}"
            f"  │  altloc: {getattr(sanitize_result, 'selected_altloc_count', '')}  occ_norm: {getattr(sanitize_result, 'normalized_occupancy_count', '')}",
            f"  heavy atoms  :  missing: {getattr(sanitize_result, 'missing_heavy_atom_count', '')}  repaired: {getattr(sanitize_result, 'repaired_heavy_atom_count', 0)}  blocked: {getattr(sanitize_result, 'blocked_repair_count', 0)}",
            f"  nonstandard  :  {', '.join(getattr(sanitize_result, 'nonstandard_residue_names', ())) or '(none)'}",
            f"  issues       :  {_summarize_sanitize_issues(sanitize_result)}",
            f"  input        :  {input_pdb}",
            f"  output       :  {getattr(sanitize_result, 'output_pdb_path', '')}",
        ],
    )

    return sanitize_result


def _run_sanitize_for_prepared_variant(
    *,
    pdb_id: str,
    pdb_dir: Path,
    variant_label: str,
    prepared_output_path: Path,
) -> Any:
    output_pdb_path = _build_sanitized_prepared_variant_path(
        pdb_id=pdb_id,
        pdb_dir=pdb_dir,
        variant_label=variant_label,
    )
    log_path = _build_prepared_variant_sanitize_log_path(
        pdb_id=pdb_id,
        pdb_dir=pdb_dir,
        variant_label=variant_label,
    )

    sanitize_result = sanitize_pdb_for_gromacs(
        input_pdb_path=prepared_output_path,
        output_pdb_path=output_pdb_path,
        pdb_id=pdb_id,
        variant_label=f"prepared_{variant_label}",
        force_field=_cfg._force_field,
        log_path=log_path,
        fail_on_missing_heavy_atoms=False,
        fail_on_nonstandard_residues=False,
    )

    _cfg._log(
        "sanitize:prepared_variant",
        [
            f"  {pdb_id}  /  {variant_label}",
            f"  status       :  {getattr(sanitize_result, 'status', '')}  "
            f"(gromacs_ok: {getattr(sanitize_result, 'can_run_gromacs', '')})",
            f"  atoms        :  {getattr(sanitize_result, 'atom_count_in', '')} → {getattr(sanitize_result, 'atom_count_out', '')}"
            f"  │  altloc: {getattr(sanitize_result, 'selected_altloc_count', '')}  occ_norm: {getattr(sanitize_result, 'normalized_occupancy_count', '')}",
            f"  heavy atoms  :  missing: {getattr(sanitize_result, 'missing_heavy_atom_count', '')}  repaired: {getattr(sanitize_result, 'repaired_heavy_atom_count', 0)}  blocked: {getattr(sanitize_result, 'blocked_repair_count', 0)}",
            f"  nonstandard  :  {', '.join(getattr(sanitize_result, 'nonstandard_residue_names', ())) or '(none)'}",
            f"  issues       :  {_summarize_sanitize_issues(sanitize_result)}",
            f"  input        :  {prepared_output_path}",
            f"  output       :  {getattr(sanitize_result, 'output_pdb_path', '')}",
        ],
    )

    return sanitize_result


def _append_parameter_audit_variant_summary(
    *,
    pdb_id: str,
    variant_label: str,
    audit_result: Any,
    accepted: bool,
) -> None:
    _accepted_tag = "ACCEPTED" if accepted else "REJECTED"
    _cfg._log(
        f"parameter_audit_variant:{pdb_id}:{variant_label}",
        [
            f"  {pdb_id}  /  {variant_label}  →  {_accepted_tag}",
            f"  status      :  {getattr(audit_result, 'status', '')}  "
            f"│  gromacs_probe: {getattr(audit_result, 'gromacs_probe_success', '')}  "
            f"can_use_ff: {getattr(audit_result, 'can_use_current_force_field', '')}",
            f"  flags       :  repair: {getattr(audit_result, 'requires_repair', '')}  "
            f"ext_params: {getattr(audit_result, 'requires_external_parameters', '')}  "
            f"qm: {getattr(audit_result, 'requires_qm_parameters', '')}  "
            f"metal: {getattr(audit_result, 'requires_metal_parameters', '')}",
            f"  supported   :  {', '.join(getattr(audit_result, 'supported_nonstandard_residue_names', ())) or '(none)'}",
            f"  unsupported :  {', '.join(getattr(audit_result, 'unsupported_residue_names', ())) or '(none)'}",
        ],
    )


def _prepare_crystal_water_component(
    *,
    pdb_id: str,
    pdb_dir: Path,
) -> None:
    from stack_protein_preparation.protonation import add_hydrogens_to_crystal_water_pdb

    water_path = pdb_dir / "components" / f"{pdb_id}_water.pdb"

    if not water_path.exists() or water_path.stat().st_size == 0:
        _cfg._log(
            "crystal_water_preparation:skipped",
            [
                f"pdb_id     : {pdb_id}",
                f"water_path : {water_path}",
                "reason     : water component is empty or absent",
            ],
        )
        return

    add_hydrogens_to_crystal_water_pdb(
        water_pdb_path=water_path,
        ff=_cfg._force_field,
        water_model=_cfg._water_model,
    )

    n_atom_lines = sum(
        1 for line in water_path.read_text(encoding="utf-8").splitlines()
        if line.startswith("ATOM") or line.startswith("HETATM")
    )
    _cfg._log(
        "crystal_water_preparation:success",
        [
            f"pdb_id       : {pdb_id}",
            f"water_path   : {water_path}",
            f"n_atom_lines : {n_atom_lines}",
            "SOL_format   : OW + HW1 + HW2 per molecule",
        ],
    )


def _check_water_clashes_for_prepared_variant(
    *,
    pdb_id: str,
    variant_label: str,
    prepared_output_path: Path,
    hard_clash_dist: float = 2.0,
    soft_clash_dist: float = 2.5,
) -> dict[str, object]:
    import numpy as np

    protein_coords: list[tuple[float, float, float]] = []
    water_records: list[tuple[str, float, float, float]] = []

    if not prepared_output_path.exists():
        return {"n_hard_clashes": 0, "n_soft_clashes": 0, "n_waters_checked": 0}

    with prepared_output_path.open("r", encoding="utf-8") as fh:
        for line in fh:
            try:
                record = line[:6]
                if record == "ATOM  ":
                    elem = line[76:78].strip().upper()
                    if elem not in ("H", "D", ""):
                        protein_coords.append(
                            (float(line[30:38]), float(line[38:46]), float(line[46:54]))
                        )
                elif record == "HETATM":
                    if line[17:20].strip() == "SOL" and line[12:16].strip() == "OW":
                        chain_id = line[21]
                        res_seq = int(line[22:26])
                        icode = line[26].strip()
                        water_id = f"{chain_id}/{res_seq}{icode}"
                        water_records.append(
                            (water_id, float(line[30:38]), float(line[38:46]), float(line[46:54]))
                        )
            except (ValueError, IndexError):
                continue

    if not protein_coords or not water_records:
        _cfg._log(
            f"water_clash_check:{pdb_id}:{variant_label}",
            [
                f"pdb_id             : {pdb_id}",
                f"variant_label      : {variant_label}",
                "result             : skipped (no protein atoms or no waters in prepared structure)",
            ],
        )
        return {"n_hard_clashes": 0, "n_soft_clashes": 0, "n_waters_checked": len(water_records)}

    prot = np.array(protein_coords)
    hard_clashes: list[str] = []
    soft_clashes: list[str] = []

    for water_id, wx, wy, wz in water_records:
        diffs = prot - np.array([wx, wy, wz])
        min_dist = float(np.sqrt((diffs * diffs).sum(axis=1)).min())
        if min_dist < hard_clash_dist:
            hard_clashes.append(f"{water_id}({min_dist:.2f}A)")
        elif min_dist < soft_clash_dist:
            soft_clashes.append(f"{water_id}({min_dist:.2f}A)")

    result: dict[str, object] = {
        "n_waters_checked": len(water_records),
        "n_hard_clashes": len(hard_clashes),
        "n_soft_clashes": len(soft_clashes),
        "hard_clash_examples": hard_clashes[:10],
        "soft_clash_examples": soft_clashes[:10],
    }
    _cfg._log(
        f"water_clash_check:{pdb_id}:{variant_label}",
        [
            f"pdb_id              : {pdb_id}",
            f"variant_label       : {variant_label}",
            f"n_waters_checked    : {len(water_records)}",
            f"n_hard_clashes      : {len(hard_clashes)}  (dist < {hard_clash_dist} A)",
            f"n_soft_clashes      : {len(soft_clashes)}  ({hard_clash_dist} <= dist < {soft_clash_dist} A)",
            f"hard_clash_examples : {', '.join(hard_clashes[:10]) or '(none)'}",
            f"soft_clash_examples : {', '.join(soft_clashes[:10]) or '(none)'}",
        ],
    )
    return result


def _run_metals_check_for_input_inventory(
    *,
    pdb_id: str,
    pdb_dir: Path,
    monomer_path: Path,
) -> Any:
    result = run_metals_inventory_for_structure(
        pdb_id=pdb_id,
        protein_dir=pdb_dir,
        model_pdb_path=monomer_path,
        metals_pdb_path=pdb_dir / "components" / f"{pdb_id}_metals.pdb",
        water_pdb_path=pdb_dir / "components" / f"{pdb_id}_water.pdb",
        output_dir=pdb_dir / "components" / "metals_check" / "input_inventory",
        use_chimera=False,
    )
    _cfg._log(
        "metalls_check:input_inventory",
        [
            f"pdb_id                    : {pdb_id}",
            f"monomer_path              : {monomer_path}",
            f"status                    : {getattr(result, 'status', '')}",
            f"metal_class               : {getattr(result, 'metal_class', '')}",
            f"ion_type_text             : {getattr(result, 'ion_type_text', '')}",
            f"parameter_reference_text  : {getattr(result, 'parameter_reference_text', '')}",
            f"found_geometry_text       : {getattr(result, 'found_geometry_text', '')}",
            f"probable_geometry_text    : {getattr(result, 'probable_geometry_text', '')}",
            f"geometry_match            : {getattr(result, 'geometry_match', '')}",
            f"model_metal_ready         : {getattr(result, 'model_metal_ready', '')}",
            f"contacts_tsv_path         : {getattr(result, 'contacts_tsv_path', '')}",
            f"manifest_path             : {getattr(result, 'manifest_path', '')}",
            f"log_path                  : {getattr(result, 'log_path', '')}",
        ],
    )
    return result


def _run_metals_check_for_prepared_variant(
    *,
    pdb_id: str,
    pdb_dir: Path,
    variant_label: str,
    prepared_output_path: Path,
) -> Any:
    result = run_metals_check_for_model(
        pdb_id=pdb_id,
        protein_dir=pdb_dir,
        variant_label=variant_label,
        model_pdb_path=prepared_output_path,
        metals_pdb_path=pdb_dir / "components" / f"{pdb_id}_metals.pdb",
        water_pdb_path=pdb_dir / "components" / f"{pdb_id}_water.pdb",
        output_dir=pdb_dir / "components" / "metals_check" / variant_label,
        use_chimera=True,
    )
    _cfg._log(
        "metalls_check:prepared_variant",
        [
            f"pdb_id                    : {pdb_id}",
            f"variant_label             : {variant_label}",
            f"prepared_output_path      : {prepared_output_path}",
            f"status                    : {getattr(result, 'status', '')}",
            f"metal_class               : {getattr(result, 'metal_class', '')}",
            f"ion_type_text             : {getattr(result, 'ion_type_text', '')}",
            f"parameter_reference_text  : {getattr(result, 'parameter_reference_text', '')}",
            f"found_geometry_text       : {getattr(result, 'found_geometry_text', '')}",
            f"probable_geometry_text    : {getattr(result, 'probable_geometry_text', '')}",
            f"geometry_match            : {getattr(result, 'geometry_match', '')}",
            f"model_metal_ready         : {getattr(result, 'model_metal_ready', '')}",
            f"contacts_tsv_path         : {getattr(result, 'contacts_tsv_path', '')}",
            f"manifest_path             : {getattr(result, 'manifest_path', '')}",
            f"log_path                  : {getattr(result, 'log_path', '')}",
        ],
    )
    return result


def _run_metall_params_for_protein(
    *,
    pdb_id: str,
    pdb_dir: Path,
) -> dict:
    from stack_protein_preparation.metall_params import run_metal_parametrization_for_protein_dir

    result = run_metal_parametrization_for_protein_dir(protein_dir=pdb_dir)
    _cfg._log(
        "metall_params",
        [
            f"pdb_id                     : {pdb_id}",
            f"status                     : {result.get('status', '')}",
            f"transition_metal_site_count: {result.get('transition_metal_site_count', 0)}",
            f"manifest_path              : {result.get('manifest_path', '')}",
            f"message                    : {result.get('message', '')}",
        ],
    )
    return result


def _run_nonstd_residue_params_for_protein(
    *,
    pdb_id: str,
    pdb_dir: Path,
) -> dict:
    from stack_protein_preparation.nonstd_residue_params import run_nonstd_residue_params

    protein_pdb = pdb_dir / "components" / f"{pdb_id}_protein.pdb"
    result = run_nonstd_residue_params(
        protein_dir=pdb_dir,
        pdb_id=pdb_id,
        protein_pdb=protein_pdb,
    )
    _cfg._log(
        "nonstd_residue_params",
        [
            f"pdb_id         : {pdb_id}",
            f"status         : {result.get('status', '')}",
            f"n_residues     : {result.get('n_residues', 0)}",
            f"manifest_path  : {result.get('manifest_path', '')}",
            f"message        : {result.get('message', '')}",
        ],
    )
    return result


def _run_fruton_protonation_for_variant(
    *,
    pdb_id: str,
    pdb_dir: Path,
    variant_label: str,
    input_pdb: Path,
    input_source: str,
) -> dict[str, str | bool | float | int]:
    return protonate_variant_structure(
        pdb_id=pdb_id,
        protein_dir=pdb_dir,
        variant_label=variant_label,
        input_pdb=input_pdb,
        input_source=input_source,  # type: ignore[arg-type]
        ph=_cfg._ph,
        ff=_cfg._force_field,
        water_model=_cfg._water_model,
    )


def _run_protonation_route_for_variant(
    pdb_id: str,
    pdb_dir: Path,
    variant_label: str,
    model_source: str,
    model_path: Path | None,
) -> dict[str, str]:
    if variant_label == "gaps":
        return _run_gaps_fragment_route(
            pdb_id=pdb_id,
            pdb_dir=pdb_dir,
            variant_label=variant_label,
        )
    return _run_linear_protonation_route_for_variant(
        pdb_id=pdb_id,
        pdb_dir=pdb_dir,
        variant_label=variant_label,
        model_source=model_source,
        model_path=model_path,
    )


def _run_linear_protonation_route_for_variant(
    pdb_id: str,
    pdb_dir: Path,
    variant_label: str,
    model_source: str,
    model_path: Path | None,
) -> dict[str, str]:
    _cfg._log(
        "protonation_route_linear:start",
        [
            f"pdb_id             : {pdb_id}",
            f"variant_label      : {variant_label}",
            f"model_source       : {model_source}",
            f"model_path         : {model_path}",
            "internal_capping    : skipped",
            f"protonation_method : gmx pdb2gmx -ignh",
            f"protonation_ph     : {_cfg._ph} compatibility metadata",
            f"gromacs_ff         : {_cfg._force_field}",
            f"gromacs_water      : {_cfg._water_model}",
            "route              : monomer/model -> sanitize -> gmx pdb2gmx protonation -> prepared_structure",
        ],
    )

    selected_input_path, input_source_label = _choose_variant_structure_input_path(
        pdb_id=pdb_id,
        pdb_dir=pdb_dir,
        model_source=model_source,
        model_path=model_path,
    )

    try:
        sanitize_result = _run_sanitize_for_protonation_input(
            pdb_id=pdb_id,
            pdb_dir=pdb_dir,
            variant_label=variant_label,
            input_pdb=selected_input_path,
        )
    except Exception as error:
        _cfg._log(
            f"protonation_route_linear:sanitize_exception:{pdb_id}:{variant_label}",
            [
                f"exception_type: {type(error).__name__}",
                f"exception_repr: {error!r}",
                "traceback:",
                traceback.format_exc().rstrip(),
            ],
        )
        return {
            "status": STATUS_FAILED,
            "variant_label": variant_label,
            "structure_input_path": str(selected_input_path),
            "structure_input_source": input_source_label,
            "internal_capping_status": STATUS_SKIPPED,
            "internal_capping_input_path": "",
            "internal_capping_output_path": "",
            "protonation_input_source": "sanitize_exception",
            "protonation_input_path": str(selected_input_path),
            "protonation_output_path": "",
            "final_protein_input_paths": "",
        }

    sanitized_input_path = Path(str(getattr(sanitize_result, "output_pdb_path", "")))
    if not bool(getattr(sanitize_result, "can_run_gromacs", False)):
        _cfg._log(
            "protonation_route_linear:sanitize_failed",
            [
                f"pdb_id                 : {pdb_id}",
                f"variant_label          : {variant_label}",
                f"input_path             : {selected_input_path}",
                f"sanitized_input_path   : {sanitized_input_path}",
                f"sanitize_log_path      : {getattr(sanitize_result, 'log_path', '')}",
                f"status                 : {getattr(sanitize_result, 'status', '')}",
                f"issues                 : {_summarize_sanitize_issues(sanitize_result)}",
            ],
        )
        return {
            "status": STATUS_FAILED,
            "variant_label": variant_label,
            "structure_input_path": str(selected_input_path),
            "structure_input_source": input_source_label,
            "internal_capping_status": STATUS_SKIPPED,
            "internal_capping_input_path": "",
            "internal_capping_output_path": "",
            "protonation_input_source": "sanitize_failed",
            "protonation_input_path": str(sanitized_input_path),
            "protonation_output_path": "",
            "final_protein_input_paths": "",
        }

    protonation_result = _run_fruton_protonation_for_variant(
        pdb_id=pdb_id,
        pdb_dir=pdb_dir,
        variant_label=variant_label,
        input_pdb=sanitized_input_path,
        input_source=input_source_label,
    )

    if not protonation_result.get("protonation_success", False):
        _cfg._log(
            "protonation_route_linear:protonation_failed",
            [
                f"pdb_id                    : {pdb_id}",
                f"variant_label             : {variant_label}",
                f"input_path                : {selected_input_path}",
                f"sanitized_input_path      : {sanitized_input_path}",
                f"input_source              : {input_source_label}",
                f"protonation_method        : {protonation_result.get('protonation_method', '')}",
                f"protonation_output_path   : {protonation_result.get('protonation_output_path', '')}",
                f"protonation_topology_path : {protonation_result.get('protonation_topology_path', '')}",
                f"protonation_error_message : {protonation_result.get('protonation_error_message', '')}",
            ],
        )
        return {
            "status": STATUS_FAILED,
            "variant_label": variant_label,
            "structure_input_path": str(selected_input_path),
            "structure_input_source": input_source_label,
            "internal_capping_status": STATUS_SKIPPED,
            "internal_capping_input_path": "",
            "internal_capping_output_path": "",
            "protonation_input_source": str(protonation_result.get("protonation_input_source", "")),
            "protonation_input_path": str(protonation_result.get("protonation_input_path", "")),
            "protonation_output_path": str(protonation_result.get("protonation_output_path", "")),
            "final_protein_input_paths": "",
        }

    protonation_output_path = Path(str(protonation_result["protonation_output_path"]))

    return {
        "status": STATUS_SUCCESS,
        "variant_label": variant_label,
        "structure_input_path": str(selected_input_path),
        "structure_input_source": input_source_label,
        "internal_capping_status": STATUS_SKIPPED,
        "internal_capping_input_path": "",
        "internal_capping_output_path": "",
        "protonation_input_source": str(protonation_result["protonation_input_source"]),
        "protonation_input_path": str(protonation_result["protonation_input_path"]),
        "protonation_output_path": str(protonation_output_path),
        "protonation_propka_his_assignments": str(
            protonation_result.get("protonation_propka_his_assignments", "")
        ),
        "final_protein_input_paths": str(protonation_output_path),
    }


def _run_gaps_fragment_route(
    pdb_id: str,
    pdb_dir: Path,
    variant_label: str,
) -> dict[str, str]:
    from stack_protein_preparation.pipeline_variants import _select_default_protein_input_path
    import shutil

    if variant_label != "gaps":
        raise ValueError(
            f"_run_gaps_fragment_route only supports variant_label='gaps', got {variant_label!r}"
        )

    monomer_protein_input_path = _select_default_protein_input_path(
        pdb_id=pdb_id,
        pdb_dir=pdb_dir,
    )
    components_dir = pdb_dir / "components"
    temporary_fragment_dir = components_dir / ".tmp_gaps_fragments"

    if temporary_fragment_dir.exists():
        shutil.rmtree(temporary_fragment_dir)
    temporary_fragment_dir.mkdir(parents=True, exist_ok=True)

    _cfg._log(
        "protonation_route_gaps:start",
        [
            f"pdb_id                     : {pdb_id}",
            f"variant_label              : {variant_label}",
            f"monomer_protein_input_path : {monomer_protein_input_path}",
            f"temporary_fragment_dir     : {temporary_fragment_dir}",
            "internal_capping            : skipped",
            "internal_capping_backend    : none",
            f"protonation_method         : gmx pdb2gmx -ignh",
            f"gromacs_ff                 : {_cfg._force_field}",
            f"gromacs_water              : {_cfg._water_model}",
            "route                      : monomer protein -> temporary fragment split -> per-fragment sanitize -> per-fragment gmx pdb2gmx -> prepared_structure merge",
        ],
    )

    try:
        split_result = split_gap_variant_structure(
            pdb_id=pdb_id,
            protein_dir=pdb_dir,
            variant_label=variant_label,
            input_pdb_path=monomer_protein_input_path,
            output_dir=temporary_fragment_dir,
            verbose=False,
        )
    except Exception as error:
        _cfg._log(
            f"protonation_route_gaps:fragment_split_failed:{pdb_id}",
            [
                f"exception_type: {type(error).__name__}",
                f"exception_repr: {error!r}",
                "traceback:",
                traceback.format_exc().rstrip(),
            ],
        )
        return {
            "status": STATUS_FAILED,
            "variant_label": variant_label,
            "structure_input_path": str(monomer_protein_input_path),
            "structure_input_source": "protein_monomer",
            "internal_capping_status": STATUS_SKIPPED,
            "internal_capping_input_path": "",
            "internal_capping_output_path": "",
            "protonation_input_source": "",
            "protonation_input_path": "",
            "protonation_output_path": "",
            "final_protein_input_paths": "",
        }

    fragment_path_list = _extract_fragment_paths_from_split_result(split_result)
    if not fragment_path_list:
        _cfg._log(
            "protonation_route_gaps:fragment_split_empty",
            [
                f"pdb_id                 : {pdb_id}",
                f"variant_label          : {variant_label}",
                f"temporary_fragment_dir : {temporary_fragment_dir}",
                f"split_result           : {split_result}",
            ],
        )
        return {
            "status": STATUS_FAILED,
            "variant_label": variant_label,
            "structure_input_path": str(monomer_protein_input_path),
            "structure_input_source": "protein_monomer",
            "internal_capping_status": STATUS_SKIPPED,
            "internal_capping_input_path": "",
            "internal_capping_output_path": "",
            "protonation_input_source": "",
            "protonation_input_path": "",
            "protonation_output_path": "",
            "final_protein_input_paths": "",
        }

    protonated_fragment_path_list: list[Path] = []
    protonation_input_path_list: list[Path] = []
    protonation_output_path_list: list[Path] = []

    for fragment_index, fragment_path in enumerate(fragment_path_list, start=1):
        fragment_variant_label = f"{variant_label}_fragment_{fragment_index:02d}"

        try:
            sanitize_result = _run_sanitize_for_protonation_input(
                pdb_id=pdb_id,
                pdb_dir=pdb_dir,
                variant_label=fragment_variant_label,
                input_pdb=fragment_path,
            )
        except Exception as error:
            protonation_input_path_list.append(fragment_path)
            _cfg._log(
                f"protonation_route_gaps:fragment_sanitize_exception:{pdb_id}:{fragment_variant_label}",
                [
                    f"exception_type: {type(error).__name__}",
                    f"exception_repr: {error!r}",
                    "traceback:",
                    traceback.format_exc().rstrip(),
                ],
            )
            continue

        sanitized_fragment_path = Path(str(getattr(sanitize_result, "output_pdb_path", "")))
        protonation_input_path_list.append(sanitized_fragment_path)

        if not bool(getattr(sanitize_result, "can_run_gromacs", False)):
            _cfg._log(
                "protonation_route_gaps:fragment_sanitize_failed",
                [
                    f"pdb_id                  : {pdb_id}",
                    f"variant_label           : {variant_label}",
                    f"fragment_variant_label  : {fragment_variant_label}",
                    f"raw_fragment_path       : {fragment_path}",
                    f"sanitized_fragment_path : {sanitized_fragment_path}",
                    f"sanitize_log_path       : {getattr(sanitize_result, 'log_path', '')}",
                    f"status                  : {getattr(sanitize_result, 'status', '')}",
                    f"issues                  : {_summarize_sanitize_issues(sanitize_result)}",
                ],
            )
            continue

        try:
            protonation_result = _run_fruton_protonation_for_variant(
                pdb_id=pdb_id,
                pdb_dir=pdb_dir,
                variant_label=fragment_variant_label,
                input_pdb=sanitized_fragment_path,
                input_source="protein",
            )
        except Exception as error:
            _cfg._log(
                f"protonation_route_gaps:fragment_protonation_exception:{pdb_id}:{fragment_variant_label}",
                [
                    f"exception_type: {type(error).__name__}",
                    f"exception_repr: {error!r}",
                    "traceback:",
                    traceback.format_exc().rstrip(),
                ],
            )
            continue

        output_path_text = str(protonation_result.get("protonation_output_path", "")).strip()
        output_path = Path(output_path_text) if output_path_text else None
        if (
            protonation_result.get("protonation_success", False)
            and output_path is not None
            and output_path.exists()
            and output_path.stat().st_size > 0
        ):
            protonated_fragment_path_list.append(output_path)
            protonation_output_path_list.append(output_path)
            continue

        _cfg._log(
            "protonation_route_gaps:fragment_protonation_failed",
            [
                f"pdb_id                    : {pdb_id}",
                f"variant_label             : {variant_label}",
                f"fragment_variant_label    : {fragment_variant_label}",
                f"raw_fragment_path         : {fragment_path}",
                f"sanitized_fragment_path   : {sanitized_fragment_path}",
                f"protonation_output_path   : {output_path_text}",
                f"protonation_error_message : {protonation_result.get('protonation_error_message', '')}",
            ],
        )

    if len(protonated_fragment_path_list) != len(fragment_path_list):
        return {
            "status": STATUS_FAILED,
            "variant_label": variant_label,
            "structure_input_path": str(monomer_protein_input_path),
            "structure_input_source": "protein_monomer",
            "internal_capping_status": STATUS_SKIPPED,
            "internal_capping_input_path": "",
            "internal_capping_output_path": "",
            "protonation_input_source": "temporary_fragments",
            "protonation_input_path": _path_list_to_text(protonation_input_path_list),
            "protonation_output_path": _path_list_to_text(protonation_output_path_list),
            "final_protein_input_paths": "",
        }

    _cfg._log(
        "protonation_route_gaps:success",
        [
            f"pdb_id                         : {pdb_id}",
            f"variant_label                  : {variant_label}",
            f"n_raw_fragments                : {len(fragment_path_list)}",
            f"n_protonated_fragments         : {len(protonated_fragment_path_list)}",
            f"raw_fragment_paths             : {_path_list_to_text(fragment_path_list)}",
            f"sanitized_fragment_paths       : {_path_list_to_text(protonation_input_path_list)}",
            f"protonated_fragment_paths      : {_path_list_to_text(protonated_fragment_path_list)}",
            "prepared_structure_input_policy : merge protonated fragment PDB files with TER separation",
        ],
    )

    return {
        "status": STATUS_SUCCESS,
        "variant_label": variant_label,
        "structure_input_path": str(monomer_protein_input_path),
        "structure_input_source": "protein_monomer",
        "internal_capping_status": STATUS_SKIPPED,
        "internal_capping_input_path": "",
        "internal_capping_output_path": "",
        "protonation_input_source": "temporary_fragments",
        "protonation_input_path": _path_list_to_text(protonation_input_path_list),
        "protonation_output_path": _path_list_to_text(protonated_fragment_path_list),
        "final_protein_input_paths": _path_list_to_text(protonated_fragment_path_list),
    }
