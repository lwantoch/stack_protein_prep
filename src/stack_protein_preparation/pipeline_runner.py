"""Stage-execution helpers for the FRUTON pipeline.

Functions that formerly called ``_append_fruton_log`` now call the module-level
``_log`` callable.  Fruton.py calls ``init_runner(_append_fruton_log, ...)``
immediately after defining ``_append_fruton_log`` so the real log function is
wired in before any pipeline stage runs.
"""

from __future__ import annotations

from stack_protein_preparation._pipeline_config import (  # noqa: F401  (re-exported)
    DEFAULT_FORCE_FIELD,
    DEFAULT_PH,
    DEFAULT_WATER_MODEL,
    _force_field,
    _log,
    _ph,
    _water_model,
    init_runner,
)
from stack_protein_preparation._pipeline_helpers import (  # noqa: F401  (re-exported)
    _apply_accepted_variant_audit_decision,
    _build_prepared_cofactor_structure_for_variant,
    _build_prepared_structure_for_variant,
    _build_prepared_wat_structure_for_variant,
    _clear_component_output_files,
    _compute_actual_protein_range,
    _copy_representative_protein_component_to_monomer_path,
    _extract_fragment_paths_from_split_result,
    _extract_monomer_output_path,
    _find_template_pdb_for_filler,
    _find_uniprot_id_for_protein,
    _get_sanitize_result_output_path,
    _path_list_to_text,
    _run_parameter_audit_for_protein,
    _run_parameter_audit_for_variant,
    _store_metal_check_result_in_record,
    _summarize_sanitize_issues,
    _variant_audit_accepts_end_model,
    _variant_metal_check_accepts_end_model,
)
from stack_protein_preparation._pipeline_stages import (  # noqa: F401  (re-exported)
    _append_parameter_audit_summary,
    _append_parameter_audit_variant_summary,
    _check_water_clashes_for_prepared_variant,
    _prepare_crystal_water_component,
    _run_fruton_protonation_for_variant,
    _run_gaps_fragment_route,
    _run_linear_protonation_route_for_variant,
    _run_metals_check_for_input_inventory,
    _run_metals_check_for_prepared_variant,
    _run_metall_params_for_protein,
    _run_nonstd_residue_params_for_protein,
    _run_protonation_route_for_variant,
    _run_sanitize_for_prepared_variant,
    _run_sanitize_for_protonation_input,
    _run_sanitize_for_representative_protein,
    _run_standard_residue_repair_for_protonation_input,
    _write_representative_monomer_for_protein,
)
