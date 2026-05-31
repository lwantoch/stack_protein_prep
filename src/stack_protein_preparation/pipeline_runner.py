"""Stage-execution helpers for the FRUTON pipeline.

Functions that formerly called ``_append_fruton_log`` now call the module-level
``_log`` callable.  Fruton.py calls ``init_runner(_append_fruton_log, ...)``
immediately after defining ``_append_fruton_log`` so the real log function is
wired in before any pipeline stage runs.
"""

from __future__ import annotations

import re
from pathlib import Path
from typing import Any, Callable

from stack_protein_preparation.metalls_check import (
    run_metals_check_for_model,
    run_metals_inventory_for_structure,
)
from stack_protein_preparation.monomer import write_single_representative_monomer_unit
from stack_protein_preparation.parameter_audit import audit_protein_residue_parameters
from stack_protein_preparation.pdb_components import split_pdb_components
from stack_protein_preparation.pipeline_paths import (
    _build_monomer_protein_path,
    _build_prepared_variant_sanitize_log_path,
    _build_raw_protein_path,
    _build_representative_sanitize_log_path,
    _build_representative_unit_path,
    _build_sanitize_log_path,
    _build_sanitized_monomer_protein_path,
    _build_sanitized_prepared_variant_path,
    _build_sanitized_protonation_path,
    _build_standard_residue_repair_log_path,
    _build_standard_residue_repair_path,
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
    STATUS_FAILED,
    STATUS_SKIPPED,
    STATUS_SUCCESS,
)
from stack_protein_preparation.pipeline_variants import (
    _choose_representative_successful_variant,
    _update_available_models_field,
    _choose_variant_structure_input_path,
    _select_default_protein_input_path,
)
from stack_protein_preparation.prepared_structure import build_prepared_structure_for_variant
from stack_protein_preparation.protonation import protonate_variant_structure
from stack_protein_preparation.sanitize import sanitize_pdb_for_gromacs
from stack_protein_preparation.standard_residue_repair import (
    repair_standard_residue_heavy_atoms_for_gromacs,
)
from stack_protein_preparation.fragment_split import split_gap_variant_structure

# ---------------------------------------------------------------------------
# Module-level injectable logger
# ---------------------------------------------------------------------------

DEFAULT_FORCE_FIELD = "amber99sb-ildn"
DEFAULT_WATER_MODEL = "tip3p"
DEFAULT_PH = 7.0

_log: Callable[[str, list[str]], None] = lambda _tag, _lines: None
_force_field: str = DEFAULT_FORCE_FIELD
_water_model: str = DEFAULT_WATER_MODEL
_ph: float = DEFAULT_PH


def init_runner(
    log_fn: Callable[[str, list[str]], None],
    force_field: str = DEFAULT_FORCE_FIELD,
    water_model: str = DEFAULT_WATER_MODEL,
    ph: float = DEFAULT_PH,
) -> None:
    global _log, _force_field, _water_model, _ph
    _log = log_fn
    _force_field = force_field
    _water_model = water_model
    _ph = ph


# ---------------------------------------------------------------------------
# Pure helpers (no log needed)
# ---------------------------------------------------------------------------

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
    """Remove stale component files before splitting the representative unit.

    The component splitter writes category files such as protein, ligand, water,
    and metal outputs, but an old output can survive when a new representative
    unit no longer contains that category. That is dangerous here because a
    ligand file from the full multimer can be reintroduced later by prepared
    structure assembly. Removing known component outputs before splitting keeps
    all component files synchronized with the representative unit.
    """

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
    """Copy the split representative protein component to the monomer path.

    Downstream stages already use ``components/<PDBID>_protein_monomer.pdb`` as
    the strict default protein input. After Step 8 is changed to split the full
    representative unit, ``components/<PDBID>_protein.pdb`` already contains only
    the representative protein chains. This helper mirrors that file to the
    historical monomer path so the rest of the runner does not need a broader
    path-interface change.
    """

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
    """Return a usable sanitizer output path or the original fallback path."""

    raw_output_path = str(getattr(sanitize_result, "output_pdb_path", "")).strip()
    if raw_output_path:
        output_path = Path(raw_output_path)
        if output_path.exists() and output_path.stat().st_size > 0:
            return output_path

    return fallback_path


def _summarize_sanitize_issues(sanitize_result: Any) -> str:
    """Return a compact one-line summary of sanitizer issues, with counts."""

    issue_list = list(getattr(sanitize_result, "issues", ()) or ())
    if not issue_list:
        return "<none>"

    # Count occurrences of each severity:code pair, preserving first-seen order.
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
    """Run the residue parameter audit for one representative protein.

    The audit writes a per-protein JSON manifest and a per-protein text log
    under the normal FRUTON output tree. It does not change the input PDB and it
    does not generate Gaussian files. Its role is only to decide whether the
    current GROMACS force field accepts the concrete representative protein,
    whether the input needs repair, or whether unsupported residues should be
    sent to the separate residue-parameter workflow.
    """

    audit_manifest_path = pdb_dir / "parameter_audit.json"
    audit_log_path = pdb_dir.parent / "logs" / pdb_id / "parameter_audit.log"
    audit_log_path.parent.mkdir(parents=True, exist_ok=True)

    return audit_protein_residue_parameters(
        input_pdb_path=input_pdb_path,
        pdb_id=pdb_id,
        force_field=_force_field,
        water_model=_water_model,
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
    """Run parameter audit on the sanitized protein-only copy of a final variant.

    This audit is the final end-model decision gate. Unlike the earlier
    representative-protein audit, it examines the exact produced model variant
    after FRUTON has completed gap handling, filler selection, protonation, and
    prepared-structure assembly.
    """

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
        force_field=_force_field,
        water_model=_water_model,
        manifest_path=audit_manifest_path,
        log_path=audit_log_path,
    )


def _store_metal_check_result_in_record(
    *,
    pipeline_record: dict[str, str],
    metal_check_result: Any,
) -> None:
    """Store the compact metal check summary in the flat pipeline state."""

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
    """Return whether the variant metal check allows final-model acceptance."""

    value = str(variant_result.get("metals.model_ready", "")).strip().lower()
    return value in {"", "yes", "true", "1"}


def _variant_audit_accepts_end_model(audit_result: Any) -> bool:
    """Return whether a sanitized prepared variant can remain an end model."""

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
    """Return 'first-last' residue range from ATOM records in *pdb_path*.

    Reads only standard ATOM records (not HETATM) to determine the min and max
    residue sequence numbers present in the prepared structure.  Returns an
    empty string if the file cannot be read or contains no ATOM records.
    """
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


def _build_prepared_structure_for_variant(
    pdb_id: str,
    pdb_dir: Path,
    variant_label: str,
    final_protein_input_path: Path | None = None,
    final_protein_input_paths: list[Path] | None = None,
) -> Path:
    ligand_input_path = pdb_dir / "components" / f"{pdb_id}_ligand.pdb"
    metals_input_path = pdb_dir / "components" / f"{pdb_id}_metals.pdb"

    protein_component_path = pdb_dir / "components" / f"{pdb_id}_protein.pdb"
    kwargs: dict[str, object] = {
        "pdb_directory": pdb_dir,
        "pdb_id": pdb_id,
        "structure_variant": variant_label,
        "water_input_path": None,
        "ligand_input_path": ligand_input_path if ligand_input_path.exists() else None,
        "metals_input_path": metals_input_path if metals_input_path.exists() else None,
        "backbone_nonstd_input_path": protein_component_path if protein_component_path.exists() else None,
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
    return summary.output_pdb_path


def _build_prepared_wat_structure_for_variant(
    pdb_id: str,
    pdb_dir: Path,
    variant_label: str,
    final_protein_input_path: Path | None = None,
    final_protein_input_paths: list[Path] | None = None,
) -> Path | None:
    """Build the protein + crystal-water companion model (<PDBID>_WAT.pdb).

    Returns the output path when the water component is present and non-empty,
    or None when no crystal waters are available.
    """
    from stack_protein_preparation.prepared_structure import (
        build_prepared_structure,
        get_prepared_structure_output_path,
        sanitize_variant_label,
    )

    water_path = pdb_dir / "components" / f"{pdb_id}_water.pdb"
    if not water_path.exists() or water_path.stat().st_size == 0:
        return None

    ligand_input_path = pdb_dir / "components" / f"{pdb_id}_ligand.pdb"
    metals_input_path = pdb_dir / "components" / f"{pdb_id}_metals.pdb"
    protein_component_path = pdb_dir / "components" / f"{pdb_id}_protein.pdb"

    safe_variant = sanitize_variant_label(variant_label)
    wat_output_path = pdb_dir / "prepared" / safe_variant / f"{pdb_id}_WAT.pdb"

    kwargs: dict[str, object] = {
        "output_pdb_path": wat_output_path,
        "water_input_path": water_path,
        "ligand_input_path": ligand_input_path if ligand_input_path.exists() else None,
        "metals_input_path": metals_input_path if metals_input_path.exists() else None,
        "backbone_nonstd_input_path": protein_component_path if protein_component_path.exists() else None,
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
    """Build the protein + cofactor companion model (<PDBID>_cofactor.pdb).

    Uses the auto-detected cofactor component.  Returns None when no cofactors
    were extracted.  Manual inspection of components/per_hetatm/ is recommended
    to verify the auto-classification.
    """
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
        "backbone_nonstd_input_path": protein_component_path if protein_component_path.exists() else None,
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


def _run_fruton_protonation_for_variant(
    *,
    pdb_id: str,
    pdb_dir: Path,
    variant_label: str,
    input_pdb: Path,
    input_source: str,
) -> dict[str, str | bool | float | int]:
    """Run the FRUTON-default GROMACS protonation route for one variant.

    The runner treats ``gmx pdb2gmx -ignh`` as the protein protonation and
    GROMACS naming step. The output is passed directly into prepared-structure
    assembly because it should remain in the naming convention produced by the
    selected GROMACS force field. The historical pH value is still passed as
    compatibility metadata because the protonation module records it, but it is
    not applied directly by GROMACS.
    """

    return protonate_variant_structure(
        pdb_id=pdb_id,
        protein_dir=pdb_dir,
        variant_label=variant_label,
        input_pdb=input_pdb,
        input_source=input_source,  # type: ignore[arg-type]
        ph=_ph,
        ff=_force_field,
        water_model=_water_model,
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


# ---------------------------------------------------------------------------
# Functions that use _log
# ---------------------------------------------------------------------------

def _write_representative_monomer_for_protein(
    pdb_id: str,
    pdb_dir: Path,
    input_pdb_path: Path,
) -> Path:
    """Write the representative unit from the full processed PDB input.

    The monomer selector must run on ``<PDBID>_delins.pdb`` rather than on
    ``components/<PDBID>_protein.pdb``. The protein component file is already
    stripped of ligands, metals, cofactors, and other HETATM records, so a
    representative monomer written from that file can never retain ligand records
    per selected chain. Running on the full processed input lets ``monomer.py``
    keep non-water HETATM records assigned to the selected representative chains.

    This helper only writes ``components/<PDBID>_representative_unit.pdb``. Step 8
    then splits that representative unit into protein, ligand, water, and metal
    component files, so downstream prepared-structure assembly cannot reintroduce
    ligand records from discarded duplicate chains.
    """

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

    result = write_single_representative_monomer_unit(
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

    _log(
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
    """Write the compact parameter-audit decision into the FRUTON log."""

    _log(
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
    """Run minimal MODELLER missing-heavy-atom repair before sanitizer."""

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
        force_field=_force_field,
        allow_nonstandard_residues=False,
    )
    _log(
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
    """Sanitize one structure before passing it to ``pdb2gmx``.

    FRUTON treats sanitizing as a deterministic cleanup gate immediately before
    protonation. The gate selects alternate locations, normalizes occupancies,
    removes non-polymer heterogens for the protein-only route, and writes a
    clean PDB copy. It does not build atoms, replace PTMs, add loops, add
    hydrogens, or generate non-standard residue parameters; unresolved chemistry
    is left visible for parameter audit and the final GROMACS check.
    """

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
        force_field=_force_field,
        log_path=log_path,
        fail_on_missing_heavy_atoms=False,
        fail_on_nonstandard_residues=False,
    )

    _atom_in = getattr(sanitize_result, 'atom_count_in', '')
    _atom_out = getattr(sanitize_result, 'atom_count_out', '')
    _altloc = getattr(sanitize_result, 'selected_altloc_count', '')
    _occ = getattr(sanitize_result, 'normalized_occupancy_count', '')
    _missing = getattr(sanitize_result, 'missing_heavy_atom_count', '')
    _repaired = getattr(sanitize_result, 'repaired_heavy_atom_count', 0)
    _blocked = getattr(sanitize_result, 'blocked_repair_count', 0)
    _nonstd = ', '.join(getattr(sanitize_result, 'nonstandard_residue_names', ())) or '(none)'
    _log(
        "sanitize:protonation_input",
        [
            f"  {pdb_id}  /  {variant_label}",
            f"  status       :  {getattr(sanitize_result, 'status', '')}  "
            f"(gromacs_ok: {getattr(sanitize_result, 'can_run_gromacs', '')})",
            f"  atoms        :  {_atom_in} → {_atom_out}"
            f"  │  altloc: {_altloc}  occ_norm: {_occ}",
            f"  heavy atoms  :  missing: {_missing}  repaired: {_repaired}  blocked: {_blocked}",
            f"  nonstandard  :  {_nonstd}",
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
    """Sanitize the representative monomer protein before gap/filler/protonation.

    This visible sanitize gate runs after gap detection and filler decisions but
    before parameter audit. It writes a cleaned protein-only copy for audit and
    later GROMACS-facing steps without replacing the original monomer file used
    for gap detection and fragment splitting.
    """

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
        force_field=_force_field,
        log_path=log_path,
        fail_on_missing_heavy_atoms=False,
        fail_on_nonstandard_residues=False,
    )

    _atom_in = getattr(sanitize_result, 'atom_count_in', '')
    _atom_out = getattr(sanitize_result, 'atom_count_out', '')
    _altloc = getattr(sanitize_result, 'selected_altloc_count', '')
    _occ = getattr(sanitize_result, 'normalized_occupancy_count', '')
    _missing = getattr(sanitize_result, 'missing_heavy_atom_count', '')
    _repaired = getattr(sanitize_result, 'repaired_heavy_atom_count', 0)
    _blocked = getattr(sanitize_result, 'blocked_repair_count', 0)
    _nonstd = ', '.join(getattr(sanitize_result, 'nonstandard_residue_names', ())) or '(none)'
    _log(
        "sanitize:representative_protein",
        [
            f"  {pdb_id}  /  representative_protein",
            f"  status       :  {getattr(sanitize_result, 'status', '')}  "
            f"(gromacs_ok: {getattr(sanitize_result, 'can_run_gromacs', '')})",
            f"  atoms        :  {_atom_in} → {_atom_out}"
            f"  │  altloc: {_altloc}  occ_norm: {_occ}",
            f"  heavy atoms  :  missing: {_missing}  repaired: {_repaired}  blocked: {_blocked}",
            f"  nonstandard  :  {_nonstd}",
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
    """Sanitize one produced prepared variant before final parameter audit.

    This is intentionally later than gap detection, filler, and protonation. The
    sanitizer therefore tests the model FRUTON actually produced, not the raw
    representative monomer. It is still sanitize-only: it selects alternate
    locations, normalizes occupancies, drops non-polymer heterogens for the
    protein-only audit route, and writes a cleaned audit copy without repairing
    missing atoms or replacing non-standard residues.
    """

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
        force_field=_force_field,
        log_path=log_path,
        fail_on_missing_heavy_atoms=False,
        fail_on_nonstandard_residues=False,
    )

    _atom_in = getattr(sanitize_result, 'atom_count_in', '')
    _atom_out = getattr(sanitize_result, 'atom_count_out', '')
    _altloc = getattr(sanitize_result, 'selected_altloc_count', '')
    _occ = getattr(sanitize_result, 'normalized_occupancy_count', '')
    _missing = getattr(sanitize_result, 'missing_heavy_atom_count', '')
    _repaired = getattr(sanitize_result, 'repaired_heavy_atom_count', 0)
    _blocked = getattr(sanitize_result, 'blocked_repair_count', 0)
    _nonstd = ', '.join(getattr(sanitize_result, 'nonstandard_residue_names', ())) or '(none)'
    _log(
        "sanitize:prepared_variant",
        [
            f"  {pdb_id}  /  {variant_label}",
            f"  status       :  {getattr(sanitize_result, 'status', '')}  "
            f"(gromacs_ok: {getattr(sanitize_result, 'can_run_gromacs', '')})",
            f"  atoms        :  {_atom_in} → {_atom_out}"
            f"  │  altloc: {_altloc}  occ_norm: {_occ}",
            f"  heavy atoms  :  missing: {_missing}  repaired: {_repaired}  blocked: {_blocked}",
            f"  nonstandard  :  {_nonstd}",
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
    """Write compact final variant-audit decision into the FRUTON log."""

    _accepted_tag = "ACCEPTED" if accepted else "REJECTED"
    _log(
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
    """Add TIP3P hydrogens to the crystal water component and rename to GROMACS SOL.

    Calls ``add_hydrogens_to_crystal_water_pdb`` which renames the water residues
    to SOL (OW + HW1 + HW2) using ``gmx pdb2gmx -ignh`` with the configured
    force field and water model. Skipped silently when the water component is
    absent or empty.
    """
    from stack_protein_preparation.protonation import add_hydrogens_to_crystal_water_pdb

    water_path = pdb_dir / "components" / f"{pdb_id}_water.pdb"

    if not water_path.exists() or water_path.stat().st_size == 0:
        _log(
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
        ff=_force_field,
        water_model=_water_model,
    )

    n_atom_lines = sum(
        1 for line in water_path.read_text(encoding="utf-8").splitlines()
        if line.startswith("ATOM") or line.startswith("HETATM")
    )
    _log(
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
    """Report steric clashes between crystal waters (SOL OW) and protein atoms.

    Reads the prepared structure (which already contains protein + SOL water),
    extracts protein heavy-atom coordinates (ATOM records, non-H element) and
    water OW coordinates (HETATM SOL OW), then for each OW finds the minimum
    distance to any protein heavy atom. Results are written to the pipeline log;
    clashes do not affect variant acceptance.
    """
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
        _log(
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
    _log(
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
    """Record metal identity from the component-split input structure.

    This inventory runs before any prepared variant exists. It is deliberately
    independent from final-model acceptance: it fills the workbook with metal
    identity and standard-ion policy even for structures that later fail
    protonation or parameter audit. Crystal waters and the separated metals
    component are included so coordination analysis has the same context that
    the prepared-variant check will later use.
    """

    result = run_metals_inventory_for_structure(
        pdb_id=pdb_id,
        protein_dir=pdb_dir,
        model_pdb_path=monomer_path,
        metals_pdb_path=pdb_dir / "components" / f"{pdb_id}_metals.pdb",
        water_pdb_path=pdb_dir / "components" / f"{pdb_id}_water.pdb",
        output_dir=pdb_dir / "components" / "metals_check" / "input_inventory",
        use_chimera=False,
    )
    _log(
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
    """Run metal classification on one produced prepared variant."""

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
    _log(
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
    """Run per-site MCPB scaffold preparation for one protein directory.

    Calls :func:`~stack_protein_preparation.metall_params.run_metal_parametrization_for_protein_dir`
    which internally:
    - Enumerates transition-metal atoms in ``components/{pdb_id}_metal.pdb``
    - Runs deterministic H cleanup on each site (geometry-aware)
    - Generates the MCPB scaffold only when site geometry matches the standard
      coordination geometry for that element
    - Writes ``manifest.json`` and ``all_sites_summary.tsv`` under ``metall_params/``
    """
    from stack_protein_preparation.metall_params import run_metal_parametrization_for_protein_dir

    result = run_metal_parametrization_for_protein_dir(protein_dir=pdb_dir)
    _log(
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
    """Run non-standard residue RESP scaffold generation for one protein directory."""
    from stack_protein_preparation.nonstd_residue_params import run_nonstd_residue_params

    protein_pdb = pdb_dir / "components" / f"{pdb_id}_protein.pdb"
    result = run_nonstd_residue_params(
        protein_dir=pdb_dir,
        pdb_id=pdb_id,
        protein_pdb=protein_pdb,
    )
    _log(
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


def _apply_accepted_variant_audit_decision(
    *,
    pipeline_record: dict[str, str],
    accepted_variant_result_list: list[dict[str, str]],
) -> None:
    """Update final prepared fields after variant-level parameter audit.

    Step 11 may produce several tentative prepared variants. Step 13 is the
    final end-model decision: only variants whose sanitized prepared model passes
    parameter audit remain visible as final prepared outputs in the pipeline
    state and workbook.
    """
    from stack_protein_preparation.pipeline_records import _clear_prepared_variant_output_paths
    from stack_protein_preparation.pipeline_state import (
        PREPARED_ALPHAFOLD_OUTPUT_PATH_COLUMN_NAME,
        PREPARED_GAPS_OUTPUT_PATH_COLUMN_NAME,
        PREPARED_MODELLER_OUTPUT_PATH_COLUMN_NAME,
        PREPARED_STRUCTURE_ACTUAL_RANGE_COLUMN_NAME,
        PREPARED_STRUCTURE_OUTPUT_PATH_COLUMN_NAME,
        PREPARED_STRUCTURE_PROTEIN_INPUT_PATH_COLUMN_NAME,
        PREPARED_STRUCTURE_STATUS_COLUMN_NAME,
        PREPARED_STRUCTURE_VARIANT_COLUMN_NAME,
    )

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


def _run_linear_protonation_route_for_variant(
    pdb_id: str,
    pdb_dir: Path,
    variant_label: str,
    model_source: str,
    model_path: Path | None,
) -> dict[str, str]:
    _log(
        "protonation_route_linear:start",
        [
            f"pdb_id             : {pdb_id}",
            f"variant_label      : {variant_label}",
            f"model_source       : {model_source}",
            f"model_path         : {model_path}",
            "internal_capping    : skipped",
            f"protonation_method : gmx pdb2gmx -ignh",
            f"protonation_ph     : {_ph} compatibility metadata",
            f"gromacs_ff         : {_force_field}",
            f"gromacs_water      : {_water_model}",
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
        import traceback as _tb
        _log(
            f"protonation_route_linear:sanitize_exception:{pdb_id}:{variant_label}",
            [
                f"exception_type: {type(error).__name__}",
                f"exception_repr: {error!r}",
                "traceback:",
                _tb.format_exc().rstrip(),
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
        _log(
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
        _log(
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
            "protonation_input_source": str(
                protonation_result.get("protonation_input_source", "")
            ),
            "protonation_input_path": str(
                protonation_result.get("protonation_input_path", "")
            ),
            "protonation_output_path": str(
                protonation_result.get("protonation_output_path", "")
            ),
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
    """Run the retained gaps variant as temporary physical fragments.

    The conservative gaps variant should preserve the observed PDB fragments
    rather than reconstructing missing residues. For this route, FRUTON splits
    the representative monomer protein into connected fragments, lets GROMACS
    treat every fragment as its own molecule with proper termini, and then asks
    prepared_structure to merge the protonated fragment coordinate outputs. The
    raw split fragments are temporary implementation artefacts; the protonated
    fragment outputs remain visible until prepared assembly so failures stay
    inspectable.
    """
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
        import shutil

        shutil.rmtree(temporary_fragment_dir)
    temporary_fragment_dir.mkdir(parents=True, exist_ok=True)

    _log(
        "protonation_route_gaps:start",
        [
            f"pdb_id                     : {pdb_id}",
            f"variant_label              : {variant_label}",
            f"monomer_protein_input_path : {monomer_protein_input_path}",
            f"temporary_fragment_dir     : {temporary_fragment_dir}",
            "internal_capping            : skipped",
            "internal_capping_backend    : none",
            f"protonation_method         : gmx pdb2gmx -ignh",
            f"gromacs_ff                 : {_force_field}",
            f"gromacs_water              : {_water_model}",
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
        import traceback
        _log(
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
        _log(
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
            import traceback
            protonation_input_path_list.append(fragment_path)
            _log(
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
            _log(
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
            import traceback
            _log(
                f"protonation_route_gaps:fragment_protonation_exception:{pdb_id}:{fragment_variant_label}",
                [
                    f"exception_type: {type(error).__name__}",
                    f"exception_repr: {error!r}",
                    "traceback:",
                    traceback.format_exc().rstrip(),
                ],
            )
            continue

        output_path_text = str(
            protonation_result.get("protonation_output_path", "")
        ).strip()
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

        _log(
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

    _log(
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
