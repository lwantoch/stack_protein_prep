"""Audit whether protein residues are compatible with the active GROMACS force field.

This module keeps the chemistry-specific logic deliberately small. It uses
BioPython to inspect which residue names occur in a coordinate file, but it does
not try to reimplement the GROMACS residue-template machinery. The authoritative
check is a temporary ``gmx pdb2gmx`` probe run with the same force field, water
model, chain-separation mode, and merge mode that the FRUTON protonation stage
uses.

The audit is meant to sit after component splitting and before reconstruction or
protonation decisions. A residue can be non-standard from the PDB/BioPython point
of view while still being usable by the selected force field, for example a
force-field alias or a known protonation-state name. Conversely, a standard
residue can still fail if the input coordinate file is incomplete, for example a
proline missing the heavy atom ``CD``. The result therefore separates unsupported
chemistry from repairable input-quality problems.

The public entry point is ``audit_protein_residue_parameters``. It returns a
``ParameterAuditResult`` dataclass and can optionally write a JSON manifest and a
plain text probe log. The module performs no mutation of the input PDB and does
not modify FRUTON pipeline state by itself.
"""

from __future__ import annotations

from pathlib import Path

# Re-export everything from _parameter_audit_core so existing imports keep working
from stack_protein_preparation._parameter_audit_core import (  # noqa: F401 (re-exported)
    MODULE_NAME,
    AuditStatus,
    IssueCategory,
    ResidueCategory,
    DEFAULT_GROMACS_FORCE_FIELD,
    DEFAULT_GROMACS_WATER_MODEL,
    DEFAULT_GROMACS_CHAIN_SEPARATION,
    DEFAULT_GROMACS_MERGE_MODE,
    DEFAULT_GROMACS_TIMEOUT_SECONDS,
    STANDARD_AMINO_ACID_RESNAMES,
    FORCE_FIELD_ALIAS_RESNAMES,
    CANONICALIZABLE_RESNAMES,
    COMMON_PROTEIN_LIKE_NONSTANDARD_RESNAMES,
    METAL_RESNAMES,
    WATER_RESNAMES,
    ResidueSite,
    ParameterAuditIssue,
    GromacsProbeResult,
    ParameterAuditResult,
    TleapProbeResult,
    RespParameterProbeResult,
    collect_residue_sites,
    classify_residue_site,
    run_gromacs_parameter_probe,
    parse_gromacs_probe_issues,
    write_parameter_audit_manifest,
    write_parameter_audit_log,
    probe_amber_resp_parameters,
    result_to_pipeline_fields,
    _resolve_pdb_id,
    _find_gromacs_executable,
    _supported_nonstandard_residue_names,
    _unsupported_residue_names,
    _audit_status,
    _normalize_chain_id,
    _compact_message,
    _stderr_tail,
    _deduplicate_issues,
    _yes_no,
    _tleap_has_fatal,
    _tleap_error_summary,
)


def audit_protein_residue_parameters(
    input_pdb_path: str | Path,
    *,
    pdb_id: str | None = None,
    force_field: str = DEFAULT_GROMACS_FORCE_FIELD,
    water_model: str = DEFAULT_GROMACS_WATER_MODEL,
    chain_separation: str = DEFAULT_GROMACS_CHAIN_SEPARATION,
    merge: str = DEFAULT_GROMACS_MERGE_MODE,
    gromacs_executable: str | None = None,
    timeout_seconds: int = DEFAULT_GROMACS_TIMEOUT_SECONDS,
    manifest_path: str | Path | None = None,
    log_path: str | Path | None = None,
    keep_probe_directory: bool = False,
) -> ParameterAuditResult:
    """Audit one protein PDB against the selected GROMACS force field.

    The function first collects residue sites from the input file with
    BioPython. It then performs a temporary ``gmx pdb2gmx`` run using the same
    basic options as the current FRUTON protonation route. The result combines
    the local residue classification with parsed GROMACS failure categories so
    the runner can decide whether to continue, repair the input, load external
    parameters, or branch to a QM/metal parametrization workflow.

    Parameters are intentionally close to the protonation module defaults. This
    keeps the audit predictive for the actual route that will run later. The
    probe directory is deleted by default, because all relevant stdout/stderr and
    command metadata are preserved in the returned result and optional log file.
    """

    resolved_input_pdb_path = Path(input_pdb_path).resolve()
    if not resolved_input_pdb_path.is_file():
        raise FileNotFoundError(f"Input PDB not found: {resolved_input_pdb_path}")

    resolved_pdb_id = _resolve_pdb_id(
        pdb_id=pdb_id,
        input_pdb_path=resolved_input_pdb_path,
    )
    resolved_gromacs_executable = _find_gromacs_executable(gromacs_executable)

    residue_sites = collect_residue_sites(resolved_input_pdb_path)
    nonstandard_sites = tuple(
        site for site in residue_sites if site.category not in {"standard", "forcefield_alias"}
    )
    metal_sites = tuple(site for site in residue_sites if site.category == "metal")

    probe = run_gromacs_parameter_probe(
        input_pdb_path=resolved_input_pdb_path,
        gromacs_executable=resolved_gromacs_executable,
        force_field=force_field,
        water_model=water_model,
        chain_separation=chain_separation,
        merge=merge,
        timeout_seconds=timeout_seconds,
        keep_probe_directory=keep_probe_directory,
    )

    issues = parse_gromacs_probe_issues(probe.stderr, probe.error_message)
    supported_nonstandard_names = _supported_nonstandard_residue_names(
        nonstandard_sites=nonstandard_sites,
        probe_success=probe.success,
    )
    unsupported_names = _unsupported_residue_names(
        nonstandard_sites=nonstandard_sites,
        issues=issues,
        probe_success=probe.success,
    )

    requires_repair = any(issue.category == "missing_heavy_atom" for issue in issues)
    requires_metal_parameters = bool(metal_sites)
    requires_external_parameters = bool(unsupported_names) and not probe.success
    requires_qm_parameters = requires_external_parameters and not requires_metal_parameters

    status = _audit_status(
        probe_success=probe.success,
        requires_repair=requires_repair,
        requires_external_parameters=requires_external_parameters,
        requires_metal_parameters=requires_metal_parameters,
    )

    result = ParameterAuditResult(
        pdb_id=resolved_pdb_id,
        input_pdb_path=resolved_input_pdb_path,
        force_field=force_field,
        water_model=water_model,
        gromacs_executable=resolved_gromacs_executable,
        status=status,
        gromacs_probe_success=probe.success,
        can_use_current_force_field=probe.success,
        requires_repair=requires_repair,
        requires_external_parameters=requires_external_parameters,
        requires_qm_parameters=requires_qm_parameters,
        requires_metal_parameters=requires_metal_parameters,
        residue_sites=residue_sites,
        nonstandard_residue_sites=nonstandard_sites,
        metal_sites=metal_sites,
        supported_nonstandard_residue_names=supported_nonstandard_names,
        unsupported_residue_names=unsupported_names,
        issues=issues,
        probe=probe,
        manifest_path=Path(manifest_path).resolve() if manifest_path is not None else None,
        log_path=Path(log_path).resolve() if log_path is not None else None,
    )

    if result.log_path is not None:
        write_parameter_audit_log(result=result, log_path=result.log_path)

    if result.manifest_path is not None:
        write_parameter_audit_manifest(result=result, manifest_path=result.manifest_path)

    return result
