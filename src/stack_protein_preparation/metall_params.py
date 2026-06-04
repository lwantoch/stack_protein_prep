"""Metal-center parametrization preparation module.

Purpose
-------
Prepare per-site MCPB scaffold folders for every transition-metal atom found
in a FRUTON protein directory.  For each site the module:

1. Detects transition-metal atoms using the same deterministic logic as
   :mod:`stack_protein_preparation.metalls_check`.
2. Creates a numbered site folder under ``metall_params/`` with four
   sub-directories: ``00_input/``, ``01_hydrogen_cleanup/``, ``02_contacts/``,
   and ``03_mcpb/``.
3. Runs :func:`~stack_protein_preparation.metal_hydrogen_cleanup.remove_metal_coordinating_hydrogens`
   to strip spurious hydrogens from the environment before parametrization.
4. Writes a deterministic ``deterministic_metal_contacts.tsv`` (Python only,
   no Chimera dependency).
5. Writes an MCPB.py input scaffold with ``TODO`` placeholders for the user to
   fill in before running Gaussian.
6. Optionally writes and runs a ChimeraX inspection script (never used for
   decisions; for visual inspection only).
7. Writes a top-level ``manifest.json`` and ``all_sites_summary.tsv``.

Public API
----------
:func:`run_metal_parametrization_for_protein_dir`
"""
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

# Re-export everything from _metall_params_helpers so existing imports keep working
from stack_protein_preparation._metall_params_helpers import (  # noqa: F401 (re-exported)
    DEFAULT_CHIMERAX_EXECUTABLE_CANDIDATES,
    DEFAULT_CHIMERA_EXECUTABLE_CANDIDATES,
    _METAL_RESNAMES,
    _NON_NAA_RESNAMES,
    _GROMACS_TO_AMBER_ATOM,
    _choose_best_protein_input,
    _get_optional_component_paths,
    _has_metal_file,
    _read_pdb_payload_lines,
    _find_naa_resnames_in_contacts,
    _extract_residue_pdb,
    _generate_naa_mol2,
    _is_metal_hetatm_line,
    _write_combined_tmp_param_pdb,
    _copy_metal_only_pdb,
    _write_single_metal_pdb,
    _combined_mcpb_pdb,
    _renumber_pdb_atoms,
    _rename_gromacs_to_amber_atoms,
    _rename_his_by_protonation,
    _generate_water_mol2,
    _resolve_chimera_executable,
    _write_chimerax_cxc_script,
    _run_chimerax_script,
    _write_deterministic_contacts_tsv,
    _process_one_site,
)

# Re-export _metall_params_mcpb symbols (originally re-exported from metall_params)
from stack_protein_preparation._metall_params_mcpb import (  # noqa: F401 (re-exported)
    _UNAMBIGUOUS_SPIN,
    _write_mcpb_in,
    _write_commands_sh,
    _write_submit_gaussian_sh,
    _write_run_after_gaussian_sh,
    _write_metal_mol2,
    _find_metal_serial_in_pdb,
    _write_mcpb_readme,
)

from stack_protein_preparation.metalls_check import (  # noqa: F401 (re-exported)
    DONOR_ELEMENTS,
    METAL_RESIDUE_IDENTITY,
    TRANSITION_ION_GEOMETRY_PREFERENCES,
    TRANSITION_METAL_ELEMENTS,
    classify_coordination_geometry,
    distance_between_atoms,
    find_metal_contacts,
    is_true_metal_atom,
    read_pdb_atom_records,
)


# ---------------------------------------------------------------------------
# Public entry point
# ---------------------------------------------------------------------------

def run_metal_parametrization_for_protein_dir(
    protein_dir: str | Path,
    chimera_executable: str | None = None,
) -> dict[str, Any]:
    """Run the metal parametrization preparation step for one protein directory.

    For every transition-metal atom found in ``<protein_dir>/components/``
    this function creates a numbered site folder under
    ``<protein_dir>/metall_params/`` containing:

    - ``00_input/``         environment PDB and metal-only PDB
    - ``01_hydrogen_cleanup/`` environment with bad H removed
    - ``02_contacts/``      deterministic contacts TSV and optional ChimeraX script
    - ``03_mcpb/``          MCPB.py scaffold (input file, commands, README)
    - ``logs/``             per-site log

    A top-level ``manifest.json`` and ``all_sites_summary.tsv`` are written to
    ``<protein_dir>/metall_params/``.

    Parameters
    ----------
    protein_dir:
        Path to the FRUTON protein directory, e.g. ``data/proteins/2AFX``.
    chimera_executable:
        Optional explicit ChimeraX/Chimera executable for optional visual
        inspection scripts.  If ``None``, the function auto-discovers ChimeraX
        from well-known locations.  Chimera output is never used for decisions.

    Returns
    -------
    dict
        Structured result for pipeline-state integration with keys:
        ``status``, ``pdb_id``, ``metall_params_dir``,
        ``transition_metal_site_count``, ``manifest_path``,
        ``site_results``, ``message``.
    """
    protein_dir = Path(protein_dir).resolve()
    pdb_id = protein_dir.name
    components_dir = protein_dir / "components"
    metall_params_dir = protein_dir / "metall_params"

    result: dict[str, Any] = {
        "status": "failed",
        "pdb_id": pdb_id,
        "metall_params_dir": str(metall_params_dir),
        "transition_metal_site_count": 0,
        "manifest_path": None,
        "site_results": [],
        "message": "",
    }

    if not components_dir.is_dir():
        result["message"] = f"Missing components directory: {components_dir}"
        return result

    if not _has_metal_file(components_dir, pdb_id):
        result["status"] = "skipped"
        result["message"] = "No metal file found. metall_params step skipped."
        return result

    protein_input = _choose_best_protein_input(components_dir, pdb_id)
    optional = _get_optional_component_paths(components_dir, pdb_id)
    water_input = optional["water"]
    ligand_input = optional["ligand"]
    metal_input = optional["metal"]

    if protein_input is None:
        result["message"] = (
            "No suitable protein input found. Expected one of: "
            f"{pdb_id}_proteinH_single.pdb, {pdb_id}_proteinH_best_complete.pdb, "
            f"{pdb_id}_proteinH_gaps.pdb, {pdb_id}_protein_monomer.pdb, "
            f"{pdb_id}_protein.pdb, {pdb_id}_protein_final.pdb, "
            f"{pdb_id}_protein_as_Amber.pdb, {pdb_id}_proteinH.pdb"
        )
        return result

    if metal_input is None:
        result["message"] = "Metal file missing despite positive metal detection."
        return result

    # Read metal records to enumerate individual transition-metal atoms
    try:
        metal_records = read_pdb_atom_records(metal_input, source_role="metal_detection")
    except Exception as exc:
        result["message"] = f"Could not read metal PDB: {exc!r}"
        return result

    transition_metals = [
        a for a in metal_records
        if is_true_metal_atom(a) and a.element in TRANSITION_METAL_ELEMENTS
    ]

    if not transition_metals:
        result["status"] = "skipped"
        result["message"] = (
            "Metal file present but contains no transition-metal atoms. "
            "metall_params step skipped."
        )
        return result

    metall_params_dir.mkdir(parents=True, exist_ok=True)

    # Also read environment for contact analysis
    env_records_for_analysis: list[Any] = []
    try:
        env_records_for_analysis = read_pdb_atom_records(
            protein_input, source_role="env_analysis"
        )
        if water_input:
            env_records_for_analysis += read_pdb_atom_records(
                water_input, source_role="env_analysis"
            )
        if ligand_input:
            env_records_for_analysis += read_pdb_atom_records(
                ligand_input, source_role="env_analysis"
            )
    except Exception:
        pass  # Non-fatal; contacts step will re-read from file

    site_results: list[dict[str, Any]] = []
    contact_cutoff = 3.5  # fixed default; could be parameterised later

    for idx, metal_atom in enumerate(transition_metals, start=1):
        try:
            site_res = _process_one_site(
                site_index=idx,
                metal_atom=metal_atom,
                all_records=env_records_for_analysis + metal_records,
                protein_dir=protein_dir,
                pdb_id=pdb_id,
                protein_input=protein_input,
                water_input=water_input,
                ligand_input=ligand_input,
                metal_input=metal_input,
                metall_params_dir=metall_params_dir,
                chimera_executable=chimera_executable,
                contact_cutoff_angstrom=contact_cutoff,
            )
        except Exception as exc:
            element = getattr(metal_atom, "element", "?")
            chain = getattr(metal_atom, "chain_id", "?") or "?"
            resseq = getattr(metal_atom, "residue_number", 0)
            site_res = {
                "site_index": idx,
                "site_folder": f"site_{idx:03d}_{element}_{chain}_{resseq}",
                "element": element,
                "chain_id": chain,
                "residue_number": resseq,
                "status": "failed",
                "hydrogen_cleanup_status": "failed",
                "removed_h_count": 0,
                "kept_h_count": 0,
                "manual_review_count": 0,
                "coordination_number": 0,
                "mcpb_folder": f"site_{idx:03d}_{element}_{chain}_{resseq}/03_mcpb",
                "message": f"Exception: {exc!r}",
            }
        site_results.append(site_res)

    result["transition_metal_site_count"] = len(transition_metals)
    result["site_results"] = site_results

    # ------------------------------------------------------------------
    # manifest.json
    # ------------------------------------------------------------------
    manifest_sites = [
        {
            "site_index": s["site_index"],
            "site_folder": s["site_folder"],
            "element": s["element"],
            "chain_id": s["chain_id"],
            "residue_number": s["residue_number"],
            "coordination_number": s["coordination_number"],
            "found_geometry": s["found_geometry"],
            "geometry_ok": s["geometry_ok"],
            "hydrogen_cleanup_status": s["hydrogen_cleanup_status"],
            "removed_h_count": s["removed_h_count"],
            "mcpb_scaffold_generated": s["mcpb_scaffold_generated"],
            "mcpb_folder": s["mcpb_folder"],
        }
        for s in site_results
    ]

    all_statuses = {s["status"] for s in site_results}
    if all_statuses == {"success"}:
        overall_status = "success"
    elif "success" in all_statuses:
        overall_status = "partial"
    elif all_statuses == {"skipped"}:
        overall_status = "skipped"
    else:
        overall_status = "failed"

    manifest = {
        "pdb_id": pdb_id,
        "status": overall_status,
        "transition_metal_site_count": len(transition_metals),
        "sites": manifest_sites,
    }
    manifest_path = metall_params_dir / "manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2), encoding="utf-8")

    # ------------------------------------------------------------------
    # all_sites_summary.tsv
    # ------------------------------------------------------------------
    summary_tsv = metall_params_dir / "all_sites_summary.tsv"
    summary_header = "\t".join([
        "site_index", "element", "chain_id", "residue_number",
        "coordination_number", "found_geometry", "geometry_ok",
        "h_removed", "h_kept", "h_manual_review",
        "mcpb_scaffold_generated", "mcpb_folder", "status",
    ])
    summary_rows = []
    for s in site_results:
        summary_rows.append("\t".join([
            str(s["site_index"]),
            s["element"],
            s["chain_id"],
            str(s["residue_number"]),
            str(s["coordination_number"]),
            s.get("found_geometry", ""),
            str(s.get("geometry_ok", False)),
            str(s["removed_h_count"]),
            str(s["kept_h_count"]),
            str(s["manual_review_count"]),
            str(s.get("mcpb_scaffold_generated", False)),
            s["mcpb_folder"],
            s["status"],
        ]))
    summary_tsv.write_text(
        summary_header + "\n" + "\n".join(summary_rows) + "\n",
        encoding="utf-8",
    )

    result["status"] = overall_status
    result["manifest_path"] = str(manifest_path)
    result["message"] = (
        f"{overall_status.capitalize()}: processed {len(transition_metals)} "
        f"transition-metal site(s) for {pdb_id}."
    )
    return result


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description=(
            "Prepare per-site MCPB scaffold folders for metal-containing "
            "protein directories."
        )
    )
    parser.add_argument(
        "protein_dir",
        type=str,
        help="Path to one protein directory, e.g. data/proteins/2AFX",
    )
    parser.add_argument(
        "--chimera",
        dest="chimera_executable",
        type=str,
        default=None,
        help="Optional explicit ChimeraX executable path or name.",
    )

    args = parser.parse_args()
    cli_result = run_metal_parametrization_for_protein_dir(
        protein_dir=args.protein_dir,
        chimera_executable=args.chimera_executable,
    )
    import json as _json
    print(_json.dumps(cli_result, indent=2))
