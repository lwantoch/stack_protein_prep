"""PROPKA-guided, GROMACS-based protonation for FRUTON protein structures.

Protonation proceeds in two stages:

1. **PROPKA pKa prediction** — PROPKA 3 is run on the input structure at the
   requested pH. Histidine residues whose predicted pKa exceeds the pH are
   renamed HIP (doubly protonated) in a temporary copy of the input file.
   All other HIS residues are left as HIE, which is the AMBER99SB-ILDN
   default for neutral histidine. If PROPKA is unavailable or the run fails
   the stage is skipped silently and the original input is passed unchanged.

2. **GROMACS hydrogen placement** — ``gmx pdb2gmx -ignh`` rebuilds all
   hydrogens from the force-field residue templates, using the PROPKA-renamed
   residue names to select the correct HIS variant templates.

The ``ph`` parameter controls the PROPKA decision boundary and is now
functionally meaningful; it is not passed to pdb2gmx (GROMACS has no
``--with-ph`` equivalent).
"""

from __future__ import annotations

import json
import shlex
import tempfile
from pathlib import Path

from stack_protein_preparation._protonation_core import (  # noqa: F401  (re-exported)
    DEFAULT_FILLER_FINAL_FILENAME,
    DEFAULT_GROMACS_CHAIN_SEPARATION,
    DEFAULT_GROMACS_FORCE_FIELD,
    DEFAULT_GROMACS_MERGE_MODE,
    DEFAULT_GROMACS_TIMEOUT_SECONDS,
    DEFAULT_GROMACS_WATER_MODEL,
    GromacsMergeMode,
    GromacsChainSeparation,
    GromacsWaterModel,
    MODULE_NAME,
    ProtonationInputSource,
    _append_log_block,
    _append_text_block,
    _apply_protonation_renames,
    _basic_structure_diagnostics,
    _fill_pdb_element_column,
    _find_gmx_executable,
    _get_module_log_path,
    _get_propka_his_assignments,
    _infer_pdb_id_from_path,
    _infer_protein_dir_from_path,
    _log_exception,
    _preview_text,
    _rebuild_missing_sidechain_atoms,
    _reinject_phospho_residues,
    _remove_backbone_incomplete_residues,
    _remove_standalone_residue_chains,
    _rename_water_pdb_to_sol,
    _run_propka,
    _strip_phospho_residues,
    _screen_header,
    _screen_item,
    _screen_result,
    _timestamp,
    _write_text_log,
    add_hydrogens_to_crystal_water_pdb,
    count_atoms_in_pdb,
    count_atoms_in_structure_file,
    predict_protonation_states,
    run_gmx_pdb2gmx_protonation,
    select_protonation_input,
)
from stack_protein_preparation._protonation_paths import (  # noqa: F401  (re-exported)
    build_gromacs_position_restraints_output_path,
    build_gromacs_topology_output_path,
    build_protonation_stderr_log_path,
    build_protonation_stdout_log_path,
    build_standard_protonation_output_path,
    build_variant_protonation_output_path,
    find_default_filler_final_model_path,
    sanitize_variant_label,
)


# ============================================================================
# main protonation entry points
# ============================================================================


def protonate_selected_structure(
    *,
    input_pdb: str | Path,
    output_pdb: str | Path,
    input_source: ProtonationInputSource,
    pdb_id: str | None = None,
    protein_dir: str | Path | None = None,
    ph: float = 7.4,
    ff: str = DEFAULT_GROMACS_FORCE_FIELD,
    water_model: GromacsWaterModel = DEFAULT_GROMACS_WATER_MODEL,
    chain_separation: GromacsChainSeparation = DEFAULT_GROMACS_CHAIN_SEPARATION,
    merge: GromacsMergeMode = DEFAULT_GROMACS_MERGE_MODE,
    ignore_input_hydrogens: bool = True,
    variant_label: str | None = None,
    timeout_seconds: int = DEFAULT_GROMACS_TIMEOUT_SECONDS,
) -> dict[str, str | bool | float | int]:
    """Protonate one explicitly selected structure with GROMACS."""

    input_pdb = Path(input_pdb)
    output_pdb = Path(output_pdb)
    topology_output_path = build_gromacs_topology_output_path(output_pdb)
    position_restraints_output_path = build_gromacs_position_restraints_output_path(output_pdb)

    if not input_pdb.is_file():
        raise FileNotFoundError(f"Structure file not found: {input_pdb}")

    pdb_id_guess = (pdb_id or _infer_pdb_id_from_path(input_pdb)).upper()
    protein_dir_guess = (
        Path(protein_dir)
        if protein_dir is not None
        else _infer_protein_dir_from_path(input_pdb)
    )

    input_atom_count = count_atoms_in_structure_file(input_pdb)
    diagnostics = _basic_structure_diagnostics(input_pdb)

    stdout_log_path = build_protonation_stdout_log_path(
        pdb_id=pdb_id_guess,
        protein_dir=protein_dir_guess,
        variant_label=variant_label,
    )
    stderr_log_path = build_protonation_stderr_log_path(
        pdb_id=pdb_id_guess,
        protein_dir=protein_dir_guess,
        variant_label=variant_label,
    )
    module_log_path = _get_module_log_path(
        pdb_id=pdb_id_guess,
        protein_dir=protein_dir_guess,
    )

    try:
        executable = _find_gmx_executable()
    except Exception as exc:
        _log_exception(module_log_path, "protonate_selected_structure:executable_error", exc)
        raise

    command_preview = [
        executable, "pdb2gmx",
        "-f", str(input_pdb.resolve()),
        "-o", output_pdb.name,
        "-p", topology_output_path.name,
        "-i", position_restraints_output_path.name,
        "-ff", ff,
        "-water", water_model,
        "-chainsep", chain_separation,
        "-merge", merge,
    ]
    if ignore_input_hydrogens:
        command_preview.append("-ignh")

    _screen_header(pdb_id_guess, "starting PROPKA + GROMACS protonation")
    _screen_item("method", "PROPKA pKa prediction → gmx pdb2gmx -ignh")
    _screen_item("input_source", str(input_source))
    _screen_item("input", str(input_pdb))
    _screen_item("output", str(output_pdb))
    _screen_item("topology", str(topology_output_path))
    _screen_item("posre", str(position_restraints_output_path))
    if variant_label is not None:
        _screen_item("variant", variant_label)
    _screen_item("ff", ff)
    _screen_item("water_model", water_model)
    _screen_item("ph", str(ph))
    _screen_item("stdout_log", str(stdout_log_path))
    _screen_item("stderr_log", str(stderr_log_path))
    _screen_item("module_log", str(module_log_path))

    _append_log_block(
        module_log_path,
        "protonate_selected_structure:start",
        [
            "protonation_method          : PROPKA pKa prediction + gmx pdb2gmx -ignh",
            f"input_source                : {input_source}",
            f"input_pdb                   : {input_pdb}",
            f"output_pdb                  : {output_pdb}",
            f"topology_output_path        : {topology_output_path}",
            f"position_restraints_path    : {position_restraints_output_path}",
            f"variant_label               : {variant_label!r}",
            f"ph                          : {ph}",
            f"ff                          : {ff}",
            f"water_model                 : {water_model}",
            f"chain_separation            : {chain_separation}",
            f"merge                       : {merge}",
            f"ignore_input_hydrogens      : {ignore_input_hydrogens}",
            f"timeout_seconds             : {timeout_seconds}",
            f"input_atom_count            : {input_atom_count}",
            f"input_hetatm_count          : {diagnostics['hetatm_count']}",
            f"input_ter_count             : {diagnostics['ter_count']}",
            f"input_end_count             : {diagnostics['end_count']}",
            f"input_has_any_atoms         : {diagnostics['has_any_atoms']}",
            f"protonation_stdout_log_path : {stdout_log_path}",
            f"protonation_stderr_log_path : {stderr_log_path}",
            f"command                     : {shlex.join(command_preview)}",
            f"cwd                         : {output_pdb.parent}",
        ],
    )

    propka_mol = _run_propka(input_pdb, output_dir=protein_dir_guess / "protonation")
    propka_ran = propka_mol is not None
    his_assignments = _get_propka_his_assignments(propka_mol, ph) if propka_ran else []
    propka_renames: dict[tuple[str, int, str], str] = {
        (str(d["chain"]), int(d["res_num"]), str(d["icode"])): str(d["assigned"])
        for d in his_assignments
        if d["assigned"] == "HIP"
    }
    propka_him_count = len(propka_renames)

    # Override propka for any His whose ring nitrogen coordinates a metal
    # (PDB REMARK 620).  propka's pKa-based decision does not know about metal
    # coordination, and frequently picks HIP for a His that is actually a
    # neutral HID (NE2 → metal lone-pair donor) or HIE (ND1 → metal).  Letting
    # the HIP label through means pdb2gmx hands MCPB.py a doubly-protonated
    # ring on a metal donor, the small-model extractor strips both ring Hs
    # without updating the charge, and the resulting QM region has an odd
    # electron count that forces an unphysical doublet Gaussian run.
    from stack_protein_preparation._protonation_core import (
        parse_metal_coordinating_his_overrides,
    )
    # REMARK 620 lives in the original PDB file (component_split strips it).
    # Try the canonical "<protein_dir>/<PDB_ID>.pdb" first; fall back to the
    # protonation input (mostly for tests that pass an already-prepared PDB
    # with REMARK records preserved).
    _remark_source = protein_dir_guess / f"{pdb_id_guess}.pdb"
    if not _remark_source.is_file():
        _remark_source = input_pdb
    metal_his_overrides = parse_metal_coordinating_his_overrides(_remark_source)
    if metal_his_overrides:
        propka_renames.update(metal_his_overrides)
    metal_his_override_count = len(metal_his_overrides)

    # NEW: user-supplied active-site overrides from
    # <pdb_dir>/active_site_overrides.json.  Highest priority: applied AFTER
    # PROPKA + metal-based overrides so a literature-driven catalytic HIP or
    # CYM assignment always wins.
    from stack_protein_preparation._active_site_overrides import (
        load_active_site_overrides,
        load_active_site_overrides_with_reasons,
    )
    user_overrides = load_active_site_overrides(protein_dir_guess)
    n_user_overrides = 0
    if user_overrides:
        propka_renames.update(user_overrides)
        n_user_overrides = len(user_overrides)
        _screen_item(
            "user_active_site_overrides_applied", str(n_user_overrides)
        )
        for entry in load_active_site_overrides_with_reasons(protein_dir_guess):
            _screen_item(
                f"  {entry.get('chain')}/{entry.get('res_num')} -> {entry.get('to')}",
                str(entry.get("reason", "")),
            )
    propka_applied = bool(propka_renames)

    _screen_item("propka_ran", str(propka_ran))
    _screen_item("propka_his_total", str(len(his_assignments)))
    _screen_item("propka_hip_count", str(propka_him_count))
    _screen_item("metal_his_override_count", str(metal_his_override_count))
    _append_log_block(
        module_log_path,
        "protonate_selected_structure:propka",
        [
            f"propka_ran               : {propka_ran}",
            f"propka_ph                : {ph}",
            f"propka_his_total         : {len(his_assignments)}",
            f"propka_hip_count         : {propka_him_count}",
            f"metal_his_override_count : {metal_his_override_count}",
            f"metal_his_overrides      : {metal_his_overrides}",
            f"propka_applied           : {propka_applied}",
            f"propka_assignments       : {his_assignments}",
        ],
    )

    if propka_applied:
        tmp_propka_dir = tempfile.mkdtemp()
        propka_input = Path(tmp_propka_dir) / "propka_renamed.pdb"
        _apply_protonation_renames(input_pdb, propka_input, propka_renames)
        pdb2gmx_input: Path = propka_input
    else:
        tmp_propka_dir = None
        pdb2gmx_input = input_pdb

    _cleaned_input = pdb2gmx_input.parent / (pdb2gmx_input.stem + "_backbone_cleaned.pdb")
    _removed_residues = _remove_backbone_incomplete_residues(
        pdb2gmx_input, _cleaned_input, log_path=module_log_path
    )
    if _removed_residues:
        _append_log_block(
            module_log_path,
            "protonate_selected_structure:backbone_cleanup",
            [f"Removed {len(_removed_residues)} residue(s) with incomplete backbone before pdb2gmx:"]
            + _removed_residues,
        )
        pdb2gmx_input = _cleaned_input

    # Drop chains containing only a single standard residue.  pdb2gmx cannot
    # cap a residue as both N- and C-terminus simultaneously, so 1-residue
    # chains (substrate mimetics, crystal peptide fragments) abort the run.
    _no_singletons = pdb2gmx_input.parent / (pdb2gmx_input.stem + "_no_singletons.pdb")
    _dropped_chains = _remove_standalone_residue_chains(
        pdb2gmx_input, _no_singletons, log_path=module_log_path
    )
    if _dropped_chains:
        _append_log_block(
            module_log_path,
            "protonate_selected_structure:standalone_chain_drop",
            [f"Dropped {len(_dropped_chains)} single-residue chain(s) before pdb2gmx:"]
            + _dropped_chains,
        )
        # Guard: if singleton removal left the file with zero ATOM records
        # (input was a lone 1-residue chain), the "cleaned" file is useless
        # for pdb2gmx.  Fall back to the input path so downstream handlers
        # see the original and can emit a proper error / skip.
        try:
            _no_singletons_atoms = sum(
                1 for _ln in _no_singletons.read_text(
                    encoding="utf-8", errors="replace"
                ).splitlines() if _ln.startswith("ATOM  ")
            )
        except Exception:  # noqa: BLE001
            _no_singletons_atoms = -1
        if _no_singletons_atoms > 0:
            pdb2gmx_input = _no_singletons

    # Rebuild missing sidechain heavy atoms via PDBFixer.  Crystal structures
    # routinely omit atoms with poor density (LYS-CG, GLU-CG, SER-OG, ARG-CG)
    # which causes pdb2gmx to abort with "atom X not found in input file".
    _sidechain_rebuilt = pdb2gmx_input.parent / (pdb2gmx_input.stem + "_sidechain_rebuilt.pdb")
    _rebuilt_atoms = _rebuild_missing_sidechain_atoms(
        pdb2gmx_input, _sidechain_rebuilt, log_path=module_log_path
    )
    if _rebuilt_atoms:
        _append_log_block(
            module_log_path,
            "protonate_selected_structure:sidechain_rebuild",
            [f"Rebuilt sidechain atoms on {len(_rebuilt_atoms)} residue(s) before pdb2gmx:"]
            + _rebuilt_atoms,
        )
        pdb2gmx_input = _sidechain_rebuilt
    elif _sidechain_rebuilt.is_file():
        # PDBFixer ran but found nothing to add — still adopt its normalised
        # PDB (element column, atom ordering) as it is what pdb2gmx will see.
        pdb2gmx_input = _sidechain_rebuilt

    # Strip phospho residues (PTR/TPO/SEP/...) — amber99sb-ildn has no entry
    # for them and pdb2gmx bails out with "chain does not have consistent
    # type".  We reinject them into the output so downstream tleap
    # (leaprc.phosaa14SB) parametrises them.
    _phospho_stripped = pdb2gmx_input.parent / (pdb2gmx_input.stem + "_phospho_stripped.pdb")
    _phospho_extract = pdb2gmx_input.parent / (pdb2gmx_input.stem + "_phospho_extract.pdb")
    _phospho_keys = _strip_phospho_residues(
        pdb2gmx_input,
        _phospho_stripped,
        _phospho_extract,
        log_path=module_log_path,
    )
    if _phospho_keys:
        _append_log_block(
            module_log_path,
            "protonate_selected_structure:phospho_extract",
            [
                f"Extracted {len(_phospho_keys)} phospho residue(s) before pdb2gmx:",
                *[
                    f"  {resname} {resnum.strip()}{icode.strip()} chain {chain}"
                    for chain, resnum, icode, resname in _phospho_keys
                ],
            ],
        )
        pdb2gmx_input = _phospho_stripped

    try:
        result = run_gmx_pdb2gmx_protonation(
            input_pdb=pdb2gmx_input,
            output_pdb=output_pdb,
            topology_output_path=topology_output_path,
            position_restraints_output_path=position_restraints_output_path,
            ff=ff,
            water_model=water_model,
            chain_separation=chain_separation,
            merge=merge,
            ignore_input_hydrogens=ignore_input_hydrogens,
            timeout_seconds=timeout_seconds,
        )
    except Exception as exc:
        _log_exception(module_log_path, "protonate_selected_structure:exception", exc)
        raise
    finally:
        if tmp_propka_dir is not None:
            import shutil as _shutil
            _shutil.rmtree(tmp_propka_dir, ignore_errors=True)

    _write_text_log(stdout_log_path, result.stdout)
    _write_text_log(stderr_log_path, result.stderr)

    output_exists = output_pdb.is_file()
    output_nonempty = output_exists and output_pdb.stat().st_size > 0
    topology_exists = topology_output_path.is_file()
    topology_nonempty = topology_exists and topology_output_path.stat().st_size > 0
    position_restraints_exists = position_restraints_output_path.is_file()
    position_restraints_nonempty = (
        position_restraints_exists and position_restraints_output_path.stat().st_size > 0
    )

    if output_nonempty:
        _fill_pdb_element_column(output_pdb)
        if _phospho_keys and _phospho_extract.is_file():
            _reinjected = _reinject_phospho_residues(output_pdb, _phospho_extract)
            if _reinjected:
                _append_log_block(
                    module_log_path,
                    "protonate_selected_structure:phospho_reinject",
                    [f"Reinjected {_reinjected} phospho-residue atom record(s) into protonated output"],
                )
        output_atom_count = count_atoms_in_structure_file(output_pdb)
        atom_count_increased = output_atom_count > input_atom_count
    else:
        output_atom_count = 0
        atom_count_increased = False

    protonation_success = result.returncode == 0 and output_nonempty

    _PDBFIXER_TRIGGERS = (
        "incomplete ring",
        "missing atom",
        "not found in residue",
        "not found in the input file",
    )
    stderr_lower = result.stderr.lower()
    if not protonation_success and any(kw in stderr_lower for kw in _PDBFIXER_TRIGGERS):
        _append_log_block(
            module_log_path,
            "protonate_selected_structure:modeller_retry_start",
            [
                "pdb2gmx failed with structural defect; attempting MODELLER repair",
                f"trigger keywords matched: {[kw for kw in _PDBFIXER_TRIGGERS if kw in stderr_lower]}",
                f"input: {pdb2gmx_input}",
            ],
        )
        try:
            import tempfile as _tempfile
            from modeller import Environ as _MEnviron
            from modeller.scripts import complete_pdb as _mcomplete_pdb

            _fixed_dir = _tempfile.mkdtemp(prefix="modeller_retry_")
            _fixed_path = Path(_fixed_dir) / "modeller_fixed.pdb"

            _menv = _MEnviron()
            _menv.libs.topology.read(file="$(LIB)/top_heav.lib")
            _menv.libs.parameters.read(file="$(LIB)/par.lib")
            _mdl = _mcomplete_pdb(_menv, str(pdb2gmx_input))
            _mdl.write(file=str(_fixed_path))

            # Preserve HETATMs from input (MODELLER only handles standard AAs)
            _het_lines = [
                _ln for _ln in Path(pdb2gmx_input).read_text(
                    encoding="utf-8", errors="replace"
                ).splitlines()
                if _ln.startswith("HETATM") or _ln.startswith("CONECT")
            ]
            if _het_lines:
                _existing = _fixed_path.read_text(encoding="utf-8")
                if "END\n" in _existing:
                    _existing = _existing.replace(
                        "END\n", "\n".join(_het_lines) + "\nEND\n", 1
                    )
                else:
                    _existing = _existing.rstrip() + "\n" + "\n".join(_het_lines) + "\nEND\n"
                _fixed_path.write_text(_existing, encoding="utf-8")

            _append_log_block(
                module_log_path,
                "protonate_selected_structure:modeller_repair_done",
                [f"fixed pdb written to: {_fixed_path}"],
            )

            _retry_result = run_gmx_pdb2gmx_protonation(
                input_pdb=_fixed_path,
                output_pdb=output_pdb,
                topology_output_path=topology_output_path,
                position_restraints_output_path=position_restraints_output_path,
                ff=ff,
                water_model=water_model,
                chain_separation=chain_separation,
                merge=merge,
                ignore_input_hydrogens=ignore_input_hydrogens,
                timeout_seconds=timeout_seconds,
            )
            result = _retry_result
            _write_text_log(stdout_log_path, result.stdout)
            _write_text_log(stderr_log_path, result.stderr)

            output_exists = output_pdb.is_file()
            output_nonempty = output_exists and output_pdb.stat().st_size > 0
            topology_exists = topology_output_path.is_file()
            topology_nonempty = topology_exists and topology_output_path.stat().st_size > 0
            position_restraints_exists = position_restraints_output_path.is_file()
            position_restraints_nonempty = (
                position_restraints_exists and position_restraints_output_path.stat().st_size > 0
            )
            if output_nonempty:
                _fill_pdb_element_column(output_pdb)
                output_atom_count = count_atoms_in_structure_file(output_pdb)
                atom_count_increased = output_atom_count > input_atom_count
            else:
                output_atom_count = 0
                atom_count_increased = False
            protonation_success = result.returncode == 0 and output_nonempty

            _append_log_block(
                module_log_path,
                "protonate_selected_structure:modeller_retry_result",
                [
                    f"returncode       : {result.returncode}",
                    f"protonation_success: {protonation_success}",
                ],
            )
        except Exception as _modeller_exc:
            _append_log_block(
                module_log_path,
                "protonate_selected_structure:modeller_retry_exception",
                [str(_modeller_exc)],
            )

    _append_log_block(
        module_log_path,
        "protonate_selected_structure:subprocess_result",
        [
            f"returncode                  : {result.returncode}",
            f"stdout_length               : {len(result.stdout)}",
            f"stderr_length               : {len(result.stderr)}",
            f"stdout_log_path             : {stdout_log_path}",
            f"stderr_log_path             : {stderr_log_path}",
            f"output_exists               : {output_exists}",
            f"output_nonempty             : {output_nonempty}",
            f"topology_exists             : {topology_exists}",
            f"topology_nonempty           : {topology_nonempty}",
            f"position_restraints_exists  : {position_restraints_exists}",
            f"position_restraints_nonempty: {position_restraints_nonempty}",
        ],
    )
    _append_text_block(module_log_path, "protonate_selected_structure:stdout", result.stdout)
    _append_text_block(module_log_path, "protonate_selected_structure:stderr", result.stderr)

    # Restore residue numbering against the ORIGINAL caller-supplied input.
    # The in-function restore inside ``run_gmx_pdb2gmx_protonation`` uses
    # whichever file was last passed to pdb2gmx as template — and on the
    # pdbfixer-retry path that is the pdbfixer-fixed coordinates, which
    # pdbfixer itself can renumber.  Re-applying the restore here against the
    # caller's original ``input_pdb`` makes the final output's numbering
    # consistent with what the user / upstream stages saw, regardless of
    # which intermediate transformations were used.
    if protonation_success and output_pdb.is_file():
        try:
            from stack_protein_preparation._protonation_core import (
                restore_residue_numbering_from_template,
            )
            restore_residue_numbering_from_template(
                template_pdb=input_pdb,
                target_pdb=output_pdb,
            )
        except Exception as _restore_exc:
            _append_log_block(
                module_log_path,
                "protonate_selected_structure:restore_residue_numbering_failed",
                [str(_restore_exc)],
            )

    error_message = ""
    if not protonation_success:
        stderr_preview = _preview_text(result.stderr)
        stdout_preview = _preview_text(result.stdout)
        error_parts = [
            f"gmx pdb2gmx failed for input {input_pdb}",
            f"returncode={result.returncode}",
            f"output_exists={output_exists}",
            f"output_nonempty={output_nonempty}",
            f"topology_exists={topology_exists}",
            f"topology_nonempty={topology_nonempty}",
            f"position_restraints_exists={position_restraints_exists}",
            f"position_restraints_nonempty={position_restraints_nonempty}",
            f"input_atom_count={input_atom_count}",
            f"input_hetatm_count={diagnostics['hetatm_count']}",
            f"input_TER_count={diagnostics['ter_count']}",
            f"stdout_log={stdout_log_path}",
            f"stderr_log={stderr_log_path}",
            f"module_log={module_log_path}",
        ]
        if stderr_preview:
            error_parts.append(f"stderr_preview={stderr_preview}")
        elif stdout_preview:
            error_parts.append(f"stdout_preview={stdout_preview}")
        else:
            error_parts.append("No stdout/stderr produced by gmx pdb2gmx.")
        error_message = " | ".join(error_parts)

    _append_log_block(
        module_log_path,
        "protonate_selected_structure:final_status",
        [
            f"protonation_success         : {protonation_success}",
            f"output_exists               : {output_exists}",
            f"output_nonempty             : {output_nonempty}",
            f"output_atom_count           : {output_atom_count}",
            f"atom_count_increased        : {atom_count_increased}",
            f"topology_exists             : {topology_exists}",
            f"topology_nonempty           : {topology_nonempty}",
            f"position_restraints_exists  : {position_restraints_exists}",
            f"position_restraints_nonempty: {position_restraints_nonempty}",
            f"protonation_error_message   : {error_message or '<none>'}",
        ],
    )

    if protonation_success:
        _screen_result("success", "gmx pdb2gmx completed")
    else:
        _screen_result("failed", f"returncode={result.returncode}")

    output_dict: dict[str, str | bool | float | int] = {
        "protonation_success": protonation_success,
        "protonation_method": "gmx pdb2gmx -ignh",
        "protonation_input_source": input_source,
        "protonation_input_path": str(input_pdb),
        "protonation_output_path": str(output_pdb),
        "protonation_topology_path": str(topology_output_path),
        "protonation_position_restraints_path": str(position_restraints_output_path),
        "protonation_ph": ph,
        "protonation_ph_applied_by_gromacs": False,
        "protonation_propka_ran": propka_ran,
        "protonation_propka_applied": propka_applied,
        "protonation_propka_his_total": len(his_assignments),
        "protonation_propka_hip_count": propka_him_count,
        "protonation_propka_his_assignments": json.dumps(his_assignments) if his_assignments else "",
        "protonation_state_policy": (
            f"PROPKA-guided HIS states at pH {ph}; gmx pdb2gmx -ignh"
            if propka_ran
            else "PROPKA unavailable; gmx pdb2gmx -ignh with GROMACS defaults"
        ),
        "gromacs_force_field": ff,
        "gromacs_water_model": water_model,
        "gromacs_chain_separation": chain_separation,
        "gromacs_merge_mode": merge,
        "gromacs_ignore_input_hydrogens": ignore_input_hydrogens,
        "input_atom_count": input_atom_count,
        "output_atom_count": output_atom_count,
        "atom_count_increased": atom_count_increased,
        "input_hetatm_count": int(diagnostics["hetatm_count"]),
        "input_ter_count": int(diagnostics["ter_count"]),
        "gmx_pdb2gmx_returncode": int(result.returncode),
        "gmx_pdb2gmx_stdout": result.stdout,
        "gmx_pdb2gmx_stderr": result.stderr,
        "protonation_stdout_log_path": str(stdout_log_path),
        "protonation_stderr_log_path": str(stderr_log_path),
        "protonation_error_message": error_message,
        "topology_exists": topology_exists,
        "topology_nonempty": topology_nonempty,
        "position_restraints_exists": position_restraints_exists,
        "position_restraints_nonempty": position_restraints_nonempty,
    }

    if variant_label is not None:
        output_dict["protonation_variant"] = variant_label

    return output_dict


def protonate_variant_structure(
    pdb_id: str,
    protein_dir: str | Path,
    *,
    variant_label: str,
    input_pdb: str | Path,
    input_source: ProtonationInputSource,
    ph: float = 7.4,
    ff: str = DEFAULT_GROMACS_FORCE_FIELD,
    water_model: GromacsWaterModel = DEFAULT_GROMACS_WATER_MODEL,
) -> dict[str, str | bool | float | int]:
    """Protonate one named variant with GROMACS."""

    output_pdb = build_variant_protonation_output_path(
        pdb_id=pdb_id,
        protein_dir=protein_dir,
        variant_label=variant_label,
    )

    module_log_path = _get_module_log_path(pdb_id=pdb_id, protein_dir=protein_dir)
    _append_log_block(
        module_log_path,
        "protonate_variant_structure",
        [
            f"pdb_id        : {pdb_id}",
            f"protein_dir   : {protein_dir}",
            f"variant_label : {variant_label}",
            f"input_pdb     : {input_pdb}",
            f"output_pdb    : {output_pdb}",
            f"input_source  : {input_source}",
            f"ph            : {ph}",
            f"ff            : {ff}",
            f"water_model   : {water_model}",
        ],
    )

    return protonate_selected_structure(
        input_pdb=input_pdb,
        output_pdb=output_pdb,
        input_source=input_source,
        pdb_id=pdb_id,
        protein_dir=protein_dir,
        ph=ph,
        ff=ff,
        water_model=water_model,
        chain_separation=DEFAULT_GROMACS_CHAIN_SEPARATION,
        merge=DEFAULT_GROMACS_MERGE_MODE,
        ignore_input_hydrogens=True,
        variant_label=variant_label,
    )


def protonate_protein_structure(
    pdb_id: str,
    protein_dir: str | Path,
    *,
    filler_final_model_path: str | Path | None = None,
    filler_output_dir: str | Path | None = None,
    modeller_model_path: str | Path | None = None,
    alphafold_model_path: str | Path | None = None,
    ph: float = 7.4,
    ff: str = DEFAULT_GROMACS_FORCE_FIELD,
    water_model: GromacsWaterModel = DEFAULT_GROMACS_WATER_MODEL,
) -> dict[str, str | bool | float | int]:
    """Backward-compatible single-path protonation workflow."""

    protein_dir = Path(protein_dir)
    output_pdb = build_standard_protonation_output_path(
        pdb_id=pdb_id,
        protein_dir=protein_dir,
    )
    module_log_path = _get_module_log_path(pdb_id=pdb_id, protein_dir=protein_dir)

    _append_log_block(
        module_log_path,
        "protonate_protein_structure:selection_inputs",
        [
            f"pdb_id                   : {pdb_id}",
            f"protein_dir              : {protein_dir}",
            f"filler_final_model_path  : {filler_final_model_path}",
            f"filler_output_dir        : {filler_output_dir}",
            f"modeller_model_path      : {modeller_model_path}",
            f"alphafold_model_path     : {alphafold_model_path}",
            f"ph                       : {ph}",
            f"ff                       : {ff}",
            f"water_model              : {water_model}",
            f"output_pdb               : {output_pdb}",
        ],
    )

    input_pdb, input_source = select_protonation_input(
        pdb_id=pdb_id,
        protein_dir=protein_dir,
        filler_final_model_path=filler_final_model_path,
        filler_output_dir=filler_output_dir,
        modeller_model_path=modeller_model_path,
        alphafold_model_path=alphafold_model_path,
    )

    _append_log_block(
        module_log_path,
        "protonate_protein_structure:selected_input",
        [
            f"selected_input_pdb    : {input_pdb}",
            f"selected_input_source : {input_source}",
        ],
    )

    return protonate_selected_structure(
        input_pdb=input_pdb,
        output_pdb=output_pdb,
        input_source=input_source,
        pdb_id=pdb_id,
        protein_dir=protein_dir,
        ph=ph,
        ff=ff,
        water_model=water_model,
        chain_separation=DEFAULT_GROMACS_CHAIN_SEPARATION,
        merge=DEFAULT_GROMACS_MERGE_MODE,
        ignore_input_hydrogens=True,
    )
