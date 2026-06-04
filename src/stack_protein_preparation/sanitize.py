"""PDB sanitizing gates for GROMACS-facing FRUTON stages.

The module prepares protein-like PDB inputs for ``gmx pdb2gmx`` without trying
to solve scientific modeling problems inside the sanitizer. It uses Bio.PDB for
PDB parsing, disorder handling, residue traversal, and PDB writing rather than
maintaining a hand-written PDB object model. The only intentionally custom part
is the GROMACS force-field template check: the module reads the installed
``aminoacids.rtp`` file so missing-heavy-atom diagnostics match the force field
that ``pdb2gmx`` will use.

The sanitizer is meant to sit before protonation of full protein variants and
temporary gap fragments. It selects one alternate location per disordered atom,
normalizes selected occupancies to 1.00, removes non-polymer heterogens from the
protein-only route, writes a clean PDB, and reports what cleanup and diagnostics
were observed before ``pdb2gmx``.  It does not build missing residues, it does
not protonate, and it does not generate parameters for non-standard residues or
metals.

Optionally (``repair_missing_sidechain_heavy_atoms=True``) it can fill in
missing standard-residue side-chain heavy atoms using PDBFixer in controlled
mode.  Backbone heavy atoms (N, CA, C, O), nonstandard residues, and any
residue whose atoms collide with nearby heavy atoms are blocked and reported
rather than silently fixed.  The repair runs *after* the BioPython cleanup
write so altloc selection, occupancy normalisation, and heterogen removal are
already applied before PDBFixer sees the file.
"""

from __future__ import annotations

from pathlib import Path

from Bio.PDB import PDBParser

# Re-export everything from _sanitize_core so existing imports keep working
from stack_protein_preparation._sanitize_core import (  # noqa: F401 (re-exported)
    IssueSeverity,
    SanitizeIssue,
    HeavyAtomTemplateLookup,
    SanitizeResult,
    GromacsProteinSelect,
    STANDARD_RESIDUE_TEMPLATE_ALIASES,
    AltlocSelectionSummary,
    select_altlocs_by_highest_occupancy,
    count_nonunit_occupancies,
    find_nonstandard_residue_names,
    count_written_residues,
    load_gromacs_heavy_atom_templates,
    find_missing_standard_heavy_atoms,
    find_gromacs_aminoacids_rtp_path,
    parse_gromacs_rtp_heavy_atoms,
    _BACKBONE_HEAVY,
    _CLASH_WARN_DIST,
    _RepairResult,
    _repair_sidechain_heavy_atoms,
    _write_structure_with_biopython,
    _iter_residue_child_atoms,
    _altloc_sort_key,
    _resolve_expected_heavy_atoms,
    _is_amino_acid_like,
    _is_standard_amino_acid_name,
    _is_hydrogen_atom,
    _resolve_gromacs_executable,
    _read_gromacs_data_prefix,
    _result_for_missing_input,
    _result_for_parse_failure,
    _write_sanitize_log,
)


def sanitize_pdb_for_gromacs(
    input_pdb_path: str | Path,
    output_pdb_path: str | Path,
    *,
    pdb_id: str | None = None,
    variant_label: str | None = None,
    force_field: str = "amber99sb-ildn",
    gmx_executable: str | None = None,
    log_path: str | Path | None = None,
    drop_nonpolymer_heterogens: bool = True,
    normalize_occupancy: bool = True,
    fail_on_missing_heavy_atoms: bool = False,
    fail_on_nonstandard_residues: bool = False,
    repair_missing_sidechain_heavy_atoms: bool = False,
    min_residue_count_warning: int = 3,
) -> SanitizeResult:
    """Clean and classify one PDB before the FRUTON GROMACS protonation route.

    The function uses Bio.PDB to parse the structure, select alternate
    locations, and write the sanitized output. It then classifies the result
    against GROMACS force-field residue templates so FRUTON can decide whether
    to call ``pdb2gmx``. The design deliberately separates simple formatting
    cleanup from repair: missing heavy atoms are reported, not constructed, and
    non-standard residues are reported, not parameterized.

    Parameters
    ----------
    input_pdb_path
        PDB file to sanitize.
    output_pdb_path
        Cleaned PDB path. The file is written even when the result is not
        ``pdb2gmx``-ready, because that file is useful for inspection.
    pdb_id
        Optional PDB identifier used only for log context.
    variant_label
        Optional FRUTON variant label, for example ``single`` or
        ``gaps_fragment_02``.
    force_field
        GROMACS force-field name without the ``.ff`` suffix.
    gmx_executable
        Optional GROMACS command. If omitted, ``gmx`` and then ``gmx_mpi`` are
        searched.
    log_path
        Optional plain-text log path.
    drop_nonpolymer_heterogens
        Remove waters, ligands, ions, and metals from the protein-only route.
        Amino-acid-like hetero residues are retained so the parameter audit can
        still detect PTMs such as PTR.
    normalize_occupancy
        Write selected atoms with occupancy 1.00.
    fail_on_missing_heavy_atoms
        Treat missing standard heavy atoms as hard blockers when explicitly requested.
    fail_on_nonstandard_residues
        Treat non-standard residues as hard blockers. This usually stays false
        because FRUTON should let the parameter-audit policy decide PTM handling.
    repair_missing_sidechain_heavy_atoms
        When True, attempt to fill in missing standard-residue side-chain heavy
        atoms using PDBFixer after the BioPython cleanup write.  Backbone atoms
        (N, CA, C, O), nonstandard residues, and clash-producing additions are
        blocked and recorded as issues rather than silently applied.  Original
        atom positions must be unchanged after repair; any drift is treated as
        an error and the repaired file is discarded.  Requires PDBFixer and
        OpenMM to be installed.
    min_residue_count_warning
        Emit a warning for very small fragments. This does not block by itself.
    """

    input_path = Path(input_pdb_path).resolve()
    output_path = Path(output_pdb_path).resolve()
    resolved_log_path = Path(log_path).resolve() if log_path is not None else None

    issues: list[SanitizeIssue] = []

    if not input_path.exists():
        result = _result_for_missing_input(
            input_path=input_path,
            output_path=output_path,
            log_path=resolved_log_path,
        )
        _write_sanitize_log(result, pdb_id=pdb_id, variant_label=variant_label)
        return result

    parser = PDBParser(QUIET=True, PERMISSIVE=True)

    try:
        structure = parser.get_structure(pdb_id or input_path.stem, str(input_path))
    except Exception as error:
        issues.append(
            SanitizeIssue(
                severity="error",
                code="pdb_parse_failed",
                message=f"Bio.PDB could not parse input PDB: {error!r}",
            )
        )
        result = _result_for_parse_failure(
            input_path=input_path,
            output_path=output_path,
            log_path=resolved_log_path,
            issues=issues,
        )
        _write_sanitize_log(result, pdb_id=pdb_id, variant_label=variant_label)
        return result

    atom_count_in = sum(1 for _atom in structure.get_atoms())
    if atom_count_in == 0:
        issues.append(
            SanitizeIssue(
                severity="error",
                code="no_atoms",
                message="Input PDB contains no atoms after Bio.PDB parsing.",
            )
        )

    altloc_summary = select_altlocs_by_highest_occupancy(structure)
    if altloc_summary.selected_altloc_count:
        issues.append(
            SanitizeIssue(
                severity="warning",
                code="altloc_selected",
                message=(
                    f"Selected one alternate location for "
                    f"{altloc_summary.selected_altloc_count} disordered atoms; "
                    f"dropped {altloc_summary.dropped_altloc_atom_count} duplicate positions."
                ),
            )
        )

    nonunit_occupancy_count = count_nonunit_occupancies(structure)
    if nonunit_occupancy_count:
        issues.append(
            SanitizeIssue(
                severity="warning",
                code="occupancy_normalized",
                message=(
                    f"{nonunit_occupancy_count} selected atoms had occupancy other "
                    "than 1.00 and will be written as 1.00."
                ),
            )
        )

    nonstandard_residue_names = find_nonstandard_residue_names(structure)
    if nonstandard_residue_names:
        issues.append(
            SanitizeIssue(
                severity="warning",
                code="nonstandard_residues_present",
                message=(
                    "Non-standard amino-acid-like residues are present: "
                    + ", ".join(nonstandard_residue_names)
                ),
            )
        )

    if fail_on_nonstandard_residues and nonstandard_residue_names:
        issues.append(
            SanitizeIssue(
                severity="error",
                code="nonstandard_residues_block_gromacs",
                message="Current sanitizer policy blocks non-standard residues before pdb2gmx.",
            )
        )

    template_lookup = load_gromacs_heavy_atom_templates(
        force_field=force_field,
        gmx_executable=gmx_executable,
    )
    if template_lookup.rtp_path is None:
        issues.append(
            SanitizeIssue(
                severity="warning",
                code="force_field_template_missing",
                message=(
                    f"Could not find {force_field}.ff/aminoacids.rtp; "
                    "standard-residue heavy-atom completeness was not checked."
                ),
            )
        )

    missing_heavy_atom_messages = find_missing_standard_heavy_atoms(
        structure,
        template_lookup=template_lookup,
    )
    for message in missing_heavy_atom_messages:
        issues.append(
            SanitizeIssue(
                severity="error" if fail_on_missing_heavy_atoms else "warning",
                code="missing_standard_heavy_atoms",
                message=message,
            )
        )

    residue_count = count_written_residues(
        structure,
        drop_nonpolymer_heterogens=drop_nonpolymer_heterogens,
    )
    if 0 < residue_count < min_residue_count_warning:
        issues.append(
            SanitizeIssue(
                severity="warning",
                code="tiny_fragment",
                message=(
                    f"Sanitized structure has only {residue_count} polymer-like residues; "
                    "this may be an unintended fragment."
                ),
            )
        )

    output_path.parent.mkdir(parents=True, exist_ok=True)
    selector = GromacsProteinSelect(
        drop_nonpolymer_heterogens=drop_nonpolymer_heterogens,
        normalize_occupancy=normalize_occupancy,
    )
    _write_structure_with_biopython(structure, output_path=output_path, selector=selector)

    if selector.dropped_heterogen_atom_count:
        issues.append(
            SanitizeIssue(
                severity="info",
                code="heterogens_dropped",
                message=(
                    f"Dropped {selector.dropped_heterogen_atom_count} non-polymer "
                    "HETATM atoms from the protein-only route."
                ),
            )
        )

    # --- Optional PDBFixer sidechain heavy-atom repair ----------------------
    repaired_heavy_atom_count = 0
    blocked_repair_count = 0
    if repair_missing_sidechain_heavy_atoms and missing_heavy_atom_messages:
        repair = _repair_sidechain_heavy_atoms(
            pdb_path=output_path,
            template_lookup=template_lookup,
        )
        repaired_heavy_atom_count = repair.n_repaired
        blocked_repair_count = repair.n_blocked
        for msg in repair.repaired:
            issues.append(SanitizeIssue(
                severity="info",
                code="heavy_atom_repaired",
                message=f"PDBFixer added missing sidechain atoms: {msg}",
            ))
        for msg in repair.blocked:
            issues.append(SanitizeIssue(
                severity="warning",
                code="heavy_atom_repair_blocked",
                message=f"Repair blocked (not applied): {msg}",
            ))
        if repair.error:
            issues.append(SanitizeIssue(
                severity="error",
                code="heavy_atom_repair_failed",
                message=f"PDBFixer repair failed, original file kept: {repair.error}",
            ))

    hard_error_present = any(issue.severity == "error" for issue in issues)

    result = SanitizeResult(
        input_pdb_path=input_path,
        output_pdb_path=output_path,
        log_path=resolved_log_path,
        output_written=output_path.exists() and output_path.stat().st_size > 0,
        can_run_gromacs=not hard_error_present,
        atom_count_in=atom_count_in,
        atom_count_out=selector.written_atom_count,
        residue_count_out=residue_count,
        selected_altloc_count=altloc_summary.selected_altloc_count,
        dropped_altloc_atom_count=altloc_summary.dropped_altloc_atom_count,
        normalized_occupancy_count=nonunit_occupancy_count,
        dropped_heterogen_atom_count=selector.dropped_heterogen_atom_count,
        nonstandard_residue_names=tuple(nonstandard_residue_names),
        missing_heavy_atom_count=len(missing_heavy_atom_messages),
        repaired_heavy_atom_count=repaired_heavy_atom_count,
        blocked_repair_count=blocked_repair_count,
        force_field_template_path=template_lookup.rtp_path,
        issues=tuple(issues),
    )
    _write_sanitize_log(result, pdb_id=pdb_id, variant_label=variant_label)
    return result
