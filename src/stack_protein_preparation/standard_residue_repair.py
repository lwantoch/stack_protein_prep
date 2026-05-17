"""Minimal standard-residue heavy-atom repair for FRUTON protonation inputs.

This module is intentionally narrower than a general PDB fixer. It only attempts
one repair class: missing heavy atoms in existing standard protein residues. It
uses MODELLER's ``complete_pdb`` helper, which builds missing atoms from internal
coordinates, and it refuses to silently handle non-standard amino-acid-like
residues by default. Missing loops, missing residues, PTMs, metals, ligands, and
water placement remain outside this module.

The repair function is meant to run immediately before the sanitizer and
``gmx pdb2gmx``. Gap detection and filler decisions should still use the original
representative monomer so that residue-number discontinuities are not hidden by
any downstream cleanup or atom rebuilding.
"""

from __future__ import annotations

import contextlib
import io
from dataclasses import dataclass
from pathlib import Path

from Bio.PDB import PDBParser

from stack_protein_preparation.sanitize import (
    find_missing_standard_heavy_atoms,
    find_nonstandard_residue_names,
    load_gromacs_heavy_atom_templates,
)


@dataclass(frozen=True)
class StandardResidueRepairResult:
    """Result of the optional MODELLER heavy-atom repair gate."""

    input_pdb_path: Path
    output_pdb_path: Path
    log_path: Path
    repair_attempted: bool
    repair_success: bool
    can_continue: bool
    used_output_path: Path
    missing_heavy_atom_count_before: int
    missing_heavy_atom_count_after: int
    nonstandard_residue_names: tuple[str, ...]
    message: str

    @property
    def status(self) -> str:
        if self.repair_success:
            return "success"
        if self.repair_attempted:
            return "failed"
        if self.missing_heavy_atom_count_before:
            return "warning"
        return "skipped"


def repair_standard_residue_heavy_atoms_for_gromacs(
    *,
    input_pdb_path: str | Path,
    output_pdb_path: str | Path,
    log_path: str | Path,
    force_field: str = "amber99sb-ildn",
    allow_nonstandard_residues: bool = False,
) -> StandardResidueRepairResult:
    """Repair missing heavy atoms in standard residues with MODELLER if needed.

    The function first checks the input against the installed GROMACS RTP
    templates through the sanitizer's existing template reader. If no standard
    heavy atoms are missing, no repair is attempted and the original input is
    returned as the path to use. If non-standard residues are present, repair is
    skipped by default because MODELLER may not understand the residue chemistry
    and FRUTON should route those cases through parameter generation instead.
    """

    input_path = Path(input_pdb_path).resolve()
    output_path = Path(output_pdb_path).resolve()
    resolved_log_path = Path(log_path).resolve()
    resolved_log_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    structure = _parse_structure(input_path)
    if structure is None:
        result = StandardResidueRepairResult(
            input_pdb_path=input_path,
            output_pdb_path=output_path,
            log_path=resolved_log_path,
            repair_attempted=False,
            repair_success=False,
            can_continue=True,
            used_output_path=input_path,
            missing_heavy_atom_count_before=0,
            missing_heavy_atom_count_after=0,
            nonstandard_residue_names=(),
            message="Bio.PDB could not parse input; repair skipped and original input retained.",
        )
        _write_log(result)
        return result

    template_lookup = load_gromacs_heavy_atom_templates(force_field=force_field)
    missing_before = find_missing_standard_heavy_atoms(
        structure,
        template_lookup=template_lookup,
    )
    nonstandard_names = tuple(find_nonstandard_residue_names(structure))

    if not missing_before:
        result = StandardResidueRepairResult(
            input_pdb_path=input_path,
            output_pdb_path=output_path,
            log_path=resolved_log_path,
            repair_attempted=False,
            repair_success=False,
            can_continue=True,
            used_output_path=input_path,
            missing_heavy_atom_count_before=0,
            missing_heavy_atom_count_after=0,
            nonstandard_residue_names=nonstandard_names,
            message="No missing standard heavy atoms detected; repair skipped.",
        )
        _write_log(result)
        return result

    if nonstandard_names and not allow_nonstandard_residues:
        result = StandardResidueRepairResult(
            input_pdb_path=input_path,
            output_pdb_path=output_path,
            log_path=resolved_log_path,
            repair_attempted=False,
            repair_success=False,
            can_continue=True,
            used_output_path=input_path,
            missing_heavy_atom_count_before=len(missing_before),
            missing_heavy_atom_count_after=len(missing_before),
            nonstandard_residue_names=nonstandard_names,
            message=(
                "Missing standard heavy atoms detected, but repair skipped because "
                "non-standard residues are present: " + ", ".join(nonstandard_names)
            ),
        )
        _write_log(result, missing_before=missing_before)
        return result

    try:
        modeller_stdout = _run_modeller_complete_pdb(
            input_pdb_path=input_path,
            output_pdb_path=output_path,
        )
    except Exception as error:
        result = StandardResidueRepairResult(
            input_pdb_path=input_path,
            output_pdb_path=output_path,
            log_path=resolved_log_path,
            repair_attempted=True,
            repair_success=False,
            can_continue=True,
            used_output_path=input_path,
            missing_heavy_atom_count_before=len(missing_before),
            missing_heavy_atom_count_after=len(missing_before),
            nonstandard_residue_names=nonstandard_names,
            message=f"MODELLER complete_pdb repair failed: {error!r}; original input retained.",
        )
        _write_log(result, missing_before=missing_before)
        return result

    repaired_structure = _parse_structure(output_path)
    if repaired_structure is None:
        result = StandardResidueRepairResult(
            input_pdb_path=input_path,
            output_pdb_path=output_path,
            log_path=resolved_log_path,
            repair_attempted=True,
            repair_success=False,
            can_continue=True,
            used_output_path=input_path,
            missing_heavy_atom_count_before=len(missing_before),
            missing_heavy_atom_count_after=len(missing_before),
            nonstandard_residue_names=nonstandard_names,
            message="MODELLER wrote output, but Bio.PDB could not parse it; original input retained.",
        )
        _write_log(result, missing_before=missing_before, modeller_stdout=modeller_stdout)
        return result

    missing_after = find_missing_standard_heavy_atoms(
        repaired_structure,
        template_lookup=template_lookup,
    )
    repair_success = output_path.exists() and output_path.stat().st_size > 0 and len(missing_after) < len(missing_before)
    used_output_path = output_path if repair_success else input_path
    result = StandardResidueRepairResult(
        input_pdb_path=input_path,
        output_pdb_path=output_path,
        log_path=resolved_log_path,
        repair_attempted=True,
        repair_success=repair_success,
        can_continue=True,
        used_output_path=used_output_path,
        missing_heavy_atom_count_before=len(missing_before),
        missing_heavy_atom_count_after=len(missing_after),
        nonstandard_residue_names=nonstandard_names,
        message=(
            "MODELLER complete_pdb repair applied."
            if repair_success
            else "MODELLER complete_pdb did not reduce missing-heavy-atom diagnostics; original input retained."
        ),
    )
    _write_log(
        result,
        missing_before=missing_before,
        missing_after=missing_after,
        modeller_stdout=modeller_stdout,
    )
    return result


def _run_modeller_complete_pdb(*, input_pdb_path: Path, output_pdb_path: Path) -> str:
    """Run MODELLER complete_pdb and write a completed PDB."""

    try:
        from modeller import environ
        from modeller.scripts import complete_pdb
    except Exception as error:
        raise ImportError(
            "MODELLER is required for minimal standard-heavy-atom repair."
        ) from error

    buffer = io.StringIO()
    with contextlib.redirect_stdout(buffer), contextlib.redirect_stderr(buffer):
        env = environ()
        model = complete_pdb(env, str(input_pdb_path))
        model.write(file=str(output_pdb_path))
    return buffer.getvalue()


def _parse_structure(pdb_path: Path):
    parser = PDBParser(QUIET=True, PERMISSIVE=True)
    try:
        return parser.get_structure(pdb_path.stem, str(pdb_path))
    except Exception:
        return None


def _write_log(
    result: StandardResidueRepairResult,
    *,
    missing_before: list[str] | None = None,
    missing_after: list[str] | None = None,
    modeller_stdout: str = "",
) -> None:
    lines = [
        "standard_residue_repair",
        "=======================",
        "",
        f"input_pdb_path                   : {result.input_pdb_path}",
        f"output_pdb_path                  : {result.output_pdb_path}",
        f"used_output_path                 : {result.used_output_path}",
        f"repair_attempted                 : {result.repair_attempted}",
        f"repair_success                   : {result.repair_success}",
        f"can_continue                     : {result.can_continue}",
        f"missing_heavy_atom_count_before  : {result.missing_heavy_atom_count_before}",
        f"missing_heavy_atom_count_after   : {result.missing_heavy_atom_count_after}",
        f"nonstandard_residue_names        : {', '.join(result.nonstandard_residue_names) or '<none>'}",
        f"status                           : {result.status}",
        f"message                          : {result.message}",
        "",
        "missing_before",
        "--------------",
    ]
    if missing_before:
        lines.extend(missing_before)
    else:
        lines.append("<none>")
    lines.extend(["", "missing_after", "-------------"])
    if missing_after:
        lines.extend(missing_after)
    else:
        lines.append("<none>")
    if modeller_stdout:
        lines.extend(["", "modeller_output", "---------------", modeller_stdout])
    result.log_path.write_text("\n".join(lines) + "\n", encoding="utf-8")