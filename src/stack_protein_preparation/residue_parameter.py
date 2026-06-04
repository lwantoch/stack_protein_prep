"""Build AMBER residue templates for non-standard amino-acid residues.

This module prepares a capped model peptide for one non-standard residue and
uses AmberTools to derive a reusable residue template for insertion into a
protein. The chemical model used for charge derivation is ACE-RES-NME, because
an internal amino-acid residue should be parametrized in a peptide-like
environment rather than as an isolated molecule. The final insertion artefact is
not the capped model itself: prepgen receives a main-chain definition file that
removes ACE/NME atoms and writes an uncapped RES template with N/C peptide
connection atoms.

The module intentionally orchestrates established external tools rather than
implementing force-field chemistry locally. BioPython is used only for coordinate
selection and PDB writing, PTM-Psi may optionally build the ACE/NME capped model,
Gaussian input is written as a normal text input file, and AmberTools programs
perform atom typing, charge consumption, template generation, parameter checking,
and LEaP validation. A Gaussian log is optional at the first stage because the
normal workflow is two-pass: first write the Gaussian input, then rerun the same
function after Gaussian has produced the output log.
"""

from __future__ import annotations

import shutil
from pathlib import Path

# Re-export everything from _residue_parameter_core so existing imports keep working
from stack_protein_preparation._residue_parameter_core import (  # noqa: F401 (re-exported)
    ChargeMethod,
    AtomTypeSet,
    CapBackend,
    CAP_RESIDUE_NAMES,
    ResidueSpec,
    GaussianSettings,
    CommandResult,
    ResidueParameterResult,
    _SingleResidueSelect,
    _normalize_chain_id,
    _residue_name,
    _safe_label,
    _write_single_residue_pdb,
    _prepare_capped_model,
    _normalize_capped_model_atom_names,
    _contains_cap_residue,
    _replace_pdb_atom_name,
    _write_gaussian_input,
    _read_gaussian_atom_rows,
    _infer_element_from_atom_name,
    _write_prepgen_mainchain_file,
    _run_antechamber_outputs,
    _parmchk2_atom_type_flag,
    _write_tleap_validation_input,
    _run_command,
    _require_executable,
    _raise_if_last_commands_failed,
    _raise_if_command_failed,
    _frcmod_needs_revision,
    _validate_with_biosimspace,
    _make_result,
    _write_manifest,
    _path_to_str,
    _result_to_json,
)


def build_residue_parameter_candidate(
    *,
    pdb_id: str,
    input_pdb_path: str | Path,
    residue_spec: ResidueSpec,
    output_root_dir: str | Path,
    model_charge: int,
    residue_charge: float,
    multiplicity: int = 1,
    cap_backend: CapBackend = "pre_capped",
    pre_capped_pdb_path: str | Path | None = None,
    gaussian_output_path: str | Path | None = None,
    gaussian_settings: GaussianSettings | None = None,
    charge_method: ChargeMethod = "resp",
    atom_type_set: AtomTypeSet = "gaff2",
    run_tleap: bool = True,
    validate_with_biosimspace: bool = False,
    overwrite: bool = True,
) -> ResidueParameterResult:
    """Build a parameter candidate for one non-standard residue.

    The first pass always writes a residue extract, an ACE-RES-NME capped model
    or a copy of the caller-supplied capped model, a normalized capped PDB, a
    Gaussian input file, and a prepgen main-chain definition file. If
    ``charge_method`` is ``"resp"`` and ``gaussian_output_path`` is not supplied,
    the function stops with status ``"requires_gaussian"`` because the Gaussian
    output log is the input for RESP charge extraction through Antechamber. If a
    Gaussian log is supplied, the function continues through Antechamber,
    prepgen, parmchk2, and optionally tleap. The function returns a structured
    result and always writes a JSON manifest, including failed and partial runs.
    """

    normalized_pdb_id = pdb_id.strip().upper()
    input_pdb_path = Path(input_pdb_path).resolve()
    output_root_dir = Path(output_root_dir).resolve()
    gaussian_settings = gaussian_settings or GaussianSettings()
    gaussian_output = (
        Path(gaussian_output_path).resolve()
        if gaussian_output_path is not None
        else None
    )

    label = _safe_label(normalized_pdb_id, residue_spec)
    output_dir = output_root_dir / label

    if overwrite and output_dir.exists():
        shutil.rmtree(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    extracted_residue_pdb_path = output_dir / "extracted_residue.pdb"
    capped_model_pdb_path = output_dir / "ace_res_nme_model.pdb"
    normalized_capped_model_pdb_path = output_dir / "ace_res_nme_model.named.pdb"
    gaussian_input_path = output_dir / "ace_res_nme_gaussian.gjf"
    mainchain_file_path = output_dir / f"{_residue_name(residue_spec)}.mc"
    manifest_path = output_dir / "manifest.json"

    ac_path = output_dir / f"{_residue_name(residue_spec)}.ac"
    mol2_path = output_dir / f"{_residue_name(residue_spec)}.mol2"
    prepin_path = output_dir / f"{_residue_name(residue_spec)}.prepin"
    frcmod_path = output_dir / f"{_residue_name(residue_spec)}.frcmod"
    leap_input_path = output_dir / "tleap_validate_residue.in"
    test_prm7_path = output_dir / "test_ace_res_nme.prm7"
    test_rst7_path = output_dir / "test_ace_res_nme.rst7"

    command_results: list[CommandResult] = []
    status = "success"
    error_message = ""
    needs_revision = False
    biosimspace_validation_result = "not_requested"

    try:
        _write_single_residue_pdb(
            input_pdb_path=input_pdb_path,
            output_pdb_path=extracted_residue_pdb_path,
            residue_spec=residue_spec,
        )

        _prepare_capped_model(
            extracted_residue_pdb_path=extracted_residue_pdb_path,
            capped_model_pdb_path=capped_model_pdb_path,
            residue_spec=residue_spec,
            cap_backend=cap_backend,
            pre_capped_pdb_path=pre_capped_pdb_path,
        )

        cap_atom_names = _normalize_capped_model_atom_names(
            input_pdb_path=capped_model_pdb_path,
            output_pdb_path=normalized_capped_model_pdb_path,
            residue_spec=residue_spec,
        )

        _write_gaussian_input(
            capped_pdb_path=normalized_capped_model_pdb_path,
            gaussian_input_path=gaussian_input_path,
            settings=gaussian_settings,
            charge=model_charge,
            multiplicity=multiplicity,
        )

        _write_prepgen_mainchain_file(
            output_path=mainchain_file_path,
            cap_atom_names=cap_atom_names,
            residue_charge=residue_charge,
        )

        if charge_method == "resp" and gaussian_output is None:
            status = "requires_gaussian"
            error_message = (
                "Gaussian input was written. Provide gaussian_output_path on the "
                "next run to continue with RESP/Antechamber."
            )
            result = _make_result(
                status=status,
                pdb_id=normalized_pdb_id,
                residue_spec=residue_spec,
                output_dir=output_dir,
                extracted_residue_pdb_path=extracted_residue_pdb_path,
                capped_model_pdb_path=capped_model_pdb_path,
                normalized_capped_model_pdb_path=normalized_capped_model_pdb_path,
                gaussian_input_path=gaussian_input_path,
                gaussian_output_path=gaussian_output,
                mainchain_file_path=mainchain_file_path,
                ac_path=None,
                mol2_path=None,
                prepin_path=None,
                frcmod_path=None,
                leap_input_path=None,
                test_prm7_path=None,
                test_rst7_path=None,
                manifest_path=manifest_path,
                model_charge=model_charge,
                residue_charge=residue_charge,
                charge_method=charge_method,
                atom_type_set=atom_type_set,
                needs_revision=False,
                biosimspace_validation=biosimspace_validation_result,
                command_results=tuple(command_results),
                error_message=error_message,
            )
            _write_manifest(result)
            return result

        command_results.extend(
            _run_antechamber_outputs(
                working_dir=output_dir,
                residue_spec=residue_spec,
                normalized_capped_model_pdb_path=normalized_capped_model_pdb_path,
                gaussian_output_path=gaussian_output,
                ac_path=ac_path,
                mol2_path=mol2_path,
                model_charge=model_charge,
                charge_method=charge_method,
                atom_type_set=atom_type_set,
            )
        )

        _raise_if_last_commands_failed(command_results, "antechamber")

        command_results.append(
            _run_command(
                command=(
                    _require_executable("prepgen"),
                    "-i",
                    ac_path.name,
                    "-o",
                    prepin_path.name,
                    "-m",
                    mainchain_file_path.name,
                    "-rn",
                    _residue_name(residue_spec),
                ),
                cwd=output_dir,
                stdout_path=output_dir / "prepgen_stdout.log",
                stderr_path=output_dir / "prepgen_stderr.log",
            )
        )
        _raise_if_command_failed(command_results[-1], "prepgen")

        command_results.append(
            _run_command(
                command=(
                    _require_executable("parmchk2"),
                    "-i",
                    prepin_path.name,
                    "-f",
                    "prepi",
                    "-o",
                    frcmod_path.name,
                    "-s",
                    _parmchk2_atom_type_flag(atom_type_set),
                ),
                cwd=output_dir,
                stdout_path=output_dir / "parmchk2_stdout.log",
                stderr_path=output_dir / "parmchk2_stderr.log",
            )
        )
        _raise_if_command_failed(command_results[-1], "parmchk2")

        # Re-run parmchk2 on a prepin→mol2 conversion so that bonds at the
        # ff14SB/GAFF2 interface (e.g. 2C-s4 after prepgen renames backbone atoms)
        # get parameters. parmchk2 -f prepi only queries the GAFF2 database and
        # silently skips ff14SB atom types, leaving cross-type bonds unparametrized.
        _prepin_as_mol2_path = output_dir / f"{_residue_name(residue_spec)}_from_prepin.mol2"
        command_results.append(
            _run_command(
                command=(
                    _require_executable("antechamber"),
                    "-i",
                    prepin_path.name,
                    "-fi",
                    "prepi",
                    "-o",
                    _prepin_as_mol2_path.name,
                    "-fo",
                    "mol2",
                    "-at",
                    atom_type_set,
                ),
                cwd=output_dir,
                stdout_path=output_dir / "antechamber_prepin_mol2_stdout.log",
                stderr_path=output_dir / "antechamber_prepin_mol2_stderr.log",
            )
        )
        if command_results[-1].returncode == 0 and _prepin_as_mol2_path.exists():
            command_results.append(
                _run_command(
                    command=(
                        _require_executable("parmchk2"),
                        "-i",
                        _prepin_as_mol2_path.name,
                        "-f",
                        "mol2",
                        "-o",
                        frcmod_path.name,
                        "-s",
                        _parmchk2_atom_type_flag(atom_type_set),
                    ),
                    cwd=output_dir,
                    stdout_path=output_dir / "parmchk2_mol2_stdout.log",
                    stderr_path=output_dir / "parmchk2_mol2_stderr.log",
                )
            )
            _raise_if_command_failed(command_results[-1], "parmchk2 (mol2 re-run)")

        needs_revision = _frcmod_needs_revision(frcmod_path)
        if needs_revision:
            status = "warning"

        if run_tleap:
            _write_tleap_validation_input(
                output_path=leap_input_path,
                residue_spec=residue_spec,
                prepin_path=prepin_path,
                frcmod_path=frcmod_path,
                test_prm7_path=test_prm7_path,
                test_rst7_path=test_rst7_path,
                atom_type_set=atom_type_set,
            )
            command_results.append(
                _run_command(
                    command=(
                        _require_executable("tleap"),
                        "-f",
                        leap_input_path.name,
                    ),
                    cwd=output_dir,
                    stdout_path=output_dir / "tleap_stdout.log",
                    stderr_path=output_dir / "tleap_stderr.log",
                )
            )
            _raise_if_command_failed(command_results[-1], "tleap")

        if validate_with_biosimspace and test_prm7_path.exists() and test_rst7_path.exists():
            biosimspace_validation_result = _validate_with_biosimspace(
                prm7_path=test_prm7_path,
                rst7_path=test_rst7_path,
            )

    except Exception as exc:
        status = "failed"
        error_message = f"{type(exc).__name__}: {exc}"

    result = _make_result(
        status=status,
        pdb_id=normalized_pdb_id,
        residue_spec=residue_spec,
        output_dir=output_dir,
        extracted_residue_pdb_path=extracted_residue_pdb_path,
        capped_model_pdb_path=capped_model_pdb_path,
        normalized_capped_model_pdb_path=normalized_capped_model_pdb_path,
        gaussian_input_path=gaussian_input_path,
        gaussian_output_path=gaussian_output,
        mainchain_file_path=mainchain_file_path,
        ac_path=ac_path if ac_path.exists() else None,
        mol2_path=mol2_path if mol2_path.exists() else None,
        prepin_path=prepin_path if prepin_path.exists() else None,
        frcmod_path=frcmod_path if frcmod_path.exists() else None,
        leap_input_path=leap_input_path if leap_input_path.exists() else None,
        test_prm7_path=test_prm7_path if test_prm7_path.exists() else None,
        test_rst7_path=test_rst7_path if test_rst7_path.exists() else None,
        manifest_path=manifest_path,
        model_charge=model_charge,
        residue_charge=residue_charge,
        charge_method=charge_method,
        atom_type_set=atom_type_set,
        needs_revision=needs_revision,
        biosimspace_validation=biosimspace_validation_result,
        command_results=tuple(command_results),
        error_message=error_message,
    )
    _write_manifest(result)
    return result


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description=(
            "Prepare Gaussian input and AmberTools artefacts for one "
            "ACE-RES-NME capped non-standard residue."
        )
    )
    parser.add_argument("pdb_id")
    parser.add_argument("input_pdb", type=Path)
    parser.add_argument("chain_id")
    parser.add_argument("residue_number", type=int)
    parser.add_argument("residue_name")
    parser.add_argument("output_root_dir", type=Path)
    parser.add_argument("--insertion-code", default="")
    parser.add_argument("--model-charge", type=int, required=True)
    parser.add_argument("--residue-charge", type=float, required=True)
    parser.add_argument("--multiplicity", type=int, default=1)
    parser.add_argument(
        "--cap-backend",
        choices=["pre_capped", "ptmpsi", "none"],
        default="pre_capped",
    )
    parser.add_argument("--pre-capped-pdb", type=Path)
    parser.add_argument("--gaussian-output", type=Path)
    parser.add_argument("--charge-method", choices=["resp", "bcc"], default="resp")
    parser.add_argument("--atom-type-set", choices=["gaff", "gaff2"], default="gaff2")
    parser.add_argument("--no-tleap", action="store_true")
    args = parser.parse_args()

    result = build_residue_parameter_candidate(
        pdb_id=args.pdb_id,
        input_pdb_path=args.input_pdb,
        residue_spec=ResidueSpec(
            chain_id=args.chain_id,
            residue_number=args.residue_number,
            residue_name=args.residue_name,
            insertion_code=args.insertion_code,
        ),
        output_root_dir=args.output_root_dir,
        model_charge=args.model_charge,
        residue_charge=args.residue_charge,
        multiplicity=args.multiplicity,
        cap_backend=args.cap_backend,
        pre_capped_pdb_path=args.pre_capped_pdb,
        gaussian_output_path=args.gaussian_output,
        charge_method=args.charge_method,
        atom_type_set=args.atom_type_set,
        run_tleap=not args.no_tleap,
    )

    print(f"status={result.status}")
    print(f"gaussian_input={result.gaussian_input_path}")
    print(f"manifest={result.manifest_path}")
    if result.error_message:
        print(f"error={result.error_message}")
