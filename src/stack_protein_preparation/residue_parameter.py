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

import json
import re
import shutil
import subprocess
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any, Literal

from Bio.PDB import PDBIO, PDBParser
from Bio.PDB.PDBIO import Select

ChargeMethod = Literal["resp", "bcc"]
AtomTypeSet = Literal["gaff", "gaff2"]
CapBackend = Literal["pre_capped", "ptmpsi", "none"]

CAP_RESIDUE_NAMES = {"ACE", "NME", "NMA"}


@dataclass(frozen=True, slots=True)
class ResidueSpec:
    """Identify one residue in a PDB structure.

    ``chain_id`` is normalized so blank chain identifiers are represented as
    ``"_"`` in file names and logs. ``residue_number`` is the author residue
    number read from the PDB coordinate record. ``residue_name`` should be the
    final residue name that Amber should learn, such as ``PTR``, ``SEP``, or a
    project-specific three-letter code. ``insertion_code`` is optional and is
    needed only when the same residue number appears more than once.
    """

    chain_id: str
    residue_number: int
    residue_name: str
    insertion_code: str = ""


@dataclass(frozen=True, slots=True)
class GaussianSettings:
    """Configuration for the Gaussian input written for the capped model.

    The default route follows the traditional AMBER RESP-style model commonly
    used for capped residues: HF/6-31G* with electrostatic-potential output. The
    module does not run Gaussian because that is normally a licensed external
    step and may require a cluster scheduler. The generated input is meant to be
    reviewable before submission, so memory, processor count, route line, title,
    charge, and multiplicity are all explicit in the written file. The same
    capped PDB coordinates are used for the coordinate block.
    """

    route_line: str = "#P HF/6-31G* Opt Pop=MK IOp(6/33=2) SCF=Tight"
    nprocshared: int = 8
    mem: str = "8GB"
    title: str = "ACE-RES-NME capped residue model for RESP charge derivation"
    checkpoint_filename: str = "ace_res_nme.chk"


@dataclass(frozen=True, slots=True)
class CommandResult:
    """Store one external command execution result.

    The full command and working directory are recorded so the generated
    manifest can be used for reproducibility. Standard output and standard error
    are written to files and only the paths are kept here. This avoids storing
    long AmberTools logs directly in the JSON manifest while still preserving
    all information needed to debug a failed run. A non-zero return code is left
    for the caller to interpret because the module may want to return a partial
    result with all generated artefacts.
    """

    command: tuple[str, ...]
    cwd: Path
    returncode: int
    stdout_path: Path
    stderr_path: Path


@dataclass(frozen=True, slots=True)
class ResidueParameterResult:
    """Describe all artefacts produced for one residue parameter candidate.

    The result separates the capped model used for QM/charge derivation from the
    uncapped residue template used for protein insertion. ``gaussian_input_path``
    is produced before any AmberTools step and remains the primary output when
    ``charge_method='resp'`` and no Gaussian log is supplied. ``ac_path``,
    ``prepin_path``, ``frcmod_path``, and optional test topology files are
    produced only when the required external tools complete. ``needs_revision``
    reports whether ``parmchk2`` emitted ``ATTN``/``NEEDS REVISION`` warnings.
    """

    status: str
    pdb_id: str
    residue_spec: ResidueSpec
    output_dir: Path
    extracted_residue_pdb_path: Path
    capped_model_pdb_path: Path
    normalized_capped_model_pdb_path: Path
    gaussian_input_path: Path
    gaussian_output_path: Path | None
    mainchain_file_path: Path
    ac_path: Path | None
    mol2_path: Path | None
    prepin_path: Path | None
    frcmod_path: Path | None
    leap_input_path: Path | None
    test_prm7_path: Path | None
    test_rst7_path: Path | None
    manifest_path: Path
    model_charge: int
    residue_charge: float
    charge_method: str
    atom_type_set: str
    needs_revision: bool
    biosimspace_validation: str
    command_results: tuple[CommandResult, ...]
    error_message: str = ""


class _SingleResidueSelect(Select):
    """BioPython selector for one residue from one chain in the first model."""

    def __init__(
        self,
        *,
        model_id: Any,
        chain_id: str,
        residue_number: int,
        insertion_code: str,
    ) -> None:
        self._model_id = model_id
        self._chain_id = _normalize_chain_id(chain_id)
        self._residue_number = residue_number
        self._insertion_code = insertion_code.strip()

    def accept_model(self, model: Any) -> bool:
        return model.id == self._model_id

    def accept_chain(self, chain: Any) -> bool:
        return _normalize_chain_id(str(chain.id)) == self._chain_id

    def accept_residue(self, residue: Any) -> bool:
        residue_number = int(residue.id[1])
        insertion_code = str(residue.id[2]).strip()
        return (
            residue_number == self._residue_number
            and insertion_code == self._insertion_code
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
    biosimspace_validation = "not_requested"

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
                biosimspace_validation=biosimspace_validation,
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
            biosimspace_validation = _validate_with_biosimspace(
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
        biosimspace_validation=biosimspace_validation,
        command_results=tuple(command_results),
        error_message=error_message,
    )
    _write_manifest(result)
    return result


def _normalize_chain_id(chain_id: str | None) -> str:
    return (chain_id or "").strip() or "_"


def _residue_name(residue_spec: ResidueSpec) -> str:
    return residue_spec.residue_name.strip().upper()


def _safe_label(pdb_id: str, residue_spec: ResidueSpec) -> str:
    chain = _normalize_chain_id(residue_spec.chain_id)
    icode = residue_spec.insertion_code.strip()
    suffix = f"{residue_spec.residue_number}{icode}" if icode else str(residue_spec.residue_number)
    raw = f"{pdb_id}_{_residue_name(residue_spec)}_chain_{chain}_{suffix}"
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", raw)


def _write_single_residue_pdb(
    *,
    input_pdb_path: Path,
    output_pdb_path: Path,
    residue_spec: ResidueSpec,
) -> None:
    """Write the selected residue as a small PDB file.

    The selected residue is not the final parametrization model. It is an
    inspectable starting point and a fallback input for the capping backend. The
    function uses BioPython selection rather than manual PDB line slicing so
    insertion codes, chain IDs, and first-model semantics remain explicit. A
    missing or empty output is treated as an error because all downstream steps
    require coordinates.
    """

    if not input_pdb_path.exists():
        raise FileNotFoundError(f"Input PDB not found: {input_pdb_path}")

    parser = PDBParser(QUIET=True)
    structure = parser.get_structure(input_pdb_path.stem, str(input_pdb_path))
    models = list(structure.get_models())
    if not models:
        raise ValueError(f"No coordinate model found in {input_pdb_path}")

    output_pdb_path.parent.mkdir(parents=True, exist_ok=True)
    io = PDBIO()
    io.set_structure(structure)
    io.save(
        str(output_pdb_path),
        select=_SingleResidueSelect(
            model_id=models[0].id,
            chain_id=residue_spec.chain_id,
            residue_number=residue_spec.residue_number,
            insertion_code=residue_spec.insertion_code,
        ),
        preserve_atom_numbering=True,
    )

    if not output_pdb_path.exists() or output_pdb_path.stat().st_size == 0:
        raise RuntimeError(f"Could not write selected residue PDB: {output_pdb_path}")


def _prepare_capped_model(
    *,
    extracted_residue_pdb_path: Path,
    capped_model_pdb_path: Path,
    residue_spec: ResidueSpec,
    cap_backend: CapBackend,
    pre_capped_pdb_path: str | Path | None,
) -> None:
    """Prepare the ACE-RES-NME model used for Gaussian and AmberTools.

    ``pre_capped`` is the safest scientific route because the caller can inspect
    the ACE/NME geometry before quantum chemistry. ``ptmpsi`` is provided as an
    automation backend when PTM-Psi is installed and the single-residue extract
    is sufficient for capping. ``none`` copies the residue without caps and is
    intended only for debugging or cases where the caller has a non-peptide
    residue model. The function writes one canonical capped-model path so later
    steps do not need to know which backend was used.
    """

    capped_model_pdb_path.parent.mkdir(parents=True, exist_ok=True)

    if cap_backend == "pre_capped":
        if pre_capped_pdb_path is None:
            raise ValueError("pre_capped_pdb_path is required for cap_backend='pre_capped'.")
        source = Path(pre_capped_pdb_path).resolve()
        if not source.exists():
            raise FileNotFoundError(f"Pre-capped PDB not found: {source}")
        shutil.copyfile(source, capped_model_pdb_path)
        return

    if cap_backend == "ptmpsi":
        try:
            from ptmpsi.protein import Protein
        except ImportError as exc:
            raise ImportError(
                "PTM-Psi is required for cap_backend='ptmpsi'. Install PTM-Psi "
                "or provide a pre-capped ACE-RES-NME PDB."
            ) from exc

        protein = Protein(
            filename=str(extracted_residue_pdb_path),
            interactive=False,
            delwat=True,
            delhet=False,
        )
        chain_id = _normalize_chain_id(residue_spec.chain_id)
        protein.prepend(chain_id, "ACE")
        protein.append(chain_id, "NME")
        protein.write_pdb(str(capped_model_pdb_path))
        return

    if cap_backend == "none":
        shutil.copyfile(extracted_residue_pdb_path, capped_model_pdb_path)
        return

    raise ValueError(f"Unsupported cap_backend: {cap_backend!r}")


def _normalize_capped_model_atom_names(
    *,
    input_pdb_path: Path,
    output_pdb_path: Path,
    residue_spec: ResidueSpec,
) -> tuple[str, ...]:
    """Rename cap atoms conservatively and return atom names to omit in prepgen.

    prepgen's ``OMIT_NAME`` entries refer to atom names in the Antechamber AC
    file. Atom names in capped peptides can collide between ACE, RES, and NME,
    for example multiple atoms named ``C`` or ``O``. This helper rewrites cap
    atom names to short, predictable names so the generated main-chain file can
    omit the caps without accidentally omitting central-residue atoms. Central
    residue atom names are left unchanged because they must match the protein
    PDB as closely as possible.
    """

    lines = input_pdb_path.read_text(encoding="utf-8").splitlines()
    output_lines: list[str] = []
    omit_names: list[str] = []
    cap_counter = 0

    for line in lines:
        if not line.startswith(("ATOM  ", "HETATM")):
            output_lines.append(line)
            continue

        resname = line[17:20].strip().upper()

        if resname in CAP_RESIDUE_NAMES:
            cap_counter += 1
            prefix = "A" if resname == "ACE" else "M"
            new_atom_name = f"{prefix}{cap_counter:03d}"[-4:]
            omit_names.append(new_atom_name)
            output_lines.append(_replace_pdb_atom_name(line, new_atom_name))
        else:
            output_lines.append(line)

    if not omit_names and _contains_cap_residue(lines):
        raise RuntimeError("Cap residues were detected, but no cap atom names were collected.")

    output_pdb_path.write_text("\n".join(output_lines).rstrip() + "\n", encoding="utf-8")
    return tuple(omit_names)


def _contains_cap_residue(lines: list[str]) -> bool:
    for line in lines:
        if line.startswith(("ATOM  ", "HETATM")) and line[17:20].strip().upper() in CAP_RESIDUE_NAMES:
            return True
    return False


def _replace_pdb_atom_name(line: str, atom_name: str) -> str:
    return f"{line[:12]}{atom_name:>4}{line[16:]}"


def _write_gaussian_input(
    *,
    capped_pdb_path: Path,
    gaussian_input_path: Path,
    settings: GaussianSettings,
    charge: int,
    multiplicity: int,
) -> None:
    """Write a Gaussian input file from the capped model coordinates.

    The coordinate block is generated from PDB ``ATOM`` and ``HETATM`` records
    so the file can be submitted directly after manual inspection. Elements are
    read from the PDB element column when present and otherwise inferred from the
    atom name. The function does not attempt to add hydrogens, optimize
    geometry, or assign protonation states; those choices must be reflected in
    the capped input model. This keeps the module deterministic and avoids
    hidden chemical edits.
    """

    atom_rows = _read_gaussian_atom_rows(capped_pdb_path)
    if not atom_rows:
        raise ValueError(f"No atoms found for Gaussian input: {capped_pdb_path}")

    lines = [
        f"%chk={settings.checkpoint_filename}",
        f"%nprocshared={settings.nprocshared}",
        f"%mem={settings.mem}",
        settings.route_line,
        "",
        settings.title,
        "",
        f"{charge} {multiplicity}",
    ]

    for element, x, y, z in atom_rows:
        lines.append(f"{element:<2} {x:>14.6f} {y:>14.6f} {z:>14.6f}")

    lines.extend(["", ""])
    gaussian_input_path.write_text("\n".join(lines), encoding="utf-8")


def _read_gaussian_atom_rows(pdb_path: Path) -> list[tuple[str, float, float, float]]:
    rows: list[tuple[str, float, float, float]] = []
    for line in pdb_path.read_text(encoding="utf-8").splitlines():
        if not line.startswith(("ATOM  ", "HETATM")):
            continue

        element = line[76:78].strip()
        if not element:
            element = _infer_element_from_atom_name(line[12:16].strip())

        x = float(line[30:38])
        y = float(line[38:46])
        z = float(line[46:54])
        rows.append((element.capitalize(), x, y, z))

    return rows


def _infer_element_from_atom_name(atom_name: str) -> str:
    cleaned = re.sub(r"[^A-Za-z]", "", atom_name).upper()
    if not cleaned:
        return "X"

    two_letter = cleaned[:2].capitalize()
    if two_letter in {
        "Cl",
        "Br",
        "Na",
        "Mg",
        "Al",
        "Si",
        "Ca",
        "Fe",
        "Zn",
        "Cu",
        "Mn",
        "Co",
        "Ni",
        "Se",
        "Li",
        "Be",
    }:
        return two_letter

    return cleaned[0].upper()


def _write_prepgen_mainchain_file(
    *,
    output_path: Path,
    cap_atom_names: tuple[str, ...],
    residue_charge: float,
) -> None:
    """Write the prepgen main-chain definition for cap removal.

    The generated file defines ``N`` and ``C`` as the peptide connection atoms
    for the uncapped residue template and ``CA`` as the main-chain atom. It then
    omits every cap atom name collected from the normalized capped model. prepgen
    redistributes charges from omitted atoms so that the resulting residue has
    the requested total residue charge. The defaults are appropriate for a
    normal internal alpha-amino-acid residue; unusual backbones should be handled
    by providing or editing this file manually.
    """

    lines = [
        "HEAD_NAME N",
        "TAIL_NAME C",
        "MAIN_CHAIN CA",
    ]

    for atom_name in cap_atom_names:
        lines.append(f"OMIT_NAME {atom_name}")

    lines.extend(
        [
            "PRE_HEAD_TYPE C",
            "POST_TAIL_TYPE N",
            f"CHARGE {residue_charge:.6f}",
            "",
        ]
    )

    output_path.write_text("\n".join(lines), encoding="utf-8")


def _run_antechamber_outputs(
    *,
    working_dir: Path,
    residue_spec: ResidueSpec,
    normalized_capped_model_pdb_path: Path,
    gaussian_output_path: Path | None,
    ac_path: Path,
    mol2_path: Path,
    model_charge: int,
    charge_method: ChargeMethod,
    atom_type_set: AtomTypeSet,
) -> tuple[CommandResult, ...]:
    results: list[CommandResult] = []
    residue_name = _residue_name(residue_spec)

    if charge_method == "resp":
        if gaussian_output_path is None:
            raise ValueError("gaussian_output_path is required for charge_method='resp'.")
        if not gaussian_output_path.exists():
            raise FileNotFoundError(f"Gaussian output log not found: {gaussian_output_path}")
        input_path = gaussian_output_path
        input_format = "gout"
        charge_flag = "resp"
    else:
        input_path = normalized_capped_model_pdb_path
        input_format = "pdb"
        charge_flag = "bcc"

    for output_path, output_format, tag in (
        (ac_path, "ac", "ac"),
        (mol2_path, "mol2", "mol2"),
    ):
        results.append(
            _run_command(
                command=(
                    _require_executable("antechamber"),
                    "-i",
                    str(input_path),
                    "-fi",
                    input_format,
                    "-o",
                    output_path.name,
                    "-fo",
                    output_format,
                    "-rn",
                    residue_name,
                    "-c",
                    charge_flag,
                    "-nc",
                    str(model_charge),
                    "-at",
                    atom_type_set,
                ),
                cwd=working_dir,
                stdout_path=working_dir / f"antechamber_{tag}_stdout.log",
                stderr_path=working_dir / f"antechamber_{tag}_stderr.log",
            )
        )

    return tuple(results)


def _parmchk2_atom_type_flag(atom_type_set: AtomTypeSet) -> str:
    if atom_type_set == "gaff2":
        return "2"
    return "1"


def _write_tleap_validation_input(
    *,
    output_path: Path,
    residue_spec: ResidueSpec,
    prepin_path: Path,
    frcmod_path: Path,
    test_prm7_path: Path,
    test_rst7_path: Path,
    atom_type_set: AtomTypeSet,
) -> None:
    """Write a small LEaP validation script for the uncapped residue template.

    The validation system is ``sequence { ACE RES NME }`` because that directly
    tests whether the generated residue template links as an internal peptide
    residue. This is not the final protein system and should not be interpreted
    as the production topology. Its purpose is to verify that the prepin and
    frcmod can be loaded and that the residue can be bonded to normal peptide
    caps. A later full-protein LEaP stage should use the same prepin/frcmod
    artefacts in the real structure.
    """

    residue_name = _residue_name(residue_spec)
    leaprc_gaff = "leaprc.gaff2" if atom_type_set == "gaff2" else "leaprc.gaff"

    lines = [
        "source leaprc.protein.ff14SB",
        f"source {leaprc_gaff}",
        f"loadamberprep {prepin_path.name}",
        f"loadamberparams {frcmod_path.name}",
        f"model = sequence {{ ACE {residue_name} NME }}",
        f"saveamberparm model {test_prm7_path.name} {test_rst7_path.name}",
        "quit",
        "",
    ]

    output_path.write_text("\n".join(lines), encoding="utf-8")


def _run_command(
    *,
    command: tuple[str, ...],
    cwd: Path,
    stdout_path: Path,
    stderr_path: Path,
) -> CommandResult:
    stdout_path.parent.mkdir(parents=True, exist_ok=True)
    stderr_path.parent.mkdir(parents=True, exist_ok=True)

    completed = subprocess.run(
        list(command),
        cwd=str(cwd),
        text=True,
        capture_output=True,
        check=False,
    )

    stdout_path.write_text(completed.stdout or "", encoding="utf-8")
    stderr_path.write_text(completed.stderr or "", encoding="utf-8")

    return CommandResult(
        command=tuple(command),
        cwd=cwd,
        returncode=completed.returncode,
        stdout_path=stdout_path,
        stderr_path=stderr_path,
    )


def _require_executable(name: str) -> str:
    executable = shutil.which(name)
    if executable is None:
        raise FileNotFoundError(f"Required executable not found in PATH: {name}")
    return executable


def _raise_if_last_commands_failed(
    command_results: list[CommandResult],
    label: str,
) -> None:
    for result in command_results[-2:]:
        _raise_if_command_failed(result, label)


def _raise_if_command_failed(command_result: CommandResult, label: str) -> None:
    if command_result.returncode == 0:
        return

    stderr_preview = ""
    if command_result.stderr_path.exists():
        stderr_preview = command_result.stderr_path.read_text(encoding="utf-8")[-1000:]

    raise RuntimeError(
        f"{label} failed with return code {command_result.returncode}. "
        f"stderr_log={command_result.stderr_path}. stderr_tail={stderr_preview}"
    )


def _frcmod_needs_revision(frcmod_path: Path) -> bool:
    if not frcmod_path.exists():
        return False

    text = frcmod_path.read_text(encoding="utf-8", errors="replace").upper()
    return "ATTN" in text or "NEEDS REVISION" in text


def _validate_with_biosimspace(*, prm7_path: Path, rst7_path: Path) -> str:
    try:
        import BioSimSpace as BSS  # type: ignore
    except Exception as exc:
        return f"not_available: {type(exc).__name__}: {exc}"

    try:
        _ = BSS.IO.readMolecules([str(rst7_path), str(prm7_path)])
    except Exception as exc:
        return f"failed: {type(exc).__name__}: {exc}"

    return "success"


def _make_result(
    *,
    status: str,
    pdb_id: str,
    residue_spec: ResidueSpec,
    output_dir: Path,
    extracted_residue_pdb_path: Path,
    capped_model_pdb_path: Path,
    normalized_capped_model_pdb_path: Path,
    gaussian_input_path: Path,
    gaussian_output_path: Path | None,
    mainchain_file_path: Path,
    ac_path: Path | None,
    mol2_path: Path | None,
    prepin_path: Path | None,
    frcmod_path: Path | None,
    leap_input_path: Path | None,
    test_prm7_path: Path | None,
    test_rst7_path: Path | None,
    manifest_path: Path,
    model_charge: int,
    residue_charge: float,
    charge_method: str,
    atom_type_set: str,
    needs_revision: bool,
    biosimspace_validation: str,
    command_results: tuple[CommandResult, ...],
    error_message: str,
) -> ResidueParameterResult:
    return ResidueParameterResult(
        status=status,
        pdb_id=pdb_id,
        residue_spec=residue_spec,
        output_dir=output_dir,
        extracted_residue_pdb_path=extracted_residue_pdb_path,
        capped_model_pdb_path=capped_model_pdb_path,
        normalized_capped_model_pdb_path=normalized_capped_model_pdb_path,
        gaussian_input_path=gaussian_input_path,
        gaussian_output_path=gaussian_output_path,
        mainchain_file_path=mainchain_file_path,
        ac_path=ac_path,
        mol2_path=mol2_path,
        prepin_path=prepin_path,
        frcmod_path=frcmod_path,
        leap_input_path=leap_input_path,
        test_prm7_path=test_prm7_path,
        test_rst7_path=test_rst7_path,
        manifest_path=manifest_path,
        model_charge=model_charge,
        residue_charge=residue_charge,
        charge_method=charge_method,
        atom_type_set=atom_type_set,
        needs_revision=needs_revision,
        biosimspace_validation=biosimspace_validation,
        command_results=command_results,
        error_message=error_message,
    )


def _write_manifest(result: ResidueParameterResult) -> None:
    payload = _result_to_json(result)
    result.manifest_path.write_text(
        json.dumps(payload, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )


def _path_to_str(path: Path | None) -> str | None:
    if path is None:
        return None
    return str(path)


def _result_to_json(result: ResidueParameterResult) -> dict[str, Any]:
    return {
        "status": result.status,
        "pdb_id": result.pdb_id,
        "residue_spec": asdict(result.residue_spec),
        "output_dir": str(result.output_dir),
        "extracted_residue_pdb_path": str(result.extracted_residue_pdb_path),
        "capped_model_pdb_path": str(result.capped_model_pdb_path),
        "normalized_capped_model_pdb_path": str(result.normalized_capped_model_pdb_path),
        "gaussian_input_path": str(result.gaussian_input_path),
        "gaussian_output_path": _path_to_str(result.gaussian_output_path),
        "mainchain_file_path": str(result.mainchain_file_path),
        "ac_path": _path_to_str(result.ac_path),
        "mol2_path": _path_to_str(result.mol2_path),
        "prepin_path": _path_to_str(result.prepin_path),
        "frcmod_path": _path_to_str(result.frcmod_path),
        "leap_input_path": _path_to_str(result.leap_input_path),
        "test_prm7_path": _path_to_str(result.test_prm7_path),
        "test_rst7_path": _path_to_str(result.test_rst7_path),
        "manifest_path": str(result.manifest_path),
        "model_charge": result.model_charge,
        "residue_charge": result.residue_charge,
        "charge_method": result.charge_method,
        "atom_type_set": result.atom_type_set,
        "needs_revision": result.needs_revision,
        "biosimspace_validation": result.biosimspace_validation,
        "command_results": [
            {
                "command": list(command_result.command),
                "cwd": str(command_result.cwd),
                "returncode": command_result.returncode,
                "stdout_path": str(command_result.stdout_path),
                "stderr_path": str(command_result.stderr_path),
            }
            for command_result in result.command_results
        ],
        "error_message": result.error_message,
    }


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
