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

import json
import re
import shutil
import subprocess
import tempfile
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any, Literal

from Bio.PDB import PDBParser
from Bio.PDB.Polypeptide import is_aa

MODULE_NAME = "parameter_audit"

AuditStatus = Literal["success", "warning", "failed"]
IssueCategory = Literal[
    "unsupported_residue",
    "missing_heavy_atom",
    "atom_name_mismatch",
    "terminal_patch_problem",
    "metal_or_special_bond_problem",
    "gromacs_runtime_error",
    "environment_error",
]
ResidueCategory = Literal[
    "standard",
    "forcefield_alias",
    "canonicalizable",
    "protein_like_nonstandard",
    "metal",
    "other_hetero",
]

DEFAULT_GROMACS_FORCE_FIELD = "amber99sb-ildn"
DEFAULT_GROMACS_WATER_MODEL = "tip3p"
DEFAULT_GROMACS_CHAIN_SEPARATION = "id_or_ter"
DEFAULT_GROMACS_MERGE_MODE = "no"
DEFAULT_GROMACS_TIMEOUT_SECONDS = 600

STANDARD_AMINO_ACID_RESNAMES = {
    "ALA",
    "ARG",
    "ASN",
    "ASP",
    "CYS",
    "GLN",
    "GLU",
    "GLY",
    "HIS",
    "ILE",
    "LEU",
    "LYS",
    "MET",
    "PHE",
    "PRO",
    "SER",
    "THR",
    "TRP",
    "TYR",
    "VAL",
}

FORCE_FIELD_ALIAS_RESNAMES = {
    "ASH",
    "GLH",
    "HID",
    "HIE",
    "HIP",
    "CYM",
    "CYX",
    "LYN",
    "NALA",
    "NARG",
    "NASN",
    "NASP",
    "NASH",
    "NCYS",
    "NCYM",
    "NCYX",
    "NGLN",
    "NGLU",
    "NGLH",
    "NGLY",
    "NHIS",
    "NHID",
    "NHIE",
    "NHIP",
    "NILE",
    "NLEU",
    "NLYS",
    "NLYN",
    "NMET",
    "NPHE",
    "NPRO",
    "NSER",
    "NTHR",
    "NTRP",
    "NTYR",
    "NVAL",
    "CALA",
    "CARG",
    "CASN",
    "CASP",
    "CASH",
    "CCYS",
    "CCYM",
    "CCYX",
    "CGLN",
    "CGLU",
    "CGLH",
    "CGLY",
    "CHIS",
    "CHID",
    "CHIE",
    "CHIP",
    "CILE",
    "CLEU",
    "CLYS",
    "CLYN",
    "CMET",
    "CPHE",
    "CPRO",
    "CSER",
    "CTHR",
    "CTRP",
    "CTYR",
    "CVAL",
}

CANONICALIZABLE_RESNAMES = {
    "MSE": "MET",
}

COMMON_PROTEIN_LIKE_NONSTANDARD_RESNAMES = {
    "SEC",
    "PYL",
    "SEP",
    "TPO",
    "PTR",
    "CSO",
    "CSD",
    "OCS",
    "KCX",
    "LLP",
    "MLY",
    "HYP",
}

METAL_RESNAMES = {
    "LI",
    "NA",
    "K",
    "RB",
    "CS",
    "MG",
    "CA",
    "SR",
    "BA",
    "ZN",
    "FE",
    "CU",
    "MN",
    "CO",
    "NI",
    "CR",
    "AL",
    "AG",
    "AU",
    "HG",
    "CD",
}

WATER_RESNAMES = {
    "HOH",
    "WAT",
    "H2O",
    "TIP",
    "TIP3",
    "TIP3P",
    "SOL",
}


@dataclass(frozen=True, slots=True)
class ResidueSite:
    """Describe one residue site observed in the audited coordinate file.

    The site stores only stable residue metadata needed for reporting and policy
    decisions. It intentionally does not expose BioPython residue objects because
    those are parser-specific and awkward to serialize. ``category`` is a local
    FRUTON classification, not a guarantee that GROMACS can parameterize the
    residue. The GROMACS probe result remains the authoritative compatibility
    decision.
    """

    chain_id: str
    resname: str
    resseq: int
    icode: str
    category: ResidueCategory
    atom_names: tuple[str, ...]

    @property
    def label(self) -> str:
        insertion = self.icode.strip()
        residue_number = f"{self.resseq}{insertion}" if insertion else str(self.resseq)
        return f"{self.chain_id}:{self.resname}{residue_number}"


@dataclass(frozen=True, slots=True)
class ParameterAuditIssue:
    """Represent one force-field or input-quality issue found by the audit.

    Issues are mostly parsed from the GROMACS ``pdb2gmx`` stderr text. The parser
    is deliberately conservative and captures only high-value categories that are
    useful for branching the pipeline. The original message is retained so the
    caller can inspect the exact GROMACS wording in the JSON manifest or log.
    """

    category: IssueCategory
    resname: str | None
    residue_number: int | None
    atom_name: str | None
    message: str


@dataclass(frozen=True, slots=True)
class GromacsProbeResult:
    """Store the raw outcome of the temporary ``gmx pdb2gmx`` compatibility run.

    The probe is executed in a temporary directory and should not be used as a
    production topology. Successful probe output only means the current GROMACS
    force field accepted the input structure with the requested options. The real
    protonation/preparation stage should still run separately in the normal
    FRUTON output directory.
    """

    command: tuple[str, ...]
    cwd: Path
    returncode: int | None
    success: bool
    stdout: str
    stderr: str
    timeout_seconds: int
    error_message: str


@dataclass(frozen=True, slots=True)
class ParameterAuditResult:
    """Complete parameter-audit result for one protein coordinate file.

    ``gromacs_probe_success`` is the main compatibility signal. Non-standard
    residue lists are reported separately so FRUTON can distinguish harmless
    aliases or canonicalizable residues from residues that need external
    parameters. Missing heavy atoms are treated as repair problems rather than
    as immediate Gaussian/QM-parameterization problems. Metal residues are listed
    separately because they usually require a dedicated metal-center workflow
    instead of generic small-molecule parametrization.
    """

    pdb_id: str
    input_pdb_path: Path
    force_field: str
    water_model: str
    gromacs_executable: str
    status: AuditStatus
    gromacs_probe_success: bool
    can_use_current_force_field: bool
    requires_repair: bool
    requires_external_parameters: bool
    requires_qm_parameters: bool
    requires_metal_parameters: bool
    residue_sites: tuple[ResidueSite, ...]
    nonstandard_residue_sites: tuple[ResidueSite, ...]
    metal_sites: tuple[ResidueSite, ...]
    supported_nonstandard_residue_names: tuple[str, ...]
    unsupported_residue_names: tuple[str, ...]
    issues: tuple[ParameterAuditIssue, ...]
    probe: GromacsProbeResult
    manifest_path: Path | None
    log_path: Path | None

    def to_json_dict(self) -> dict[str, Any]:
        """Return a JSON-serializable representation of the audit result."""

        payload = asdict(self)
        payload["input_pdb_path"] = str(self.input_pdb_path)
        payload["manifest_path"] = str(self.manifest_path) if self.manifest_path else None
        payload["log_path"] = str(self.log_path) if self.log_path else None
        payload["probe"]["cwd"] = str(self.probe.cwd)
        payload["probe"]["command"] = list(self.probe.command)
        return payload


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


def collect_residue_sites(input_pdb_path: str | Path) -> tuple[ResidueSite, ...]:
    """Collect residue sites from the first model of a PDB file.

    BioPython is used for parsing so chain, residue, insertion-code, and atom
    information come from a tested structural parser. Water residues are skipped
    because they are not audited as protein residue templates. Hetero residues
    are retained when they look metal-like or protein-like, since those are the
    cases that can affect the downstream protein force-field route.
    """

    input_pdb_path = Path(input_pdb_path).resolve()
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure(input_pdb_path.stem, str(input_pdb_path))
    models = list(structure.get_models())

    if not models:
        raise ValueError(f"No coordinate model found in input PDB: {input_pdb_path}")

    residue_sites: list[ResidueSite] = []

    for chain in models[0].get_chains():
        chain_id = _normalize_chain_id(str(chain.id))

        for residue in chain.get_residues():
            resname = residue.get_resname().strip().upper()
            if resname in WATER_RESNAMES:
                continue

            residue_number = int(residue.id[1])
            insertion_code = str(residue.id[2]).strip()
            atom_names = tuple(
                sorted({str(atom.id).strip().upper() for atom in residue.get_atoms()})
            )
            category = classify_residue_site(residue)

            if category == "other_hetero":
                continue

            residue_sites.append(
                ResidueSite(
                    chain_id=chain_id,
                    resname=resname,
                    resseq=residue_number,
                    icode=insertion_code,
                    category=category,
                    atom_names=atom_names,
                )
            )

    return tuple(residue_sites)


def classify_residue_site(residue: Any) -> ResidueCategory:
    """Classify one BioPython residue for reporting and pipeline branching.

    This classification is intentionally not a force-field compatibility test.
    It only assigns a useful human-readable category before the GROMACS probe is
    run. The probe remains the authoritative test because residue names, atom
    names, termini, missing heavy atoms, and special-bond logic interact inside
    GROMACS.
    """

    resname = residue.get_resname().strip().upper()
    element_names = {
        str(getattr(atom, "element", "")).strip().upper()
        for atom in residue.get_atoms()
    }

    if resname in METAL_RESNAMES or any(element in METAL_RESNAMES for element in element_names):
        return "metal"

    if resname in STANDARD_AMINO_ACID_RESNAMES:
        return "standard"

    if resname in FORCE_FIELD_ALIAS_RESNAMES:
        return "forcefield_alias"

    if resname in CANONICALIZABLE_RESNAMES:
        return "canonicalizable"

    if resname in COMMON_PROTEIN_LIKE_NONSTANDARD_RESNAMES:
        return "protein_like_nonstandard"

    if is_aa(residue, standard=False):
        return "protein_like_nonstandard"

    return "other_hetero"


def run_gromacs_parameter_probe(
    *,
    input_pdb_path: Path,
    gromacs_executable: str,
    force_field: str,
    water_model: str,
    chain_separation: str,
    merge: str,
    timeout_seconds: int,
    keep_probe_directory: bool = False,
) -> GromacsProbeResult:
    """Run a temporary ``gmx pdb2gmx`` probe and return its raw result.

    The probe is a compatibility check only. It runs in a temporary directory so
    generated topologies and coordinate files cannot collide with normal FRUTON
    artefacts. ``GMX_MAXBACKUP=-1`` avoids backup-file prompts and keeps the
    command deterministic in repeated audit runs.
    """

    tmp_parent = None if keep_probe_directory else tempfile.TemporaryDirectory(
        prefix="fruton_parameter_audit_"
    )

    if tmp_parent is None:
        tmp_dir = Path(tempfile.mkdtemp(prefix="fruton_parameter_audit_"))
    else:
        tmp_dir = Path(tmp_parent.name)

    command = (
        gromacs_executable,
        "pdb2gmx",
        "-f",
        str(input_pdb_path),
        "-o",
        "probe_output.pdb",
        "-p",
        "probe_topol.top",
        "-i",
        "probe_posre.itp",
        "-ff",
        force_field,
        "-water",
        water_model,
        "-chainsep",
        chain_separation,
        "-merge",
        merge,
        "-ignh",
    )

    env = dict(**{key: value for key, value in __import__("os").environ.items()})
    env["GMX_MAXBACKUP"] = "-1"

    try:
        completed = subprocess.run(
            command,
            cwd=tmp_dir,
            input="",
            text=True,
            capture_output=True,
            timeout=timeout_seconds,
            check=False,
            env=env,
        )
        result = GromacsProbeResult(
            command=command,
            cwd=tmp_dir,
            returncode=completed.returncode,
            success=completed.returncode == 0,
            stdout=completed.stdout,
            stderr=completed.stderr,
            timeout_seconds=timeout_seconds,
            error_message="" if completed.returncode == 0 else _stderr_tail(completed.stderr),
        )
    except subprocess.TimeoutExpired as exc:
        result = GromacsProbeResult(
            command=command,
            cwd=tmp_dir,
            returncode=None,
            success=False,
            stdout=str(exc.stdout or ""),
            stderr=str(exc.stderr or ""),
            timeout_seconds=timeout_seconds,
            error_message=f"gmx pdb2gmx timed out after {timeout_seconds} seconds",
        )
    except OSError as exc:
        result = GromacsProbeResult(
            command=command,
            cwd=tmp_dir,
            returncode=None,
            success=False,
            stdout="",
            stderr="",
            timeout_seconds=timeout_seconds,
            error_message=f"{type(exc).__name__}: {exc}",
        )

    if tmp_parent is not None:
        tmp_parent.cleanup()

    return result


def parse_gromacs_probe_issues(
    stderr_text: str,
    fallback_error_message: str = "",
) -> tuple[ParameterAuditIssue, ...]:
    """Parse high-value failure categories from GROMACS stderr text.

    The parser intentionally avoids covering every possible GROMACS error. It
    only captures categories that are useful for FRUTON branching: unsupported
    residue templates, missing heavy atoms, atom-name mismatches, terminal patch
    problems, and special-bond/metal warnings. Unknown failures are kept as a
    generic runtime issue so the full stderr remains available in the audit log.
    """

    text = stderr_text or ""
    issues: list[ParameterAuditIssue] = []

    missing_atom_pattern = re.compile(
        r"Residue\s+(?P<residue_number>\d+)\s+named\s+(?P<resname>[A-Za-z0-9]+).*?"
        r"atom\s+(?P<atom_name>[A-Za-z0-9]+)\s+used.*?not found",
        re.IGNORECASE | re.DOTALL,
    )
    for match in missing_atom_pattern.finditer(text):
        issues.append(
            ParameterAuditIssue(
                category="missing_heavy_atom",
                resname=match.group("resname").upper(),
                residue_number=int(match.group("residue_number")),
                atom_name=match.group("atom_name").upper(),
                message=_compact_message(match.group(0)),
            )
        )

    unsupported_patterns = [
        re.compile(
            r"Residue\s+'?(?P<resname>[A-Za-z0-9]+)'?\s+not\s+found\s+in\s+residue\s+topology\s+database",
            re.IGNORECASE,
        ),
        re.compile(
            r"There\s+is\s+no\s+residue\s+type\s+for\s+'?(?P<resname>[A-Za-z0-9]+)'?",
            re.IGNORECASE,
        ),
        re.compile(
            r"could\s+not\s+find\s+residue\s+'?(?P<resname>[A-Za-z0-9]+)'?",
            re.IGNORECASE,
        ),
    ]
    for pattern in unsupported_patterns:
        for match in pattern.finditer(text):
            issues.append(
                ParameterAuditIssue(
                    category="unsupported_residue",
                    resname=match.group("resname").upper(),
                    residue_number=None,
                    atom_name=None,
                    message=_compact_message(match.group(0)),
                )
            )

    if "atom" in text.lower() and "not found in" in text.lower() and not issues:
        issues.append(
            ParameterAuditIssue(
                category="atom_name_mismatch",
                resname=None,
                residue_number=None,
                atom_name=None,
                message=_stderr_tail(text),
            )
        )

    lower_text = text.lower()
    if "termin" in lower_text and "fatal error" in lower_text and not issues:
        issues.append(
            ParameterAuditIssue(
                category="terminal_patch_problem",
                resname=None,
                residue_number=None,
                atom_name=None,
                message=_stderr_tail(text),
            )
        )

    if any(token in lower_text for token in ("specbond", "metal", "coordination")) and "fatal error" in lower_text and not issues:
        issues.append(
            ParameterAuditIssue(
                category="metal_or_special_bond_problem",
                resname=None,
                residue_number=None,
                atom_name=None,
                message=_stderr_tail(text),
            )
        )

    if not issues and ("fatal error" in lower_text or fallback_error_message):
        issues.append(
            ParameterAuditIssue(
                category="gromacs_runtime_error",
                resname=None,
                residue_number=None,
                atom_name=None,
                message=_stderr_tail(text) or fallback_error_message,
            )
        )

    return tuple(_deduplicate_issues(issues))


def write_parameter_audit_manifest(
    *,
    result: ParameterAuditResult,
    manifest_path: str | Path,
) -> Path:
    """Write a JSON manifest for one audit result."""

    manifest_path = Path(manifest_path)
    manifest_path.parent.mkdir(parents=True, exist_ok=True)
    manifest_path.write_text(
        json.dumps(result.to_json_dict(), indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    return manifest_path


def write_parameter_audit_log(
    *,
    result: ParameterAuditResult,
    log_path: str | Path,
) -> Path:
    """Write a compact text log containing command, summary, stdout, and stderr."""

    log_path = Path(log_path)
    log_path.parent.mkdir(parents=True, exist_ok=True)

    lines = [
        "=" * 100,
        f"module                  : {MODULE_NAME}",
        f"pdb_id                  : {result.pdb_id}",
        f"input_pdb_path          : {result.input_pdb_path}",
        f"force_field             : {result.force_field}",
        f"water_model             : {result.water_model}",
        f"gromacs_executable      : {result.gromacs_executable}",
        f"status                  : {result.status}",
        f"gromacs_probe_success   : {result.gromacs_probe_success}",
        f"requires_repair         : {result.requires_repair}",
        f"requires_external_params: {result.requires_external_parameters}",
        f"requires_qm_parameters  : {result.requires_qm_parameters}",
        f"requires_metal_params   : {result.requires_metal_parameters}",
        f"nonstandard_count       : {len(result.nonstandard_residue_sites)}",
        f"metal_count             : {len(result.metal_sites)}",
        f"issues_count            : {len(result.issues)}",
        "-" * 100,
        "[COMMAND]",
        " ".join(result.probe.command),
        "-" * 100,
        "[NONSTANDARD RESIDUES]",
    ]

    if result.nonstandard_residue_sites:
        for site in result.nonstandard_residue_sites:
            lines.append(f"{site.label:<20} {site.category}")
    else:
        lines.append("<none>")

    lines.extend(["-" * 100, "[ISSUES]"])
    if result.issues:
        for issue in result.issues:
            lines.append(
                f"{issue.category}: resname={issue.resname or '-'} "
                f"residue_number={issue.residue_number or '-'} "
                f"atom={issue.atom_name or '-'} | {issue.message}"
            )
    else:
        lines.append("<none>")

    lines.extend(
        [
            "-" * 100,
            "[STDOUT]",
            result.probe.stdout.strip() or "<empty>",
            "-" * 100,
            "[STDERR]",
            result.probe.stderr.strip() or result.probe.error_message or "<empty>",
            "",
        ]
    )

    log_path.write_text("\n".join(lines), encoding="utf-8")
    return log_path


def result_to_pipeline_fields(result: ParameterAuditResult) -> dict[str, str]:
    """Return flat string fields suitable for future FRUTON pipeline records.

    The current pipeline state module does not yet define these columns. This
    helper makes the intended integration explicit without importing or mutating
    pipeline state here. A runner can merge the returned mapping only after the
    corresponding columns have been added to ``pipeline_state.py`` and the XLSX
    export.
    """

    return {
        "parameter_audit.status": result.status,
        "parameter_audit.gromacs_probe_success": _yes_no(result.gromacs_probe_success),
        "parameter_audit.requires_repair": _yes_no(result.requires_repair),
        "parameter_audit.requires_external_parameters": _yes_no(
            result.requires_external_parameters
        ),
        "parameter_audit.requires_qm_parameters": _yes_no(result.requires_qm_parameters),
        "parameter_audit.requires_metal_parameters": _yes_no(
            result.requires_metal_parameters
        ),
        "parameter_audit.supported_nonstandard_residues": "|".join(
            result.supported_nonstandard_residue_names
        ),
        "parameter_audit.unsupported_residues": "|".join(
            result.unsupported_residue_names
        ),
        "parameter_audit.issue_categories": "|".join(
            issue.category for issue in result.issues
        ),
    }


def _resolve_pdb_id(*, pdb_id: str | None, input_pdb_path: Path) -> str:
    if pdb_id is not None and str(pdb_id).strip():
        return str(pdb_id).strip().upper()

    for parent in input_pdb_path.parents:
        if parent.name == "components" and parent.parent.name:
            return parent.parent.name.upper()

    stem = input_pdb_path.stem.strip()
    if stem:
        return stem.split("_", maxsplit=1)[0].upper()

    return "UNKNOWN"


def _find_gromacs_executable(explicit_executable: str | None) -> str:
    if explicit_executable is not None and str(explicit_executable).strip():
        executable = str(explicit_executable).strip()
        if shutil.which(executable) or Path(executable).is_file():
            return executable
        raise FileNotFoundError(f"GROMACS executable not found: {executable}")

    for candidate in ("gmx", "gmx_mpi"):
        resolved = shutil.which(candidate)
        if resolved:
            return candidate

    raise FileNotFoundError("Could not find GROMACS executable in PATH: gmx or gmx_mpi")


def _supported_nonstandard_residue_names(
    *,
    nonstandard_sites: tuple[ResidueSite, ...],
    probe_success: bool,
) -> tuple[str, ...]:
    names = {site.resname for site in nonstandard_sites}
    if probe_success:
        return tuple(sorted(names))

    supported = {
        site.resname
        for site in nonstandard_sites
        if site.category in {"canonicalizable", "forcefield_alias"}
    }
    return tuple(sorted(supported))


def _unsupported_residue_names(
    *,
    nonstandard_sites: tuple[ResidueSite, ...],
    issues: tuple[ParameterAuditIssue, ...],
    probe_success: bool,
) -> tuple[str, ...]:
    if probe_success:
        return ()

    issue_names = {
        str(issue.resname).upper()
        for issue in issues
        if issue.category == "unsupported_residue" and issue.resname
    }
    if issue_names:
        return tuple(sorted(issue_names))

    unknown_nonstandard = {
        site.resname
        for site in nonstandard_sites
        if site.category == "protein_like_nonstandard"
    }
    return tuple(sorted(unknown_nonstandard))


def _audit_status(
    *,
    probe_success: bool,
    requires_repair: bool,
    requires_external_parameters: bool,
    requires_metal_parameters: bool,
) -> AuditStatus:
    if probe_success and not requires_metal_parameters:
        return "success"

    if probe_success and requires_metal_parameters:
        return "warning"

    if requires_repair:
        return "warning"

    if requires_external_parameters or requires_metal_parameters:
        return "failed"

    return "failed"


def _normalize_chain_id(chain_id: str | None) -> str:
    return (chain_id or "").strip() or "_"


def _compact_message(message: str) -> str:
    return " ".join(str(message).split())


def _stderr_tail(stderr_text: str, max_lines: int = 20) -> str:
    lines = [line.rstrip() for line in str(stderr_text or "").splitlines() if line.strip()]
    if not lines:
        return ""
    return "\n".join(lines[-max_lines:])


def _deduplicate_issues(
    issues: list[ParameterAuditIssue],
) -> list[ParameterAuditIssue]:
    seen: set[tuple[Any, ...]] = set()
    deduplicated: list[ParameterAuditIssue] = []

    for issue in issues:
        key = (
            issue.category,
            issue.resname,
            issue.residue_number,
            issue.atom_name,
            issue.message,
        )
        if key in seen:
            continue
        seen.add(key)
        deduplicated.append(issue)

    return deduplicated


def _yes_no(value: bool) -> str:
    return "yes" if value else "no"
