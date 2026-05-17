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
protein-only route, writes a clean PDB, and reports what cleanup and diagnostics were observed before ``pdb2gmx``. It does not build missing residues, it does
not add missing heavy atoms, it does not protonate, and it does not generate
parameters for non-standard residues or metals.
"""

from __future__ import annotations

import os
import shutil
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Literal

from Bio.PDB import PDBIO, PDBParser, Select
from Bio.PDB.Polypeptide import is_aa

IssueSeverity = Literal["info", "warning", "error"]


@dataclass(frozen=True)
class SanitizeIssue:
    """One inspectable sanitizer diagnostic.

    ``severity`` decides whether the issue is only informative, needs review, or
    blocks automatic ``pdb2gmx``. ``code`` is intentionally stable and compact so
    it can later be stored in JSON or displayed in the FRUTON workbook.
    ``message`` is explicit enough for a user to understand the failing residue
    or cleanup action without opening the PDB file first.
    """

    severity: IssueSeverity
    code: str
    message: str


@dataclass(frozen=True)
class HeavyAtomTemplateLookup:
    """Expected heavy atoms loaded from one GROMACS force-field RTP file."""

    force_field: str
    rtp_path: Path | None
    heavy_atoms_by_residue: dict[str, frozenset[str]]


@dataclass(frozen=True)
class SanitizeResult:
    """Summary returned by ``sanitize_pdb_for_gromacs``.

    ``output_written`` only says that a cleaned PDB was written. It is not the
    same as ``can_run_gromacs``. A file with missing standard side-chain heavy
    atoms is still useful for inspection, but it should not be passed blindly to
    ``pdb2gmx`` because GROMACS will usually fail later with a less local error.
    """

    input_pdb_path: Path
    output_pdb_path: Path
    log_path: Path | None
    output_written: bool
    can_run_gromacs: bool
    atom_count_in: int
    atom_count_out: int
    residue_count_out: int
    selected_altloc_count: int
    dropped_altloc_atom_count: int
    normalized_occupancy_count: int
    dropped_heterogen_atom_count: int
    nonstandard_residue_names: tuple[str, ...]
    missing_heavy_atom_count: int
    force_field_template_path: Path | None
    issues: tuple[SanitizeIssue, ...]

    @property
    def status(self) -> str:
        """Return ``success``, ``warning``, or ``failed`` for pipeline tables."""

        if not self.can_run_gromacs:
            return "failed"
        if any(issue.severity == "warning" for issue in self.issues):
            return "warning"
        return "success"


class GromacsProteinSelect(Select):
    """Bio.PDB writer selector for the FRUTON protein-only route.

    The selector keeps polymer-like amino-acid residues, including modified
    amino acids that Bio.PDB recognizes as amino-acid-like with
    ``standard=False``. It drops waters, ligands, ions, and other non-polymer
    heterogens when ``drop_nonpolymer_heterogens`` is true. Selected atom
    occupancies are normalized and selected alternate-location labels are
    cleared at write time so the output is less noisy for ``pdb2gmx``.
    """

    def __init__(
        self,
        *,
        drop_nonpolymer_heterogens: bool,
        normalize_occupancy: bool,
    ) -> None:
        self.drop_nonpolymer_heterogens = drop_nonpolymer_heterogens
        self.normalize_occupancy = normalize_occupancy
        self.dropped_heterogen_atom_count = 0
        self.written_atom_count = 0

    def accept_residue(self, residue) -> int:  # noqa: ANN001
        """Return whether a residue should be written."""

        if not self.drop_nonpolymer_heterogens:
            return 1

        hetfield = residue.id[0]
        if hetfield == " ":
            return 1

        if _is_amino_acid_like(residue.get_resname()):
            return 1

        self.dropped_heterogen_atom_count += sum(1 for _atom in residue.get_atoms())
        return 0

    def accept_atom(self, atom) -> int:  # noqa: ANN001
        """Normalize selected atoms immediately before PDBIO writes them."""

        if self.normalize_occupancy:
            atom.set_occupancy(1.0)

        if atom.get_altloc() not in ("", " "):
            atom.set_altloc(" ")

        self.written_atom_count += 1
        return 1


STANDARD_RESIDUE_TEMPLATE_ALIASES: dict[str, tuple[str, ...]] = {
    "HIS": ("HIS", "HISE", "HISD", "HISH", "HIE", "HID", "HIP", "HSD", "HSE", "HSP"),
    "HIE": ("HIE", "HISE", "HIS"),
    "HID": ("HID", "HISD", "HIS"),
    "HIP": ("HIP", "HISH", "HIS"),
    "HSE": ("HSE", "HISE", "HIS"),
    "HSD": ("HSD", "HISD", "HIS"),
    "HSP": ("HSP", "HISH", "HIS"),
    "CYS": ("CYS", "CYX", "CYM"),
}


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
        force_field_template_path=template_lookup.rtp_path,
        issues=tuple(issues),
    )
    _write_sanitize_log(result, pdb_id=pdb_id, variant_label=variant_label)
    return result


@dataclass(frozen=True)
class AltlocSelectionSummary:
    """Summary of Bio.PDB disordered-atom selection."""

    selected_altloc_count: int
    dropped_altloc_atom_count: int


def select_altlocs_by_highest_occupancy(structure) -> AltlocSelectionSummary:  # noqa: ANN001
    """Select the best alternate location for every Bio.PDB disordered atom.

    The policy follows the same common rule used by pdb-tools: prefer the
    location with the highest occupancy, then blank altloc, then altloc ``A``.
    The structure is modified in place because Bio.PDB stores selected disordered
    atoms on the parsed object before PDBIO writes it.
    """

    selected_altloc_count = 0
    dropped_altloc_atom_count = 0

    for atom in _iter_residue_child_atoms(structure):
        child_dict = getattr(atom, "child_dict", None)
        disordered_select = getattr(atom, "disordered_select", None)

        if not child_dict or disordered_select is None:
            continue

        children = list(child_dict.values())
        if len(children) <= 1:
            continue

        selected_child = sorted(children, key=_altloc_sort_key)[0]
        disordered_select(selected_child.get_altloc())

        selected_altloc_count += 1
        dropped_altloc_atom_count += len(children) - 1

    return AltlocSelectionSummary(
        selected_altloc_count=selected_altloc_count,
        dropped_altloc_atom_count=dropped_altloc_atom_count,
    )


def count_nonunit_occupancies(structure) -> int:  # noqa: ANN001
    """Count selected atoms whose occupancy is absent or not equal to 1.00."""

    count = 0

    for atom in structure.get_atoms():
        occupancy = atom.get_occupancy()
        if occupancy is None or abs(float(occupancy) - 1.0) > 1e-6:
            count += 1

    return count


def find_nonstandard_residue_names(structure) -> list[str]:  # noqa: ANN001
    """Return amino-acid-like residue names that are not standard amino acids."""

    residue_names: set[str] = set()

    for residue in structure.get_residues():
        resname = residue.get_resname().strip().upper()
        if not _is_amino_acid_like(resname):
            continue
        if _is_standard_amino_acid_name(resname):
            continue
        residue_names.add(resname)

    return sorted(residue_names)


def count_written_residues(
    structure,  # noqa: ANN001
    *,
    drop_nonpolymer_heterogens: bool,
) -> int:
    """Count residues that will be written by the protein-only selector."""

    count = 0

    for residue in structure.get_residues():
        hetfield = residue.id[0]
        if not drop_nonpolymer_heterogens:
            count += 1
            continue
        if hetfield == " " or _is_amino_acid_like(residue.get_resname()):
            count += 1

    return count


def load_gromacs_heavy_atom_templates(
    *,
    force_field: str,
    gmx_executable: str | None = None,
) -> HeavyAtomTemplateLookup:
    """Load residue heavy atoms from the installed GROMACS RTP file."""

    rtp_path = find_gromacs_aminoacids_rtp_path(
        force_field=force_field,
        gmx_executable=gmx_executable,
    )

    if rtp_path is None:
        return HeavyAtomTemplateLookup(
            force_field=force_field,
            rtp_path=None,
            heavy_atoms_by_residue={},
        )

    return HeavyAtomTemplateLookup(
        force_field=force_field,
        rtp_path=rtp_path,
        heavy_atoms_by_residue=parse_gromacs_rtp_heavy_atoms(rtp_path),
    )


def find_missing_standard_heavy_atoms(
    structure,  # noqa: ANN001
    *,
    template_lookup: HeavyAtomTemplateLookup,
) -> list[str]:
    """Return missing-heavy-atom diagnostics for standard residues."""

    if not template_lookup.heavy_atoms_by_residue:
        return []

    messages: list[str] = []

    for residue in structure.get_residues():
        resname = residue.get_resname().strip().upper()
        if not _is_standard_amino_acid_name(resname):
            continue

        expected_atoms = _resolve_expected_heavy_atoms(
            resname=resname,
            template_lookup=template_lookup,
        )
        if not expected_atoms:
            continue

        observed_atoms = {
            atom.get_name().strip()
            for atom in residue.get_atoms()
            if not _is_hydrogen_atom(atom_name=atom.get_name(), element=atom.element)
        }

        missing_atoms = sorted(expected_atoms - observed_atoms)
        if not missing_atoms:
            continue

        chain = residue.get_parent()
        chain_id = chain.id if chain is not None else ""
        hetfield, resseq, icode = residue.id
        insertion_code = "" if icode == " " else str(icode)
        messages.append(
            f"Residue {resname} {chain_id or '<blank>'} {resseq}{insertion_code} "
            "is missing expected heavy atoms: "
            + ", ".join(missing_atoms)
        )

    return messages


def find_gromacs_aminoacids_rtp_path(
    *,
    force_field: str,
    gmx_executable: str | None = None,
) -> Path | None:
    """Find ``aminoacids.rtp`` for the requested installed GROMACS force field."""

    force_field_dir = f"{force_field}.ff"
    candidate_paths: list[Path] = []

    gmxlibrary = os.environ.get("GMXLIB")
    if gmxlibrary:
        candidate_paths.append(Path(gmxlibrary) / force_field_dir / "aminoacids.rtp")

    conda_prefix = os.environ.get("CONDA_PREFIX")
    if conda_prefix:
        candidate_paths.append(
            Path(conda_prefix)
            / "share"
            / "gromacs"
            / "top"
            / force_field_dir
            / "aminoacids.rtp"
        )

    candidate_paths.append(
        Path(sys.prefix)
        / "share"
        / "gromacs"
        / "top"
        / force_field_dir
        / "aminoacids.rtp"
    )

    executable = _resolve_gromacs_executable(gmx_executable)
    data_prefix = _read_gromacs_data_prefix(executable) if executable else None
    if data_prefix:
        candidate_paths.append(
            data_prefix / "share" / "gromacs" / "top" / force_field_dir / "aminoacids.rtp"
        )

    for candidate_path in candidate_paths:
        if candidate_path.exists():
            return candidate_path.resolve()

    return None


def parse_gromacs_rtp_heavy_atoms(rtp_path: str | Path) -> dict[str, frozenset[str]]:
    """Parse heavy atom names from GROMACS ``aminoacids.rtp``.

    The RTP parser is intentionally narrow. It only reads residue ``[ atoms ]``
    blocks and ignores hydrogens, terminal-only OXT, bonds, impropers, and
    charge-group details. That is sufficient for the sanitizer gate because we
    only need to know whether a standard residue lacks a heavy atom expected by
    the same force field that ``pdb2gmx`` will use.
    """

    path = Path(rtp_path)
    heavy_atoms_by_residue: dict[str, set[str]] = {}
    current_residue_name: str | None = None
    in_atoms_block = False

    for raw_line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        line = raw_line.split(";", maxsplit=1)[0].strip()
        if not line:
            continue

        if line.startswith("[") and line.endswith("]"):
            block_name = line.strip("[]").strip().upper()

            if block_name == "ATOMS":
                in_atoms_block = current_residue_name is not None
                continue

            current_residue_name = block_name
            heavy_atoms_by_residue.setdefault(current_residue_name, set())
            in_atoms_block = False
            continue

        if current_residue_name is None or not in_atoms_block:
            continue

        parts = line.split()
        if not parts:
            continue

        atom_name = parts[0].strip()
        if not atom_name or atom_name.startswith(("-", "+")):
            continue
        if atom_name == "OXT":
            continue
        if _is_hydrogen_atom(atom_name=atom_name, element=""):
            continue

        heavy_atoms_by_residue[current_residue_name].add(atom_name)

    return {
        residue_name: frozenset(heavy_atoms)
        for residue_name, heavy_atoms in heavy_atoms_by_residue.items()
        if heavy_atoms
    }


def _write_structure_with_biopython(
    structure,  # noqa: ANN001
    *,
    output_path: Path,
    selector: GromacsProteinSelect,
) -> None:
    """Write a structure with Bio.PDB while tolerating older PDBIO versions."""

    try:
        io = PDBIO(use_model_flag=0)
    except TypeError:
        io = PDBIO()

    io.set_structure(structure)
    io.save(str(output_path), select=selector, preserve_atom_numbering=False)


def _iter_residue_child_atoms(structure) -> Iterable[object]:  # noqa: ANN001
    """Yield direct residue atom children, including DisorderedAtom wrappers."""

    for residue in structure.get_residues():
        for atom in residue:
            yield atom


def _altloc_sort_key(atom) -> tuple[float, int, str]:  # noqa: ANN001
    """Return deterministic sort key for one alternate-position child atom."""

    occupancy = atom.get_occupancy()
    occupancy_value = float(occupancy) if occupancy is not None else 0.0
    altloc = (atom.get_altloc() or " ").strip().upper()

    if altloc == "":
        altloc_rank = 0
    elif altloc == "A":
        altloc_rank = 1
    else:
        altloc_rank = 2

    return (-occupancy_value, altloc_rank, altloc)


def _resolve_expected_heavy_atoms(
    *,
    resname: str,
    template_lookup: HeavyAtomTemplateLookup,
) -> frozenset[str]:
    """Return expected heavy atoms for one residue if available."""

    candidate_names = STANDARD_RESIDUE_TEMPLATE_ALIASES.get(resname, (resname,))

    for candidate_name in candidate_names:
        expected_atoms = template_lookup.heavy_atoms_by_residue.get(candidate_name)
        if expected_atoms:
            return expected_atoms

    return frozenset()


def _is_amino_acid_like(resname: str) -> bool:
    """Return whether Bio.PDB treats ``resname`` as amino-acid-like."""

    clean_resname = str(resname).strip().upper()

    try:
        return bool(is_aa(clean_resname, standard=False))
    except Exception:
        return _is_standard_amino_acid_name(clean_resname)


def _is_standard_amino_acid_name(resname: str) -> bool:
    """Return whether ``resname`` is one of the standard amino acids."""

    clean_resname = str(resname).strip().upper()

    try:
        return bool(is_aa(clean_resname, standard=True))
    except Exception:
        return clean_resname in {
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


def _is_hydrogen_atom(*, atom_name: str, element: str | None) -> bool:
    """Return whether an atom should be treated as hydrogen."""

    clean_atom_name = str(atom_name).strip().upper()
    clean_element = "" if element is None else str(element).strip().upper()

    if clean_element == "H":
        return True
    if clean_atom_name.startswith("H"):
        return True
    return len(clean_atom_name) > 1 and clean_atom_name[0].isdigit() and clean_atom_name[1] == "H"


def _resolve_gromacs_executable(gmx_executable: str | None) -> str | None:
    """Return a usable GROMACS executable path or command name."""

    if gmx_executable:
        return gmx_executable

    for executable_name in ("gmx", "gmx_mpi"):
        executable_path = shutil.which(executable_name)
        if executable_path:
            return executable_path

    return None


def _read_gromacs_data_prefix(executable: str | None) -> Path | None:
    """Read the GROMACS data prefix from ``gmx --version`` output."""

    if executable is None:
        return None

    try:
        completed_process = subprocess.run(
            [executable, "--version"],
            check=False,
            capture_output=True,
            text=True,
            timeout=30,
        )
    except Exception:
        return None

    combined_output = completed_process.stdout + "\n" + completed_process.stderr

    for line in combined_output.splitlines():
        if "Data prefix:" not in line:
            continue
        _label, value = line.split("Data prefix:", maxsplit=1)
        prefix = Path(value.strip())
        if prefix.exists():
            return prefix.resolve()

    return None


def _result_for_missing_input(
    *,
    input_path: Path,
    output_path: Path,
    log_path: Path | None,
) -> SanitizeResult:
    """Return a failed result for a missing input file."""

    issue = SanitizeIssue(
        severity="error",
        code="input_missing",
        message=f"Input PDB does not exist: {input_path}",
    )
    return SanitizeResult(
        input_pdb_path=input_path,
        output_pdb_path=output_path,
        log_path=log_path,
        output_written=False,
        can_run_gromacs=False,
        atom_count_in=0,
        atom_count_out=0,
        residue_count_out=0,
        selected_altloc_count=0,
        dropped_altloc_atom_count=0,
        normalized_occupancy_count=0,
        dropped_heterogen_atom_count=0,
        nonstandard_residue_names=(),
        missing_heavy_atom_count=0,
        force_field_template_path=None,
        issues=(issue,),
    )


def _result_for_parse_failure(
    *,
    input_path: Path,
    output_path: Path,
    log_path: Path | None,
    issues: list[SanitizeIssue],
) -> SanitizeResult:
    """Return a failed result for a PDB parse failure."""

    return SanitizeResult(
        input_pdb_path=input_path,
        output_pdb_path=output_path,
        log_path=log_path,
        output_written=False,
        can_run_gromacs=False,
        atom_count_in=0,
        atom_count_out=0,
        residue_count_out=0,
        selected_altloc_count=0,
        dropped_altloc_atom_count=0,
        normalized_occupancy_count=0,
        dropped_heterogen_atom_count=0,
        nonstandard_residue_names=(),
        missing_heavy_atom_count=0,
        force_field_template_path=None,
        issues=tuple(issues),
    )


def _write_sanitize_log(
    result: SanitizeResult,
    *,
    pdb_id: str | None,
    variant_label: str | None,
) -> None:
    """Write a plain-text sanitizer log."""

    if result.log_path is None:
        return

    result.log_path.parent.mkdir(parents=True, exist_ok=True)

    lines = [
        "sanitize_pdb_for_gromacs",
        "========================",
        "",
        f"pdb_id                       : {pdb_id or '<none>'}",
        f"variant_label                : {variant_label or '<none>'}",
        f"input_pdb_path               : {result.input_pdb_path}",
        f"output_pdb_path              : {result.output_pdb_path}",
        f"output_written               : {result.output_written}",
        f"can_run_gromacs              : {result.can_run_gromacs}",
        f"status                       : {result.status}",
        f"atom_count_in                : {result.atom_count_in}",
        f"atom_count_out               : {result.atom_count_out}",
        f"residue_count_out            : {result.residue_count_out}",
        f"selected_altloc_count        : {result.selected_altloc_count}",
        f"dropped_altloc_atom_count    : {result.dropped_altloc_atom_count}",
        f"normalized_occupancy_count   : {result.normalized_occupancy_count}",
        f"dropped_heterogen_atom_count : {result.dropped_heterogen_atom_count}",
        f"nonstandard_residue_names    : {', '.join(result.nonstandard_residue_names) or '<none>'}",
        f"missing_heavy_atom_count     : {result.missing_heavy_atom_count}",
        f"force_field_template_path    : {result.force_field_template_path or '<not found>'}",
        "",
        "issues",
        "------",
    ]

    if result.issues:
        for issue in result.issues:
            lines.append(f"[{issue.severity}] {issue.code}: {issue.message}")
    else:
        lines.append("<none>")

    result.log_path.write_text("\n".join(lines) + "\n", encoding="utf-8")