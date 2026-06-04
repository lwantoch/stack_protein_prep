"""Internal helpers for sanitize.py.

All private/helper functions, dataclasses, and helper classes are kept here so
sanitize.py can be a thin orchestrator.
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
    repaired_heavy_atom_count: int
    blocked_repair_count: int
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


# ---------------------------------------------------------------------------
# PDBFixer-based sidechain heavy-atom repair
# ---------------------------------------------------------------------------

_BACKBONE_HEAVY: frozenset[str] = frozenset({"N", "CA", "C", "O"})
# Minimum distance (Å) from a newly placed atom to any pre-existing heavy
# atom before a clash warning is raised.
_CLASH_WARN_DIST: float = 1.6


@dataclass(frozen=True)
class _RepairResult:
    n_repaired: int
    n_blocked: int
    repaired: tuple[str, ...]
    blocked: tuple[str, ...]
    error: str | None


def _repair_sidechain_heavy_atoms(
    pdb_path: Path,
    template_lookup: HeavyAtomTemplateLookup,
) -> _RepairResult:
    """Fill missing standard side-chain heavy atoms with PDBFixer (controlled).

    Rules:
    - ``missingResidues`` is emptied before ``addMissingAtoms`` — no loop building.
    - ``missingTerminals`` is emptied — no OXT / terminal patching.
    - Backbone atoms (N, CA, C, O) in ``missingAtoms`` are blocked, not repaired.
    - Non-standard residues are blocked.
    - ``addMissingHydrogens``, ``addSolvent``, ``replaceNonstandardResidues``,
      and ``removeHeterogens`` are never called.
    - After repair, original atom coordinates must be unchanged within 0.01 Å;
      any drift reverts the file to the pre-repair version.
    - Newly placed atoms within ``_CLASH_WARN_DIST`` of existing heavy atoms are
      flagged as warnings but not blocked (PDBFixer's template placement is
      usually sane; a downstream clash check is the right gate).
    """

    try:
        from pdbfixer import PDBFixer
        from openmm.app import PDBFile as OpenMMPDBFile
    except ImportError:
        return _RepairResult(
            n_repaired=0, n_blocked=0,
            repaired=(), blocked=(),
            error="PDBFixer / OpenMM not available; skipping repair.",
        )

    import math
    import tempfile

    # ---- run PDBFixer in controlled mode ------------------------------------
    fixer = PDBFixer(filename=str(pdb_path))
    fixer.findMissingResidues()
    fixer.missingResidues = {}          # no loop / residue insertion
    fixer.findMissingAtoms()

    allowed: dict = {}
    blocked_msgs: list[str] = []
    repaired_msgs: list[str] = []

    for residue, atoms in list((fixer.missingAtoms or {}).items()):
        resname = residue.name.strip().upper()
        atom_names = [a.name for a in atoms]
        label = f"{resname} {residue.chain.id} {residue.id}"

        if not _is_standard_amino_acid_name(resname):
            blocked_msgs.append(f"{label}: {atom_names} — nonstandard residue")
            continue

        backbone_hits = [a for a in atom_names if a in _BACKBONE_HEAVY]
        if backbone_hits:
            blocked_msgs.append(f"{label}: {backbone_hits} — backbone heavy atom")
            continue

        allowed[residue] = atoms
        repaired_msgs.append(f"{label}: {atom_names}")

    if not allowed:
        return _RepairResult(
            n_repaired=0, n_blocked=len(blocked_msgs),
            repaired=(), blocked=tuple(blocked_msgs),
            error=None,
        )

    fixer.missingAtoms = allowed
    fixer.missingTerminals = {}
    fixer.addMissingAtoms()

    # ---- write to a temp file and validate ----------------------------------
    tmp_fd, tmp_path_str = tempfile.mkstemp(suffix=".pdb", prefix="fruton_repair_")
    tmp_path = Path(tmp_path_str)
    try:
        import os
        os.close(tmp_fd)
        with tmp_path.open("w") as fh:
            OpenMMPDBFile.writeFile(fixer.topology, fixer.positions, fh)

        # Strict diff: original atom coordinates must be preserved.
        parser = PDBParser(QUIET=True, PERMISSIVE=True)
        orig = {
            (ch.id, res.id, atom.get_name().strip()): atom.get_vector().get_array()
            for ch in parser.get_structure("o", str(pdb_path)).get_chains()
            for res in ch.get_residues()
            for atom in res.get_atoms()
            if not _is_hydrogen_atom(atom_name=atom.get_name(), element=atom.element)
        }
        repaired_struct = parser.get_structure("r", str(tmp_path))
        repaired_atoms = {
            (ch.id, res.id, atom.get_name().strip()): atom.get_vector().get_array()
            for ch in repaired_struct.get_chains()
            for res in ch.get_residues()
            for atom in res.get_atoms()
            if not _is_hydrogen_atom(atom_name=atom.get_name(), element=atom.element)
        }

        drifted = []
        for key, orig_xyz in orig.items():
            rep_xyz = repaired_atoms.get(key)
            if rep_xyz is None:
                drifted.append(f"{key} disappeared")
                continue
            dist = math.sqrt(sum((a - b) ** 2 for a, b in zip(orig_xyz, rep_xyz)))
            if dist > 0.01:
                drifted.append(f"{key} moved {dist:.3f} Å")

        if drifted:
            return _RepairResult(
                n_repaired=0, n_blocked=len(blocked_msgs),
                repaired=(), blocked=tuple(blocked_msgs),
                error="Original atom drift detected after PDBFixer repair — file not replaced: "
                      + "; ".join(drifted[:5]),
            )

        # Clash check on newly added atoms (warn only).
        new_keys = set(repaired_atoms) - set(orig)
        orig_xyz_list = list(orig.values())
        clash_warnings: list[str] = []
        for key in new_keys:
            new_xyz = repaired_atoms[key]
            min_dist = min(
                math.sqrt(sum((a - b) ** 2 for a, b in zip(new_xyz, o)))
                for o in orig_xyz_list
            ) if orig_xyz_list else float("inf")
            if min_dist < _CLASH_WARN_DIST:
                clash_warnings.append(
                    f"{key[2]} {key[0]}:{key[1]} — nearest existing atom {min_dist:.2f} Å"
                )

        if clash_warnings:
            blocked_msgs.extend(
                f"clash-warn: {w}" for w in clash_warnings
            )

        # Replace the original output with the repaired version.
        import shutil
        shutil.move(str(tmp_path), str(pdb_path))

    except Exception as exc:
        return _RepairResult(
            n_repaired=0, n_blocked=len(blocked_msgs),
            repaired=(), blocked=tuple(blocked_msgs),
            error=f"{type(exc).__name__}: {exc}",
        )
    finally:
        if tmp_path.exists():
            tmp_path.unlink(missing_ok=True)

    n_repaired = sum(len(atoms) for atoms in allowed.values())
    return _RepairResult(
        n_repaired=n_repaired, n_blocked=len(blocked_msgs),
        repaired=tuple(repaired_msgs), blocked=tuple(blocked_msgs),
        error=None,
    )


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
        repaired_heavy_atom_count=0,
        blocked_repair_count=0,
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
        repaired_heavy_atom_count=0,
        blocked_repair_count=0,
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
        f"repaired_heavy_atom_count    : {result.repaired_heavy_atom_count}",
        f"blocked_repair_count         : {result.blocked_repair_count}",
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
