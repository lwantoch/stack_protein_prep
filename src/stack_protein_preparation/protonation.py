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
import os
import shlex
import shutil
import subprocess
import traceback
from datetime import datetime
from pathlib import Path
from typing import Literal

from stack_protein_preparation.pdb_components import WATER_NAMES

ProtonationInputSource = Literal["protein", "filler", "modeller", "alphafold"]

GromacsWaterModel = Literal[
    "none",
    "opc",
    "opc3",
    "spc",
    "spce",
    "tip3p",
    "tip4p",
    "tip4pew",
    "tip5p",
    "tips3p",
]

GromacsChainSeparation = Literal[
    "id_or_ter",
    "id_and_ter",
    "ter",
    "id",
]

GromacsMergeMode = Literal[
    "no",
    "all",
]

DEFAULT_FILLER_FINAL_FILENAME = "final_filled_model.pdb"
DEFAULT_GROMACS_FORCE_FIELD = "amber99sb-ildn"
DEFAULT_GROMACS_WATER_MODEL: GromacsWaterModel = "tip3p"
DEFAULT_GROMACS_CHAIN_SEPARATION: GromacsChainSeparation = "id_or_ter"
DEFAULT_GROMACS_MERGE_MODE: GromacsMergeMode = "no"
DEFAULT_GROMACS_TIMEOUT_SECONDS = 600
MODULE_NAME = "protonation"


# ============================================================================
# terminal helpers
# ============================================================================


def _screen_header(pdb_id: str, message: str) -> None:
    print(f"┏━━ [{MODULE_NAME}] {pdb_id}")
    print(f"┗━ {message}")


def _screen_item(label: str, value: str) -> None:
    print(f"   ├─ {label}: {value}")


def _screen_result(status: str, message: str) -> None:
    print(f"   └─ result: {status} | {message}")


# ============================================================================
# logging helpers
# ============================================================================


def _timestamp() -> str:
    return datetime.now().strftime("%Y-%m-%d %H:%M:%S")


def _infer_pdb_id_from_path(path: str | Path) -> str:
    path = Path(path)

    for parent in path.parents:
        if parent.name == "components" and parent.parent.name:
            return parent.parent.name.upper()

    stem = path.stem.strip()
    if stem:
        return stem.split("_")[0].upper()

    return "UNKNOWN"


def _infer_protein_dir_from_path(path: str | Path) -> Path:
    path = Path(path)

    if path.parent.name == "components":
        return path.parent.parent

    for parent in path.parents:
        if parent.name == "components":
            return parent.parent

    return path.parent


def _get_module_log_path(
    *,
    pdb_id: str | None = None,
    protein_dir: str | Path | None = None,
    probe_path: str | Path | None = None,
) -> Path:
    """Return the standard FRUTON log path for this module.

    The expected convention is ``data/proteins/logs/<PDB_ID>/protonation.log``.
    ``protein_dir`` should normally be ``data/proteins/<PDB_ID>``. When no
    explicit protein directory is provided, the helper infers it from a path
    below a ``components`` directory. This keeps logs for original structures and
    variant structures in the same per-protein logging location.
    """

    if protein_dir is None:
        if probe_path is None:
            raise ValueError("Need protein_dir or probe_path to infer log path.")
        protein_dir = _infer_protein_dir_from_path(probe_path)

    protein_dir = Path(protein_dir)

    if pdb_id is None:
        if probe_path is not None:
            pdb_id = _infer_pdb_id_from_path(probe_path)
        else:
            pdb_id = protein_dir.name.upper()

    logs_dir = protein_dir.parent / "logs" / str(pdb_id).upper()
    logs_dir.mkdir(parents=True, exist_ok=True)
    return logs_dir / f"{MODULE_NAME}.log"


def _append_log_block(
    log_path: Path,
    title: str,
    lines: list[str],
) -> None:
    log_path.parent.mkdir(parents=True, exist_ok=True)

    with log_path.open("a", encoding="utf-8") as handle:
        handle.write("═" * 100 + "\n")
        handle.write(f"[{_timestamp()}] {title}\n")
        handle.write("─" * 100 + "\n")
        for line in lines:
            handle.write(f"{line}\n")
        handle.write("\n")


def _append_text_block(
    log_path: Path,
    title: str,
    text: str,
) -> None:
    safe_text = text if text.strip() else "<empty>"

    log_path.parent.mkdir(parents=True, exist_ok=True)

    with log_path.open("a", encoding="utf-8") as handle:
        handle.write("═" * 100 + "\n")
        handle.write(f"[{_timestamp()}] {title}\n")
        handle.write("─" * 100 + "\n")
        handle.write(f"{safe_text}\n\n")


def _log_exception(
    log_path: Path,
    title: str,
    exc: Exception,
) -> None:
    _append_log_block(
        log_path,
        title,
        [
            f"exception_type: {type(exc).__name__}",
            f"exception_repr: {exc!r}",
            "traceback:",
            traceback.format_exc().rstrip(),
        ],
    )


# ============================================================================
# path / naming helpers
# ============================================================================


def _find_gmx_executable() -> str:
    """Return the first available GROMACS executable.

    Local workstation installations usually expose ``gmx``. Some MPI-enabled or
    cluster installations expose ``gmx_mpi`` instead. The function returns the
    executable name found in the active ``PATH`` so the subprocess uses the same
    environment that launched the FRUTON pipeline. A missing executable is an
    environment error and is therefore raised before any protonation attempt is
    treated as a chemistry failure.
    """

    candidates = ["gmx", "gmx_mpi"]

    for candidate in candidates:
        if shutil.which(candidate):
            return candidate

    raise FileNotFoundError(
        "Could not find a GROMACS executable in PATH. "
        "Expected one of: gmx, gmx_mpi."
    )


def count_atoms_in_structure_file(pdb_path: str | Path) -> int:
    """Count coordinate records in a PDB-like text structure file.

    This intentionally counts ``ATOM`` and ``HETATM`` records directly instead
    of using a structural parser. The helper is used only for cheap diagnostics
    before and after GROMACS execution. It does not validate chemistry,
    connectivity, residue states, or force-field compatibility.
    """

    pdb_path = Path(pdb_path)

    if not pdb_path.is_file():
        raise FileNotFoundError(f"Structure file not found: {pdb_path}")

    atom_count = 0

    with pdb_path.open("r", encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("ATOM") or line.startswith("HETATM"):
                atom_count += 1

    return atom_count


def count_atoms_in_pdb(pdb_path: str | Path) -> int:
    """Backward-compatible wrapper around ``count_atoms_in_structure_file``."""

    return count_atoms_in_structure_file(pdb_path)


def _is_valid_structure_file(path: str | Path | None) -> bool:
    """Return ``True`` when a candidate structure path exists and is non-empty."""

    if path is None:
        return False

    path = Path(path)
    return path.is_file() and path.stat().st_size > 0


def sanitize_variant_label(variant_label: str) -> str:
    """Convert a variant label into a filesystem-safe token.

    Variant labels are used in output filenames and log filenames. This helper
    keeps letters, digits, underscores, and hyphens unchanged while replacing
    everything else by underscores. Empty labels are rejected so accidental blank
    variant names do not overwrite standard protonation outputs.
    """

    cleaned = "".join(
        character if character.isalnum() or character in {"_", "-"} else "_"
        for character in variant_label.strip()
    )
    cleaned = cleaned.strip("_")

    if not cleaned:
        raise ValueError(f"Invalid empty variant label derived from {variant_label!r}")

    return cleaned


def build_standard_protonation_output_path(
    pdb_id: str,
    protein_dir: str | Path,
) -> Path:
    """Return the standard single-path protonated coordinate output.

    Example
    -------
    ``data/proteins/1ABC/components/1ABC_proteinH.pdb``
    """

    protein_dir = Path(protein_dir)
    components_dir = protein_dir / "components"
    return components_dir / f"{pdb_id}_proteinH.pdb"


def build_variant_protonation_output_path(
    pdb_id: str,
    protein_dir: str | Path,
    variant_label: str,
) -> Path:
    """Return the variant-specific protonated coordinate output.

    Example
    -------
    ``components/1ABC_proteinH_gaps.pdb``
    ``components/1ABC_proteinH_best_complete.pdb``
    """

    protein_dir = Path(protein_dir)
    components_dir = protein_dir / "components"
    safe_variant_label = sanitize_variant_label(variant_label)
    return components_dir / f"{pdb_id}_proteinH_{safe_variant_label}.pdb"


def build_gromacs_topology_output_path(output_pdb: str | Path) -> Path:
    """Return the GROMACS topology path paired with a protonated coordinate file.

    ``pdb2gmx`` always generates a topology in addition to coordinates. Keeping
    the topology beside the protonated PDB makes the artefacts inspectable and
    prevents topology includes from being scattered across unrelated runtime
    directories. The topology is not used as the authoritative downstream
    topology unless later FRUTON stages explicitly choose to consume it.
    """

    output_pdb = Path(output_pdb)
    return output_pdb.with_suffix(".top")


def build_gromacs_position_restraints_output_path(output_pdb: str | Path) -> Path:
    """Return the position-restraints include path paired with the output PDB."""

    output_pdb = Path(output_pdb)
    return output_pdb.with_name(f"{output_pdb.stem}_posre.itp")


def build_protonation_stdout_log_path(
    pdb_id: str,
    protein_dir: str | Path,
    variant_label: str | None = None,
) -> Path:
    protein_dir = Path(protein_dir)
    components_dir = protein_dir / "components"

    if variant_label is None:
        return components_dir / f"{pdb_id}_protonation_stdout.log"

    safe_variant_label = sanitize_variant_label(variant_label)
    return components_dir / f"{pdb_id}_protonation_stdout_{safe_variant_label}.log"


def build_protonation_stderr_log_path(
    pdb_id: str,
    protein_dir: str | Path,
    variant_label: str | None = None,
) -> Path:
    protein_dir = Path(protein_dir)
    components_dir = protein_dir / "components"

    if variant_label is None:
        return components_dir / f"{pdb_id}_protonation_stderr.log"

    safe_variant_label = sanitize_variant_label(variant_label)
    return components_dir / f"{pdb_id}_protonation_stderr_{safe_variant_label}.log"


def find_default_filler_final_model_path(
    filler_output_dir: str | Path,
) -> Path:
    """Return the default unified filler final model path.

    Expected filler convention
    --------------------------
    ``<filler_output_dir>/final_filled_model.pdb``
    """

    filler_output_dir = Path(filler_output_dir)
    return filler_output_dir / DEFAULT_FILLER_FINAL_FILENAME


# ============================================================================
# selection helpers
# ============================================================================


def select_protonation_input(
    pdb_id: str,
    protein_dir: str | Path,
    *,
    filler_final_model_path: str | Path | None = None,
    filler_output_dir: str | Path | None = None,
    modeller_model_path: str | Path | None = None,
    alphafold_model_path: str | Path | None = None,
) -> tuple[Path, ProtonationInputSource]:
    """Select the input structure for protonation.

    Preference order mirrors the previous FRUTON protonation module. Explicit
    filler output is preferred, then the default filler final model, then
    MODELLER, then AlphaFold, and finally the original component protein. The
    function checks only whether candidate files exist and are non-empty. It
    does not inspect whether ``pdb2gmx`` will accept the residues, termini, or
    chain layout.
    """

    protein_dir = Path(protein_dir)
    components_dir = protein_dir / "components"
    default_protein_path = components_dir / f"{pdb_id}_protein.pdb"

    if _is_valid_structure_file(filler_final_model_path):
        return Path(filler_final_model_path), "filler"

    if filler_output_dir is not None:
        default_filler_final_model_path = find_default_filler_final_model_path(
            filler_output_dir=filler_output_dir
        )
        if _is_valid_structure_file(default_filler_final_model_path):
            return default_filler_final_model_path, "filler"

    if _is_valid_structure_file(modeller_model_path):
        return Path(modeller_model_path), "modeller"

    if _is_valid_structure_file(alphafold_model_path):
        return Path(alphafold_model_path), "alphafold"

    if _is_valid_structure_file(default_protein_path):
        return default_protein_path, "protein"

    raise FileNotFoundError(
        f"No valid protonation input found for {pdb_id}. "
        f"Checked filler_final_model_path={filler_final_model_path!s}, "
        f"filler_output_dir={filler_output_dir!s}, "
        f"modeller_model_path={modeller_model_path!s}, "
        f"alphafold_model_path={alphafold_model_path!s}, and "
        f"default protein path={default_protein_path!s}."
    )


def _write_text_log(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")


def _preview_text(text: str, max_chars: int = 1200) -> str:
    stripped = text.strip()
    if not stripped:
        return ""
    if len(stripped) <= max_chars:
        return stripped
    return stripped[:max_chars] + "\n...[truncated]..."


def _basic_structure_diagnostics(input_pdb: Path) -> dict[str, int | bool]:
    """Return cheap pre-run diagnostics for GROMACS failures.

    These counters help distinguish empty/malformed inputs from force-field
    interpretation errors. They are intentionally text-record diagnostics rather
    than structural validation. The external ``pdb2gmx`` process remains the
    actual source of truth for whether the input can be processed.
    """

    atom_count = 0
    hetatm_count = 0
    ter_count = 0
    end_count = 0

    with input_pdb.open("r", encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("ATOM"):
                atom_count += 1
            elif line.startswith("HETATM"):
                hetatm_count += 1
            elif line.startswith("TER"):
                ter_count += 1
            elif line.startswith("END"):
                end_count += 1

    return {
        "atom_count": atom_count,
        "hetatm_count": hetatm_count,
        "ter_count": ter_count,
        "end_count": end_count,
        "has_any_atoms": (atom_count + hetatm_count) > 0,
    }


# ============================================================================
# PDB element-column repair
# ============================================================================

_TWO_CHAR_ELEMENTS = frozenset({
    "CL", "BR", "FE", "ZN", "MG", "CA", "MN", "CU", "CO", "NI",
    "NA", "HG", "SE", "SI",
})


def _infer_element_from_atom_name(atom_name_field: str) -> str:
    """Return the element symbol inferred from a 4-char PDB atom name field.

    Hydrogen names (H, HA, HB, HG1, HW1, 1HB …) always resolve to 'H'.
    Two-character heavy elements (CL, FE, ZN …) are only matched when the
    atom name core does not begin with 'H', so 'HG1' (gamma-H) is never
    mis-classified as mercury (HG).
    """
    name = atom_name_field.strip()
    if not name:
        return ""
    i = 0
    while i < len(name) and name[i].isdigit():
        i += 1
    core = name[i:] if i < len(name) else name
    if not core:
        return ""
    if core[0].upper() == "H":
        return "H"
    two = core[:2].upper()
    if two in _TWO_CHAR_ELEMENTS:
        return two
    return core[0].upper()


def _fill_pdb_element_column(pdb_path: Path) -> None:
    """Fill the element symbol column (cols 77-78) for ATOM/HETATM records.

    gmx pdb2gmx writes hydrogen records without the element column, leaving
    lines at 66 characters.  This pass fills the column from the atom name so
    that downstream tools that rely on the element field (MDAnalysis, VMD,
    OpenMM) get correct element assignments.
    """
    raw = pdb_path.read_text(encoding="utf-8")
    lines = raw.splitlines(keepends=True)
    fixed: list[str] = []
    changed = False
    for line in lines:
        if line.startswith(("ATOM  ", "HETATM")):
            body = line.rstrip("\n\r")
            elem_field = body[76:78] if len(body) >= 78 else ""
            elem_str = elem_field.strip()
            elem_missing = not elem_str
            # A 2-char element whose first character is 'H' (e.g. "HG") was
            # mis-assigned by a previous pass that lacked the H-early-return
            # rule.  Re-infer those as well.
            has_wrong_h = (
                len(elem_str) == 2 and elem_str[0].upper() == "H"
            )
            # Lines longer than 78 chars have junk beyond the element column
            # (never written by GROMACS) that must be trimmed.
            has_junk = len(body) > 78
            if elem_missing or has_wrong_h or has_junk:
                atom_name_field = body[12:16] if len(body) >= 16 else ""
                elem = _infer_element_from_atom_name(atom_name_field)
                if not elem and not elem_missing:
                    elem = elem_str  # keep existing if inference yields nothing
                if elem:
                    new_body = body[:76].ljust(76) + elem.rjust(2)
                    if new_body != body:
                        body = new_body
                        changed = True
            eol = line[len(line.rstrip("\n\r")):]
            fixed.append(body + eol)
        else:
            fixed.append(line)
    if changed:
        pdb_path.write_text("".join(fixed), encoding="utf-8")


# ============================================================================
# crystal water preparation
# ============================================================================

def _rename_water_pdb_to_sol(input_path: Path, output_path: Path) -> None:
    """Rewrite a water PDB renaming residues to SOL and O atoms to OW.

    This is a plain-text preparation step so that gmx pdb2gmx recognises the
    water residues and applies the TIP3P template when adding hydrogens.
    """
    result_lines: list[str] = []
    with input_path.open("r", encoding="utf-8") as fh:
        for line in fh:
            if line.startswith("ATOM") or line.startswith("HETATM"):
                if line[17:20].strip() in WATER_NAMES:
                    atom_name = line[12:16].strip()
                    padded = " OW " if atom_name in ("O", "OW") else f" {atom_name:<3}"
                    line = line[:12] + padded + line[16:17] + "SOL" + line[20:]
            result_lines.append(line)
    output_path.write_text("".join(result_lines), encoding="utf-8")


def add_hydrogens_to_crystal_water_pdb(
    water_pdb_path: str | Path,
    ff: str = DEFAULT_GROMACS_FORCE_FIELD,
    water_model: GromacsWaterModel = DEFAULT_GROMACS_WATER_MODEL,
) -> None:
    """Add TIP3P hydrogens to crystal waters using gmx pdb2gmx.

    Renames water residues to SOL and O atoms to OW so that pdb2gmx recognises
    them, then runs ``gmx pdb2gmx -ignh`` with the TIP3P template to place HW1
    and HW2. The input file is overwritten in-place with the result. Empty or
    missing files are silently skipped.
    """
    import shutil
    import tempfile

    water_pdb_path = Path(water_pdb_path)
    if not water_pdb_path.exists() or water_pdb_path.stat().st_size == 0:
        return

    with tempfile.TemporaryDirectory() as tmp_dir:
        tmp = Path(tmp_dir)
        sol_input = tmp / "water_sol.pdb"
        sol_output = tmp / "water_solH.pdb"

        _rename_water_pdb_to_sol(water_pdb_path, sol_input)

        run_gmx_pdb2gmx_protonation(
            input_pdb=sol_input,
            output_pdb=sol_output,
            topology_output_path=tmp / "water_sol.top",
            position_restraints_output_path=tmp / "water_sol_posre.itp",
            ff=ff,
            water_model=water_model,
            ignore_input_hydrogens=True,
        )

        if not sol_output.exists() or sol_output.stat().st_size == 0:
            raise RuntimeError(
                f"gmx pdb2gmx produced no output for water component: {water_pdb_path}"
            )

        shutil.copy2(sol_output, water_pdb_path)
        _fill_pdb_element_column(water_pdb_path)


# ============================================================================
# PROPKA protonation state prediction
# ============================================================================


def _run_propka(pdb_path: Path) -> object | None:
    """Run PROPKA on pdb_path; return the MolecularContainer or None on error."""
    import contextlib
    import io

    try:
        import propka.run  # type: ignore[import]
    except ImportError:
        return None

    try:
        buf = io.StringIO()
        with contextlib.redirect_stdout(buf), contextlib.redirect_stderr(buf):
            mol = propka.run.single(str(pdb_path), optargs=["--quiet"])
        return mol
    except Exception:
        return None


def _get_propka_his_assignments(
    mol: object,
    ph: float,
) -> list[dict[str, str | int | float]]:
    """Return one record per titratable HIS group with pKa and assigned state.

    Each record has keys: chain, res_num, icode, pka, assigned.
    assigned is 'HIP' when pka > ph, else 'HIE'.
    Records are sorted by (chain, res_num, icode).
    """
    assignments: list[dict[str, str | int | float]] = []
    try:
        groups = mol.conformations["AVR"].groups  # type: ignore[union-attr]
    except (KeyError, AttributeError):
        return assignments
    for g in groups:
        if g.residue_type != "HIS" or not g.titratable:
            continue
        icode = (g.atom.icode or " ").strip()
        assignments.append({
            "chain": g.atom.chain_id,
            "res_num": int(g.atom.res_num),
            "icode": icode,
            "pka": round(float(g.pka_value), 2),
            "assigned": "HIP" if g.pka_value > ph else "HIE",
        })
    return sorted(assignments, key=lambda d: (str(d["chain"]), int(d["res_num"]), str(d["icode"])))


def predict_protonation_states(
    pdb_path: Path | str,
    ph: float = 7.4,
) -> dict[tuple[str, int, str], str]:
    """Return PROPKA-based HIP rename keys for residues with pKa > ph.

    Maps (chain_id, res_num, icode) → 'HIP'. Only HIP entries are returned;
    HIE is the GROMACS/AMBER default and needs no rename. Returns an empty
    dict when PROPKA is unavailable or the run fails.
    """
    pdb_path = Path(pdb_path)
    mol = _run_propka(pdb_path)
    if mol is None:
        return {}
    return {
        (str(d["chain"]), int(d["res_num"]), str(d["icode"])): d["assigned"]
        for d in _get_propka_his_assignments(mol, ph)
        if d["assigned"] == "HIP"
    }


def _apply_protonation_renames(
    input_path: Path,
    output_path: Path,
    renames: dict[tuple[str, int, str], str],
) -> None:
    """Rewrite input_path with per-residue name substitutions and copy to output_path.

    Only ATOM/HETATM records are touched; residue name is at PDB columns 17-19.
    """
    if not renames:
        import shutil

        shutil.copy2(input_path, output_path)
        return

    lines: list[str] = []
    with input_path.open("r", encoding="utf-8") as fh:
        for line in fh:
            if line.startswith(("ATOM  ", "HETATM")):
                chain = line[21] if len(line) > 21 else " "
                try:
                    res_num = int(line[22:26])
                except ValueError:
                    lines.append(line)
                    continue
                icode = line[26].strip() if len(line) > 26 else ""
                key = (chain, res_num, icode)
                if key in renames:
                    new_name = f"{renames[key]:<3}"
                    line = line[:17] + new_name + line[20:]
            lines.append(line)
    output_path.write_text("".join(lines), encoding="utf-8")


# ============================================================================
# external tool execution
# ============================================================================


def run_gmx_pdb2gmx_protonation(
    input_pdb: str | Path,
    output_pdb: str | Path,
    *,
    topology_output_path: str | Path | None = None,
    position_restraints_output_path: str | Path | None = None,
    ff: str = DEFAULT_GROMACS_FORCE_FIELD,
    water_model: GromacsWaterModel = DEFAULT_GROMACS_WATER_MODEL,
    chain_separation: GromacsChainSeparation = DEFAULT_GROMACS_CHAIN_SEPARATION,
    merge: GromacsMergeMode = DEFAULT_GROMACS_MERGE_MODE,
    ignore_input_hydrogens: bool = True,
    timeout_seconds: int = DEFAULT_GROMACS_TIMEOUT_SECONDS,
) -> subprocess.CompletedProcess[str]:
    """Run ``gmx pdb2gmx`` and write protonated coordinates.

    The command is executed in the output directory. The input path is passed as
    an absolute path, while output coordinate, topology, and position-restraint
    filenames are passed relative to the working directory. This keeps generated
    GROMACS artefacts beside the FRUTON protonated output and avoids topology
    includes that refer to unrelated temporary paths.

    The subprocess receives empty stdin so FRUTON does not hang if GROMACS asks
    an unexpected interactive question. Force field and water model are provided
    explicitly to avoid normal force-field selection prompts. The command uses
    ``-ignh`` by default because the requested route is to rebuild hydrogens via
    GROMACS rather than preserving hydrogens already present in the input file.
    """

    input_pdb = Path(input_pdb).resolve()
    output_pdb = Path(output_pdb)
    topology_output_path = (
        Path(topology_output_path)
        if topology_output_path is not None
        else build_gromacs_topology_output_path(output_pdb)
    )
    position_restraints_output_path = (
        Path(position_restraints_output_path)
        if position_restraints_output_path is not None
        else build_gromacs_position_restraints_output_path(output_pdb)
    )

    if not input_pdb.is_file():
        raise FileNotFoundError(f"Protonation input structure not found: {input_pdb}")

    output_pdb.parent.mkdir(parents=True, exist_ok=True)
    topology_output_path.parent.mkdir(parents=True, exist_ok=True)
    position_restraints_output_path.parent.mkdir(parents=True, exist_ok=True)

    if topology_output_path.parent.resolve() != output_pdb.parent.resolve():
        raise ValueError(
            "Topology output must be written beside the protonated coordinate "
            f"file. Got topology={topology_output_path}, output={output_pdb}."
        )

    if position_restraints_output_path.parent.resolve() != output_pdb.parent.resolve():
        raise ValueError(
            "Position-restraints output must be written beside the protonated "
            f"coordinate file. Got posre={position_restraints_output_path}, "
            f"output={output_pdb}."
        )

    executable = _find_gmx_executable()

    cmd = [
        executable,
        "pdb2gmx",
        "-f",
        str(input_pdb),
        "-o",
        output_pdb.name,
        "-p",
        topology_output_path.name,
        "-i",
        position_restraints_output_path.name,
        "-ff",
        ff,
        "-water",
        water_model,
        "-chainsep",
        chain_separation,
        "-merge",
        merge,
    ]

    if ignore_input_hydrogens:
        cmd.append("-ignh")

    env = os.environ.copy()
    env.setdefault("GMX_MAXBACKUP", "-1")

    return subprocess.run(
        cmd,
        check=False,
        capture_output=True,
        text=True,
        input="",
        cwd=str(output_pdb.parent),
        timeout=timeout_seconds,
        env=env,
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
    """Protonate one explicitly selected structure with GROMACS.

    This is the variant-friendly core function. It does not decide which input
    should be used; it only protonates the path it receives. The return
    dictionary keeps the old FRUTON style and adds GROMACS-specific artefact
    paths. ``ph`` is kept to preserve the previous API and to record the intended
    historical protonation setting, but it is not applied by ``pdb2gmx``.
    """

    input_pdb = Path(input_pdb)
    output_pdb = Path(output_pdb)
    topology_output_path = build_gromacs_topology_output_path(output_pdb)
    position_restraints_output_path = build_gromacs_position_restraints_output_path(
        output_pdb
    )

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
        executable,
        "pdb2gmx",
        "-f",
        str(input_pdb.resolve()),
        "-o",
        output_pdb.name,
        "-p",
        topology_output_path.name,
        "-i",
        position_restraints_output_path.name,
        "-ff",
        ff,
        "-water",
        water_model,
        "-chainsep",
        chain_separation,
        "-merge",
        merge,
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

    import tempfile

    propka_mol = _run_propka(input_pdb)
    propka_ran = propka_mol is not None
    his_assignments = _get_propka_his_assignments(propka_mol, ph) if propka_ran else []
    propka_renames: dict[tuple[str, int, str], str] = {
        (str(d["chain"]), int(d["res_num"]), str(d["icode"])): str(d["assigned"])
        for d in his_assignments
        if d["assigned"] == "HIP"
    }
    propka_him_count = len(propka_renames)
    propka_applied = bool(propka_renames)

    _screen_item("propka_ran", str(propka_ran))
    _screen_item("propka_his_total", str(len(his_assignments)))
    _screen_item("propka_hip_count", str(propka_him_count))
    _append_log_block(
        module_log_path,
        "protonate_selected_structure:propka",
        [
            f"propka_ran              : {propka_ran}",
            f"propka_ph               : {ph}",
            f"propka_his_total        : {len(his_assignments)}",
            f"propka_hip_count        : {propka_him_count}",
            f"propka_applied          : {propka_applied}",
            f"propka_assignments      : {his_assignments}",
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
        output_atom_count = count_atoms_in_structure_file(output_pdb)
        atom_count_increased = output_atom_count > input_atom_count
    else:
        output_atom_count = 0
        atom_count_increased = False

    protonation_success = result.returncode == 0 and output_nonempty

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
    _append_text_block(
        module_log_path,
        "protonate_selected_structure:stdout",
        result.stdout,
    )
    _append_text_block(
        module_log_path,
        "protonate_selected_structure:stderr",
        result.stderr,
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
        "protonation_method": "PROPKA pKa prediction + gmx pdb2gmx -ignh",
        "protonation_input_source": input_source,
        "protonation_input_path": str(input_pdb),
        "protonation_output_path": str(output_pdb),
        "protonation_topology_path": str(topology_output_path),
        "protonation_position_restraints_path": str(position_restraints_output_path),
        "protonation_ph": ph,
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
    """Protonate one named FRUTON variant with GROMACS.

    The output path follows the previous variant-specific convention, e.g.
    ``components/<PDBID>_proteinH_gaps.pdb`` or
    ``components/<PDBID>_proteinH_best_complete.pdb``. The paired GROMACS
    topology and position-restraints files are written beside that coordinate
    output. ``ph`` is accepted for compatibility with the old API but is only
    recorded, not applied by GROMACS.
    """

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
    """Backward-compatible single-path FRUTON protonation workflow.

    Selection priority is unchanged: filler final model, MODELLER model,
    AlphaFold model, and original component protein. The coordinate output path
    remains ``components/<PDBID>_proteinH.pdb`` so downstream FRUTON stages do
    not need to change their expected filename. The paired GROMACS topology and
    position-restraints files are written beside the protonated coordinate file.
    """

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