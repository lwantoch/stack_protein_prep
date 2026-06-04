"""Core helpers for the protonation pipeline — shared by protonation.py."""
from __future__ import annotations

import contextlib
import io
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
from stack_protein_preparation._protonation_paths import (
    build_gromacs_position_restraints_output_path,
    build_gromacs_topology_output_path,
    build_protonation_stderr_log_path,
    build_protonation_stdout_log_path,
    find_default_filler_final_model_path,
)

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


def _append_log_block(log_path: Path, title: str, lines: list[str]) -> None:
    log_path.parent.mkdir(parents=True, exist_ok=True)
    with log_path.open("a", encoding="utf-8") as handle:
        handle.write("═" * 100 + "\n")
        handle.write(f"[{_timestamp()}] {title}\n")
        handle.write("─" * 100 + "\n")
        for line in lines:
            handle.write(f"{line}\n")
        handle.write("\n")


def _append_text_block(log_path: Path, title: str, text: str) -> None:
    safe_text = text if text.strip() else "<empty>"
    log_path.parent.mkdir(parents=True, exist_ok=True)
    with log_path.open("a", encoding="utf-8") as handle:
        handle.write("═" * 100 + "\n")
        handle.write(f"[{_timestamp()}] {title}\n")
        handle.write("─" * 100 + "\n")
        handle.write(f"{safe_text}\n\n")


def _log_exception(log_path: Path, title: str, exc: Exception) -> None:
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
    candidates = ["gmx", "gmx_mpi"]
    for candidate in candidates:
        if shutil.which(candidate):
            return candidate
    raise FileNotFoundError(
        "Could not find a GROMACS executable in PATH. "
        "Expected one of: gmx, gmx_mpi."
    )


def count_atoms_in_structure_file(pdb_path: str | Path) -> int:
    """Count ATOM and HETATM records in a PDB-like structure file."""
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
    if path is None:
        return False
    path = Path(path)
    return path.is_file() and path.stat().st_size > 0


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

    Priority: explicit filler path → default filler final model → MODELLER →
    AlphaFold → original component protein.
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


# ============================================================================
# PDB editing helpers
# ============================================================================


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


_TWO_CHAR_ELEMENTS = frozenset({
    "CL", "BR", "FE", "ZN", "MG", "CA", "MN", "CU", "CO", "NI",
    "NA", "HG", "SE", "SI",
})


def _infer_element_from_atom_name(atom_name_field: str) -> str:
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


def _remove_backbone_incomplete_residues(
    input_pdb: Path,
    output_pdb: Path,
    log_path: Path | None = None,
) -> list[str]:
    """Remove residues missing backbone atoms (N, CA, C, O) from a PDB file."""
    _BACKBONE = {"N", "CA", "C", "O"}
    _SKIP_RESNAMES = {
        "HOH", "WAT", "SOL",
        "NA", "K", "CL", "MG", "CA", "ZN", "FE", "MN", "CU", "CO", "NI",
    }

    lines = input_pdb.read_text(encoding="utf-8", errors="replace").splitlines(keepends=True)

    res_atoms: dict[tuple[str, str, str, str], set[str]] = {}
    res_order: list[tuple[str, str, str, str]] = []
    for line in lines:
        if not (line.startswith("ATOM  ") or line.startswith("HETATM")):
            continue
        if len(line) < 27:
            continue
        chain = line[21]
        resnum = line[22:26]
        icode = line[26]
        resname = line[17:20].strip().upper()
        atom_name = line[12:16].strip().upper()
        key = (chain, resnum, icode, resname)
        if key not in res_atoms:
            res_atoms[key] = set()
            res_order.append(key)
        res_atoms[key].add(atom_name)

    removed: list[str] = []
    bad_keys: set[tuple[str, str, str, str]] = set()
    for key in res_order:
        chain, resnum, icode, resname = key
        if resname in _SKIP_RESNAMES:
            continue
        present = res_atoms[key]
        missing_backbone = _BACKBONE - present
        if missing_backbone:
            label = f"{resname} {resnum.strip()}{icode.strip()} chain {chain}"
            removed.append(
                f"removed incomplete residue {label} "
                f"(missing {', '.join(sorted(missing_backbone))})"
            )
            bad_keys.add(key)

    if not removed:
        if input_pdb != output_pdb:
            output_pdb.write_text("".join(lines), encoding="utf-8")
        return []

    cleaned: list[str] = []
    for line in lines:
        if (line.startswith("ATOM  ") or line.startswith("HETATM")) and len(line) >= 27:
            chain = line[21]
            resnum = line[22:26]
            icode = line[26]
            resname = line[17:20].strip().upper()
            if (chain, resnum, icode, resname) in bad_keys:
                continue
        cleaned.append(line)

    output_pdb.parent.mkdir(parents=True, exist_ok=True)
    output_pdb.write_text("".join(cleaned), encoding="utf-8")

    if log_path is not None:
        try:
            with log_path.open("a", encoding="utf-8") as fh:
                for msg in removed:
                    fh.write(f"[pre-pdb2gmx cleanup] {msg}\n")
        except Exception:
            pass

    return removed


def _fill_pdb_element_column(pdb_path: Path) -> None:
    """Fill element symbol column (cols 77-78) for ATOM/HETATM records."""
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
            has_wrong_h = len(elem_str) == 2 and elem_str[0].upper() == "H"
            has_junk = len(body) > 78
            if elem_missing or has_wrong_h or has_junk:
                atom_name_field = body[12:16] if len(body) >= 16 else ""
                elem = _infer_element_from_atom_name(atom_name_field)
                if not elem and not elem_missing:
                    elem = elem_str
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
    """Add TIP3P hydrogens to crystal waters using gmx pdb2gmx."""
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


def _run_propka(pdb_path: Path, output_dir: Path | None = None) -> object | None:
    """Run PROPKA on pdb_path; return the MolecularContainer or None on error."""
    try:
        import propka.run  # type: ignore[import]
    except ImportError:
        return None

    if output_dir is not None:
        output_dir.mkdir(parents=True, exist_ok=True)
        output_pka = str(output_dir / (pdb_path.stem + ".pka"))
    else:
        output_pka = str(pdb_path.with_suffix(".pka"))

    try:
        buf = io.StringIO()
        with contextlib.redirect_stdout(buf), contextlib.redirect_stderr(buf):
            mol = propka.run.single(
                str(pdb_path),
                optargs=["--quiet", "-o", output_pka],
            )
        return mol
    except SystemExit:
        pass
    except Exception:
        return None

    target_dir = output_dir if output_dir is not None else pdb_path.parent
    target_dir.mkdir(parents=True, exist_ok=True)
    _orig_cwd = os.getcwd()
    try:
        os.chdir(str(target_dir))
        buf = io.StringIO()
        with contextlib.redirect_stdout(buf), contextlib.redirect_stderr(buf):
            mol = propka.run.single(str(pdb_path), optargs=["--quiet"])
        return mol
    except (Exception, SystemExit):
        return None
    finally:
        try:
            os.chdir(_orig_cwd)
        except Exception:
            pass


def _get_propka_his_assignments(
    mol: object,
    ph: float,
) -> list[dict[str, str | int | float]]:
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
    return sorted(
        assignments,
        key=lambda d: (str(d["chain"]), int(d["res_num"]), str(d["icode"])),
    )


def predict_protonation_states(
    pdb_path: Path | str,
    ph: float = 7.4,
) -> dict[tuple[str, int, str], str]:
    """Return PROPKA-based HIP rename keys for residues with pKa > ph."""
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
    if not renames:
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
    """Run ``gmx pdb2gmx`` and write protonated coordinates."""
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
        "-f", str(input_pdb),
        "-o", output_pdb.name,
        "-p", topology_output_path.name,
        "-i", position_restraints_output_path.name,
        "-ff", ff,
        "-water", water_model,
        "-chainsep", chain_separation,
        "-merge", merge,
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
