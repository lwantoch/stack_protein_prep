"""
/home/grheco/repositorios/stack_protein_prep/src/stack_protein_preparation/metall_params.py

Metal-center parametrization preparation module.

Purpose
-------
Prepare a reduced PDB around metal-containing systems for downstream
MCPB.py / Gaussian parametrization.

Current responsibilities
------------------------
1. Detect whether a protein directory contains metals
2. Create the directory:
       data/proteins/<PDB_ID>/metall_params/
   only when metals are present
3. Build:
       metall_params/tmp_param.pdb
   by concatenating the best available protein structure with optional
   water / ligand / metal component files
4. Write a headless UCSF Chimera script that runs a contact / clash search
   and saves the result to:
       metall_params/contacts.data
5. Execute Chimera in terminal mode
6. Return a structured result dict for pipeline-state integration
7. Print extremely verbose debug information to screen for every step

Notes
-----
- This version is intentionally noisy for debugging.
- It prints paths, file-existence checks, file sizes, selected inputs,
  generated files, Chimera command, return code, and final outcome.
- This module targets classic UCSF Chimera, not ChimeraX.
- Current metal filename convention in this project is:
      <PDB_ID>_metal.pdb
"""

from __future__ import annotations

import os
import shutil
import subprocess
from pathlib import Path
from typing import Any

# ============================================================================
# USER-ADJUSTABLE CONSTANTS
# ============================================================================

CHIMERA_FINDCLASH_COMMAND = (
    "findclash #0 test other overlap -0.4 hbond 0.0 bondSeparation 2 "
    "saveFile contacts.data"
)

CHIMERA_PRE_COMMANDS: list[str] = []
CHIMERA_POST_COMMANDS: list[str] = []

DEFAULT_CHIMERA_EXECUTABLE_CANDIDATES = [
    "chimera",
    "/usr/bin/chimera",
    "/home/grheco/.local/bin/chimera",
    "/opt/UCSF/Chimera64-1.17/bin/chimera",
    "/opt/UCSF/Chimera64-1.16/bin/chimera",
]

# ============================================================================
# DEBUG PRINTS
# ============================================================================


def _p(label: str, value: Any) -> None:
    print(f"[metall_params][DEBUG] {label}: {value}")


def _section(title: str) -> None:
    print(f"\n[metall_params] ===== {title} =====")


# ============================================================================
# FILE SELECTION
# ============================================================================


def _choose_best_protein_input(components_dir: Path, pdb_id: str) -> Path | None:
    """
    Return the preferred protein structure for metal parametrization.

    Priority:
    1. <PDB_ID>_protein_final.pdb
    2. <PDB_ID>_protein_as_Amber.pdb
    3. <PDB_ID>_proteinH.pdb
    4. <PDB_ID>_protein.pdb
    """
    _section("_choose_best_protein_input")

    candidates = [
        components_dir / f"{pdb_id}_protein_final.pdb",
        components_dir / f"{pdb_id}_protein_as_Amber.pdb",
        components_dir / f"{pdb_id}_proteinH.pdb",
        components_dir / f"{pdb_id}_protein.pdb",
    ]

    for path in candidates:
        _p("candidate_protein_input", path)
        _p("candidate_exists", path.is_file())
        if path.is_file():
            _p("selected_protein_input", path)
            return path

    _p("selected_protein_input", None)
    return None


def _get_optional_component_paths(
    components_dir: Path,
    pdb_id: str,
) -> dict[str, Path | None]:
    """
    Return optional component files if present.
    """
    _section("_get_optional_component_paths")

    water_path = components_dir / f"{pdb_id}_water.pdb"
    ligand_path = components_dir / f"{pdb_id}_ligand.pdb"
    metal_path = components_dir / f"{pdb_id}_metal.pdb"

    _p("water_path", water_path)
    _p("water_exists", water_path.is_file())
    _p("ligand_path", ligand_path)
    _p("ligand_exists", ligand_path.is_file())
    _p("metal_path", metal_path)
    _p("metal_exists", metal_path.is_file())

    result = {
        "water": water_path if water_path.is_file() else None,
        "ligand": ligand_path if ligand_path.is_file() else None,
        "metal": metal_path if metal_path.is_file() else None,
    }

    _p("optional_component_result", result)
    return result


def _has_metal_file(components_dir: Path, pdb_id: str) -> bool:
    """
    Decide whether the current protein directory contains metals.

    Current rule:
    - metal file exists and is non-empty
    """
    _section("_has_metal_file")

    metal_path = components_dir / f"{pdb_id}_metal.pdb"
    exists = metal_path.is_file()
    size = metal_path.stat().st_size if exists else 0

    _p("metal_path_checked", metal_path)
    _p("metal_exists", exists)
    _p("metal_size", size)

    decision = exists and size > 0
    _p("has_metal_file_decision", decision)
    return decision


# ============================================================================
# PDB BUILDING
# ============================================================================


def _read_pdb_payload_lines(pdb_path: Path) -> list[str]:
    """
    Read only payload lines that should be kept in a merged PDB.

    Dropped:
    - END
    - ENDMDL
    """
    _section("_read_pdb_payload_lines")
    _p("pdb_path", pdb_path)
    _p("exists", pdb_path.exists())

    payload_lines: list[str] = []

    with pdb_path.open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            stripped = line.strip()
            if stripped in {"END", "ENDMDL"}:
                continue
            payload_lines.append(line)

    _p("payload_line_count", len(payload_lines))
    return payload_lines


def _write_combined_tmp_param_pdb(
    output_path: Path,
    protein_path: Path,
    water_path: Path | None,
    ligand_path: Path | None,
    metal_path: Path | None,
) -> dict[str, Any]:
    """
    Concatenate the selected component PDB files into tmp_param.pdb.
    """
    _section("_write_combined_tmp_param_pdb")
    _p("output_path", output_path)
    _p("protein_path", protein_path)
    _p("water_path", water_path)
    _p("ligand_path", ligand_path)
    _p("metal_path", metal_path)

    source_paths: list[Path] = [protein_path]
    if water_path is not None:
        source_paths.append(water_path)
    if ligand_path is not None:
        source_paths.append(ligand_path)
    if metal_path is not None:
        source_paths.append(metal_path)

    _p("source_paths", source_paths)

    line_count = 0

    with output_path.open("w", encoding="utf-8") as out_handle:
        for source_path in source_paths:
            _p("adding_source_path", source_path)
            lines = _read_pdb_payload_lines(source_path)

            for line in lines:
                out_handle.write(line)
                line_count += 1

            if lines and not lines[-1].startswith("TER"):
                out_handle.write("TER\n")
                line_count += 1
                _p("added_TER_after", source_path)

        out_handle.write("END\n")
        line_count += 1

    _p("line_count_written", line_count)
    _p("tmp_param_exists_after_write", output_path.exists())
    _p(
        "tmp_param_size_after_write",
        output_path.stat().st_size if output_path.exists() else 0,
    )

    return {
        "source_paths": [str(path) for path in source_paths],
        "line_count_written": line_count,
    }


# ============================================================================
# CHIMERA
# ============================================================================


def _resolve_chimera_executable(chimera_executable: str | None = None) -> str | None:
    """
    Resolve Chimera executable.

    Priority:
    1. explicit function argument
    2. CHIMERA_EXECUTABLE environment variable
    3. fallback candidate list
    """
    _section("_resolve_chimera_executable")

    if chimera_executable:
        _p("explicit_arg", chimera_executable)

        explicit_resolved = shutil.which(chimera_executable)
        _p("explicit_arg_shutil_which", explicit_resolved)
        if explicit_resolved:
            return explicit_resolved

        explicit_path = Path(chimera_executable)
        _p("explicit_arg_is_file", explicit_path.is_file())
        if explicit_path.is_file():
            return str(explicit_path)

    env_path = os.environ.get("CHIMERA_EXECUTABLE")
    _p("env_CHIMERA_EXECUTABLE", env_path)

    if env_path:
        env_resolved = shutil.which(env_path)
        _p("env_path_shutil_which", env_resolved)
        if env_resolved:
            return env_resolved

        env_file = Path(env_path)
        _p("env_path_is_file", env_file.is_file())
        if env_file.is_file():
            return str(env_file)

    for candidate in DEFAULT_CHIMERA_EXECUTABLE_CANDIDATES:
        _p("candidate_chimera", candidate)

        resolved = shutil.which(candidate)
        _p("candidate_shutil_which", resolved)
        if resolved:
            _p("using_candidate_from_PATH", resolved)
            return resolved

        candidate_path = Path(candidate)
        _p("candidate_is_file", candidate_path.is_file())
        if candidate_path.is_file():
            _p("using_candidate_as_file", candidate_path)
            return str(candidate_path)

    _p("resolved_chimera", None)
    return None


def _write_chimera_python_script(
    script_path: Path,
    tmp_param_pdb_name: str = "tmp_param.pdb",
    contacts_file_name: str = "contacts.data",
) -> None:
    """
    Write a headless Chimera Python script.
    """
    _section("_write_chimera_python_script")
    _p("script_path", script_path)
    _p("tmp_param_pdb_name", tmp_param_pdb_name)
    _p("contacts_file_name", contacts_file_name)
    _p("findclash_command", CHIMERA_FINDCLASH_COMMAND)
    _p("pre_commands", CHIMERA_PRE_COMMANDS)
    _p("post_commands", CHIMERA_POST_COMMANDS)

    pre_commands_block = "\n".join(f'rc("{cmd}")' for cmd in CHIMERA_PRE_COMMANDS)
    post_commands_block = "\n".join(f'rc("{cmd}")' for cmd in CHIMERA_POST_COMMANDS)

    script_text = f"""# Auto-generated by metall_params.py
from chimera import runCommand as rc

rc("open {tmp_param_pdb_name}")
{pre_commands_block}
rc("{CHIMERA_FINDCLASH_COMMAND}")
{post_commands_block}
rc("stop now")
"""

    script_path.write_text(script_text, encoding="utf-8")

    _p("chimera_script_written", script_path.exists())
    _p(
        "chimera_script_size",
        script_path.stat().st_size if script_path.exists() else 0,
    )
    _p("chimera_script_text", script_text)


def _run_chimera_headless(
    metall_params_dir: Path,
    chimera_executable: str,
    script_path: Path,
    log_path: Path,
) -> subprocess.CompletedProcess[str]:
    """
    Run UCSF Chimera in headless mode inside the metall_params directory.
    """
    _section("_run_chimera_headless")
    _p("metall_params_dir", metall_params_dir)
    _p("chimera_executable", chimera_executable)
    _p("script_path", script_path)
    _p("log_path", log_path)

    command = [chimera_executable, "--nogui", str(script_path.name)]
    _p("chimera_command", command)

    result = subprocess.run(
        command,
        cwd=metall_params_dir,
        capture_output=True,
        text=True,
        check=False,
    )

    _p("chimera_returncode", result.returncode)
    _p("chimera_stdout", result.stdout)
    _p("chimera_stderr", result.stderr)

    log_text = [
        "=== COMMAND ===",
        " ".join(command),
        "",
        "=== STDOUT ===",
        result.stdout or "",
        "",
        "=== STDERR ===",
        result.stderr or "",
        "",
        f"=== RETURN CODE ===\n{result.returncode}\n",
    ]
    log_path.write_text("\n".join(log_text), encoding="utf-8")

    _p("chimera_log_written", log_path.exists())
    _p("chimera_log_size", log_path.stat().st_size if log_path.exists() else 0)

    return result


# ============================================================================
# PUBLIC API
# ============================================================================


def run_metal_parametrization_for_protein_dir(
    protein_dir: str | Path,
    chimera_executable: str | None = None,
) -> dict[str, Any]:
    """
    Run the metal-preparation step for one protein directory.
    """
    _section("run_metal_parametrization_for_protein_dir START")

    protein_dir = Path(protein_dir).resolve()
    pdb_id = protein_dir.name
    components_dir = protein_dir / "components"
    metall_params_dir = protein_dir / "metall_params"

    _p("protein_dir", protein_dir)
    _p("pdb_id", pdb_id)
    _p("components_dir", components_dir)
    _p("components_dir_exists", components_dir.is_dir())
    _p("metall_params_dir", metall_params_dir)
    _p("requested_chimera_executable", chimera_executable)

    result: dict[str, Any] = {
        "status": "required",
        "pdb_id": pdb_id,
        "protein_dir": str(protein_dir),
        "components_dir": str(components_dir),
        "metall_params_dir": str(metall_params_dir),
        "used_protein_input": None,
        "used_water_input": None,
        "used_ligand_input": None,
        "used_metal_input": None,
        "tmp_param_pdb": None,
        "chimera_script": None,
        "contacts_file": None,
        "chimera_log": None,
        "chimera_executable": None,
        "message": "",
        "merge_source_paths": [],
        "tmp_param_line_count": 0,
        "debug": {
            "components_dir_exists": False,
            "has_metal_file": False,
            "protein_input_found": False,
            "water_input_found": False,
            "ligand_input_found": False,
            "metal_input_found": False,
            "tmp_param_written": False,
            "chimera_script_written": False,
            "chimera_resolved": False,
            "chimera_returncode_zero": False,
            "contacts_file_exists": False,
            "contacts_file_nonempty": False,
            "findclash_command": CHIMERA_FINDCLASH_COMMAND,
            "tmp_param_size": 0,
            "contacts_file_size": 0,
        },
    }

    result["debug"]["components_dir_exists"] = components_dir.is_dir()

    if not components_dir.is_dir():
        result["status"] = "failed"
        result["message"] = f"Missing components directory: {components_dir}"
        _p("FINAL_STATUS", result["status"])
        _p("FINAL_MESSAGE", result["message"])
        return result

    result["debug"]["has_metal_file"] = _has_metal_file(components_dir, pdb_id)
    _p("debug.has_metal_file", result["debug"]["has_metal_file"])

    if not result["debug"]["has_metal_file"]:
        result["status"] = "skipped"
        result["message"] = "No metal file found. metall_params step skipped."
        _p("FINAL_STATUS", result["status"])
        _p("FINAL_MESSAGE", result["message"])
        return result

    protein_input = _choose_best_protein_input(components_dir, pdb_id)
    optional_components = _get_optional_component_paths(components_dir, pdb_id)

    water_input = optional_components["water"]
    ligand_input = optional_components["ligand"]
    metal_input = optional_components["metal"]

    result["debug"]["protein_input_found"] = protein_input is not None
    result["debug"]["water_input_found"] = water_input is not None
    result["debug"]["ligand_input_found"] = ligand_input is not None
    result["debug"]["metal_input_found"] = metal_input is not None

    _p("protein_input", protein_input)
    _p("water_input", water_input)
    _p("ligand_input", ligand_input)
    _p("metal_input", metal_input)
    _p("debug.protein_input_found", result["debug"]["protein_input_found"])
    _p("debug.water_input_found", result["debug"]["water_input_found"])
    _p("debug.ligand_input_found", result["debug"]["ligand_input_found"])
    _p("debug.metal_input_found", result["debug"]["metal_input_found"])

    if protein_input is None:
        result["status"] = "failed"
        result["message"] = (
            "No suitable protein input found. Expected one of: "
            f"{pdb_id}_protein_final.pdb, "
            f"{pdb_id}_protein_as_Amber.pdb, "
            f"{pdb_id}_proteinH.pdb, "
            f"{pdb_id}_protein.pdb"
        )
        _p("FINAL_STATUS", result["status"])
        _p("FINAL_MESSAGE", result["message"])
        return result

    if metal_input is None:
        result["status"] = "failed"
        result["message"] = "Metal file missing despite positive metal detection."
        _p("FINAL_STATUS", result["status"])
        _p("FINAL_MESSAGE", result["message"])
        return result

    _section("create metall_params directory")
    metall_params_dir.mkdir(parents=True, exist_ok=True)
    _p("metall_params_dir_exists_after_mkdir", metall_params_dir.exists())

    tmp_param_pdb = metall_params_dir / "tmp_param.pdb"
    chimera_script = metall_params_dir / "chimera_contacts.py"
    contacts_file = metall_params_dir / "contacts.data"
    chimera_log = metall_params_dir / "chimera_run.log"

    _p("tmp_param_pdb", tmp_param_pdb)
    _p("chimera_script", chimera_script)
    _p("contacts_file", contacts_file)
    _p("chimera_log", chimera_log)

    merge_info = _write_combined_tmp_param_pdb(
        output_path=tmp_param_pdb,
        protein_path=protein_input,
        water_path=water_input,
        ligand_path=ligand_input,
        metal_path=metal_input,
    )

    result["debug"]["tmp_param_written"] = (
        tmp_param_pdb.is_file() and tmp_param_pdb.stat().st_size > 0
    )
    result["debug"]["tmp_param_size"] = (
        tmp_param_pdb.stat().st_size if tmp_param_pdb.exists() else 0
    )

    _p("debug.tmp_param_written", result["debug"]["tmp_param_written"])
    _p("debug.tmp_param_size", result["debug"]["tmp_param_size"])
    _p("merge_info", merge_info)

    _write_chimera_python_script(
        script_path=chimera_script,
        tmp_param_pdb_name=tmp_param_pdb.name,
        contacts_file_name=contacts_file.name,
    )

    result["debug"]["chimera_script_written"] = (
        chimera_script.is_file() and chimera_script.stat().st_size > 0
    )
    _p("debug.chimera_script_written", result["debug"]["chimera_script_written"])

    resolved_chimera = _resolve_chimera_executable(
        chimera_executable=chimera_executable
    )

    result["used_protein_input"] = str(protein_input)
    result["used_water_input"] = str(water_input) if water_input else None
    result["used_ligand_input"] = str(ligand_input) if ligand_input else None
    result["used_metal_input"] = str(metal_input) if metal_input else None
    result["tmp_param_pdb"] = str(tmp_param_pdb)
    result["chimera_script"] = str(chimera_script)
    result["contacts_file"] = str(contacts_file)
    result["chimera_log"] = str(chimera_log)
    result["merge_source_paths"] = merge_info["source_paths"]
    result["tmp_param_line_count"] = merge_info["line_count_written"]

    result["debug"]["chimera_resolved"] = resolved_chimera is not None
    _p("resolved_chimera", resolved_chimera)
    _p("debug.chimera_resolved", result["debug"]["chimera_resolved"])

    if resolved_chimera is None:
        result["status"] = "failed"
        result["message"] = (
            "Could not resolve UCSF Chimera executable. "
            "Pass `chimera_executable=` explicitly or set CHIMERA_EXECUTABLE."
        )
        _p("FINAL_STATUS", result["status"])
        _p("FINAL_MESSAGE", result["message"])
        return result

    result["chimera_executable"] = resolved_chimera

    run_result = _run_chimera_headless(
        metall_params_dir=metall_params_dir,
        chimera_executable=resolved_chimera,
        script_path=chimera_script,
        log_path=chimera_log,
    )

    result["debug"]["chimera_returncode_zero"] = run_result.returncode == 0
    _p("debug.chimera_returncode_zero", result["debug"]["chimera_returncode_zero"])

    if run_result.returncode != 0:
        result["status"] = "failed"
        result["message"] = (
            f"Chimera returned non-zero exit code {run_result.returncode}. "
            f"tmp_param_written={result['debug']['tmp_param_written']}, "
            f"chimera_script_written={result['debug']['chimera_script_written']}, "
            f"chimera_resolved={result['debug']['chimera_resolved']}. "
            f"See log: {chimera_log}"
        )
        _p("FINAL_STATUS", result["status"])
        _p("FINAL_MESSAGE", result["message"])
        return result

    result["debug"]["contacts_file_exists"] = contacts_file.is_file()
    result["debug"]["contacts_file_nonempty"] = (
        contacts_file.is_file() and contacts_file.stat().st_size > 0
    )
    result["debug"]["contacts_file_size"] = (
        contacts_file.stat().st_size if contacts_file.exists() else 0
    )

    _p("debug.contacts_file_exists", result["debug"]["contacts_file_exists"])
    _p("debug.contacts_file_nonempty", result["debug"]["contacts_file_nonempty"])
    _p("debug.contacts_file_size", result["debug"]["contacts_file_size"])

    if not contacts_file.is_file() or contacts_file.stat().st_size == 0:
        result["status"] = "warning"
        result["message"] = (
            "Chimera finished but contacts.data is missing or empty. "
            f"contacts_file_exists={result['debug']['contacts_file_exists']}, "
            f"contacts_file_nonempty={result['debug']['contacts_file_nonempty']}. "
            f"See log: {chimera_log}"
        )
        _p("FINAL_STATUS", result["status"])
        _p("FINAL_MESSAGE", result["message"])
        return result

    result["status"] = "success"
    result["message"] = "Metal parametrization preparation finished successfully."

    _p("FINAL_STATUS", result["status"])
    _p("FINAL_MESSAGE", result["message"])
    _p("FINAL_RESULT", result)

    _section("run_metal_parametrization_for_protein_dir END")
    return result


# ============================================================================
# OPTIONAL CLI
# ============================================================================

if __name__ == "__main__":
    import argparse
    import json

    parser = argparse.ArgumentParser(
        description="Prepare tmp_param.pdb and contacts.data for metal-containing protein systems."
    )
    parser.add_argument(
        "protein_dir",
        type=str,
        help="Path to one protein directory, e.g. data/proteins/1W4R",
    )
    parser.add_argument(
        "--chimera",
        dest="chimera_executable",
        type=str,
        default=None,
        help="Optional explicit UCSF Chimera executable path or name.",
    )

    args = parser.parse_args()

    cli_result = run_metal_parametrization_for_protein_dir(
        protein_dir=args.protein_dir,
        chimera_executable=args.chimera_executable,
    )
    print(json.dumps(cli_result, indent=2))
