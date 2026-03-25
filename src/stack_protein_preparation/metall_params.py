"""
/home/grheco/repositorios/stack_protein_prep/src/stack_protein_preparation/metall_params.py

Metal-center parametrization preparation module.

Purpose
-------
Prepare a reduced structural input around metal-containing protein systems for
later MCPB.py / Gaussian / metal-center parametrization workflows.

Responsibilities
----------------
1. Detect whether a protein directory contains a metal component
2. Create:
       data/proteins/<PDB_ID>/metall_params/
   only when a metal file is present
3. Build:
       metall_params/tmp_param.pdb
   by concatenating the best available protein structure with optional
   water / ligand component files
4. Copy the metal component separately to:
       metall_params/metal_only.pdb
5. Write a headless UCSF Chimera script that:
   - opens the full environment system as model #0
   - opens the metal-only file as model #1
   - runs a contact / clash search from metal to environment
   - saves the result to:
         metall_params/contacts.data
6. Execute Chimera in terminal mode
7. Return a structured result dict for pipeline-state integration

Important filename convention
-----------------------------
This module assumes the metal component file is named:

    <PDB_ID>_metal.pdb

Example:
    2AFX_metal.pdb

Important model separation
--------------------------
tmp_param.pdb intentionally does NOT contain the metal atoms.

Model usage in Chimera:
- #0 -> tmp_param.pdb = protein + optional water + optional ligand
- #1 -> metal_only.pdb = metal-only component

This avoids false self-contacts such as:
- ZN in #1 against the same ZN already present in #0

Notes
-----
- This module targets classic UCSF Chimera, not ChimeraX.
- The merged tmp_param.pdb may still contain duplicate atom serial numbers if
  the component files were previously generated independently.
"""

from __future__ import annotations

import os
import shutil
import subprocess
from pathlib import Path
from typing import Any

CHIMERA_FINDCLASH_COMMAND = (
    "findclash #1 test #0 overlap -0.4 hbond 0.0 bondSeparation 2 "
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


def _choose_best_protein_input(components_dir: Path, pdb_id: str) -> Path | None:
    """
    Return the preferred protein structure for metal parametrization.

    Priority:
    1. <PDB_ID>_protein_final.pdb
    2. <PDB_ID>_protein_internal_capped.pdb
    3. <PDB_ID>_protein_amber_termini.pdb
    4. <PDB_ID>_protein_as_Amber.pdb
    5. <PDB_ID>_proteinH.pdb
    6. <PDB_ID>_protein.pdb
    """
    candidates = [
        components_dir / f"{pdb_id}_protein_final.pdb",
        components_dir / f"{pdb_id}_protein_internal_capped.pdb",
        components_dir / f"{pdb_id}_protein_amber_termini.pdb",
        components_dir / f"{pdb_id}_protein_as_Amber.pdb",
        components_dir / f"{pdb_id}_proteinH.pdb",
        components_dir / f"{pdb_id}_protein.pdb",
    ]

    for path in candidates:
        if path.is_file():
            return path

    return None


def _get_optional_component_paths(
    components_dir: Path,
    pdb_id: str,
) -> dict[str, Path | None]:
    """
    Return optional component files if present.
    """
    water_path = components_dir / f"{pdb_id}_water.pdb"
    ligand_path = components_dir / f"{pdb_id}_ligand.pdb"
    metal_path = components_dir / f"{pdb_id}_metal.pdb"

    return {
        "water": water_path if water_path.is_file() else None,
        "ligand": ligand_path if ligand_path.is_file() else None,
        "metal": metal_path if metal_path.is_file() else None,
    }


def _has_metal_file(components_dir: Path, pdb_id: str) -> bool:
    """
    Decide whether the current protein directory contains a metal file.

    Current rule:
    - metal file exists and is non-empty
    """
    metal_path = components_dir / f"{pdb_id}_metal.pdb"
    return metal_path.is_file() and metal_path.stat().st_size > 0


def _read_pdb_payload_lines(pdb_path: Path) -> list[str]:
    """
    Read only payload lines that should be kept in a merged PDB.

    Dropped:
    - END
    - ENDMDL

    Kept:
    - everything else, including ATOM / HETATM / TER / REMARK, etc.
    """
    payload_lines: list[str] = []

    with pdb_path.open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            stripped = line.strip()
            if stripped in {"END", "ENDMDL"}:
                continue
            payload_lines.append(line)

    return payload_lines


def _write_combined_tmp_param_pdb(
    output_path: Path,
    protein_path: Path,
    water_path: Path | None,
    ligand_path: Path | None,
) -> dict[str, Any]:
    """
    Concatenate the selected component PDB files into tmp_param.pdb.

    Order:
    1. protein
    2. water (optional)
    3. ligand (optional)

    Important
    ---------
    The metal component is intentionally NOT included here.
    It is written separately to metal_only.pdb so that Chimera can compare:

    - #1 = metal-only
    - #0 = environment without metal
    """
    source_paths: list[Path] = [protein_path]
    if water_path is not None:
        source_paths.append(water_path)
    if ligand_path is not None:
        source_paths.append(ligand_path)

    line_count = 0

    with output_path.open("w", encoding="utf-8") as out_handle:
        for source_path in source_paths:
            lines = _read_pdb_payload_lines(source_path)

            for line in lines:
                out_handle.write(line)
                line_count += 1

            if lines and not lines[-1].startswith("TER"):
                out_handle.write("TER\n")
                line_count += 1

        out_handle.write("END\n")
        line_count += 1

    return {
        "source_paths": [str(path) for path in source_paths],
        "line_count_written": line_count,
    }


def _copy_metal_only_pdb(
    metal_input_path: Path,
    metal_only_output_path: Path,
) -> None:
    """
    Copy the metal component to a dedicated local file for separate opening
    in Chimera as model #1.
    """
    text = metal_input_path.read_text(encoding="utf-8", errors="replace")
    metal_only_output_path.write_text(text, encoding="utf-8")


def _resolve_chimera_executable(chimera_executable: str | None = None) -> str | None:
    """
    Resolve Chimera executable.

    Priority:
    1. explicit function argument
    2. CHIMERA_EXECUTABLE environment variable
    3. fallback candidate list
    """
    if chimera_executable:
        explicit_resolved = shutil.which(chimera_executable)
        if explicit_resolved:
            return explicit_resolved

        explicit_path = Path(chimera_executable)
        if explicit_path.is_file():
            return str(explicit_path)

    env_path = os.environ.get("CHIMERA_EXECUTABLE")
    if env_path:
        env_resolved = shutil.which(env_path)
        if env_resolved:
            return env_resolved

        env_file = Path(env_path)
        if env_file.is_file():
            return str(env_file)

    for candidate in DEFAULT_CHIMERA_EXECUTABLE_CANDIDATES:
        resolved = shutil.which(candidate)
        if resolved:
            return resolved

        candidate_path = Path(candidate)
        if candidate_path.is_file():
            return str(candidate_path)

    return None


def _write_chimera_python_script(
    script_path: Path,
    tmp_param_pdb_name: str = "tmp_param.pdb",
    metal_only_pdb_name: str = "metal_only.pdb",
    contacts_file_name: str = "contacts.data",
) -> None:
    """
    Write a headless Chimera Python script.

    Model layout
    ------------
    #0 -> full environment system (tmp_param.pdb)
    #1 -> metal-only PDB
    """
    pre_commands_block = "\n".join(f'rc("{cmd}")' for cmd in CHIMERA_PRE_COMMANDS)
    post_commands_block = "\n".join(f'rc("{cmd}")' for cmd in CHIMERA_POST_COMMANDS)

    script_text = f"""# Auto-generated by metall_params.py
from chimera import runCommand as rc

rc("open {tmp_param_pdb_name}")
rc("open {metal_only_pdb_name}")
{pre_commands_block}
rc("{CHIMERA_FINDCLASH_COMMAND}")
{post_commands_block}
rc("stop now")
"""

    script_path.write_text(script_text, encoding="utf-8")


def _run_chimera_headless(
    metall_params_dir: Path,
    chimera_executable: str,
    script_path: Path,
    log_path: Path,
) -> subprocess.CompletedProcess[str]:
    """
    Run UCSF Chimera in headless mode inside the metall_params directory.
    """
    command = [chimera_executable, "--nogui", str(script_path.name)]

    result = subprocess.run(
        command,
        cwd=metall_params_dir,
        capture_output=True,
        text=True,
        check=False,
    )

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

    return result


def run_metal_parametrization_for_protein_dir(
    protein_dir: str | Path,
    chimera_executable: str | None = None,
) -> dict[str, Any]:
    """
    Run the metal-preparation step for one protein directory.

    Parameters
    ----------
    protein_dir:
        Path like:
            data/proteins/<PDB_ID>
    chimera_executable:
        Optional explicit Chimera executable path or name.

    Returns
    -------
    dict
        Structured result suitable for pipeline state logging.
    """
    protein_dir = Path(protein_dir).resolve()
    pdb_id = protein_dir.name
    components_dir = protein_dir / "components"
    metall_params_dir = protein_dir / "metall_params"

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
        "metal_only_pdb": None,
        "chimera_script": None,
        "contacts_file": None,
        "chimera_log": None,
        "chimera_executable": None,
        "message": "",
        "merge_source_paths": [],
        "tmp_param_line_count": 0,
    }

    if not components_dir.is_dir():
        result["status"] = "failed"
        result["message"] = f"Missing components directory: {components_dir}"
        return result

    if not _has_metal_file(components_dir, pdb_id):
        result["status"] = "skipped"
        result["message"] = "No metal file found. metall_params step skipped."
        return result

    protein_input = _choose_best_protein_input(components_dir, pdb_id)
    optional_components = _get_optional_component_paths(components_dir, pdb_id)

    water_input = optional_components["water"]
    ligand_input = optional_components["ligand"]
    metal_input = optional_components["metal"]

    if protein_input is None:
        result["status"] = "failed"
        result["message"] = (
            "No suitable protein input found. Expected one of: "
            f"{pdb_id}_protein_final.pdb, "
            f"{pdb_id}_protein_internal_capped.pdb, "
            f"{pdb_id}_protein_amber_termini.pdb, "
            f"{pdb_id}_protein_as_Amber.pdb, "
            f"{pdb_id}_proteinH.pdb, "
            f"{pdb_id}_protein.pdb"
        )
        return result

    if metal_input is None:
        result["status"] = "failed"
        result["message"] = "Metal file missing despite positive metal detection."
        return result

    metall_params_dir.mkdir(parents=True, exist_ok=True)

    tmp_param_pdb = metall_params_dir / "tmp_param.pdb"
    metal_only_pdb = metall_params_dir / "metal_only.pdb"
    chimera_script = metall_params_dir / "chimera_contacts.py"
    contacts_file = metall_params_dir / "contacts.data"
    chimera_log = metall_params_dir / "chimera_run.log"

    merge_info = _write_combined_tmp_param_pdb(
        output_path=tmp_param_pdb,
        protein_path=protein_input,
        water_path=water_input,
        ligand_path=ligand_input,
    )

    _copy_metal_only_pdb(
        metal_input_path=metal_input,
        metal_only_output_path=metal_only_pdb,
    )

    _write_chimera_python_script(
        script_path=chimera_script,
        tmp_param_pdb_name=tmp_param_pdb.name,
        metal_only_pdb_name=metal_only_pdb.name,
        contacts_file_name=contacts_file.name,
    )

    resolved_chimera = _resolve_chimera_executable(
        chimera_executable=chimera_executable
    )

    result["used_protein_input"] = str(protein_input)
    result["used_water_input"] = str(water_input) if water_input else None
    result["used_ligand_input"] = str(ligand_input) if ligand_input else None
    result["used_metal_input"] = str(metal_input)
    result["tmp_param_pdb"] = str(tmp_param_pdb)
    result["metal_only_pdb"] = str(metal_only_pdb)
    result["chimera_script"] = str(chimera_script)
    result["contacts_file"] = str(contacts_file)
    result["chimera_log"] = str(chimera_log)
    result["merge_source_paths"] = merge_info["source_paths"]
    result["tmp_param_line_count"] = merge_info["line_count_written"]

    if resolved_chimera is None:
        result["status"] = "failed"
        result["message"] = (
            "Could not resolve UCSF Chimera executable. "
            "Pass `chimera_executable=` explicitly or set CHIMERA_EXECUTABLE."
        )
        return result

    result["chimera_executable"] = resolved_chimera

    run_result = _run_chimera_headless(
        metall_params_dir=metall_params_dir,
        chimera_executable=resolved_chimera,
        script_path=chimera_script,
        log_path=chimera_log,
    )

    if run_result.returncode != 0:
        result["status"] = "failed"
        result["message"] = (
            f"Chimera returned non-zero exit code {run_result.returncode}. "
            f"See log: {chimera_log}"
        )
        return result

    if not contacts_file.is_file() or contacts_file.stat().st_size == 0:
        result["status"] = "warning"
        result["message"] = (
            "Chimera finished but contacts.data is missing or empty. "
            f"See log: {chimera_log}"
        )
        return result

    result["status"] = "success"
    result["message"] = "Metal parametrization preparation finished successfully."
    return result


if __name__ == "__main__":
    import argparse
    import json

    parser = argparse.ArgumentParser(
        description=(
            "Prepare tmp_param.pdb and contacts.data for metal-containing "
            "protein systems."
        )
    )
    parser.add_argument(
        "protein_dir",
        type=str,
        help="Path to one protein directory, e.g. data/proteins/2AFX",
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
