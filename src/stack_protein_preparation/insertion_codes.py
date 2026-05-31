# /home/grheco/repositorios/stack_protein_prep/src/stack_protein_preparation/insertion_codes.py

"""
insertion_codes.py

Handle insertion codes in PDB files.

Responsibilities
----------------
- detect insertion codes
- remove them using pdb_delinsertion
- provide standardized status output

Design
------
- no structural modification beyond insertion handling
- returns structured result dict for pipeline integration
- preserves existing public function interfaces
- adds structured per-protein logging and compact terminal output

Important
---------
- does NOT perform renumbering
- does NOT depend on alignment/mapping
"""

from __future__ import annotations

import logging
import shutil
import subprocess
from pathlib import Path

STATUS_NONE = "none"
STATUS_SUCCESS = "success"
STATUS_FAILED = "failed"

COLOR_GREEN = "green"
COLOR_RED = "red"

MODULE_NAME = "insertion_codes"


# ============================================================================
# logging helpers
# ============================================================================


def _infer_pdb_id_from_path(pdb_path: Path) -> str:
    """
    Infer the protein/PDB identifier from a pipeline-style path.

    Expected typical layouts:
    - data/proteins/<PDB_ID>/<PDB_ID>.pdb
    - data/proteins/<PDB_ID>/something/<file>.pdb

    Fallback:
    - use pdb_path.stem uppercased
    """
    resolved_parts = pdb_path.resolve().parts
    if "proteins" in resolved_parts:
        proteins_index = resolved_parts.index("proteins")
        if proteins_index + 1 < len(resolved_parts):
            return resolved_parts[proteins_index + 1].upper()

    parent_name = pdb_path.parent.name.strip()
    if parent_name:
        return parent_name.upper()

    return pdb_path.stem.upper()


def _infer_logs_root_from_path(pdb_path: Path) -> Path:
    """
    Infer the shared logs root.

    Preferred:
        .../data/proteins/logs/

    If the file path is somewhere below a 'proteins' directory, use that.
    Otherwise fall back to:
        <input parent>/logs
    """
    resolved = pdb_path.resolve()
    parts = resolved.parts

    if "proteins" in parts:
        proteins_index = parts.index("proteins")
        proteins_root = Path(*parts[: proteins_index + 1])
        return proteins_root / "logs"

    return pdb_path.parent / "logs"


def _get_module_log_path(pdb_path: Path) -> Path:
    """
    Build the standard module log path for this protein and module.

    Example:
        data/proteins/logs/1ABC/insertion_codes.log
    """
    pdb_id = _infer_pdb_id_from_path(pdb_path)
    logs_root = _infer_logs_root_from_path(pdb_path)
    return logs_root / pdb_id / f"{MODULE_NAME}.log"


def _get_logger_for_pdb(pdb_path: Path) -> logging.Logger:
    """
    Return a file-backed logger dedicated to one protein/module log file.
    """
    log_path = _get_module_log_path(pdb_path)
    log_path.parent.mkdir(parents=True, exist_ok=True)

    logger_name = f"{MODULE_NAME}:{log_path}"
    logger = logging.getLogger(logger_name)
    logger.setLevel(logging.INFO)
    logger.propagate = False

    if not logger.handlers:
        handler = logging.FileHandler(log_path, encoding="utf-8")
        formatter = logging.Formatter(
            fmt="%(asctime)s | %(levelname)-7s | %(message)s",
            datefmt="%Y-%m-%d %H:%M:%S",
        )
        handler.setFormatter(formatter)
        logger.addHandler(handler)

    return logger


def _log_block(
    logger: logging.Logger,
    title: str,
    content: str,
) -> None:
    """
    Write a clearly separated multi-line block into the log file.
    """
    safe_content = content if content.strip() else "<empty>"
    logger.info("┏━ %s", title)
    for line in safe_content.splitlines():
        logger.info("┃ %s", line)
    logger.info("┗━ end %s", title)


def _print_header(pdb_id: str, message: str) -> None:
    """
    Compact but more structured terminal header.
    """
    print(f"┏━━ [{MODULE_NAME}] {pdb_id}")
    print(f"┗━ {message}")


def _print_item(label: str, value: str) -> None:
    """
    Structured terminal body line.
    """
    print(f"   ├─ {label}: {value}")


def _print_result(status: str, message: str) -> None:
    """
    Terminal result line.
    """
    print(f"   └─ result: {status} | {message}")


# ============================================================================
# core logic
# ============================================================================


def pdb_has_insertion_codes(pdb_path: Path) -> bool:
    """
    Detect insertion codes directly from a PDB file.

    In the PDB format, the insertion code is column 27 (1-based),
    i.e. line[26] in Python 0-based indexing.
    """
    with pdb_path.open("r", encoding="utf-8") as handle:
        for line in handle:
            if not (line.startswith("ATOM") or line.startswith("HETATM")):
                continue

            if len(line) < 27:
                continue

            insertion_code = line[26].strip()
            if insertion_code:
                return True

    return False


def _python_delinsertion(pdb_text: str) -> str:
    """Pure-Python replacement for pdb_delinsertion from pdb-tools.

    For each chain, maps every unique (resnum_str, icode) tuple to a new
    sequential integer starting from the first residue number in that chain.
    Insertion-code letters are removed (column 27 set to space).
    """
    lines = pdb_text.splitlines(keepends=True)

    chain_order: dict[str, list[tuple[str, str]]] = {}
    for line in lines:
        if not (line.startswith("ATOM  ") or line.startswith("HETATM")):
            continue
        if len(line) < 27:
            continue
        chain = line[21]
        resnum_str = line[22:26]
        icode = line[26]
        key = (resnum_str, icode)
        if chain not in chain_order:
            chain_order[chain] = []
        if key not in chain_order[chain]:
            chain_order[chain].append(key)

    chain_renum: dict[str, dict[tuple[str, str], str]] = {}
    for chain, keys in chain_order.items():
        chain_renum[chain] = {}
        try:
            counter = int(keys[0][0].strip())
        except (ValueError, IndexError):
            counter = 1
        for key in keys:
            chain_renum[chain][key] = f"{counter:4d}"
            counter += 1

    result: list[str] = []
    for line in lines:
        if (line.startswith("ATOM  ") or line.startswith("HETATM")) and len(line) >= 27:
            chain = line[21]
            resnum_str = line[22:26]
            icode = line[26]
            key = (resnum_str, icode)
            if chain in chain_renum and key in chain_renum[chain]:
                new_resnum = chain_renum[chain][key]
                line = line[:22] + new_resnum + " " + line[27:]
        result.append(line)

    return "".join(result)


def run_pdb_delinsertion(
    input_pdb_path: Path,
    output_pdb_path: Path,
) -> tuple[bool, str]:
    """
    Run pdb_delinsertion on one file and write output to output_pdb_path.

    Returns
    -------
    tuple[bool, str]
        (success, message)
    """
    logger = _get_logger_for_pdb(input_pdb_path)
    pdb_id = _infer_pdb_id_from_path(input_pdb_path)

    command = ["pdb_delinsertion", str(input_pdb_path)]

    _print_header(pdb_id, "running pdb_delinsertion")
    _print_item("input", str(input_pdb_path))
    _print_item("output", str(output_pdb_path))
    _print_item("command", " ".join(command))

    logger.info("Starting pdb_delinsertion run")
    logger.info("Input PDB: %s", input_pdb_path)
    logger.info("Output PDB: %s", output_pdb_path)
    logger.info("Command: %s", " ".join(command))

    try:
        result = subprocess.run(
            command,
            capture_output=True,
            text=True,
            check=False,
        )
    except FileNotFoundError:
        logger.info("pdb_delinsertion not found in PATH; using Python fallback")
        pdb_text = input_pdb_path.read_text(encoding="utf-8", errors="replace")
        delinsertion_output = _python_delinsertion(pdb_text)
        if not delinsertion_output.strip():
            message = "Python delinsertion produced empty output"
            logger.error(message)
            _print_result(STATUS_FAILED, message)
            return False, message
        output_pdb_path.parent.mkdir(parents=True, exist_ok=True)
        output_pdb_path.write_text(delinsertion_output, encoding="utf-8")
        message = "pdb_delinsertion completed (Python fallback)"
        logger.info(message)
        _print_result(STATUS_SUCCESS, message)
        return True, message

    logger.info("Return code: %s", result.returncode)
    _log_block(logger, "stdout", result.stdout)
    _log_block(logger, "stderr", result.stderr)

    if result.returncode != 0:
        message = result.stderr.strip() or "pdb_delinsertion failed"
        logger.error("pdb_delinsertion failed: %s", message)
        _print_result(STATUS_FAILED, message)
        return False, message

    stdout = result.stdout
    if not stdout.strip():
        message = "pdb_delinsertion returned empty stdout"
        logger.error(message)
        _print_result(STATUS_FAILED, message)
        return False, message

    output_pdb_path.parent.mkdir(parents=True, exist_ok=True)
    output_pdb_path.write_text(stdout, encoding="utf-8")

    message = "pdb_delinsertion completed"
    logger.info(message)
    _print_result(STATUS_SUCCESS, message)
    return True, message


def copy_pdb_to_output(
    input_pdb_path: Path,
    output_pdb_path: Path,
) -> tuple[bool, str]:
    """
    Copy input PDB to output location.
    """
    logger = _get_logger_for_pdb(input_pdb_path)
    pdb_id = _infer_pdb_id_from_path(input_pdb_path)

    _print_header(pdb_id, "no insertion codes detected -> copying input")
    _print_item("input", str(input_pdb_path))
    _print_item("output", str(output_pdb_path))

    logger.info("No insertion codes detected; copying input PDB unchanged")
    logger.info("Input PDB: %s", input_pdb_path)
    logger.info("Output PDB: %s", output_pdb_path)

    output_pdb_path.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(input_pdb_path, output_pdb_path)

    message = "copied original pdb (no insertion codes)"
    logger.info(message)
    _print_result(STATUS_NONE, message)
    return True, message


def process_pdb_for_delinsertion(
    input_pdb_path: Path,
    output_pdb_path: Path,
) -> dict[str, str | bool]:
    """
    Process one PDB:
    - if insertion codes exist -> run pdb_delinsertion
    - else -> copy file unchanged

    Status convention
    -----------------
    - none    -> green
    - success -> green
    - failed  -> red

    Returns
    -------
    dict[str, str | bool]
        Structured result dictionary.
    """
    logger = _get_logger_for_pdb(input_pdb_path)
    pdb_id = _infer_pdb_id_from_path(input_pdb_path)
    log_path = _get_module_log_path(input_pdb_path)

    _print_header(pdb_id, "checking insertion codes")
    _print_item("input", str(input_pdb_path))
    _print_item("output", str(output_pdb_path))
    _print_item("log", str(log_path))

    logger.info("============================================================")
    logger.info("Starting insertion code processing")
    logger.info("Input PDB: %s", input_pdb_path)
    logger.info("Output PDB: %s", output_pdb_path)

    had_insertion_codes = pdb_has_insertion_codes(input_pdb_path)

    logger.info("Insertion codes detected: %s", had_insertion_codes)
    _print_item("had_insertion_codes", str(had_insertion_codes))

    if had_insertion_codes:
        success, message = run_pdb_delinsertion(
            input_pdb_path=input_pdb_path,
            output_pdb_path=output_pdb_path,
        )
        if success:
            status = STATUS_SUCCESS
            color = COLOR_GREEN
        else:
            status = STATUS_FAILED
            color = COLOR_RED
    else:
        success, message = copy_pdb_to_output(
            input_pdb_path=input_pdb_path,
            output_pdb_path=output_pdb_path,
        )
        status = STATUS_NONE
        color = COLOR_GREEN

    logger.info("Final status: %s", status)
    logger.info("Color: %s", color)
    logger.info("Success: %s", success)
    logger.info("Message: %s", message)
    logger.info("Finished insertion code processing")

    return {
        "input_pdb": str(input_pdb_path),
        "output_pdb": str(output_pdb_path),
        "had_insertion_codes": had_insertion_codes,
        "status": status,
        "color": color,
        "success": success,
        "message": message,
    }


def _preferred_pdb_candidates_for_protein_dir(protein_dir: Path) -> list[Path]:
    """
    Return preferred PDB candidates in robust priority order.

    Priority
    --------
    1. <PDB_ID>.pdb directly inside protein_dir
    2. any other direct .pdb inside protein_dir
    3. first-level subdirectory PDBs
    """
    pdb_id = protein_dir.name.upper()
    preferred_candidates: list[Path] = []

    exact_top_level = protein_dir / f"{pdb_id}.pdb"
    if exact_top_level.exists():
        preferred_candidates.append(exact_top_level)

    direct_pdb_files = sorted(
        path for path in protein_dir.glob("*.pdb") if path != exact_top_level
    )
    preferred_candidates.extend(direct_pdb_files)

    subdirs = sorted(path for path in protein_dir.iterdir() if path.is_dir())
    for subdir in subdirs:
        exact_subdir = subdir / f"{pdb_id}.pdb"
        if exact_subdir.exists():
            preferred_candidates.append(exact_subdir)

        other_subdir_pdbs = sorted(
            path for path in subdir.glob("*.pdb") if path != exact_subdir
        )
        preferred_candidates.extend(other_subdir_pdbs)

    return preferred_candidates


def find_input_pdb_for_protein(protein_dir: Path) -> Path | None:
    """
    Find the most likely original input PDB for one protein directory.

    Search order:
    1. exact top-level file: data/proteins/<PDB_ID>/<PDB_ID>.pdb
    2. any other top-level .pdb
    3. exact <PDB_ID>.pdb in first-level subdirectories
    4. any other first-level subdirectory .pdb

    This is more robust than simply taking the first alphabetical .pdb file,
    because later pipeline stages may create additional derived PDB files.
    """
    pdb_id = protein_dir.name.upper()
    log_probe_path = protein_dir / f"{pdb_id}.pdb"
    logger = _get_logger_for_pdb(log_probe_path)

    logger.info("Searching preferred input PDB for protein_dir=%s", protein_dir)

    if not protein_dir.exists():
        logger.warning("Protein directory does not exist: %s", protein_dir)
        return None

    candidate_paths = _preferred_pdb_candidates_for_protein_dir(protein_dir)

    if not candidate_paths:
        logger.warning("No PDB candidates found for protein_dir=%s", protein_dir)
        return None

    logger.info("Candidate count: %d", len(candidate_paths))
    for index, candidate_path in enumerate(candidate_paths, start=1):
        logger.info("Candidate %d: %s", index, candidate_path)

    chosen_path = candidate_paths[0]
    logger.info("Chosen input PDB: %s", chosen_path)
    return chosen_path
