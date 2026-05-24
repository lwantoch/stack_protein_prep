"""
gaussian_hpc.py — SLURM-based Gaussian submission and waiting for FRUTON pipeline.

Handles the full QM computation loop for both parametrization branches:

MCPB (metal centers):
    commands.sh  (MCPB.py -s 1  →  .com files)
    ↓  sbatch 3 independent jobs in parallel
    ↓  wait (squeue poll)
    run_after_gaussian.sh  (formchk + MCPB.py -s 2,3,4 + tleap)

RESP (non-standard residues):
    commands.sh  (antechamber  →  opt.com + esp.com)
    ↓  sbatch opt
    ↓  sbatch esp  (with --dependency=afterok:<opt_job>)
    ↓  wait
    run_after_gaussian.sh  (antechamber RESP + parmchk2 + cap stripping)

Public API
----------
run_gaussian_parametrization_for_protein(protein_dir, pdb_id, slurm_config)
"""

from __future__ import annotations

import json
import shutil
import subprocess
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Any


# ---------------------------------------------------------------------------
# SLURM configuration
# ---------------------------------------------------------------------------

@dataclass
class SlurmConfig:
    partition: str = "compute"
    time: str = "24:00:00"
    mem: str = "32G"
    cpus: int = 16

    def resource_args(self) -> list[str]:
        return [
            f"--partition={self.partition}",
            "--nodes=1",
            "--ntasks=1",
            f"--cpus-per-task={self.cpus}",
            f"--mem={self.mem}",
            f"--time={self.time}",
            "--requeue",
        ]


# ---------------------------------------------------------------------------
# SLURM helpers
# ---------------------------------------------------------------------------

def _sbatch_available() -> bool:
    return shutil.which("sbatch") is not None


def _find_central_worker() -> Path | None:
    """Locate run_gaussian.sh in scripts/CESGA_SLURM/ relative to this package."""
    here = Path(__file__).resolve().parent          # .../src/stack_protein_preparation/
    repo_root = here.parent.parent                  # .../stack_protein_prep/
    candidate = repo_root / "scripts" / "CESGA_SLURM" / "run_gaussian.sh"
    return candidate if candidate.is_file() else None


def _sbatch_single(
    com_path: Path,
    worker: Path,
    slurm_config: SlurmConfig,
    log_dir: Path,
    job_name: str,
    dependency: str | None = None,
) -> str:
    """Submit one Gaussian .com file. Returns SLURM job ID string."""
    cmd = ["sbatch", "--parsable"]
    cmd += slurm_config.resource_args()
    cmd += [
        f"--job-name={job_name}",
        f"--output={log_dir}/{job_name}_%j.out",
        f"--error={log_dir}/{job_name}_%j.err",
    ]
    if dependency:
        cmd.append(f"--dependency={dependency}")
    cmd += ["--wrap", f"bash {worker} {com_path}"]
    result = subprocess.run(cmd, capture_output=True, text=True, check=True)
    return result.stdout.strip()


def _wait_for_slurm_jobs(
    job_ids: list[str],
    poll_interval_s: int = 60,
) -> dict[str, str]:
    """Block until all job_ids leave squeue. Returns {job_id: final_sacct_state}.

    Uses squeue polling rather than sbatch --wait so failures in any job
    do not silently block the caller.
    """
    pending = set(job_ids)
    while pending:
        result = subprocess.run(
            ["squeue", "--noheader", "--format=%i", "-j", ",".join(pending)],
            capture_output=True, text=True,
        )
        active = {line.strip() for line in result.stdout.splitlines() if line.strip()}
        pending &= active
        if pending:
            time.sleep(poll_interval_s)

    # sacct for final states
    sacct = subprocess.run(
        ["sacct", "-j", ",".join(job_ids),
         "--format=JobID,State", "--noheader", "--parsable2"],
        capture_output=True, text=True,
    )
    final: dict[str, str] = {}
    for line in sacct.stdout.splitlines():
        parts = line.split("|")
        if len(parts) < 2:
            continue
        base_id = parts[0].strip().split(".")[0]  # strip .batch / .extern
        state = parts[1].strip()
        if base_id in job_ids:
            if base_id not in final or state == "COMPLETED":
                final[base_id] = state
    return final


def _gaussian_log_ok(log_path: Path) -> bool:
    """Return True if log_path exists and contains Normal termination."""
    if not log_path.is_file():
        return False
    try:
        return "Normal termination of Gaussian" in log_path.read_text(errors="replace")
    except Exception:
        return False


def _run_sh(script_name: str, work_dir: Path, timeout_s: int = 900) -> None:
    """Run a bash script in work_dir; raise CalledProcessError on failure."""
    script = work_dir / script_name
    if not script.is_file():
        raise FileNotFoundError(f"{script_name} not found in {work_dir}")
    result = subprocess.run(
        ["bash", script_name],
        cwd=work_dir,
        capture_output=True,
        text=True,
        timeout=timeout_s,
    )
    if result.returncode != 0:
        raise subprocess.CalledProcessError(
            result.returncode, script_name,
            output=result.stdout, stderr=result.stderr,
        )


# ---------------------------------------------------------------------------
# Group-name extraction helpers
# ---------------------------------------------------------------------------

def _read_mcpb_group_name(mcpb_dir: Path, pdb_id: str) -> str | None:
    """Read group_name from {pdb_id}.in inside mcpb_dir."""
    in_file = mcpb_dir / f"{pdb_id}.in"
    if not in_file.is_file():
        candidates = list(mcpb_dir.glob("*.in"))
        if not candidates:
            return None
        in_file = candidates[0]
    for line in in_file.read_text(errors="replace").splitlines():
        stripped = line.strip()
        if stripped.startswith("group_name") and "=" in stripped:
            return stripped.split("=", 1)[1].strip()
    return None


def _resname_from_resp_label(label: str) -> str:
    """Extract residue name from directory label like 'A_PTR_488' → 'PTR'."""
    parts = label.split("_")
    return parts[1] if len(parts) >= 3 else parts[-1]


# ---------------------------------------------------------------------------
# MCPB Gaussian workflow
# ---------------------------------------------------------------------------

def _find_mcpb_site_dirs(protein_dir: Path) -> list[Path]:
    """Return sorted list of 03_mcpb/ dirs that have a commands.sh."""
    metall_dir = protein_dir / "metall_params"
    if not metall_dir.is_dir():
        return []
    result = []
    for site_dir in sorted(metall_dir.iterdir()):
        if not site_dir.is_dir():
            continue
        mcpb_dir = site_dir / "03_mcpb"
        if (mcpb_dir / "commands.sh").is_file():
            result.append(mcpb_dir)
    return result


def run_mcpb_gaussian_for_protein(
    protein_dir: Path,
    pdb_id: str,
    slurm_config: SlurmConfig | None = None,
    poll_interval_s: int = 60,
) -> dict[str, Any]:
    """Run the full MCPB Gaussian workflow for all transition-metal sites.

    Steps per site:
    1. bash commands.sh   (MCPB.py -s 1  — generates .com files)
    2. sbatch 3 jobs      (small_opt, small_fc, large_mk — parallel)
    3. squeue poll        (wait for all 3)
    4. bash run_after_gaussian.sh  (formchk + MCPB.py -s 2,3,4 + tleap)

    All sites' Gaussian jobs are submitted before waiting, so parallel HPC
    execution is maximised across sites of the same protein.
    """
    if slurm_config is None:
        slurm_config = SlurmConfig()

    result: dict[str, Any] = {
        "pdb_id": pdb_id,
        "type": "mcpb",
        "status": "skipped",
        "message": "",
        "n_sites": 0,
        "n_success": 0,
        "n_failed": 0,
        "site_results": [],
    }

    site_dirs = _find_mcpb_site_dirs(protein_dir)
    result["n_sites"] = len(site_dirs)
    if not site_dirs:
        result["message"] = "No MCPB sites with commands.sh found."
        return result

    hpc_available = _sbatch_available()

    if not hpc_available:
        # Local machine — run commands.sh to generate .com files, then stop.
        for mcpb_dir in site_dirs:
            step02 = mcpb_dir / "step02_gaussian"
            if step02.is_dir() and list(step02.glob("*.com")):
                continue  # already generated
            try:
                _run_sh("commands.sh", mcpb_dir)
            except Exception:
                pass  # best-effort; MCPB.py may not be installed locally
        result["status"] = "required"
        result["message"] = (
            f"{len(site_dirs)} MCPB site(s): Gaussian input files prepared. "
            "Submit on CESGA to complete parametrization."
        )
        return result

    worker = _find_central_worker()
    if worker is None:
        result["status"] = "failed"
        result["message"] = "scripts/CESGA_SLURM/run_gaussian.sh not found."
        return result

    log_dir = protein_dir / "metall_params" / "logs"
    log_dir.mkdir(parents=True, exist_ok=True)

    site_infos: list[dict] = []

    # ---- Phase 1: commands.sh (serial, fast) --------------------------------
    for mcpb_dir in site_dirs:
        group_name = _read_mcpb_group_name(mcpb_dir, pdb_id) or mcpb_dir.parent.name
        step02 = mcpb_dir / "step02_gaussian"

        site_info: dict[str, Any] = {
            "mcpb_dir": mcpb_dir,
            "group_name": group_name,
            "job_ids": [],
            "skip": False,
            "error": None,
        }

        # Skip if already done
        frcmod = mcpb_dir / "step03_amber_params" / f"{group_name}.frcmod"
        if frcmod.is_file():
            site_info["skip"] = True
            site_info["error"] = None
            result["site_results"].append({
                "group_name": group_name,
                "status": "already_done",
                "message": "frcmod already exists.",
            })
            site_infos.append(site_info)
            continue

        if not step02.is_dir() or not list(step02.glob("*.com")):
            try:
                _run_sh("commands.sh", mcpb_dir)
            except Exception as exc:
                stderr = getattr(exc, "stderr", None) or ""
                site_info["error"] = f"commands.sh rc={getattr(exc, 'returncode', '?')}: {stderr[-800:] or exc!r}"
                result["site_results"].append({
                    "group_name": group_name,
                    "status": "failed",
                    "message": site_info["error"],
                })
                result["n_failed"] += 1
                site_infos.append(site_info)
                continue

        site_infos.append(site_info)

    # ---- Phase 2: submit ALL Gaussian jobs ----------------------------------
    all_job_ids: list[str] = []
    for info in site_infos:
        if info["skip"] or info["error"]:
            continue
        group_name = info["group_name"]
        step02 = info["mcpb_dir"] / "step02_gaussian"
        for suffix in ("small_opt", "small_fc", "large_mk"):
            com = step02 / f"{group_name}_{suffix}.com"
            log = step02 / f"{group_name}_{suffix}.log"
            if not com.is_file() or _gaussian_log_ok(log):
                continue
            try:
                job_id = _sbatch_single(
                    com_path=com,
                    worker=worker,
                    slurm_config=slurm_config,
                    log_dir=log_dir,
                    job_name=f"{pdb_id}_{suffix}",
                )
                info["job_ids"].append((suffix, job_id))
                all_job_ids.append(job_id)
            except subprocess.CalledProcessError as exc:
                info["error"] = f"sbatch {suffix}: {exc.stderr}"
                result["n_failed"] += 1
                break

    # ---- Phase 3: wait ------------------------------------------------------
    if all_job_ids:
        _wait_for_slurm_jobs(all_job_ids, poll_interval_s=poll_interval_s)

    # ---- Phase 4: post-processing per site ----------------------------------
    for info in site_infos:
        if info["skip"]:
            continue
        if info["error"]:
            continue  # already counted as failed above

        group_name = info["group_name"]
        mcpb_dir = info["mcpb_dir"]
        step02 = mcpb_dir / "step02_gaussian"

        # Verify Gaussian logs
        failures = [
            suf for suf in ("small_opt", "small_fc", "large_mk")
            if not _gaussian_log_ok(step02 / f"{group_name}_{suf}.log")
        ]
        if failures:
            result["site_results"].append({
                "group_name": group_name,
                "status": "failed",
                "message": f"Gaussian failed/missing: {failures}",
            })
            result["n_failed"] += 1
            continue

        # Post-processing
        try:
            _run_sh("run_after_gaussian.sh", mcpb_dir, timeout_s=1800)
        except Exception as exc:
            msg = str(exc)
            if hasattr(exc, "stderr"):
                msg = (exc.stderr or "")[-600:]
            result["site_results"].append({
                "group_name": group_name,
                "status": "failed",
                "message": f"run_after_gaussian.sh: {msg}",
            })
            result["n_failed"] += 1
            continue

        frcmod = mcpb_dir / "step03_amber_params" / f"{group_name}.frcmod"
        result["site_results"].append({
            "group_name": group_name,
            "status": "success",
            "message": "MCPB parametrization complete.",
            "frcmod_path": str(frcmod) if frcmod.is_file() else "",
        })
        result["n_success"] += 1

    n_s, n_f = result["n_success"], result["n_failed"]
    n_total = result["n_sites"]
    if n_f == 0 and n_s > 0:
        result["status"] = "success"
    elif n_f == 0 and n_s == 0:
        result["status"] = "skipped"
    elif n_s > 0:
        result["status"] = "partial"
    else:
        result["status"] = "failed"
    result["message"] = f"{n_s}/{n_total} MCPB site(s) complete; {n_f} failed."
    return result


# ---------------------------------------------------------------------------
# RESP Gaussian workflow
# ---------------------------------------------------------------------------

def _find_resp_dirs(protein_dir: Path) -> list[Path]:
    """Return sorted list of resp/ dirs that have commands.sh."""
    nonstd_dir = protein_dir / "nonstandard_params"
    if not nonstd_dir.is_dir():
        return []
    result = []
    for label_dir in sorted(nonstd_dir.iterdir()):
        if not label_dir.is_dir():
            continue
        resp_dir = label_dir / "resp"
        if (resp_dir / "commands.sh").is_file():
            result.append(resp_dir)
    return result


def run_resp_gaussian_for_protein(
    protein_dir: Path,
    pdb_id: str,
    slurm_config: SlurmConfig | None = None,
    poll_interval_s: int = 60,
) -> dict[str, Any]:
    """Run the full RESP Gaussian workflow for all non-standard residues.

    Steps per residue:
    1. bash commands.sh   (antechamber  —  generates opt.com + esp.com)
    2. sbatch opt job
    3. sbatch esp job     (with --dependency=afterok:<opt_job_id>)
    4. squeue poll        (wait for esp, which implies opt done)
    5. bash run_after_gaussian.sh  (antechamber RESP + parmchk2 + cap strip)

    Jobs from all residues are submitted before the first wait, maximising
    parallel HPC utilisation.
    """
    if slurm_config is None:
        slurm_config = SlurmConfig()

    result: dict[str, Any] = {
        "pdb_id": pdb_id,
        "type": "resp",
        "status": "skipped",
        "message": "",
        "n_residues": 0,
        "n_success": 0,
        "n_failed": 0,
        "residue_results": [],
    }

    resp_dirs = _find_resp_dirs(protein_dir)
    result["n_residues"] = len(resp_dirs)
    if not resp_dirs:
        result["message"] = "No RESP directories with commands.sh found."
        return result

    hpc_available = _sbatch_available()

    if not hpc_available:
        # Local machine — run commands.sh to generate .com files, then stop.
        for resp_dir in resp_dirs:
            step02 = resp_dir / "step02_gaussian"
            if step02.is_dir() and list(step02.glob("*.com")):
                continue  # already generated
            try:
                _run_sh("commands.sh", resp_dir)
            except Exception:
                pass  # best-effort; antechamber may not be installed locally
        result["status"] = "required"
        result["message"] = (
            f"{len(resp_dirs)} RESP residue(s): Gaussian input files prepared. "
            "Submit on CESGA to complete parametrization."
        )
        return result

    worker = _find_central_worker()
    if worker is None:
        result["status"] = "failed"
        result["message"] = "scripts/CESGA_SLURM/run_gaussian.sh not found."
        return result

    log_dir = protein_dir / "nonstandard_params" / "logs"
    log_dir.mkdir(parents=True, exist_ok=True)

    residue_infos: list[dict] = []

    # ---- Phase 1: commands.sh + submit (interleaved; esp depends on opt) ----
    all_esp_job_ids: list[str] = []

    for resp_dir in resp_dirs:
        label = resp_dir.parent.name         # e.g. "A_PTR_488"
        resname = _resname_from_resp_label(label)
        step02 = resp_dir / "step02_gaussian"

        info: dict[str, Any] = {
            "resp_dir": resp_dir,
            "resname": resname,
            "label": label,
            "esp_job_id": None,
            "skip": False,
            "error": None,
        }

        # Skip if already done
        frcmod = resp_dir / "step03_resp_params" / f"{resname}.frcmod"
        if frcmod.is_file():
            info["skip"] = True
            result["residue_results"].append({
                "resname": resname,
                "status": "already_done",
                "message": "frcmod already exists.",
            })
            residue_infos.append(info)
            continue

        # Run commands.sh if needed
        if not step02.is_dir() or not list(step02.glob("*.com")):
            try:
                _run_sh("commands.sh", resp_dir)
            except Exception as exc:
                stderr = getattr(exc, "stderr", None) or ""
                info["error"] = f"commands.sh rc={getattr(exc, 'returncode', '?')}: {stderr[-800:] or exc!r}"
                result["residue_results"].append({
                    "resname": resname,
                    "status": "failed",
                    "message": info["error"],
                })
                result["n_failed"] += 1
                residue_infos.append(info)
                continue

        # Submit opt
        opt_com = step02 / f"{resname}_opt.com"
        opt_log = step02 / f"{resname}_opt.log"
        esp_com = step02 / f"{resname}_esp.com"
        esp_log = step02 / f"{resname}_esp.log"

        if not opt_com.is_file() or not esp_com.is_file():
            info["error"] = f"Expected .com files missing in {step02}"
            result["residue_results"].append({
                "resname": resname,
                "status": "failed",
                "message": info["error"],
            })
            result["n_failed"] += 1
            residue_infos.append(info)
            continue

        try:
            if _gaussian_log_ok(opt_log):
                opt_job_id = None  # already done
            else:
                opt_job_id = _sbatch_single(
                    com_path=opt_com,
                    worker=worker,
                    slurm_config=slurm_config,
                    log_dir=log_dir,
                    job_name=f"{pdb_id}_{resname}_opt",
                )

            if _gaussian_log_ok(esp_log):
                esp_job_id = None  # already done
            else:
                dep = f"afterok:{opt_job_id}" if opt_job_id else None
                esp_job_id = _sbatch_single(
                    com_path=esp_com,
                    worker=worker,
                    slurm_config=slurm_config,
                    log_dir=log_dir,
                    job_name=f"{pdb_id}_{resname}_esp",
                    dependency=dep,
                )

            if esp_job_id:
                info["esp_job_id"] = esp_job_id
                all_esp_job_ids.append(esp_job_id)
            if opt_job_id and not esp_job_id:
                # opt submitted but esp already done — wait for opt via its ID
                all_esp_job_ids.append(opt_job_id)

        except subprocess.CalledProcessError as exc:
            info["error"] = f"sbatch: {exc.stderr}"
            result["residue_results"].append({
                "resname": resname,
                "status": "failed",
                "message": info["error"],
            })
            result["n_failed"] += 1

        residue_infos.append(info)

    # ---- Phase 2: wait for all esp (and opt via dependency chain) -----------
    if all_esp_job_ids:
        _wait_for_slurm_jobs(all_esp_job_ids, poll_interval_s=poll_interval_s)

    # ---- Phase 3: post-processing per residue --------------------------------
    for info in residue_infos:
        if info["skip"] or info["error"]:
            continue

        resname = info["resname"]
        resp_dir = info["resp_dir"]
        step02 = resp_dir / "step02_gaussian"

        # Verify both Gaussian logs
        failures = [
            s for s in ("opt", "esp")
            if not _gaussian_log_ok(step02 / f"{resname}_{s}.log")
        ]
        if failures:
            result["residue_results"].append({
                "resname": resname,
                "status": "failed",
                "message": f"Gaussian failed/missing: {failures}",
            })
            result["n_failed"] += 1
            continue

        try:
            _run_sh("run_after_gaussian.sh", resp_dir, timeout_s=900)
        except Exception as exc:
            msg = str(exc)
            if hasattr(exc, "stderr"):
                msg = (exc.stderr or "")[-600:]
            result["residue_results"].append({
                "resname": resname,
                "status": "failed",
                "message": f"run_after_gaussian.sh: {msg}",
            })
            result["n_failed"] += 1
            continue

        frcmod = resp_dir / "step03_resp_params" / f"{resname}.frcmod"
        mol2 = resp_dir / "step03_resp_params" / f"{resname}.mol2"
        result["residue_results"].append({
            "resname": resname,
            "status": "success",
            "message": "RESP parametrization complete.",
            "frcmod_path": str(frcmod) if frcmod.is_file() else "",
            "mol2_path": str(mol2) if mol2.is_file() else "",
        })
        result["n_success"] += 1

    n_s, n_f = result["n_success"], result["n_failed"]
    n_total = result["n_residues"]
    if n_f == 0 and n_s > 0:
        result["status"] = "success"
    elif n_f == 0 and n_s == 0:
        result["status"] = "skipped"
    elif n_s > 0:
        result["status"] = "partial"
    else:
        result["status"] = "failed"
    result["message"] = f"{n_s}/{n_total} RESP residue(s) complete; {n_f} failed."
    return result


# ---------------------------------------------------------------------------
# Combined public entry point
# ---------------------------------------------------------------------------

def run_gaussian_parametrization_for_protein(
    protein_dir: Path,
    pdb_id: str,
    slurm_config: SlurmConfig | None = None,
    poll_interval_s: int = 60,
) -> dict[str, Any]:
    """Run all Gaussian parametrization workflows for one protein.

    Runs MCPB (metal centers) and RESP (non-standard residues) sequentially.
    Within each, all jobs are submitted before waiting so the HPC scheduler
    can parallelise them.

    Returns
    -------
    dict with keys: status, message, mcpb_result, resp_result, manifest_path.
    """
    result: dict[str, Any] = {
        "pdb_id": pdb_id,
        "status": "skipped",
        "message": "",
        "mcpb_result": None,
        "resp_result": None,
        "manifest_path": "",
    }

    mcpb_result = run_mcpb_gaussian_for_protein(
        protein_dir=protein_dir,
        pdb_id=pdb_id,
        slurm_config=slurm_config,
        poll_interval_s=poll_interval_s,
    )
    resp_result = run_resp_gaussian_for_protein(
        protein_dir=protein_dir,
        pdb_id=pdb_id,
        slurm_config=slurm_config,
        poll_interval_s=poll_interval_s,
    )

    result["mcpb_result"] = mcpb_result
    result["resp_result"] = resp_result

    statuses = {mcpb_result["status"], resp_result["status"]}
    if "failed" in statuses:
        result["status"] = "failed"
    elif "partial" in statuses:
        result["status"] = "partial"
    elif "required" in statuses or "pending_hpc" in statuses:
        result["status"] = "required"
    elif statuses == {"skipped"}:
        result["status"] = "skipped"
    elif "success" in statuses:
        result["status"] = "success"
    else:
        result["status"] = "skipped"

    result["message"] = (
        f"MCPB: {mcpb_result['message']} | RESP: {resp_result['message']}"
    )

    # Write manifest
    manifest_path = protein_dir / "gaussian_params_manifest.json"
    manifest_path.write_text(json.dumps(result, indent=2, default=str), encoding="utf-8")
    result["manifest_path"] = str(manifest_path)

    return result
