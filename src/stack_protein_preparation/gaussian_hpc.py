"""
gaussian_hpc.py — SLURM-based Gaussian submission and waiting for FRUTON pipeline.

Two-pass workflow across ALL proteins simultaneously:

    Pass 1  commands.sh for every protein (MCPB + RESP) — fast, login-node
    Pass 2  Submit ALL Gaussian jobs across ALL proteins at once
    Wait    Single squeue poll until every job is done
    Pass 3  run_after_gaussian.sh for every protein — produces frcmods

MCPB dependency chain per site:
    small_opt  ──┐ (parallel)
    large_mk     │
    small_fc  ←──┘ afterok:small_opt

RESP dependency chain per residue:
    opt ──→ esp  (afterok:opt)

Public API
----------
run_gaussian_parametrization_for_all_proteins(proteins, slurm_config)
run_gaussian_parametrization_for_protein(protein_dir, pdb_id, slurm_config)
"""

from __future__ import annotations

import json
import os
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
    partition: str = "short"
    time: str = "06:00:00"
    mem: str = "32G"
    cpus: int = 16
    gpus: int = 0

    def resource_args(self) -> list[str]:
        args = [
            f"--partition={self.partition}",
            "--nodes=1",
            "--ntasks=1",
            f"--cpus-per-task={self.cpus}",
            f"--mem={self.mem}",
            f"--time={self.time}",
        ]
        if self.gpus > 0:
            args.append(f"--gres=gpu:a100:{self.gpus}")
        return args


# ---------------------------------------------------------------------------
# SLURM helpers
# ---------------------------------------------------------------------------

def _sbatch_available() -> bool:
    return shutil.which("sbatch") is not None


def _find_central_worker() -> Path | None:
    """Locate run_gaussian.sh.

    Search order:
    1. $FRUTON_SCRIPTS env var  (set this in .bashrc on CESGA when working
       from $LUSTRE where the git repo is not present)
    2. Relative to this package (covers repo-local installs on $STORE)
    """
    env_scripts = os.environ.get("FRUTON_SCRIPTS", "").strip()
    if env_scripts:
        candidate = Path(env_scripts) / "run_gaussian.sh"
        if candidate.is_file():
            return candidate

    here = Path(__file__).resolve().parent
    repo_root = here.parent.parent
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
    """Block until all job_ids leave squeue. Returns {job_id: final_sacct_state}."""
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
        base_id = parts[0].strip().split(".")[0]
        state = parts[1].strip()
        if base_id in job_ids:
            if base_id not in final or state == "COMPLETED":
                final[base_id] = state
    return final


def _gaussian_log_ok(log_path: Path) -> bool:
    if not log_path.is_file():
        return False
    try:
        return "Normal termination of Gaussian" in log_path.read_text(errors="replace")
    except Exception:
        return False


def _run_sh(script_name: str, work_dir: Path, timeout_s: int = 900) -> None:
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
# Discovery helpers
# ---------------------------------------------------------------------------

def _find_mcpb_site_dirs(protein_dir: Path) -> list[Path]:
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


def _find_resp_dirs(protein_dir: Path) -> list[Path]:
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


def _read_mcpb_group_name(mcpb_dir: Path, pdb_id: str) -> str | None:
    in_file = mcpb_dir / f"{pdb_id}.in"
    if not in_file.is_file():
        candidates = list(mcpb_dir.glob("*.in"))
        if not candidates:
            return None
        in_file = candidates[0]
    for line in in_file.read_text(errors="replace").splitlines():
        stripped = line.strip()
        if stripped.startswith("group_name"):
            parts = stripped.split("=", 1) if "=" in stripped else stripped.split(None, 1)
            return parts[1].strip() if len(parts) > 1 else None
    return None


def _read_mcpb_charge_spin(mcpb_dir: Path, pdb_id: str) -> tuple[int | None, int | None]:
    """Parse charge and spin from charge_m_sm line in the MCPB.py .in file."""
    in_file = mcpb_dir / f"{pdb_id}.in"
    if not in_file.is_file():
        candidates = list(mcpb_dir.glob("*.in"))
        if not candidates:
            return None, None
        in_file = candidates[0]
    for line in in_file.read_text(errors="replace").splitlines():
        stripped = line.strip()
        if stripped.startswith("charge_m_sm"):
            parts = stripped.split()
            try:
                charge = int(parts[1]) if parts[1] != "TODO_charge" else None
                spin = int(parts[2]) if len(parts) > 2 and parts[2] != "TODO_spin" else None
                return charge, spin
            except (IndexError, ValueError):
                return None, None
    return None, None


def _read_resp_charge(resp_dir: Path, resname: str) -> int | None:
    """Parse the formal charge for a RESP residue.

    Tries (in order):
    1. nonstd_params_manifest.json two levels up (authoritative, always -formal charge)
    2. CHARGE= bash variable in commands.sh (current template)
    3. -nc X in antechamber line in commands.sh (old template)
    """
    import re as _re

    # 1. nonstd_params_manifest.json
    manifest = resp_dir.parent.parent / "nonstd_params_manifest.json"
    if manifest.is_file():
        try:
            data = json.loads(manifest.read_text(errors="replace"))
            for res in data.get("residues", []):
                if res.get("resname", "").upper() == resname.upper():
                    ch = res.get("formal_charge")
                    if ch is not None:
                        return int(ch)
        except Exception:
            pass

    # 2 & 3. commands.sh
    cmd_sh = resp_dir / "commands.sh"
    if not cmd_sh.is_file():
        return None
    text = cmd_sh.read_text(errors="replace")
    for line in text.splitlines():
        stripped = line.strip()
        # Current template: bare CHARGE=-2 (before parity script substitutes it)
        if _re.match(r'^CHARGE=(-?\d+)$', stripped):
            try:
                return int(stripped.split("=", 1)[1])
            except ValueError:
                pass
        # Old template: antechamber -nc -1
        m = _re.search(r'\s-nc\s+(-?\d+)', stripped)
        if m:
            try:
                return int(m.group(1))
            except ValueError:
                pass
    return None


def _resname_from_resp_label(label: str) -> str:
    parts = label.split("_")
    return parts[1] if len(parts) >= 3 else parts[-1]


# ---------------------------------------------------------------------------
# Phase 1: prepare inputs (commands.sh)
# ---------------------------------------------------------------------------

def _prepare_mcpb_sites(protein_dir: Path, pdb_id: str) -> tuple[list[dict], Path]:
    """Run commands.sh for all MCPB sites. Returns (site_infos, log_dir)."""
    site_dirs = _find_mcpb_site_dirs(protein_dir)
    log_dir = protein_dir / "metall_params" / "logs"
    log_dir.mkdir(parents=True, exist_ok=True)
    site_infos = []
    for mcpb_dir in site_dirs:
        group_name = _read_mcpb_group_name(mcpb_dir, pdb_id) or mcpb_dir.parent.name
        charge, spin = _read_mcpb_charge_spin(mcpb_dir, pdb_id)
        info: dict[str, Any] = {
            "mcpb_dir": mcpb_dir,
            "group_name": group_name,
            "pdb_id": pdb_id,
            "log_dir": log_dir,
            "job_ids": [],
            "skip": False,
            "error": None,
            "metal_charge": charge,
            "metal_spin": spin,
        }
        frcmod = mcpb_dir / "step03_amber_params" / f"{group_name}.frcmod"
        if frcmod.is_file():
            info["skip"] = True
            site_infos.append(info)
            continue
        step02 = mcpb_dir / "step02_gaussian"
        if not step02.is_dir() or not list(step02.glob("*.com")):
            try:
                _run_sh("commands.sh", mcpb_dir)
            except Exception as exc:
                stderr = getattr(exc, "stderr", None) or ""
                info["error"] = f"commands.sh rc={getattr(exc, 'returncode', '?')}: {stderr[-800:] or exc!r}"
        site_infos.append(info)
    return site_infos, log_dir


def _prepare_resp_residues(protein_dir: Path, pdb_id: str) -> tuple[list[dict], Path]:
    """Run commands.sh for all RESP residues. Returns (residue_infos, log_dir)."""
    resp_dirs = _find_resp_dirs(protein_dir)
    log_dir = protein_dir / "nonstandard_params" / "logs"
    log_dir.mkdir(parents=True, exist_ok=True)
    residue_infos = []
    for resp_dir in resp_dirs:
        label = resp_dir.parent.name
        resname = _resname_from_resp_label(label)
        info: dict[str, Any] = {
            "resp_dir": resp_dir,
            "resname": resname,
            "label": label,
            "pdb_id": pdb_id,
            "log_dir": log_dir,
            "esp_job_id": None,
            "skip": False,
            "error": None,
            "resp_charge": _read_resp_charge(resp_dir, resname),
        }
        frcmod = resp_dir / "step03_resp_params" / f"{resname}.frcmod"
        if frcmod.is_file():
            info["skip"] = True
            residue_infos.append(info)
            continue
        step02 = resp_dir / "step02_gaussian"
        if not step02.is_dir() or not list(step02.glob("*.com")):
            try:
                _run_sh("commands.sh", resp_dir)
            except Exception as exc:
                stderr = getattr(exc, "stderr", None) or ""
                info["error"] = f"commands.sh rc={getattr(exc, 'returncode', '?')}: {stderr[-800:] or exc!r}"
        residue_infos.append(info)
    return residue_infos, log_dir


# ---------------------------------------------------------------------------
# Phase 2: submit jobs
# ---------------------------------------------------------------------------

def _submit_mcpb_jobs(
    site_infos: list[dict],
    worker: Path,
    slurm_config: SlurmConfig,
) -> list[str]:
    """Submit Gaussian jobs for all MCPB sites. Returns submitted job IDs."""
    all_job_ids = []
    for info in site_infos:
        if info["skip"] or info["error"]:
            continue
        group_name = info["group_name"]
        pdb_id = info["pdb_id"]
        step02 = info["mcpb_dir"] / "step02_gaussian"
        opt_job_id: str | None = None
        for suffix in ("small_opt", "large_mk", "small_fc"):
            com = step02 / f"{group_name}_{suffix}.com"
            log = step02 / f"{group_name}_{suffix}.log"
            if not com.is_file() or _gaussian_log_ok(log):
                continue
            dep = f"afterok:{opt_job_id}" if suffix == "small_fc" and opt_job_id else None
            try:
                job_id = _sbatch_single(
                    com, worker, slurm_config, info["log_dir"],
                    f"{pdb_id}_{suffix}", dep,
                )
                if suffix == "small_opt":
                    opt_job_id = job_id
                info["job_ids"].append((suffix, job_id))
                all_job_ids.append(job_id)
            except subprocess.CalledProcessError as exc:
                info["error"] = f"sbatch {suffix}: {exc.stderr}"
                break
    return all_job_ids


def _submit_resp_jobs(
    residue_infos: list[dict],
    worker: Path,
    slurm_config: SlurmConfig,
) -> list[str]:
    """Submit Gaussian jobs for all RESP residues. Returns submitted job IDs."""
    all_job_ids = []
    for info in residue_infos:
        if info["skip"] or info["error"]:
            continue
        resname = info["resname"]
        pdb_id = info["pdb_id"]
        resp_dir = info["resp_dir"]
        step02 = resp_dir / "step02_gaussian"
        opt_com = step02 / f"{resname}_opt.com"
        opt_log = step02 / f"{resname}_opt.log"
        esp_com = step02 / f"{resname}_esp.com"
        esp_log = step02 / f"{resname}_esp.log"
        if not opt_com.is_file() or not esp_com.is_file():
            info["error"] = f"Expected .com files missing in {step02}"
            continue
        try:
            if _gaussian_log_ok(opt_log):
                opt_job_id = None
            else:
                opt_job_id = _sbatch_single(
                    opt_com, worker, slurm_config, info["log_dir"],
                    f"{pdb_id}_{resname}_opt",
                )
            dep = f"afterok:{opt_job_id}" if opt_job_id else None
            if _gaussian_log_ok(esp_log):
                esp_job_id = None
            else:
                esp_job_id = _sbatch_single(
                    esp_com, worker, slurm_config, info["log_dir"],
                    f"{pdb_id}_{resname}_esp", dep,
                )
            info["esp_job_id"] = esp_job_id
            if esp_job_id:
                all_job_ids.append(esp_job_id)
            elif opt_job_id:
                all_job_ids.append(opt_job_id)
        except subprocess.CalledProcessError as exc:
            info["error"] = f"sbatch: {exc.stderr}"
    return all_job_ids


# ---------------------------------------------------------------------------
# Phase 3: post-processing (run_after_gaussian.sh)
# ---------------------------------------------------------------------------

def _postprocess_mcpb_sites(site_infos: list[dict]) -> dict[str, Any]:
    """Run run_after_gaussian.sh for all MCPB sites. Returns result dict."""
    n_success = n_failed = 0
    site_results = []
    for info in site_infos:
        group_name = info["group_name"]
        if info["skip"]:
            frcmod = info["mcpb_dir"] / "step03_amber_params" / f"{group_name}.frcmod"
            site_results.append({
                "group_name": group_name, "status": "already_done",
                "message": "frcmod already exists.",
                "frcmod_path": str(frcmod) if frcmod.is_file() else "",
            })
            n_success += 1
            continue
        if info["error"]:
            site_results.append({
                "group_name": group_name, "status": "failed",
                "message": info["error"],
            })
            n_failed += 1
            continue
        mcpb_dir = info["mcpb_dir"]
        step02 = mcpb_dir / "step02_gaussian"
        failures = [
            suf for suf in ("small_opt", "small_fc", "large_mk")
            if not _gaussian_log_ok(step02 / f"{group_name}_{suf}.log")
        ]
        if failures:
            site_results.append({
                "group_name": group_name, "status": "failed",
                "message": f"Gaussian failed/missing: {failures}",
            })
            n_failed += 1
            continue
        try:
            _run_sh("run_after_gaussian.sh", mcpb_dir, timeout_s=1800)
        except Exception as exc:
            msg = (getattr(exc, "stderr", None) or "")[-600:] or str(exc)
            site_results.append({
                "group_name": group_name, "status": "failed",
                "message": f"run_after_gaussian.sh: {msg}",
            })
            n_failed += 1
            continue
        frcmod = mcpb_dir / "step03_amber_params" / f"{group_name}.frcmod"
        site_results.append({
            "group_name": group_name, "status": "success",
            "message": "MCPB parametrization complete.",
            "frcmod_path": str(frcmod) if frcmod.is_file() else "",
        })
        n_success += 1

    n_total = len(site_infos)
    if n_failed == 0 and n_success > 0:
        status = "success"
    elif n_failed == 0 and n_success == 0:
        status = "skipped"
    elif n_success > 0:
        status = "partial"
    else:
        status = "failed"
    return {
        "type": "mcpb",
        "n_sites": n_total,
        "n_success": n_success,
        "n_failed": n_failed,
        "site_results": site_results,
        "status": status,
        "message": f"{n_success}/{n_total} MCPB site(s) complete; {n_failed} failed.",
    }


def _postprocess_resp_residues(residue_infos: list[dict]) -> dict[str, Any]:
    """Run run_after_gaussian.sh for all RESP residues. Returns result dict."""
    n_success = n_failed = 0
    residue_results = []
    for info in residue_infos:
        resname = info["resname"]
        if info["skip"]:
            frcmod = info["resp_dir"] / "step03_resp_params" / f"{resname}.frcmod"
            mol2 = info["resp_dir"] / "step03_resp_params" / f"{resname}.mol2"
            residue_results.append({
                "resname": resname, "status": "already_done",
                "message": "frcmod already exists.",
                "frcmod_path": str(frcmod) if frcmod.is_file() else "",
                "mol2_path": str(mol2) if mol2.is_file() else "",
            })
            n_success += 1
            continue
        if info["error"]:
            residue_results.append({
                "resname": resname, "status": "failed",
                "message": info["error"],
            })
            n_failed += 1
            continue
        resp_dir = info["resp_dir"]
        step02 = resp_dir / "step02_gaussian"
        failures = [
            s for s in ("opt", "esp")
            if not _gaussian_log_ok(step02 / f"{resname}_{s}.log")
        ]
        if failures:
            residue_results.append({
                "resname": resname, "status": "failed",
                "message": f"Gaussian failed/missing: {failures}",
            })
            n_failed += 1
            continue
        try:
            _run_sh("run_after_gaussian.sh", resp_dir, timeout_s=900)
        except Exception as exc:
            msg = (getattr(exc, "stderr", None) or "")[-600:] or str(exc)
            residue_results.append({
                "resname": resname, "status": "failed",
                "message": f"run_after_gaussian.sh: {msg}",
            })
            n_failed += 1
            continue
        frcmod = resp_dir / "step03_resp_params" / f"{resname}.frcmod"
        mol2 = resp_dir / "step03_resp_params" / f"{resname}.mol2"
        residue_results.append({
            "resname": resname, "status": "success",
            "message": "RESP parametrization complete.",
            "frcmod_path": str(frcmod) if frcmod.is_file() else "",
            "mol2_path": str(mol2) if mol2.is_file() else "",
        })
        n_success += 1

    n_total = len(residue_infos)
    if n_failed == 0 and n_success > 0:
        status = "success"
    elif n_failed == 0 and n_success == 0:
        status = "skipped"
    elif n_success > 0:
        status = "partial"
    else:
        status = "failed"
    return {
        "type": "resp",
        "n_residues": n_total,
        "n_success": n_success,
        "n_failed": n_failed,
        "residue_results": residue_results,
        "status": status,
        "message": f"{n_success}/{n_total} RESP residue(s) complete; {n_failed} failed.",
    }


# ---------------------------------------------------------------------------
# Combined public entry point — all proteins at once
# ---------------------------------------------------------------------------

def run_gaussian_parametrization_for_all_proteins(
    proteins: list[tuple[Path, str]],
    slurm_config: SlurmConfig | None = None,
    poll_interval_s: int = 60,
) -> dict[str, dict[str, Any]]:
    """Two-pass Gaussian parametrization across ALL proteins simultaneously.

    Pass 1  Run commands.sh for every protein (MCPB + RESP) — fast, login-node.
    Pass 2  Submit ALL Gaussian jobs for all proteins in one batch.
    Wait    Single blocking poll until every job is done.
    Pass 3  Run run_after_gaussian.sh for every protein.

    Returns
    -------
    dict mapping pdb_id → result dict with keys:
        status, message, mcpb_result, resp_result, manifest_path.
    """
    if slurm_config is None:
        slurm_config = SlurmConfig()

    hpc_available = _sbatch_available()

    worker: Path | None = None
    if hpc_available:
        worker = _find_central_worker()
        if worker is None:
            return {
                pdb_id: {
                    "pdb_id": pdb_id, "status": "failed",
                    "message": "scripts/CESGA_SLURM/run_gaussian.sh not found.",
                    "mcpb_result": None, "resp_result": None, "manifest_path": "",
                }
                for _, pdb_id in proteins
            }

    # ---- Phase 1: commands.sh for all proteins (serial, fast) ---------------
    all_mcpb: dict[str, list[dict]] = {}
    all_resp: dict[str, list[dict]] = {}
    protein_dir_map: dict[str, Path] = {}

    for protein_dir, pdb_id in proteins:
        protein_dir_map[pdb_id] = protein_dir
        site_infos, _ = _prepare_mcpb_sites(protein_dir, pdb_id)
        residue_infos, _ = _prepare_resp_residues(protein_dir, pdb_id)
        all_mcpb[pdb_id] = site_infos
        all_resp[pdb_id] = residue_infos

    if not hpc_available:
        results = {}
        for protein_dir, pdb_id in proteins:
            site_infos = all_mcpb.get(pdb_id, [])
            residue_infos = all_resp.get(pdb_id, [])
            n_mcpb = len(site_infos)
            n_resp = len(residue_infos)
            mcpb_charges = [i["metal_charge"] for i in site_infos if i.get("metal_charge") is not None]
            mcpb_spins = [i["metal_spin"] for i in site_infos if i.get("metal_spin") is not None]
            resp_charges = [i["resp_charge"] for i in residue_infos if i.get("resp_charge") is not None]
            results[pdb_id] = {
                "pdb_id": pdb_id, "status": "required",
                "message": (
                    f"{n_mcpb} MCPB site(s), {n_resp} RESP residue(s): "
                    "Gaussian inputs prepared. Submit on CESGA to complete."
                ),
                "mcpb_result": {"status": "required", "n_sites": n_mcpb,
                                "message": "Gaussian input files prepared."},
                "resp_result": {"status": "required", "n_residues": n_resp,
                                "message": "Gaussian input files prepared."},
                "mcpb_metal_charge": mcpb_charges[0] if len(mcpb_charges) == 1 else (mcpb_charges or None),
                "mcpb_metal_spin": mcpb_spins[0] if len(mcpb_spins) == 1 else (mcpb_spins or None),
                "resp_charge": resp_charges[0] if len(resp_charges) == 1 else (resp_charges or None),
                "frcmod_paths": "",
                "manifest_path": "",
            }
        return results

    # ---- Phase 2: submit ALL jobs for ALL proteins --------------------------
    all_job_ids: list[str] = []
    for pdb_id, site_infos in all_mcpb.items():
        all_job_ids.extend(_submit_mcpb_jobs(site_infos, worker, slurm_config))
    for pdb_id, residue_infos in all_resp.items():
        all_job_ids.extend(_submit_resp_jobs(residue_infos, worker, slurm_config))

    # ---- Wait for all jobs --------------------------------------------------
    if all_job_ids:
        _wait_for_slurm_jobs(all_job_ids, poll_interval_s=poll_interval_s)

    # ---- Phase 3: post-processing for all proteins --------------------------
    results: dict[str, dict[str, Any]] = {}
    for pdb_id, protein_dir in protein_dir_map.items():
        mcpb_result = _postprocess_mcpb_sites(all_mcpb.get(pdb_id, []))
        resp_result = _postprocess_resp_residues(all_resp.get(pdb_id, []))

        statuses = {mcpb_result["status"], resp_result["status"]}
        if "failed" in statuses:
            status = "failed"
        elif "partial" in statuses:
            status = "partial"
        elif "required" in statuses:
            status = "required"
        elif statuses == {"skipped"}:
            status = "skipped"
        elif "success" in statuses:
            status = "success"
        else:
            status = "skipped"

        # Collect charge/spin/frcmod metadata for pipeline record fields
        mcpb_charges = [
            info["metal_charge"] for info in all_mcpb.get(pdb_id, [])
            if info.get("metal_charge") is not None
        ]
        mcpb_spins = [
            info["metal_spin"] for info in all_mcpb.get(pdb_id, [])
            if info.get("metal_spin") is not None
        ]
        resp_charges = [
            info["resp_charge"] for info in all_resp.get(pdb_id, [])
            if info.get("resp_charge") is not None
        ]
        frcmod_paths = []
        for sr in mcpb_result.get("site_results", []):
            if sr.get("frcmod_path"):
                frcmod_paths.append(sr["frcmod_path"])
        for rr in resp_result.get("residue_results", []):
            if rr.get("frcmod_path"):
                frcmod_paths.append(rr["frcmod_path"])

        result: dict[str, Any] = {
            "pdb_id": pdb_id,
            "status": status,
            "message": f"MCPB: {mcpb_result['message']} | RESP: {resp_result['message']}",
            "mcpb_result": mcpb_result,
            "resp_result": resp_result,
            "mcpb_metal_charge": mcpb_charges[0] if len(mcpb_charges) == 1 else (mcpb_charges or None),
            "mcpb_metal_spin": mcpb_spins[0] if len(mcpb_spins) == 1 else (mcpb_spins or None),
            "resp_charge": resp_charges[0] if len(resp_charges) == 1 else (resp_charges or None),
            "frcmod_paths": ";".join(frcmod_paths),
            "manifest_path": "",
        }
        manifest_path = protein_dir / "gaussian_params_manifest.json"
        manifest_path.write_text(json.dumps(result, indent=2, default=str), encoding="utf-8")
        result["manifest_path"] = str(manifest_path)
        results[pdb_id] = result

    return results


def run_gaussian_parametrization_for_protein(
    protein_dir: Path,
    pdb_id: str,
    slurm_config: SlurmConfig | None = None,
    poll_interval_s: int = 60,
) -> dict[str, Any]:
    """Single-protein convenience wrapper — delegates to the all-proteins function."""
    results = run_gaussian_parametrization_for_all_proteins(
        proteins=[(protein_dir, pdb_id)],
        slurm_config=slurm_config,
        poll_interval_s=poll_interval_s,
    )
    return results.get(pdb_id, {
        "pdb_id": pdb_id, "status": "skipped", "message": "",
        "mcpb_result": None, "resp_result": None, "manifest_path": "",
    })
