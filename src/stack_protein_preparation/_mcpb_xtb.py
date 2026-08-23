"""xtb (GFN2-xTB) wrapper for MCPB bonded parametrization.

Scaffold module — subprocess plumbing + output-file discovery + fail-open.
The scientifically-heavy Seminario force-constant extraction is left as
a well-marked follow-up (see ``TODO(seminario)`` markers below).  This
module gives the pipeline enough hooks that a follow-up refactor can
drop in the real physics without touching the orchestration layer.

Public API:

    run_xtb_optimize_hessian(cluster_xyz, charge, multiplicity, workdir) -> XtbRunResult
    emit_frcmod_scaffold(xtb_result, output_frcmod) -> Path
    xtb_binary_available() -> bool

Fail-open pattern (matches _uniprot_idr, _cofactor_params, _dssp_ss3):
    when the xtb binary is missing from PATH, ``run_xtb_optimize_hessian``
    returns ran=False with a helpful fallback_reason instead of raising.
    The dispatch layer downgrades that to ``manual`` or falls back to
    the nonbonded Li-Merz tier.

Reference: Bannwarth, Ehlert, Grimme (2019) JCTC 15:1652 — GFN2-xTB.
"""
from __future__ import annotations

import shutil
import subprocess
from dataclasses import dataclass, field
from pathlib import Path


DEFAULT_XTB_TIMEOUT_SECONDS = 900   # 15 min per cluster; scale up if needed
DEFAULT_ACCURACY = 0.5              # xtb --acc; 0.5 = balanced accuracy/speed
DEFAULT_METHOD = "gfn2"


def xtb_binary_available(xtb_bin: str | None = None) -> str | None:
    """Return the resolved xtb binary path, or None when not on PATH."""
    return shutil.which(xtb_bin or "xtb")


@dataclass
class XtbRunResult:
    ran: bool
    passed: bool
    fallback_reason: str | None = None
    # Optimised geometry (xyz file path once xtb finishes)
    opt_geometry_xyz: Path | None = None
    # Raw output files xtb produces — parsed by the Seminario follow-up
    hessian_path: Path | None = None
    charges_path: Path | None = None       # CM5 or Mulliken partial charges
    wbo_path: Path | None = None           # Wiberg bond orders
    xtb_log_path: Path | None = None
    workdir: Path | None = None
    exit_code: int | None = None
    diagnostics: list[str] = field(default_factory=list)

    def to_dict(self) -> dict:
        return {
            "ran": self.ran,
            "passed": self.passed,
            "fallback_reason": self.fallback_reason,
            "opt_geometry_xyz": str(self.opt_geometry_xyz) if self.opt_geometry_xyz else None,
            "hessian_path": str(self.hessian_path) if self.hessian_path else None,
            "charges_path": str(self.charges_path) if self.charges_path else None,
            "wbo_path": str(self.wbo_path) if self.wbo_path else None,
            "exit_code": self.exit_code,
            "diagnostics": list(self.diagnostics),
        }


def write_cluster_xyz(atoms: list[tuple[str, float, float, float]], output_path: Path) -> Path:
    """Write a cluster (element, x, y, z) list to a standard XYZ file."""
    output_path.parent.mkdir(parents=True, exist_ok=True)
    lines = [str(len(atoms)), "FRUTON MCPB cluster"]
    for elem, x, y, z in atoms:
        lines.append(f"{elem:<3s} {x:>14.6f} {y:>14.6f} {z:>14.6f}")
    output_path.write_text("\n".join(lines) + "\n")
    return output_path


def run_xtb_optimize_hessian(
    cluster_xyz: str | Path,
    charge: int,
    multiplicity: int,
    workdir: str | Path,
    method: str = DEFAULT_METHOD,
    accuracy: float = DEFAULT_ACCURACY,
    xtb_bin: str | None = None,
    timeout_seconds: int = DEFAULT_XTB_TIMEOUT_SECONDS,
) -> XtbRunResult:
    """Run xtb opt + hess on the cluster; return an XtbRunResult.

    Emits everything xtb needs into ``workdir`` and invokes:

        xtb <cluster.xyz> --opt --hess --gfn <method> --chrg <c> --uhf <n> --acc <a>

    ``multiplicity`` is 2S+1; xtb wants ``--uhf N`` where N = 2S (i.e.
    N = multiplicity − 1).
    """
    cluster_xyz = Path(cluster_xyz)
    workdir = Path(workdir)

    if not cluster_xyz.is_file():
        return XtbRunResult(
            ran=False, passed=False,
            fallback_reason=f"cluster xyz not found: {cluster_xyz}",
        )

    resolved = xtb_binary_available(xtb_bin)
    if resolved is None:
        return XtbRunResult(
            ran=False, passed=False,
            fallback_reason=(
                "xtb binary not on PATH; install via 'module load xtb/6.4.0' "
                "on CESGA or 'conda install -c conda-forge xtb'.  Falling "
                "back so the pipeline can continue (metal will be marked "
                "manual by the dispatch layer)."
            ),
        )

    workdir.mkdir(parents=True, exist_ok=True)
    log_path = workdir / "xtb.log"

    uhf = max(0, int(multiplicity) - 1)
    cmd = [
        resolved,
        str(cluster_xyz),
        "--opt",
        "--hess",
        "--gfn", "2" if method == "gfn2" else "1",
        "--chrg", str(int(charge)),
        "--uhf", str(uhf),
        "--acc", f"{accuracy:.3f}",
    ]

    diagnostics = [f"xtb_bin={resolved}", f"cmd={' '.join(cmd)}"]

    try:
        proc = subprocess.run(
            cmd,
            cwd=str(workdir),
            capture_output=True,
            text=True,
            timeout=timeout_seconds,
        )
    except subprocess.TimeoutExpired:
        return XtbRunResult(
            ran=True, passed=False, workdir=workdir,
            fallback_reason=f"xtb timed out after {timeout_seconds}s",
            xtb_log_path=log_path,
            diagnostics=diagnostics,
        )
    except Exception as exc:  # noqa: BLE001
        return XtbRunResult(
            ran=True, passed=False, workdir=workdir,
            fallback_reason=f"xtb subprocess failed: {exc!r}",
            diagnostics=diagnostics,
        )

    log_path.write_text((proc.stdout or "") + "\n---STDERR---\n" + (proc.stderr or ""))
    diagnostics.append(f"returncode={proc.returncode}")

    # xtb writes several files in the cwd it was called from.
    opt_xyz = workdir / "xtbopt.xyz"
    hess = workdir / "hessian"
    charges = workdir / "charges"
    wbo = workdir / "wbo"

    if proc.returncode != 0:
        return XtbRunResult(
            ran=True, passed=False, exit_code=proc.returncode,
            workdir=workdir,
            fallback_reason=f"xtb exited nonzero ({proc.returncode}); see log",
            xtb_log_path=log_path,
            opt_geometry_xyz=opt_xyz if opt_xyz.is_file() else None,
            diagnostics=diagnostics,
        )

    # Success path: check that at least the optimised geometry landed.
    if not opt_xyz.is_file():
        return XtbRunResult(
            ran=True, passed=False, exit_code=0, workdir=workdir,
            fallback_reason="xtb finished but xtbopt.xyz not written",
            xtb_log_path=log_path,
            diagnostics=diagnostics,
        )

    return XtbRunResult(
        ran=True,
        passed=True,
        exit_code=0,
        workdir=workdir,
        opt_geometry_xyz=opt_xyz,
        hessian_path=hess if hess.is_file() else None,
        charges_path=charges if charges.is_file() else None,
        wbo_path=wbo if wbo.is_file() else None,
        xtb_log_path=log_path,
        diagnostics=diagnostics,
    )


# ---------------------------------------------------------------------------
# Seminario force-constant extraction — TODO(seminario)
# ---------------------------------------------------------------------------


def emit_frcmod_scaffold(
    xtb_result: XtbRunResult,
    metal_element: str,
    metal_atom_type: str,
    donor_atom_types: list[str],
    output_frcmod: str | Path,
) -> Path:
    """Emit a MINIMAL well-formed frcmod stub from an xtb run.

    IMPORTANT — this is the scaffold that unblocks the rest of the
    pipeline (tleap will load this file).  The force constants written
    here are conservative reference values, NOT extracted from the
    xtb Hessian yet.

    TODO(seminario): replace the placeholder K_R / K_theta constants
    below with values derived from the xtb Hessian via the Seminario
    method (Seminario 1996 Int. J. Quantum Chem. 60:1271):

        K_R(A-B) = -u_AB^T · H_AB · u_AB
        where u_AB = unit vector A→B, H_AB = 3×3 Hessian block for (A,B)

    The Hessian file is at ``xtb_result.hessian_path`` (xtb writes it in
    a flat Fortran-format text file with the full 3N × 3N matrix).  A
    numpy loader + Seminario projector + frcmod re-emitter belongs here.

    Until then: this scaffold writes a frcmod that tleap will accept
    with generic parameters (K_R = 200 kcal/mol/Å², eq = 2.1 Å for
    M-N/O donors), so downstream tleap.in generation + sander sanity
    check don't fall over — the metal center is present, but a
    scientific review would flag the placeholder force constants.
    """
    output_frcmod = Path(output_frcmod)
    output_frcmod.parent.mkdir(parents=True, exist_ok=True)

    header = "FRUTON MCPB scaffold — xtb-based bonded metal parameters"
    body: list[str] = [header, "", "MASS"]

    # Placeholder atomic mass for the metal.  In a real emit, the
    # frcmod would look up the mass from the periodic table + xtb output.
    body.append(f"{metal_atom_type} 65.38  0.000  ! {metal_element} (placeholder mass)")
    body.append("")
    body.append("BOND")
    # TODO(seminario): replace K_R = 200.0, eq from xtb geometry
    for donor in donor_atom_types:
        body.append(
            f"{metal_atom_type}-{donor}   200.0    2.100    ! "
            f"placeholder — TODO(seminario) Seminario from xtb Hessian"
        )
    body.append("")
    body.append("ANGLE")
    # TODO(seminario): replace K_theta = 50.0, eq from xtb geometry per donor pair
    if len(donor_atom_types) >= 2:
        seen: set[tuple[str, str]] = set()
        for i, a in enumerate(donor_atom_types):
            for b in donor_atom_types[i + 1:]:
                pair = tuple(sorted((a, b)))
                if pair in seen:
                    continue
                seen.add(pair)
                body.append(
                    f"{a}-{metal_atom_type}-{b}   50.0   109.5   ! "
                    f"placeholder — TODO(seminario) Seminario from xtb Hessian"
                )
    body.append("")
    body.append("DIHE")
    body.append("")
    body.append("IMPROPER")
    body.append("")
    body.append("NONBON")
    # TODO(seminario): replace with LJ params derived from xtb (or fall
    # back to Li-Merz 12-6-4 entry from data/ion_12_6_4_reference.csv).
    body.append(
        f"  {metal_atom_type}          1.100  0.0125  ! placeholder LJ "
        f"— TODO(seminario) replace with LiMerz or xtb-derived"
    )
    body.append("")

    output_frcmod.write_text("\n".join(body) + "\n")
    return output_frcmod
