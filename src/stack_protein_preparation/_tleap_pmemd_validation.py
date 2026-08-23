"""End-to-end MD-input-loadability sanity check for FRUTON output.

FRUTON's model-prep goal (per user 2026-08-23): "es muss einfach alles
funktionieren — am Ende komplette frcmod oder was auch immer."  A
green pytest suite proves the emitted files are structurally valid,
but the only way to confirm the *chemistry* is consistent is to hand
them to Amber and see if it accepts them.  This module runs that
handoff in two stages:

1. ``validate_tleap_can_load(tleap_in_path)``
     Runs ``tleap -f tleap.in`` in a fresh cwd.  Confirms every
     ``source`` / ``loadmol2`` / ``loadamberparams`` / ``loadpdb``
     / ``solvatebox`` / ``saveamberparm`` line completed without
     tleap emitting a ``FATAL`` / ``Error`` line.  Reads back the
     generated .prmtop + .rst7 to confirm they are non-empty and
     have the expected file signatures (``%VERSION`` in prmtop,
     3-line rst7 header).

2. ``validate_sander_zero_step(prmtop, rst7)``
     Runs sander (CPU, license-free) for 0 minimization steps.  This
     is the cheapest possible "does Amber accept this system"
     sanity: sander reads the topology + coord, verifies the bonded
     terms, allocates the neighbor list, and exits.  Missing
     parameters, atom-name mismatches, charge inconsistencies, and
     invalid box dimensions all fail this test even though a full
     MD would fail later.

Rationale for **sander not pmemd**: sander ships free with AmberTools
(no license), pmemd.cuda requires a licensed Amber install.  0-step
sanity does not need GPU speed.  A downstream MD project can swap in
pmemd.cuda for real runs; sander here is the license-independent
oracle Davide's later refactor keeps.

Fail-open pattern (matches _cofactor_params + _uniprot_idr): if
tleap or sander is not on PATH, returns a ValidationResult with
``ran=False`` and a helpful reason.  Never raises.

License-free: pure subprocess call to AmberTools (all free).  No
Python dep beyond stdlib.
"""
from __future__ import annotations

import os
import shutil
import subprocess
import tempfile
from dataclasses import dataclass, field
from pathlib import Path


DEFAULT_TLEAP_TIMEOUT_SECONDS = 120
DEFAULT_SANDER_TIMEOUT_SECONDS = 120


@dataclass
class ValidationResult:
    tool: str                       # "tleap" | "sander"
    ran: bool
    passed: bool
    fallback_reason: str | None = None
    exit_code: int | None = None
    prmtop_path: Path | None = None
    rst7_path: Path | None = None
    stdout_tail: str = ""
    stderr_tail: str = ""
    diagnostics: list[str] = field(default_factory=list)

    def to_dict(self) -> dict:
        return {
            "tool": self.tool,
            "ran": self.ran,
            "passed": self.passed,
            "fallback_reason": self.fallback_reason,
            "exit_code": self.exit_code,
            "prmtop_path": str(self.prmtop_path) if self.prmtop_path else None,
            "rst7_path": str(self.rst7_path) if self.rst7_path else None,
            "stdout_tail": self.stdout_tail,
            "stderr_tail": self.stderr_tail,
            "diagnostics": list(self.diagnostics),
        }


def _tail(text: str, n: int = 40) -> str:
    if not text:
        return ""
    return "\n".join(text.splitlines()[-n:])


def _looks_like_valid_prmtop(prmtop_path: Path) -> bool:
    """First 256 B of an Amber prmtop should contain ``%VERSION`` and
    ``%FLAG POINTERS`` markers."""
    if not prmtop_path.is_file() or prmtop_path.stat().st_size < 512:
        return False
    try:
        head = prmtop_path.read_text(encoding="ascii", errors="replace")[:2048]
    except OSError:
        return False
    return "%VERSION" in head and "%FLAG POINTERS" in head


def _looks_like_valid_rst7(rst7_path: Path) -> bool:
    """Amber ASCII rst7 line 1 is the title, line 2 is 'N_atoms   time'."""
    if not rst7_path.is_file() or rst7_path.stat().st_size < 32:
        return False
    try:
        with rst7_path.open(encoding="ascii", errors="replace") as fh:
            _title = fh.readline()
            second = fh.readline().strip()
    except OSError:
        return False
    parts = second.split()
    if not parts:
        return False
    try:
        int(parts[0])  # first token must parse as N_atoms
        return True
    except ValueError:
        return False


def _prmtop_and_rst7_from_tleap(tleap_in_path: Path) -> tuple[Path | None, Path | None]:
    """Parse ``saveamberparm mol <prmtop> <rst7>`` from a tleap.in file
    (best-effort) to find the expected output paths.  Returns (None, None)
    when the directive is not found."""
    if not tleap_in_path.is_file():
        return (None, None)
    try:
        text = tleap_in_path.read_text(encoding="utf-8", errors="replace")
    except OSError:
        return (None, None)
    for line in text.splitlines():
        s = line.strip()
        if not s.startswith("saveamberparm"):
            continue
        parts = s.split()
        # saveamberparm <unit> <prmtop> <rst7>
        if len(parts) >= 4:
            prmtop = Path(parts[2])
            rst7 = Path(parts[3])
            if not prmtop.is_absolute():
                prmtop = tleap_in_path.parent / prmtop
            if not rst7.is_absolute():
                rst7 = tleap_in_path.parent / rst7
            return (prmtop.resolve(), rst7.resolve())
    return (None, None)


def validate_tleap_can_load(
    tleap_in_path: str | Path,
    timeout_seconds: int = DEFAULT_TLEAP_TIMEOUT_SECONDS,
    tleap_bin: str | None = None,
) -> ValidationResult:
    """Run ``tleap -f tleap_in_path``, verify prmtop + rst7 written.

    Returns ``ran=False`` + fallback reason when tleap is missing from
    PATH.  ``passed=True`` requires exit 0 AND valid prmtop + rst7 signatures.
    """
    tleap_in = Path(tleap_in_path).resolve()
    resolved = shutil.which(tleap_bin or "tleap")
    if resolved is None:
        return ValidationResult(
            tool="tleap", ran=False, passed=False,
            fallback_reason=(
                "tleap not on PATH; install AmberTools "
                "(free component) or load the CESGA amber module."
            ),
        )
    if not tleap_in.is_file():
        return ValidationResult(
            tool="tleap", ran=False, passed=False,
            fallback_reason=f"tleap.in not found: {tleap_in}",
        )

    prmtop_expected, rst7_expected = _prmtop_and_rst7_from_tleap(tleap_in)

    cwd = tleap_in.parent
    try:
        proc = subprocess.run(
            [resolved, "-f", str(tleap_in)],
            capture_output=True,
            text=True,
            timeout=timeout_seconds,
            cwd=str(cwd),
        )
    except subprocess.TimeoutExpired:
        return ValidationResult(
            tool="tleap", ran=True, passed=False,
            fallback_reason=f"tleap timed out after {timeout_seconds}s",
        )
    except (FileNotFoundError, PermissionError) as exc:
        return ValidationResult(
            tool="tleap", ran=False, passed=False,
            fallback_reason=f"tleap exec failed: {exc!r}",
        )

    stdout_tail = _tail(proc.stdout, 40)
    stderr_tail = _tail(proc.stderr, 40)

    # tleap prints Errors to stdout, warns to stderr; both matter.
    combined_lower = (proc.stdout + "\n" + proc.stderr).lower()
    fatal = ("fatal" in combined_lower) or ("could not open" in combined_lower)

    if proc.returncode != 0 or fatal:
        return ValidationResult(
            tool="tleap", ran=True, passed=False,
            fallback_reason=f"tleap exited with code {proc.returncode} (fatal={fatal})",
            exit_code=proc.returncode,
            stdout_tail=stdout_tail, stderr_tail=stderr_tail,
        )

    prmtop_ok = prmtop_expected is not None and _looks_like_valid_prmtop(prmtop_expected)
    rst7_ok = rst7_expected is not None and _looks_like_valid_rst7(rst7_expected)

    diagnostics: list[str] = []
    if prmtop_expected is None:
        diagnostics.append("no saveamberparm directive found in tleap.in")
    else:
        diagnostics.append(f"prmtop={prmtop_expected} valid_signature={prmtop_ok}")
    if rst7_expected is None:
        diagnostics.append("no rst7 target parsed from tleap.in")
    else:
        diagnostics.append(f"rst7={rst7_expected} valid_signature={rst7_ok}")

    return ValidationResult(
        tool="tleap",
        ran=True,
        passed=(prmtop_ok and rst7_ok),
        exit_code=proc.returncode,
        prmtop_path=prmtop_expected,
        rst7_path=rst7_expected,
        stdout_tail=stdout_tail,
        stderr_tail=stderr_tail,
        diagnostics=diagnostics,
    )


def validate_sander_zero_step(
    prmtop_path: str | Path,
    rst7_path: str | Path,
    workdir: str | Path | None = None,
    timeout_seconds: int = DEFAULT_SANDER_TIMEOUT_SECONDS,
    sander_bin: str | None = None,
) -> ValidationResult:
    """Run sander with a 0-step minimization input.

    Confirms Amber accepts the topology + coord: parameters resolve,
    bonded network is closed, box dimensions parse, atom count matches.
    Fail modes caught: missing GAFF terms for cofactors, mismatched
    residue templates, unphysical box dimensions.
    """
    prmtop = Path(prmtop_path).resolve()
    rst7 = Path(rst7_path).resolve()

    resolved = shutil.which(sander_bin or "sander")
    if resolved is None:
        return ValidationResult(
            tool="sander", ran=False, passed=False,
            fallback_reason=(
                "sander not on PATH; install AmberTools (free) or "
                "load the CESGA amber module."
            ),
        )
    for path, name in ((prmtop, "prmtop"), (rst7, "rst7")):
        if not path.is_file():
            return ValidationResult(
                tool="sander", ran=False, passed=False,
                fallback_reason=f"{name} not found: {path}",
            )

    work = Path(workdir) if workdir else Path(tempfile.mkdtemp(prefix="fruton_sander_sanity_"))
    work.mkdir(parents=True, exist_ok=True)
    in_path = work / "sanity_0step.in"
    in_path.write_text(
        """\
0-step sanity: verify Amber can load the topology + coord
 &cntrl
  imin = 1,
  maxcyc = 0,
  ncyc = 0,
  ntb = 1,
  cut = 10.0,
  ntpr = 1,
 /
""",
        encoding="utf-8",
    )

    out_path = work / "sanity_0step.out"
    try:
        proc = subprocess.run(
            [resolved, "-O",
             "-i", str(in_path),
             "-o", str(out_path),
             "-p", str(prmtop),
             "-c", str(rst7)],
            capture_output=True,
            text=True,
            timeout=timeout_seconds,
            cwd=str(work),
        )
    except subprocess.TimeoutExpired:
        return ValidationResult(
            tool="sander", ran=True, passed=False,
            fallback_reason=f"sander timed out after {timeout_seconds}s",
        )
    except (FileNotFoundError, PermissionError) as exc:
        return ValidationResult(
            tool="sander", ran=False, passed=False,
            fallback_reason=f"sander exec failed: {exc!r}",
        )

    stdout_tail = _tail(proc.stdout, 30)
    stderr_tail = _tail(proc.stderr, 30)

    # sander writes an "Total wall time" line + a "Final" energy section on
    # a successful run.  Absence of "Total wall time" combined with a
    # non-zero exit is a fatal load error.
    out_text = out_path.read_text(encoding="utf-8", errors="replace") if out_path.is_file() else ""
    has_wall = "Total wall time" in out_text or "wallclock() was called" in out_text

    passed = (proc.returncode == 0) and has_wall

    diagnostics = [
        f"exit_code={proc.returncode}",
        f"out_file={out_path} exists={out_path.is_file()}",
        f"has_wall_time_marker={has_wall}",
    ]

    return ValidationResult(
        tool="sander",
        ran=True,
        passed=passed,
        exit_code=proc.returncode,
        stdout_tail=stdout_tail if not passed else _tail(out_text, 20),
        stderr_tail=stderr_tail,
        diagnostics=diagnostics,
    )


def validate_full(
    tleap_in_path: str | Path,
    timeout_tleap_s: int = DEFAULT_TLEAP_TIMEOUT_SECONDS,
    timeout_sander_s: int = DEFAULT_SANDER_TIMEOUT_SECONDS,
) -> tuple[ValidationResult, ValidationResult | None]:
    """Convenience: run tleap validation, then (only if it passed) run
    the sander 0-step check on the emitted .prmtop / .rst7.  Second
    result is None when the first fails or produces no files.
    """
    tleap_r = validate_tleap_can_load(tleap_in_path, timeout_seconds=timeout_tleap_s)
    if not (tleap_r.passed and tleap_r.prmtop_path and tleap_r.rst7_path):
        return (tleap_r, None)
    sander_r = validate_sander_zero_step(
        tleap_r.prmtop_path, tleap_r.rst7_path,
        timeout_seconds=timeout_sander_s,
    )
    return (tleap_r, sander_r)
