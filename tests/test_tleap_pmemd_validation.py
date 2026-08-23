"""Tests for stack_protein_preparation._tleap_pmemd_validation.

Never launches real tleap / sander — subprocess.run is monkey-patched
so the tests exercise the wrapper logic (fail-open, exit-code parsing,
prmtop / rst7 signature checks, saveamberparm directive parsing)
independent of AmberTools availability.
"""
from __future__ import annotations

import subprocess
from pathlib import Path

import pytest

from stack_protein_preparation import _tleap_pmemd_validation as tpv


def _write_prmtop(path: Path, valid: bool = True) -> Path:
    """Emit a minimal file that passes / fails our signature check."""
    if valid:
        # Padded to > 512 bytes with %VERSION + %FLAG POINTERS markers
        header = (
            "%VERSION  VERSION_STAMP = V0001.000  DATE = 08/23/26  12:34:56\n"
            "%FLAG POINTERS\n"
            "%FORMAT(10I8)\n"
            + " " * 600
        )
        path.write_text(header, encoding="utf-8")
    else:
        path.write_text("not a real prmtop\n", encoding="utf-8")
    return path


def _write_rst7(path: Path, valid: bool = True) -> Path:
    if valid:
        path.write_text("title line\n     100    0.0000000E+00\n", encoding="utf-8")
    else:
        path.write_text("junk\n", encoding="utf-8")
    return path


def _write_tleap_in(path: Path, prmtop: Path, rst7: Path) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        f"source leaprc.protein.ff14SB\n"
        f"source leaprc.water.tip3p\n"
        f"mol = loadpdb /tmp/dummy.pdb\n"
        f"saveamberparm mol {prmtop} {rst7}\n"
        f"quit\n",
        encoding="utf-8",
    )
    return path


class _FakeProc:
    def __init__(self, returncode: int, stdout: str = "", stderr: str = ""):
        self.returncode = returncode
        self.stdout = stdout
        self.stderr = stderr


@pytest.fixture
def _no_amber_tools(monkeypatch: pytest.MonkeyPatch):
    """Force shutil.which to return None for tleap and sander."""
    import shutil as _sh
    monkeypatch.setattr(_sh, "which", lambda name: None)


# ---------------------------------------------------------------------------
# Fail-open when AmberTools is not on PATH
# ---------------------------------------------------------------------------


def test_tleap_missing_returns_fail_open(_no_amber_tools, tmp_path: Path):
    r = tpv.validate_tleap_can_load(tmp_path / "tleap.in")
    assert r.ran is False
    assert r.passed is False
    assert r.tool == "tleap"
    assert "not on PATH" in r.fallback_reason


def test_sander_missing_returns_fail_open(_no_amber_tools, tmp_path: Path):
    _write_prmtop(tmp_path / "p.prmtop")
    _write_rst7(tmp_path / "r.rst7")
    r = tpv.validate_sander_zero_step(tmp_path / "p.prmtop", tmp_path / "r.rst7")
    assert r.ran is False
    assert r.passed is False
    assert "not on PATH" in r.fallback_reason


def test_tleap_missing_input_file(tmp_path: Path, monkeypatch: pytest.MonkeyPatch):
    """tleap on PATH but tleap.in missing -> ran=False + helpful reason."""
    import shutil as _sh
    monkeypatch.setattr(_sh, "which", lambda name: "/opt/amber/bin/tleap")
    r = tpv.validate_tleap_can_load(tmp_path / "nope.in")
    assert r.ran is False
    assert "not found" in r.fallback_reason


def test_sander_missing_prmtop_or_rst7(tmp_path: Path, monkeypatch: pytest.MonkeyPatch):
    import shutil as _sh
    monkeypatch.setattr(_sh, "which", lambda name: "/opt/amber/bin/sander")
    r = tpv.validate_sander_zero_step(tmp_path / "nope.prmtop", tmp_path / "nope.rst7")
    assert r.ran is False
    assert "not found" in r.fallback_reason


# ---------------------------------------------------------------------------
# Fake subprocess: tleap paths
# ---------------------------------------------------------------------------


def test_tleap_success_with_valid_output_files(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch,
):
    import shutil as _sh
    monkeypatch.setattr(_sh, "which", lambda name: "/opt/amber/bin/tleap")
    prmtop = _write_prmtop(tmp_path / "system.prmtop", valid=True)
    rst7 = _write_rst7(tmp_path / "system.rst7", valid=True)
    tleap_in = _write_tleap_in(tmp_path / "tleap.in", prmtop, rst7)

    def _fake_run(cmd, **kw):
        return _FakeProc(returncode=0, stdout="tleap finished OK\n")
    monkeypatch.setattr(subprocess, "run", _fake_run)

    r = tpv.validate_tleap_can_load(tleap_in)
    assert r.passed is True
    assert r.exit_code == 0
    assert r.prmtop_path == prmtop.resolve()
    assert r.rst7_path == rst7.resolve()


def test_tleap_nonzero_exit_fails(tmp_path: Path, monkeypatch: pytest.MonkeyPatch):
    import shutil as _sh
    monkeypatch.setattr(_sh, "which", lambda name: "/opt/amber/bin/tleap")
    prmtop = tmp_path / "system.prmtop"
    rst7 = tmp_path / "system.rst7"
    tleap_in = _write_tleap_in(tmp_path / "tleap.in", prmtop, rst7)

    def _fake_run(cmd, **kw):
        return _FakeProc(returncode=1, stdout="FATAL: parameter missing for CX\n")
    monkeypatch.setattr(subprocess, "run", _fake_run)

    r = tpv.validate_tleap_can_load(tleap_in)
    assert r.passed is False
    assert r.exit_code == 1
    assert "fatal=True" in r.fallback_reason


def test_tleap_zero_exit_but_fatal_in_stdout(tmp_path: Path, monkeypatch: pytest.MonkeyPatch):
    """Some tleap versions return 0 even after a FATAL -- catch that."""
    import shutil as _sh
    monkeypatch.setattr(_sh, "which", lambda name: "/opt/amber/bin/tleap")
    prmtop = _write_prmtop(tmp_path / "system.prmtop")
    rst7 = _write_rst7(tmp_path / "system.rst7")
    tleap_in = _write_tleap_in(tmp_path / "tleap.in", prmtop, rst7)

    def _fake_run(cmd, **kw):
        return _FakeProc(returncode=0, stdout="FATAL:  could not open XYZ.frcmod\n")
    monkeypatch.setattr(subprocess, "run", _fake_run)

    r = tpv.validate_tleap_can_load(tleap_in)
    assert r.passed is False
    assert "fatal=True" in r.fallback_reason


def test_tleap_success_but_invalid_prmtop_signature_fails(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch,
):
    import shutil as _sh
    monkeypatch.setattr(_sh, "which", lambda name: "/opt/amber/bin/tleap")
    prmtop = _write_prmtop(tmp_path / "system.prmtop", valid=False)  # bad signature
    rst7 = _write_rst7(tmp_path / "system.rst7", valid=True)
    tleap_in = _write_tleap_in(tmp_path / "tleap.in", prmtop, rst7)

    def _fake_run(cmd, **kw):
        return _FakeProc(returncode=0, stdout="tleap finished OK\n")
    monkeypatch.setattr(subprocess, "run", _fake_run)

    r = tpv.validate_tleap_can_load(tleap_in)
    assert r.passed is False  # prmtop signature invalid


def test_tleap_timeout(tmp_path: Path, monkeypatch: pytest.MonkeyPatch):
    import shutil as _sh
    monkeypatch.setattr(_sh, "which", lambda name: "/opt/amber/bin/tleap")
    tleap_in = _write_tleap_in(
        tmp_path / "tleap.in", tmp_path / "p.prmtop", tmp_path / "r.rst7",
    )

    def _fake_run(cmd, **kw):
        raise subprocess.TimeoutExpired(cmd=cmd, timeout=kw.get("timeout", 0))
    monkeypatch.setattr(subprocess, "run", _fake_run)

    r = tpv.validate_tleap_can_load(tleap_in, timeout_seconds=1)
    assert r.ran is True
    assert r.passed is False
    assert "timed out" in r.fallback_reason


# ---------------------------------------------------------------------------
# Prmtop / rst7 signature helpers
# ---------------------------------------------------------------------------


def test_valid_prmtop_signature_check(tmp_path: Path):
    p = _write_prmtop(tmp_path / "p.prmtop", valid=True)
    assert tpv._looks_like_valid_prmtop(p) is True


def test_invalid_prmtop_signature_check(tmp_path: Path):
    p = _write_prmtop(tmp_path / "p.prmtop", valid=False)
    assert tpv._looks_like_valid_prmtop(p) is False


def test_valid_rst7_signature_check(tmp_path: Path):
    r = _write_rst7(tmp_path / "r.rst7", valid=True)
    assert tpv._looks_like_valid_rst7(r) is True


def test_invalid_rst7_signature_check(tmp_path: Path):
    r = _write_rst7(tmp_path / "r.rst7", valid=False)
    assert tpv._looks_like_valid_rst7(r) is False


def test_saveamberparm_parser(tmp_path: Path):
    prmtop = tmp_path / "system.prmtop"
    rst7 = tmp_path / "system.rst7"
    tleap_in = _write_tleap_in(tmp_path / "tleap.in", prmtop, rst7)
    p, r = tpv._prmtop_and_rst7_from_tleap(tleap_in)
    assert p == prmtop.resolve()
    assert r == rst7.resolve()


def test_saveamberparm_parser_missing_directive(tmp_path: Path):
    """tleap.in without saveamberparm -> (None, None)."""
    (tmp_path / "tleap.in").write_text("source leaprc.protein.ff14SB\nquit\n")
    p, r = tpv._prmtop_and_rst7_from_tleap(tmp_path / "tleap.in")
    assert p is None and r is None


# ---------------------------------------------------------------------------
# validate_full convenience wrapper
# ---------------------------------------------------------------------------


def test_validate_full_skips_sander_when_tleap_fails(
    _no_amber_tools, tmp_path: Path,
):
    tleap_in = _write_tleap_in(
        tmp_path / "tleap.in", tmp_path / "p.prmtop", tmp_path / "r.rst7",
    )
    tleap_r, sander_r = tpv.validate_full(tleap_in)
    assert tleap_r.passed is False
    assert sander_r is None
