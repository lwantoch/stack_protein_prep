"""Tests for stack_protein_preparation._external_generator_runner.

Exercises the shared subprocess-cascade helper directly (independent
of the two scaffold callers) so refactors that touch this file don't
have to rely on the scaffold-level tests to catch regressions.
"""
from __future__ import annotations

import subprocess
from pathlib import Path

import pytest

from stack_protein_preparation import _external_generator_runner as egr


class _FakeProc:
    def __init__(self, returncode: int, stdout: str = "", stderr: str = ""):
        self.returncode = returncode
        self.stdout = stdout
        self.stderr = stderr


# ---------------------------------------------------------------------------
# resolve_binary_from_candidates
# ---------------------------------------------------------------------------


def test_resolve_returns_first_hit(monkeypatch: pytest.MonkeyPatch):
    import shutil as _sh
    monkeypatch.setattr(
        _sh, "which",
        lambda name: "/opt/" + name if name == "second" else None,
    )
    assert egr.resolve_binary_from_candidates(["first", "second"]) == "/opt/second"


def test_resolve_returns_none_when_all_missing(monkeypatch: pytest.MonkeyPatch):
    import shutil as _sh
    monkeypatch.setattr(_sh, "which", lambda name: None)
    assert egr.resolve_binary_from_candidates(["a", "b", "c"]) is None


def test_resolve_short_circuits_on_first_hit(monkeypatch: pytest.MonkeyPatch):
    """Once a binary is found we do not keep scanning."""
    import shutil as _sh
    checked: list[str] = []

    def _which(name: str):
        checked.append(name)
        return "/opt/a" if name == "a" else None

    monkeypatch.setattr(_sh, "which", _which)
    assert egr.resolve_binary_from_candidates(["a", "b", "c"]) == "/opt/a"
    assert checked == ["a"]  # b and c never scanned


# ---------------------------------------------------------------------------
# execute_generator_subprocess — success path
# ---------------------------------------------------------------------------


def test_execute_success_collects_candidates(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch,
):
    out = tmp_path / "out"
    out.mkdir()

    def _fake_run(cmd, **kw):
        assert "--manifest" in cmd
        outdir_arg = Path(cmd[cmd.index("--output-dir") + 1])
        (outdir_arg / "candidate_1.pdb").write_text("ATOM\n")
        (outdir_arg / "candidate_2.pdb").write_text("ATOM\n")
        return _FakeProc(0, stdout="ok\n")

    monkeypatch.setattr(subprocess, "run", _fake_run)
    r = egr.execute_generator_subprocess(
        binary="/opt/tool",
        manifest_path=tmp_path / "manifest.json",
        output_dir=out,
        n_candidates=2,
        timeout_seconds=60,
        tool_name_for_messages="tool",
    )
    assert r.ran is True
    assert r.accepted is True
    assert len(r.candidate_pdbs) == 2
    assert r.used_binary == "/opt/tool"


def test_execute_success_but_no_candidates_rejects(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch,
):
    out = tmp_path / "out"
    out.mkdir()

    def _fake_run(cmd, **kw):
        return _FakeProc(0, stdout="nothing produced\n")

    monkeypatch.setattr(subprocess, "run", _fake_run)
    r = egr.execute_generator_subprocess(
        binary="/opt/tool",
        manifest_path=tmp_path / "manifest.json",
        output_dir=out,
        n_candidates=2,
        timeout_seconds=60,
        tool_name_for_messages="tool",
    )
    assert r.ran is True
    assert r.accepted is False
    assert "tool ran but produced no" in r.fallback_reason


# ---------------------------------------------------------------------------
# execute_generator_subprocess — failure paths
# ---------------------------------------------------------------------------


def test_execute_nonzero_exit_reports_stderr(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch,
):
    def _fake_run(cmd, **kw):
        return _FakeProc(2, stderr="CUDA OOM: not enough VRAM\n")

    monkeypatch.setattr(subprocess, "run", _fake_run)
    r = egr.execute_generator_subprocess(
        binary="/opt/tool",
        manifest_path=tmp_path / "manifest.json",
        output_dir=tmp_path,
        n_candidates=2,
        timeout_seconds=60,
        tool_name_for_messages="my_generator",
    )
    assert r.ran is True
    assert r.accepted is False
    assert "my_generator exit=2" in r.fallback_reason
    assert "CUDA OOM" in r.fallback_reason


def test_execute_timeout_returns_ran_true(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch,
):
    def _fake_run(cmd, **kw):
        raise subprocess.TimeoutExpired(cmd, kw.get("timeout", 0))

    monkeypatch.setattr(subprocess, "run", _fake_run)
    r = egr.execute_generator_subprocess(
        binary="/opt/tool",
        manifest_path=tmp_path / "manifest.json",
        output_dir=tmp_path,
        n_candidates=2,
        timeout_seconds=5,
        tool_name_for_messages="tool",
    )
    assert r.ran is True
    assert r.accepted is False
    assert "tool timed out after 5s" in r.fallback_reason


def test_execute_os_error_returns_ran_false(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch,
):
    def _fake_run(cmd, **kw):
        raise OSError("cannot exec")

    monkeypatch.setattr(subprocess, "run", _fake_run)
    r = egr.execute_generator_subprocess(
        binary="/opt/tool",
        manifest_path=tmp_path / "manifest.json",
        output_dir=tmp_path,
        n_candidates=2,
        timeout_seconds=60,
        tool_name_for_messages="tool",
    )
    assert r.ran is False
    assert r.accepted is False
    assert "failed to spawn /opt/tool" in r.fallback_reason


# ---------------------------------------------------------------------------
# custom candidate glob
# ---------------------------------------------------------------------------


def test_execute_custom_candidate_glob(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch,
):
    out = tmp_path / "out"
    out.mkdir()

    def _fake_run(cmd, **kw):
        outdir = Path(cmd[cmd.index("--output-dir") + 1])
        (outdir / "pred_0.pdb").write_text("ATOM\n")
        (outdir / "pred_1.pdb").write_text("ATOM\n")
        return _FakeProc(0)

    monkeypatch.setattr(subprocess, "run", _fake_run)
    r = egr.execute_generator_subprocess(
        binary="/opt/tool",
        manifest_path=tmp_path / "manifest.json",
        output_dir=out,
        n_candidates=2,
        timeout_seconds=60,
        tool_name_for_messages="tool",
        candidate_glob="pred_*.pdb",
    )
    assert r.accepted is True
    assert len(r.candidate_pdbs) == 2


def test_execute_used_binary_set_even_on_failure(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch,
):
    def _fake_run(cmd, **kw):
        return _FakeProc(1)

    monkeypatch.setattr(subprocess, "run", _fake_run)
    r = egr.execute_generator_subprocess(
        binary="/opt/x",
        manifest_path=tmp_path / "m.json",
        output_dir=tmp_path,
        n_candidates=2,
        timeout_seconds=60,
        tool_name_for_messages="tool",
    )
    assert r.used_binary == "/opt/x"
