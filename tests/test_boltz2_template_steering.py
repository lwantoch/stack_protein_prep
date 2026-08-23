"""Tests for stack_protein_preparation._boltz2_template_steering.

Never launches a real Boltz-2 binary — shutil.which is monkey-patched
and subprocess.run is replaced with a fake for the happy path.
"""
from __future__ import annotations

import json
import subprocess
from pathlib import Path

import pytest

from stack_protein_preparation import _boltz2_template_steering as bt


def _write_ctx(path: Path) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        "ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00  0.00           C\n"
        "ATOM      2  CA  ALA A  10      15.000   0.000   0.000  1.00  0.00           C\n"
        "TER\nEND\n"
    )
    return path


def _write_tpl(path: Path) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        "ATOM      1  CA  ALA B   1       0.000   0.000   0.000  1.00  0.00           C\n"
        "TER\nEND\n"
    )
    return path


def _spec_ok(tmp_path: Path) -> bt.BoltzTemplateSteeringSpec:
    tpl = _write_tpl(tmp_path / "tpl.pdb")
    return bt.BoltzTemplateSteeringSpec(
        chain="A",
        first_resnum=2, last_resnum=9,
        sequence="AAAAAAAA",
        template_pdb=tpl,
    )


class _FakeProc:
    def __init__(self, returncode: int, stdout: str = "", stderr: str = ""):
        self.returncode = returncode
        self.stdout = stdout
        self.stderr = stderr


# ---------------------------------------------------------------------------
# BoltzTemplateSteeringSpec
# ---------------------------------------------------------------------------


def test_spec_n_residues(tmp_path: Path):
    assert _spec_ok(tmp_path).n_residues() == 8


def test_spec_sanity_ok(tmp_path: Path):
    assert _spec_ok(tmp_path).sanity_check() == []


def test_spec_sanity_bad_seq_length(tmp_path: Path):
    tpl = _write_tpl(tmp_path / "tpl.pdb")
    bad = bt.BoltzTemplateSteeringSpec(
        chain="A", first_resnum=2, last_resnum=9,
        sequence="AA", template_pdb=tpl,
    )
    assert any("sequence length" in p for p in bad.sanity_check())


def test_spec_sanity_empty_gap(tmp_path: Path):
    tpl = _write_tpl(tmp_path / "tpl.pdb")
    bad = bt.BoltzTemplateSteeringSpec(
        chain="A", first_resnum=9, last_resnum=2,  # inverted
        sequence="", template_pdb=tpl,
    )
    assert any("empty gap" in p for p in bad.sanity_check())


def test_spec_sanity_missing_template(tmp_path: Path):
    bad = bt.BoltzTemplateSteeringSpec(
        chain="A", first_resnum=2, last_resnum=9,
        sequence="AAAAAAAA",
        template_pdb=tmp_path / "nope.pdb",
    )
    assert any("template pdb not found" in p for p in bad.sanity_check())


def test_spec_default_plddt_threshold_is_50(tmp_path: Path):
    assert _spec_ok(tmp_path).pLDDT_threshold == 50.0


# ---------------------------------------------------------------------------
# _resolve_binary
# ---------------------------------------------------------------------------


def test_resolve_binary_returns_first_hit(monkeypatch: pytest.MonkeyPatch):
    import shutil as _sh
    monkeypatch.setattr(
        _sh, "which",
        lambda name: "/opt/x/boltz2" if name == "boltz2" else None,
    )
    assert bt._resolve_binary() == "/opt/x/boltz2"


def test_resolve_binary_returns_none_when_all_missing(
    monkeypatch: pytest.MonkeyPatch,
):
    import shutil as _sh
    monkeypatch.setattr(_sh, "which", lambda name: None)
    assert bt._resolve_binary() is None


def test_resolve_binary_honours_override(monkeypatch: pytest.MonkeyPatch):
    import shutil as _sh
    monkeypatch.setattr(_sh, "which",
                        lambda name: "/x/" + name if name == "custom" else None)
    assert bt._resolve_binary(["custom"]) == "/x/custom"


# ---------------------------------------------------------------------------
# _write_steering_manifest
# ---------------------------------------------------------------------------


def test_manifest_writes_expected_fields(tmp_path: Path):
    ctx = _write_ctx(tmp_path / "ctx.pdb")
    spec = _spec_ok(tmp_path)
    path = bt._write_steering_manifest(tmp_path, spec, ctx)
    payload = json.loads(path.read_text())
    assert payload["generator"] == "boltz2_template_steering"
    assert payload["template_pdb"].endswith("tpl.pdb")
    assert payload["gap"]["chain"] == "A"
    assert payload["gap"]["first_resnum"] == 2
    assert payload["gap"]["last_resnum"] == 9
    assert payload["steering"]["mode"] == "template_bias"


# ---------------------------------------------------------------------------
# attempt_template_steering: fail-open paths
# ---------------------------------------------------------------------------


def test_missing_context_pdb(tmp_path: Path):
    r = bt.attempt_template_steering(
        tmp_path / "nope.pdb", _spec_ok(tmp_path), tmp_path / "out",
    )
    assert r.ran is False
    assert "context pdb not found" in r.fallback_reason


def test_invalid_spec(tmp_path: Path):
    ctx = _write_ctx(tmp_path / "ctx.pdb")
    tpl = _write_tpl(tmp_path / "tpl.pdb")
    bad = bt.BoltzTemplateSteeringSpec(
        chain="A", first_resnum=2, last_resnum=9,
        sequence="A", template_pdb=tpl,
    )
    r = bt.attempt_template_steering(ctx, bad, tmp_path / "out")
    assert r.ran is False
    assert "invalid" in r.fallback_reason


def test_no_binary_on_path(tmp_path: Path, monkeypatch: pytest.MonkeyPatch):
    import shutil as _sh
    monkeypatch.setattr(_sh, "which", lambda name: None)
    ctx = _write_ctx(tmp_path / "ctx.pdb")
    r = bt.attempt_template_steering(ctx, _spec_ok(tmp_path), tmp_path / "out")
    assert r.ran is False
    assert "no boltz binary" in r.fallback_reason
    assert "MODELLER LoopModel" in r.fallback_reason  # advises fallback


# ---------------------------------------------------------------------------
# attempt_template_steering: subprocess paths (monkey-patched)
# ---------------------------------------------------------------------------


def _bin(monkeypatch: pytest.MonkeyPatch, name: str = "boltz"):
    import shutil as _sh
    monkeypatch.setattr(
        _sh, "which",
        lambda n: "/opt/x/" + n if n == name else None,
    )


def test_subprocess_success_produces_candidates(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch,
):
    _bin(monkeypatch)
    ctx = _write_ctx(tmp_path / "ctx.pdb")
    outdir = tmp_path / "out"

    def _fake_run(cmd, **kw):
        assert "--manifest" in cmd
        out = Path(cmd[cmd.index("--output-dir") + 1])
        out.mkdir(parents=True, exist_ok=True)
        (out / "candidate_1.pdb").write_text("ATOM\nEND\n")
        (out / "candidate_2.pdb").write_text("ATOM\nEND\n")
        (out / "candidate_3.pdb").write_text("ATOM\nEND\n")
        return _FakeProc(0, stdout="boltz ok\n")

    monkeypatch.setattr(subprocess, "run", _fake_run)
    r = bt.attempt_template_steering(ctx, _spec_ok(tmp_path), outdir)
    assert r.ran is True
    assert r.accepted is True
    assert len(r.candidate_pdbs) == 3
    assert r.used_binary == "/opt/x/boltz"


def test_subprocess_nonzero_exit_falls_back(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch,
):
    _bin(monkeypatch)
    ctx = _write_ctx(tmp_path / "ctx.pdb")

    def _fake_run(cmd, **kw):
        return _FakeProc(1, stderr="CUDA OOM\n")
    monkeypatch.setattr(subprocess, "run", _fake_run)

    r = bt.attempt_template_steering(ctx, _spec_ok(tmp_path), tmp_path / "out")
    assert r.ran is True
    assert r.accepted is False
    assert "exit=1" in r.fallback_reason
    assert "CUDA OOM" in r.fallback_reason


def test_subprocess_no_candidates_produced_rejected(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch,
):
    _bin(monkeypatch)
    ctx = _write_ctx(tmp_path / "ctx.pdb")

    def _fake_run(cmd, **kw):
        return _FakeProc(0, stdout="nothing\n")
    monkeypatch.setattr(subprocess, "run", _fake_run)

    r = bt.attempt_template_steering(ctx, _spec_ok(tmp_path), tmp_path / "out")
    assert r.ran is True
    assert r.accepted is False
    assert "no candidate" in r.fallback_reason


def test_subprocess_timeout_falls_back(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch,
):
    _bin(monkeypatch)
    ctx = _write_ctx(tmp_path / "ctx.pdb")

    def _fake_run(cmd, **kw):
        raise subprocess.TimeoutExpired(cmd, kw.get("timeout", 0))
    monkeypatch.setattr(subprocess, "run", _fake_run)

    r = bt.attempt_template_steering(
        ctx, _spec_ok(tmp_path), tmp_path / "out", timeout_seconds=1,
    )
    assert r.ran is True
    assert r.accepted is False
    assert "timed out" in r.fallback_reason


def test_subprocess_os_error_falls_back(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch,
):
    _bin(monkeypatch)
    ctx = _write_ctx(tmp_path / "ctx.pdb")

    def _fake_run(cmd, **kw):
        raise OSError("cannot exec")
    monkeypatch.setattr(subprocess, "run", _fake_run)

    r = bt.attempt_template_steering(ctx, _spec_ok(tmp_path), tmp_path / "out")
    assert r.ran is False
    assert "failed to spawn" in r.fallback_reason


# ---------------------------------------------------------------------------
# to_dict / summarise
# ---------------------------------------------------------------------------


def test_to_dict_json_serialisable():
    r = bt.BoltzSteeringAttempt(
        ran=True, accepted=True,
        candidate_pdbs=[Path("/tmp/candidate_1.pdb")],
        used_binary="/opt/x/boltz",
        diagnostics=["n=1"],
    )
    assert "candidate_pdbs" in json.dumps(r.to_dict())


def test_summarise_skipped():
    r = bt.BoltzSteeringAttempt(ran=False, accepted=False,
                                fallback_reason="not on PATH")
    assert "skipped" in bt.summarise(r)


def test_summarise_rejected():
    r = bt.BoltzSteeringAttempt(ran=True, accepted=False,
                                fallback_reason="exit=1")
    assert "rejected" in bt.summarise(r)


def test_summarise_accepted():
    r = bt.BoltzSteeringAttempt(
        ran=True, accepted=True,
        candidate_pdbs=[Path("/tmp/c_1.pdb"), Path("/tmp/c_2.pdb")],
        used_binary="/opt/x/boltz",
    )
    s = bt.summarise(r)
    assert "2 candidates" in s
    assert "boltz" in s
