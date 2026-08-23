"""Tests for stack_protein_preparation._rfdiffusion2_gap.

Never launches a real diffusion binary — shutil.which is monkey-patched
and subprocess.run is replaced with a fake for the happy path.
"""
from __future__ import annotations

import json
import subprocess
from pathlib import Path

import pytest

from stack_protein_preparation import _rfdiffusion2_gap as rd


def _write_context(path: Path) -> Path:
    """Tiny 2-residue PDB used as the 'existing structure' context."""
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        "ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00  0.00           C\n"
        "ATOM      2  CA  ALA A  10      15.000   0.000   0.000  1.00  0.00           C\n"
        "TER\nEND\n"
    )
    return path


def _spec_ok() -> rd.GapConstraintSpec:
    return rd.GapConstraintSpec(
        chain="A",
        first_resnum=2,
        last_resnum=9,
        sequence="AAAAAAAA",  # 8 residues
        anchor_ca_by_resnum={
            1: (0.0, 0.0, 0.0),
            10: (15.0, 0.0, 0.0),
        },
    )


class _FakeProc:
    def __init__(self, returncode: int, stdout: str = "", stderr: str = ""):
        self.returncode = returncode
        self.stdout = stdout
        self.stderr = stderr


# ---------------------------------------------------------------------------
# GapConstraintSpec
# ---------------------------------------------------------------------------


def test_gap_spec_n_residues():
    spec = _spec_ok()
    assert spec.n_residues() == 8


def test_gap_spec_sanity_check_ok():
    assert _spec_ok().sanity_check() == []


def test_gap_spec_sanity_check_length_mismatch():
    spec = rd.GapConstraintSpec(
        chain="A", first_resnum=2, last_resnum=9,
        sequence="AAA",  # only 3 vs 8 residues
        anchor_ca_by_resnum={1: (0, 0, 0)},
    )
    problems = spec.sanity_check()
    assert any("sequence length" in p for p in problems)


def test_gap_spec_sanity_check_empty_range():
    spec = rd.GapConstraintSpec(
        chain="A", first_resnum=5, last_resnum=4,  # inverted
        sequence="",
        anchor_ca_by_resnum={5: (0, 0, 0)},
    )
    problems = spec.sanity_check()
    assert any("empty gap" in p for p in problems)


def test_gap_spec_sanity_check_no_anchors():
    spec = rd.GapConstraintSpec(
        chain="A", first_resnum=2, last_resnum=9,
        sequence="AAAAAAAA",
        anchor_ca_by_resnum={},
    )
    problems = spec.sanity_check()
    assert any("anchor" in p for p in problems)


# ---------------------------------------------------------------------------
# _resolve_binary
# ---------------------------------------------------------------------------


def test_resolve_binary_returns_first_hit(monkeypatch: pytest.MonkeyPatch):
    import shutil as _sh
    seen: list[str] = []

    def _fake_which(name):
        seen.append(name)
        return "/opt/tools/" + name if name == "rfdiffusion" else None

    monkeypatch.setattr(_sh, "which", _fake_which)
    result = rd._resolve_binary()
    assert result == "/opt/tools/rfdiffusion"
    # rfdiffusion2 attempted first, then rfdiffusion
    assert seen[:2] == ["rfdiffusion2", "rfdiffusion"]


def test_resolve_binary_returns_none_when_all_missing(
    monkeypatch: pytest.MonkeyPatch,
):
    import shutil as _sh
    monkeypatch.setattr(_sh, "which", lambda name: None)
    assert rd._resolve_binary() is None


def test_resolve_binary_honours_override(monkeypatch: pytest.MonkeyPatch):
    import shutil as _sh
    monkeypatch.setattr(_sh, "which",
                        lambda name: "/x/" + name if name == "my_diffusion" else None)
    assert rd._resolve_binary(["my_diffusion"]) == "/x/my_diffusion"


# ---------------------------------------------------------------------------
# _write_constraint_manifest
# ---------------------------------------------------------------------------


def test_manifest_writes_expected_fields(tmp_path: Path):
    ctx = _write_context(tmp_path / "ctx.pdb")
    spec = _spec_ok()
    path = rd._write_constraint_manifest(tmp_path, spec, ctx)
    payload = json.loads(path.read_text())
    assert payload["generator"] == "constrained_diffusion"
    assert payload["context_pdb"].endswith("ctx.pdb")
    assert payload["gap"]["chain"] == "A"
    assert payload["gap"]["first_resnum"] == 2
    assert payload["gap"]["last_resnum"] == 9
    assert payload["gap"]["sequence"] == "AAAAAAAA"
    assert set(payload["anchor_constraints"]["ca_xyz_by_resnum"]) == {"1", "10"}


# ---------------------------------------------------------------------------
# attempt_gap_fill: fail-open paths
# ---------------------------------------------------------------------------


def test_attempt_missing_context_pdb(tmp_path: Path):
    r = rd.attempt_gap_fill(
        tmp_path / "nope.pdb", _spec_ok(), tmp_path / "out",
    )
    assert r.ran is False
    assert r.accepted is False
    assert "context pdb not found" in r.fallback_reason


def test_attempt_invalid_gap_spec(tmp_path: Path):
    ctx = _write_context(tmp_path / "ctx.pdb")
    bad_spec = rd.GapConstraintSpec(
        chain="A", first_resnum=2, last_resnum=9,
        sequence="A", anchor_ca_by_resnum={1: (0, 0, 0)},
    )
    r = rd.attempt_gap_fill(ctx, bad_spec, tmp_path / "out")
    assert r.ran is False
    assert "invalid" in r.fallback_reason


def test_attempt_no_binary_on_path(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch,
):
    import shutil as _sh
    monkeypatch.setattr(_sh, "which", lambda name: None)
    ctx = _write_context(tmp_path / "ctx.pdb")
    r = rd.attempt_gap_fill(ctx, _spec_ok(), tmp_path / "out")
    assert r.ran is False
    assert "no constrained-diffusion binary" in r.fallback_reason


# ---------------------------------------------------------------------------
# attempt_gap_fill: subprocess paths (monkey-patched)
# ---------------------------------------------------------------------------


def _monkey_bin(monkeypatch: pytest.MonkeyPatch, name: str = "rfdiffusion2"):
    import shutil as _sh
    monkeypatch.setattr(
        _sh, "which",
        lambda n: "/opt/x/" + n if n == name else None,
    )


def test_attempt_subprocess_success_produces_candidates(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch,
):
    _monkey_bin(monkeypatch)
    ctx = _write_context(tmp_path / "ctx.pdb")
    outdir = tmp_path / "out"

    def _fake_run(cmd, **kw):
        # Simulate the CLI writing two candidate PDBs.
        assert "--manifest" in cmd
        assert "--output-dir" in cmd
        out = Path(cmd[cmd.index("--output-dir") + 1])
        out.mkdir(parents=True, exist_ok=True)
        (out / "candidate_1.pdb").write_text("ATOM\nEND\n")
        (out / "candidate_2.pdb").write_text("ATOM\nEND\n")
        return _FakeProc(0, stdout="diffusion ok\n")

    monkeypatch.setattr(subprocess, "run", _fake_run)
    r = rd.attempt_gap_fill(ctx, _spec_ok(), outdir)
    assert r.ran is True
    assert r.accepted is True
    assert len(r.candidate_pdbs) == 2
    assert r.used_binary == "/opt/x/rfdiffusion2"


def test_attempt_subprocess_nonzero_exit_falls_back(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch,
):
    _monkey_bin(monkeypatch)
    ctx = _write_context(tmp_path / "ctx.pdb")

    def _fake_run(cmd, **kw):
        return _FakeProc(1, stderr="CUDA OOM\n")
    monkeypatch.setattr(subprocess, "run", _fake_run)

    r = rd.attempt_gap_fill(ctx, _spec_ok(), tmp_path / "out")
    assert r.ran is True
    assert r.accepted is False
    assert "exit=1" in r.fallback_reason
    assert "CUDA OOM" in r.fallback_reason


def test_attempt_subprocess_success_but_no_candidates_produced(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch,
):
    _monkey_bin(monkeypatch)
    ctx = _write_context(tmp_path / "ctx.pdb")

    def _fake_run(cmd, **kw):
        # Returncode 0 but writes nothing → treated as rejection.
        return _FakeProc(0, stdout="no output\n")
    monkeypatch.setattr(subprocess, "run", _fake_run)

    r = rd.attempt_gap_fill(ctx, _spec_ok(), tmp_path / "out")
    assert r.ran is True
    assert r.accepted is False
    assert "no candidate" in r.fallback_reason


def test_attempt_subprocess_timeout_falls_back(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch,
):
    _monkey_bin(monkeypatch)
    ctx = _write_context(tmp_path / "ctx.pdb")

    def _fake_run(cmd, **kw):
        raise subprocess.TimeoutExpired(cmd, kw.get("timeout", 0))
    monkeypatch.setattr(subprocess, "run", _fake_run)

    r = rd.attempt_gap_fill(ctx, _spec_ok(), tmp_path / "out", timeout_seconds=1)
    assert r.ran is True
    assert r.accepted is False
    assert "timed out" in r.fallback_reason


def test_attempt_subprocess_os_error_falls_back(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch,
):
    _monkey_bin(monkeypatch)
    ctx = _write_context(tmp_path / "ctx.pdb")

    def _fake_run(cmd, **kw):
        raise OSError("cannot exec")
    monkeypatch.setattr(subprocess, "run", _fake_run)

    r = rd.attempt_gap_fill(ctx, _spec_ok(), tmp_path / "out")
    assert r.ran is False
    assert "failed to spawn" in r.fallback_reason


# ---------------------------------------------------------------------------
# to_dict / summarise
# ---------------------------------------------------------------------------


def test_to_dict_json_serialisable():
    r = rd.RFDiffusionAttempt(
        ran=True, accepted=True,
        candidate_pdbs=[Path("/tmp/candidate_1.pdb")],
        used_binary="/opt/x/rfdiffusion2",
        diagnostics=["n=1"],
    )
    d = r.to_dict()
    s = json.dumps(d)
    assert "candidate_pdbs" in s
    assert "used_binary" in s


def test_summarise_skipped():
    r = rd.RFDiffusionAttempt(ran=False, accepted=False,
                              fallback_reason="not on PATH")
    assert "skipped" in rd.summarise(r)


def test_summarise_rejected():
    r = rd.RFDiffusionAttempt(ran=True, accepted=False,
                              fallback_reason="exit=1")
    assert "rejected" in rd.summarise(r)


def test_summarise_accepted():
    r = rd.RFDiffusionAttempt(
        ran=True, accepted=True,
        candidate_pdbs=[Path("/tmp/c_1.pdb"), Path("/tmp/c_2.pdb")],
        used_binary="/opt/x/rfdiffusion2",
    )
    s = rd.summarise(r)
    assert "2 candidates" in s
    assert "rfdiffusion2" in s
