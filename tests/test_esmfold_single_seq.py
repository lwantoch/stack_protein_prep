"""Tests for stack_protein_preparation._esmfold_single_seq.

Never launches a real esmfold binary — shutil.which is monkey-patched
and subprocess.run is faked when needed.
"""
from __future__ import annotations

import json
import subprocess
from pathlib import Path

import pytest

from stack_protein_preparation import _esmfold_single_seq as ef


class _FakeProc:
    def __init__(self, returncode: int, stdout: str = "", stderr: str = ""):
        self.returncode = returncode
        self.stdout = stdout
        self.stderr = stderr


def _spec_ok() -> ef.EsmFoldSequenceSpec:
    return ef.EsmFoldSequenceSpec(
        chain_id="A",
        sequence="MVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHF",
    )


# ---------------------------------------------------------------------------
# EsmFoldSequenceSpec
# ---------------------------------------------------------------------------


def test_spec_ok_passes_sanity():
    assert _spec_ok().sanity_check() == []


def test_spec_empty_sequence_flagged():
    bad = ef.EsmFoldSequenceSpec(chain_id="A", sequence="")
    assert any("empty sequence" in p for p in bad.sanity_check())


def test_spec_non_standard_aa_flagged():
    bad = ef.EsmFoldSequenceSpec(chain_id="A", sequence="MVLZUB")
    problems = bad.sanity_check()
    assert any("non-standard amino-acid" in p for p in problems)


def test_spec_wildcard_X_allowed():
    """X is the standard wildcard for unknown residue; must not be flagged."""
    ok = ef.EsmFoldSequenceSpec(chain_id="A", sequence="MVLXY")
    assert ok.sanity_check() == []


def test_spec_max_recycles_out_of_range_flagged():
    bad = ef.EsmFoldSequenceSpec(chain_id="A", sequence="AA", max_recycles=0)
    assert any("max_recycles" in p for p in bad.sanity_check())
    bad_high = ef.EsmFoldSequenceSpec(chain_id="A", sequence="AA", max_recycles=99)
    assert any("max_recycles" in p for p in bad_high.sanity_check())


# ---------------------------------------------------------------------------
# _resolve_binary
# ---------------------------------------------------------------------------


def test_resolve_returns_first_hit(monkeypatch: pytest.MonkeyPatch):
    import shutil as _sh
    monkeypatch.setattr(
        _sh, "which",
        lambda name: "/opt/x/esm-fold" if name == "esm-fold" else None,
    )
    assert ef._resolve_binary() == "/opt/x/esm-fold"


def test_resolve_returns_none(monkeypatch: pytest.MonkeyPatch):
    import shutil as _sh
    monkeypatch.setattr(_sh, "which", lambda name: None)
    assert ef._resolve_binary() is None


# ---------------------------------------------------------------------------
# _write_esmfold_manifest
# ---------------------------------------------------------------------------


def test_manifest_writes_expected_fields(tmp_path: Path):
    spec = _spec_ok()
    path = ef._write_esmfold_manifest(tmp_path, spec)
    payload = json.loads(path.read_text())
    assert payload["generator"] == "esmfold_single_sequence"
    assert payload["sequence"]["chain_id"] == "A"
    assert payload["sequence"]["length"] == len(spec.sequence)
    assert payload["inference"]["max_recycles"] == 3


# ---------------------------------------------------------------------------
# attempt_single_sequence_fold: fail-open paths
# ---------------------------------------------------------------------------


def test_attempt_invalid_spec(tmp_path: Path):
    r = ef.attempt_single_sequence_fold(
        ef.EsmFoldSequenceSpec(chain_id="A", sequence=""),
        tmp_path / "out",
    )
    assert r.ran is False
    assert "invalid" in r.fallback_reason


def test_attempt_no_binary_on_path(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch,
):
    import shutil as _sh
    monkeypatch.setattr(_sh, "which", lambda name: None)
    r = ef.attempt_single_sequence_fold(_spec_ok(), tmp_path / "out")
    assert r.ran is False
    assert "no esmfold binary" in r.fallback_reason
    assert "AF-splice" in r.fallback_reason  # advises fallback path


# ---------------------------------------------------------------------------
# attempt_single_sequence_fold: subprocess paths (via shared runner)
# ---------------------------------------------------------------------------


def _bin(monkeypatch: pytest.MonkeyPatch, name: str = "esmfold"):
    import shutil as _sh
    monkeypatch.setattr(
        _sh, "which",
        lambda n: "/opt/x/" + n if n == name else None,
    )


def test_subprocess_success_produces_candidate(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch,
):
    _bin(monkeypatch)
    outdir = tmp_path / "out"

    def _fake_run(cmd, **kw):
        assert "--manifest" in cmd
        out = Path(cmd[cmd.index("--output-dir") + 1])
        out.mkdir(parents=True, exist_ok=True)
        (out / "candidate_1.pdb").write_text("ATOM\n")
        return _FakeProc(0, stdout="esmfold ok\n")

    monkeypatch.setattr(subprocess, "run", _fake_run)
    r = ef.attempt_single_sequence_fold(_spec_ok(), outdir)
    assert r.ran is True
    assert r.accepted is True
    assert len(r.candidate_pdbs) == 1
    assert r.used_binary == "/opt/x/esmfold"


def test_subprocess_nonzero_exit_falls_back(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch,
):
    _bin(monkeypatch)

    def _fake_run(cmd, **kw):
        return _FakeProc(1, stderr="CUDA OOM\n")

    monkeypatch.setattr(subprocess, "run", _fake_run)
    r = ef.attempt_single_sequence_fold(_spec_ok(), tmp_path / "out")
    assert r.ran is True
    assert r.accepted is False
    assert "esmfold exit=1" in r.fallback_reason


def test_subprocess_timeout_falls_back(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch,
):
    _bin(monkeypatch)

    def _fake_run(cmd, **kw):
        raise subprocess.TimeoutExpired(cmd, kw.get("timeout", 0))

    monkeypatch.setattr(subprocess, "run", _fake_run)
    r = ef.attempt_single_sequence_fold(
        _spec_ok(), tmp_path / "out", timeout_seconds=1,
    )
    assert r.ran is True
    assert r.accepted is False
    assert "esmfold timed out" in r.fallback_reason


# ---------------------------------------------------------------------------
# to_dict / summarise
# ---------------------------------------------------------------------------


def test_to_dict_json_serialisable():
    r = ef.EsmFoldAttempt(
        ran=True, accepted=True,
        candidate_pdbs=[Path("/tmp/c_1.pdb")],
        used_binary="/opt/x/esmfold",
    )
    json.dumps(r.to_dict())


def test_summarise_skipped():
    r = ef.EsmFoldAttempt(ran=False, accepted=False,
                          fallback_reason="not on PATH")
    assert "skipped" in ef.summarise(r)


def test_summarise_accepted():
    r = ef.EsmFoldAttempt(
        ran=True, accepted=True,
        candidate_pdbs=[Path("/tmp/c_1.pdb")],
        used_binary="/opt/x/esmfold",
    )
    s = ef.summarise(r)
    assert "1 candidates" in s
    assert "esmfold" in s
