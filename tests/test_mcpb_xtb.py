"""Tests for stack_protein_preparation._mcpb_xtb.

Never launches real xtb — subprocess.run is monkey-patched so the
tests exercise the wrapper logic (fail-open, exit-code parsing,
output-file discovery) independent of xtb availability.
"""
from __future__ import annotations

import subprocess
from pathlib import Path

import pytest

from stack_protein_preparation import _mcpb_xtb as xtb


class _FakeProc:
    def __init__(self, returncode: int, stdout: str = "", stderr: str = ""):
        self.returncode = returncode
        self.stdout = stdout
        self.stderr = stderr


# ---------------------------------------------------------------------------
# xtb_binary_available
# ---------------------------------------------------------------------------


def test_binary_available_returns_path_when_present(monkeypatch):
    import shutil as _sh
    monkeypatch.setattr(_sh, "which", lambda name: "/opt/xtb/bin/xtb")
    assert xtb.xtb_binary_available() == "/opt/xtb/bin/xtb"


def test_binary_available_returns_none_when_missing(monkeypatch):
    import shutil as _sh
    monkeypatch.setattr(_sh, "which", lambda name: None)
    assert xtb.xtb_binary_available() is None


# ---------------------------------------------------------------------------
# write_cluster_xyz
# ---------------------------------------------------------------------------


def test_write_cluster_xyz_shape(tmp_path):
    atoms = [("Zn", 0.0, 0.0, 0.0), ("N", 2.1, 0.0, 0.0), ("O", 0.0, 2.2, 0.0)]
    out = xtb.write_cluster_xyz(atoms, tmp_path / "cluster.xyz")
    lines = out.read_text().splitlines()
    assert lines[0] == "3"
    assert lines[1].startswith("FRUTON")
    assert lines[2].strip().startswith("Zn")
    assert lines[3].strip().startswith("N")
    assert lines[4].strip().startswith("O")


# ---------------------------------------------------------------------------
# run_xtb_optimize_hessian — fail-open paths
# ---------------------------------------------------------------------------


def test_missing_xtb_binary_fail_open(tmp_path, monkeypatch):
    import shutil as _sh
    monkeypatch.setattr(_sh, "which", lambda name: None)
    xyz = xtb.write_cluster_xyz([("Zn", 0, 0, 0)], tmp_path / "c.xyz")
    r = xtb.run_xtb_optimize_hessian(xyz, 2, 1, tmp_path / "work")
    assert r.ran is False
    assert r.passed is False
    assert "not on PATH" in r.fallback_reason


def test_missing_cluster_xyz_fail_open(tmp_path, monkeypatch):
    import shutil as _sh
    monkeypatch.setattr(_sh, "which", lambda name: "/opt/xtb/bin/xtb")
    r = xtb.run_xtb_optimize_hessian(tmp_path / "nope.xyz", 2, 1, tmp_path / "work")
    assert r.ran is False
    assert "cluster xyz not found" in r.fallback_reason


# ---------------------------------------------------------------------------
# run_xtb_optimize_hessian — happy path via fake subprocess
# ---------------------------------------------------------------------------


def test_successful_run_finds_output_files(tmp_path, monkeypatch):
    import shutil as _sh
    monkeypatch.setattr(_sh, "which", lambda name: "/opt/xtb/bin/xtb")
    xyz = xtb.write_cluster_xyz([("Zn", 0, 0, 0)], tmp_path / "c.xyz")
    workdir = tmp_path / "work"

    def _fake_run(cmd, **kw):
        # Simulate xtb writing its output files in cwd
        cwd = Path(kw["cwd"])
        (cwd / "xtbopt.xyz").write_text("1\nfake\nZn 0 0 0\n")
        (cwd / "hessian").write_text("$hessian\n1.0 0.0 0.0\n")
        (cwd / "charges").write_text("1.234\n")
        (cwd / "wbo").write_text("1 2 0.98\n")
        return _FakeProc(returncode=0, stdout="normal termination of xtb\n")
    monkeypatch.setattr(subprocess, "run", _fake_run)

    r = xtb.run_xtb_optimize_hessian(xyz, 2, 1, workdir)
    assert r.ran is True
    assert r.passed is True
    assert r.exit_code == 0
    assert r.opt_geometry_xyz.name == "xtbopt.xyz"
    assert r.hessian_path is not None
    assert r.charges_path is not None
    assert r.wbo_path is not None


def test_nonzero_exit_marks_failed(tmp_path, monkeypatch):
    import shutil as _sh
    monkeypatch.setattr(_sh, "which", lambda name: "/opt/xtb/bin/xtb")
    xyz = xtb.write_cluster_xyz([("Zn", 0, 0, 0)], tmp_path / "c.xyz")

    def _fake_run(cmd, **kw):
        return _FakeProc(returncode=1, stdout="", stderr="scf did not converge")
    monkeypatch.setattr(subprocess, "run", _fake_run)

    r = xtb.run_xtb_optimize_hessian(xyz, 2, 1, tmp_path / "work")
    assert r.ran is True
    assert r.passed is False
    assert r.exit_code == 1
    assert "nonzero" in r.fallback_reason


def test_success_but_no_output_file_marks_failed(tmp_path, monkeypatch):
    """xtb sometimes prints success but writes nothing (disk full etc.)."""
    import shutil as _sh
    monkeypatch.setattr(_sh, "which", lambda name: "/opt/xtb/bin/xtb")
    xyz = xtb.write_cluster_xyz([("Zn", 0, 0, 0)], tmp_path / "c.xyz")

    def _fake_run(cmd, **kw):
        return _FakeProc(returncode=0, stdout="normal termination\n")
    monkeypatch.setattr(subprocess, "run", _fake_run)

    r = xtb.run_xtb_optimize_hessian(xyz, 2, 1, tmp_path / "work")
    assert r.passed is False
    assert "xtbopt.xyz not written" in r.fallback_reason


def test_timeout_returns_ran_but_not_passed(tmp_path, monkeypatch):
    import shutil as _sh
    monkeypatch.setattr(_sh, "which", lambda name: "/opt/xtb/bin/xtb")
    xyz = xtb.write_cluster_xyz([("Zn", 0, 0, 0)], tmp_path / "c.xyz")

    def _fake_run(cmd, **kw):
        raise subprocess.TimeoutExpired(cmd=cmd, timeout=kw.get("timeout", 0))
    monkeypatch.setattr(subprocess, "run", _fake_run)

    r = xtb.run_xtb_optimize_hessian(xyz, 2, 1, tmp_path / "work", timeout_seconds=5)
    assert r.ran is True
    assert r.passed is False
    assert "timed out" in r.fallback_reason


def test_multiplicity_maps_to_uhf_correctly(tmp_path, monkeypatch):
    """Multiplicity 2S+1 → xtb --uhf = multiplicity - 1 = 2S."""
    import shutil as _sh
    monkeypatch.setattr(_sh, "which", lambda name: "/opt/xtb/bin/xtb")
    xyz = xtb.write_cluster_xyz([("Fe", 0, 0, 0)], tmp_path / "c.xyz")
    captured = {}

    def _fake_run(cmd, **kw):
        captured["cmd"] = cmd
        cwd = Path(kw["cwd"])
        (cwd / "xtbopt.xyz").write_text("1\nfake\nFe 0 0 0\n")
        return _FakeProc(returncode=0)
    monkeypatch.setattr(subprocess, "run", _fake_run)

    # High-spin Fe(II) = multiplicity 5, so uhf = 4
    xtb.run_xtb_optimize_hessian(xyz, 2, 5, tmp_path / "work")
    assert "--uhf" in captured["cmd"]
    idx = captured["cmd"].index("--uhf")
    assert captured["cmd"][idx + 1] == "4"


# ---------------------------------------------------------------------------
# emit_frcmod_scaffold — placeholder frcmod is well-formed enough for tleap
# ---------------------------------------------------------------------------


def test_frcmod_scaffold_writes_expected_sections(tmp_path):
    result = xtb.XtbRunResult(ran=True, passed=True)
    out = xtb.emit_frcmod_scaffold(
        result,
        metal_element="Zn",
        metal_atom_type="M1",
        donor_atom_types=["N1", "O1", "O2"],
        output_frcmod=tmp_path / "metal.frcmod",
    )
    text = out.read_text()
    for section in ("MASS", "BOND", "ANGLE", "DIHE", "IMPROPER", "NONBON"):
        assert section in text, f"frcmod missing section {section}"


def test_frcmod_scaffold_contains_todo_marker(tmp_path):
    """Reviewer transparency: the scaffold must announce itself as such."""
    result = xtb.XtbRunResult(ran=True, passed=True)
    out = xtb.emit_frcmod_scaffold(
        result, "Zn", "M1", ["N1", "O1"], tmp_path / "metal.frcmod",
    )
    text = out.read_text()
    assert "TODO(seminario)" in text
    assert "placeholder" in text.lower()


def test_frcmod_scaffold_emits_bond_per_donor(tmp_path):
    result = xtb.XtbRunResult(ran=True, passed=True)
    out = xtb.emit_frcmod_scaffold(
        result, "Zn", "M1", ["N1", "O1", "S1"], tmp_path / "metal.frcmod",
    )
    text = out.read_text()
    # Three BOND lines, one per donor
    bond_section = text.split("BOND", 1)[1].split("ANGLE", 1)[0]
    for donor in ("N1", "O1", "S1"):
        assert f"M1-{donor}" in bond_section


def test_frcmod_scaffold_emits_unique_angle_pairs(tmp_path):
    """N donors → C(N,2) angle terms, no duplicates."""
    result = xtb.XtbRunResult(ran=True, passed=True)
    out = xtb.emit_frcmod_scaffold(
        result, "Zn", "M1", ["N1", "O1", "S1"], tmp_path / "metal.frcmod",
    )
    text = out.read_text()
    angle_section = text.split("ANGLE", 1)[1].split("DIHE", 1)[0]
    angle_lines = [ln for ln in angle_section.splitlines() if "-M1-" in ln]
    # 3 donors → 3 unique unordered pairs
    assert len(angle_lines) == 3


def test_frcmod_scaffold_single_donor_has_no_angles(tmp_path):
    result = xtb.XtbRunResult(ran=True, passed=True)
    out = xtb.emit_frcmod_scaffold(
        result, "Zn", "M1", ["N1"], tmp_path / "metal.frcmod",
    )
    text = out.read_text()
    angle_section = text.split("ANGLE", 1)[1].split("DIHE", 1)[0]
    angle_lines = [ln for ln in angle_section.splitlines() if "-M1-" in ln]
    assert angle_lines == []


# ---------------------------------------------------------------------------
# XtbRunResult.to_dict round-trip
# ---------------------------------------------------------------------------


def test_run_result_to_dict_json_safe():
    import json
    r = xtb.XtbRunResult(
        ran=True, passed=True, exit_code=0,
        opt_geometry_xyz=Path("/tmp/xtbopt.xyz"),
        diagnostics=["cmd=xtb ..."],
    )
    payload = json.dumps(r.to_dict())
    reload = json.loads(payload)
    assert reload["ran"] is True
    assert reload["opt_geometry_xyz"] == "/tmp/xtbopt.xyz"
