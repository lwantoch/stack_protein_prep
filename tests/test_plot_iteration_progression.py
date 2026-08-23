"""Tests for scripts/plot_iteration_progression.py."""
from __future__ import annotations

import importlib.util
import json
import sys
from pathlib import Path

import pytest


_SCRIPT = (
    Path(__file__).resolve().parents[1]
    / "scripts" / "plot_iteration_progression.py"
)


def _load():
    spec = importlib.util.spec_from_file_location("_plip", _SCRIPT)
    mod = importlib.util.module_from_spec(spec)
    sys.modules["_plip"] = mod
    spec.loader.exec_module(mod)
    return mod


def _write(path: Path, obj) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(obj))
    return path


# ---------------------------------------------------------------------------
# load_history
# ---------------------------------------------------------------------------


def test_load_history_reads_list(tmp_path: Path):
    mod = _load()
    p = _write(tmp_path / "h.json", [
        {"iter_label": "iter-2", "n_pdbs": 48, "gate_pass_pct": 89.6},
        {"iter_label": "iter-3", "n_pdbs": 46, "gate_pass_pct": 91.3,
         "n_omega_non_planar": 6},
    ])
    r = mod.load_history(p)
    assert len(r) == 2
    assert r[0]["iter_label"] == "iter-2"


def test_load_history_rejects_non_list(tmp_path: Path):
    mod = _load()
    p = tmp_path / "h.json"
    p.write_text(json.dumps({"iter_label": "iter-2"}))
    with pytest.raises(ValueError, match="not a JSON list"):
        mod.load_history(p)


def test_load_history_rejects_missing_iter_label(tmp_path: Path):
    mod = _load()
    p = _write(tmp_path / "h.json", [{"gate_pass_pct": 89.6}])
    with pytest.raises(ValueError, match="missing iter_label"):
        mod.load_history(p)


def test_load_history_rejects_missing_gate_pass_pct(tmp_path: Path):
    mod = _load()
    p = _write(tmp_path / "h.json", [{"iter_label": "iter-2"}])
    with pytest.raises(ValueError, match="missing gate_pass_pct"):
        mod.load_history(p)


# ---------------------------------------------------------------------------
# CLI end-to-end
# ---------------------------------------------------------------------------


def test_end_to_end_writes_png_csv(tmp_path: Path):
    mod = _load()
    p = _write(tmp_path / "h.json", [
        {"iter_label": "iter-2", "n_pdbs": 48, "gate_pass_pct": 89.6,
         "n_omega_non_planar": None, "notes": "baseline"},
        {"iter_label": "iter-3", "n_pdbs": 46, "gate_pass_pct": 91.3,
         "n_omega_non_planar": 6, "notes": "adaptive fast->slow"},
        {"iter_label": "iter-4", "n_pdbs": 46, "gate_pass_pct": 97.8,
         "n_omega_non_planar": 3, "notes": "ceiling guard"},
    ])
    outdir = tmp_path / "out"
    rc = mod.main([
        "--history", str(p),
        "--outdir", str(outdir),
    ])
    assert rc == 0
    assert (outdir / "iteration_progression.png").is_file()
    assert (outdir / "iteration_progression.csv").is_file()

    csv_text = (outdir / "iteration_progression.csv").read_text()
    assert "iter-2" in csv_text
    assert "iter-3" in csv_text
    assert "iter-4" in csv_text
    assert "89.6" in csv_text
    assert "97.8" in csv_text


def test_end_to_end_missing_history_returns_nonzero(tmp_path: Path):
    mod = _load()
    rc = mod.main([
        "--history", str(tmp_path / "nope.json"),
        "--outdir", str(tmp_path / "out"),
    ])
    assert rc == 2


def test_end_to_end_empty_history_returns_nonzero(tmp_path: Path):
    mod = _load()
    p = _write(tmp_path / "h.json", [])
    rc = mod.main([
        "--history", str(p),
        "--outdir", str(tmp_path / "out"),
    ])
    assert rc == 2


def test_end_to_end_null_omega_still_plots(tmp_path: Path):
    """Iterations with null n_omega_non_planar must not break the plot."""
    mod = _load()
    p = _write(tmp_path / "h.json", [
        {"iter_label": "A", "gate_pass_pct": 50.0, "n_omega_non_planar": None},
        {"iter_label": "B", "gate_pass_pct": 80.0, "n_omega_non_planar": 10},
        {"iter_label": "C", "gate_pass_pct": 95.0, "n_omega_non_planar": None},
    ])
    outdir = tmp_path / "out"
    rc = mod.main([
        "--history", str(p),
        "--outdir", str(outdir),
    ])
    assert rc == 0
    assert (outdir / "iteration_progression.png").is_file()


def test_end_to_end_all_null_omega(tmp_path: Path):
    """No ω data at all → right panel shows placeholder, left works."""
    mod = _load()
    p = _write(tmp_path / "h.json", [
        {"iter_label": "A", "gate_pass_pct": 50.0, "n_omega_non_planar": None},
        {"iter_label": "B", "gate_pass_pct": 80.0, "n_omega_non_planar": None},
    ])
    outdir = tmp_path / "out"
    rc = mod.main([
        "--history", str(p),
        "--outdir", str(outdir),
    ])
    assert rc == 0
    assert (outdir / "iteration_progression.png").is_file()


# ---------------------------------------------------------------------------
# Ships iteration_history_20260823.json — must be valid + parseable.
# ---------------------------------------------------------------------------


def test_shipped_history_json_loads(tmp_path: Path):
    """The shipped iteration_history_20260823.json must be valid."""
    mod = _load()
    shipped = (
        Path(__file__).resolve().parents[1]
        / "src" / "stack_protein_preparation" / "data"
        / "iteration_history_20260823.json"
    )
    assert shipped.is_file(), f"shipped history missing: {shipped}"
    history = mod.load_history(shipped)
    assert len(history) >= 3, "expect at least iter-2/3/4"
    labels = [e["iter_label"] for e in history]
    assert "iter-3" in labels
    assert "iter-4" in labels
    # gate_pass_pct rises monotonically across iter-2 -> iter-4
    numeric = [e["gate_pass_pct"] for e in history if e.get("n_pdbs")]
    assert numeric[-1] > numeric[0], (
        "gate_pass_pct should rise across iterations; got " + str(numeric)
    )


def test_shipped_history_json_produces_figure(tmp_path: Path):
    """End-to-end: render the killer paper figure from the shipped history."""
    mod = _load()
    shipped = (
        Path(__file__).resolve().parents[1]
        / "src" / "stack_protein_preparation" / "data"
        / "iteration_history_20260823.json"
    )
    outdir = tmp_path / "figures"
    rc = mod.main([
        "--history", str(shipped),
        "--outdir", str(outdir),
    ])
    assert rc == 0
    assert (outdir / "iteration_progression.png").is_file()
