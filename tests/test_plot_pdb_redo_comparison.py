"""Smoke tests for scripts/plot_pdb_redo_comparison.py."""
from __future__ import annotations

import importlib.util
import json
import sys
from pathlib import Path


_SCRIPT = (
    Path(__file__).resolve().parents[1]
    / "scripts" / "plot_pdb_redo_comparison.py"
)


def _load_mod():
    spec = importlib.util.spec_from_file_location("_plot_pdb_redo", _SCRIPT)
    mod = importlib.util.module_from_spec(spec)
    sys.modules["_plot_pdb_redo"] = mod
    spec.loader.exec_module(mod)
    return mod


def _write(path: Path, obj: dict) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(obj))
    return path


def test_end_to_end_writes_png_md_json(tmp_path: Path):
    frdir = tmp_path / "fruton"
    prdir = tmp_path / "pdb_redo"
    _write(frdir / "8ABC.json", {"n_rama_outlier": 3, "clashscore": 5.0})
    _write(frdir / "8XYZ.json", {"n_rama_outlier": 8, "clashscore": 10.0})
    _write(prdir / "8ABC.json", {"n_rama_outlier": 7, "clashscore": 12.0})
    _write(prdir / "8XYZ.json", {"n_rama_outlier": 4, "clashscore": 3.0})

    mod = _load_mod()
    outdir = tmp_path / "figs"
    rc = mod.main([
        "--fruton-dir", str(frdir),
        "--pdb-redo-dir", str(prdir),
        "--outdir", str(outdir),
        "--metric", "n_rama_outlier",
        "--metric", "clashscore",
    ])
    assert rc == 0
    assert (outdir / "pdb_redo_comparison.png").is_file()
    assert (outdir / "pdb_redo_comparison.md").is_file()
    assert (outdir / "pdb_redo_comparison.json").is_file()

    md = (outdir / "pdb_redo_comparison.md").read_text()
    assert "Ramachandran outliers" in md
    assert "MolProbity clashscore" in md


def test_default_metric_selection_uses_all(tmp_path: Path):
    frdir = tmp_path / "fruton"
    prdir = tmp_path / "pdb_redo"
    _write(frdir / "8ABC.json", {"n_rama_outlier": 1})
    _write(prdir / "8ABC.json", {"n_rama_outlier": 1})
    mod = _load_mod()
    outdir = tmp_path / "figs"
    rc = mod.main([
        "--fruton-dir", str(frdir),
        "--pdb-redo-dir", str(prdir),
        "--outdir", str(outdir),
    ])
    assert rc == 0
    js = json.loads((outdir / "pdb_redo_comparison.json").read_text())
    # Nine metrics in the catalogue.
    assert len(js["per_metric"]) == 9


def test_missing_directory_returns_nonzero(tmp_path: Path):
    mod = _load_mod()
    outdir = tmp_path / "figs"
    rc = mod.main([
        "--fruton-dir", str(tmp_path / "nope"),
        "--pdb-redo-dir", str(tmp_path),
        "--outdir", str(outdir),
    ])
    assert rc == 2


def test_disjoint_pdbs_produces_empty_plot(tmp_path: Path):
    frdir = tmp_path / "f"
    prdir = tmp_path / "p"
    _write(frdir / "8ABC.json", {"n_rama_outlier": 3})
    _write(prdir / "8XYZ.json", {"n_rama_outlier": 7})
    mod = _load_mod()
    outdir = tmp_path / "figs"
    rc = mod.main([
        "--fruton-dir", str(frdir),
        "--pdb-redo-dir", str(prdir),
        "--outdir", str(outdir),
        "--metric", "n_rama_outlier",
    ])
    assert rc == 0
    js = json.loads((outdir / "pdb_redo_comparison.json").read_text())
    assert js["per_pdb_deltas"] == {}
    # The metric row exists but n_matched=0.
    assert js["per_metric"][0]["n_matched"] == 0
    assert (outdir / "pdb_redo_comparison.png").is_file()


def test_unknown_metric_returns_nonzero(tmp_path: Path):
    frdir = tmp_path / "f"
    prdir = tmp_path / "p"
    _write(frdir / "8ABC.json", {})
    _write(prdir / "8ABC.json", {})
    mod = _load_mod()
    outdir = tmp_path / "figs"
    rc = mod.main([
        "--fruton-dir", str(frdir),
        "--pdb-redo-dir", str(prdir),
        "--outdir", str(outdir),
        "--metric", "not_a_real_metric",
    ])
    assert rc == 2
