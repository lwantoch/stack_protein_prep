"""Smoke tests for scripts/plot_omega_planarity_distribution.py.

Verifies the CLI runs end-to-end on a tiny toy PDB set: emits PNG,
CSV, and per_pdb_summary.txt.  Matplotlib is used with the Agg backend
so no display is needed.
"""
from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import pytest


_SCRIPT = (
    Path(__file__).resolve().parents[1]
    / "scripts" / "plot_omega_planarity_distribution.py"
)


def _load_script_module():
    spec = importlib.util.spec_from_file_location("_plot_omega_cli", _SCRIPT)
    mod = importlib.util.module_from_spec(spec)
    sys.modules["_plot_omega_cli"] = mod
    spec.loader.exec_module(mod)
    return mod


def _write_trans_pdb(path: Path, resname_j: str = "ALA"):
    path.parent.mkdir(parents=True, exist_ok=True)
    text = (
        "ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00  0.00           N\n"
        "ATOM      2  CA  ALA A   1       1.458   0.000   0.000  1.00  0.00           C\n"
        "ATOM      3  C   ALA A   1       2.009   1.420   0.000  1.00  0.00           C\n"
        "ATOM      4  O   ALA A   1       1.251   2.390   0.000  1.00  0.00           O\n"
        f"ATOM      5  N   {resname_j} A   2       3.332   1.540   0.000  1.00  0.00           N\n"
        f"ATOM      6  CA  {resname_j} A   2       4.000   2.828   0.000  1.00  0.00           C\n"
        f"ATOM      7  C   {resname_j} A   2       5.500   2.700   0.000  1.00  0.00           C\n"
        f"ATOM      8  O   {resname_j} A   2       6.000   1.600   0.000  1.00  0.00           O\n"
        "TER\nEND\n"
    )
    path.write_text(text)


def test_end_to_end_writes_all_outputs(tmp_path: Path):
    bench = tmp_path / "bench"
    (bench / "8ABC").mkdir(parents=True)
    (bench / "8XYZ").mkdir(parents=True)
    _write_trans_pdb(bench / "8ABC" / "final_model.pdb")
    _write_trans_pdb(bench / "8XYZ" / "final_model.pdb", resname_j="GLY")

    mod = _load_script_module()
    outdir = tmp_path / "figs"
    argv = [
        "--pdb-glob", str(bench / "*/final_model.pdb"),
        "--outdir", str(outdir),
    ]
    rc = mod.main(argv)
    assert rc == 0
    assert (outdir / "omega_distribution.csv").is_file()
    assert (outdir / "omega_distribution.png").is_file()
    assert (outdir / "per_pdb_summary.txt").is_file()

    # CSV has 2 rows (one pair per PDB) + header
    lines = (outdir / "omega_distribution.csv").read_text().strip().splitlines()
    assert len(lines) == 3  # header + 2 data rows
    assert "pdb_id" in lines[0]
    assert "8ABC" in lines[1] or "8ABC" in lines[2]
    assert "8XYZ" in lines[1] or "8XYZ" in lines[2]


def test_no_matches_returns_nonzero(tmp_path: Path):
    mod = _load_script_module()
    outdir = tmp_path / "figs"
    argv = [
        "--pdb-glob", str(tmp_path / "nonexistent_*.pdb"),
        "--outdir", str(outdir),
    ]
    rc = mod.main(argv)
    assert rc == 2


def test_pdb_id_from_stem_mode(tmp_path: Path):
    bench = tmp_path / "bench"
    bench.mkdir()
    _write_trans_pdb(bench / "9WXY.pdb")
    mod = _load_script_module()
    outdir = tmp_path / "figs"
    argv = [
        "--pdb-glob", str(bench / "*.pdb"),
        "--outdir", str(outdir),
        "--pdb-id-from", "stem",
    ]
    rc = mod.main(argv)
    assert rc == 0
    csv_text = (outdir / "omega_distribution.csv").read_text()
    assert "9WXY" in csv_text


def test_log_y_option(tmp_path: Path):
    """--log-y must not crash even when a bin is zero."""
    bench = tmp_path / "bench"
    (bench / "8ABC").mkdir(parents=True)
    _write_trans_pdb(bench / "8ABC" / "final_model.pdb")
    mod = _load_script_module()
    outdir = tmp_path / "figs"
    argv = [
        "--pdb-glob", str(bench / "*/final_model.pdb"),
        "--outdir", str(outdir),
        "--log-y",
    ]
    rc = mod.main(argv)
    assert rc == 0
    assert (outdir / "omega_distribution.png").is_file()


def test_bin_width_argument(tmp_path: Path):
    bench = tmp_path / "bench"
    (bench / "8ABC").mkdir(parents=True)
    _write_trans_pdb(bench / "8ABC" / "final_model.pdb")
    mod = _load_script_module()
    outdir = tmp_path / "figs"
    argv = [
        "--pdb-glob", str(bench / "*/final_model.pdb"),
        "--outdir", str(outdir),
        "--bin-width", "5.0",
    ]
    rc = mod.main(argv)
    assert rc == 0
