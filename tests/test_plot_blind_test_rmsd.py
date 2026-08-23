"""Smoke tests for scripts/plot_blind_test_rmsd.py."""
from __future__ import annotations

import importlib.util
import json
import sys
from pathlib import Path


_SCRIPT = (
    Path(__file__).resolve().parents[1]
    / "scripts" / "plot_blind_test_rmsd.py"
)


def _load_mod():
    spec = importlib.util.spec_from_file_location("_plot_blind_test", _SCRIPT)
    mod = importlib.util.module_from_spec(spec)
    sys.modules["_plot_blind_test"] = mod
    spec.loader.exec_module(mod)
    return mod


def _write_crystal(path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        "ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00  0.00           N\n"
        "ATOM      2  CA  ALA A   1       1.458   0.000   0.000  1.00  0.00           C\n"
        "ATOM      3  C   ALA A   1       2.009   1.420   0.000  1.00  0.00           C\n"
        "ATOM      4  O   ALA A   1       1.251   2.390   0.000  1.00  0.00           O\n"
        "ATOM      5  N   ALA A   2       3.332   1.540   0.000  1.00  0.00           N\n"
        "ATOM      6  CA  ALA A   2       4.000   2.828   0.000  1.00  0.00           C\n"
        "ATOM      7  C   ALA A   2       5.500   2.700   0.000  1.00  0.00           C\n"
        "ATOM      8  O   ALA A   2       6.000   1.600   0.000  1.00  0.00           O\n"
        "TER\nEND\n"
    )


def _write_filled(path: Path, ca2_dx: float = 0.0) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    ca2_x = 4.000 + ca2_dx
    path.write_text(
        "ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00  0.00           N\n"
        "ATOM      2  CA  ALA A   1       1.458   0.000   0.000  1.00  0.00           C\n"
        "ATOM      3  C   ALA A   1       2.009   1.420   0.000  1.00  0.00           C\n"
        "ATOM      4  O   ALA A   1       1.251   2.390   0.000  1.00  0.00           O\n"
        "ATOM      5  N   ALA A   2       3.332   1.540   0.000  1.00  0.00           N\n"
        f"ATOM      6  CA  ALA A   2     {ca2_x:7.3f}   2.828   0.000  1.00  0.00           C\n"
        "ATOM      7  C   ALA A   2       5.500   2.700   0.000  1.00  0.00           C\n"
        "ATOM      8  O   ALA A   2       6.000   1.600   0.000  1.00  0.00           O\n"
        "TER\nEND\n"
    )


def test_end_to_end_writes_png_csv_summary(tmp_path: Path):
    bench = tmp_path / "bench"
    (bench / "8ABC").mkdir(parents=True)
    (bench / "8XYZ").mkdir(parents=True)
    _write_crystal(bench / "8ABC" / "crystal.pdb")
    _write_filled(bench / "8ABC" / "filled.pdb", ca2_dx=0.4)
    _write_crystal(bench / "8XYZ" / "crystal.pdb")
    _write_filled(bench / "8XYZ" / "filled.pdb", ca2_dx=1.2)

    spec = tmp_path / "spec.json"
    spec.write_text(json.dumps([
        {
            "pdb_id": "8ABC",
            "crystal_pdb": str(bench / "8ABC" / "crystal.pdb"),
            "filled_pdb": str(bench / "8ABC" / "filled.pdb"),
            "held_out": [{"chain": "A", "first_resnum": 2, "last_resnum": 2}],
        },
        {
            "pdb_id": "8XYZ",
            "crystal_pdb": str(bench / "8XYZ" / "crystal.pdb"),
            "filled_pdb": str(bench / "8XYZ" / "filled.pdb"),
            "held_out": [{"chain": "A", "first_resnum": 2, "last_resnum": 2}],
        },
    ]))

    mod = _load_mod()
    outdir = tmp_path / "figs"
    rc = mod.main([
        "--bench-spec", str(spec),
        "--outdir", str(outdir),
    ])
    assert rc == 0
    assert (outdir / "blind_test_rmsd.png").is_file()
    assert (outdir / "blind_test_rmsd.csv").is_file()
    assert (outdir / "per_pdb_summary.txt").is_file()

    csv_lines = (outdir / "blind_test_rmsd.csv").read_text().strip().splitlines()
    assert len(csv_lines) == 3  # header + 2 rows
    assert "pdb_id" in csv_lines[0]


def test_missing_spec_returns_nonzero(tmp_path: Path):
    mod = _load_mod()
    outdir = tmp_path / "figs"
    rc = mod.main([
        "--bench-spec", str(tmp_path / "nope.json"),
        "--outdir", str(outdir),
    ])
    assert rc == 2


def test_spec_not_a_list_returns_nonzero(tmp_path: Path):
    spec = tmp_path / "s.json"
    spec.write_text(json.dumps({"pdb_id": "8ABC"}))
    mod = _load_mod()
    outdir = tmp_path / "figs"
    rc = mod.main([
        "--bench-spec", str(spec),
        "--outdir", str(outdir),
    ])
    assert rc == 2


def test_entry_missing_held_out_is_skipped(tmp_path: Path):
    bench = tmp_path / "bench"
    (bench / "8ABC").mkdir(parents=True)
    _write_crystal(bench / "8ABC" / "crystal.pdb")
    _write_filled(bench / "8ABC" / "filled.pdb")
    spec = tmp_path / "spec.json"
    spec.write_text(json.dumps([
        {
            "pdb_id": "8ABC",
            "crystal_pdb": str(bench / "8ABC" / "crystal.pdb"),
            "filled_pdb": str(bench / "8ABC" / "filled.pdb"),
            "held_out": [],
        },
    ]))
    mod = _load_mod()
    outdir = tmp_path / "figs"
    rc = mod.main([
        "--bench-spec", str(spec),
        "--outdir", str(outdir),
    ])
    assert rc == 0
    summary = (outdir / "per_pdb_summary.txt").read_text()
    assert "SKIP" in summary
    assert "no held_out ranges" in summary


def test_missing_files_in_entry_are_skipped_gracefully(tmp_path: Path):
    spec = tmp_path / "spec.json"
    spec.write_text(json.dumps([
        {
            "pdb_id": "8ABC",
            "crystal_pdb": str(tmp_path / "nope1.pdb"),
            "filled_pdb": str(tmp_path / "nope2.pdb"),
            "held_out": [{"chain": "A", "first_resnum": 1, "last_resnum": 5}],
        },
    ]))
    mod = _load_mod()
    outdir = tmp_path / "figs"
    rc = mod.main([
        "--bench-spec", str(spec),
        "--outdir", str(outdir),
    ])
    assert rc == 0  # graceful skip, not error
    summary = (outdir / "per_pdb_summary.txt").read_text()
    assert "SKIP" in summary


def test_log_y_and_bin_width_smoke(tmp_path: Path):
    bench = tmp_path / "bench"
    (bench / "8ABC").mkdir(parents=True)
    _write_crystal(bench / "8ABC" / "crystal.pdb")
    _write_filled(bench / "8ABC" / "filled.pdb", ca2_dx=0.3)
    spec = tmp_path / "spec.json"
    spec.write_text(json.dumps([
        {
            "pdb_id": "8ABC",
            "crystal_pdb": str(bench / "8ABC" / "crystal.pdb"),
            "filled_pdb": str(bench / "8ABC" / "filled.pdb"),
            "held_out": [{"chain": "A", "first_resnum": 2, "last_resnum": 2}],
        },
    ]))
    mod = _load_mod()
    outdir = tmp_path / "figs"
    rc = mod.main([
        "--bench-spec", str(spec),
        "--outdir", str(outdir),
        "--log-y",
        "--bin-width", "0.5",
        "--x-max", "5.0",
    ])
    assert rc == 0
    assert (outdir / "blind_test_rmsd.png").is_file()
