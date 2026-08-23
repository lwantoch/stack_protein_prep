"""Smoke tests for scripts/plot_metal_coord_preservation.py."""
from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import pytest


_SCRIPT = (
    Path(__file__).resolve().parents[1]
    / "scripts" / "plot_metal_coord_preservation.py"
)


def _load_mod():
    spec = importlib.util.spec_from_file_location("_plot_metal_cli", _SCRIPT)
    mod = importlib.util.module_from_spec(spec)
    sys.modules["_plot_metal_cli"] = mod
    spec.loader.exec_module(mod)
    return mod


def _write_zn_his(path: Path, distance: float) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        "ATOM      1  N   HIS A   1       0.000   0.000   3.000  1.00  0.00           N\n"
        "ATOM      2  CA  HIS A   1       0.000   0.000   4.000  1.00  0.00           C\n"
        "ATOM      3  C   HIS A   1       0.000   0.000   5.000  1.00  0.00           C\n"
        "ATOM      4  O   HIS A   1       1.000   0.000   5.000  1.00  0.00           O\n"
        "ATOM      5  CB  HIS A   1       0.000   1.000   4.000  1.00  0.00           C\n"
        "ATOM      6  ND1 HIS A   1       0.000   2.000   4.000  1.00  0.00           N\n"
        "ATOM      7  NE2 HIS A   1       0.000   0.000   0.000  1.00  0.00           N\n"
        f"HETATM    8 ZN    ZN B 101    {distance: >8.3f}   0.000   0.000  1.00  0.00          ZN\n"
        "END\n"
    )
    return path


def test_pair_parser_ok():
    mod = _load_mod()
    label, cry, fru = mod._parse_pair("8ABC=/tmp/c.pdb:/tmp/f.pdb")
    assert label == "8ABC"
    assert cry == Path("/tmp/c.pdb")
    assert fru == Path("/tmp/f.pdb")


def test_pair_parser_rejects_missing_colon():
    mod = _load_mod()
    with pytest.raises(Exception):
        mod._parse_pair("8ABC=/tmp/onlyone.pdb")


def test_end_to_end_writes_png_csv_summary(tmp_path: Path):
    crystal = _write_zn_his(tmp_path / "crystal" / "8ABC.pdb", distance=2.05)
    fruton = _write_zn_his(tmp_path / "fruton" / "8ABC.pdb", distance=2.30)
    mod = _load_mod()
    outdir = tmp_path / "figs"
    rc = mod.main([
        "--pair", f"8ABC={crystal}:{fruton}",
        "--outdir", str(outdir),
    ])
    assert rc == 0
    assert (outdir / "metal_coord_preservation.png").is_file()
    assert (outdir / "metal_coord_preservation.csv").is_file()
    assert (outdir / "per_pdb_summary.txt").is_file()

    csv_text = (outdir / "metal_coord_preservation.csv").read_text()
    assert "8ABC" in csv_text
    assert "NE2" in csv_text
    # delta = 2.30 - 2.05 ≈ +0.25
    assert "+0.25" in csv_text or "+0.250" in csv_text


def test_lost_donor_flagged(tmp_path: Path):
    """FRUTON model has Zn moved beyond cutoff → donor 'lost'."""
    crystal = _write_zn_his(tmp_path / "c.pdb", distance=2.05)
    fruton = _write_zn_his(tmp_path / "f.pdb", distance=5.00)
    mod = _load_mod()
    outdir = tmp_path / "figs"
    rc = mod.main([
        "--pair", f"8ABC={crystal}:{fruton}",
        "--outdir", str(outdir),
    ])
    assert rc == 0
    csv_text = (outdir / "metal_coord_preservation.csv").read_text()
    assert "lost" in csv_text


def test_missing_crystal_recorded_as_skip(tmp_path: Path):
    fruton = _write_zn_his(tmp_path / "f.pdb", distance=2.10)
    mod = _load_mod()
    outdir = tmp_path / "figs"
    rc = mod.main([
        "--pair", f"8ABC={tmp_path/'nope.pdb'}:{fruton}",
        "--outdir", str(outdir),
    ])
    assert rc == 0
    summary = (outdir / "per_pdb_summary.txt").read_text()
    assert "SKIP" in summary


def test_multiple_pairs_all_written(tmp_path: Path):
    c1 = _write_zn_his(tmp_path / "c1.pdb", 2.05)
    f1 = _write_zn_his(tmp_path / "f1.pdb", 2.15)
    c2 = _write_zn_his(tmp_path / "c2.pdb", 2.10)
    f2 = _write_zn_his(tmp_path / "f2.pdb", 2.35)
    mod = _load_mod()
    outdir = tmp_path / "figs"
    rc = mod.main([
        "--pair", f"A={c1}:{f1}",
        "--pair", f"B={c2}:{f2}",
        "--outdir", str(outdir),
    ])
    assert rc == 0
    csv_text = (outdir / "metal_coord_preservation.csv").read_text()
    assert "A" in csv_text and "B" in csv_text
    assert csv_text.count("preserved") >= 2
