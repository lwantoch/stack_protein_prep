"""Tests for scripts/plot_clashscore_histogram.py."""
from __future__ import annotations

import importlib.util
import json
import sys
from pathlib import Path

import pytest


_SCRIPT = (
    Path(__file__).resolve().parents[1]
    / "scripts" / "plot_clashscore_histogram.py"
)


def _load_mod():
    spec = importlib.util.spec_from_file_location("_plot_clash", _SCRIPT)
    mod = importlib.util.module_from_spec(spec)
    sys.modules["_plot_clash"] = mod
    spec.loader.exec_module(mod)
    return mod


def _write(path: Path, obj) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(obj))
    return path


def _rec(pdb: str, n_clash: int = 0, n_vdw: int = 0) -> dict:
    return {"pdb": pdb, "n_clash_pairs": n_clash, "n_vdw_clashes": n_vdw}


# ---------------------------------------------------------------------------
# _percentile
# ---------------------------------------------------------------------------


def test_percentile_boundaries():
    mod = _load_mod()
    assert mod._percentile([], 50) == 0.0
    assert mod._percentile([5], 50) == 5.0
    assert mod._percentile([1, 2, 3, 4, 5], 0) == 1.0
    assert mod._percentile([1, 2, 3, 4, 5], 100) == 5.0
    assert mod._percentile([1, 2, 3, 4, 5], 50) == 3.0
    # p90 of 10 items: 90th percentile lies between index 8 and 9 -> ~9.1
    vals = list(range(1, 11))
    assert 8.9 < mod._percentile(vals, 90) < 9.5


def test_percentile_matches_expected_p50_p90_p99():
    mod = _load_mod()
    vals = list(range(1, 101))  # 1..100
    assert mod._percentile(vals, 50) == pytest.approx(50.5, abs=0.5)
    assert mod._percentile(vals, 90) == pytest.approx(90.9, abs=1.0)
    assert mod._percentile(vals, 99) == pytest.approx(99.01, abs=1.0)


# ---------------------------------------------------------------------------
# _to_int helper
# ---------------------------------------------------------------------------


def test_to_int_handles_none_and_bad():
    mod = _load_mod()
    assert mod._to_int(None) is None
    assert mod._to_int(None, default=0) == 0
    assert mod._to_int("bad", default=0) == 0
    assert mod._to_int("5") == 5
    assert mod._to_int(3.7) == 3


# ---------------------------------------------------------------------------
# build_summary
# ---------------------------------------------------------------------------


def test_build_summary_reports_counts_and_percentiles():
    mod = _load_mod()
    cp = list(range(1, 11))
    vdw = [0, 1, 2, 3, 4]
    s = mod.build_summary(cp, vdw)
    assert "n_clash_pairs" in s
    assert "n_vdw_clashes" in s
    assert "n=10" in s
    assert "p50=" in s
    assert "p90=" in s


def test_build_summary_empty_data():
    mod = _load_mod()
    s = mod.build_summary([], [])
    assert "no data" in s


# ---------------------------------------------------------------------------
# _load_from_bench_json / _load_from_metrics_dir
# ---------------------------------------------------------------------------


def test_load_from_bench_json_reads_list(tmp_path: Path):
    mod = _load_mod()
    p = _write(tmp_path / "b.json", [_rec("A", 3, 1), _rec("B", 5, 2)])
    r = mod._load_from_bench_json(p)
    assert len(r) == 2
    assert r[0]["pdb"] == "A"


def test_load_from_bench_json_rejects_non_list(tmp_path: Path):
    mod = _load_mod()
    p = tmp_path / "b.json"
    p.write_text(json.dumps({"pdb": "A"}))
    with pytest.raises(ValueError, match="not a JSON list"):
        mod._load_from_bench_json(p)


def test_load_from_metrics_dir_pulls_json_files(tmp_path: Path):
    mod = _load_mod()
    _write(tmp_path / "8ABC.json", {"n_clash_pairs": 3})
    _write(tmp_path / "8XYZ.json", {"n_clash_pairs": 5, "n_vdw_clashes": 1})
    r = mod._load_from_metrics_dir(tmp_path)
    assert len(r) == 2
    # pdb backfilled from stem
    pdbs = {rec["pdb"] for rec in r}
    assert pdbs == {"8ABC", "8XYZ"}


def test_load_from_metrics_dir_raises_when_missing(tmp_path: Path):
    mod = _load_mod()
    with pytest.raises(FileNotFoundError):
        mod._load_from_metrics_dir(tmp_path / "nope")


# ---------------------------------------------------------------------------
# CLI end-to-end
# ---------------------------------------------------------------------------


def test_end_to_end_bench_json_writes_png_csv_summary(tmp_path: Path):
    mod = _load_mod()
    p = _write(tmp_path / "b.json", [
        _rec("A", 2, 0), _rec("B", 10, 3), _rec("C", 50, 20),
    ])
    outdir = tmp_path / "figs"
    rc = mod.main([
        "--bench-json", str(p),
        "--outdir", str(outdir),
    ])
    assert rc == 0
    assert (outdir / "clashscore_histogram.png").is_file()
    assert (outdir / "clashscore_histogram.csv").is_file()
    assert (outdir / "summary.txt").is_file()

    csv_text = (outdir / "clashscore_histogram.csv").read_text()
    assert "pdb_id" in csv_text.splitlines()[0]
    assert csv_text.count("\n") >= 4  # header + 3 data rows


def test_end_to_end_metrics_dir_reads_all(tmp_path: Path):
    mod = _load_mod()
    mdir = tmp_path / "metrics"
    _write(mdir / "8ABC.json", {"n_clash_pairs": 2, "n_vdw_clashes": 0})
    _write(mdir / "8XYZ.json", {"n_clash_pairs": 15, "n_vdw_clashes": 4})
    outdir = tmp_path / "figs"
    rc = mod.main([
        "--pdb-metrics-dir", str(mdir),
        "--outdir", str(outdir),
    ])
    assert rc == 0
    text = (outdir / "clashscore_histogram.csv").read_text()
    assert "8ABC" in text
    assert "8XYZ" in text


def test_missing_bench_json_returns_nonzero(tmp_path: Path):
    mod = _load_mod()
    outdir = tmp_path / "figs"
    rc = mod.main([
        "--bench-json", str(tmp_path / "nope.json"),
        "--outdir", str(outdir),
    ])
    assert rc == 2


def test_missing_metrics_dir_returns_nonzero(tmp_path: Path):
    mod = _load_mod()
    outdir = tmp_path / "figs"
    rc = mod.main([
        "--pdb-metrics-dir", str(tmp_path / "nope"),
        "--outdir", str(outdir),
    ])
    assert rc == 2


def test_bench_json_and_metrics_dir_are_mutually_exclusive(tmp_path: Path):
    mod = _load_mod()
    p = _write(tmp_path / "b.json", [_rec("A")])
    with pytest.raises(SystemExit):  # argparse aborts
        mod.main([
            "--bench-json", str(p),
            "--pdb-metrics-dir", str(tmp_path),
            "--outdir", str(tmp_path / "figs"),
        ])


def test_log_y_and_bin_width(tmp_path: Path):
    mod = _load_mod()
    p = _write(tmp_path / "b.json", [
        _rec("A", 2, 0), _rec("B", 8, 2),
    ])
    outdir = tmp_path / "figs"
    rc = mod.main([
        "--bench-json", str(p),
        "--outdir", str(outdir),
        "--log-y",
        "--bin-width", "2",
        "--x-max", "30",
    ])
    assert rc == 0
    assert (outdir / "clashscore_histogram.png").is_file()


def test_falls_back_to_clash_field_when_n_clash_pairs_missing(tmp_path: Path):
    """Legacy bench JSONs use 'clash' — the loader should accept that."""
    mod = _load_mod()
    p = _write(tmp_path / "b.json", [
        {"pdb": "A", "clash": 7},
        {"pdb": "B", "clash": 22},
    ])
    outdir = tmp_path / "figs"
    rc = mod.main([
        "--bench-json", str(p),
        "--outdir", str(outdir),
    ])
    assert rc == 0
    summary = (outdir / "summary.txt").read_text()
    assert "n=2" in summary
