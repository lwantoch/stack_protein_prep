"""Smoke tests for scripts/plot_family_stratification.py."""
from __future__ import annotations

import importlib.util
import json
import sys
from pathlib import Path


_SCRIPT = (
    Path(__file__).resolve().parents[1]
    / "scripts" / "plot_family_stratification.py"
)


def _load_mod():
    spec = importlib.util.spec_from_file_location("_plot_family", _SCRIPT)
    mod = importlib.util.module_from_spec(spec)
    sys.modules["_plot_family"] = mod
    spec.loader.exec_module(mod)
    return mod


def _write(path: Path, obj) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(obj))
    return path


def _rec(pdb: str, gate_pass: bool = True, delta_n: int = 0) -> dict:
    return {"pdb": pdb, "gate_pass": gate_pass, "delta_n": delta_n,
            "clash": 0, "brk": 0}


def test_end_to_end_writes_png_md_json(tmp_path: Path):
    bench = _write(tmp_path / "bench.json", [
        _rec("8ABC", True, 10),
        _rec("8XYZ", True, 5),
        _rec("8QQQ", False, 0),
    ])
    fam = _write(tmp_path / "family.json", {
        "8ABC": "kinase", "8XYZ": "kinase", "8QQQ": "gpcr",
    })
    mod = _load_mod()
    outdir = tmp_path / "figs"
    rc = mod.main([
        "--bench-json", str(bench),
        "--family-map", str(fam),
        "--outdir", str(outdir),
    ])
    assert rc == 0
    assert (outdir / "family_stratification.png").is_file()
    assert (outdir / "family_stratification.md").is_file()
    assert (outdir / "family_stratification.json").is_file()

    md = (outdir / "family_stratification.md").read_text()
    assert "kinase" in md
    assert "gpcr" in md


def test_missing_bench_file_returns_nonzero(tmp_path: Path):
    fam = _write(tmp_path / "f.json", {"8ABC": "kinase"})
    mod = _load_mod()
    outdir = tmp_path / "figs"
    rc = mod.main([
        "--bench-json", str(tmp_path / "nope.json"),
        "--family-map", str(fam),
        "--outdir", str(outdir),
    ])
    assert rc == 2


def test_missing_family_map_returns_nonzero(tmp_path: Path):
    bench = _write(tmp_path / "bench.json", [_rec("8ABC")])
    mod = _load_mod()
    outdir = tmp_path / "figs"
    rc = mod.main([
        "--bench-json", str(bench),
        "--family-map", str(tmp_path / "nope.json"),
        "--outdir", str(outdir),
    ])
    assert rc == 2


def test_bench_json_not_a_list_returns_nonzero(tmp_path: Path):
    bench = _write(tmp_path / "bench.json", {"pdb": "8ABC"})  # dict, not list
    fam = _write(tmp_path / "f.json", {"8ABC": "kinase"})
    mod = _load_mod()
    outdir = tmp_path / "figs"
    rc = mod.main([
        "--bench-json", str(bench),
        "--family-map", str(fam),
        "--outdir", str(outdir),
    ])
    assert rc == 2


def test_unassigned_pdbs_appear_in_output(tmp_path: Path):
    bench = _write(tmp_path / "bench.json", [
        _rec("8ABC", True), _rec("8XYZ", False),
    ])
    fam = _write(tmp_path / "f.json", {"8ABC": "kinase"})  # 8XYZ unmapped
    mod = _load_mod()
    outdir = tmp_path / "figs"
    rc = mod.main([
        "--bench-json", str(bench),
        "--family-map", str(fam),
        "--outdir", str(outdir),
    ])
    assert rc == 0
    md = (outdir / "family_stratification.md").read_text()
    assert "__unassigned__" in md


def test_output_json_contains_overall_aggregate(tmp_path: Path):
    bench = _write(tmp_path / "bench.json", [
        _rec("A", True, 10), _rec("B", False, 0),
    ])
    fam = _write(tmp_path / "f.json", {"A": "kinase", "B": "gpcr"})
    mod = _load_mod()
    outdir = tmp_path / "figs"
    rc = mod.main([
        "--bench-json", str(bench),
        "--family-map", str(fam),
        "--outdir", str(outdir),
    ])
    assert rc == 0
    js = json.loads((outdir / "family_stratification.json").read_text())
    assert "overall" in js
    assert js["overall"]["n_pdbs"] == 2
    assert js["overall"]["n_passed"] == 1
