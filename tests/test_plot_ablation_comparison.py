"""Smoke tests for scripts/plot_ablation_comparison.py."""
from __future__ import annotations

import importlib.util
import json
import sys
from pathlib import Path

import pytest


_SCRIPT = (
    Path(__file__).resolve().parents[1]
    / "scripts" / "plot_ablation_comparison.py"
)


def _load_mod():
    spec = importlib.util.spec_from_file_location("_plot_ablation_cli", _SCRIPT)
    mod = importlib.util.module_from_spec(spec)
    sys.modules["_plot_ablation_cli"] = mod
    spec.loader.exec_module(mod)
    return mod


def _write(path: Path, records: list[dict]) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(records))
    return path


def _rec(pdb: str, gate_pass: bool = True, clash: int = 0, delta_n: int = 0) -> dict:
    return {
        "pdb": pdb, "gate_pass": gate_pass, "clash": clash,
        "brk": 0, "delta_n": delta_n, "refine_seconds": 0.0,
    }


def test_variant_spec_parser():
    mod = _load_mod()
    assert mod._parse_variant("v=/tmp/x.json:plddt") == ("v", "/tmp/x.json", "plddt")
    assert mod._parse_variant("v=/tmp/x.json") == ("v", "/tmp/x.json", None)


def test_variant_spec_parser_rejects_bare_path():
    mod = _load_mod()
    with pytest.raises(Exception):
        mod._parse_variant("/tmp/x.json")


def test_end_to_end_writes_png_md_json(tmp_path: Path):
    baseline = _write(tmp_path / "b.json", [
        _rec("A", gate_pass=True, clash=3, delta_n=10),
        _rec("B", gate_pass=False, clash=50, delta_n=0),
    ])
    var = _write(tmp_path / "v.json", [
        _rec("A", gate_pass=True, clash=3, delta_n=10),
        _rec("B", gate_pass=True, clash=50, delta_n=4),
    ])
    mod = _load_mod()
    outdir = tmp_path / "figs"
    rc = mod.main([
        "--baseline", str(baseline),
        "--variant", f"no_clash={var}:clash",
        "--outdir", str(outdir),
    ])
    assert rc == 0
    assert (outdir / "ablation_comparison.png").is_file()
    assert (outdir / "ablation_comparison.md").is_file()
    assert (outdir / "ablation_comparison.json").is_file()

    md = (outdir / "ablation_comparison.md").read_text()
    assert "baseline" in md
    assert "no_clash" in md

    js = json.loads((outdir / "ablation_comparison.json").read_text())
    assert "baseline" in js
    assert js["variants"][0]["variant"] == "no_clash"
    assert js["variants"][0]["disabled_gate"] == "clash"


def test_end_to_end_variant_without_gate_label(tmp_path: Path):
    baseline = _write(tmp_path / "b.json", [_rec("A", gate_pass=True)])
    var = _write(tmp_path / "v.json", [_rec("A", gate_pass=True)])
    mod = _load_mod()
    outdir = tmp_path / "figs"
    rc = mod.main([
        "--baseline", str(baseline),
        "--variant", f"v_only={var}",
        "--outdir", str(outdir),
    ])
    assert rc == 0
    js = json.loads((outdir / "ablation_comparison.json").read_text())
    assert js["variants"][0]["disabled_gate"] is None


def test_end_to_end_missing_baseline_returns_nonzero(tmp_path: Path):
    mod = _load_mod()
    var = _write(tmp_path / "v.json", [_rec("A")])
    outdir = tmp_path / "figs"
    rc = mod.main([
        "--baseline", str(tmp_path / "nope.json"),
        "--variant", f"v={var}",
        "--outdir", str(outdir),
    ])
    assert rc == 2


def test_end_to_end_unknown_gate_returns_nonzero(tmp_path: Path):
    mod = _load_mod()
    baseline = _write(tmp_path / "b.json", [_rec("A")])
    var = _write(tmp_path / "v.json", [_rec("A")])
    outdir = tmp_path / "figs"
    rc = mod.main([
        "--baseline", str(baseline),
        "--variant", f"weird={var}:not_a_gate",
        "--outdir", str(outdir),
    ])
    assert rc == 2
