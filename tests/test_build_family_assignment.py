"""Tests for scripts/build_family_assignment.py."""
from __future__ import annotations

import importlib.util
import io
import json
import sys
from contextlib import redirect_stdout
from pathlib import Path

import pytest


_SCRIPT = (
    Path(__file__).resolve().parents[1]
    / "scripts" / "build_family_assignment.py"
)


def _load_mod():
    spec = importlib.util.spec_from_file_location("_bfa", _SCRIPT)
    mod = importlib.util.module_from_spec(spec)
    sys.modules["_bfa"] = mod
    spec.loader.exec_module(mod)
    return mod


def _write(path: Path, obj) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(obj))
    return path


def _rec(pdb: str) -> dict:
    return {"pdb": pdb}


# ---------------------------------------------------------------------------
# _load_bench_pdb_ids
# ---------------------------------------------------------------------------


def test_load_bench_uppercases_and_dedupes(tmp_path: Path):
    mod = _load_mod()
    p = _write(tmp_path / "b.json", [_rec("8abc"), _rec("8ABC"), _rec("8xyz")])
    r = mod._load_bench_pdb_ids(p)
    assert r == ["8ABC", "8XYZ"]


def test_load_bench_skips_empty_pdb(tmp_path: Path):
    mod = _load_mod()
    p = _write(tmp_path / "b.json", [{"pdb": ""}, _rec("A"), {"other": 1}])
    r = mod._load_bench_pdb_ids(p)
    assert r == ["A"]


def test_load_bench_rejects_non_list(tmp_path: Path):
    mod = _load_mod()
    p = _write(tmp_path / "b.json", {"pdb": "A"})
    with pytest.raises(ValueError, match="not a list"):
        mod._load_bench_pdb_ids(p)


# ---------------------------------------------------------------------------
# _merge_maps
# ---------------------------------------------------------------------------


def test_merge_maps_later_overrides():
    mod = _load_mod()
    seed = {"A": "kinase", "B": "gpcr"}
    extra = {"B": "protease", "C": "hydrolase"}
    r = mod._merge_maps(seed, [extra])
    assert r == {"A": "kinase", "B": "protease", "C": "hydrolase"}


def test_merge_maps_uppercases_extras():
    mod = _load_mod()
    seed = {"A": "kinase"}
    extra = {"b": "gpcr"}
    r = mod._merge_maps(seed, [extra])
    assert "B" in r
    assert "b" not in r


# ---------------------------------------------------------------------------
# build_coverage_report
# ---------------------------------------------------------------------------


def test_coverage_report_counts():
    mod = _load_mod()
    r = mod.build_coverage_report(
        ["A", "B", "C", "D"],
        {"A": "kinase", "B": "kinase", "C": "gpcr"},
    )
    assert r["n_bench_pdbs"] == 4
    assert r["n_covered"] == 3
    assert r["n_unassigned"] == 1
    assert r["unassigned"] == ["D"]
    assert r["by_family"]["kinase"] == ["A", "B"]
    assert r["by_family"]["gpcr"] == ["C"]
    assert r["coverage_pct"] == pytest.approx(75.0)


def test_coverage_report_all_unassigned():
    mod = _load_mod()
    r = mod.build_coverage_report(["X", "Y"], {})
    assert r["n_covered"] == 0
    assert r["n_unassigned"] == 2
    assert r["by_family"] == {}


def test_coverage_report_empty_bench():
    mod = _load_mod()
    r = mod.build_coverage_report([], {"A": "kinase"})
    assert r["n_bench_pdbs"] == 0
    assert r["coverage_pct"] == 0.0


# ---------------------------------------------------------------------------
# format_report
# ---------------------------------------------------------------------------


def test_format_report_includes_headline_line():
    mod = _load_mod()
    r = mod.build_coverage_report(["A", "B"], {"A": "kinase"})
    t = mod.format_report(r)
    assert "bench_pdbs=2" in t
    assert "covered=1" in t
    assert "unassigned=1" in t
    assert "kinase" in t


def test_format_report_omits_family_block_when_empty():
    mod = _load_mod()
    r = mod.build_coverage_report(["A"], {})
    t = mod.format_report(r)
    assert "Per-family coverage" not in t
    assert "Unassigned" in t


# ---------------------------------------------------------------------------
# CLI end-to-end
# ---------------------------------------------------------------------------


def test_end_to_end_prints_summary(tmp_path: Path):
    mod = _load_mod()
    bench = _write(tmp_path / "b.json", [_rec("1ATP"), _rec("XXXX")])
    buf = io.StringIO()
    with redirect_stdout(buf):
        rc = mod.main(["--bench-json", str(bench)])
    assert rc == 0
    out = buf.getvalue()
    assert "bench_pdbs=2" in out
    # 1ATP is in the shipped seed; XXXX is not.
    assert "kinase" in out
    assert "XXXX" in out  # in the unassigned list


def test_end_to_end_emit_writes_merged_map(tmp_path: Path):
    mod = _load_mod()
    bench = _write(tmp_path / "b.json", [_rec("1ATP"), _rec("XXXX")])
    emitted = tmp_path / "merged.json"
    buf = io.StringIO()
    with redirect_stdout(buf):
        rc = mod.main([
            "--bench-json", str(bench),
            "--emit", str(emitted),
        ])
    assert rc == 0
    payload = json.loads(emitted.read_text())
    assert payload["1ATP"] == "kinase"
    # Without --fill-unassigned, XXXX is omitted
    assert "XXXX" not in payload


def test_end_to_end_fill_unassigned_includes_missing_pdbs(tmp_path: Path):
    mod = _load_mod()
    bench = _write(tmp_path / "b.json", [_rec("1ATP"), _rec("XXXX")])
    emitted = tmp_path / "merged.json"
    buf = io.StringIO()
    with redirect_stdout(buf):
        rc = mod.main([
            "--bench-json", str(bench),
            "--emit", str(emitted),
            "--fill-unassigned",
        ])
    assert rc == 0
    payload = json.loads(emitted.read_text())
    assert payload["1ATP"] == "kinase"
    assert payload["XXXX"] == "__unassigned__"


def test_end_to_end_extra_map_overrides_seed(tmp_path: Path):
    mod = _load_mod()
    bench = _write(tmp_path / "b.json", [_rec("1ATP")])
    extra = _write(tmp_path / "extra.json", {"1ATP": "phosphatase"})
    emitted = tmp_path / "merged.json"
    buf = io.StringIO()
    with redirect_stdout(buf):
        rc = mod.main([
            "--bench-json", str(bench),
            "--extra-map", str(extra),
            "--emit", str(emitted),
        ])
    assert rc == 0
    payload = json.loads(emitted.read_text())
    assert payload["1ATP"] == "phosphatase"


def test_end_to_end_non_canonical_label_prints_warning(tmp_path: Path, capsys):
    mod = _load_mod()
    bench = _write(tmp_path / "b.json", [_rec("ZZZZ")])
    extra = _write(tmp_path / "extra.json", {"ZZZZ": "not_a_family"})
    rc = mod.main([
        "--bench-json", str(bench),
        "--extra-map", str(extra),
    ])
    captured = capsys.readouterr()
    assert "non-canonical labels" in captured.err
    assert rc == 0


def test_missing_bench_json_returns_nonzero(tmp_path: Path):
    mod = _load_mod()
    rc = mod.main(["--bench-json", str(tmp_path / "nope.json")])
    assert rc == 2


def test_missing_seed_returns_nonzero(tmp_path: Path):
    mod = _load_mod()
    bench = _write(tmp_path / "b.json", [_rec("A")])
    rc = mod.main([
        "--bench-json", str(bench),
        "--seed", str(tmp_path / "nope.json"),
    ])
    assert rc == 2
