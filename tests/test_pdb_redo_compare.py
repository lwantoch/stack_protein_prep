"""Tests for stack_protein_preparation._pdb_redo_compare."""
from __future__ import annotations

import json
from pathlib import Path

import pytest

from stack_protein_preparation import _pdb_redo_compare as prc


def _write(path: Path, obj: dict) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(obj))
    return path


# ---------------------------------------------------------------------------
# METRIC_NAMES catalogue
# ---------------------------------------------------------------------------


def test_metric_names_cover_expected_set():
    expected = {
        "n_rama_outlier", "n_clash_pairs", "n_omega_non_planar",
        "n_omega_cis_nonpro", "clashscore", "rama_favoured_pct",
        "rotamer_outlier_pct", "bond_rms_z", "angle_rms_z",
    }
    assert set(prc.METRIC_NAMES.keys()) == expected


def test_every_metric_has_direction_and_label():
    for k, spec in prc.METRIC_NAMES.items():
        assert "lower_is_better" in spec
        assert isinstance(spec["reviewer_label"], str)
        assert isinstance(spec["units"], str)


# ---------------------------------------------------------------------------
# load_metrics_dir
# ---------------------------------------------------------------------------


def test_load_metrics_dir_reads_json_files(tmp_path: Path):
    _write(tmp_path / "8ABC.json", {"n_rama_outlier": 3})
    _write(tmp_path / "8xyz.json", {"n_rama_outlier": 7})
    out = prc.load_metrics_dir(tmp_path)
    assert set(out.keys()) == {"8ABC", "8XYZ"}
    assert out["8ABC"]["n_rama_outlier"] == 3


def test_load_metrics_dir_raises_on_missing_dir(tmp_path: Path):
    with pytest.raises(FileNotFoundError):
        prc.load_metrics_dir(tmp_path / "nope")


def test_load_metrics_dir_raises_on_list_payload(tmp_path: Path):
    _write(tmp_path / "8ABC.json", {})
    (tmp_path / "8XYZ.json").write_text("[1,2,3]")
    with pytest.raises(ValueError, match="not a JSON object"):
        prc.load_metrics_dir(tmp_path)


# ---------------------------------------------------------------------------
# PerPdbMetricDelta
# ---------------------------------------------------------------------------


def test_delta_negative_for_lower_is_better_favours_fruton():
    e = prc.PerPdbMetricDelta(
        pdb="8ABC", metric="n_rama_outlier",
        fruton_value=3, pdb_redo_value=7,
    )
    assert e.delta() == -4.0
    assert e.better_side(lower_is_better=True) == "fruton"


def test_delta_positive_for_lower_is_better_favours_pdb_redo():
    e = prc.PerPdbMetricDelta(
        pdb="8ABC", metric="clashscore",
        fruton_value=10, pdb_redo_value=5,
    )
    assert e.delta() == 5.0
    assert e.better_side(lower_is_better=True) == "pdb_redo"


def test_delta_positive_for_higher_is_better_favours_fruton():
    e = prc.PerPdbMetricDelta(
        pdb="8ABC", metric="rama_favoured_pct",
        fruton_value=98.0, pdb_redo_value=95.0,
    )
    assert e.delta() == 3.0
    assert e.better_side(lower_is_better=False) == "fruton"


def test_delta_returns_none_when_either_side_missing():
    e = prc.PerPdbMetricDelta(
        pdb="8ABC", metric="clashscore",
        fruton_value=None, pdb_redo_value=5,
    )
    assert e.delta() is None
    assert e.better_side(lower_is_better=True) == "unavailable"


def test_delta_tie_when_equal():
    e = prc.PerPdbMetricDelta(
        pdb="8ABC", metric="n_rama_outlier",
        fruton_value=3, pdb_redo_value=3,
    )
    assert e.better_side(lower_is_better=True) == "tie"


# ---------------------------------------------------------------------------
# build_comparison
# ---------------------------------------------------------------------------


def test_build_comparison_end_to_end_matches_shared_pdbs_only(tmp_path: Path):
    frdir = tmp_path / "fruton"
    prdir = tmp_path / "pdb_redo"
    _write(frdir / "8ABC.json", {"n_rama_outlier": 3, "clashscore": 5.0})
    _write(frdir / "8XYZ.json", {"n_rama_outlier": 8, "clashscore": 10.0})
    _write(frdir / "9ONLY_FRUTON.json", {"n_rama_outlier": 1})
    _write(prdir / "8ABC.json", {"n_rama_outlier": 7, "clashscore": 12.0})
    _write(prdir / "8XYZ.json", {"n_rama_outlier": 4, "clashscore": 3.0})
    _write(prdir / "9ONLY_REDO.json", {"n_rama_outlier": 9})

    comp = prc.build_comparison(frdir, prdir, metrics=["n_rama_outlier", "clashscore"])

    # only 8ABC + 8XYZ appear (intersection)
    assert set(comp.per_pdb_deltas) == {"8ABC", "8XYZ"}

    rama = next(m for m in comp.per_metric if m.metric == "n_rama_outlier")
    # 8ABC: fruton 3 < pdb_redo 7 → fruton better
    # 8XYZ: fruton 8 > pdb_redo 4 → pdb_redo better
    assert rama.n_fruton_better == 1
    assert rama.n_pdb_redo_better == 1
    assert rama.n_tie == 0
    assert rama.mean_delta == pytest.approx(0.0)  # (-4 + 4)/2

    clash = next(m for m in comp.per_metric if m.metric == "clashscore")
    # 8ABC: -7 (fruton better), 8XYZ: +7 (pdb_redo better)
    assert clash.mean_delta == pytest.approx(0.0)


def test_build_comparison_rejects_unknown_metric(tmp_path: Path):
    frdir = tmp_path / "f"
    prdir = tmp_path / "p"
    _write(frdir / "8ABC.json", {"n_rama_outlier": 3})
    _write(prdir / "8ABC.json", {"n_rama_outlier": 7})
    with pytest.raises(ValueError, match="unknown metric"):
        prc.build_comparison(frdir, prdir, metrics=["not_a_metric"])


def test_build_comparison_handles_all_default_metrics(tmp_path: Path):
    frdir = tmp_path / "f"
    prdir = tmp_path / "p"
    all_metrics = {k: 1.0 for k in prc.METRIC_NAMES}
    _write(frdir / "8ABC.json", all_metrics)
    _write(prdir / "8ABC.json", all_metrics)
    comp = prc.build_comparison(frdir, prdir)
    assert len(comp.per_metric) == len(prc.METRIC_NAMES)


def test_build_comparison_counts_unavailable_when_missing_metric(tmp_path: Path):
    frdir = tmp_path / "f"
    prdir = tmp_path / "p"
    _write(frdir / "8ABC.json", {})  # no metrics at all
    _write(prdir / "8ABC.json", {"clashscore": 5.0})
    comp = prc.build_comparison(frdir, prdir, metrics=["clashscore"])
    m = comp.per_metric[0]
    assert m.n_unavailable == 1
    assert m.n_matched == 1
    assert m.mean_delta is None


# ---------------------------------------------------------------------------
# format_comparison_table + to_dict
# ---------------------------------------------------------------------------


def test_format_comparison_table_has_headers_and_rows(tmp_path: Path):
    frdir = tmp_path / "f"
    prdir = tmp_path / "p"
    _write(frdir / "8ABC.json", {"n_rama_outlier": 3, "rama_favoured_pct": 98.0})
    _write(prdir / "8ABC.json", {"n_rama_outlier": 7, "rama_favoured_pct": 95.0})
    comp = prc.build_comparison(
        frdir, prdir,
        metrics=["n_rama_outlier", "rama_favoured_pct"],
    )
    md = prc.format_comparison_table(comp)
    assert "Ramachandran outliers" in md
    assert "Ramachandran favoured %" in md
    assert "↓ better" in md
    assert "↑ better" in md
    assert "-4.000" in md  # fruton 3 vs redo 7
    assert "+3.000" in md  # 98 vs 95


def test_to_dict_is_json_serialisable(tmp_path: Path):
    frdir = tmp_path / "f"
    prdir = tmp_path / "p"
    _write(frdir / "8ABC.json", {"n_rama_outlier": 3})
    _write(prdir / "8ABC.json", {"n_rama_outlier": 7})
    comp = prc.build_comparison(frdir, prdir, metrics=["n_rama_outlier"])
    d = comp.to_dict()
    s = json.dumps(d)
    assert "per_metric" in s
    assert "per_pdb_deltas" in s
