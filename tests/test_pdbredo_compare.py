"""Tests for stack_protein_preparation._pdbredo_compare."""
from __future__ import annotations

import json
from pathlib import Path

import pytest

from stack_protein_preparation import _pdbredo_compare as prc


# ---------------------------------------------------------------------------
# constants
# ---------------------------------------------------------------------------


def test_lower_and_higher_is_better_are_disjoint():
    assert prc.LOWER_IS_BETTER.isdisjoint(prc.HIGHER_IS_BETTER)


def test_fruton_to_canonical_covers_common_metrics():
    for canonical in ("clashscore", "rama_outlier", "rama_favoured"):
        assert canonical in prc.FRUTON_TO_CANONICAL.values()


# ---------------------------------------------------------------------------
# load_pdbredo_metrics_json
# ---------------------------------------------------------------------------


def test_load_pdbredo_json_happy_path(tmp_path: Path):
    p = tmp_path / "redo.json"
    p.write_text(json.dumps({
        "8abc": {"clashscore": 3.5, "rama_favoured": 97.2, "rama_outlier": 0.1},
        "8XYZ": {"clashscore": 5.1},
    }))
    m = prc.load_pdbredo_metrics_json(p)
    assert set(m.keys()) == {"8ABC", "8XYZ"}
    assert m["8ABC"].clashscore == 3.5
    assert m["8ABC"].rama_favoured == 97.2
    assert m["8XYZ"].clashscore == 5.1


def test_load_pdbredo_json_drops_unknown_fields(tmp_path: Path):
    """Unknown PDB-REDO metric fields must be silently dropped."""
    p = tmp_path / "redo.json"
    p.write_text(json.dumps({
        "8ABC": {"clashscore": 3.5, "future_metric_2027": 42.0},
    }))
    m = prc.load_pdbredo_metrics_json(p)
    assert m["8ABC"].clashscore == 3.5


def test_load_pdbredo_json_rejects_non_dict(tmp_path: Path):
    p = tmp_path / "redo.json"
    p.write_text(json.dumps([{"pdb_id": "8ABC"}]))
    with pytest.raises(ValueError, match="not a dict"):
        prc.load_pdbredo_metrics_json(p)


def test_load_pdbredo_json_ignores_non_dict_values(tmp_path: Path):
    p = tmp_path / "redo.json"
    p.write_text(json.dumps({"8ABC": "not a dict"}))
    m = prc.load_pdbredo_metrics_json(p)
    assert m == {}


def test_load_pdbredo_json_drops_null_values(tmp_path: Path):
    p = tmp_path / "redo.json"
    p.write_text(json.dumps({
        "8ABC": {"clashscore": 3.5, "rama_favoured": None},
    }))
    m = prc.load_pdbredo_metrics_json(p)
    assert m["8ABC"].clashscore == 3.5
    assert m["8ABC"].rama_favoured is None


# ---------------------------------------------------------------------------
# _canonical_fruton_metrics
# ---------------------------------------------------------------------------


def test_canonical_fruton_metrics_translates_field_names():
    r = {"pdb": "8ABC", "clash": 5, "n_rama_outlier_pct": 0.3}
    out = prc._canonical_fruton_metrics(r)
    assert out["clashscore"] == 5.0
    assert out["rama_outlier"] == 0.3


def test_canonical_fruton_metrics_skips_non_numeric():
    r = {"pdb": "8ABC", "clash": "not-a-number"}
    out = prc._canonical_fruton_metrics(r)
    assert "clashscore" not in out


# ---------------------------------------------------------------------------
# compare_bench_to_pdbredo
# ---------------------------------------------------------------------------


def _pdbredo(pdb_id: str, **kw) -> prc.PdbRedoMetrics:
    return prc.PdbRedoMetrics(pdb_id=pdb_id, **kw)


def test_compare_emits_one_delta_per_shared_metric():
    fruton = [{"pdb": "8ABC", "clash": 5, "n_rama_outlier_pct": 0.5}]
    pdbredo = {"8ABC": _pdbredo("8ABC", clashscore=3.5, rama_outlier=0.2)}
    ds = prc.compare_bench_to_pdbredo(fruton, pdbredo)
    metrics = {d.metric for d in ds}
    assert metrics == {"clashscore", "rama_outlier"}


def test_compare_skips_pdbs_missing_in_pdbredo():
    fruton = [{"pdb": "8ABC", "clash": 5}, {"pdb": "8XYZ", "clash": 3}]
    pdbredo = {"8ABC": _pdbredo("8ABC", clashscore=3.5)}
    ds = prc.compare_bench_to_pdbredo(fruton, pdbredo)
    assert {d.pdb_id for d in ds} == {"8ABC"}


def test_compare_upcases_pdb_id():
    fruton = [{"pdb": "8abc", "clash": 5}]
    pdbredo = {"8ABC": _pdbredo("8ABC", clashscore=3.5)}
    ds = prc.compare_bench_to_pdbredo(fruton, pdbredo)
    assert len(ds) == 1


def test_compare_skips_metric_absent_in_pdbredo():
    fruton = [{"pdb": "8ABC", "clash": 5}]
    pdbredo = {"8ABC": _pdbredo("8ABC", rama_favoured=97.0)}  # no clashscore
    ds = prc.compare_bench_to_pdbredo(fruton, pdbredo)
    assert ds == []


# ---------------------------------------------------------------------------
# MetricDelta.fruton_better
# ---------------------------------------------------------------------------


def test_fruton_better_for_lower_is_better_metric():
    d = prc.MetricDelta("8ABC", "clashscore", fruton_value=2.0, pdbredo_value=5.0)
    assert d.delta() == -3.0
    assert d.fruton_better() is True  # lower clashscore = better


def test_fruton_worse_for_lower_is_better_metric():
    d = prc.MetricDelta("8ABC", "clashscore", fruton_value=8.0, pdbredo_value=5.0)
    assert d.fruton_better() is False


def test_fruton_better_for_higher_is_better_metric():
    d = prc.MetricDelta("8ABC", "rama_favoured", fruton_value=98.0, pdbredo_value=95.0)
    assert d.fruton_better() is True  # more favoured = better


def test_fruton_better_returns_none_on_missing_value():
    d = prc.MetricDelta("8ABC", "clashscore", fruton_value=None, pdbredo_value=5.0)
    assert d.fruton_better() is None


def test_fruton_better_returns_none_on_unclassified_metric():
    d = prc.MetricDelta("8ABC", "custom_metric", fruton_value=2.0, pdbredo_value=5.0)
    assert d.fruton_better() is None


# ---------------------------------------------------------------------------
# aggregate_deltas_by_metric
# ---------------------------------------------------------------------------


def test_aggregate_computes_mean_median_and_win_counts():
    deltas = [
        prc.MetricDelta("A", "clashscore", 2.0, 5.0),   # FRUTON better (-3)
        prc.MetricDelta("B", "clashscore", 3.0, 4.0),   # FRUTON better (-1)
        prc.MetricDelta("C", "clashscore", 6.0, 5.0),   # PDB-REDO better (+1)
    ]
    stats = prc.aggregate_deltas_by_metric(deltas)
    s = stats["clashscore"]
    assert s.n_pairs == 3
    assert s.n_fruton_better == 2
    assert s.n_pdbredo_better == 1
    assert s.mean_delta == pytest.approx(-1.0)
    assert s.median_delta == -1.0


def test_aggregate_groups_by_metric():
    deltas = [
        prc.MetricDelta("A", "clashscore", 2.0, 5.0),
        prc.MetricDelta("A", "rama_favoured", 96.0, 95.0),
    ]
    stats = prc.aggregate_deltas_by_metric(deltas)
    assert set(stats.keys()) == {"clashscore", "rama_favoured"}
    assert stats["rama_favoured"].n_fruton_better == 1


def test_aggregate_ignores_deltas_with_missing_value():
    deltas = [
        prc.MetricDelta("A", "clashscore", None, 5.0),
        prc.MetricDelta("B", "clashscore", 2.0, 5.0),
    ]
    stats = prc.aggregate_deltas_by_metric(deltas)
    assert stats["clashscore"].n_pairs == 1


# ---------------------------------------------------------------------------
# format_comparison_table + summarise_delta_direction
# ---------------------------------------------------------------------------


def test_format_comparison_table_lists_all_metrics_sorted():
    deltas = [
        prc.MetricDelta("A", "clashscore", 2.0, 5.0),
        prc.MetricDelta("A", "rama_favoured", 96.0, 95.0),
    ]
    stats = prc.aggregate_deltas_by_metric(deltas)
    md = prc.format_comparison_table(stats)
    # sorted alphabetically -> clashscore before rama_favoured
    idx_clash = md.index("clashscore")
    idx_rama = md.index("rama_favoured")
    assert idx_clash < idx_rama


def test_format_comparison_table_empty_case():
    md = prc.format_comparison_table({})
    assert "no overlapping metrics" in md


def test_summarise_delta_direction_reports_winner():
    deltas = [
        prc.MetricDelta("A", "clashscore", 2.0, 5.0),  # FRUTON
        prc.MetricDelta("B", "clashscore", 3.0, 4.0),  # FRUTON
    ]
    stats = prc.aggregate_deltas_by_metric(deltas)
    lines = prc.summarise_delta_direction(stats)
    assert len(lines) == 1
    assert "FRUTON wins on 2/2" in lines[0]


# ---------------------------------------------------------------------------
# end-to-end via JSON files
# ---------------------------------------------------------------------------


def test_end_to_end_with_json_files(tmp_path: Path):
    fruton_json = tmp_path / "bench.json"
    fruton_json.write_text(json.dumps([
        {"pdb": "8ABC", "clash": 3, "n_rama_favoured_pct": 96.5},
        {"pdb": "8XYZ", "clash": 4},
    ]))
    pdbredo_json = tmp_path / "redo.json"
    pdbredo_json.write_text(json.dumps({
        "8ABC": {"clashscore": 5.0, "rama_favoured": 95.0},
        "8XYZ": {"clashscore": 3.0},
    }))

    fruton = json.loads(fruton_json.read_text())
    pdbredo = prc.load_pdbredo_metrics_json(pdbredo_json)
    deltas = prc.compare_bench_to_pdbredo(fruton, pdbredo)
    stats = prc.aggregate_deltas_by_metric(deltas)

    assert stats["clashscore"].n_pairs == 2
    assert stats["rama_favoured"].n_pairs == 1
    assert stats["clashscore"].n_fruton_better == 1  # 8ABC: 3<5
    assert stats["clashscore"].n_pdbredo_better == 1  # 8XYZ: 4>3
