"""Tests for stack_protein_preparation._family_stratification."""
from __future__ import annotations

import json
from pathlib import Path

import pytest

from stack_protein_preparation import _family_stratification as fs


def _rec(pdb: str, *, gate_pass: bool = True, clash: int = 0,
         brk: int = 0, delta_n: int = 0) -> dict:
    return {
        "pdb": pdb, "gate_pass": gate_pass, "clash": clash,
        "brk": brk, "delta_n": delta_n,
    }


# ---------------------------------------------------------------------------
# FAMILY_LABELS
# ---------------------------------------------------------------------------


def test_family_labels_contains_expected_canonical_names():
    canonical = {
        "kinase", "gpcr", "protease", "metalloenzyme", "phosphatase",
        "nuclear_receptor", "hydrolase", "transferase", "oxidoreductase",
        "cofactor_dependent", "multi_domain", "other", "__unassigned__",
    }
    assert canonical.issubset(set(fs.FAMILY_LABELS))


# ---------------------------------------------------------------------------
# load_family_mapping
# ---------------------------------------------------------------------------


def test_load_family_mapping_uppercases_pdb_ids(tmp_path: Path):
    p = tmp_path / "map.json"
    p.write_text(json.dumps({"8abc": "kinase", "8XYZ": "gpcr"}))
    m = fs.load_family_mapping(p)
    assert m == {"8ABC": "kinase", "8XYZ": "gpcr"}


def test_load_family_mapping_skips_non_string_values(tmp_path: Path):
    p = tmp_path / "map.json"
    p.write_text(json.dumps({"8ABC": "kinase", "8XYZ": None, "8YYY": 42}))
    m = fs.load_family_mapping(p)
    assert m == {"8ABC": "kinase"}


def test_load_family_mapping_rejects_top_level_list(tmp_path: Path):
    p = tmp_path / "map.json"
    p.write_text(json.dumps(["8ABC", "kinase"]))
    with pytest.raises(ValueError, match="not a JSON object"):
        fs.load_family_mapping(p)


# ---------------------------------------------------------------------------
# stratify_bench_by_family
# ---------------------------------------------------------------------------


def test_stratify_groups_matching_pdbs():
    bench = [
        _rec("A", gate_pass=True, delta_n=10),
        _rec("B", gate_pass=True, delta_n=5),
        _rec("C", gate_pass=False, delta_n=0),
    ]
    mapping = {"A": "kinase", "B": "kinase", "C": "gpcr"}
    strat = fs.stratify_bench_by_family(bench, mapping)
    by_family = {a.family: a for a in strat.per_family}
    assert by_family["kinase"].n_pdbs == 2
    assert by_family["kinase"].n_passed == 2
    assert by_family["kinase"].total_delta_n == 15
    assert by_family["gpcr"].n_pdbs == 1
    assert by_family["gpcr"].n_passed == 0


def test_stratify_unassigned_pdbs_land_in_bucket():
    bench = [
        _rec("A", gate_pass=True),
        _rec("B", gate_pass=False),
    ]
    mapping = {"A": "kinase"}  # B unmapped
    strat = fs.stratify_bench_by_family(bench, mapping)
    labels = {a.family for a in strat.per_family}
    assert "kinase" in labels
    assert "__unassigned__" in labels


def test_stratify_overall_aggregates_all_records():
    bench = [
        _rec("A", gate_pass=True, delta_n=10),
        _rec("B", gate_pass=False, delta_n=0),
    ]
    strat = fs.stratify_bench_by_family(bench, {})
    assert strat.overall.n_pdbs == 2
    assert strat.overall.n_passed == 1
    assert strat.overall.total_delta_n == 10


def test_stratify_orders_canonical_first_unassigned_last():
    """kinase (canonical) < unknown label < __unassigned__."""
    bench = [
        _rec("A"),
        _rec("B"),
        _rec("C"),
    ]
    mapping = {"A": "unknown_odd_family", "B": "kinase"}  # C unmapped
    strat = fs.stratify_bench_by_family(bench, mapping)
    fams = [a.family for a in strat.per_family]
    assert fams.index("kinase") < fams.index("unknown_odd_family")
    assert fams.index("unknown_odd_family") < fams.index("__unassigned__")


def test_stratify_drops_empty_families():
    """An empty family (no matching pdbs) must not appear."""
    bench = [_rec("A")]
    mapping = {"A": "kinase"}
    strat = fs.stratify_bench_by_family(bench, mapping)
    fams = {a.family for a in strat.per_family}
    assert "gpcr" not in fams


def test_stratify_skips_records_without_pdb_field():
    bench = [_rec("A", gate_pass=True), {"gate_pass": True}]
    strat = fs.stratify_bench_by_family(bench, {"A": "kinase"})
    kinase = next(a for a in strat.per_family if a.family == "kinase")
    assert kinase.n_pdbs == 1


# ---------------------------------------------------------------------------
# format_stratification_table + summary
# ---------------------------------------------------------------------------


def test_format_stratification_table_contains_families_and_overall():
    bench = [
        _rec("A", gate_pass=True, delta_n=10),
        _rec("B", gate_pass=False, delta_n=0),
    ]
    mapping = {"A": "kinase", "B": "gpcr"}
    strat = fs.stratify_bench_by_family(bench, mapping)
    md = fs.format_stratification_table(strat)
    assert "| kinase " in md
    assert "| gpcr " in md
    assert "| overall " in md
    assert "100.0%" in md  # kinase 1/1
    assert "0.0%" in md    # gpcr 0/1
    assert "50.0%" in md   # overall 1/2


def test_imbalance_summary_names_best_and_worst():
    bench = [
        _rec("A", gate_pass=True),
        _rec("B", gate_pass=True),
        _rec("C", gate_pass=False),
        _rec("D", gate_pass=False),
    ]
    mapping = {"A": "kinase", "B": "kinase", "C": "gpcr", "D": "gpcr"}
    strat = fs.stratify_bench_by_family(bench, mapping)
    s = fs.imbalance_summary(strat)
    assert "kinase" in s and "100.0%" in s
    assert "gpcr" in s and "0.0%" in s


def test_imbalance_summary_when_empty():
    strat = fs.stratify_bench_by_family([], {})
    assert "no populated families" in fs.imbalance_summary(strat)


def test_to_dict_round_trip_json_serialisable():
    bench = [_rec("A", gate_pass=True, delta_n=10)]
    strat = fs.stratify_bench_by_family(bench, {"A": "kinase"})
    d = strat.to_dict()
    s = json.dumps(d)
    assert "per_family" in s
    assert "overall" in s
