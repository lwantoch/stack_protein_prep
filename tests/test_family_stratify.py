"""Tests for stack_protein_preparation._family_stratify."""
from __future__ import annotations

from pathlib import Path

import pytest

from stack_protein_preparation import _family_stratify as fs


def _rec(pdb: str, gate_pass=True, clash=0, brk=0, delta_n=0) -> dict:
    return {"pdb": pdb, "gate_pass": gate_pass, "clash": clash,
            "brk": brk, "delta_n": delta_n}


def _write_mapping_csv(path: Path, rows: list[tuple[str, str]],
                       header=("pdb_id", "family")) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    lines = [",".join(header)] + [f"{p},{f}" for p, f in rows]
    path.write_text("\n".join(lines) + "\n")
    return path


# ---------------------------------------------------------------------------
# load_family_mapping
# ---------------------------------------------------------------------------


def test_load_mapping_happy_path(tmp_path: Path):
    p = _write_mapping_csv(tmp_path / "m.csv", [
        ("8ABC", "kinase"),
        ("8XYZ", "gpcr"),
    ])
    m = fs.load_family_mapping(p)
    assert m == {"8ABC": "kinase", "8XYZ": "gpcr"}


def test_load_mapping_upcases_pdb_ids(tmp_path: Path):
    p = _write_mapping_csv(tmp_path / "m.csv", [("8abc", "kinase")])
    m = fs.load_family_mapping(p)
    assert "8ABC" in m


def test_load_mapping_accepts_alt_column_names(tmp_path: Path):
    p = _write_mapping_csv(tmp_path / "m.csv",
                           [("8ABC", "kinase")],
                           header=("pdb", "class"))
    m = fs.load_family_mapping(p)
    assert m == {"8ABC": "kinase"}


def test_load_mapping_rejects_missing_columns(tmp_path: Path):
    p = tmp_path / "m.csv"
    p.write_text("randomcol,othercol\n8ABC,kinase\n")
    with pytest.raises(ValueError, match="expected columns"):
        fs.load_family_mapping(p)


def test_load_mapping_skips_rows_with_missing_fields(tmp_path: Path):
    p = tmp_path / "m.csv"
    p.write_text("pdb_id,family\n8ABC,kinase\n,gpcr\n8XYZ,\n8WWW,phospholipase\n")
    m = fs.load_family_mapping(p)
    assert m == {"8ABC": "kinase", "8WWW": "phospholipase"}


def test_load_mapping_empty_file_returns_empty(tmp_path: Path):
    p = tmp_path / "m.csv"
    p.write_text("")
    m = fs.load_family_mapping(p)
    assert m == {}


# ---------------------------------------------------------------------------
# stratify_bench_by_family
# ---------------------------------------------------------------------------


def test_stratify_groups_by_family():
    bench = [
        _rec("A1", gate_pass=True, clash=2, delta_n=8),
        _rec("A2", gate_pass=True, clash=3, delta_n=6),
        _rec("B1", gate_pass=False, clash=50, delta_n=0),
    ]
    mapping = {"A1": "kinase", "A2": "kinase", "B1": "gpcr"}
    aggs = fs.stratify_bench_by_family(bench, mapping)
    by_family = {a.family: a for a in aggs}
    assert by_family["kinase"].n_pdbs == 2
    assert by_family["kinase"].pass_rate == 1.0
    assert by_family["gpcr"].n_pdbs == 1
    assert by_family["gpcr"].pass_rate == 0.0


def test_stratify_computes_correct_aggregates():
    bench = [
        _rec("A1", gate_pass=True, clash=2, brk=0, delta_n=10),
        _rec("A2", gate_pass=True, clash=4, brk=0, delta_n=6),
        _rec("A3", gate_pass=False, clash=30, brk=2, delta_n=0),
    ]
    mapping = {"A1": "kinase", "A2": "kinase", "A3": "kinase"}
    aggs = fs.stratify_bench_by_family(bench, mapping)
    assert len(aggs) == 1
    a = aggs[0]
    assert a.n_pdbs == 3
    assert a.n_passed == 2
    assert a.pass_rate == pytest.approx(2/3)
    assert a.mean_clash == pytest.approx((2 + 4 + 30) / 3)
    assert a.total_delta_n == 16


def test_stratify_upcases_bench_pdb_ids():
    bench = [_rec("a1", gate_pass=True)]
    mapping = {"A1": "kinase"}
    aggs = fs.stratify_bench_by_family(bench, mapping)
    assert aggs[0].family == "kinase"


def test_stratify_include_unassigned_true():
    bench = [_rec("A1", gate_pass=True), _rec("MISS", gate_pass=False)]
    mapping = {"A1": "kinase"}
    aggs = fs.stratify_bench_by_family(bench, mapping, include_unassigned=True)
    families = {a.family for a in aggs}
    assert fs.UNKNOWN_FAMILY in families


def test_stratify_include_unassigned_false_drops():
    bench = [_rec("A1", gate_pass=True), _rec("MISS", gate_pass=False)]
    mapping = {"A1": "kinase"}
    aggs = fs.stratify_bench_by_family(bench, mapping, include_unassigned=False)
    families = {a.family for a in aggs}
    assert fs.UNKNOWN_FAMILY not in families
    assert families == {"kinase"}


def test_stratify_sorts_unassigned_last():
    bench = [_rec("ZZZ", gate_pass=True), _rec("A1", gate_pass=True)]
    mapping = {"A1": "kinase"}
    aggs = fs.stratify_bench_by_family(bench, mapping)
    assert aggs[-1].family == fs.UNKNOWN_FAMILY


def test_stratify_member_pdbs_populated():
    bench = [_rec("A1", gate_pass=True), _rec("A2", gate_pass=True)]
    mapping = {"A1": "kinase", "A2": "kinase"}
    aggs = fs.stratify_bench_by_family(bench, mapping)
    assert aggs[0].member_pdbs == ["A1", "A2"]


def test_stratify_empty_bench_returns_empty_list():
    aggs = fs.stratify_bench_by_family([], {})
    assert aggs == []


def test_stratify_skips_records_with_no_pdb():
    bench = [_rec("A1", gate_pass=True), {"gate_pass": True}]  # missing pdb
    mapping = {"A1": "kinase"}
    aggs = fs.stratify_bench_by_family(bench, mapping)
    assert sum(a.n_pdbs for a in aggs) == 1


# ---------------------------------------------------------------------------
# reporting helpers
# ---------------------------------------------------------------------------


def test_format_family_table_has_expected_headers():
    aggs = fs.stratify_bench_by_family(
        [_rec("A1", gate_pass=True, clash=2, delta_n=5)],
        {"A1": "kinase"},
    )
    md = fs.format_family_table(aggs)
    assert "Family" in md
    assert "kinase" in md
    assert "pass %" in md


def test_format_family_table_empty():
    md = fs.format_family_table([])
    assert "no family aggregates" in md


def test_hardest_families_orders_by_pass_rate():
    bench = [
        _rec("A1", gate_pass=True), _rec("A2", gate_pass=True), _rec("A3", gate_pass=True),
        _rec("B1", gate_pass=False), _rec("B2", gate_pass=False), _rec("B3", gate_pass=True),
        _rec("C1", gate_pass=True), _rec("C2", gate_pass=False), _rec("C3", gate_pass=True),
    ]
    mapping = {
        "A1": "kinase", "A2": "kinase", "A3": "kinase",
        "B1": "gpcr", "B2": "gpcr", "B3": "gpcr",
        "C1": "phospholipase", "C2": "phospholipase", "C3": "phospholipase",
    }
    aggs = fs.stratify_bench_by_family(bench, mapping)
    hardest = fs.hardest_families(aggs, min_family_size=3, top_n=2)
    assert hardest[0].family == "gpcr"        # 1/3 pass
    assert hardest[1].family == "phospholipase"  # 2/3 pass


def test_hardest_families_respects_min_size_gate():
    """A family with fewer than min_family_size members is excluded."""
    bench = [
        _rec("A1", gate_pass=True), _rec("A2", gate_pass=True), _rec("A3", gate_pass=True),
        _rec("B1", gate_pass=False),  # single-member family
    ]
    mapping = {"A1": "k", "A2": "k", "A3": "k", "B1": "solo"}
    aggs = fs.stratify_bench_by_family(bench, mapping)
    hardest = fs.hardest_families(aggs, min_family_size=3, top_n=5)
    assert {a.family for a in hardest} == {"k"}


def test_hardest_families_excludes_unassigned():
    bench = [
        _rec("A1", gate_pass=True), _rec("A2", gate_pass=True), _rec("A3", gate_pass=True),
        _rec("Z1", gate_pass=False), _rec("Z2", gate_pass=False), _rec("Z3", gate_pass=False),
    ]
    mapping = {"A1": "kinase", "A2": "kinase", "A3": "kinase"}  # Z* unassigned
    aggs = fs.stratify_bench_by_family(bench, mapping)
    hardest = fs.hardest_families(aggs, min_family_size=3, top_n=5)
    assert {a.family for a in hardest} == {"kinase"}


def test_aggregate_to_dict_round_trip():
    aggs = fs.stratify_bench_by_family(
        [_rec("A1", gate_pass=True, clash=2, delta_n=5)],
        {"A1": "kinase"},
    )
    d = aggs[0].to_dict()
    assert d["family"] == "kinase"
    assert d["n_pdbs"] == 1
    assert d["mean_clash"] == 2.0
    assert "member_pdbs" in d
