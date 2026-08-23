"""Tests for stack_protein_preparation._ablation."""
from __future__ import annotations

import json
from pathlib import Path

import pytest

from stack_protein_preparation import _ablation as ab


def _write_json(path: Path, records: list[dict]) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(records))
    return path


def _fake_record(
    pdb: str, *, gate_pass: bool = True, clash: int = 0,
    brk: int = 0, delta_n: int = 0, refine_seconds: float = 0.0,
) -> dict:
    return {
        "pdb": pdb, "gate_pass": gate_pass, "clash": clash,
        "brk": brk, "delta_n": delta_n, "refine_seconds": refine_seconds,
    }


# ---------------------------------------------------------------------------
# GATE_NAMES catalogue
# ---------------------------------------------------------------------------


def test_gate_names_covers_expected_gates():
    """Sanity: nine canonical toggle-able gates."""
    expected = {
        "plddt", "idr", "clash", "chirality_d", "chirality_improper",
        "omega", "rollback", "loop_refine", "reviewer_gate",
    }
    assert set(ab.GATE_NAMES.keys()) == expected


def test_gate_names_have_nonempty_descriptions():
    for name, desc in ab.GATE_NAMES.items():
        assert isinstance(desc, str) and len(desc) > 20


# ---------------------------------------------------------------------------
# load_bench_json
# ---------------------------------------------------------------------------


def test_load_bench_json_happy_path(tmp_path: Path):
    p = _write_json(tmp_path / "bench.json", [_fake_record("8ABC")])
    r = ab.load_bench_json(p)
    assert len(r) == 1
    assert r[0]["pdb"] == "8ABC"


def test_load_bench_json_rejects_non_list(tmp_path: Path):
    p = tmp_path / "bench.json"
    p.write_text('{"pdb": "8ABC"}')  # dict, not list
    with pytest.raises(ValueError, match="not a list"):
        ab.load_bench_json(p)


# ---------------------------------------------------------------------------
# compute_variant_aggregate
# ---------------------------------------------------------------------------


def test_aggregate_counts_pass_and_totals():
    results = [
        _fake_record("A", gate_pass=True, clash=2, brk=0, delta_n=10),
        _fake_record("B", gate_pass=True, clash=3, brk=1, delta_n=5),
        _fake_record("C", gate_pass=False, clash=50, brk=8, delta_n=0),
    ]
    agg = ab.compute_variant_aggregate("baseline", results)
    assert agg.n_pdbs == 3
    assert agg.n_passed == 2
    assert agg.pass_rate == pytest.approx(2 / 3)
    assert agg.total_delta_n == 15
    assert agg.mean_delta_n == pytest.approx(5.0)
    # (2 + 3 + 50) / 3
    assert agg.mean_clash == pytest.approx(55 / 3)


def test_aggregate_empty_bench():
    agg = ab.compute_variant_aggregate("empty", [])
    assert agg.n_pdbs == 0
    assert agg.pass_rate == 0.0
    assert agg.mean_clash == 0.0
    assert agg.total_delta_n == 0


def test_aggregate_missing_fields_default_to_zero():
    """Records without 'clash', 'brk', 'delta_n' are still counted."""
    results = [{"pdb": "A", "gate_pass": True}]
    agg = ab.compute_variant_aggregate("sparse", results)
    assert agg.mean_clash == 0.0
    assert agg.mean_brk == 0.0
    assert agg.total_delta_n == 0
    assert agg.n_passed == 1


def test_aggregate_records_disabled_gate_label():
    agg = ab.compute_variant_aggregate("no_plddt", [], disabled_gate="plddt")
    assert agg.disabled_gate == "plddt"


def test_aggregate_to_dict_round_trip():
    agg = ab.compute_variant_aggregate(
        "v1", [_fake_record("A", clash=3, brk=1, delta_n=5)],
        disabled_gate="clash",
    )
    d = agg.to_dict()
    assert d["variant"] == "v1"
    assert d["disabled_gate"] == "clash"
    assert d["n_pdbs"] == 1
    assert d["mean_clash"] == 3.0


# ---------------------------------------------------------------------------
# diff_per_pdb
# ---------------------------------------------------------------------------


def test_diff_per_pdb_identifies_rescued_and_lost():
    baseline = [
        _fake_record("A", gate_pass=False, clash=50, delta_n=0),  # fails baseline
        _fake_record("B", gate_pass=True, clash=2, delta_n=8),   # passes baseline
    ]
    variant = [
        _fake_record("A", gate_pass=True, clash=45, delta_n=6),  # now passes
        _fake_record("B", gate_pass=False, clash=60, delta_n=8), # now fails
    ]
    deltas = ab.diff_per_pdb(baseline, variant)
    by_pdb = {d.pdb: d for d in deltas}
    assert by_pdb["A"].sign() == "rescued"
    assert by_pdb["B"].sign() == "lost"


def test_diff_per_pdb_skips_pdbs_only_in_variant():
    baseline = [_fake_record("A")]
    variant = [_fake_record("A"), _fake_record("Z_new")]
    deltas = ab.diff_per_pdb(baseline, variant)
    assert {d.pdb for d in deltas} == {"A"}


def test_diff_per_pdb_marks_unchanged():
    baseline = [_fake_record("A", gate_pass=True, delta_n=5, clash=3)]
    variant = [_fake_record("A", gate_pass=True, delta_n=5, clash=3)]
    deltas = ab.diff_per_pdb(baseline, variant)
    assert deltas[0].sign() == "unchanged"


# ---------------------------------------------------------------------------
# build_comparison — end-to-end
# ---------------------------------------------------------------------------


def test_build_comparison_end_to_end(tmp_path: Path):
    base = _write_json(tmp_path / "base.json", [
        _fake_record("A", gate_pass=True, clash=3, delta_n=10),
        _fake_record("B", gate_pass=False, clash=50, delta_n=0),
    ])
    no_clash_gate = _write_json(tmp_path / "no_clash.json", [
        _fake_record("A", gate_pass=True, clash=3, delta_n=10),
        _fake_record("B", gate_pass=True, clash=50, delta_n=4),  # rescued by disabling clash gate
    ])
    comp = ab.build_comparison(
        base,
        {"no_clash": no_clash_gate},
        disabled_gate_by_variant={"no_clash": "clash"},
    )
    assert comp.baseline.n_passed == 1
    assert comp.variants[0].n_passed == 2
    assert comp.variants[0].disabled_gate == "clash"
    deltas = comp.per_pdb_deltas["no_clash"]
    assert any(d.sign() == "rescued" and d.pdb == "B" for d in deltas)


def test_build_comparison_rejects_unknown_gate(tmp_path: Path):
    base = _write_json(tmp_path / "base.json", [_fake_record("A")])
    var = _write_json(tmp_path / "v.json", [_fake_record("A")])
    with pytest.raises(ValueError, match="not in GATE_NAMES"):
        ab.build_comparison(
            base, {"weird": var},
            disabled_gate_by_variant={"weird": "not_a_real_gate"},
        )


# ---------------------------------------------------------------------------
# reporting helpers
# ---------------------------------------------------------------------------


def test_format_comparison_table_contains_baseline_and_variant(tmp_path: Path):
    base = _write_json(tmp_path / "base.json", [
        _fake_record("A", gate_pass=True, clash=3, delta_n=10),
        _fake_record("B", gate_pass=False, clash=50, delta_n=0),
    ])
    var = _write_json(tmp_path / "no_omega.json", [
        _fake_record("A", gate_pass=True, clash=3, delta_n=10),
        _fake_record("B", gate_pass=True, clash=10, delta_n=6),
    ])
    comp = ab.build_comparison(
        base, {"no_omega": var},
        disabled_gate_by_variant={"no_omega": "omega"},
    )
    md = ab.format_comparison_table(comp)
    assert "| baseline |" in md
    assert "| no_omega |" in md
    assert "omega" in md
    assert "50.0%" in md  # baseline pass 1/2 = 50%
    assert "100.0%" in md  # variant 2/2 = 100%


def test_rescued_and_lost_counts(tmp_path: Path):
    base = _write_json(tmp_path / "base.json", [
        _fake_record("A", gate_pass=False),
        _fake_record("B", gate_pass=True),
        _fake_record("C", gate_pass=True),
    ])
    var = _write_json(tmp_path / "v.json", [
        _fake_record("A", gate_pass=True),   # rescued
        _fake_record("B", gate_pass=False),  # lost
        _fake_record("C", gate_pass=True),   # unchanged
    ])
    comp = ab.build_comparison(base, {"v1": var})
    counts = ab.rescued_and_lost_counts(comp)
    assert counts["v1"] == (1, 1)


def test_comparison_to_dict_is_json_serialisable(tmp_path: Path):
    base = _write_json(tmp_path / "base.json", [_fake_record("A")])
    var = _write_json(tmp_path / "v.json", [_fake_record("A")])
    comp = ab.build_comparison(base, {"v": var})
    d = comp.to_dict()
    s = json.dumps(d)  # must not raise
    assert "baseline" in s
    assert "per_pdb_deltas_by_variant" in s
