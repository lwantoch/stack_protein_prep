"""End-to-end integration smoke test for the reviewer-figure scripts.

Composes plot_family_stratification + plot_clashscore_histogram +
plot_ablation_comparison + build_family_assignment against a single
tiny synthetic bench dataset in a single fixture directory, exercising
the full reviewer-facing pipeline.  If any script's CLI or module
contract regresses this test catches it in one shot.

Uses the shipped family_by_pdb_seed.json so the family script picks up
canonical labels.  All matplotlib output uses the Agg backend so no
display is needed.
"""
from __future__ import annotations

import importlib.util
import json
import sys
from pathlib import Path

import pytest


_SCRIPTS_DIR = Path(__file__).resolve().parents[1] / "scripts"


def _load(script_name: str):
    spec = importlib.util.spec_from_file_location(
        f"_{script_name}", _SCRIPTS_DIR / f"{script_name}.py",
    )
    mod = importlib.util.module_from_spec(spec)
    sys.modules[f"_{script_name}"] = mod
    spec.loader.exec_module(mod)
    return mod


def _rec(pdb: str, *, gate_pass: bool, n_clash: int,
         delta_n: int, brk: int = 0) -> dict:
    return {
        "pdb": pdb,
        "gate_pass": gate_pass,
        "clash": n_clash,
        "n_clash_pairs": n_clash,
        "n_vdw_clashes": n_clash // 2,
        "brk": brk,
        "delta_n": delta_n,
        "refine_seconds": 0.0,
    }


@pytest.fixture
def synthetic_bench(tmp_path: Path) -> dict:
    """Emit a tiny synthetic bench dir + baseline + variant JSONs."""
    bench = tmp_path / "bench"
    bench.mkdir()

    # 6 PDBs from 4 families in the shipped seed
    baseline = [
        _rec("1ATP", gate_pass=True, n_clash=2, delta_n=10),      # kinase
        _rec("2SRC", gate_pass=True, n_clash=1, delta_n=5),       # kinase
        _rec("1F88", gate_pass=False, n_clash=45, delta_n=0),     # gpcr
        _rec("1CA2", gate_pass=True, n_clash=3, delta_n=8),       # metalloenzyme
        _rec("1TRN", gate_pass=True, n_clash=1, delta_n=4),       # protease
        _rec("XXXX", gate_pass=False, n_clash=60, delta_n=0),     # unassigned
    ]
    (bench / "baseline.json").write_text(json.dumps(baseline))

    # Variant: disabling clash gate rescues 1F88 and XXXX
    no_clash = list(baseline)
    no_clash[2] = _rec("1F88", gate_pass=True, n_clash=45, delta_n=3)
    no_clash[5] = _rec("XXXX", gate_pass=True, n_clash=60, delta_n=2)
    (bench / "no_clash.json").write_text(json.dumps(no_clash))

    return {
        "bench_dir": bench,
        "baseline_json": bench / "baseline.json",
        "no_clash_json": bench / "no_clash.json",
    }


# ---------------------------------------------------------------------------
# Individual scripts
# ---------------------------------------------------------------------------


def test_build_family_assignment_covers_seed_entries(
    synthetic_bench, tmp_path: Path,
):
    mod = _load("build_family_assignment")
    emitted = tmp_path / "merged_family.json"
    rc = mod.main([
        "--bench-json", str(synthetic_bench["baseline_json"]),
        "--emit", str(emitted),
        "--fill-unassigned",
    ])
    assert rc == 0
    payload = json.loads(emitted.read_text())
    # 5 of 6 covered by the shipped seed; XXXX -> __unassigned__
    assert payload["1ATP"] == "kinase"
    assert payload["1F88"] == "gpcr"
    assert payload["1CA2"] == "metalloenzyme"
    assert payload["1TRN"] == "protease"
    assert payload["XXXX"] == "__unassigned__"


def test_plot_family_stratification_end_to_end(synthetic_bench, tmp_path: Path):
    """Uses the merged family map from build_family_assignment."""
    bfa = _load("build_family_assignment")
    plot = _load("plot_family_stratification")

    merged = tmp_path / "merged_family.json"
    bfa.main([
        "--bench-json", str(synthetic_bench["baseline_json"]),
        "--emit", str(merged),
        "--fill-unassigned",
    ])

    outdir = tmp_path / "family_fig"
    rc = plot.main([
        "--bench-json", str(synthetic_bench["baseline_json"]),
        "--family-map", str(merged),
        "--outdir", str(outdir),
    ])
    assert rc == 0
    for name in ("family_stratification.png", "family_stratification.md",
                 "family_stratification.json"):
        assert (outdir / name).is_file(), f"missing {name}"

    strat = json.loads((outdir / "family_stratification.json").read_text())
    families = {f["family"] for f in strat["per_family"]}
    assert {"kinase", "gpcr", "metalloenzyme", "protease"}.issubset(families)
    # Both kinases pass -> kinase pass_rate == 1.0
    kinase = next(f for f in strat["per_family"] if f["family"] == "kinase")
    assert kinase["n_pdbs"] == 2
    assert kinase["pass_rate"] == 1.0


def test_plot_clashscore_histogram_end_to_end(synthetic_bench, tmp_path: Path):
    plot = _load("plot_clashscore_histogram")
    outdir = tmp_path / "clash_fig"
    rc = plot.main([
        "--bench-json", str(synthetic_bench["baseline_json"]),
        "--outdir", str(outdir),
    ])
    assert rc == 0
    for name in ("clashscore_histogram.png", "clashscore_histogram.csv",
                 "summary.txt"):
        assert (outdir / name).is_file(), f"missing {name}"

    summary = (outdir / "summary.txt").read_text()
    # n=6 records total
    assert "n=6" in summary
    # p90 above 40 (1F88=45 and XXXX=60 dominate)
    csv_text = (outdir / "clashscore_histogram.csv").read_text()
    assert "1ATP" in csv_text
    assert "XXXX" in csv_text


def test_plot_ablation_comparison_end_to_end(synthetic_bench, tmp_path: Path):
    plot = _load("plot_ablation_comparison")
    outdir = tmp_path / "abl_fig"
    rc = plot.main([
        "--baseline", str(synthetic_bench["baseline_json"]),
        "--variant", f"no_clash={synthetic_bench['no_clash_json']}:clash",
        "--outdir", str(outdir),
    ])
    assert rc == 0
    for name in ("ablation_comparison.png", "ablation_comparison.md",
                 "ablation_comparison.json"):
        assert (outdir / name).is_file(), f"missing {name}"

    payload = json.loads((outdir / "ablation_comparison.json").read_text())
    # baseline: 4 pass / 6; variant: 6 pass / 6
    assert payload["baseline"]["n_passed"] == 4
    assert payload["variants"][0]["n_passed"] == 6
    assert payload["variants"][0]["disabled_gate"] == "clash"
    # Deltas should include two 'rescued' entries (1F88 + XXXX)
    deltas = payload["per_pdb_deltas_by_variant"]["no_clash"]
    signs = [d["sign"] for d in deltas]
    assert signs.count("rescued") == 2


# ---------------------------------------------------------------------------
# Composed pipeline: family + clash + ablation from the same bench
# ---------------------------------------------------------------------------


def test_full_reviewer_bundle_composition(synthetic_bench, tmp_path: Path):
    """Assert that four reviewer scripts all run on the same bench
    without conflicts and each produces its expected artefacts."""
    bfa = _load("build_family_assignment")
    fam = _load("plot_family_stratification")
    clash = _load("plot_clashscore_histogram")
    abl = _load("plot_ablation_comparison")

    outdir = tmp_path / "bundle"
    (outdir / "family").mkdir(parents=True)
    (outdir / "clash").mkdir(parents=True)
    (outdir / "ablation").mkdir(parents=True)

    merged = outdir / "merged_family.json"
    assert bfa.main([
        "--bench-json", str(synthetic_bench["baseline_json"]),
        "--emit", str(merged),
        "--fill-unassigned",
    ]) == 0
    assert fam.main([
        "--bench-json", str(synthetic_bench["baseline_json"]),
        "--family-map", str(merged),
        "--outdir", str(outdir / "family"),
    ]) == 0
    assert clash.main([
        "--bench-json", str(synthetic_bench["baseline_json"]),
        "--outdir", str(outdir / "clash"),
    ]) == 0
    assert abl.main([
        "--baseline", str(synthetic_bench["baseline_json"]),
        "--variant", f"no_clash={synthetic_bench['no_clash_json']}:clash",
        "--outdir", str(outdir / "ablation"),
    ]) == 0

    # Each script wrote its expected trio (png/md/json or png/csv/summary).
    assert (outdir / "family" / "family_stratification.png").is_file()
    assert (outdir / "family" / "family_stratification.md").is_file()
    assert (outdir / "family" / "family_stratification.json").is_file()
    assert (outdir / "clash" / "clashscore_histogram.png").is_file()
    assert (outdir / "clash" / "clashscore_histogram.csv").is_file()
    assert (outdir / "clash" / "summary.txt").is_file()
    assert (outdir / "ablation" / "ablation_comparison.png").is_file()
    assert (outdir / "ablation" / "ablation_comparison.md").is_file()
    assert (outdir / "ablation" / "ablation_comparison.json").is_file()
