"""Tests for scripts/aggregate_reviewer_bench_output.py."""
from __future__ import annotations

import importlib.util
import json
import sys
from pathlib import Path

import pytest


_SCRIPT = (
    Path(__file__).resolve().parents[1]
    / "scripts" / "aggregate_reviewer_bench_output.py"
)


def _load():
    spec = importlib.util.spec_from_file_location("_arb", _SCRIPT)
    mod = importlib.util.module_from_spec(spec)
    sys.modules["_arb"] = mod
    spec.loader.exec_module(mod)
    return mod


def _write_provenance(path: Path, headline_bullets: list[str]) -> None:
    """Emit a minimal PROVENANCE.md with a Headline numbers section."""
    path.parent.mkdir(parents=True, exist_ok=True)
    body = ["# Figure title", "", "Some prose.", "", "## Headline numbers", ""]
    body += [f"- {b}" for b in headline_bullets]
    body += ["", "## Some other section", "", "Ignored bullet:", "- ignored"]
    path.write_text("\n".join(body) + "\n")


# ---------------------------------------------------------------------------
# _extract_headline_bullets
# ---------------------------------------------------------------------------


def test_extract_finds_bullets_in_headline_section(tmp_path: Path):
    mod = _load()
    p = tmp_path / "PROVENANCE.md"
    _write_provenance(p, ["n = 199", "mean = 2.60", "p90 = 6.2"])
    bullets = mod._extract_headline_bullets(p)
    assert bullets == ["n = 199", "mean = 2.60", "p90 = 6.2"]


def test_extract_stops_at_next_heading(tmp_path: Path):
    mod = _load()
    p = tmp_path / "PROVENANCE.md"
    _write_provenance(p, ["n = 199"])
    bullets = mod._extract_headline_bullets(p)
    assert bullets == ["n = 199"]
    assert "ignored" not in " ".join(bullets)


def test_extract_returns_empty_when_no_headline_section(tmp_path: Path):
    mod = _load()
    p = tmp_path / "PROVENANCE.md"
    p.write_text("# Title\n\nJust prose.\n")
    assert mod._extract_headline_bullets(p) == []


def test_extract_headline_case_insensitive(tmp_path: Path):
    mod = _load()
    p = tmp_path / "PROVENANCE.md"
    p.write_text("## HEADLINE NUMBERS\n\n- a = 1\n")
    assert mod._extract_headline_bullets(p) == ["a = 1"]


def test_extract_headline_with_parenthetical(tmp_path: Path):
    mod = _load()
    p = tmp_path / "PROVENANCE.md"
    p.write_text("## Headline numbers (post-FRUTON)\n\n- x = 1\n")
    assert mod._extract_headline_bullets(p) == ["x = 1"]


# ---------------------------------------------------------------------------
# _bullet_to_kv
# ---------------------------------------------------------------------------


def test_bullet_to_kv_equals():
    mod = _load()
    assert mod._bullet_to_kv("n = 199") == ("n", "199")


def test_bullet_to_kv_colon():
    mod = _load()
    assert mod._bullet_to_kv("mean: 2.60") == ("mean", "2.60")


def test_bullet_to_kv_strips_bold():
    mod = _load()
    r = mod._bullet_to_kv("**n = 199**")
    assert r == ("n", "199")


def test_bullet_to_kv_returns_none_for_freeform():
    mod = _load()
    assert mod._bullet_to_kv("this is a sentence") is None


def test_bullet_to_kv_handles_units():
    mod = _load()
    r = mod._bullet_to_kv("mean clashscore = 2.60 clash pairs / crystal")
    assert r == ("mean clashscore", "2.60 clash pairs / crystal")


# ---------------------------------------------------------------------------
# collect_figures + summarise_figure
# ---------------------------------------------------------------------------


def test_summarise_figure_returns_expected_shape(tmp_path: Path):
    mod = _load()
    p = tmp_path / "figA" / "PROVENANCE.md"
    _write_provenance(p, ["n = 46", "gate PASS = 42/46 (91.3%)"])
    fig = mod.summarise_figure("figA", p)
    assert fig["figure_dir"] == "figA"
    assert fig["headline_kv"] == {
        "n": "46",
        "gate PASS": "42/46 (91.3%)",
    }
    assert fig["freeform_bullets"] == []


def test_summarise_figure_captures_freeform(tmp_path: Path):
    mod = _load()
    p = tmp_path / "figB" / "PROVENANCE.md"
    _write_provenance(p, [
        "n = 46",
        "this is a note without kv",
    ])
    fig = mod.summarise_figure("figB", p)
    assert fig["headline_kv"] == {"n": "46"}
    assert fig["freeform_bullets"] == ["this is a note without kv"]


def test_collect_figures_walks_multiple_subdirs(tmp_path: Path):
    mod = _load()
    for fig_dir in ("fig1", "fig2", "fig3"):
        _write_provenance(tmp_path / fig_dir / "PROVENANCE.md", ["n = 1"])
    (tmp_path / "no_prov_dir").mkdir()  # skipped — no PROVENANCE.md
    figures = mod.collect_figures(tmp_path)
    assert {f["figure_dir"] for f in figures} == {"fig1", "fig2", "fig3"}


def test_collect_figures_sorted_alphabetically(tmp_path: Path):
    mod = _load()
    for name in ("z_fig", "a_fig", "m_fig"):
        _write_provenance(tmp_path / name / "PROVENANCE.md", ["n = 1"])
    figures = mod.collect_figures(tmp_path)
    assert [f["figure_dir"] for f in figures] == ["a_fig", "m_fig", "z_fig"]


# ---------------------------------------------------------------------------
# format_markdown_summary
# ---------------------------------------------------------------------------


def test_markdown_summary_lists_every_figure(tmp_path: Path):
    mod = _load()
    for name in ("figA", "figB"):
        _write_provenance(tmp_path / name / "PROVENANCE.md",
                          [f"n = {name}"])
    figures = mod.collect_figures(tmp_path)
    md = mod.format_markdown_summary(figures)
    assert "## figA" in md
    assert "## figB" in md
    assert "n | figA" in md
    assert "n | figB" in md


def test_markdown_summary_escapes_pipe_in_values():
    mod = _load()
    figures = [{
        "figure_dir": "f",
        "provenance_path": "x/PROVENANCE.md",
        "headline_kv": {"k": "a | b"},
        "freeform_bullets": [],
    }]
    md = mod.format_markdown_summary(figures)
    assert "a \\| b" in md


def test_markdown_summary_empty_input():
    mod = _load()
    md = mod.format_markdown_summary([])
    assert "**0 figure(s)**" in md


# ---------------------------------------------------------------------------
# CLI end-to-end
# ---------------------------------------------------------------------------


def test_end_to_end_writes_md_and_json(tmp_path: Path):
    mod = _load()
    _write_provenance(
        tmp_path / "bo" / "figA" / "PROVENANCE.md",
        ["n = 199", "mean = 2.60"],
    )
    outdir = tmp_path / "out"
    rc = mod.main([
        "--bench-output-dir", str(tmp_path / "bo"),
        "--outdir", str(outdir),
    ])
    assert rc == 0
    assert (outdir / "reviewer_summary.md").is_file()
    assert (outdir / "reviewer_summary.json").is_file()
    payload = json.loads((outdir / "reviewer_summary.json").read_text())
    assert len(payload) == 1
    assert payload[0]["headline_kv"]["n"] == "199"


def test_missing_bench_output_dir_returns_nonzero(tmp_path: Path):
    mod = _load()
    rc = mod.main(["--bench-output-dir", str(tmp_path / "nope")])
    assert rc == 2


def test_empty_bench_output_dir_still_writes(tmp_path: Path):
    mod = _load()
    (tmp_path / "bo").mkdir()
    outdir = tmp_path / "out"
    rc = mod.main([
        "--bench-output-dir", str(tmp_path / "bo"),
        "--outdir", str(outdir),
    ])
    assert rc == 0
    payload = json.loads((outdir / "reviewer_summary.json").read_text())
    assert payload == []
