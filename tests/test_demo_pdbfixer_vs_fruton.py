"""Tests for scripts/demo_pdbfixer_vs_fruton_protonation.py."""
from __future__ import annotations

import importlib.util
import io
import sys
from contextlib import redirect_stdout
from pathlib import Path


_SCRIPT = (
    Path(__file__).resolve().parents[1]
    / "scripts" / "demo_pdbfixer_vs_fruton_protonation.py"
)


def _load_mod():
    spec = importlib.util.spec_from_file_location("_demo_ppp", _SCRIPT)
    mod = importlib.util.module_from_spec(spec)
    sys.modules["_demo_ppp"] = mod
    spec.loader.exec_module(mod)
    return mod


# ---------------------------------------------------------------------------
# naive_ph_rule
# ---------------------------------------------------------------------------


def test_naive_rule_his_ph7_is_hie():
    mod = _load_mod()
    assert mod.naive_ph_rule("HIS", pH=7.0) == "HIE"


def test_naive_rule_his_low_ph_is_hip():
    mod = _load_mod()
    assert mod.naive_ph_rule("HIS", pH=5.5) == "HIP"


def test_naive_rule_asp_glu_stay_deprotonated():
    mod = _load_mod()
    assert mod.naive_ph_rule("ASP", pH=7.0) == "ASP"
    assert mod.naive_ph_rule("GLU", pH=7.0) == "GLU"


def test_naive_rule_cys_stays_protonated():
    mod = _load_mod()
    assert mod.naive_ph_rule("CYS", pH=7.0) == "CYS"


# ---------------------------------------------------------------------------
# _scenarios_as_paper_evidence_md
# ---------------------------------------------------------------------------


def test_scenarios_to_evidence_emits_headings_and_quotes():
    mod = _load_mod()
    md = mod._scenarios_as_paper_evidence_md([
        ("demo", "HIS57", "H57 is the catalytic base", "HIP"),
    ])
    assert "### HIS57" in md
    assert "> H57 is the catalytic base" in md


def test_scenarios_to_evidence_skips_unparseable_labels():
    mod = _load_mod()
    md = mod._scenarios_as_paper_evidence_md([
        ("bad", "not_a_label", "some text", "HIP"),
    ])
    assert "### " not in md


# ---------------------------------------------------------------------------
# build_verdicts
# ---------------------------------------------------------------------------


def test_build_verdicts_catches_catalytic_base_as_hip(tmp_path: Path):
    mod = _load_mod()
    ev_path = tmp_path / "e.md"
    ev_path.write_text("### HIS57\n> H57 is the catalytic base in the triad.\n")
    verdicts = mod.build_verdicts(
        [("t", "HIS57", "H57 is the catalytic base in the triad.", "HIP")],
        ev_path,
    )
    assert len(verdicts) == 1
    v = verdicts[0]
    assert v.naive_default_ph7 == "HIE"
    assert v.fruton_override == "HIP"
    assert v.literature_recommended == "HIP"
    assert v.fruton_correct is True
    assert v.naive_wrong is True


def test_build_verdicts_catches_thiolate_nucleophile_as_cym(tmp_path: Path):
    mod = _load_mod()
    ev_path = tmp_path / "e.md"
    ev_path.write_text(
        "### CYS25\n> C25 is the catalytic nucleophile (thiolate form).\n"
    )
    verdicts = mod.build_verdicts(
        [("t", "CYS25", "C25 is the catalytic nucleophile (thiolate form).", "CYM")],
        ev_path,
    )
    assert verdicts[0].fruton_override == "CYM"
    assert verdicts[0].fruton_correct is True


def test_build_verdicts_needs_review_returns_none(tmp_path: Path):
    """Vague evidence with no protonation keyword → FRUTON returns None."""
    mod = _load_mod()
    ev_path = tmp_path / "e.md"
    ev_path.write_text(
        "### ASP162\n> D162 sits at the entrance to the channel.\n"
    )
    verdicts = mod.build_verdicts(
        [("t", "ASP162", "D162 sits at the entrance to the channel.", "ASP")],
        ev_path,
    )
    assert verdicts[0].fruton_override is None
    # naive default matches literature ASP, so this is tie(naive)
    assert verdicts[0].naive_default_ph7 == "ASP"


def test_load_evidence_path_falls_back_to_tmp(tmp_path: Path):
    mod = _load_mod()
    path = mod._load_evidence_path(None, [
        ("t", "HIS57", "H57 is the catalytic base.", "HIP"),
    ])
    assert path.is_file()
    assert "### HIS57" in path.read_text()


def test_residue_number_extractor():
    mod = _load_mod()
    assert mod._residue_number("HIS57") == 57
    assert mod._residue_number("CYS25") == 25
    assert mod._residue_number("not_a_label") is None


def test_resname_from_label():
    mod = _load_mod()
    assert mod._resname_from_label("HIS57") == "HIS"
    assert mod._resname_from_label("CYS25") == "CYS"
    assert mod._resname_from_label("garbage") == ""


# ---------------------------------------------------------------------------
# CLI table + main
# ---------------------------------------------------------------------------


def test_cli_table_contains_headers_and_row():
    mod = _load_mod()
    v = mod.Verdict(
        scenario="s", residue="HIS57",
        naive_default_ph7="HIE",
        fruton_override="HIP",
        literature_recommended="HIP",
    )
    t = mod._cli_table([v])
    assert "Scenario" in t
    assert "HIS57" in t
    assert "FRUTON" in t or "tie" in t


def test_main_runs_end_to_end():
    mod = _load_mod()
    buf = io.StringIO()
    with redirect_stdout(buf):
        rc = mod.main([])
    assert rc == 0
    output = buf.getvalue()
    assert "summary" in output
    assert "FRUTON" in output


def test_main_with_evidence_path(tmp_path: Path):
    mod = _load_mod()
    ev = tmp_path / "e.md"
    ev.write_text("### HIS57\n> H57 is the catalytic base.\n")
    buf = io.StringIO()
    with redirect_stdout(buf):
        rc = mod.main(["--evidence", str(ev)])
    assert rc == 0
