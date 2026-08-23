"""Tests for stack_protein_preparation._md_input_generator."""
from __future__ import annotations

from pathlib import Path

import pytest

from stack_protein_preparation import _md_input_generator as mig


def test_writes_all_five_stages_plus_readme(tmp_path: Path):
    r = mig.write_md_input_decks(tmp_path / "md")
    for stage in ("min1", "min2", "heat", "eq", "prod", "readme"):
        assert stage in r.files
        assert r.files[stage].is_file()


def test_min1_has_restraint_wt_10(tmp_path: Path):
    r = mig.write_md_input_decks(tmp_path / "md")
    content = r.files["min1"].read_text()
    assert "restraint_wt = 10.0" in content
    assert "restraintmask = '!:WAT & !@H='" in content
    assert "imin = 1" in content


def test_min2_unrestrained(tmp_path: Path):
    r = mig.write_md_input_decks(tmp_path / "md")
    content = r.files["min2"].read_text()
    assert "ntr = 0" in content
    assert "imin = 1" in content


def test_heat_uses_langevin_and_ca_restraint(tmp_path: Path):
    r = mig.write_md_input_decks(tmp_path / "md")
    content = r.files["heat"].read_text()
    assert "ntt = 3" in content
    assert f"gamma_ln = {mig.DEFAULT_LANGEVIN_GAMMA_PS:.1f}" in content
    assert "restraintmask = '@CA'" in content
    assert "tempi = 100.0" in content
    assert "TEMP0" in content


def test_heat_step_count_matches_duration(tmp_path: Path):
    """20 ps at 2 fs step = 10000 steps."""
    r = mig.write_md_input_decks(tmp_path / "md", heat_duration_ps=20.0, dt_fs=2.0)
    content = r.files["heat"].read_text()
    assert "nstlim = 10000" in content


def test_eq_uses_berendsen_barostat(tmp_path: Path):
    """Berendsen (barostat=1) during eq for fast relaxation."""
    r = mig.write_md_input_decks(tmp_path / "md")
    content = r.files["eq"].read_text()
    assert "barostat = 1" in content
    assert "ntp = 1" in content
    assert "irest = 1" in content
    assert "ntx = 5" in content


def test_prod_uses_mc_barostat_and_iwrap(tmp_path: Path):
    """Monte-Carlo barostat for prod per Roe & Brooks JCTC 2020."""
    r = mig.write_md_input_decks(tmp_path / "md")
    content = r.files["prod"].read_text()
    assert "barostat = 2" in content
    assert "mcbarint = 100" in content
    assert "iwrap = 1" in content


def test_prod_length_10ns_default(tmp_path: Path):
    """10 ns * 500000 steps/ns @ 2 fs = 5000000 steps."""
    r = mig.write_md_input_decks(tmp_path / "md")
    content = r.files["prod"].read_text()
    assert "nstlim = 5000000" in content


def test_custom_prod_duration(tmp_path: Path):
    """1 ns at 2 fs step = 500000 steps."""
    r = mig.write_md_input_decks(tmp_path / "md", prod_duration_ns=1.0)
    content = r.files["prod"].read_text()
    assert "nstlim = 500000" in content


def test_custom_temperature(tmp_path: Path):
    r = mig.write_md_input_decks(tmp_path / "md", target_temperature_K=310.0)
    heat = r.files["heat"].read_text()
    eq = r.files["eq"].read_text()
    prod = r.files["prod"].read_text()
    assert "temp0 = 310.0" in heat
    assert "temp0 = 310.0" in eq
    assert "temp0 = 310.0" in prod


def test_aggressive_timestep_produces_warning(tmp_path: Path):
    """dt > 2.5 fs without HMR is unstable -- must warn."""
    r = mig.write_md_input_decks(tmp_path / "md", dt_fs=4.0)
    assert any("2.5" in w or "HMR" in w for w in r.warnings)


def test_long_prod_produces_warning(tmp_path: Path):
    r = mig.write_md_input_decks(tmp_path / "md", prod_duration_ns=2000.0)
    assert any("us" in w or "disk" in w for w in r.warnings)


def test_shake_settings_consistent_across_dynamic_stages(tmp_path: Path):
    """ntc=2 + ntf=2 in heat, eq, prod (allows 2 fs step)."""
    r = mig.write_md_input_decks(tmp_path / "md")
    for stage in ("heat", "eq", "prod"):
        content = r.files[stage].read_text()
        assert "ntc = 2" in content, f"{stage}.in missing ntc=2"
        assert "ntf = 2" in content, f"{stage}.in missing ntf=2"


def test_readme_documents_pmemd_command_chain(tmp_path: Path):
    r = mig.write_md_input_decks(tmp_path / "md")
    readme = r.files["readme"].read_text()
    assert "pmemd.cuda" in readme
    assert "-i min1.in" in readme
    assert "-i prod.in" in readme
    assert "system.prmtop" in readme


def test_write_readme_can_be_disabled(tmp_path: Path):
    r = mig.write_md_input_decks(tmp_path / "md", write_readme=False)
    assert "readme" not in r.files


def test_custom_prefix_flows_into_titles(tmp_path: Path):
    r = mig.write_md_input_decks(tmp_path / "md", prefix="fruton_kin1")
    content = r.files["min1"].read_text()
    assert "fruton_kin1: min1" in content


def test_result_to_dict_includes_all_files(tmp_path: Path):
    r = mig.write_md_input_decks(tmp_path / "md")
    d = r.to_dict()
    assert d["output_dir"].endswith("md")
    assert set(d["files"].keys()) == {"min1", "min2", "heat", "eq", "prod", "readme"}


def test_output_dir_is_created_if_missing(tmp_path: Path):
    nested = tmp_path / "deep" / "nested" / "md"
    r = mig.write_md_input_decks(nested)
    assert nested.is_dir()
    assert r.files["min1"].is_file()


def test_all_files_are_nonempty(tmp_path: Path):
    r = mig.write_md_input_decks(tmp_path / "md")
    for name, path in r.files.items():
        assert path.stat().st_size > 0, f"{name} is empty"
