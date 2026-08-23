"""Tests for stack_protein_preparation._tleap_generator.

Verifies the emitted tleap.in has the right sections for common
combinations (protein-only, protein+ions, +cofactor, +nonstandard,
+lanthanide), that ion source lines dedup correctly, that missing
model paths still produce a script + a warning rather than raising,
and that param-directory scans skip orphan files.
"""
from __future__ import annotations

from pathlib import Path

import pytest

from stack_protein_preparation import _tleap_generator as tg


def _make_pdb(path: Path) -> Path:
    path.write_text(
        "ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00  0.00           N  \n"
        "END\n",
        encoding="utf-8",
    )
    return path


def _make_param_pair(dir_: Path, resname: str) -> tuple[Path, Path]:
    dir_.mkdir(parents=True, exist_ok=True)
    mol2 = dir_ / f"{resname}.mol2"
    frcmod = dir_ / f"{resname}.frcmod"
    mol2.write_text(f"@<TRIPOS>MOLECULE\n{resname}\n", encoding="utf-8")
    frcmod.write_text(f"remark {resname}\nBOND\n", encoding="utf-8")
    return mol2, frcmod


def test_protein_only_emits_minimal_script(tmp_path: Path):
    pdb = _make_pdb(tmp_path / "prot.pdb")
    r = tg.write_tleap_script(pdb, tmp_path / "out")
    assert r.tleap_in_path.is_file()
    content = r.tleap_in_path.read_text()
    assert "source leaprc.protein.ff14SB" in content
    assert "source leaprc.water.tip3p" in content
    assert f"loadpdb {pdb.resolve()}" in content
    assert "solvatebox mol TIP3PBOX 10.0" in content
    assert "addions mol Na+ 0" in content
    assert "saveamberparm mol system.prmtop system.rst7" in content
    assert "quit" in content
    assert r.n_cofactors_loaded == 0
    assert r.n_nonstandard_loaded == 0


def test_missing_model_still_emits_script_with_warning(tmp_path: Path):
    r = tg.write_tleap_script(tmp_path / "does_not_exist.pdb", tmp_path / "out")
    assert r.tleap_in_path.is_file()
    assert any("does not exist" in w for w in r.warnings)


def test_ion_resnames_source_dedup(tmp_path: Path):
    pdb = _make_pdb(tmp_path / "prot.pdb")
    r = tg.write_tleap_script(
        pdb, tmp_path / "out",
        ion_resnames=["ZN", "MG", "CA", "FE2"],
    )
    content = r.tleap_in_path.read_text()
    assert content.count("source leaprc.water.tip3p_ion12-6-4") == 1
    assert r.n_ion_leaprc_lines == 1


def test_lanthanide_adds_hfe_leaprc(tmp_path: Path):
    pdb = _make_pdb(tmp_path / "prot.pdb")
    r = tg.write_tleap_script(pdb, tmp_path / "out", ion_resnames=["ZN", "GD"])
    content = r.tleap_in_path.read_text()
    assert "source leaprc.water.tip3p_ion12-6-4" in content
    assert "HFE" in content
    assert r.n_ion_leaprc_lines == 2


def test_cofactor_params_loaded(tmp_path: Path):
    pdb = _make_pdb(tmp_path / "prot.pdb")
    cof = tmp_path / "cof"
    _make_param_pair(cof, "NAD")
    _make_param_pair(cof, "XYZ")
    r = tg.write_tleap_script(pdb, tmp_path / "out", cofactor_param_dir=cof)
    content = r.tleap_in_path.read_text()
    assert "NAD = loadmol2" in content
    assert "XYZ = loadmol2" in content
    assert content.count("loadamberparams") == 2
    assert "source leaprc.gaff2" in content
    assert r.n_cofactors_loaded == 2


def test_orphan_param_files_are_skipped(tmp_path: Path):
    pdb = _make_pdb(tmp_path / "prot.pdb")
    cof = tmp_path / "cof"
    cof.mkdir()
    # Orphan mol2 with no frcmod -> skipped
    (cof / "ORPHAN.mol2").write_text("@<TRIPOS>MOLECULE\nORPHAN\n")
    _make_param_pair(cof, "COMPLETE")
    r = tg.write_tleap_script(pdb, tmp_path / "out", cofactor_param_dir=cof)
    content = r.tleap_in_path.read_text()
    assert "COMPLETE = loadmol2" in content
    assert "ORPHAN = loadmol2" not in content
    assert r.n_cofactors_loaded == 1


def test_nonstandard_and_cofactor_both_loaded(tmp_path: Path):
    pdb = _make_pdb(tmp_path / "prot.pdb")
    cof = tmp_path / "cof"; _make_param_pair(cof, "NAD")
    nst = tmp_path / "nst"; _make_param_pair(nst, "SEP")
    r = tg.write_tleap_script(
        pdb, tmp_path / "out",
        cofactor_param_dir=cof, nonstandard_param_dir=nst,
    )
    content = r.tleap_in_path.read_text()
    assert "NAD = loadmol2" in content
    assert "SEP = loadmol2" in content
    assert r.n_cofactors_loaded == 1
    assert r.n_nonstandard_loaded == 1


def test_neutralise_can_be_disabled(tmp_path: Path):
    pdb = _make_pdb(tmp_path / "prot.pdb")
    r = tg.write_tleap_script(pdb, tmp_path / "out", neutralise=False)
    content = r.tleap_in_path.read_text()
    assert "addions" not in content
    assert "solvatebox mol TIP3PBOX 10.0" in content


def test_custom_topology_basename(tmp_path: Path):
    pdb = _make_pdb(tmp_path / "prot.pdb")
    r = tg.write_tleap_script(pdb, tmp_path / "out", topology_basename="fruton_run1")
    content = r.tleap_in_path.read_text()
    assert "saveamberparm mol fruton_run1.prmtop fruton_run1.rst7" in content
    assert "savepdb mol fruton_run1.pdb" in content


def test_unknown_ion_resname_is_dropped_silently(tmp_path: Path):
    """UNKNOWN metals are absent from _ion_params -> no source line added."""
    pdb = _make_pdb(tmp_path / "prot.pdb")
    r = tg.write_tleap_script(pdb, tmp_path / "out", ion_resnames=["XXX"])
    content = r.tleap_in_path.read_text()
    assert "tip3p_ion12-6-4" not in content
    assert r.n_ion_leaprc_lines == 0


def test_output_directory_is_created_if_missing(tmp_path: Path):
    pdb = _make_pdb(tmp_path / "prot.pdb")
    nested = tmp_path / "deep" / "nested" / "out"
    r = tg.write_tleap_script(pdb, nested)
    assert nested.is_dir()
    assert r.tleap_in_path.is_file()


def test_tleap_in_uses_absolute_paths(tmp_path: Path):
    pdb = _make_pdb(tmp_path / "prot.pdb")
    r = tg.write_tleap_script(pdb, tmp_path / "out")
    content = r.tleap_in_path.read_text()
    # loadpdb line must contain the absolute path
    assert str(pdb.resolve()) in content


def test_result_to_dict_contains_metrics(tmp_path: Path):
    pdb = _make_pdb(tmp_path / "prot.pdb")
    r = tg.write_tleap_script(pdb, tmp_path / "out", ion_resnames=["ZN"])
    d = r.to_dict()
    assert d["tleap_in_path"].endswith("tleap.in")
    assert d["n_cofactors_loaded"] == 0
    assert d["n_ion_leaprc_lines"] == 1


def test_header_comments_reflect_inputs(tmp_path: Path):
    """Reviewer-friendly header lists what was fed in."""
    pdb = _make_pdb(tmp_path / "prot.pdb")
    r = tg.write_tleap_script(
        pdb, tmp_path / "out",
        ion_resnames=["ZN", "MG"],
        solvate_margin_angstrom=12.0,
    )
    content = r.tleap_in_path.read_text()
    assert "solvate margin:   12.0 A" in content
    assert "['MG', 'ZN']" in content  # sorted set
