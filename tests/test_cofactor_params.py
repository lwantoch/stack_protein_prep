"""Tests for stack_protein_preparation._cofactor_params.

Exercises the pure-Python helpers (resname extraction, net-charge lookup,
known-library skip logic) and fail-open behaviour when antechamber /
parmchk2 are not on PATH.  Never launches a real antechamber process --
the subprocess call is monkey-patched.
"""
from __future__ import annotations

from pathlib import Path

import pytest

from stack_protein_preparation import _cofactor_params


def _write_pdb(path: Path, lines: list[str]) -> Path:
    path.write_text("\n".join(lines) + "\nEND\n", encoding="utf-8")
    return path


def test_read_pdb_lines_missing_file_returns_empty(tmp_path: Path):
    assert _cofactor_params._read_pdb_lines(tmp_path / "nope.pdb") == []


def test_unique_cofactor_resnames_extracts_only_hetatm_and_atom():
    lines = [
        "HETATM    1  N1  NAD A 501       0.0  0.0  0.0",
        "HETATM    2  N2  NAD A 501       1.0  0.0  0.0",
        "HETATM    3  O   HOH A 502       2.0  0.0  0.0",  # HOH must be skipped
        "REMARK   1                                       ",  # non-atom line
        "HETATM    4  FE  HEM A 503       3.0  0.0  0.0",
        "ATOM      5  CA  ALA A   1       0.0  0.0  0.0",  # protein AA passes through
    ]
    per_res = _cofactor_params._unique_cofactor_resnames(lines)
    assert list(per_res.keys()) == ["NAD", "HEM", "ALA"]
    assert len(per_res["NAD"]) == 2
    assert len(per_res["HEM"]) == 1


def test_guess_net_charge_for_common_cofactors():
    assert _cofactor_params._guess_net_charge_for_cofactor("NAD") == -1
    assert _cofactor_params._guess_net_charge_for_cofactor("ATP") == -4
    assert _cofactor_params._guess_net_charge_for_cofactor("SAM") == 1
    assert _cofactor_params._guess_net_charge_for_cofactor("SAH") == 0
    assert _cofactor_params._guess_net_charge_for_cofactor("UNKNOWN_XYZ") == 0


def test_write_single_cofactor_pdb_appends_end(tmp_path: Path):
    p = tmp_path / "one.pdb"
    _cofactor_params._write_single_cofactor_pdb(
        ["HETATM    1  N1  NAD A 501       0.0  0.0  0.0"], p,
    )
    text = p.read_text()
    assert text.endswith("END\n")


def test_parametrize_no_input_returns_skipped(tmp_path: Path):
    out = tmp_path / "out"
    result = _cofactor_params.parametrize_cofactors(
        tmp_path / "does_not_exist.pdb", out,
    )
    assert result["status"] == "skipped"
    assert "no cofactor PDB present" in result["message"]


def test_parametrize_empty_file_returns_skipped(tmp_path: Path):
    empty = tmp_path / "empty.pdb"
    empty.write_text("END\n")
    result = _cofactor_params.parametrize_cofactors(empty, tmp_path / "out")
    # Empty file has no HETATMs -> "no cofactors detected in file"
    assert result["status"] == "skipped"


def test_parametrize_known_library_cofactor_is_skipped_library(tmp_path: Path):
    inp = _write_pdb(tmp_path / "cof.pdb", [
        "HETATM    1  N1  NAD A 501       0.0  0.0  0.0",
        "HETATM    2  N2  NAD A 501       1.0  0.0  0.0",
    ])
    result = _cofactor_params.parametrize_cofactors(inp, tmp_path / "out")
    assert result["status"] == "skipped"  # all library-known
    assert len(result["cofactors"]) == 1
    assert result["cofactors"][0]["resname"] == "NAD"
    assert result["cofactors"][0]["status"] == "skipped_library"
    assert "prefer that over AM1-BCC" in result["cofactors"][0]["message"]


def test_parametrize_novel_cofactor_fails_open_without_antechamber(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch,
):
    """When antechamber is not on PATH the wrapper reports 'failed' with a
    helpful message but never raises."""
    # Force which() to return None so the pipeline believes antechamber missing.
    import shutil as _sh
    monkeypatch.setattr(_sh, "which", lambda name: None)

    inp = _write_pdb(tmp_path / "cof.pdb", [
        "HETATM    1  C1  XYZ A 501       0.0  0.0  0.0",
        "HETATM    2  C2  XYZ A 501       1.5  0.0  0.0",
    ])
    result = _cofactor_params.parametrize_cofactors(inp, tmp_path / "out")
    assert result["status"] == "failed"
    assert len(result["cofactors"]) == 1
    cof = result["cofactors"][0]
    assert cof["resname"] == "XYZ"
    assert cof["status"] == "failed"
    assert "antechamber not on PATH" in cof["message"]


def test_parametrize_mixed_library_and_novel(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch,
):
    """One known-library + one novel; expected status 'failed' because the
    novel one cannot be parametrized without antechamber, but the summary
    counts each category correctly."""
    import shutil as _sh
    monkeypatch.setattr(_sh, "which", lambda name: None)

    inp = _write_pdb(tmp_path / "mix.pdb", [
        "HETATM    1  N1  NAD A 501       0.0  0.0  0.0",  # library
        "HETATM    2  C1  XYZ B 601       0.0  0.0  0.0",  # novel
    ])
    result = _cofactor_params.parametrize_cofactors(inp, tmp_path / "out")
    assert len(result["cofactors"]) == 2
    library = next(c for c in result["cofactors"] if c["resname"] == "NAD")
    novel = next(c for c in result["cofactors"] if c["resname"] == "XYZ")
    assert library["status"] == "skipped_library"
    assert novel["status"] == "failed"
    assert result["summary"]["n_skipped_library"] == 1
    assert result["summary"]["n_failed"] == 1
    assert result["summary"]["n_success"] == 0


def test_known_library_table_covers_common_cofactors():
    """Sanity: the library table should cover the most-cited cofactors so
    users don't accidentally re-parametrize the standard ones."""
    must_have = ("NAD", "NAP", "FAD", "HEM", "SAM", "ATP")
    for cofactor in must_have:
        assert cofactor in _cofactor_params._KNOWN_LIBRARY_COFACTORS, (
            f"expected {cofactor} in known-library table"
        )
