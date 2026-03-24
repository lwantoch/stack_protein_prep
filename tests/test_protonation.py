# /home/grheco/repositorios/stack_protein_prep/tests/test_protonation.py

from __future__ import annotations

from pathlib import Path

import pytest

from stack_protein_preparation.protonation import (
    _find_pdb2pqr_executable,
    protonate_protein_structure,
    run_pdb2pqr_protonation,
    select_protonation_input,
)


def _write_text_file(path: Path, content: str = "ATOM\n") -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(content, encoding="utf-8")
    return path


# =========================
# selection logic
# =========================


def test_select_protonation_input_prefers_modeller_over_alphafold_and_protein(
    tmp_path: Path,
) -> None:
    pdb_id = "1ABC"
    protein_dir = tmp_path / pdb_id

    protein_path = _write_text_file(
        protein_dir / "components" / f"{pdb_id}_protein.pdb"
    )
    alphafold_path = _write_text_file(
        protein_dir / "alphafold" / f"{pdb_id}_alphafold_model.pdb"
    )
    modeller_path = _write_text_file(
        protein_dir / "MODELLER" / "chain_A" / f"{pdb_id}_A_target_model.pdb"
    )

    selected_path, source = select_protonation_input(
        pdb_id=pdb_id,
        protein_dir=protein_dir,
        modeller_model_path=modeller_path,
        alphafold_model_path=alphafold_path,
    )

    assert selected_path == modeller_path
    assert source == "modeller"
    assert protein_path.is_file()


def test_select_protonation_input_uses_alphafold_when_modeller_missing(
    tmp_path: Path,
) -> None:
    pdb_id = "2DEF"
    protein_dir = tmp_path / pdb_id

    _write_text_file(protein_dir / "components" / f"{pdb_id}_protein.pdb")
    alphafold_path = _write_text_file(
        protein_dir / "alphafold" / f"{pdb_id}_alphafold_model.pdb"
    )

    selected_path, source = select_protonation_input(
        pdb_id=pdb_id,
        protein_dir=protein_dir,
        modeller_model_path=protein_dir / "MODELLER" / "missing_model.pdb",
        alphafold_model_path=alphafold_path,
    )

    assert selected_path == alphafold_path
    assert source == "alphafold"


def test_select_protonation_input_uses_original_protein_when_no_models_exist(
    tmp_path: Path,
) -> None:
    pdb_id = "3GHI"
    protein_dir = tmp_path / pdb_id

    protein_path = _write_text_file(
        protein_dir / "components" / f"{pdb_id}_protein.pdb"
    )

    selected_path, source = select_protonation_input(
        pdb_id=pdb_id,
        protein_dir=protein_dir,
        modeller_model_path=None,
        alphafold_model_path=None,
    )

    assert selected_path == protein_path
    assert source == "protein"


def test_select_protonation_input_ignores_empty_modeller_and_uses_alphafold(
    tmp_path: Path,
) -> None:
    pdb_id = "4JKL"
    protein_dir = tmp_path / pdb_id

    empty_modeller = protein_dir / "MODELLER" / "chain_A" / f"{pdb_id}_empty_model.pdb"
    empty_modeller.parent.mkdir(parents=True, exist_ok=True)
    empty_modeller.write_text("", encoding="utf-8")

    alphafold_path = _write_text_file(
        protein_dir / "alphafold" / f"{pdb_id}_alphafold_model.pdb"
    )
    _write_text_file(protein_dir / "components" / f"{pdb_id}_protein.pdb")

    selected_path, source = select_protonation_input(
        pdb_id=pdb_id,
        protein_dir=protein_dir,
        modeller_model_path=empty_modeller,
        alphafold_model_path=alphafold_path,
    )

    assert selected_path == alphafold_path
    assert source == "alphafold"


def test_select_protonation_input_raises_if_no_valid_input_exists(
    tmp_path: Path,
) -> None:
    pdb_id = "5MNO"
    protein_dir = tmp_path / pdb_id

    with pytest.raises(FileNotFoundError, match="No valid protonation input found"):
        select_protonation_input(
            pdb_id=pdb_id,
            protein_dir=protein_dir,
            modeller_model_path=None,
            alphafold_model_path=None,
        )


# =========================
# pdb2pqr executable
# =========================


def test_find_pdb2pqr_executable_prefers_first_available(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    def fake_which(name: str) -> str | None:
        if name == "pdb2pqr":
            return "/usr/bin/pdb2pqr"
        if name == "pdb2pqr30":
            return "/usr/bin/pdb2pqr30"
        return None

    monkeypatch.setattr(
        "stack_protein_preparation.protonation.shutil.which", fake_which
    )

    executable = _find_pdb2pqr_executable()

    assert executable == "pdb2pqr"


def test_find_pdb2pqr_executable_raises_if_not_found(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    monkeypatch.setattr(
        "stack_protein_preparation.protonation.shutil.which",
        lambda _name: None,
    )

    with pytest.raises(FileNotFoundError, match="Could not find a pdb2pqr executable"):
        _find_pdb2pqr_executable()


# =========================
# command construction
# =========================


def test_run_pdb2pqr_protonation_builds_expected_command(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    input_pdb = _write_text_file(tmp_path / "input.pdb")
    output_pdb = tmp_path / "nested" / "output.pdb"

    monkeypatch.setattr(
        "stack_protein_preparation.protonation._find_pdb2pqr_executable",
        lambda: "pdb2pqr",
    )

    captured: dict[str, object] = {}

    class DummyCompletedProcess:
        def __init__(self) -> None:
            self.stdout = "stdout text"
            self.stderr = "stderr text"

    def fake_run(
        cmd: list[str],
        check: bool,
        capture_output: bool,
        text: bool,
    ) -> DummyCompletedProcess:
        captured["cmd"] = cmd
        captured["check"] = check
        captured["capture_output"] = capture_output
        captured["text"] = text
        return DummyCompletedProcess()

    monkeypatch.setattr(
        "stack_protein_preparation.protonation.subprocess.run", fake_run
    )

    result = run_pdb2pqr_protonation(
        input_pdb=input_pdb,
        output_pdb=output_pdb,
        ph=6.8,
        ff="AMBER",
        keep_chain=True,
        keep_heterogens=False,
    )

    assert output_pdb.parent.is_dir()
    assert captured["cmd"] == [
        "pdb2pqr",
        "--ff=AMBER",
        "--titration-state-method=propka",
        "--with-ph=6.8",
        "--keep-chain",
        str(input_pdb),
        str(output_pdb),
    ]
    assert result.stdout == "stdout text"
    assert result.stderr == "stderr text"


def test_run_pdb2pqr_protonation_can_keep_heterogens(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    input_pdb = _write_text_file(tmp_path / "input.pdb")
    output_pdb = tmp_path / "output.pdb"

    monkeypatch.setattr(
        "stack_protein_preparation.protonation._find_pdb2pqr_executable",
        lambda: "pdb2pqr30",
    )

    captured: dict[str, list[str]] = {}

    class DummyCompletedProcess:
        def __init__(self) -> None:
            self.stdout = ""
            self.stderr = ""

    def fake_run(
        cmd: list[str],
        check: bool,
        capture_output: bool,
        text: bool,
    ) -> DummyCompletedProcess:
        captured["cmd"] = cmd
        return DummyCompletedProcess()

    monkeypatch.setattr(
        "stack_protein_preparation.protonation.subprocess.run", fake_run
    )

    run_pdb2pqr_protonation(
        input_pdb=input_pdb,
        output_pdb=output_pdb,
        ph=7.4,
        ff="AMBER",
        keep_chain=False,
        keep_heterogens=True,
    )

    assert captured["cmd"] == [
        "pdb2pqr30",
        "--ff=AMBER",
        "--titration-state-method=propka",
        "--with-ph=7.4",
        "--keep-heterogens",
        str(input_pdb),
        str(output_pdb),
    ]


# =========================
# protonation main logic
# =========================


def test_protonate_protein_structure_returns_expected_metadata_for_modeller(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    pdb_id = "6PQR"
    protein_dir = tmp_path / pdb_id

    modeller_path = _write_text_file(
        protein_dir / "MODELLER" / "chain_A" / f"{pdb_id}_A_target_model.pdb"
    )

    def fake_run_pdb2pqr_protonation(*args, **kwargs):
        output_pdb = Path(kwargs["output_pdb"])
        output_pdb.parent.mkdir(parents=True, exist_ok=True)
        output_pdb.write_text("ATOM\nATOM\n", encoding="utf-8")

        class Dummy:
            stdout = "ok"
            stderr = ""

        return Dummy()

    monkeypatch.setattr(
        "stack_protein_preparation.protonation.run_pdb2pqr_protonation",
        fake_run_pdb2pqr_protonation,
    )

    monkeypatch.setattr(
        "stack_protein_preparation.protonation.count_atoms_in_pdb",
        lambda p: 2,
    )

    result = protonate_protein_structure(
        pdb_id=pdb_id,
        protein_dir=protein_dir,
        modeller_model_path=modeller_path,
        alphafold_model_path=None,
    )

    assert result["protonation_input_source"] == "modeller"
    assert result["protonation_success"] is True


def test_protonate_protein_structure_atom_count_increases(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    pdb_id = "TEST"
    protein_dir = tmp_path / pdb_id

    _write_text_file(
        protein_dir / "components" / f"{pdb_id}_protein.pdb",
        content="ATOM\nATOM\nATOM\n",
    )

    def fake_run_pdb2pqr_protonation(*args, **kwargs):
        output_pdb = Path(kwargs["output_pdb"])
        output_pdb.parent.mkdir(parents=True, exist_ok=True)
        output_pdb.write_text(
            "ATOM\nATOM\nATOM\nATOM\nATOM\n",
            encoding="utf-8",
        )

        class Dummy:
            stdout = ""
            stderr = ""

        return Dummy()

    monkeypatch.setattr(
        "stack_protein_preparation.protonation.run_pdb2pqr_protonation",
        fake_run_pdb2pqr_protonation,
    )

    counts = {"input": 3, "output": 5}

    def fake_count_atoms(pdb_path):
        return counts["input"] if "protein.pdb" in str(pdb_path) else counts["output"]

    monkeypatch.setattr(
        "stack_protein_preparation.protonation.count_atoms_in_pdb",
        fake_count_atoms,
    )

    result = protonate_protein_structure(
        pdb_id=pdb_id,
        protein_dir=protein_dir,
    )

    assert result["input_atom_count"] == 3
    assert result["output_atom_count"] == 5
    assert result["atom_count_increased"] is True
    assert result["protonation_success"] is True
