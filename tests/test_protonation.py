from __future__ import annotations

import subprocess
from pathlib import Path
from typing import Any

import pytest

from stack_protein_preparation.protonation import (
    build_gromacs_position_restraints_output_path,
    build_gromacs_topology_output_path,
    build_standard_protonation_output_path,
    build_variant_protonation_output_path,
    count_atoms_in_pdb,
    count_atoms_in_structure_file,
    protonate_protein_structure,
    protonate_selected_structure,
    protonate_variant_structure,
    run_gmx_pdb2gmx_protonation,
    sanitize_variant_label,
    select_protonation_input,
)


def _make_atom_line(
    serial: int,
    atom_name: str,
    residue_name: str,
    chain_id: str,
    residue_number: int,
    x: float,
    y: float,
    z: float,
    element: str,
    record_name: str = "ATOM",
) -> str:
    return (
        f"{record_name:<6}"
        f"{serial:5d} "
        f"{atom_name:>4}"
        f" "
        f"{residue_name:>3} "
        f"{chain_id:1}"
        f"{residue_number:4d}    "
        f"{x:8.3f}{y:8.3f}{z:8.3f}"
        f"{1.00:6.2f}{20.00:6.2f}          "
        f"{element:>2}"
    )


def _write_structure(path: Path, lines: list[str]) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def _write_text_file(path: str | Path, text: str) -> Path:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")
    return path


def _build_minimal_protein(path: Path, residue_name: str = "ALA") -> Path:
    lines = [
        _make_atom_line(1, "N", residue_name, "A", 1, 0.000, 0.000, 0.000, "N"),
        _make_atom_line(2, "CA", residue_name, "A", 1, 1.460, 0.000, 0.000, "C"),
        _make_atom_line(3, "C", residue_name, "A", 1, 2.900, 0.000, 0.000, "C"),
        _make_atom_line(4, "O", residue_name, "A", 1, 3.900, 0.000, 0.000, "O"),
        "TER",
        "END",
    ]
    return _write_structure(path, lines)


def _build_protonated_like_output(path: str | Path) -> Path:
    path = Path(path)
    lines = [
        _make_atom_line(1, "N", "ALA", "A", 1, 0.000, 0.000, 0.000, "N"),
        _make_atom_line(2, "H", "ALA", "A", 1, -0.700, 0.000, 0.000, "H"),
        _make_atom_line(3, "CA", "ALA", "A", 1, 1.460, 0.000, 0.000, "C"),
        _make_atom_line(4, "HA", "ALA", "A", 1, 1.460, 0.900, 0.000, "H"),
        _make_atom_line(5, "C", "ALA", "A", 1, 2.900, 0.000, 0.000, "C"),
        _make_atom_line(6, "O", "ALA", "A", 1, 3.900, 0.000, 0.000, "O"),
        "TER",
        "END",
    ]
    return _write_structure(path, lines)


def _write_fake_gromacs_outputs(
    *,
    output_pdb: str | Path,
    topology_output_path: str | Path | None,
    position_restraints_output_path: str | Path | None,
) -> None:
    _build_protonated_like_output(Path(output_pdb))

    if topology_output_path is not None:
        _write_text_file(topology_output_path, "; fake topology\n")

    if position_restraints_output_path is not None:
        _write_text_file(position_restraints_output_path, "; fake position restraints\n")


def test_count_atoms_in_structure_file_counts_atom_and_hetatm(tmp_path: Path) -> None:
    pdb_path = tmp_path / "count_test.pdb"
    lines = [
        _make_atom_line(1, "N", "ALA", "A", 1, 0, 0, 0, "N", "ATOM"),
        _make_atom_line(2, "CA", "ALA", "A", 1, 1, 0, 0, "C", "ATOM"),
        _make_atom_line(3, "O", "HOH", "A", 2, 2, 0, 0, "O", "HETATM"),
        "TER",
        "END",
    ]
    _write_structure(pdb_path, lines)

    assert count_atoms_in_structure_file(pdb_path) == 3
    assert count_atoms_in_pdb(pdb_path) == 3


def test_count_atoms_in_structure_file_raises_for_missing_file(tmp_path: Path) -> None:
    with pytest.raises(FileNotFoundError):
        count_atoms_in_structure_file(tmp_path / "missing.pdb")


def test_select_protonation_input_prefers_explicit_filler_over_other_inputs(
    tmp_path: Path,
) -> None:
    protein_dir = tmp_path / "1ABC"
    components_dir = protein_dir / "components"

    _build_minimal_protein(components_dir / "1ABC_protein.pdb")
    alphafold_model = _build_minimal_protein(tmp_path / "af_model.pdb")
    modeller_model = _build_minimal_protein(tmp_path / "modeller_model.pdb")
    filler_model = _build_minimal_protein(tmp_path / "filler_final.pdb")

    selected_path, source = select_protonation_input(
        pdb_id="1ABC",
        protein_dir=protein_dir,
        filler_final_model_path=filler_model,
        modeller_model_path=modeller_model,
        alphafold_model_path=alphafold_model,
    )

    assert selected_path == filler_model
    assert source == "filler"


def test_select_protonation_input_uses_default_filler_final_model(
    tmp_path: Path,
) -> None:
    protein_dir = tmp_path / "1ABC"
    components_dir = protein_dir / "components"
    filler_output_dir = tmp_path / "filler"

    _build_minimal_protein(components_dir / "1ABC_protein.pdb")
    default_filler_model = _build_minimal_protein(
        filler_output_dir / "final_filled_model.pdb"
    )

    selected_path, source = select_protonation_input(
        pdb_id="1ABC",
        protein_dir=protein_dir,
        filler_output_dir=filler_output_dir,
    )

    assert selected_path == default_filler_model
    assert source == "filler"


def test_select_protonation_input_prefers_modeller_over_alphafold_and_protein(
    tmp_path: Path,
) -> None:
    protein_dir = tmp_path / "1ABC"
    components_dir = protein_dir / "components"
    default_protein = _build_minimal_protein(components_dir / "1ABC_protein.pdb")
    alphafold_model = _build_minimal_protein(tmp_path / "af_model.pdb")
    modeller_model = _build_minimal_protein(tmp_path / "modeller_model.pdb")

    selected_path, source = select_protonation_input(
        pdb_id="1ABC",
        protein_dir=protein_dir,
        modeller_model_path=modeller_model,
        alphafold_model_path=alphafold_model,
    )

    assert selected_path == modeller_model
    assert source == "modeller"
    assert default_protein.exists()


def test_select_protonation_input_prefers_alphafold_over_default_when_no_modeller(
    tmp_path: Path,
) -> None:
    protein_dir = tmp_path / "1ABC"
    components_dir = protein_dir / "components"
    _build_minimal_protein(components_dir / "1ABC_protein.pdb")
    alphafold_model = _build_minimal_protein(tmp_path / "af_model.pdb")

    selected_path, source = select_protonation_input(
        pdb_id="1ABC",
        protein_dir=protein_dir,
        modeller_model_path=None,
        alphafold_model_path=alphafold_model,
    )

    assert selected_path == alphafold_model
    assert source == "alphafold"


def test_select_protonation_input_falls_back_to_default_protein(tmp_path: Path) -> None:
    protein_dir = tmp_path / "1ABC"
    components_dir = protein_dir / "components"
    default_protein = _build_minimal_protein(components_dir / "1ABC_protein.pdb")

    selected_path, source = select_protonation_input(
        pdb_id="1ABC",
        protein_dir=protein_dir,
        modeller_model_path=None,
        alphafold_model_path=None,
    )

    assert selected_path == default_protein
    assert source == "protein"


def test_select_protonation_input_raises_if_no_valid_input_exists(
    tmp_path: Path,
) -> None:
    with pytest.raises(FileNotFoundError):
        select_protonation_input(
            pdb_id="1ABC",
            protein_dir=tmp_path / "1ABC",
            modeller_model_path=None,
            alphafold_model_path=None,
        )


def test_sanitize_variant_label_normalizes_symbols() -> None:
    assert sanitize_variant_label("modeller complete") == "modeller_complete"
    assert sanitize_variant_label("alpha/fold") == "alpha_fold"
    assert sanitize_variant_label("gaps") == "gaps"


def test_sanitize_variant_label_raises_for_empty_result() -> None:
    with pytest.raises(ValueError):
        sanitize_variant_label("///")


def test_build_standard_protonation_output_path_returns_expected_path(
    tmp_path: Path,
) -> None:
    output_path = build_standard_protonation_output_path(
        pdb_id="1ABC",
        protein_dir=tmp_path / "1ABC",
    )

    assert output_path == tmp_path / "1ABC" / "components" / "1ABC_proteinH.pdb"


def test_build_variant_protonation_output_path_returns_expected_path(
    tmp_path: Path,
) -> None:
    output_path = build_variant_protonation_output_path(
        pdb_id="1ABC",
        protein_dir=tmp_path / "1ABC",
        variant_label="modeller complete",
    )

    assert output_path == (
        tmp_path / "1ABC" / "components" / "1ABC_proteinH_modeller_complete.pdb"
    )


def test_build_gromacs_output_paths_are_paired_with_coordinate_output(
    tmp_path: Path,
) -> None:
    output_pdb = tmp_path / "1ABC" / "components" / "1ABC_proteinH.pdb"

    assert build_gromacs_topology_output_path(output_pdb) == (
        tmp_path / "1ABC" / "components" / "1ABC_proteinH.top"
    )
    assert build_gromacs_position_restraints_output_path(output_pdb) == (
        tmp_path / "1ABC" / "components" / "1ABC_proteinH_posre.itp"
    )


def test_run_gmx_pdb2gmx_protonation_constructs_expected_command(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    input_pdb = _build_minimal_protein(tmp_path / "input.pdb")
    output_pdb = tmp_path / "outputs" / "proteinH.pdb"

    observed: dict[str, Any] = {}

    def fake_run(
        cmd: list[str],
        check: bool,
        capture_output: bool,
        text: bool,
        input: str,
        cwd: str,
        timeout: int,
        env: dict[str, str],
    ) -> subprocess.CompletedProcess[str]:
        observed["cmd"] = cmd
        observed["check"] = check
        observed["capture_output"] = capture_output
        observed["text"] = text
        observed["input"] = input
        observed["cwd"] = cwd
        observed["timeout"] = timeout
        observed["env"] = env
        return subprocess.CompletedProcess(
            args=cmd,
            returncode=0,
            stdout="gmx stdout",
            stderr="gmx stderr",
        )

    monkeypatch.setattr(
        "stack_protein_preparation._protonation_core._find_gmx_executable",
        lambda: "gmx",
    )
    monkeypatch.setattr(
        "stack_protein_preparation._protonation_core.subprocess.run",
        fake_run,
    )

    result = run_gmx_pdb2gmx_protonation(
        input_pdb=input_pdb,
        output_pdb=output_pdb,
        ff="amber99sb-ildn",
        water_model="tip3p",
        timeout_seconds=123,
    )

    assert result.returncode == 0
    assert observed["cmd"] == [
        "gmx",
        "pdb2gmx",
        "-f",
        str(input_pdb.resolve()),
        "-o",
        "proteinH.pdb",
        "-p",
        "proteinH.top",
        "-i",
        "proteinH_posre.itp",
        "-ff",
        "amber99sb-ildn",
        "-water",
        "tip3p",
        "-chainsep",
        "id_or_ter",
        "-merge",
        "no",
        "-ignh",
    ]
    assert observed["cwd"] == str(output_pdb.parent)
    assert observed["input"] == ""
    assert observed["timeout"] == 123
    assert observed["env"]["GMX_MAXBACKUP"] == "-1"


def test_protonate_selected_structure_runs_on_explicit_input_and_returns_stats(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    input_pdb = _build_minimal_protein(tmp_path / "input.pdb")
    output_pdb = tmp_path / "outputH.pdb"

    def fake_run_gmx_pdb2gmx_protonation(
        input_pdb: str | Path,
        output_pdb: str | Path,
        *,
        topology_output_path: str | Path | None = None,
        position_restraints_output_path: str | Path | None = None,
        ff: str = "amber99sb-ildn",
        water_model: str = "tip3p",
        chain_separation: str = "id_or_ter",
        merge: str = "no",
        ignore_input_hydrogens: bool = True,
        timeout_seconds: int = 600,
    ) -> subprocess.CompletedProcess[str]:
        assert ff == "amber99sb-ildn"
        assert water_model == "tip3p"
        assert chain_separation == "id_or_ter"
        assert merge == "no"
        assert ignore_input_hydrogens is True
        _write_fake_gromacs_outputs(
            output_pdb=output_pdb,
            topology_output_path=topology_output_path,
            position_restraints_output_path=position_restraints_output_path,
        )
        return subprocess.CompletedProcess(
            args=["gmx", "pdb2gmx"],
            returncode=0,
            stdout="fake gmx stdout",
            stderr="fake gmx stderr",
        )

    monkeypatch.setattr(
        "stack_protein_preparation.protonation._find_gmx_executable",
        lambda: "gmx",
    )
    monkeypatch.setattr(
        "stack_protein_preparation.protonation.run_gmx_pdb2gmx_protonation",
        fake_run_gmx_pdb2gmx_protonation,
    )

    result = protonate_selected_structure(
        input_pdb=input_pdb,
        output_pdb=output_pdb,
        input_source="modeller",
        ph=6.5,
        variant_label="modeller_complete",
    )

    expected_topology = output_pdb.with_suffix(".top")
    expected_posre = output_pdb.with_name("outputH_posre.itp")

    assert output_pdb.exists()
    assert expected_topology.exists()
    assert expected_posre.exists()
    assert result["protonation_success"] is True
    assert result["protonation_method"] == "gmx pdb2gmx -ignh"
    assert result["protonation_input_source"] == "modeller"
    assert result["protonation_input_path"] == str(input_pdb)
    assert result["protonation_output_path"] == str(output_pdb)
    assert result["protonation_topology_path"] == str(expected_topology)
    assert result["protonation_position_restraints_path"] == str(expected_posre)
    assert result["protonation_ph"] == 6.5
    assert result["protonation_ph_applied_by_gromacs"] is False
    assert result["gromacs_force_field"] == "amber99sb-ildn"
    assert result["gromacs_water_model"] == "tip3p"
    assert result["gromacs_ignore_input_hydrogens"] is True
    assert result["input_atom_count"] == 4
    assert result["output_atom_count"] == 6
    assert result["atom_count_increased"] is True
    assert result["gmx_pdb2gmx_returncode"] == 0
    assert result["gmx_pdb2gmx_stdout"] == "fake gmx stdout"
    assert result["gmx_pdb2gmx_stderr"] == "fake gmx stderr"
    assert result["topology_exists"] is True
    assert result["topology_nonempty"] is True
    assert result["position_restraints_exists"] is True
    assert result["position_restraints_nonempty"] is True
    assert result["protonation_variant"] == "modeller_complete"


def test_protonate_selected_structure_returns_structured_failure_for_gmx_failure(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    input_pdb = _build_minimal_protein(tmp_path / "input.pdb")
    output_pdb = tmp_path / "outputH.pdb"

    def fake_run_gmx_pdb2gmx_protonation(
        input_pdb: str | Path,
        output_pdb: str | Path,
        *,
        topology_output_path: str | Path | None = None,
        position_restraints_output_path: str | Path | None = None,
        ff: str = "amber99sb-ildn",
        water_model: str = "tip3p",
        chain_separation: str = "id_or_ter",
        merge: str = "no",
        ignore_input_hydrogens: bool = True,
        timeout_seconds: int = 600,
    ) -> subprocess.CompletedProcess[str]:
        return subprocess.CompletedProcess(
            args=["gmx", "pdb2gmx"],
            returncode=1,
            stdout="gmx failed stdout",
            stderr="Fatal error: residue not found",
        )

    monkeypatch.setattr(
        "stack_protein_preparation.protonation._find_gmx_executable",
        lambda: "gmx",
    )
    monkeypatch.setattr(
        "stack_protein_preparation.protonation.run_gmx_pdb2gmx_protonation",
        fake_run_gmx_pdb2gmx_protonation,
    )

    result = protonate_selected_structure(
        input_pdb=input_pdb,
        output_pdb=output_pdb,
        input_source="protein",
    )

    assert result["protonation_success"] is False
    assert result["gmx_pdb2gmx_returncode"] == 1
    assert result["output_atom_count"] == 0
    assert result["atom_count_increased"] is False
    assert "returncode=1" in str(result["protonation_error_message"])
    assert "Fatal error: residue not found" in str(result["protonation_error_message"])


def test_protonate_selected_structure_raises_for_missing_explicit_input(
    tmp_path: Path,
) -> None:
    with pytest.raises(FileNotFoundError):
        protonate_selected_structure(
            input_pdb=tmp_path / "missing.pdb",
            output_pdb=tmp_path / "out.pdb",
            input_source="protein",
        )


def test_protonate_variant_structure_writes_variant_specific_output(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    protein_dir = tmp_path / "1ABC"
    input_pdb = _build_minimal_protein(tmp_path / "variant_input.pdb")

    def fake_run_gmx_pdb2gmx_protonation(
        input_pdb: str | Path,
        output_pdb: str | Path,
        *,
        topology_output_path: str | Path | None = None,
        position_restraints_output_path: str | Path | None = None,
        ff: str = "amber99sb-ildn",
        water_model: str = "tip3p",
        chain_separation: str = "id_or_ter",
        merge: str = "no",
        ignore_input_hydrogens: bool = True,
        timeout_seconds: int = 600,
    ) -> subprocess.CompletedProcess[str]:
        _write_fake_gromacs_outputs(
            output_pdb=output_pdb,
            topology_output_path=topology_output_path,
            position_restraints_output_path=position_restraints_output_path,
        )
        return subprocess.CompletedProcess(
            args=["gmx", "pdb2gmx"],
            returncode=0,
            stdout="variant stdout",
            stderr="",
        )

    monkeypatch.setattr(
        "stack_protein_preparation.protonation._find_gmx_executable",
        lambda: "gmx",
    )
    monkeypatch.setattr(
        "stack_protein_preparation.protonation.run_gmx_pdb2gmx_protonation",
        fake_run_gmx_pdb2gmx_protonation,
    )

    result = protonate_variant_structure(
        pdb_id="1ABC",
        protein_dir=protein_dir,
        variant_label="alphafold_complete",
        input_pdb=input_pdb,
        input_source="alphafold",
        ph=7.4,
    )

    expected_output = (
        protein_dir / "components" / "1ABC_proteinH_alphafold_complete.pdb"
    )

    assert expected_output.exists()
    assert expected_output.with_suffix(".top").exists()
    assert expected_output.with_name("1ABC_proteinH_alphafold_complete_posre.itp").exists()
    assert result["protonation_success"] is True
    assert result["protonation_input_source"] == "alphafold"
    assert result["protonation_output_path"] == str(expected_output)
    assert result["protonation_variant"] == "alphafold_complete"


def test_protonate_protein_structure_keeps_backward_compatible_single_output(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    pdb_id = "1ABC"
    protein_dir = tmp_path / pdb_id
    components_dir = protein_dir / "components"
    default_protein = _build_minimal_protein(components_dir / f"{pdb_id}_protein.pdb")

    def fake_run_gmx_pdb2gmx_protonation(
        input_pdb: str | Path,
        output_pdb: str | Path,
        *,
        topology_output_path: str | Path | None = None,
        position_restraints_output_path: str | Path | None = None,
        ff: str = "amber99sb-ildn",
        water_model: str = "tip3p",
        chain_separation: str = "id_or_ter",
        merge: str = "no",
        ignore_input_hydrogens: bool = True,
        timeout_seconds: int = 600,
    ) -> subprocess.CompletedProcess[str]:
        assert Path(input_pdb) == default_protein
        _write_fake_gromacs_outputs(
            output_pdb=output_pdb,
            topology_output_path=topology_output_path,
            position_restraints_output_path=position_restraints_output_path,
        )
        return subprocess.CompletedProcess(
            args=["gmx", "pdb2gmx"],
            returncode=0,
            stdout="backward stdout",
            stderr="backward stderr",
        )

    monkeypatch.setattr(
        "stack_protein_preparation.protonation._find_gmx_executable",
        lambda: "gmx",
    )
    monkeypatch.setattr(
        "stack_protein_preparation.protonation.run_gmx_pdb2gmx_protonation",
        fake_run_gmx_pdb2gmx_protonation,
    )

    result = protonate_protein_structure(
        pdb_id=pdb_id,
        protein_dir=protein_dir,
        modeller_model_path=None,
        alphafold_model_path=None,
        ph=7.2,
    )

    expected_output = components_dir / f"{pdb_id}_proteinH.pdb"

    assert expected_output.exists()
    assert expected_output.with_suffix(".top").exists()
    assert expected_output.with_name(f"{pdb_id}_proteinH_posre.itp").exists()
    assert result["protonation_success"] is True
    assert result["protonation_input_source"] == "protein"
    assert result["protonation_input_path"] == str(default_protein)
    assert result["protonation_output_path"] == str(expected_output)
    assert result["protonation_ph"] == 7.2
    assert result["protonation_ph_applied_by_gromacs"] is False
    assert result["input_atom_count"] == 4
    assert result["output_atom_count"] == 6
    assert result["atom_count_increased"] is True