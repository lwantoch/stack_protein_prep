"""Tests for metal_hydrogen_cleanup and metall_params using synthetic PDB strings.

All tests use in-memory PDB content written to tmp_path — no external tools
(Chimera, Gaussian, MCPB.py, AmberTools, GROMACS) are required.
"""
from __future__ import annotations

import csv
import json
from pathlib import Path

import pytest

from stack_protein_preparation.metal_hydrogen_cleanup import (
    remove_metal_coordinating_hydrogens,
)


# ---------------------------------------------------------------------------
# PDB line helper
# ---------------------------------------------------------------------------

def _pdb_atom_line(
    serial: int,
    name: str,
    resname: str,
    chain: str,
    resseq: int,
    x: float,
    y: float,
    z: float,
    element: str,
    record: str = "ATOM",
) -> str:
    """Return a correctly formatted 80-column PDB ATOM/HETATM line."""
    rec = record.ljust(6)[:6]
    # 4-char atom name: single-char elements are padded with a leading space
    if len(element) == 1:
        atom_name = f" {name:<3s}"
    else:
        atom_name = f"{name:<4s}"
    elem = element.rjust(2)[:2]
    return (
        f"{rec}{serial:5d} {atom_name} {resname:>3s} {chain}{resseq:4d}    "
        f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00"
        f"          {elem}  \n"
    )


def _write_pdb(tmp_path: Path, name: str, lines: list[str]) -> Path:
    p = tmp_path / name
    p.write_text("".join(lines) + "END\n", encoding="utf-8")
    return p


# ---------------------------------------------------------------------------
# Test 1: HIS NE2-H pointing toward Zn is removed
# ---------------------------------------------------------------------------

def test_his_ne2_h_toward_zn_is_removed(tmp_path: Path) -> None:
    """HE2 on HIS NE2 pointing directly toward Zn must be removed.

    Geometry:
        Zn at (0,0,0), NE2 at (2,0,0), HE2 at (1.2,0,0)
        angle M-NE2-HE2 = 0° ≤ 60° → removal criterion met
        dist(Zn,HE2) = 1.2 Å ≤ 2.10 Å → also meets absolute dist criterion
    """
    lines = [
        _pdb_atom_line(1, "NE2", "HIS", "A", 1, 2.0, 0.0, 0.0, "N"),
        _pdb_atom_line(2, "HE2", "HIS", "A", 1, 1.2, 0.0, 0.0, "H"),
        _pdb_atom_line(100, "ZN", "ZN", "A", 301, 0.0, 0.0, 0.0, "ZN", record="HETATM"),
    ]
    pdb = _write_pdb(tmp_path, "test1.pdb", lines)
    out = tmp_path / "out1"

    result = remove_metal_coordinating_hydrogens(
        input_pdb_path=pdb,
        output_dir=out,
        pdb_id="TEST",
        variant_label="test1",
    )

    assert result.status == "success"
    assert result.removed_count == 1
    assert result.transition_metal_count == 1

    # HE2 (serial 2) must be absent from the cleaned PDB
    cleaned_lines = result.cleaned_pdb_path.read_text().splitlines()
    atom_serials_in_cleaned = {
        int(line[6:11]) for line in cleaned_lines if line.startswith(("ATOM", "HETATM"))
    }
    assert 2 not in atom_serials_in_cleaned

    # Removed TSV must record the removal
    with result.removed_tsv_path.open(newline="") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    assert len(rows) == 1
    assert rows[0]["metal_element"] == "ZN"
    assert rows[0]["donor_atom_name"] == "NE2"
    assert rows[0]["h_atom_name"] == "HE2"


# ---------------------------------------------------------------------------
# Test 2: HIS hydrogen pointing away from Zn is retained
# ---------------------------------------------------------------------------

def test_his_h_pointing_away_from_zn_is_retained(tmp_path: Path) -> None:
    """HD1 on HIS ND1 that points directly away from Zn must be kept.

    Geometry:
        Zn at (0,0,0), ND1 at (2,0,0), HD1 at (3,0,0)
        angle M-ND1-HD1 = 180° > 60° → angle criterion NOT met
        dist(Zn,HD1) = 3.0 Å > 2.10 Å → absolute dist NOT met
        slack: 3.0 ≥ 2.0 + 0.35 = 2.35 → slack criterion NOT met
        → HD1 retained
    """
    lines = [
        _pdb_atom_line(1, "ND1", "HIS", "A", 1, 2.0, 0.0, 0.0, "N"),
        _pdb_atom_line(2, "HD1", "HIS", "A", 1, 3.0, 0.0, 0.0, "H"),
        _pdb_atom_line(100, "ZN", "ZN", "A", 301, 0.0, 0.0, 0.0, "ZN", record="HETATM"),
    ]
    pdb = _write_pdb(tmp_path, "test2.pdb", lines)
    out = tmp_path / "out2"

    result = remove_metal_coordinating_hydrogens(
        input_pdb_path=pdb,
        output_dir=out,
        pdb_id="TEST",
        variant_label="test2",
    )

    assert result.status == "success"
    assert result.removed_count == 0
    assert result.kept_count >= 1

    # HD1 (serial 2) must still be present
    cleaned_lines = result.cleaned_pdb_path.read_text().splitlines()
    atom_serials_in_cleaned = {
        int(line[6:11]) for line in cleaned_lines if line.startswith(("ATOM", "HETATM"))
    }
    assert 2 in atom_serials_in_cleaned


# ---------------------------------------------------------------------------
# Test 3: Water hydrogens near Zn are retained by default
# ---------------------------------------------------------------------------

def test_water_hydrogens_near_zn_retained_by_default(tmp_path: Path) -> None:
    """HOH hydrogens within Zn range must not be removed when keep_water_hydrogens=True.

    H1 at (1.2,0,0) passes the angle and distance removal criteria, but the
    water guard skips the whole HOH residue before hydrogen evaluation.
    """
    lines = [
        _pdb_atom_line(1, "O", "HOH", "A", 1, 2.0, 0.0, 0.0, "O", record="HETATM"),
        _pdb_atom_line(2, "H1", "HOH", "A", 1, 1.2, 0.0, 0.0, "H", record="HETATM"),
        _pdb_atom_line(3, "H2", "HOH", "A", 1, 2.5, 0.8, 0.0, "H", record="HETATM"),
        _pdb_atom_line(100, "ZN", "ZN", "A", 301, 0.0, 0.0, 0.0, "ZN", record="HETATM"),
    ]
    pdb = _write_pdb(tmp_path, "test3.pdb", lines)
    out = tmp_path / "out3"

    result = remove_metal_coordinating_hydrogens(
        input_pdb_path=pdb,
        output_dir=out,
        pdb_id="TEST",
        variant_label="test3",
        keep_water_hydrogens=True,
    )

    assert result.removed_count == 0

    # Both water H atoms must still be present
    cleaned_lines = result.cleaned_pdb_path.read_text().splitlines()
    atom_serials_in_cleaned = {
        int(line[6:11]) for line in cleaned_lines if line.startswith(("ATOM", "HETATM"))
    }
    assert 2 in atom_serials_in_cleaned
    assert 3 in atom_serials_in_cleaned


# ---------------------------------------------------------------------------
# Test 4: Far donor hydrogen is retained
# ---------------------------------------------------------------------------

def test_far_donor_hydrogen_retained(tmp_path: Path) -> None:
    """CYS SG-HG pointing away from Zn at safe distance must not be removed.

    Geometry:
        Zn at (0,0,0), SG at (2,0,0), HG at (3,0,0)
        angle M-SG-HG = 180° > 60° → not removed
        dist(Zn,HG) = 3.0 Å ≥ 2.35 Å (2.0+0.35) → slack not met
        → HG retained
    """
    lines = [
        _pdb_atom_line(1, "SG", "CYS", "A", 1, 2.0, 0.0, 0.0, "S"),
        _pdb_atom_line(2, "HG", "CYS", "A", 1, 3.0, 0.0, 0.0, "H"),
        _pdb_atom_line(100, "ZN", "ZN", "A", 301, 0.0, 0.0, 0.0, "ZN", record="HETATM"),
    ]
    pdb = _write_pdb(tmp_path, "test4.pdb", lines)
    out = tmp_path / "out4"

    result = remove_metal_coordinating_hydrogens(
        input_pdb_path=pdb,
        output_dir=out,
        pdb_id="TEST",
        variant_label="test4",
    )

    assert result.status == "success"
    assert result.removed_count == 0
    assert result.kept_count >= 1

    # HG (serial 2) must be present
    cleaned_lines = result.cleaned_pdb_path.read_text().splitlines()
    serials = {
        int(line[6:11]) for line in cleaned_lines if line.startswith(("ATOM", "HETATM"))
    }
    assert 2 in serials


# ---------------------------------------------------------------------------
# Test 5: Non-standard ligand H facing metal → removed + manual_review.tsv
# ---------------------------------------------------------------------------

def test_ligand_h_removed_and_appears_in_manual_review(tmp_path: Path) -> None:
    """Ligand (non-standard) N-H pointing toward Zn is removed and flagged.

    The donor residue "LIG" is not in STANDARD_AMINO_ACID_NAMES; geometry
    rules still apply, but the decision also lands in manual_review.tsv.
    """
    lines = [
        _pdb_atom_line(1, "N1", "LIG", "A", 1, 2.0, 0.0, 0.0, "N", record="HETATM"),
        _pdb_atom_line(2, "HN1", "LIG", "A", 1, 1.2, 0.0, 0.0, "H", record="HETATM"),
        _pdb_atom_line(100, "ZN", "ZN", "A", 301, 0.0, 0.0, 0.0, "ZN", record="HETATM"),
    ]
    pdb = _write_pdb(tmp_path, "test5.pdb", lines)
    out = tmp_path / "out5"

    result = remove_metal_coordinating_hydrogens(
        input_pdb_path=pdb,
        output_dir=out,
        pdb_id="TEST",
        variant_label="test5",
    )

    assert result.status == "success"
    assert result.removed_count >= 1
    assert result.manual_review_count >= 1

    with result.manual_review_tsv_path.open(newline="") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    assert any(r["removed"] == "yes" for r in rows)
    # Verify the ligand residue is recorded
    assert any(r["donor_residue_name"] == "LIG" for r in rows)


# ---------------------------------------------------------------------------
# Test 6: metall_params creates expected per-site folder and scaffold files
# ---------------------------------------------------------------------------

def test_metall_params_creates_expected_site_structure(tmp_path: Path) -> None:
    """run_metal_parametrization_for_protein_dir creates a complete site folder.

    Checks:
      - site_001_ZN_A_301/ folder is created
      - 00_input/, 01_hydrogen_cleanup/, 02_contacts/, 03_mcpb/ sub-folders exist
      - MCPB scaffold files are present ({pdb_id}.in, commands.sh, README.md, {pdb_id}_mcpb.pdb)
      - {pdb_id}.in contains TODO placeholders and correct original_pdb line
      - commands.sh contains the MCPB.py step-1 invocation
      - manifest.json is written with correct pdb_id and site count
    """
    from stack_protein_preparation.metall_params import run_metal_parametrization_for_protein_dir

    pdb_id = "TZNX"
    protein_dir = tmp_path / pdb_id
    components_dir = protein_dir / "components"
    components_dir.mkdir(parents=True)

    # Protein: CYS SG coordinating Zn, HG pointing away (should be kept)
    protein_lines = [
        _pdb_atom_line(1, "SG", "CYS", "A", 1, 2.0, 0.0, 0.0, "S"),
        _pdb_atom_line(2, "HG", "CYS", "A", 1, 3.0, 0.0, 0.0, "H"),
    ]
    protein_pdb = components_dir / f"{pdb_id}_protein.pdb"
    protein_pdb.write_text("".join(protein_lines) + "END\n", encoding="utf-8")

    # Metal file: single Zn
    metal_lines = [
        _pdb_atom_line(100, "ZN", "ZN", "A", 301, 0.0, 0.0, 0.0, "ZN", record="HETATM"),
    ]
    metal_pdb = components_dir / f"{pdb_id}_metal.pdb"
    metal_pdb.write_text("".join(metal_lines) + "END\n", encoding="utf-8")

    result = run_metal_parametrization_for_protein_dir(
        protein_dir=protein_dir,
        chimera_executable=None,
    )

    assert result["status"] == "success", result.get("message", "")
    assert result["transition_metal_site_count"] == 1

    metall_params_dir = protein_dir / "metall_params"
    assert metall_params_dir.is_dir()

    # Locate site folder
    site_dirs = [d for d in metall_params_dir.iterdir() if d.name.startswith("site_001_ZN")]
    assert len(site_dirs) == 1, f"Expected one site_001_ZN_* folder, found: {[d.name for d in metall_params_dir.iterdir()]}"
    site_dir = site_dirs[0]

    # Sub-folders
    for sub in ("00_input", "01_hydrogen_cleanup", "02_contacts", "03_mcpb"):
        assert (site_dir / sub).is_dir(), f"Missing sub-folder: {sub}"

    # MCPB scaffold files
    mcpb_dir = site_dir / "03_mcpb"
    assert (mcpb_dir / f"{pdb_id}.in").is_file()
    assert (mcpb_dir / "commands.sh").is_file()
    assert (mcpb_dir / "README.md").is_file()
    assert (mcpb_dir / f"{pdb_id}_mcpb.pdb").is_file()

    # .in has TODO placeholders and correct original_pdb declaration
    mcpb_in_text = (mcpb_dir / f"{pdb_id}.in").read_text()
    assert "TODO" in mcpb_in_text
    assert f"original_pdb = {pdb_id}_mcpb.pdb" in mcpb_in_text

    # commands.sh has MCPB.py step-1 invocation
    cmd_text = (mcpb_dir / "commands.sh").read_text()
    assert "MCPB.py" in cmd_text
    assert f"{pdb_id}.in" in cmd_text

    # manifest.json
    manifest_path = metall_params_dir / "manifest.json"
    assert manifest_path.is_file()
    manifest = json.loads(manifest_path.read_text())
    assert manifest["pdb_id"] == pdb_id
    assert manifest["transition_metal_site_count"] == 1
    assert len(manifest["sites"]) == 1
