"""Tests for pdb_components — residue sets, line parsers, chain selection."""
from __future__ import annotations

from pathlib import Path

import pytest

from stack_protein_preparation.pdb_components import (
    COFACTORS,
    CRYSTALLIZATION_ARTIFACTS,
    METALS,
    PROTEIN_LIKE_NONSTANDARD,
    STANDARD_AMINO_ACIDS,
    WATER_NAMES,
    _chain_id,
    _element,
    _icode,
    _resname,
    _resseq,
    choose_relevant_chain_ids_from_range,
    split_pdb_components,
)


# ---------------------------------------------------------------------------
# Residue-set membership
# ---------------------------------------------------------------------------

class TestStandardAminoAcids:
    def test_all_twenty_present(self) -> None:
        canonical = {
            "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY",
            "HIS", "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER",
            "THR", "TRP", "TYR", "VAL",
        }
        assert canonical == STANDARD_AMINO_ACIDS

    def test_does_not_contain_mse(self) -> None:
        assert "MSE" not in STANDARD_AMINO_ACIDS


class TestProteinLikeNonstandard:
    def test_contains_phosphotyrosine(self) -> None:
        assert "PTR" in PROTEIN_LIKE_NONSTANDARD

    def test_contains_cysteic_acid(self) -> None:
        assert "CSD" in PROTEIN_LIKE_NONSTANDARD

    def test_contains_selenomethionine(self) -> None:
        assert "MSE" in PROTEIN_LIKE_NONSTANDARD

    def test_contains_phosphoserine(self) -> None:
        assert "SEP" in PROTEIN_LIKE_NONSTANDARD

    def test_contains_phosphothreonine(self) -> None:
        assert "TPO" in PROTEIN_LIKE_NONSTANDARD


class TestCrystallizationArtifacts:
    def test_sulfate_is_artifact(self) -> None:
        assert "SO4" in CRYSTALLIZATION_ARTIFACTS

    def test_no_overlap_with_metals(self) -> None:
        assert METALS.isdisjoint(CRYSTALLIZATION_ARTIFACTS)

    def test_no_overlap_with_standard_aa(self) -> None:
        assert STANDARD_AMINO_ACIDS.isdisjoint(CRYSTALLIZATION_ARTIFACTS)


class TestCofactors:
    def test_nad_is_cofactor(self) -> None:
        assert "NAD" in COFACTORS

    def test_heme_is_cofactor(self) -> None:
        assert "HEM" in COFACTORS

    def test_atp_is_cofactor(self) -> None:
        assert "ATP" in COFACTORS

    def test_cofactors_not_artifacts(self) -> None:
        assert COFACTORS.isdisjoint(CRYSTALLIZATION_ARTIFACTS)


# ---------------------------------------------------------------------------
# PDB line field parsers
# ---------------------------------------------------------------------------

def _atom_line(
    resname: str = "ALA",
    chain: str = "A",
    resseq: int = 1,
    icode: str = " ",
    element: str = "C",
    record: str = "ATOM",
) -> str:
    return (
        f"{record:<6}{1:5d}  CA  {resname:3s} {chain}{resseq:4d}{icode}   "
        f"   1.000   2.000   3.000  1.00  0.00          {element:>2s}  \n"
    )


def test_resname_extraction() -> None:
    assert _resname(_atom_line(resname="GLY")) == "GLY"


def test_chain_id_extraction() -> None:
    assert _chain_id(_atom_line(chain="B")) == "B"


def test_resseq_extraction() -> None:
    assert _resseq(_atom_line(resseq=42)) == "42"


def test_icode_extraction_blank() -> None:
    assert _icode(_atom_line(icode=" ")) == ""


def test_icode_extraction_letter() -> None:
    assert _icode(_atom_line(icode="A")) == "A"


def test_element_extraction() -> None:
    assert _element(_atom_line(element="N")) == "N"


# ---------------------------------------------------------------------------
# choose_relevant_chain_ids_from_range
# ---------------------------------------------------------------------------

def _pdb_lines_with_chains(chain_residues: dict[str, range]) -> list[str]:
    lines = []
    serial = 1
    for chain, rng in chain_residues.items():
        for resseq in rng:
            lines.append(
                f"ATOM  {serial:5d}  CA  ALA {chain}{resseq:4d}    "
                f"   1.000   2.000   3.000  1.00  0.00           C  \n"
            )
            serial += 1
    return lines


def test_single_chain_chosen_for_range() -> None:
    lines = _pdb_lines_with_chains({"A": range(1, 51)})
    result = choose_relevant_chain_ids_from_range(lines, "10-40")
    assert result == ("A",)


def test_best_overlap_chain_chosen() -> None:
    lines = _pdb_lines_with_chains({
        "A": range(1, 101),
        "B": range(200, 301),
    })
    result = choose_relevant_chain_ids_from_range(lines, "1-100")
    assert result == ("A",)


def test_no_overlap_returns_all_chains() -> None:
    lines = _pdb_lines_with_chains({
        "A": range(1, 51),
        "B": range(100, 151),
    })
    result = choose_relevant_chain_ids_from_range(lines, "500-600")
    assert set(result) == {"A", "B"}


def test_chain_qualified_range_respected() -> None:
    lines = _pdb_lines_with_chains({
        "A": range(1, 51),
        "B": range(1, 51),
    })
    result = choose_relevant_chain_ids_from_range(lines, "A1-A50")
    assert result == ("A",)


# ---------------------------------------------------------------------------
# split_pdb_components — smoke test
# ---------------------------------------------------------------------------

def _write_pdb(tmp_path: Path, name: str, lines: list[str]) -> Path:
    p = tmp_path / name
    p.write_text("".join(lines) + "END\n", encoding="utf-8")
    return p


def test_split_writes_protein_file(tmp_path: Path) -> None:
    lines = [
        f"ATOM  {i:5d}  CA  ALA A{i:4d}       1.000   2.000   3.000  1.00  0.00           C  \n"
        for i in range(1, 6)
    ]
    pdb = _write_pdb(tmp_path, "TEST.pdb", lines)
    out_dir = tmp_path / "components"
    split_pdb_components(pdb, out_dir, protein_stem="TEST")
    protein_file = out_dir / "TEST_protein.pdb"
    assert protein_file.is_file()
    content = protein_file.read_text()
    assert "ALA" in content


def test_split_writes_water_file(tmp_path: Path) -> None:
    lines = [
        f"ATOM  {1:5d}  CA  ALA A{1:4d}       1.000   2.000   3.000  1.00  0.00           C  \n",
        f"HETATM{2:5d}  O   HOH A{2:4d}       4.000   5.000   6.000  1.00  0.00           O  \n",
    ]
    pdb = _write_pdb(tmp_path, "TEST.pdb", lines)
    out_dir = tmp_path / "components"
    split_pdb_components(pdb, out_dir, protein_stem="TEST")
    assert (out_dir / "TEST_water.pdb").is_file()
