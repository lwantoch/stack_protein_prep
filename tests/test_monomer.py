"""Unit tests for BioPython-based representative monomer extraction.

These tests use small synthetic PDB files because the monomer module should be
tested locally without depending on large structure fixtures. The tests focus on
the module's production API: chain-sequence analysis, conservative grouping, and
representative PDB write-out. They do not test biological assembly inference,
RMSD comparison, docking behavior, parametrization behavior, or MD preparation.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from stack_protein_preparation.monomer import (
    analyze_monomer_units,
    extract_chain_sequences,
    write_representative_monomer_units,
    write_single_representative_monomer_unit,
)


def test_analyze_monomer_units_groups_equivalent_chains_and_separates_partner(
    tmp_path: Path,
) -> None:
    """Equivalent copies should collapse into one group while partners stay separate."""

    input_pdb = tmp_path / "heteromer_with_copy.pdb"
    _write_pdb(
        input_pdb,
        [
            *_chain_records(chain_id="A", residue_names=["ALA", "GLY", "SER"]),
            *_chain_records(chain_id="B", residue_names=["ALA", "GLY", "SER"]),
            *_chain_records(chain_id="C", residue_names=["LYS", "TYR", "ASP"]),
        ],
    )

    result = analyze_monomer_units(input_pdb)

    assert [group.chain_ids for group in result.groups] == [("A", "B"), ("C",)]
    assert [group.representative_chain_id for group in result.groups] == ["A", "C"]
    assert not result.is_single_monomer_type

    comparison_ab = _find_comparison(result.comparisons, "A", "B")
    assert comparison_ab.sequence_identity == pytest.approx(1.0)
    assert comparison_ab.coverage_a == pytest.approx(1.0)
    assert comparison_ab.coverage_b == pytest.approx(1.0)


def test_analyze_monomer_units_does_not_group_heavily_truncated_chain(
    tmp_path: Path,
) -> None:
    """High identity over a short overlap should not be enough for equivalence."""

    input_pdb = tmp_path / "truncated_copy.pdb"
    _write_pdb(
        input_pdb,
        [
            *_chain_records(chain_id="A", residue_names=["ALA", "GLY", "SER", "THR"]),
            *_chain_records(chain_id="B", residue_names=["ALA", "GLY"]),
        ],
    )

    result = analyze_monomer_units(
        input_pdb,
        identity_threshold=1.0,
        coverage_threshold=0.90,
    )

    assert [group.chain_ids for group in result.groups] == [("A",), ("B",)]

    comparison_ab = _find_comparison(result.comparisons, "A", "B")
    assert comparison_ab.sequence_identity == pytest.approx(1.0)
    assert comparison_ab.coverage_a == pytest.approx(0.5)
    assert comparison_ab.coverage_b == pytest.approx(1.0)


def test_write_representative_monomer_units_writes_one_chain_per_group(
    tmp_path: Path,
) -> None:
    """Representative write-out should not write every equivalent chain copy."""

    input_pdb = tmp_path / "heteromer_with_copy.pdb"
    output_dir = tmp_path / "monomers"

    _write_pdb(
        input_pdb,
        [
            *_chain_records(chain_id="A", residue_names=["ALA", "GLY", "SER"]),
            *_chain_records(chain_id="B", residue_names=["ALA", "GLY", "SER"]),
            *_chain_records(chain_id="C", residue_names=["LYS", "TYR", "ASP"]),
        ],
    )

    results = write_representative_monomer_units(input_pdb, output_dir)

    assert [result.representative_chain_id for result in results] == ["A", "C"]
    assert [result.represented_chain_ids for result in results] == [("A", "B"), ("C",)]
    assert sorted(path.name for path in output_dir.iterdir()) == [
        "heteromer_with_copy_monomer_chain_A.pdb",
        "heteromer_with_copy_monomer_chain_C.pdb",
    ]

    chain_a_output = output_dir / "heteromer_with_copy_monomer_chain_A.pdb"
    chain_c_output = output_dir / "heteromer_with_copy_monomer_chain_C.pdb"

    assert _sequence_by_chain(chain_a_output) == {"A": "AGS"}
    assert _sequence_by_chain(chain_c_output) == {"C": "KYD"}

    assert all(result.coordinate_records_written > 0 for result in results)


def test_write_single_representative_monomer_unit_raises_for_heteromer(
    tmp_path: Path,
) -> None:
    """Strict single-output mode should fail when multiple chain groups exist."""

    input_pdb = tmp_path / "heteromer.pdb"
    output_pdb = tmp_path / "single_monomer.pdb"

    _write_pdb(
        input_pdb,
        [
            *_chain_records(chain_id="A", residue_names=["ALA", "GLY", "SER"]),
            *_chain_records(chain_id="B", residue_names=["LYS", "TYR", "ASP"]),
        ],
    )

    with pytest.raises(ValueError, match="multiple non-equivalent chain groups"):
        write_single_representative_monomer_unit(input_pdb, output_pdb)

    assert not output_pdb.exists()


def test_write_single_representative_monomer_unit_writes_homomer_representative(
    tmp_path: Path,
) -> None:
    """Strict single-output mode should write one representative for a homomer."""

    input_pdb = tmp_path / "homomer.pdb"
    output_pdb = tmp_path / "single_monomer.pdb"

    _write_pdb(
        input_pdb,
        [
            *_chain_records(chain_id="A", residue_names=["ALA", "GLY", "SER"]),
            *_chain_records(chain_id="B", residue_names=["ALA", "GLY", "SER"]),
        ],
    )

    result = write_single_representative_monomer_unit(input_pdb, output_pdb)

    assert result.representative_chain_id == "A"
    assert result.represented_chain_ids == ("A", "B")
    assert output_pdb.exists()
    assert _sequence_by_chain(output_pdb) == {"A": "AGS"}


def test_write_chain_keeps_modified_residue_but_not_generic_hetero_by_default(
    tmp_path: Path,
) -> None:
    """Default write-out should retain modified amino acids but drop ligand/water."""

    input_pdb = tmp_path / "modified_residue.pdb"
    output_pdb = tmp_path / "modified_residue_monomer.pdb"

    _write_pdb(
        input_pdb,
        [
            _pdb_record(
                record_name="ATOM",
                serial=1,
                atom_name="CA",
                residue_name="ALA",
                chain_id="A",
                residue_number=1,
                element="C",
            ),
            _pdb_record(
                record_name="HETATM",
                serial=2,
                atom_name="SE",
                residue_name="MSE",
                chain_id="A",
                residue_number=2,
                element="SE",
            ),
            _pdb_record(
                record_name="HETATM",
                serial=3,
                atom_name="C1",
                residue_name="LIG",
                chain_id="A",
                residue_number=100,
                element="C",
            ),
            _pdb_record(
                record_name="HETATM",
                serial=4,
                atom_name="O",
                residue_name="HOH",
                chain_id="A",
                residue_number=200,
                element="O",
            ),
        ],
    )

    result = write_single_representative_monomer_unit(input_pdb, output_pdb)
    output_text = output_pdb.read_text(encoding="utf-8")

    assert result.coordinate_records_written == 2
    assert "ALA" in output_text
    assert "MSE" in output_text
    assert "LIG" not in output_text
    assert "HOH" not in output_text


def _sequence_by_chain(input_pdb: Path) -> dict[str, str]:
    """Return extracted test sequences keyed by chain ID."""

    return {
        chain.chain_id: chain.sequence
        for chain in extract_chain_sequences(input_pdb)
    }


def _find_comparison(comparisons, chain_id_a: str, chain_id_b: str):
    """Return one pairwise comparison independent of stored chain order."""

    requested = {chain_id_a, chain_id_b}

    for comparison in comparisons:
        observed = {comparison.chain_id_a, comparison.chain_id_b}
        if observed == requested:
            return comparison

    raise AssertionError(f"No comparison found for chains {chain_id_a!r}, {chain_id_b!r}.")


def _write_pdb(output_pdb: Path, records: list[str]) -> None:
    """Write a minimal synthetic PDB file for parser-based tests."""

    output_pdb.write_text("\n".join([*records, "END"]) + "\n", encoding="utf-8")


def _chain_records(chain_id: str, residue_names: list[str]) -> list[str]:
    """Return one CA atom per residue for a synthetic protein chain."""

    records: list[str] = []

    for index, residue_name in enumerate(residue_names, start=1):
        records.append(
            _pdb_record(
                record_name="ATOM",
                serial=_serial_from_chain_and_residue(chain_id, index),
                atom_name="CA",
                residue_name=residue_name,
                chain_id=chain_id,
                residue_number=index,
                element="C",
            )
        )

    return records


def _serial_from_chain_and_residue(chain_id: str, residue_number: int) -> int:
    """Return deterministic atom serials for synthetic PDB records."""

    return (ord(chain_id) - ord("A") + 1) * 100 + residue_number


def _pdb_record(
    *,
    record_name: str,
    serial: int,
    atom_name: str,
    residue_name: str,
    chain_id: str,
    residue_number: int,
    element: str,
) -> str:
    """Return one fixed-column PDB coordinate record.

    The production module deliberately uses BioPython instead of local column
    parsing, but the tests still need compact valid PDB snippets. This helper
    writes the small subset of fixed-column fields that BioPython's PDBParser
    needs for atoms, residues, chains, coordinates, occupancy, B factor, and
    element names. Coordinates are deterministic but chemically irrelevant for
    these sequence/grouping tests.
    """

    x_coord = float(serial)
    y_coord = float(residue_number)
    z_coord = 0.0

    return (
        f"{record_name:<6}{serial:5d} {atom_name:^4s} "
        f"{residue_name:>3s} {chain_id:1s}{residue_number:4d}    "
        f"{x_coord:8.3f}{y_coord:8.3f}{z_coord:8.3f}"
        f"{1.00:6.2f}{20.00:6.2f}          {element:>2s}"
    )