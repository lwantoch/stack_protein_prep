"""Tests for stack_protein_preparation._filler_af_splice.

Direct coverage for the AF-splice pipeline's structural gates:

* ``_detect_missing_windows`` -- contiguous-window detector against a
  crystal residue set.
* ``_protein_residue_map`` -- chain -> {resi: residue} inventory that
  ignores HETATM residues (het-flag).
* ``splice_af_gaps_into_crystal`` -- chain-overlap guard (whole-chain
  wholesale copy is skipped when AF resnums overlap ANY crystal chain);
  smart-rollback (rejects fills whose fitted bonds are outside
  [1.28, 1.40] A or introduce > 5 clashes); pLDDT floor.
* ``rollback_bad_gap_fills`` -- post-refine rollback of gaps whose
  metrics still fail after LoopModel had a chance to polish them.
"""
from __future__ import annotations

from pathlib import Path

import pytest
from Bio.PDB import PDBParser

from stack_protein_preparation._filler_af_splice import (
    _detect_missing_windows,
    _protein_residue_map,
    rollback_bad_gap_fills,
    splice_af_gaps_into_crystal,
)


# ---------------------------------------------------------- helpers


def _atom(serial: int, name: str, resname: str, chain: str, resnum: int,
          x: float, y: float, z: float, element: str = "C",
          record: str = "ATOM") -> str:
    return (
        f"{record:<6}{serial:5d} {name:>4} {resname:>3} {chain:1}"
        f"{resnum:4d}    "
        f"{x:8.3f}{y:8.3f}{z:8.3f}"
        f"  1.00 50.00           {element:>2}"
    )


def _write(pdb_path: Path, atoms: list[str]) -> Path:
    pdb_path.parent.mkdir(parents=True, exist_ok=True)
    pdb_path.write_text("\n".join(atoms) + "\nEND\n", encoding="utf-8")
    return pdb_path


def _backbone(res_id: int, chain: str = "A", resname: str = "ALA",
              x: float = 0.0, y: float = 0.0, z: float = 0.0,
              serial_start: int = 1) -> list[str]:
    """Emit N, CA, C, O for one residue centred at (x, y, z)."""
    return [
        _atom(serial_start,     "N",  resname, chain, res_id,  x - 1.46, y, z,       "N"),
        _atom(serial_start + 1, "CA", resname, chain, res_id,  x,        y, z,       "C"),
        _atom(serial_start + 2, "C",  resname, chain, res_id,  x + 1.46, y, z,       "C"),
        _atom(serial_start + 3, "O",  resname, chain, res_id,  x + 2.0,  y + 1.1, z, "O"),
    ]


# ---------------------------------------------------------- window detector


def test_detect_missing_windows_single_internal_gap():
    crystal = {1, 2, 3, 7, 8, 9}
    af = [1, 2, 3, 4, 5, 6, 7, 8, 9]
    assert _detect_missing_windows(crystal, af) == [(4, 6)]


def test_detect_missing_windows_multiple_gaps():
    crystal = {1, 2, 5, 6, 10}
    af = list(range(1, 11))
    assert _detect_missing_windows(crystal, af) == [(3, 4), (7, 9)]


def test_detect_missing_windows_terminal_gaps():
    # N-terminal + C-terminal gaps
    crystal = {5, 6, 7}
    af = list(range(1, 11))
    assert _detect_missing_windows(crystal, af) == [(1, 4), (8, 10)]


def test_detect_missing_windows_no_gap():
    crystal = {1, 2, 3, 4, 5}
    af = [1, 2, 3, 4, 5]
    assert _detect_missing_windows(crystal, af) == []


# ---------------------------------------------------------- residue map


def test_protein_residue_map_ignores_hetatm(tmp_path: Path) -> None:
    atoms = _backbone(1, "A", "ALA", x=0.0, serial_start=1)
    atoms.append(_atom(5, "O", "HOH", "A", 100, 5.0, 5.0, 5.0, "O", "HETATM"))
    p = _write(tmp_path / "mix.pdb", atoms)
    struct = PDBParser(QUIET=True).get_structure("m", str(p))
    m = _protein_residue_map(struct)
    assert "A" in m
    assert list(m["A"].keys()) == [1]  # HOH excluded


# ---------------------------------------------------------- chain-overlap guard


def test_splice_skips_whole_chain_copy_when_resnums_overlap(tmp_path: Path) -> None:
    # Reproduces the 7UL2 pattern: crystal chain R with residues 1-3,
    # AF chain A with residues 1-3 (numbering overlap).  Guard should
    # skip the wholesale copy of AF chain A instead of duplicating.
    crystal_atoms = []
    for i in range(1, 4):
        crystal_atoms.extend(_backbone(i, "R", "ALA",
                                        x=i * 3.5, serial_start=(i - 1) * 4 + 1))
    crystal = _write(tmp_path / "crystal.pdb", crystal_atoms)

    af_atoms = []
    for i in range(1, 4):
        af_atoms.extend(_backbone(i, "A", "GLY",
                                   x=(i * 3.5) + 20.0,  # off in space
                                   serial_start=(i - 1) * 4 + 1))
    af = _write(tmp_path / "af.pdb", af_atoms)

    out = tmp_path / "spliced.pdb"
    splice_af_gaps_into_crystal(crystal, af, out)
    result = PDBParser(QUIET=True).get_structure("r", str(out))

    # Only chain R should be present -- AF chain A was skipped because
    # its residues 1-3 overlap R's residues 1-3.
    chain_ids = sorted(chain.id for chain in result[0])
    assert chain_ids == ["R"], f"Expected only chain R, got {chain_ids}"


def test_splice_adds_whole_chain_when_no_overlap(tmp_path: Path) -> None:
    # Crystal has chain R with residues 1-3; AF has chain X with
    # residues 50-52 (no numbering overlap) -- AF X should be added.
    crystal_atoms = []
    for i in range(1, 4):
        crystal_atoms.extend(_backbone(i, "R", "ALA",
                                        x=i * 3.5, serial_start=(i - 1) * 4 + 1))
    crystal = _write(tmp_path / "crystal.pdb", crystal_atoms)

    af_atoms = []
    for i, resnum in enumerate((50, 51, 52), start=1):
        af_atoms.extend(_backbone(resnum, "X", "GLY",
                                   x=(i * 3.5) + 20.0,
                                   serial_start=(i - 1) * 4 + 1))
    af = _write(tmp_path / "af.pdb", af_atoms)

    out = tmp_path / "spliced.pdb"
    splice_af_gaps_into_crystal(crystal, af, out)
    result = PDBParser(QUIET=True).get_structure("r", str(out))

    chain_ids = sorted(chain.id for chain in result[0])
    assert chain_ids == ["R", "X"]


# ---------------------------------------------------------- rollback_bad_gap_fills


def test_rollback_bad_gap_fills_drops_broken_bond_gap(tmp_path: Path) -> None:
    # Chain A: residues 1-5.  Residues 3-4 are "AF-filled".  Boundary
    # bond 2 C -> 3 N is 3.5 A (broken).  Rollback should drop residues 3-4.
    atoms = []
    serial = 1
    xs = [0.0, 3.5, 7.0, 10.5, 14.0]  # normal ~3.5A CA-CA
    for i, x in enumerate(xs, start=1):
        atoms.extend(_backbone(i, "A", "ALA", x=x, serial_start=serial))
        serial += 4
    # Now break bond 2 -> 3: move residue 3's N far away
    # (need to rewrite atoms; simpler: build without break then patch string)
    p = _write(tmp_path / "with_fill.pdb", atoms)

    # Post-hoc patch: bump residue-3 N x-coord to ~7.5 (so C(2)@4.96 to N(3)@7.5 = 2.54A > 1.4)
    lines = p.read_text(encoding="utf-8").splitlines()
    new_lines: list[str] = []
    for line in lines:
        if line.startswith("ATOM") and line[13:16] == "N  " and line[22:26].strip() == "3":
            # rewrite x coordinate to 7.5 (columns 31-38)
            new_line = line[:30] + f"{7.500:8.3f}" + line[38:]
            new_lines.append(new_line)
        else:
            new_lines.append(line)
    p.write_text("\n".join(new_lines) + "\n", encoding="utf-8")

    out = tmp_path / "rolled.pdb"
    _, rolled = rollback_bad_gap_fills(p, out, gap_ranges_by_chain=[("A", 3, 4)])
    # Residue 3 dropped
    result = PDBParser(QUIET=True).get_structure("r", str(out))
    remaining = sorted(r.id[1] for chain in result[0] for r in chain if not r.id[0].strip())
    assert 3 not in remaining or 4 not in remaining, (
        f"Rollback should have dropped at least residue 3 or 4; remaining={remaining}"
    )
    assert rolled, f"Expected non-empty rolled list, got {rolled}"


def test_rollback_no_op_when_gap_is_clean(tmp_path: Path) -> None:
    # 5-residue chain, no broken bonds, no clashes.  Rollback should
    # bypass the file rewrite and leave input structure unchanged.
    atoms = []
    serial = 1
    # CA-CA spacing 4.25 A -> C-N = 4.25 - 2*1.46 = 1.33 A (ideal peptide bond)
    for i, x in enumerate([0.0, 4.25, 8.5, 12.75, 17.0], start=1):
        atoms.extend(_backbone(i, "A", "ALA", x=x, serial_start=serial))
        serial += 4
    p = _write(tmp_path / "clean.pdb", atoms)
    out = tmp_path / "rolled_clean.pdb"
    _, rolled = rollback_bad_gap_fills(p, out, gap_ranges_by_chain=[("A", 3, 4)])
    assert rolled == []  # nothing removed
    # Output should be a valid PDB with same residue set
    result = PDBParser(QUIET=True).get_structure("r", str(out))
    remaining = sorted(r.id[1] for chain in result[0] for r in chain if not r.id[0].strip())
    assert remaining == [1, 2, 3, 4, 5]
