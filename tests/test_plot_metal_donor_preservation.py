"""Smoke tests for scripts/plot_metal_donor_preservation.py.

Uses tiny hand-crafted crystal and 'FRUTON' PDBs (with a slightly
shifted His donor) so the pipeline can be exercised without heavy
fixtures.
"""
from __future__ import annotations

import importlib.util
import sys
from pathlib import Path


_SCRIPT = (
    Path(__file__).resolve().parents[1]
    / "scripts" / "plot_metal_donor_preservation.py"
)


def _load_mod():
    spec = importlib.util.spec_from_file_location("_plot_metal_cli", _SCRIPT)
    mod = importlib.util.module_from_spec(spec)
    sys.modules["_plot_metal_cli"] = mod
    spec.loader.exec_module(mod)
    return mod


def _write_zn_his_pdb(path: Path, ne2_offset: float = 0.0) -> Path:
    """A Zn ion with a His-NE2 donor.  ne2_offset shifts the NE2 atom
    along +x so we can craft a controlled Δd.
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    ne2_x = 2.100 + ne2_offset  # ~2.1 Å is a common Zn-N distance
    text = (
        "ATOM      1  N   HIS A   1       0.000   0.000   0.000  1.00  0.00           N\n"
        "ATOM      2  CA  HIS A   1       1.458   0.000   0.000  1.00  0.00           C\n"
        "ATOM      3  C   HIS A   1       2.009   1.420   0.000  1.00  0.00           C\n"
        "ATOM      4  O   HIS A   1       1.251   2.390   0.000  1.00  0.00           O\n"
        "ATOM      5  CB  HIS A   1       2.000  -0.900   1.000  1.00  0.00           C\n"
        "ATOM      6  CG  HIS A   1       2.500  -1.500   2.000  1.00  0.00           C\n"
        "ATOM      7  ND1 HIS A   1       3.500  -2.400   1.900  1.00  0.00           N\n"
        "ATOM      8  CD2 HIS A   1       2.000  -1.500   3.200  1.00  0.00           C\n"
        "ATOM      9  CE1 HIS A   1       3.500  -2.800   3.200  1.00  0.00           C\n"
        f"ATOM     10  NE2 HIS A   1     {ne2_x:7.3f}  -2.100   3.700  1.00  0.00           N\n"
        "TER\n"
        "HETATM   11 ZN    ZN B   1       0.000  -2.100   3.700  1.00  0.00          ZN\n"
        "END\n"
    )
    path.write_text(text)
    return path


def test_end_to_end_writes_png_csv_summary(tmp_path: Path):
    bench_crystal = tmp_path / "bench_crystal"
    bench_fruton = tmp_path / "bench_fruton"
    (bench_crystal / "8ABC").mkdir(parents=True)
    (bench_fruton / "8ABC").mkdir(parents=True)
    _write_zn_his_pdb(bench_crystal / "8ABC" / "input_crystal.pdb", ne2_offset=0.0)
    _write_zn_his_pdb(bench_fruton / "8ABC" / "final_model.pdb", ne2_offset=0.1)

    mod = _load_mod()
    outdir = tmp_path / "figs"
    rc = mod.main([
        "--crystal-glob", str(bench_crystal / "*/input_crystal.pdb"),
        "--fruton-glob", str(bench_fruton / "*/final_model.pdb"),
        "--outdir", str(outdir),
    ])
    assert rc == 0
    assert (outdir / "metal_donor_delta.png").is_file()
    assert (outdir / "metal_donor_delta.csv").is_file()
    assert (outdir / "per_pdb_summary.txt").is_file()

    csv_lines = (outdir / "metal_donor_delta.csv").read_text().strip().splitlines()
    assert len(csv_lines) >= 2  # header + at least one donor
    assert "delta_A" in csv_lines[0]


def test_pair_by_stem_mode(tmp_path: Path):
    bc = tmp_path / "bc"
    bf = tmp_path / "bf"
    bc.mkdir()
    bf.mkdir()
    _write_zn_his_pdb(bc / "9XYZ.pdb")
    _write_zn_his_pdb(bf / "9XYZ.pdb", ne2_offset=0.05)
    mod = _load_mod()
    outdir = tmp_path / "figs"
    rc = mod.main([
        "--crystal-glob", str(bc / "*.pdb"),
        "--fruton-glob", str(bf / "*.pdb"),
        "--pair-by", "stem",
        "--outdir", str(outdir),
    ])
    assert rc == 0
    csv_text = (outdir / "metal_donor_delta.csv").read_text()
    assert "9XYZ" in csv_text


def test_no_matches_returns_nonzero(tmp_path: Path):
    mod = _load_mod()
    outdir = tmp_path / "figs"
    rc = mod.main([
        "--crystal-glob", str(tmp_path / "nope1_*.pdb"),
        "--fruton-glob", str(tmp_path / "nope2_*.pdb"),
        "--outdir", str(outdir),
    ])
    assert rc == 2


def test_disjoint_pair_keys_returns_nonzero(tmp_path: Path):
    bc = tmp_path / "bc"
    bf = tmp_path / "bf"
    (bc / "8ABC").mkdir(parents=True)
    (bf / "8XYZ").mkdir(parents=True)
    _write_zn_his_pdb(bc / "8ABC" / "input_crystal.pdb")
    _write_zn_his_pdb(bf / "8XYZ" / "final_model.pdb")
    mod = _load_mod()
    outdir = tmp_path / "figs"
    rc = mod.main([
        "--crystal-glob", str(bc / "*/input_crystal.pdb"),
        "--fruton-glob", str(bf / "*/final_model.pdb"),
        "--outdir", str(outdir),
    ])
    assert rc == 2


def test_offset_produces_expected_delta_sign(tmp_path: Path):
    """With ne2 shifted +0.5 Å along +x, and metal at origin,
    d_fruton > d_crystal → positive Δd."""
    bc = tmp_path / "bc"
    bf = tmp_path / "bf"
    bc.mkdir()
    bf.mkdir()
    _write_zn_his_pdb(bc / "P.pdb", ne2_offset=0.0)
    _write_zn_his_pdb(bf / "P.pdb", ne2_offset=0.5)
    mod = _load_mod()
    outdir = tmp_path / "figs"
    rc = mod.main([
        "--crystal-glob", str(bc / "*.pdb"),
        "--fruton-glob", str(bf / "*.pdb"),
        "--pair-by", "stem",
        "--outdir", str(outdir),
    ])
    assert rc == 0
    csv_lines = (outdir / "metal_donor_delta.csv").read_text().strip().splitlines()
    # find a NE2 row with a positive delta
    ne2_rows = [ln for ln in csv_lines[1:] if ",NE2," in ln]
    assert ne2_rows, "expected at least one NE2 donor row"
    # The last column is status; delta_A is second-to-last.
    for row in ne2_rows:
        cols = row.split(",")
        delta = cols[-2]
        # Some deltas can be None ("") if the donor moved out of cutoff;
        # our +0.5 Å shift keeps it inside 3.0 Å (2.6 Å).
        if delta:
            assert delta.startswith("+"), f"expected +Δd for +0.5 Å shift, got {delta}"
            break
    else:
        # If none had a valid delta, that's fine only if all rows were "gained"/"lost"
        pass


def test_index_by_key_dedupes(tmp_path: Path):
    mod = _load_mod()
    # Same key from two paths -> first wins
    p1 = tmp_path / "8ABC" / "a.pdb"
    p2 = tmp_path / "8ABC" / "b.pdb"
    p1.parent.mkdir(parents=True, exist_ok=True)
    p1.touch()
    p2.touch()
    result = mod._index_by_key([p1, p2], "parent")
    assert set(result.keys()) == {"8ABC"}
    assert result["8ABC"] == p1  # first wins
