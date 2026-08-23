"""Tests for stack_protein_preparation._bench_split.

Encodes reviewer-critical invariants:
- train / val / holdout are disjoint (no data leakage)
- packaged holdout list has exactly 30 well-formed PDB IDs
- loader tolerates both CSV format + directory-glob fallback
"""
from __future__ import annotations

from pathlib import Path

import pytest

from stack_protein_preparation import _bench_split as bs


# ---------------------------------------------------------------------------
# packaged holdout list
# ---------------------------------------------------------------------------


def test_holdout_list_has_exactly_30_ids():
    ids = bs.load_holdout_ids()
    assert len(ids) == 30
    assert len(set(ids)) == 30  # uniqueness


def test_holdout_ids_are_well_formed():
    for pid in bs.load_holdout_ids():
        assert bs.PDB_ID_RE.match(pid), f"malformed holdout PDB id: {pid!r}"


def test_holdout_disjoint_from_mmbsa_200():
    """The holdout list must not overlap the deployed MMBSA_200 CSV — that
    would defeat the whole point of a blind test set."""
    import csv
    from pathlib import Path
    csv_path = Path("/mnt/lustre/scratch/nlsas/home/otras/hcx/lwa/MMBSA_200/mmbsa200_targets.csv")
    if not csv_path.is_file():
        pytest.skip("MMBSA_200 targets CSV not available on this host")
    with csv_path.open() as fh:
        reader = csv.DictReader(fh)
        m200 = {row.get("best_pdb", "").upper() for row in reader}
    m200.discard("")
    holdout = set(bs.load_holdout_ids())
    overlap = m200 & holdout
    assert not overlap, f"holdout leaked into MMBSA_200: {sorted(overlap)}"


# ---------------------------------------------------------------------------
# split loader (dir + CSV)
# ---------------------------------------------------------------------------


def test_load_ids_from_dir_with_pdb_ids_csv(tmp_path: Path):
    root = tmp_path / "MMBSA_test"
    root.mkdir()
    (root / "pdb_ids.csv").write_text("pdb_id,range\n8ABC,\n8XYZ,\n9RCI,\n")
    got = bs._load_ids_from_dir(root)
    assert got == ["8ABC", "8XYZ", "9RCI"]


def test_load_ids_from_dir_fallback_to_glob(tmp_path: Path):
    root = tmp_path / "MMBSA_test"
    root.mkdir()
    (root / "8ABC").mkdir()
    (root / "8XYZ").mkdir()
    (root / "logs").mkdir()  # non-PDB dir should be skipped
    got = bs._load_ids_from_dir(root)
    assert got == ["8ABC", "8XYZ"]


def test_load_ids_from_dir_skips_non_pdb_ids_in_csv(tmp_path: Path):
    root = tmp_path / "MMBSA_test"
    root.mkdir()
    (root / "pdb_ids.csv").write_text(
        "pdb_id,range\n8ABC,\ninvalid_id,\n8XYZ,\n"
    )
    got = bs._load_ids_from_dir(root)
    assert got == ["8ABC", "8XYZ"]


# ---------------------------------------------------------------------------
# disjoint-splits invariant (reviewer-critical)
# ---------------------------------------------------------------------------


def test_assert_disjoint_passes_when_no_overlap():
    bs.assert_splits_disjoint(["A", "B"], ["C", "D"], ["E", "F"])  # no raise


def test_assert_disjoint_flags_train_val_overlap():
    with pytest.raises(ValueError, match="train ∩ val"):
        bs.assert_splits_disjoint(["A", "B"], ["B", "C"], ["D"])


def test_assert_disjoint_flags_train_holdout_overlap():
    with pytest.raises(ValueError, match="train ∩ holdout"):
        bs.assert_splits_disjoint(["A", "B"], ["C"], ["A", "D"])


def test_assert_disjoint_flags_val_holdout_overlap():
    with pytest.raises(ValueError, match="val ∩ holdout"):
        bs.assert_splits_disjoint(["A"], ["B", "C"], ["C", "D"])


# ---------------------------------------------------------------------------
# real LUSTRE data (skipped when unavailable, e.g. CI)
# ---------------------------------------------------------------------------


def test_train_val_holdout_all_disjoint_on_real_lustre():
    import os
    lustre = os.environ.get("LUSTRE")
    if not lustre or not Path(lustre, "MMBSA_200", "MMBSA_75").is_dir():
        pytest.skip("LUSTRE not mounted or MMBSA_75/MMBSA_125 absent")
    train = bs.load_split("train")
    val = bs.load_split("val")
    holdout = bs.load_split("holdout")
    bs.assert_splits_disjoint(train, val, holdout)


def test_train_and_val_have_expected_sizes_on_real_lustre():
    import os
    lustre = os.environ.get("LUSTRE")
    if not lustre or not Path(lustre, "MMBSA_200", "MMBSA_75").is_dir():
        pytest.skip("LUSTRE not mounted")
    train = bs.load_split("train")
    val = bs.load_split("val")
    # From the actual pdb_ids.csv files: 75 + 124 = 199 (matches baseline QC)
    assert 70 <= len(train) <= 80, f"train set size unexpected: {len(train)}"
    assert 120 <= len(val) <= 130, f"val set size unexpected: {len(val)}"


def test_load_split_rejects_unknown_name():
    with pytest.raises(ValueError, match="unknown split"):
        bs.load_split("foo")
