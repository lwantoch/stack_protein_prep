"""Tests for the new geometric-break detection added to gaps.summarize_gaps.

Drop this file at ``tests/test_gaps_geometric_break.py`` in the FRUTON repo
after applying ``fruton_fix_geometric_gap_detection.patch``.
"""
from __future__ import annotations

import tempfile
from pathlib import Path

from stack_protein_preparation.gaps import (
    geometric_breaks_by_chain_for_pdb,
    summarize_gaps,
)


def _write_pdb(text: str) -> Path:
    tmp = tempfile.NamedTemporaryFile(mode="w", suffix=".pdb", delete=False)
    tmp.write(text)
    tmp.close()
    return Path(tmp.name)


# Two ALA residues on chain A with contiguous numbering (1 -> 2) but the C
# atom of residue 1 is at (0,0,0) and the N atom of residue 2 is at (0,0,12)
# -- a 12 A physical gap that FRUTON's existing detector misses (no numbering
# jump). This is the 6Z1T-style undeclared break.
_SYNTHETIC_CONTIGUOUS_GAP = """\
ATOM      1  N   ALA A   1       0.000   0.000   1.500  1.00  0.00           N
ATOM      2  CA  ALA A   1       1.500   0.000   1.500  1.00  0.00           C
ATOM      3  C   ALA A   1       0.000   0.000   0.000  1.00  0.00           C
ATOM      4  O   ALA A   1       0.000   1.200   0.000  1.00  0.00           O
ATOM      5  N   ALA A   2       0.000   0.000  12.000  1.00  0.00           N
ATOM      6  CA  ALA A   2       1.500   0.000  12.000  1.00  0.00           C
ATOM      7  C   ALA A   2       0.000   0.000  13.500  1.00  0.00           C
ATOM      8  O   ALA A   2       0.000   1.200  13.500  1.00  0.00           O
TER       9      ALA A   2
END
"""


def test_geometric_break_detected_with_contiguous_numbering():
    pdb_path = _write_pdb(_SYNTHETIC_CONTIGUOUS_GAP)
    try:
        summary = summarize_gaps(pdb_path)
        # Existing detector: no numbering jump -> no traditional gap.
        assert summary["n_gaps"] == 0
        assert summary["has_gaps"] is False
        # New detector: 12 A C-N with contiguous numbering -> one geometric break.
        assert summary["n_geometric_breaks"] == 1
        breaks = summary["geometric_breaks"]
        assert len(breaks) == 1
        b = breaks[0]
        assert b["chain"] == "A"
        assert b["prev_res"] == 1
        assert b["next_res"] == 2
        assert b["peptide_cn_distance_angstrom"] > 2.0
    finally:
        pdb_path.unlink()


def test_geometric_break_ignored_when_c_n_normal():
    """Backbone-continuous synthetic -- no geometric break should appear."""
    # Peptide bond C-N ~ 1.33 A (not 3.8 A which is CA-CA in beta strand).
    # Place N(res2) at z=1.33 (1.33 A from C(res1) at z=0), CA at same z=1.33
    # (test template mirrors N and CA z), C at z=2.83 (1.5 A from CA), O
    # follows C.
    normal_pdb = _SYNTHETIC_CONTIGUOUS_GAP.replace(
        "  12.000  1.00  0.00           N",
        "   1.330  1.00  0.00           N",
    ).replace(
        "  12.000  1.00  0.00           C",
        "   1.330  1.00  0.00           C",
    ).replace(
        "  13.500  1.00  0.00           C",
        "   2.830  1.00  0.00           C",
    ).replace(
        "  13.500  1.00  0.00           O",
        "   2.830  1.00  0.00           O",
    )
    pdb_path = _write_pdb(normal_pdb)
    try:
        summary = summarize_gaps(pdb_path)
        assert summary["n_geometric_breaks"] == 0
    finally:
        pdb_path.unlink()


def test_helper_returns_per_chain_grouping():
    pdb_path = _write_pdb(_SYNTHETIC_CONTIGUOUS_GAP)
    try:
        grouped = geometric_breaks_by_chain_for_pdb(pdb_path)
        assert "A" in grouped
        assert len(grouped["A"]) == 1
    finally:
        pdb_path.unlink()
