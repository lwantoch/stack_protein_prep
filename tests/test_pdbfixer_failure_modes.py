"""Tests for the pdbfixer failure-mode catalogue.

Chief purpose: guarantee that the machine-readable catalogue stays in
sync with docs/methods_no_pdbfixer.md.  Each of the five documented
failure modes must appear in both places.
"""
from __future__ import annotations

from pathlib import Path

import pytest

from stack_protein_preparation import _pdbfixer_failure_modes as pfm


_DOCS_PATH = (
    Path(__file__).resolve().parents[1]
    / "docs" / "methods_no_pdbfixer.md"
)


EXPECTED_MODE_IDS = (
    "protonation-generic",
    "metal-his-tautomer",
    "cofactor-deletion",
    "topology-bond-guessing",
    "gap-fill-no-template",
)


# ---------------------------------------------------------------------------
# Catalogue structure
# ---------------------------------------------------------------------------


def test_five_documented_modes():
    assert len(pfm.FAILURE_MODES) == 5


def test_expected_mode_ids_present():
    assert pfm.all_mode_ids() == EXPECTED_MODE_IDS


def test_each_mode_has_all_fields_nonempty():
    for m in pfm.FAILURE_MODES:
        assert m.mode_id and len(m.mode_id) > 5
        assert m.title and len(m.title) > 10
        assert m.chemistry and len(m.chemistry) > 40
        assert m.consequence_for_mmbsa and len(m.consequence_for_mmbsa) > 30
        assert m.fruton_alternative and len(m.fruton_alternative) > 40
        assert m.fruton_module.startswith("stack_protein_preparation.")


def test_get_failure_mode_lookup():
    m = pfm.get_failure_mode("metal-his-tautomer")
    assert m.mode_id == "metal-his-tautomer"
    assert "Zn" in m.chemistry or "metal" in m.chemistry.lower()


def test_get_failure_mode_unknown_raises():
    with pytest.raises(KeyError, match="unknown pdbfixer failure mode"):
        pfm.get_failure_mode("not-a-real-mode")


def test_format_methods_paragraph_contains_key_terms():
    p = pfm.format_methods_paragraph()
    assert "FRUTON" in p
    assert "pdbfixer" in p
    assert "MM-GBSA" in p
    assert "methods_no_pdbfixer.md" in p


# ---------------------------------------------------------------------------
# Docs synchronisation
# ---------------------------------------------------------------------------


def test_docs_file_exists():
    assert _DOCS_PATH.is_file(), (
        f"Expected reviewer-facing markdown at {_DOCS_PATH}; "
        "if you moved it, update tests/test_pdbfixer_failure_modes.py too."
    )


def test_docs_mentions_all_failure_modes_by_topic():
    """Each catalogued mode's title concept should appear in the markdown."""
    text = _DOCS_PATH.read_text().lower()
    for keyword in (
        "protonation",   # protonation-generic
        "metal",         # metal-his-tautomer
        "cofactor",      # cofactor-deletion
        "topology",      # topology-bond-guessing
        "gap-fill",      # gap-fill-no-template  (dashed or "gap fill")
    ):
        assert keyword in text, f"docs missing keyword {keyword!r}"


def test_docs_states_hard_rule_date():
    text = _DOCS_PATH.read_text()
    assert "hard rule" in text.lower()
    assert "2026-08-22" in text  # hard-rule enactment date


def test_docs_names_fruton_alternatives():
    """The doc must cite the FRUTON module that supersedes pdbfixer per mode."""
    text = _DOCS_PATH.read_text()
    for expected in (
        "_paper_override_suggest",
        "metal_hydrogen_cleanup",
        "_cofactor_params",
        "_filler_af_splice",
        "_filler_alphafold",
    ):
        assert expected in text, f"docs missing FRUTON module ref {expected!r}"


def test_docs_ends_with_fair_use_section():
    """Reviewer-honesty check: docs must acknowledge legitimate pdbfixer uses."""
    text = _DOCS_PATH.read_text().lower()
    assert "still good for" in text or "still useful" in text or "casual" in text
