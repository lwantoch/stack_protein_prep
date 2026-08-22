"""Tests for stack_protein_preparation._uniprot_idr.

Covers the MobiDB API parser, schema-tolerance, gap-overlap decision,
and fail-open behaviour on network errors -- all with a mocked urlopen
so the tests never touch the real network.
"""
from __future__ import annotations

import io
import json
from urllib.error import HTTPError, URLError

import pytest

from stack_protein_preparation import _uniprot_idr


class _FakeResponse:
    """Minimal ``urlopen`` context-manager stand-in."""

    def __init__(self, payload: bytes) -> None:
        self._payload = payload

    def __enter__(self):
        return io.BytesIO(self._payload)

    def __exit__(self, *_exc) -> None:
        return None


def _install_fake_urlopen(monkeypatch: pytest.MonkeyPatch, payload):
    def _fake(url, timeout=None):  # noqa: ARG001
        return _FakeResponse(json.dumps(payload).encode("utf-8"))
    monkeypatch.setattr(_uniprot_idr, "urlopen", _fake)


def _install_broken_urlopen(monkeypatch: pytest.MonkeyPatch, exc: Exception):
    def _fake(url, timeout=None):  # noqa: ARG001
        raise exc
    monkeypatch.setattr(_uniprot_idr, "urlopen", _fake)


# ---------------------------------------------------------------------------
# fetch_uniprot_disorder_regions -- schema tolerance
# ---------------------------------------------------------------------------


def test_fetch_regions_from_dict_response(monkeypatch: pytest.MonkeyPatch):
    """Current MobiDB shape: single dict, disorder key present."""
    payload = {
        "acc": "P04637",
        "prediction-disorder-mobidb_lite": {
            "regions": [[50, 96], [282, 325], [351, 393]],
        },
    }
    _install_fake_urlopen(monkeypatch, payload)
    result = _uniprot_idr.fetch_uniprot_disorder_regions("P04637")
    assert result == [(50, 96), (282, 325), (351, 393)]


def test_fetch_regions_from_list_response(monkeypatch: pytest.MonkeyPatch):
    """Legacy MobiDB shape: one-element list wrapping the record."""
    payload = [{
        "acc": "P04637",
        "prediction-disorder-mobidb_lite": {"regions": [[10, 20]]},
    }]
    _install_fake_urlopen(monkeypatch, payload)
    result = _uniprot_idr.fetch_uniprot_disorder_regions("P04637")
    assert result == [(10, 20)]


def test_fetch_regions_fallback_to_curated(monkeypatch: pytest.MonkeyPatch):
    """If mobidb_lite absent, curated-disorder-priority takes over."""
    payload = {
        "acc": "X99999",
        "curated-disorder-priority": {"regions": [[5, 15], [30, 40]]},
    }
    _install_fake_urlopen(monkeypatch, payload)
    result = _uniprot_idr.fetch_uniprot_disorder_regions("X99999")
    assert result == [(5, 15), (30, 40)]


def test_fetch_regions_empty_when_no_disorder(monkeypatch: pytest.MonkeyPatch):
    payload = {"acc": "Q00000"}
    _install_fake_urlopen(monkeypatch, payload)
    result = _uniprot_idr.fetch_uniprot_disorder_regions("Q00000")
    assert result == []


def test_fetch_regions_returns_none_on_http_error(monkeypatch: pytest.MonkeyPatch):
    _install_broken_urlopen(monkeypatch, HTTPError("u", 500, "s", None, None))
    result = _uniprot_idr.fetch_uniprot_disorder_regions("X")
    assert result is None


def test_fetch_regions_returns_none_on_url_error(monkeypatch: pytest.MonkeyPatch):
    _install_broken_urlopen(monkeypatch, URLError("unreachable"))
    result = _uniprot_idr.fetch_uniprot_disorder_regions("X")
    assert result is None


def test_fetch_regions_returns_none_on_timeout(monkeypatch: pytest.MonkeyPatch):
    _install_broken_urlopen(monkeypatch, TimeoutError("slow"))
    result = _uniprot_idr.fetch_uniprot_disorder_regions("X")
    assert result is None


def test_fetch_regions_skips_malformed_entries(monkeypatch: pytest.MonkeyPatch):
    payload = {
        "acc": "Q1",
        "prediction-disorder-mobidb_lite": {
            "regions": [[1, 10], "junk", [20], [30, 40], [50, 40]],
        },
    }
    _install_fake_urlopen(monkeypatch, payload)
    result = _uniprot_idr.fetch_uniprot_disorder_regions("Q1")
    # Malformed entries dropped; end<start dropped; valid entries kept sorted.
    assert result == [(1, 10), (30, 40)]


# ---------------------------------------------------------------------------
# gap_overlaps_uniprot_idr -- decision logic
# ---------------------------------------------------------------------------


def test_gap_fully_inside_idr_is_true(monkeypatch: pytest.MonkeyPatch):
    _install_fake_urlopen(monkeypatch, {
        "prediction-disorder-mobidb_lite": {"regions": [[50, 96]]}
    })
    assert _uniprot_idr.gap_overlaps_uniprot_idr("Q1", 60, 70) is True


def test_gap_disjoint_from_idr_is_false(monkeypatch: pytest.MonkeyPatch):
    _install_fake_urlopen(monkeypatch, {
        "prediction-disorder-mobidb_lite": {"regions": [[50, 96]]}
    })
    assert _uniprot_idr.gap_overlaps_uniprot_idr("Q1", 200, 210) is False


def test_gap_partial_overlap_under_threshold_is_false(monkeypatch: pytest.MonkeyPatch):
    """25% overlap fails the 50% default threshold."""
    _install_fake_urlopen(monkeypatch, {
        "prediction-disorder-mobidb_lite": {"regions": [[50, 60]]}
    })
    # gap 55-70 (16 residues) overlaps 55-60 (6 residues) = 37% -> False
    assert _uniprot_idr.gap_overlaps_uniprot_idr("Q1", 55, 70) is False


def test_gap_partial_overlap_over_threshold_is_true(monkeypatch: pytest.MonkeyPatch):
    _install_fake_urlopen(monkeypatch, {
        "prediction-disorder-mobidb_lite": {"regions": [[50, 100]]}
    })
    # gap 55-70 (16 residues) fully inside 50-100 -> 100% overlap
    assert _uniprot_idr.gap_overlaps_uniprot_idr("Q1", 55, 70) is True


def test_gap_overlap_returns_none_when_api_fails(monkeypatch: pytest.MonkeyPatch):
    _install_broken_urlopen(monkeypatch, URLError("offline"))
    assert _uniprot_idr.gap_overlaps_uniprot_idr("Q1", 60, 70) is None


def test_gap_with_end_before_start_is_false(monkeypatch: pytest.MonkeyPatch):
    _install_fake_urlopen(monkeypatch, {
        "prediction-disorder-mobidb_lite": {"regions": [[50, 96]]}
    })
    assert _uniprot_idr.gap_overlaps_uniprot_idr("Q1", 70, 60) is False
