"""UniProt-annotated intrinsically-disordered region (IDR) check via MobiDB.

Some AF3-class models confidently fold regions that UniProt / MobiDB /
DisProt annotate as intrinsically disordered (Wang et al. arXiv:
2510.15939 2025).  Because a crystal will never resolve an IDR, filling
an IDR gap with a rigid AF prediction produces a model whose geometry
looks fine locally but is scientifically wrong: MD trajectories will
just melt the region back to disorder.  A gate on top of the pLDDT
filter (Croll IUCr Acta D 2025) catches this: reject any splice window
whose residues overlap an annotated IDR region for the parent UniProt
entry.

Fail-open: if MobiDB is unreachable or the accession is unknown, we
return None (do NOT reject) rather than blocking the whole pipeline
on a network dependency.  The pLDDT gate already filters most IDR
regions incidentally (they usually have low local confidence).
"""
from __future__ import annotations

import json
from pathlib import Path
from urllib.error import HTTPError, URLError
from urllib.request import urlopen


_MOBIDB_URL = "https://mobidb.bio.unipd.it/api/download?acc={acc}&format=json"

# On-disk cache: hand-curated snapshot of MobiDB IDR regions for the
# UniProt accessions that appear in the FRUTON benchmark sets.  Purpose:
# (1) offline reproducibility -- the JCTC reviewer wants runs to be
# self-contained rather than depending on a live web API, (2) reduce
# runtime failure surface when MobiDB is briefly down, (3) shave the
# 2-3 second per-accession network call from the pipeline.  Cache misses
# still hit the live API and can be persisted with save_cache_entry().
_CACHE_PATH = Path(__file__).parent / "data" / "mobidb_snapshot.json"
_CACHE: dict[str, list[tuple[int, int]]] | None = None


def _load_cache() -> dict[str, list[tuple[int, int]]]:
    """Read ``data/mobidb_snapshot.json`` into ``{acc: [(start,end), ...]}``.

    JSON schema:
        {
          "P04637": [[50, 96], [282, 325], [351, 393]],
          "P10636": [[1, 240], ...],
          ...
        }
    Missing / malformed file -> empty dict (fail-open).
    """
    if not _CACHE_PATH.is_file():
        return {}
    try:
        raw = json.loads(_CACHE_PATH.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return {}
    out: dict[str, list[tuple[int, int]]] = {}
    if not isinstance(raw, dict):
        return {}
    for acc, regions in raw.items():
        if not isinstance(regions, list):
            continue
        parsed: list[tuple[int, int]] = []
        for r in regions:
            try:
                s = int(r[0]); e = int(r[1])
            except (TypeError, ValueError, IndexError):
                continue
            if s <= e:
                parsed.append((s, e))
        out[str(acc).upper().strip()] = sorted(parsed)
    return out


def _get_cache() -> dict[str, list[tuple[int, int]]]:
    global _CACHE
    if _CACHE is None:
        _CACHE = _load_cache()
    return _CACHE


def _reset_cache_for_tests() -> None:
    """Force the cache to reload on the next lookup.  Testing only."""
    global _CACHE
    _CACHE = None


def fetch_uniprot_disorder_regions(
    uniprot_id: str,
    timeout_seconds: float = 3.0,
    use_cache: bool = True,
) -> list[tuple[int, int]] | None:
    """Return sorted list of (start, end) 1-based residue ranges annotated
    as disordered in MobiDB for ``uniprot_id``, or ``None`` if the API is
    unavailable / the accession is unknown.

    Cache-first: when ``use_cache`` is True (default) the on-disk
    ``data/mobidb_snapshot.json`` is checked first; a hit returns
    immediately without any network call.  Misses fall through to the
    live MobiDB API.  Pass ``use_cache=False`` to force a live query
    (e.g. cache-refresh scripts).

    We use MobiDB's consensus prediction ``prediction-disorder-mobidb_lite``
    when present, falling back to any curated ``curated-disorder`` region.
    """
    acc = uniprot_id.strip().upper()
    if use_cache:
        cached = _get_cache().get(acc)
        if cached is not None:
            return list(cached)

    url = _MOBIDB_URL.format(acc=acc)
    try:
        with urlopen(url, timeout=timeout_seconds) as response:
            payload = json.load(response)
    except (HTTPError, URLError, TimeoutError, OSError, ValueError):
        return None

    # MobiDB's /api/download endpoint returns a plain dict for a single
    # accession lookup (as of 2026-08); older releases wrapped it in a
    # one-element list.  Accept both shapes.
    if isinstance(payload, list):
        if not payload:
            return []
        record = payload[0]
    elif isinstance(payload, dict):
        record = payload
    else:
        return []

    preferred_keys = (
        "prediction-disorder-mobidb_lite",
        "curated-disorder-priority",
        "curated-disorder-merge",
        "homology-disorder-merge",
    )
    for key in preferred_keys:
        entry = record.get(key)
        if entry is None:
            continue
        regions = entry.get("regions") or []
        parsed: list[tuple[int, int]] = []
        for region in regions:
            try:
                start = int(region[0])
                end = int(region[1])
            except (IndexError, TypeError, ValueError):
                continue
            if start <= end:
                parsed.append((start, end))
        if parsed:
            parsed.sort()
            return parsed
    return []


def gap_overlaps_uniprot_idr(
    uniprot_id: str,
    gap_start: int,
    gap_end: int,
    overlap_fraction_threshold: float = 0.5,
    timeout_seconds: float = 3.0,
) -> bool | None:
    """Return True if the residue range ``[gap_start, gap_end]`` overlaps
    any MobiDB disorder region for ``uniprot_id`` by at least
    ``overlap_fraction_threshold`` of the gap length.  Returns False if
    no overlap.  Returns None if MobiDB is unavailable (fail-open).
    """
    if gap_end < gap_start:
        return False
    regions = fetch_uniprot_disorder_regions(uniprot_id, timeout_seconds)
    if regions is None:
        return None
    gap_len = gap_end - gap_start + 1
    for r_start, r_end in regions:
        overlap = min(gap_end, r_end) - max(gap_start, r_start) + 1
        if overlap > 0 and overlap / gap_len >= overlap_fraction_threshold:
            return True
    return False


def save_cache_entry(uniprot_id: str, regions: list[tuple[int, int]]) -> bool:
    """Persist a single accession's IDR regions to the on-disk snapshot.

    Idempotent: reloads the whole snapshot, updates one entry, writes back.
    Returns True on success, False on I/O error.  Used by cache-refresh
    scripts (scripts/fetch_mobidb_snapshot.py) and not called from the
    hot path -- the pipeline reads the snapshot but never writes it.
    """
    try:
        _CACHE_PATH.parent.mkdir(parents=True, exist_ok=True)
        existing = _load_cache()
        existing[uniprot_id.strip().upper()] = sorted(
            (int(s), int(e)) for s, e in regions if int(s) <= int(e)
        )
        # JSON needs plain lists, not tuples
        serialisable = {
            acc: [list(pair) for pair in regs]
            for acc, regs in existing.items()
        }
        _CACHE_PATH.write_text(
            json.dumps(serialisable, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
        _reset_cache_for_tests()
        return True
    except (OSError, ValueError):
        return False


def cache_contains(uniprot_id: str) -> bool:
    """True iff the on-disk snapshot has an entry for this accession."""
    return uniprot_id.strip().upper() in _get_cache()


def known_cached_accessions() -> tuple[str, ...]:
    """Sorted tuple of accessions in the current snapshot."""
    return tuple(sorted(_get_cache().keys()))
