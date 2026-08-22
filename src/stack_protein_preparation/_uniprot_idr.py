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
from urllib.error import HTTPError, URLError
from urllib.request import urlopen


_MOBIDB_URL = "https://mobidb.bio.unipd.it/api/download?acc={acc}&format=json"


def fetch_uniprot_disorder_regions(
    uniprot_id: str, timeout_seconds: float = 3.0
) -> list[tuple[int, int]] | None:
    """Return sorted list of (start, end) 1-based residue ranges annotated
    as disordered in MobiDB for ``uniprot_id``, or ``None`` if the API is
    unavailable / the accession is unknown.

    We use MobiDB's consensus prediction ``prediction-disorder-mobidb_lite``
    when present, falling back to any curated ``curated-disorder`` region.
    """
    url = _MOBIDB_URL.format(acc=uniprot_id.strip())
    try:
        with urlopen(url, timeout=timeout_seconds) as response:
            payload = json.load(response)
    except (HTTPError, URLError, TimeoutError, OSError, ValueError):
        return None

    if not isinstance(payload, list) or not payload:
        return []
    record = payload[0]

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
