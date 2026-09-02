"""Optional user-supplied active-site protonation overrides for FRUTON.

If a per-protein file ``active_site_overrides.json`` exists in the data-dir
(``<data-dir>/<PDB_ID>/active_site_overrides.json``), FRUTON's protonation
stage will honour those manual choices AFTER PROPKA and metal-coordination
overrides.  User's active-site overrides have the highest priority.

Use cases:
  * A published paper explicitly identifies a catalytic His as HIP (positive,
    proton-donor) that PROPKA scored as HIE.
  * A crystallographic paper's mechanism figure shows a CYM thiolate that
    PROPKA missed because the S atom looked bulk-solvent-exposed.
  * A published protein-family assignment (e.g., MAP kinase catalytic Asp
    always ASH) that requires ASP -> ASH consistently.

File format (JSON, stdlib — no PyYAML dependency):

```json
{
  "overrides": [
    {"chain": "A", "res_num": 149, "to": "HIP",
     "reason": "catalytic base per Fig 3 of https://doi.org/10.xxxx/yyy"},
    {"chain": "A", "res_num": 215, "to": "CYM",
     "reason": "catalytic thiolate, JACS 2015 Fig 4b"},
    {"chain": "B", "res_num": 42,  "to": "ASH",
     "reason": "protonated Asp catalytic acid"}
  ]
}
```

Accepted ``to`` values: HID, HIE, HIP, CYM, CYX, ASH, GLH, LYN.

The `reason` field is required — a bare override without justification would
be a silent state-change and risks propagating unverified assumptions
through downstream MD.  FRUTON logs the reason with each applied override.
"""
from __future__ import annotations

import json
from pathlib import Path
from typing import Any


_ACCEPTED_TARGET_STATES = {"HID", "HIE", "HIP", "CYM", "CYX", "ASH", "GLH", "LYN"}


def load_active_site_overrides(pdb_dir: str | Path) -> dict[tuple[str, int], str]:
    """Return a mapping {(chain, res_num): target_state} loaded from
    ``<pdb_dir>/active_site_overrides.json``.  Missing file -> empty dict.

    Invalid entries (missing chain/res_num/to, unknown state, empty reason)
    are silently skipped.  Callers should validate the file exists via
    :func:`overrides_present` if they want to warn on typos.
    """
    p = Path(pdb_dir) / "active_site_overrides.json"
    if not p.is_file():
        return {}
    try:
        data = json.loads(p.read_text(encoding="utf-8"))
    except json.JSONDecodeError:
        return {}
    result: dict[tuple[str, int], str] = {}
    for entry in (data.get("overrides") or []):
        chain = str(entry.get("chain", "")).strip()
        try:
            resi = int(entry.get("res_num"))
        except (TypeError, ValueError):
            continue
        target = str(entry.get("to", "")).strip().upper()
        reason = str(entry.get("reason", "")).strip()
        if not (chain and target and reason):
            continue
        if target not in _ACCEPTED_TARGET_STATES:
            continue
        if len(chain) != 1:
            continue
        result[(chain, resi)] = target
    return result


def overrides_present(pdb_dir: str | Path) -> bool:
    return (Path(pdb_dir) / "active_site_overrides.json").is_file()


def load_active_site_overrides_with_reasons(
    pdb_dir: str | Path,
) -> list[dict[str, Any]]:
    """Return the raw override entries with reason strings, for logging.
    Missing file -> empty list."""
    p = Path(pdb_dir) / "active_site_overrides.json"
    if not p.is_file():
        return []
    try:
        data = json.loads(p.read_text(encoding="utf-8"))
    except json.JSONDecodeError:
        return []
    return list(data.get("overrides") or [])
