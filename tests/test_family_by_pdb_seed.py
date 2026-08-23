"""Verify the shipped family_by_pdb_seed.json is loadable and canonical."""
from __future__ import annotations

from pathlib import Path

from stack_protein_preparation._family_stratification import (
    FAMILY_LABELS,
    load_family_mapping,
)


_SEED = (
    Path(__file__).resolve().parents[1]
    / "src" / "stack_protein_preparation" / "data" / "family_by_pdb_seed.json"
)


def test_seed_file_exists():
    assert _SEED.is_file(), f"family_by_pdb_seed.json missing at {_SEED}"


def test_seed_loads_through_family_stratification_loader():
    mapping = load_family_mapping(_SEED)
    assert isinstance(mapping, dict)
    assert len(mapping) >= 20, "seed should carry at least 20 well-known entries"


def test_seed_labels_are_all_canonical():
    """Every family label in the seed must come from FAMILY_LABELS
    (except __unassigned__ which the seed should never emit)."""
    mapping = load_family_mapping(_SEED)
    canonical = set(FAMILY_LABELS) - {"__unassigned__"}
    seed_labels = set(mapping.values())
    assert seed_labels.issubset(canonical), (
        f"non-canonical labels in seed: {seed_labels - canonical}"
    )


def test_seed_covers_expected_families():
    mapping = load_family_mapping(_SEED)
    labels = set(mapping.values())
    # At minimum: kinase, gpcr, protease, metalloenzyme, nuclear_receptor
    for required in ("kinase", "gpcr", "protease", "metalloenzyme"):
        assert required in labels, f"seed missing family {required}"


def test_seed_pdb_ids_uppercase_and_4char():
    mapping = load_family_mapping(_SEED)
    for pdb in mapping:
        assert pdb == pdb.upper(), f"non-uppercase pdb id in seed: {pdb}"
        assert len(pdb) == 4, f"non-4-char pdb id in seed: {pdb}"


def test_seed_no_duplicate_pdb_ids():
    import json
    raw = json.loads(_SEED.read_text())
    assert len(raw) == len({k.upper() for k in raw})
