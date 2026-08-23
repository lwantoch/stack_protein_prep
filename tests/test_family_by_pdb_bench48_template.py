"""Verify bench48 family template covers every PDB in the bench48 JSON."""
from __future__ import annotations

import json
from pathlib import Path

from stack_protein_preparation._family_stratification import (
    FAMILY_LABELS,
    load_family_mapping,
)


_TEMPLATE = (
    Path(__file__).resolve().parents[1]
    / "src" / "stack_protein_preparation" / "data"
    / "family_by_pdb_bench48.template.json"
)
_BENCH = (
    Path(__file__).resolve().parents[1]
    / "scripts" / "CESGA_SLURM" / "fruton_bench48_full_results.json"
)


def test_template_exists():
    assert _TEMPLATE.is_file()


def test_template_covers_every_bench48_pdb():
    template = load_family_mapping(_TEMPLATE)
    bench = json.loads(_BENCH.read_text())
    bench_pdbs = {r["pdb"].upper() for r in bench}
    template_pdbs = set(template.keys())
    missing = bench_pdbs - template_pdbs
    extra = template_pdbs - bench_pdbs
    assert not missing, f"template missing bench48 PDBs: {sorted(missing)}"
    assert not extra, f"template has extra PDBs not in bench48: {sorted(extra)}"


def test_template_uses_canonical_labels_only():
    template = load_family_mapping(_TEMPLATE)
    valid = set(FAMILY_LABELS)
    non_canonical = {p: f for p, f in template.items() if f not in valid}
    assert not non_canonical, (
        f"template has non-canonical labels: {non_canonical}"
    )


def test_template_defaults_to_unassigned():
    """Template ships with every entry mapped to __unassigned__ so
    populators see 'this is not yet filled in' rather than a stale
    guess."""
    template = load_family_mapping(_TEMPLATE)
    non_default = [p for p, f in template.items() if f != "__unassigned__"]
    assert not non_default, (
        f"template has pre-filled labels (bad — should be all "
        f"__unassigned__ until manually populated): {non_default}"
    )


def test_readme_documents_populate_workflow():
    readme = _TEMPLATE.with_name("family_by_pdb_bench48.README.md")
    assert readme.is_file()
    text = readme.read_text()
    assert "UniProt" in text
    assert "FAMILY_LABELS" in text
    assert "plot_family_stratification.py" in text
