"""Tests for stack_protein_preparation._paper_override_suggest.

Covers the paper_evidence.md → override-suggestion pipeline: regex
heading matcher, per-residue-type hint routing, confidence tagging,
and the fail-open ".suggested.json" writer path.
"""
from __future__ import annotations

import json
from pathlib import Path

import pytest

from stack_protein_preparation import _paper_override_suggest as pos


def _write_md(path: Path, body: str) -> Path:
    path.write_text(body, encoding="utf-8")
    return path


def test_missing_paper_evidence_returns_skipped(tmp_path: Path):
    manifest = pos.suggest_overrides_from_paper_evidence(tmp_path / "nope.md")
    assert manifest["status"] == "skipped"
    assert manifest["overrides"] == []


def test_no_hints_found_when_body_empty(tmp_path: Path):
    md = _write_md(tmp_path / "empty.md", "# nothing here\n")
    manifest = pos.suggest_overrides_from_paper_evidence(md)
    assert manifest["overrides"] == []


def test_asp_catalytic_acid_maps_to_ash(tmp_path: Path):
    md = _write_md(tmp_path / "asp.md", """\
### ASP174
> D174 acts as the catalytic acid, donating a proton to the leaving group.
""")
    manifest = pos.suggest_overrides_from_paper_evidence(md)
    assert manifest["status"] == "success"
    assert len(manifest["overrides"]) == 1
    entry = manifest["overrides"][0]
    assert entry["chain"] == "A"
    assert entry["res_num"] == 174
    assert entry["to"] == "ASH"
    assert "auto-suggested" in entry["reason"]
    assert "high confidence" in entry["reason"]


def test_glu_proton_donor_maps_to_glh(tmp_path: Path):
    md = _write_md(tmp_path / "glu.md", """\
### GLU42
> Glu42 is the proton donor in the acid catalysis mechanism.
""")
    manifest = pos.suggest_overrides_from_paper_evidence(md)
    assert manifest["overrides"][0]["to"] == "GLH"


def test_cys_nucleophile_maps_to_cym(tmp_path: Path):
    md = _write_md(tmp_path / "cys.md", """\
### CYS215
> C215 acts as the nucleophile, attacking the anomeric centre.
""")
    manifest = pos.suggest_overrides_from_paper_evidence(md)
    assert manifest["overrides"][0]["to"] == "CYM"


def test_his_catalytic_base_maps_to_hip(tmp_path: Path):
    md = _write_md(tmp_path / "his.md", """\
### HIS149
> H149 serves as the catalytic base, abstracting a proton from the substrate.
""")
    manifest = pos.suggest_overrides_from_paper_evidence(md)
    assert manifest["overrides"][0]["to"] == "HIP"


def test_his_nd1_protonated_maps_to_hie(tmp_path: Path):
    md = _write_md(tmp_path / "his.md", """\
### HIS85
> H85 has ND1 protonated based on the crystal H-bond network.
""")
    manifest = pos.suggest_overrides_from_paper_evidence(md)
    assert manifest["overrides"][0]["to"] == "HIE"


def test_ambiguous_mention_goes_to_needs_review(tmp_path: Path):
    """Just 'catalytic residue D174' without acid/base attribution -> review."""
    md = _write_md(tmp_path / "ambig.md", """\
### ASP174
> D174 is one of the two catalytic residues.
""")
    manifest = pos.suggest_overrides_from_paper_evidence(md)
    assert len(manifest["overrides"]) == 0
    assert len(manifest["diagnostics"]) == 1
    assert manifest["diagnostics"][0]["status"] == "needs_review"


def test_non_titratable_residue_is_ignored(tmp_path: Path):
    """PHE/TRP/VAL etc are structural; no protonation override applies."""
    md = _write_md(tmp_path / "phe.md", """\
### PHE223
> F223 forms hydrophobic contacts with substrate.
### TRP278
> W278 sits in the -1 subsite.
""")
    manifest = pos.suggest_overrides_from_paper_evidence(md)
    assert manifest["overrides"] == []


def test_multiple_residues_extracted_independently(tmp_path: Path):
    md = _write_md(tmp_path / "multi.md", """\
### CYS215
> C215 acts as the nucleophile in the catalytic mechanism.

### ASP174
> D174 is the catalytic acid that protonates the leaving group.

### HIS149
> H149 acts as the catalytic base for proton abstraction.
""")
    manifest = pos.suggest_overrides_from_paper_evidence(md)
    assert len(manifest["overrides"]) == 3
    tos = {e["res_num"]: e["to"] for e in manifest["overrides"]}
    assert tos == {215: "CYM", 174: "ASH", 149: "HIP"}


def test_chain_default_argument(tmp_path: Path):
    md = _write_md(tmp_path / "chain.md", """\
### CYS215
> C215 is the catalytic nucleophile.
""")
    manifest = pos.suggest_overrides_from_paper_evidence(md, chain_default="B")
    assert manifest["overrides"][0]["chain"] == "B"


def test_write_suggested_overrides_writes_only_when_suggestions_present(tmp_path: Path):
    md = _write_md(tmp_path / "cys.md", """\
### CYS215
> C215 acts as the nucleophile.
""")
    out = tmp_path / "sugg.json"
    manifest = pos.write_suggested_overrides(md, out)
    assert out.is_file()
    payload = json.loads(out.read_text(encoding="utf-8"))
    assert payload["overrides"][0]["to"] == "CYM"
    assert manifest["written_to"] == str(out)


def test_write_suggested_overrides_skips_when_no_hints(tmp_path: Path):
    md = _write_md(tmp_path / "phe.md", """\
### PHE223
> Structural residue only.
""")
    out = tmp_path / "sugg.json"
    manifest = pos.write_suggested_overrides(md, out)
    assert not out.exists()
    assert "written_to" not in manifest


def test_real_5m7u_paper_evidence_style_mineable(tmp_path: Path):
    """Poster-child real-file-shape test based on 5M7U paper_evidence.md
    (O-GlcNAcase; C215 nucleophile, D174/D175 catalytic acid pair)."""
    md = _write_md(tmp_path / "5m7u.md", """\
# 5M7U — Paper Evidence for Active-Site Residues

## All residues mentioned in catalytic/binding context

### ASP174
> Cetinbas et al. Identification of Asp174 and Asp175 as the key
> catalytic acid residues of human O-GlcNAcase.

### CYS215
> C215 acts as the nucleophile in the retaining glycoside hydrolase
> mechanism.

### PHE223
> F223 forms hydrophobic contacts.
""")
    manifest = pos.suggest_overrides_from_paper_evidence(md)
    # ASP174 (catalytic acid) + CYS215 (nucleophile) both extracted, PHE ignored
    tos = {e["res_num"]: e["to"] for e in manifest["overrides"]}
    assert tos.get(174) == "ASH"
    assert tos.get(215) == "CYM"
    assert 223 not in tos
