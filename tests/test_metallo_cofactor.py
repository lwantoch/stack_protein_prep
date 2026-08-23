"""Tests for _metallo_cofactor: heme + Fe-S cluster detection & routing.

Every route documented in the module docstring gets a test hitting the
correct method label + confidence — the reviewer chain says exactly
which literature source is used per protein.
"""
from __future__ import annotations

from pathlib import Path

import pytest

from stack_protein_preparation import _metallo_cofactor as mc
from stack_protein_preparation._component_confidence import Confidence


def _write_pdb(path: Path, atoms: list[str]) -> Path:
    path.write_text("\n".join(atoms) + "\nEND\n")
    return path


def _fmt_pdb_atom(record, serial, atom_name, resname, chain, resnum,
                  x, y, z, element):
    """Standard PDB fixed-width line — columns per PDB v3.30 spec."""
    # atom name: 4 chars, right-aligned for 1-char element, left for multi
    if len(atom_name) < 4 and len(element) == 1:
        atom_str = f" {atom_name:<3s}"
    else:
        atom_str = f"{atom_name:<4s}"
    return (
        f"{record:<6s}"           # 0-5   record
        f"{serial:>5d} "          # 6-11  serial + space
        f"{atom_str}"             # 12-15 atom name
        f" "                       # 16    altloc
        f"{resname:>3s}"          # 17-19 resname
        f" "                       # 20    space
        f"{chain:1s}"             # 21    chain
        f"{resnum:>4d}"           # 22-25 resnum
        f"    "                    # 26-29 icode + spaces
        f"{x:>8.3f}"              # 30-37 x
        f"{y:>8.3f}"              # 38-45 y
        f"{z:>8.3f}"              # 46-53 z
        f"  1.00"                  # 54-59 occupancy
        f"  0.00"                  # 60-65 tempfactor
        f"          "              # 66-75 spaces
        f"{element:>2s}"          # 76-77 element
    )


# ---------------------------------------------------------------------------
# Heme axial-type classification
# ---------------------------------------------------------------------------


def _heme_pdb(tmp_path, axial_pairs: list[tuple[str, str]],
              heme_resname="HEM") -> Path:
    """axial_pairs: [(resname, donor_atom)] each placed within 2.4 Å of Fe."""
    lines = [
        _fmt_pdb_atom("HETATM", 1, "FE",  heme_resname, "A", 500,  0.0,  0.0, 0.0, "FE"),
        _fmt_pdb_atom("HETATM", 2, "NA",  heme_resname, "A", 500,  0.0,  2.0, 0.0, "N"),
        _fmt_pdb_atom("HETATM", 3, "NB",  heme_resname, "A", 500,  2.0,  0.0, 0.0, "N"),
        _fmt_pdb_atom("HETATM", 4, "NC",  heme_resname, "A", 500,  0.0, -2.0, 0.0, "N"),
        _fmt_pdb_atom("HETATM", 5, "ND",  heme_resname, "A", 500, -2.0,  0.0, 0.0, "N"),
    ]
    for i, (rn, at) in enumerate(axial_pairs):
        z = 2.4 if i == 0 else -2.4
        elem = at[0]
        lines.append(_fmt_pdb_atom("ATOM", 100+i, at, rn, "A", 50+i, 0.0, 0.0, z, elem))
    return _write_pdb(tmp_path / f"heme_{len(axial_pairs)}.pdb", lines)


def test_bis_his_heme_high_conf_hemall_bis_his(tmp_path):
    pdb = _heme_pdb(tmp_path, [("HIS", "NE2"), ("HIS", "NE2")])
    hemes = mc.detect_heme_systems(pdb)
    assert len(hemes) == 1
    cc = hemes[0].to_component_confidence()
    assert cc.confidence is Confidence.HIGH
    assert cc.method == "autenrieth_2004_hemall_bis_his"
    assert "Autenrieth" in cc.reason


def test_his_met_heme_high_conf(tmp_path):
    pdb = _heme_pdb(tmp_path, [("HIS", "NE2"), ("MET", "SD")])
    hemes = mc.detect_heme_systems(pdb)
    cc = hemes[0].to_component_confidence()
    assert cc.confidence is Confidence.HIGH
    assert cc.method == "autenrieth_2004_hemall_his_met"
    assert "cytochrome c" in cc.reason.lower()


def test_his_only_deoxy_high_conf(tmp_path):
    pdb = _heme_pdb(tmp_path, [("HIS", "NE2")])
    hemes = mc.detect_heme_systems(pdb)
    cc = hemes[0].to_component_confidence()
    assert cc.confidence is Confidence.HIGH
    assert cc.method == "giammona_1984_deoxy"


def test_p450_cys_thiolate_medium_conf(tmp_path):
    pdb = _heme_pdb(tmp_path, [("CYS", "SG")])
    hemes = mc.detect_heme_systems(pdb)
    cc = hemes[0].to_component_confidence()
    assert cc.confidence is Confidence.MEDIUM
    assert cc.method == "shahrokh_2012_p450"
    assert "Shahrokh" in cc.reason
    assert "resting" in cc.suggested_action.lower() or "state" in cc.suggested_action.lower()


def test_orphan_heme_no_axial_low_conf(tmp_path):
    pdb = _heme_pdb(tmp_path, [])
    hemes = mc.detect_heme_systems(pdb)
    cc = hemes[0].to_component_confidence()
    assert cc.confidence is Confidence.LOW
    assert cc.method == "orphan_heme"


def test_hea_variant_johansson_2016(tmp_path):
    pdb = _heme_pdb(tmp_path, [("HIS", "NE2")], heme_resname="HEA")
    hemes = mc.detect_heme_systems(pdb)
    cc = hemes[0].to_component_confidence()
    assert cc.method == "johansson_2016_heme_a"
    assert "cytochrome c oxidase" in cc.reason.lower()


# ---------------------------------------------------------------------------
# HIS override generation (heme → HID)
# ---------------------------------------------------------------------------


def test_his_ne2_axial_generates_hid_override(tmp_path):
    pdb = _heme_pdb(tmp_path, [("HIS", "NE2")])
    hemes = mc.detect_heme_systems(pdb)
    overrides = hemes[0].his_overrides()
    assert len(overrides) == 1
    chain, resnum, icode, forced = overrides[0]
    assert forced == "HID"


def test_his_nd1_axial_generates_hie_override(tmp_path):
    pdb = _heme_pdb(tmp_path, [("HIS", "ND1")])
    hemes = mc.detect_heme_systems(pdb)
    overrides = hemes[0].his_overrides()
    assert overrides[0][3] == "HIE"


# ---------------------------------------------------------------------------
# HEC covalent (CXXCH) detection
# ---------------------------------------------------------------------------


def test_hec_covalent_detected_with_2_cys_thioethers(tmp_path):
    """Fake heme c: HEC residue with CAB/CAC atoms + 2 Cys-SG within 2.5 Å."""
    lines = [
        _fmt_pdb_atom("HETATM", 1, "FE",  "HEC", "A", 500,  0.0,  0.0, 0.0, "FE"),
        _fmt_pdb_atom("HETATM", 2, "NA",  "HEC", "A", 500,  0.0,  2.0, 0.0, "N"),
        _fmt_pdb_atom("HETATM", 3, "NB",  "HEC", "A", 500,  2.0,  0.0, 0.0, "N"),
        _fmt_pdb_atom("HETATM", 4, "NC",  "HEC", "A", 500,  0.0, -2.0, 0.0, "N"),
        _fmt_pdb_atom("HETATM", 5, "ND",  "HEC", "A", 500, -2.0,  0.0, 0.0, "N"),
        _fmt_pdb_atom("HETATM", 6, "CAB", "HEC", "A", 500,  3.0,  0.0, 0.0, "C"),
        _fmt_pdb_atom("HETATM", 7, "CAC", "HEC", "A", 500,  0.0, -3.0, 0.0, "C"),
        _fmt_pdb_atom("ATOM", 100, "SG",  "CYS", "A",  14,  4.2,  0.0, 0.0, "S"),
        _fmt_pdb_atom("ATOM", 200, "SG",  "CYS", "A",  17,  0.0, -4.2, 0.0, "S"),
        _fmt_pdb_atom("ATOM", 300, "NE2", "HIS", "A",  18,  0.0,  0.0, 2.4, "N"),
    ]
    pdb = _write_pdb(tmp_path / "hec.pdb", lines)
    hemes = mc.detect_heme_systems(pdb)
    assert hemes[0].is_covalent_hec is True
    cc = hemes[0].to_component_confidence()
    assert cc.method == "hemall_hec_cxxch_patch"


# ---------------------------------------------------------------------------
# Fe-S cluster detection
# ---------------------------------------------------------------------------


def _fes_cluster_pdb(tmp_path, hetatm_resname="SF4", cys_positions=None,
                     his_positions=None):
    """cys_positions/his_positions: distances in Å (each within 2.6 → bonded)"""
    lines = [
        _fmt_pdb_atom("HETATM", 1, "FE1", hetatm_resname, "A", 500, 0.0, 0.0, 0.0, "FE"),
        _fmt_pdb_atom("HETATM", 2, "FE2", hetatm_resname, "A", 500, 2.7, 0.0, 0.0, "FE"),
        _fmt_pdb_atom("HETATM", 3, "S1",  hetatm_resname, "A", 500, 1.4, 1.4, 0.0, "S"),
        _fmt_pdb_atom("HETATM", 4, "S2",  hetatm_resname, "A", 500, 1.4, -1.4, 0.0, "S"),
    ]
    idx = 100
    for dist in (cys_positions or []):
        lines.append(_fmt_pdb_atom("ATOM", idx, "SG", "CYS", "A", idx, dist, 0.0, 0.0, "S"))
        idx += 1
    for dist in (his_positions or []):
        lines.append(_fmt_pdb_atom("ATOM", idx, "NE2", "HIS", "A", idx, dist, 0.0, 0.0, "N"))
        idx += 1
    return _write_pdb(tmp_path / f"{hetatm_resname}.pdb", lines)


def test_sf4_cluster_with_4_cys_anchors_carvalho_swart(tmp_path):
    pdb = _fes_cluster_pdb(tmp_path, "SF4", cys_positions=[2.3, 2.4, 4.3, 5.0])
    clusters = mc.detect_fes_clusters(pdb)
    assert len(clusters) == 1
    cc = clusters[0].to_component_confidence()
    assert cc.method == "carvalho_swart_2014"
    assert "Carvalho" in cc.reason
    assert cc.confidence is Confidence.MEDIUM


def test_f3s_cluster_recognises_3fe4s(tmp_path):
    pdb = _fes_cluster_pdb(tmp_path, "F3S", cys_positions=[2.3, 2.4])
    clusters = mc.detect_fes_clusters(pdb)
    cc = clusters[0].to_component_confidence()
    assert "3Fe-4S" in cc.name
    assert cc.method == "carvalho_swart_2014"


def test_rieske_2fe2s_2cys_2his_routes_molina_molina(tmp_path):
    pdb = _fes_cluster_pdb(tmp_path, "FES", cys_positions=[2.3, 2.4],
                           his_positions=[2.5, 2.5])
    clusters = mc.detect_fes_clusters(pdb)
    assert clusters[0].is_rieske() is True
    cc = clusters[0].to_component_confidence()
    assert cc.method == "molina_molina_2014_rieske"
    assert "Rieske" in cc.reason


def test_cfn_femoco_forced_low_conf_no_static_frcmod(tmp_path):
    pdb = _fes_cluster_pdb(tmp_path, "CFN", cys_positions=[2.3])
    clusters = mc.detect_fes_clusters(pdb)
    cc = clusters[0].to_component_confidence()
    assert cc.confidence is Confidence.LOW
    assert cc.method == "femoco_qmmm_only"
    assert "qm/mm" in cc.reason.lower() or "QM/MM" in cc.reason


def test_fes_orphan_no_anchors_medium(tmp_path):
    pdb = _fes_cluster_pdb(tmp_path, "FES")
    clusters = mc.detect_fes_clusters(pdb)
    cc = clusters[0].to_component_confidence()
    assert cc.method == "fes_orphan_cluster"


def test_fes_cluster_generates_cym_overrides_for_all_bridging_cys(tmp_path):
    pdb = _fes_cluster_pdb(tmp_path, "SF4", cys_positions=[2.3, 2.4])
    clusters = mc.detect_fes_clusters(pdb)
    overrides = clusters[0].cys_overrides()
    assert len(overrides) == 2
    assert all(forced == "CYM" for (_c, _rn, _ic, forced) in overrides)


def test_no_hetatm_no_clusters(tmp_path):
    (tmp_path / "empty.pdb").write_text("HEADER\nEND\n")
    assert mc.detect_fes_clusters(tmp_path / "empty.pdb") == []
    assert mc.detect_heme_systems(tmp_path / "empty.pdb") == []
