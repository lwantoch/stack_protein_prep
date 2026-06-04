"""Tests for metalls_check pure logic — distance, geometry, contact detection."""
from __future__ import annotations

import math
from pathlib import Path

import pytest

from stack_protein_preparation.metalls_check import (
    DONOR_ELEMENTS,
    METAL_RESIDUE_IDENTITY,
    STANDARD_NONTRANSITION_IONS,
    TRANSITION_METAL_ELEMENTS,
    PDBAtomRecord,
    MetalContact,
    classify_coordination_geometry,
    distance_between_atoms,
    find_metal_contacts,
    is_true_metal_atom,
)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _make_atom(
    element: str,
    residue_name: str,
    record_name: str,
    x: float = 0.0,
    y: float = 0.0,
    z: float = 0.0,
    atom_name: str | None = None,
    chain_id: str = "A",
    residue_number: int = 1,
) -> PDBAtomRecord:
    return PDBAtomRecord(
        serial=1,
        record_name=record_name,
        atom_name=atom_name or element,
        residue_name=residue_name,
        chain_id=chain_id,
        residue_number=residue_number,
        insertion_code="",
        element=element,
        formal_charge=None,
        x=x,
        y=y,
        z=z,
        source_path="",
        source_role="test",
        raw_line="",
    )


def _make_contact(vx: float, vy: float, vz: float, element: str = "N") -> MetalContact:
    return MetalContact(
        metal_label="ZN A 1 ZN",
        metal_element="ZN",
        donor_label=f"HIS A 2 {element}",
        donor_element=element,
        donor_residue_name="HIS",
        donor_chain_id="A",
        donor_residue_number=2,
        donor_atom_name=element,
        distance_angstrom=round(math.sqrt(vx**2 + vy**2 + vz**2), 3),
        donor_source_role="test",
        vector_x=vx,
        vector_y=vy,
        vector_z=vz,
    )


# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

def test_transition_metal_elements_contains_zinc() -> None:
    assert "ZN" in TRANSITION_METAL_ELEMENTS


def test_transition_metal_elements_contains_iron() -> None:
    assert "FE" in TRANSITION_METAL_ELEMENTS


def test_donor_elements_are_n_o_s_se() -> None:
    assert DONOR_ELEMENTS == {"N", "O", "S", "SE"}


def test_metal_residue_identity_zn() -> None:
    element, charge = METAL_RESIDUE_IDENTITY["ZN"]
    assert element == "ZN"
    assert charge == 2


def test_standard_nontransition_ions_includes_na_mg_ca() -> None:
    assert "NA" in STANDARD_NONTRANSITION_IONS
    assert "MG" in STANDARD_NONTRANSITION_IONS
    assert "CA" in STANDARD_NONTRANSITION_IONS


def test_zn_not_in_standard_nontransition_ions() -> None:
    assert "ZN" not in STANDARD_NONTRANSITION_IONS


# ---------------------------------------------------------------------------
# distance_between_atoms
# ---------------------------------------------------------------------------

def test_distance_same_point_is_zero() -> None:
    a = _make_atom("ZN", "ZN", "HETATM", 0, 0, 0)
    assert distance_between_atoms(a, a) == pytest.approx(0.0)


def test_distance_along_x_axis() -> None:
    a = _make_atom("ZN", "ZN", "HETATM", 0, 0, 0)
    b = _make_atom("N", "HIS", "ATOM", 3, 0, 0)
    assert distance_between_atoms(a, b) == pytest.approx(3.0)


def test_distance_3d_pythagoras() -> None:
    a = _make_atom("ZN", "ZN", "HETATM", 0, 0, 0)
    b = _make_atom("N", "HIS", "ATOM", 1, 2, 2)
    assert distance_between_atoms(a, b) == pytest.approx(3.0)


# ---------------------------------------------------------------------------
# is_true_metal_atom
# ---------------------------------------------------------------------------

def test_zn_hetatm_is_true_metal() -> None:
    atom = _make_atom("ZN", "ZN", "HETATM")
    assert is_true_metal_atom(atom) is True


def test_ca_backbone_atom_not_metal() -> None:
    atom = _make_atom("CA", "ALA", "ATOM", atom_name="CA")
    assert is_true_metal_atom(atom) is False


def test_n_donor_not_metal() -> None:
    atom = _make_atom("N", "HIS", "ATOM", atom_name="NE2")
    assert is_true_metal_atom(atom) is False


def test_fe_hetatm_is_metal() -> None:
    atom = _make_atom("FE", "FE", "HETATM")
    assert is_true_metal_atom(atom) is True


def test_cd_backbone_atom_not_metal() -> None:
    # CD = delta carbon, ATOM record — should NOT be treated as cadmium
    atom = _make_atom("CD", "LEU", "ATOM", atom_name="CD")
    assert is_true_metal_atom(atom) is False


# ---------------------------------------------------------------------------
# classify_coordination_geometry
# ---------------------------------------------------------------------------

def _metal_atom() -> PDBAtomRecord:
    return _make_atom("ZN", "ZN", "HETATM", 0, 0, 0)


def test_zero_contacts_is_uncoordinated() -> None:
    assert classify_coordination_geometry(
        metal_atom=_metal_atom(), contacts=()
    ) == "uncoordinated"


def test_one_contact_is_monodentate() -> None:
    assert classify_coordination_geometry(
        metal_atom=_metal_atom(),
        contacts=(_make_contact(2, 0, 0),),
    ) == "monodentate_contact"


def test_three_contacts_is_trigonal_planar() -> None:
    assert classify_coordination_geometry(
        metal_atom=_metal_atom(),
        contacts=(
            _make_contact(2, 0, 0),
            _make_contact(-1, 1.7, 0),
            _make_contact(-1, -1.7, 0),
        ),
    ) == "trigonal_planar"


def test_six_contacts_is_octahedral() -> None:
    contacts = tuple(
        _make_contact(*v)
        for v in [(2, 0, 0), (-2, 0, 0), (0, 2, 0), (0, -2, 0), (0, 0, 2), (0, 0, -2)]
    )
    assert classify_coordination_geometry(
        metal_atom=_metal_atom(), contacts=contacts
    ) == "octahedral"


def test_high_coordination_returns_labeled_string() -> None:
    contacts = tuple(_make_contact(i, 0, 0) for i in range(1, 9))
    result = classify_coordination_geometry(metal_atom=_metal_atom(), contacts=contacts)
    assert result.startswith("high_coordination_")


# ---------------------------------------------------------------------------
# find_metal_contacts — synthetic structure
# ---------------------------------------------------------------------------

def test_find_contacts_returns_donors_within_cutoff(tmp_path: Path) -> None:
    zn = _make_atom("ZN", "ZN", "HETATM", 0, 0, 0, residue_number=100)
    near_n = _make_atom("N", "HIS", "ATOM", 2.0, 0, 0, atom_name="NE2", residue_number=1)
    far_n = _make_atom("N", "HIS", "ATOM", 5.0, 0, 0, atom_name="NE2", residue_number=2)

    contacts = find_metal_contacts(
        metal_atom=zn,
        atom_records=[zn, near_n, far_n],
        contact_cutoff_angstrom=3.5,
    )
    assert len(contacts) == 1
    assert contacts[0].donor_element == "N"
    assert contacts[0].distance_angstrom == pytest.approx(2.0, abs=0.01)


def test_find_contacts_excludes_metal_atoms() -> None:
    zn1 = _make_atom("ZN", "ZN", "HETATM", 0, 0, 0, residue_number=100)
    zn2 = _make_atom("ZN", "ZN", "HETATM", 1, 0, 0, residue_number=101)
    contacts = find_metal_contacts(
        metal_atom=zn1,
        atom_records=[zn1, zn2],
        contact_cutoff_angstrom=3.5,
    )
    assert contacts == []


def test_find_contacts_excludes_non_donor_elements() -> None:
    zn = _make_atom("ZN", "ZN", "HETATM", 0, 0, 0, residue_number=100)
    carbon = _make_atom("C", "ALA", "ATOM", 2.0, 0, 0, atom_name="CB", residue_number=1)
    contacts = find_metal_contacts(
        metal_atom=zn,
        atom_records=[zn, carbon],
        contact_cutoff_angstrom=3.5,
    )
    assert contacts == []
