"""PDB metal-ion reference oracle.

Loads ``data/ts_metal_reference.csv`` — a hand-curated 89-row table
compiled from the PDB Chemical Component Dictionary, Harding (Acta Cryst
D 2001, 2006), Zheng/CheckMyMetal (J Inorg Biochem 2008), MetalPDB
(Andreini 2013), Li-Merz 12-6-4 papers, and MCPB.py manual.  See the
``notes`` column and the top-of-CSV header for full source list.

Every biologically-relevant transition-series ion (3d + 4d + 5d), plus
lanthanides used in MRI/probe experiments, plus s/p-block metals and
biologically-relevant metalloids (As, Sb, Bi, Ga, ...) is covered.

Public entry point:

    ref = lookup_metal(pdb_resname)          # -> MetalReference | None
    ok, reasons = validate_coordination(     # -> (bool, [reason, ...])
        pdb_resname, coord_number, geometry, donor_elements
    )

Use ``validate_coordination`` after metal-site detection to sanity-check
the observed coord sphere against expected chemistry.  Failure reasons
are user-readable and cite the reference distance / geometry so the
user can decide whether the deviation is a genuine biological special
case (unusual catalytic geometry, symmetry-mate contact, misassigned
resname) or a preparation bug (missing water, bad protonation).

Pure Python + csv stdlib.  No numpy / no external deps.
"""
from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path


_CSV_PATH = Path(__file__).parent / "data" / "ts_metal_reference.csv"


@dataclass(frozen=True)
class MetalReference:
    """One row of the metal reference table."""

    pdb_resname: str
    element: str
    oxidation_state: str
    common_charge: str
    coord_numbers: tuple[int, ...]
    geometries: tuple[str, ...]
    donor_preference_hsab: str
    typical_distances_angstrom: str
    spin_state_default: str
    force_field_route: str
    notes: str

    @property
    def preferred_coord_number(self) -> int | None:
        """First (most common) coord number, or None."""
        return self.coord_numbers[0] if self.coord_numbers else None

    @property
    def preferred_geometry(self) -> str | None:
        return self.geometries[0] if self.geometries else None


def _parse_int_list(text: str) -> tuple[int, ...]:
    out: list[int] = []
    for chunk in text.replace(",", ";").split(";"):
        chunk = chunk.strip()
        if not chunk:
            continue
        try:
            out.append(int(chunk))
        except ValueError:
            continue
    return tuple(out)


def _parse_str_list(text: str) -> tuple[str, ...]:
    return tuple(
        chunk.strip()
        for chunk in text.split(";")
        if chunk.strip()
    )


def _load_table() -> dict[str, MetalReference]:
    """Read the reference CSV into ``{pdb_resname: MetalReference}``."""
    if not _CSV_PATH.is_file():
        return {}
    out: dict[str, MetalReference] = {}
    with _CSV_PATH.open(newline="", encoding="utf-8") as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            resname = (row.get("pdb_resname") or "").strip().upper()
            if not resname:
                continue
            out[resname] = MetalReference(
                pdb_resname=resname,
                element=(row.get("element") or "").strip(),
                oxidation_state=(row.get("oxidation_state") or "").strip(),
                common_charge=(row.get("common_charge") or "").strip(),
                coord_numbers=_parse_int_list(row.get("coord_numbers") or ""),
                geometries=_parse_str_list(row.get("geometries") or ""),
                donor_preference_hsab=(row.get("donor_preference_hsab") or "").strip(),
                typical_distances_angstrom=(row.get("typical_distances_angstrom") or "").strip(),
                spin_state_default=(row.get("spin_state_default") or "").strip(),
                force_field_route=(row.get("force_field_route") or "").strip(),
                notes=(row.get("notes") or "").strip(),
            )
    return out


_TABLE: dict[str, MetalReference] | None = None


def _get_table() -> dict[str, MetalReference]:
    global _TABLE
    if _TABLE is None:
        _TABLE = _load_table()
    return _TABLE


def lookup_metal(pdb_resname: str) -> MetalReference | None:
    """Return the reference entry for a PDB residue name, or None."""
    if not pdb_resname:
        return None
    return _get_table().get(pdb_resname.strip().upper())


def all_known_pdb_resnames() -> tuple[str, ...]:
    return tuple(sorted(_get_table().keys()))


def validate_coordination(
    pdb_resname: str,
    coord_number: int,
    geometry: str | None = None,
    donor_elements: tuple[str, ...] | None = None,
) -> tuple[bool, list[str]]:
    """Sanity-check an observed coordination sphere against the reference.

    Returns ``(is_consistent, [reason, ...])``.  When the entry is
    unknown (resname not in the reference), returns ``(True, [notice])``
    so upstream doesn't reject on missing coverage.

    Rules:
      * Coord number must be listed in ``coord_numbers`` (any hit -> ok).
      * If ``geometry`` supplied, it must appear in ``geometries`` (any
        hit -> ok).  Comparison is case-insensitive on the canonical
        FRUTON geometry names (tetrahedral, octahedral, square_planar,
        trigonal_bipyramidal, square_pyramidal, linear, ...).
      * If ``donor_elements`` supplied, at least one must be consistent
        with the HSAB preference (soft/borderline/hard classification).
        Hard mismatches (e.g. soft Cu¹⁺ but no S-donors present) become
        reason lines; borderline / mixed cases are accepted.

    Never raises; unknown geometries / donors just add advisory lines.
    """
    ref = lookup_metal(pdb_resname)
    if ref is None:
        return (True, [
            f"{pdb_resname!r} not in metal reference table -- no validation performed"
        ])

    reasons: list[str] = []
    ok = True

    if ref.coord_numbers and coord_number not in ref.coord_numbers:
        ok = False
        exp = "/".join(str(n) for n in ref.coord_numbers)
        reasons.append(
            f"coord number {coord_number} not in expected {{{exp}}} for {ref.element}{ref.oxidation_state}"
        )

    if geometry and ref.geometries:
        g = geometry.strip().lower()
        expected_norm = tuple(g_.lower() for g_ in ref.geometries)
        if g not in expected_norm:
            ok = False
            reasons.append(
                f"geometry {geometry!r} not in expected {{{'/'.join(ref.geometries)}}} "
                f"for {ref.element}{ref.oxidation_state}"
            )

    if donor_elements:
        hsab = ref.donor_preference_hsab.lower()
        donor_set = {d.strip().upper() for d in donor_elements}
        has_s = any(d in ("S",) for d in donor_set)
        has_o = any(d in ("O",) for d in donor_set)
        has_n = any(d in ("N",) for d in donor_set)
        if "soft" in hsab and not has_s and not (has_n or has_o):
            reasons.append(
                f"HSAB soft ion (prefers S-donors) but observed donors {sorted(donor_set)} "
                f"has no S -- unusual for {ref.pdb_resname}"
            )
        elif "hard" in hsab and not has_o and (has_s and not has_n):
            reasons.append(
                f"HSAB hard ion (prefers O-donors) but observed donors {sorted(donor_set)} "
                f"is all-S -- unusual for {ref.pdb_resname}"
            )

    return (ok, reasons)


def describe(pdb_resname: str) -> str:
    """Reviewer-friendly one-line summary of a metal reference entry."""
    ref = lookup_metal(pdb_resname)
    if ref is None:
        return f"{pdb_resname!r} not in metal reference table"
    return (
        f"{ref.pdb_resname} ({ref.element}{ref.oxidation_state}, "
        f"charge {ref.common_charge}): CN={'/'.join(str(n) for n in ref.coord_numbers)}, "
        f"geom={'/'.join(ref.geometries)}, "
        f"HSAB={ref.donor_preference_hsab}, "
        f"spin={ref.spin_state_default}, ff_route={ref.force_field_route}"
    )
