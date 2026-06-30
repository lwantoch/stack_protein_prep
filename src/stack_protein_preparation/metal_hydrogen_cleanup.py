"""Remove or reorient hydrogen atoms that point toward coordinated metal ions."""
from __future__ import annotations

import csv
import json
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from stack_protein_preparation.metalls_check import (
    DONOR_ELEMENTS,
    TRANSITION_METAL_ELEMENTS,
    WATER_RESIDUE_NAMES,
    PDBAtomRecord,
    distance_between_atoms,
    find_metal_contacts,
    is_true_metal_atom,
    read_pdb_atom_records,
)

# ---------------------------------------------------------------------------
# Residue-level donor-atom rules
# ---------------------------------------------------------------------------

STANDARD_RESIDUE_DONOR_ATOMS: dict[str, set[str]] = {
    "HIS": {"ND1", "NE2"},
    "CYS": {"SG"},
    "SEC": {"SE"},
    "MET": {"SD"},
    "TYR": {"OH"},
    "ASP": {"OD1", "OD2"},
    "GLU": {"OE1", "OE2"},
    "SER": {"OG"},
    "THR": {"OG1"},
    "ASN": {"OD1"},
    "GLN": {"OE1"},
    "LYS": {"NZ"},
    "ARG": {"NH1", "NH2"},
}

STANDARD_AMINO_ACID_NAMES: set[str] = {
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "HIE",
    "HID", "HIP", "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR",
    "TRP", "TYR", "VAL", "SEC", "PYL",
}


# ---------------------------------------------------------------------------
# Local angle helper (reimplement: _angle_between_vectors is private)
# ---------------------------------------------------------------------------

def _calc_angle_deg(
    v1: tuple[float, float, float],
    v2: tuple[float, float, float],
) -> float:
    """Return the angle in degrees between two vectors from the origin."""
    dot = sum(a * b for a, b in zip(v1, v2))
    n1 = math.sqrt(sum(a * a for a in v1))
    n2 = math.sqrt(sum(b * b for b in v2))
    if n1 == 0.0 or n2 == 0.0:
        return 0.0
    return math.degrees(math.acos(max(-1.0, min(1.0, dot / (n1 * n2)))))


def _angle_at_vertex(
    a: PDBAtomRecord,
    b: PDBAtomRecord,
    c: PDBAtomRecord,
) -> float:
    """Return the angle A-B-C in degrees (vertex at B)."""
    v1 = (a.x - b.x, a.y - b.y, a.z - b.z)
    v2 = (c.x - b.x, c.y - b.y, c.z - b.z)
    return _calc_angle_deg(v1, v2)


# ---------------------------------------------------------------------------
# Public result dataclass
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class HydrogenCleanupResult:
    """Result of one hydrogen-cleanup pass for a metal-containing environment PDB.

    This dataclass is returned by :func:`remove_metal_coordinating_hydrogens`
    and collects paths to all output artefacts together with summary counts so
    that downstream pipeline steps can record provenance without re-reading the
    files.

    Fields
    ------
    status : str
        ``"success"`` — cleanup ran and the cleaned PDB was written.
        ``"skipped"`` — no transition-metal atoms were found; the input PDB was
        copied unchanged and all counts are zero.
        ``"failed"``  — an unexpected exception prevented completion; the
        ``message`` field carries the error text.
    message : str
        Human-readable summary, always present.
    input_pdb_path : Path
        Absolute path to the PDB file that was read.
    cleaned_pdb_path : Path
        Absolute path to the cleaned PDB (input with offending H atoms removed).
    removed_tsv_path : Path
        TSV listing every hydrogen that was removed.
    kept_tsv_path : Path
        TSV listing every hydrogen that was evaluated but kept.
    manual_review_tsv_path : Path
        TSV listing hydrogens on non-standard ligand residues that may need
        manual inspection regardless of the automated decision.
    report_json_path : Path
        JSON report with all counts and parameter echo.
    removed_count : int
        Number of hydrogen atoms removed from the cleaned PDB.
    kept_count : int
        Number of hydrogen atoms evaluated but kept.
    manual_review_count : int
        Number of hydrogen atoms flagged for manual review (subset of removed
        or kept that belong to non-standard residues).
    transition_metal_count : int
        Number of transition-metal atoms found in the input PDB.
    """

    status: str
    message: str
    input_pdb_path: Path
    cleaned_pdb_path: Path
    removed_tsv_path: Path
    kept_tsv_path: Path
    manual_review_tsv_path: Path
    report_json_path: Path
    removed_count: int
    kept_count: int
    manual_review_count: int
    transition_metal_count: int


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

def _build_residue_lookup(
    records: list[PDBAtomRecord],
) -> dict[tuple[str, int, str], list[PDBAtomRecord]]:
    """Return {(chain_id, residue_number, insertion_code): [atoms]} mapping."""
    lookup: dict[tuple[str, int, str], list[PDBAtomRecord]] = {}
    for atom in records:
        key = (atom.chain_id, atom.residue_number, atom.insertion_code)
        lookup.setdefault(key, []).append(atom)
    return lookup


def _removal_reason(
    angle_mdh: float,
    dist_m_h: float,
    dist_m_d: float,
    angle_max: float,
    h_max: float,
    slack: float,
) -> str:
    """Return a compact reason string describing why a hydrogen is being removed."""
    reasons = []
    if angle_mdh <= angle_max:
        reasons.append(f"angle_MDH={angle_mdh:.1f}<=angle_max={angle_max}")
    if dist_m_h <= h_max:
        reasons.append(f"dist_M_H={dist_m_h:.3f}<=h_max={h_max}")
    if dist_m_h < dist_m_d + slack:
        reasons.append(f"dist_M_H={dist_m_h:.3f}<dist_M_D+slack={dist_m_d + slack:.3f}")
    return "; ".join(reasons) if reasons else "geometry"


def _write_tsv(path: Path, header: list[str], rows: list[list[Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(header)
        writer.writerows(rows)


# ---------------------------------------------------------------------------
# Public function
# ---------------------------------------------------------------------------

def remove_metal_coordinating_hydrogens(
    *,
    input_pdb_path: str | Path,
    output_dir: str | Path,
    pdb_id: str,
    variant_label: str,
    contact_cutoff_angstrom: float = 3.5,
    nh_oh_cutoff_angstrom: float = 1.25,
    sh_cutoff_angstrom: float = 1.45,
    m_d_h_angle_max_deg: float = 60.0,
    m_h_distance_max_angstrom: float = 2.10,
    m_h_slack_angstrom: float = 0.35,
    keep_water_hydrogens: bool = True,
) -> HydrogenCleanupResult:
    """Remove hydrogens that are geometrically mis-placed between a transition
    metal and its donor atom.

    Background
    ----------
    After protonation programs such as GROMACS ``pdb2gmx`` or Reduce add
    hydrogen atoms to a structure, metal-coordinating residues often receive
    chemically incorrect H atoms on the donor atom (e.g., H on the SG of a
    zinc-bound CYS, or H on the ND1/NE2 of a histidine coordinated to iron).
    These spurious hydrogens distort force-field energetics and must be removed
    before MCPB parametrization or any QM/MM calculation.

    This function applies three geometric criteria to decide whether to remove
    a hydrogen ``H`` bonded to donor ``D`` that is in contact with transition
    metal ``M``:

    1. **Angle criterion** – the M-D-H angle is ≤ ``m_d_h_angle_max_deg``
       (default 60 °): the H points toward the metal.
    2. **Absolute distance criterion** – the M↔H distance is ≤
       ``m_h_distance_max_angstrom`` (default 2.10 Å): the H is inside the
       metal coordination sphere.
    3. **Slack distance criterion** – the M↔H distance is less than the
       M↔D distance plus ``m_h_slack_angstrom`` (default 0.35 Å): H is closer
       to (or nearly as close as) the donor to the metal.

    Any one of the three criteria triggers removal.

    Residue-specific safeguards
    ---------------------------
    * Only donor atoms whose ``atom_name`` is in :data:`STANDARD_RESIDUE_DONOR_ATOMS`
      for the residue are eligible for hydrogen removal (e.g., only SG for CYS,
      only ND1/NE2 for HIS).  For non-standard residues (HETATM ligands), the
      geometry rule is applied but the decision is also recorded in the
      ``manual_review.tsv`` file.
    * If ``keep_water_hydrogens`` is ``True`` (default) hydrogens on water
      molecules are never removed regardless of geometry.

    Parameters
    ----------
    input_pdb_path:
        PDB file to read.  This should be the environment PDB (protein + water +
        ligand) *without* the metal atom already removed; the metal is expected to
        still be present so that donor contacts can be computed, but the cleanup
        targets only hydrogens attached to donor atoms.
    output_dir:
        Directory where all output files are written.  Created if absent.
    pdb_id:
        Four-character PDB code used as a file-name prefix.
    variant_label:
        Short label for this preparation variant, used in file names and the
        JSON report.
    contact_cutoff_angstrom:
        Maximum metal-to-donor distance for a donor to be considered
        coordinating (default 3.5 Å).
    nh_oh_cutoff_angstrom:
        Maximum D↔H bond distance for N/O donor atoms (default 1.25 Å).
    sh_cutoff_angstrom:
        Maximum D↔H bond distance for S/Se donor atoms (default 1.45 Å).
    m_d_h_angle_max_deg:
        Maximum M-D-H angle for the angle criterion (default 60 °).
    m_h_distance_max_angstrom:
        Maximum M↔H distance for the absolute distance criterion (default
        2.10 Å).
    m_h_slack_angstrom:
        Slack added to M↔D when checking the slack distance criterion (default
        0.35 Å).
    keep_water_hydrogens:
        If ``True``, hydrogens on water residues are never removed (default
        ``True``).

    Returns
    -------
    HydrogenCleanupResult
        Frozen dataclass with paths to all written files and summary counts.
    """
    input_pdb_path = Path(input_pdb_path).resolve()
    output_dir = Path(output_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    cleaned_pdb_path = output_dir / "environment_after_cleanup.pdb"
    removed_tsv_path = output_dir / "removed_metal_facing_hydrogens.tsv"
    kept_tsv_path = output_dir / "kept_hydrogens.tsv"
    manual_review_tsv_path = output_dir / "manual_review.tsv"
    report_json_path = output_dir / "metal_hydrogen_cleanup.json"

    def _make_result(
        status: str,
        message: str,
        removed_count: int = 0,
        kept_count: int = 0,
        manual_review_count: int = 0,
        transition_metal_count: int = 0,
    ) -> HydrogenCleanupResult:
        return HydrogenCleanupResult(
            status=status,
            message=message,
            input_pdb_path=input_pdb_path,
            cleaned_pdb_path=cleaned_pdb_path,
            removed_tsv_path=removed_tsv_path,
            kept_tsv_path=kept_tsv_path,
            manual_review_tsv_path=manual_review_tsv_path,
            report_json_path=report_json_path,
            removed_count=removed_count,
            kept_count=kept_count,
            manual_review_count=manual_review_count,
            transition_metal_count=transition_metal_count,
        )

    try:
        records = read_pdb_atom_records(input_pdb_path, source_role="cleanup")
    except Exception as exc:
        _write_tsv(removed_tsv_path, _REMOVED_HEADER, [])
        _write_tsv(kept_tsv_path, _KEPT_HEADER, [])
        _write_tsv(manual_review_tsv_path, _MANUAL_HEADER, [])
        result = _make_result("failed", f"Could not read input PDB: {exc!r}")
        _write_report_json(report_json_path, result, pdb_id, variant_label, {
            "contact_cutoff_angstrom": contact_cutoff_angstrom,
            "m_d_h_angle_max_deg": m_d_h_angle_max_deg,
            "m_h_distance_max_angstrom": m_h_distance_max_angstrom,
            "m_h_slack_angstrom": m_h_slack_angstrom,
            "keep_water_hydrogens": keep_water_hydrogens,
        })
        return result

    metals = [
        a for a in records
        if is_true_metal_atom(a) and a.element in TRANSITION_METAL_ELEMENTS
    ]

    if not metals:
        # Copy input unchanged
        cleaned_pdb_path.write_bytes(input_pdb_path.read_bytes())
        _write_tsv(removed_tsv_path, _REMOVED_HEADER, [])
        _write_tsv(kept_tsv_path, _KEPT_HEADER, [])
        _write_tsv(manual_review_tsv_path, _MANUAL_HEADER, [])
        result = _make_result(
            "skipped",
            "No transition-metal atoms found; input PDB copied unchanged.",
            transition_metal_count=0,
        )
        _write_report_json(report_json_path, result, pdb_id, variant_label, {
            "contact_cutoff_angstrom": contact_cutoff_angstrom,
            "m_d_h_angle_max_deg": m_d_h_angle_max_deg,
            "m_h_distance_max_angstrom": m_h_distance_max_angstrom,
            "m_h_slack_angstrom": m_h_slack_angstrom,
            "keep_water_hydrogens": keep_water_hydrogens,
        })
        return result

    residue_lookup = _build_residue_lookup(records)

    # Tracking sets / row collectors
    decided_serials: set[int] = set()
    remove_serials: set[int] = set()
    removed_rows: list[list[Any]] = []
    kept_rows: list[list[Any]] = []
    manual_rows: list[list[Any]] = []

    for metal in metals:
        # Find coordinating donors
        contacts = find_metal_contacts(
            metal_atom=metal,
            atom_records=records,
            contact_cutoff_angstrom=contact_cutoff_angstrom,
        )

        for contact in contacts:
            # Reconstruct donor PDBAtomRecord from records matching label
            donor_candidates = [
                a for a in records
                if (
                    a.chain_id == contact.donor_chain_id
                    and a.residue_number == contact.donor_residue_number
                    and a.atom_name.strip() == contact.donor_atom_name.strip()
                    and a.residue_name == contact.donor_residue_name
                )
            ]
            if not donor_candidates:
                continue
            donor = donor_candidates[0]

            if is_true_metal_atom(donor):
                continue

            # Water guard
            if keep_water_hydrogens and donor.residue_name in WATER_RESIDUE_NAMES:
                continue

            # Residue-specific donor atom filter
            is_standard = donor.residue_name in STANDARD_AMINO_ACID_NAMES
            if is_standard:
                allowed_donor_atoms = STANDARD_RESIDUE_DONOR_ATOMS.get(donor.residue_name, set())
                if donor.atom_name.strip() not in allowed_donor_atoms:
                    # Donor atom not a canonical coordinating atom; skip entirely
                    continue

            # Determine bond-length cutoff for H detection
            if donor.element in {"S", "SE"}:
                h_bond_cutoff = sh_cutoff_angstrom
            else:
                h_bond_cutoff = nh_oh_cutoff_angstrom

            # Find hydrogens attached to this donor in the same residue
            res_key = (donor.chain_id, donor.residue_number, donor.insertion_code)
            residue_atoms = residue_lookup.get(res_key, [])
            attached_hydrogens = [
                a for a in residue_atoms
                if a.element == "H" and distance_between_atoms(donor, a) <= h_bond_cutoff
            ]

            dist_m_d = distance_between_atoms(metal, donor)

            for h_atom in attached_hydrogens:
                if h_atom.serial in decided_serials:
                    continue
                decided_serials.add(h_atom.serial)

                dist_m_h = distance_between_atoms(metal, h_atom)
                angle_mdh = _angle_at_vertex(metal, donor, h_atom)

                remove_h = (
                    angle_mdh <= m_d_h_angle_max_deg
                    or dist_m_h <= m_h_distance_max_angstrom
                    or dist_m_h < dist_m_d + m_h_slack_angstrom
                )

                is_ligand = not is_standard

                if remove_h:
                    remove_serials.add(h_atom.serial)
                    reason = _removal_reason(
                        angle_mdh, dist_m_h, dist_m_d,
                        m_d_h_angle_max_deg, m_h_distance_max_angstrom, m_h_slack_angstrom,
                    )
                    removed_rows.append([
                        metal.element,
                        metal.serial,
                        metal.atom_label,
                        donor.element,
                        donor.residue_name,
                        donor.chain_id,
                        donor.residue_number,
                        donor.atom_name.strip(),
                        h_atom.serial,
                        h_atom.atom_name.strip(),
                        f"{dist_m_d:.3f}",
                        f"{dist_m_h:.3f}",
                        f"{angle_mdh:.1f}",
                        reason,
                    ])
                    if is_ligand:
                        manual_rows.append([
                            metal.element,
                            metal.serial,
                            donor.residue_name,
                            donor.chain_id,
                            donor.residue_number,
                            donor.atom_name.strip(),
                            h_atom.serial,
                            h_atom.atom_name.strip(),
                            f"{dist_m_h:.3f}",
                            f"{angle_mdh:.1f}",
                            "yes",
                            "non-standard ligand residue; geometry criterion applied",
                        ])
                else:
                    reason_kept = "geometry_criteria_not_met"
                    kept_rows.append([
                        metal.element,
                        metal.serial,
                        donor.atom_name.strip(),
                        h_atom.serial,
                        h_atom.atom_name.strip(),
                        f"{dist_m_h:.3f}",
                        f"{angle_mdh:.1f}",
                        reason_kept,
                    ])
                    if is_ligand:
                        manual_rows.append([
                            metal.element,
                            metal.serial,
                            donor.residue_name,
                            donor.chain_id,
                            donor.residue_number,
                            donor.atom_name.strip(),
                            h_atom.serial,
                            h_atom.atom_name.strip(),
                            f"{dist_m_h:.3f}",
                            f"{angle_mdh:.1f}",
                            "no",
                            "non-standard ligand residue; geometry criteria not met",
                        ])

    # Write cleaned PDB
    _write_cleaned_pdb(input_pdb_path, cleaned_pdb_path, remove_serials)

    # Reorient metal-coordinating water hydrogens.  Removing them would leave
    # an OH- where the structure must have HOH; flipping them keeps the water
    # complete with both Hs facing away from the metal (the correct geometry).
    waters_reoriented = reorient_metal_coordinating_water_hydrogens(cleaned_pdb_path)

    # Write TSVs
    _write_tsv(removed_tsv_path, _REMOVED_HEADER, removed_rows)
    _write_tsv(kept_tsv_path, _KEPT_HEADER, kept_rows)
    _write_tsv(manual_review_tsv_path, _MANUAL_HEADER, manual_rows)

    result = _make_result(
        "success",
        (
            f"Removed {len(remove_serials)} hydrogen(s) and reoriented "
            f"{waters_reoriented} water(s) from {len(metals)} transition-metal site(s)."
        ),
        removed_count=len(remove_serials),
        kept_count=len(kept_rows),
        manual_review_count=len(manual_rows),
        transition_metal_count=len(metals),
    )
    _write_report_json(report_json_path, result, pdb_id, variant_label, {
        "contact_cutoff_angstrom": contact_cutoff_angstrom,
        "m_d_h_angle_max_deg": m_d_h_angle_max_deg,
        "m_h_distance_max_angstrom": m_h_distance_max_angstrom,
        "m_h_slack_angstrom": m_h_slack_angstrom,
        "keep_water_hydrogens": keep_water_hydrogens,
        "waters_reoriented": waters_reoriented,
    })
    return result


# ---------------------------------------------------------------------------
# TSV column headers
# ---------------------------------------------------------------------------

_REMOVED_HEADER = [
    "metal_element",
    "metal_serial",
    "metal_label",
    "donor_element",
    "donor_residue_name",
    "donor_chain_id",
    "donor_residue_number",
    "donor_atom_name",
    "h_serial",
    "h_atom_name",
    "dist_M_D_ang",
    "dist_M_H_ang",
    "angle_MDH_deg",
    "removal_reason",
]

_KEPT_HEADER = [
    "metal_element",
    "metal_serial",
    "donor_atom_name",
    "h_serial",
    "h_atom_name",
    "dist_M_H_ang",
    "angle_MDH_deg",
    "reason_kept",
]

_MANUAL_HEADER = [
    "metal_element",
    "metal_serial",
    "donor_residue_name",
    "donor_chain_id",
    "donor_residue_number",
    "donor_atom_name",
    "h_serial",
    "h_atom_name",
    "dist_M_H_ang",
    "angle_MDH_deg",
    "removed",
    "note",
]


# ---------------------------------------------------------------------------
# File-writing helpers
# ---------------------------------------------------------------------------

def _write_cleaned_pdb(
    input_path: Path,
    output_path: Path,
    remove_serials: set[int],
) -> None:
    """Write the cleaned PDB, skipping lines whose serial is in remove_serials."""
    with input_path.open("r", encoding="utf-8", errors="replace") as src, \
            output_path.open("w", encoding="utf-8") as dst:
        for line in src:
            record = line[:6]
            if record in {"ATOM  ", "HETATM"}:
                try:
                    serial = int(line[6:11].strip())
                except ValueError:
                    dst.write(line)
                    continue
                if serial in remove_serials:
                    continue
            dst.write(line)


def reorient_metal_coordinating_water_hydrogens(
    pdb_path: Path,
    *,
    metal_o_max_angstrom: float = 2.8,
    h_m_clash_max_angstrom: float = 2.5,
    oh_bond_length: float = 0.96,
    hoh_angle_deg: float = 104.5,
) -> int:
    """Rotate water hydrogens that point toward a coordinated metal so they
    point away.  Returns the number of waters reoriented (in-place).

    pdb2gmx/propka place water hydrogens based on local hydrogen-bond geometry
    and have no concept of metal coordination, so a water whose oxygen donates
    a lone pair to a transition metal often ends up with one or both H atoms
    pointing AT the metal.  That is electrostatically nonsensical and
    catastrophic for SCF convergence in the MCPB small/large model.

    The fix: for each HOH whose O is within ``metal_o_max_angstrom`` of a
    transition metal and whose nearest H is within ``h_m_clash_max_angstrom``
    of the same metal, recompute both H positions so the H–O–H bisector points
    from the metal to the oxygen (i.e. both Hs face away from the metal) and
    the H–O–H angle is the standard 104.5°.  The original H–O–H plane is
    preserved when possible (we keep the perpendicular component from the
    original H1 direction); when both Hs are essentially along the M–O axis
    (degenerate case) we fall back to a deterministic perpendicular.
    """
    import math
    from pathlib import Path as _Path

    pdb_path = _Path(pdb_path)
    if not pdb_path.is_file():
        return 0

    lines = pdb_path.read_text(encoding="utf-8", errors="replace").splitlines(keepends=True)

    # Parse heavy atoms and waters.  We can't key waters by
    # (chain, resseq, icode) alone — pdb4amber and similar preprocessors often
    # strip chain IDs from waters AND reset the resseq, so the same
    # (chain, resseq, icode) appears for multiple physically distinct waters.
    # Instead we treat each consecutive O–H–H block as one water, restarting
    # on every "O" record we see for a water residue.
    metals: list[tuple[float, float, float]] = []
    water_atoms: list[dict[str, object]] = []
    current_water: dict[str, object] | None = None
    for idx, ln in enumerate(lines):
        if not (ln.startswith("ATOM") or ln.startswith("HETATM")):
            continue
        try:
            x = float(ln[30:38]); y = float(ln[38:46]); z = float(ln[46:54])
        except ValueError:
            continue
        resname = ln[17:20].strip().upper()
        atomname = ln[12:16].strip().upper()
        element = ln[76:78].strip().upper() if len(ln) >= 78 else ""
        if resname in TRANSITION_METAL_ELEMENTS:
            metals.append((x, y, z))
            current_water = None
            continue
        if resname not in WATER_RESIDUE_NAMES:
            current_water = None
            continue
        is_oxygen = atomname in {"O", "OW", "OH2"} or element == "O"
        is_hydrogen = (atomname.startswith("H") and atomname != "HG") or element == "H"
        if is_oxygen:
            current_water = {"o": (x, y, z), "o_idx": idx}
            water_atoms.append(current_water)
        elif is_hydrogen and current_water is not None:
            if "h1" not in current_water:
                current_water["h1"] = (x, y, z); current_water["h1_idx"] = idx
            elif "h2" not in current_water:
                current_water["h2"] = (x, y, z); current_water["h2_idx"] = idx
    if not metals:
        return 0

    def _vsub(a, b): return (a[0]-b[0], a[1]-b[1], a[2]-b[2])
    def _vadd(a, b): return (a[0]+b[0], a[1]+b[1], a[2]+b[2])
    def _vscale(a, s): return (a[0]*s, a[1]*s, a[2]*s)
    def _vdot(a, b): return a[0]*b[0]+a[1]*b[1]+a[2]*b[2]
    def _vnorm(a):
        n = math.sqrt(_vdot(a, a))
        return (a[0]/n, a[1]/n, a[2]/n) if n > 1e-9 else (1.0, 0.0, 0.0)
    def _vlen(a): return math.sqrt(_vdot(a, a))

    reoriented = 0
    half_angle = math.radians(hoh_angle_deg / 2.0)
    cos_h = math.cos(half_angle)
    sin_h = math.sin(half_angle)

    for w in water_atoms:
        if "o" not in w or "h1" not in w or "h2" not in w:
            continue
        # Nearest metal to this water's O.
        o = w["o"]
        best_m = min(metals, key=lambda m: _vlen(_vsub(o, m)))
        d_om = _vlen(_vsub(o, best_m))
        if d_om > metal_o_max_angstrom:
            continue
        d_mh1 = _vlen(_vsub(w["h1"], best_m))
        d_mh2 = _vlen(_vsub(w["h2"], best_m))
        if min(d_mh1, d_mh2) > h_m_clash_max_angstrom:
            continue

        # Bisector: from metal toward oxygen (Hs should sit on this side).
        bisector = _vnorm(_vsub(o, best_m))
        # Pick the H further from the metal to preserve plane.
        far_h = w["h1"] if d_mh1 >= d_mh2 else w["h2"]
        rel = _vsub(far_h, o)
        rel_along = _vscale(bisector, _vdot(rel, bisector))
        perp_raw = _vsub(rel, rel_along)
        if _vlen(perp_raw) < 1e-4:
            # Degenerate: pick deterministic perp.
            ref = (1.0, 0.0, 0.0) if abs(bisector[0]) < 0.9 else (0.0, 1.0, 0.0)
            perp_raw = _vsub(ref, _vscale(bisector, _vdot(ref, bisector)))
        perp = _vnorm(perp_raw)

        new_h1 = _vadd(o, _vscale(_vadd(_vscale(bisector, cos_h), _vscale(perp, sin_h)), oh_bond_length))
        new_h2 = _vadd(o, _vscale(_vadd(_vscale(bisector, cos_h), _vscale(perp, -sin_h)), oh_bond_length))

        for idx_key, new_xyz in (("h1_idx", new_h1), ("h2_idx", new_h2)):
            idx = w[idx_key]
            ln = lines[idx]
            new_coord = f"{new_xyz[0]:8.3f}{new_xyz[1]:8.3f}{new_xyz[2]:8.3f}"
            lines[idx] = ln[:30] + new_coord + ln[54:]
        reoriented += 1

    if reoriented:
        pdb_path.write_text("".join(lines), encoding="utf-8")
    return reoriented


def _write_report_json(
    path: Path,
    result: HydrogenCleanupResult,
    pdb_id: str,
    variant_label: str,
    parameters: dict[str, Any],
) -> None:
    report = {
        "pdb_id": pdb_id,
        "variant_label": variant_label,
        "status": result.status,
        "transition_metal_count": result.transition_metal_count,
        "removed_count": result.removed_count,
        "kept_count": result.kept_count,
        "manual_review_count": result.manual_review_count,
        "input_pdb_path": str(result.input_pdb_path),
        "cleaned_pdb_path": str(result.cleaned_pdb_path),
        "parameters": parameters,
    }
    path.write_text(json.dumps(report, indent=2), encoding="utf-8")
