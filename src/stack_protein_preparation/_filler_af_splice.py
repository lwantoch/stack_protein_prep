"""Splice AlphaFold-modelled loops into a crystal template.

Two-layer fix (memory notes ``project_fruton_af_wholebody_bug`` and
``feedback_fruton_gap_continuity_rule``):

1. The whole-body-replace bug: FRUTON Step 9 previously returned the ENTIRE
   aligned AF model as ``final_filled_model.pdb``.  Chain B, ligands,
   cofactors, waters were silently dropped for 13 MMBSA_200 PDBs.  This
   module keeps every residue present in the crystal template and only
   inserts AF residues that fill missing-residue windows.

2. The backbone-discontinuity bug: AF is aligned globally (CA superimposition
   over all shared residues), so AF-loop endpoints float away from the
   crystal residues flanking each gap.  The user observed C-N distances of
   2.4-14.5 A at boundaries (should be ~1.33 A).  Fix: for each missing-
   residue window we re-fit the AF loop locally using K flanking anchor
   residues that are present in BOTH crystal and AF, then apply the
   resulting rigid-body transform to the loop residues before insertion.

HETATMs, waters, ions live in the crystal template only and are preserved.
"""
from __future__ import annotations

from pathlib import Path

from Bio.PDB import PDBIO, PDBParser
from Bio.PDB.Chain import Chain as BioChain
from Bio.PDB.Superimposer import Superimposer


_PARSER = PDBParser(QUIET=True)

ANCHOR_ATOMS: tuple[str, ...] = ("N", "CA", "C")
# Deliberately excludes the carbonyl O: 3 non-colinear atoms exactly
# determine a 6-DOF rigid body, so Superimposer produces RMSD=0 and the
# anchor's C atom lands *exactly* on the crystal's C atom.  The loop's
# terminal N (bonded to that C in AF at ~1.33 A) then stitches to the
# crystal at true peptide-bond geometry.  With 4 atoms the fit is
# over-determined; residuals of ~0.5 A on C propagate directly into
# boundary-bond errors.
# Number of flanking crystal residues (per side) to use as anchors for the
# local rigid-body fit.
#
# K=1 uniformly: with 3 atoms per residue (see ANCHOR_ATOMS below), a single
# anchor residue gives 3 atoms = 6 coords, which exactly determines the
# 6-DOF rigid body.  Superimposer produces RMSD~0 and the anchor residue's
# atoms overlay their crystal partners exactly, so the loop's terminal N
# lands at true peptide-bond distance from crystal C(first-1).
# For K>=2 the fit becomes over-determined and drifts (RMSD ~0.5-1 A), which
# propagates directly into 1.5-4 A boundary-bond errors.
# For very short chains the AF model is already globally aligned to the
# crystal frame, so no extra orientation constraint is needed.
ANCHOR_K_INTERNAL = 1
ANCHOR_K_TERMINAL = 1
# How far to search past the gap for viable anchor residues (must be in both
# crystal and AF, with all ANCHOR_ATOMS present).  Prevents endless scan on
# unusual chains.
ANCHOR_SEARCH_MAX = 30


def _protein_residue_map(structure) -> dict[str, dict[int, object]]:
    """{chain_id: {resi: residue}} for protein residues (het-flag empty)."""
    out: dict[str, dict[int, object]] = {}
    for model in structure:
        for chain in model:
            per_chain: dict[int, object] = {}
            for res in chain:
                hetflag, resi, _icode = res.id
                if hetflag.strip():
                    continue
                per_chain[resi] = res
            out[chain.id] = per_chain
    return out


def _detect_missing_windows(
    crystal_resnums: set[int], af_resnums_sorted: list[int]
) -> list[tuple[int, int]]:
    """Contiguous (first, last) windows present in AF but missing in crystal."""
    windows: list[tuple[int, int]] = []
    start: int | None = None
    prev: int | None = None
    for r in af_resnums_sorted:
        missing = r not in crystal_resnums
        if missing and start is None:
            start = r
        elif not missing and start is not None:
            windows.append((start, prev))  # type: ignore[arg-type]
            start = None
        prev = r
    if start is not None:
        windows.append((start, prev))  # type: ignore[arg-type]
    return windows


def _collect_anchors(
    crystal_by_resi: dict[int, object],
    af_by_resi: dict[int, object],
    start: int,
    step: int,
    stop: int,
    max_picks: int,
) -> list[tuple[object, object]]:
    """Walk resi from ``start`` toward ``stop`` in ``step`` increments and
    collect up to ``max_picks`` (crystal_res, af_res) pairs where both have
    all ANCHOR_ATOMS."""
    picks: list[tuple[object, object]] = []
    r = start
    while (step > 0 and r <= stop) or (step < 0 and r >= stop):
        if len(picks) >= max_picks:
            break
        cres = crystal_by_resi.get(r)
        ares = af_by_resi.get(r)
        if (
            cres is not None
            and ares is not None
            and all(a in cres for a in ANCHOR_ATOMS)
            and all(a in ares for a in ANCHOR_ATOMS)
        ):
            picks.append((cres, ares))
        r += step
    return picks


def _fit_with_pairs(pairs, af_by_resi, resi_range):
    """Superimpose AF backbone atoms of ``pairs`` onto their crystal partners,
    then apply the transform to every atom of each residue in ``resi_range``.

    Returns list of transformed residue copies (may be shorter than the range
    if AF happens to lack a resi).
    """
    if not pairs:
        return None, None
    fixed_atoms: list[object] = []
    mobile_atoms: list[object] = []
    for cres, ares in pairs:
        for a in ANCHOR_ATOMS:
            fixed_atoms.append(cres[a])
            mobile_atoms.append(ares[a])
    if len(fixed_atoms) < 3:  # rigid body needs 3 non-colinear atoms
        return None, None
    sup = Superimposer()
    sup.set_atoms(fixed_atoms, mobile_atoms)

    fitted: list[object] = []
    for r in resi_range:
        ares = af_by_resi.get(r)
        if ares is None:
            continue
        rcopy = ares.copy()
        rcopy.detach_parent()
        sup.apply(list(rcopy.get_atoms()))
        fitted.append(rcopy)
    return fitted, float(sup.rms)


def _fit_af_loop(
    crystal_by_resi: dict[int, object],
    af_by_resi: dict[int, object],
    window: tuple[int, int],
) -> tuple[list[object], dict[str, float]] | None:
    """Return (fitted_residues, diagnostics) or None if not fittable.

    Strategy: fit each half of the loop to its own side's anchors.
    * The N-terminal half is rigid-body-fit using left-side crystal anchors,
      so its N-terminal N atom stitches to crystal C(first-1).
    * The C-terminal half is rigid-body-fit using right-side crystal anchors,
      so its C-terminal C atom stitches to crystal N(last+1).
    * A residual break falls somewhere in the middle of the loop -- always
      the most flexible / least constrained part of a loop.  Better than the
      alternative (both endpoints wrong: 4-15 A gaps observed on 1JCN).
    * Terminal gaps (only one side has anchors) fall back to a single fit.
    """
    first, last = window

    # Peek: is this a terminal gap? (only one side has *any* anchors)
    _peek_left = _collect_anchors(
        crystal_by_resi, af_by_resi,
        start=first - 1, step=-1,
        stop=first - 1 - ANCHOR_SEARCH_MAX,
        max_picks=1,
    )
    _peek_right = _collect_anchors(
        crystal_by_resi, af_by_resi,
        start=last + 1, step=+1,
        stop=last + 1 + ANCHOR_SEARCH_MAX,
        max_picks=1,
    )
    if not _peek_left and not _peek_right:
        return None
    is_terminal = not _peek_left or not _peek_right
    k = ANCHOR_K_TERMINAL if is_terminal else ANCHOR_K_INTERNAL

    left = _collect_anchors(
        crystal_by_resi, af_by_resi,
        start=first - 1, step=-1,
        stop=first - 1 - ANCHOR_SEARCH_MAX,
        max_picks=k,
    )
    right = _collect_anchors(
        crystal_by_resi, af_by_resi,
        start=last + 1, step=+1,
        stop=last + 1 + ANCHOR_SEARCH_MAX,
        max_picks=k,
    )

    diagnostics: dict[str, float] = {
        "n_anchor_residues_left": float(len(left)),
        "n_anchor_residues_right": float(len(right)),
        "n_loop_residues": float(last - first + 1),
    }

    # Terminal gap (N- or C-term) -- single fit, whole loop
    if not left or not right:
        pairs = left or right
        fitted, rms = _fit_with_pairs(pairs, af_by_resi, range(first, last + 1))
        if fitted is None:
            return None
        diagnostics["rmsd_single"] = rms
        return fitted, diagnostics

    # Single-residue gap: split-fit would give first==mid==last so right half
    # is empty and only the left boundary stitches.  Fall back to a combined
    # both-side fit so residuals are shared between the two boundaries.
    if first == last:
        combined = left + right
        fitted, rms = _fit_with_pairs(combined, af_by_resi, range(first, last + 1))
        if fitted is None:
            return None
        diagnostics["rmsd_combined"] = rms
        return fitted, diagnostics

    # Internal gap -- split fit
    mid = (first + last) // 2  # first half: first..mid, second half: mid+1..last
    left_fitted, rms_l = _fit_with_pairs(left, af_by_resi, range(first, mid + 1))
    right_fitted, rms_r = _fit_with_pairs(right, af_by_resi, range(mid + 1, last + 1))

    if left_fitted is None and right_fitted is None:
        return None
    fitted = (left_fitted or []) + (right_fitted or [])
    if rms_l is not None:
        diagnostics["rmsd_left"] = rms_l
    if rms_r is not None:
        diagnostics["rmsd_right"] = rms_r
    return fitted, diagnostics


def splice_af_gaps_into_crystal(
    crystal_pdb_path: str | Path,
    af_aligned_pdb_path: str | Path,
    output_pdb_path: str | Path,
    min_gap_plddt: float | None = 50.0,
    enable_rollback: bool = True,
    reject_uniprot_idr_gaps_for: str | None = None,
) -> Path:
    """Splice AF-modelled residues into the missing-residue windows of a
    crystal template with per-gap local anchor alignment.

    Parameters
    ----------
    crystal_pdb_path
        Crystal-template PDB (single- or multi-chain protein file).
    af_aligned_pdb_path
        AF model already globally superimposed onto the crystal frame.
        The global alignment is used only as a starting point -- each gap
        is re-fit locally using its flanking crystal residues.
    output_pdb_path
        Destination for the spliced PDB.
    min_gap_plddt
        Minimum mean per-residue pLDDT (B-factor column of AF model, range
        [0, 100]) required to accept a gap fill.  Windows below this
        threshold are skipped -- AF confidence < 50 corresponds to Croll
        2025's "barbed-wire" regions where the local backbone is worse
        than no model at all.  Pass ``None`` to disable this filter.
    enable_rollback
        When True (default), drop gap fills whose fitted peptide bonds
        fall outside [1.28, 1.40] A or introduce > 5 heavy-atom clashes
        with the crystal environment.  Set to False when a downstream
        loop refiner (e.g. MODELLER LoopModel) will polish the fills --
        rolling back too early denies the refiner useful starting points.
    reject_uniprot_idr_gaps_for
        Optional UniProt accession.  When provided, MobiDB is queried for
        annotated intrinsically-disordered regions and gap windows that
        overlap an IDR (>= 50 % of gap length) are skipped.  Rationale:
        AF3-class models hallucinate confident conformations for genuine
        IDRs (Wang et al. arXiv 2510.15939 2025); a crystal will never
        resolve one, so shipping the fill would mislead a reviewer.
        Fail-open when MobiDB is unreachable (network unavailable in the
        pixi env for local runs).

    Returns
    -------
    Path to the written spliced PDB.
    """
    crystal = _PARSER.get_structure("crystal", str(crystal_pdb_path))
    af = _PARSER.get_structure("af", str(af_aligned_pdb_path))

    crystal_map = _protein_residue_map(crystal)
    af_map = _protein_residue_map(af)

    for cid, af_by_resi in af_map.items():
        af_resnums_sorted = sorted(af_by_resi.keys())

        if cid not in crystal[0]:
            # Whole chain missing in crystal.  Two failure modes seen on
            # newbench_27 bench:
            #   1. 7UL2: crystal = R,D; AF = A.  AF resnums 60-379 overlap
            #      crystal R:60-...  Copying AF as new chain A duplicates
            #      the entire structure and produces 5000+ clashes.
            #   2. 8Q68: crystal = A,B; AF = A.  AF's own chain A is
            #      already used in crystal.  Actually this hits the
            #      `cid in crystal` branch below, so it's the reverse
            #      variant of #1 -- caught by the standard gap loop.
            #
            # Safety: only copy the AF chain wholesale if its residue
            # numbers do NOT overlap ANY crystal chain's residue numbers.
            # Otherwise skip (leave those AF residues out -- the crystal
            # has structural evidence for those positions in another chain).
            _af_resnum_set = set(af_by_resi.keys())
            _overlaps_any_crystal_chain = False
            for _other_cid, _other_map in crystal_map.items():
                if _af_resnum_set & set(_other_map.keys()):
                    _overlaps_any_crystal_chain = True
                    break
            if _overlaps_any_crystal_chain:
                continue
            new_chain = af[0][cid].copy()
            new_chain.detach_parent()
            crystal[0].add(new_chain)
            continue

        crystal_chain = crystal[0][cid]
        crystal_by_resi = crystal_map.get(cid, {})
        windows = _detect_missing_windows(set(crystal_by_resi.keys()), af_resnums_sorted)

        for window in windows:
            # UniProt IDR gate (2026-08-22): reject fills whose window
            # overlaps an annotated intrinsically-disordered region for the
            # parent UniProt entry (MobiDB predictions).  AF3-class models
            # hallucinate confident IDR conformations; a crystal will never
            # resolve one.  Fail-open on network failure.
            if reject_uniprot_idr_gaps_for is not None:
                lo, hi = window
                try:
                    from stack_protein_preparation._uniprot_idr import (
                        gap_overlaps_uniprot_idr,
                    )
                    _is_idr = gap_overlaps_uniprot_idr(
                        reject_uniprot_idr_gaps_for, lo, hi
                    )
                except Exception:  # noqa: BLE001 -- fail-open
                    _is_idr = None
                if _is_idr is True:
                    continue

            # pLDDT gate (2026-08-22): reject AF fills whose per-window mean
            # pLDDT falls below ``min_gap_plddt``.  Croll IUCr Acta D 2025
            # shows AF regions with pLDDT < 50 are "barbed wire" — worse
            # than no fill because they poison downstream refinement.
            if min_gap_plddt is not None:
                lo, hi = window
                _gap_bfactors: list[float] = []
                for _resnum in range(lo, hi + 1):
                    _r = af_by_resi.get(_resnum)
                    if _r is None:
                        continue
                    for _a in _r:
                        _b = float(_a.get_bfactor() or 0.0)
                        if 0.0 <= _b <= 100.0:  # sanity: AF stores pLDDT here
                            _gap_bfactors.append(_b)
                if _gap_bfactors:
                    _mean_plddt = sum(_gap_bfactors) / len(_gap_bfactors)
                    if _mean_plddt < min_gap_plddt:
                        # Low-confidence AF fill; leave the gap as-crystal
                        # (declared missing via REMARK 465 upstream).
                        continue

            result = _fit_af_loop(crystal_by_resi, af_by_resi, window)
            if result is None:
                # Not enough anchors to fit safely; skip this window rather
                # than corrupt the backbone with a bad rigid-body fit.
                continue
            fitted, _diag = result

            # SMART-ROLLBACK (2026-08-22): validate the fitted loop's own
            # internal peptide-bond geometry BEFORE committing to insert it.
            # A fit can be geometrically clean at both anchors but have a
            # mid-loop break where AF-loop cartesian length ≠ crystal-gap
            # distance.  Such breaks are unrecoverable by local minimisation
            # (see project_fruton_af_wholebody_bug memory).  If we detect
            # one, DROP this whole gap-fill -- reviewer standard prefers a
            # REMARK 465 missing residue to a broken peptide bond, and the
            # downstream FRUTON stages can handle unfilled gaps.
            import math as _math
            _fitted_by_resi = {r.id[1]: r for r in fitted}
            _fitted_resnums = sorted(_fitted_by_resi.keys())
            _has_break = False
            for _i in range(1, len(_fitted_resnums)):
                _r_prev = _fitted_by_resi[_fitted_resnums[_i - 1]]
                _r_curr = _fitted_by_resi[_fitted_resnums[_i]]
                if "C" not in _r_prev or "N" not in _r_curr:
                    continue
                if (_fitted_resnums[_i] - _fitted_resnums[_i - 1]) != 1:
                    continue  # residue-number gap inside the fit -- allow
                _d = _r_prev["C"].coord - _r_curr["N"].coord
                _dist = _math.sqrt(float(_d[0]) ** 2 + float(_d[1]) ** 2 + float(_d[2]) ** 2)
                if _dist < 1.28 or _dist > 1.40:
                    _has_break = True
                    break
            # Also check the boundary peptide bonds to existing crystal
            # residues (mid-loop breaks disguised as boundary breaks).
            for _resnum, _res in _fitted_by_resi.items():
                if (_resnum - 1) in crystal_by_resi and "N" in _res:
                    _crystal_prev = crystal_by_resi[_resnum - 1]
                    if "C" in _crystal_prev:
                        _d = _crystal_prev["C"].coord - _res["N"].coord
                        _dist = _math.sqrt(float(_d[0]) ** 2 + float(_d[1]) ** 2 + float(_d[2]) ** 2)
                        if _dist < 1.28 or _dist > 1.40:
                            _has_break = True
                            break
                if (_resnum + 1) in crystal_by_resi and "C" in _res:
                    _crystal_next = crystal_by_resi[_resnum + 1]
                    if "N" in _crystal_next:
                        _d = _res["C"].coord - _crystal_next["N"].coord
                        _dist = _math.sqrt(float(_d[0]) ** 2 + float(_d[1]) ** 2 + float(_d[2]) ** 2)
                        if _dist < 1.28 or _dist > 1.40:
                            _has_break = True
                            break
            if _has_break and enable_rollback:
                # Skip this window; the gap remains as-crystal (missing
                # residues), which is what REMARK 465 already declares.
                continue

            # CLASH-CHECK (2026-08-22): reject the fill if inserting it
            # would produce > _MAX_NEW_CLASHES_PER_GAP heavy-atom overlaps
            # against the existing crystal environment.  164-residue fills
            # can carry ~1000 steric conflicts (5HJS live test) that would
            # dominate the quality gate downstream; a reviewer would rather
            # see REMARK 465 than clashed side chains.
            _MAX_NEW_CLASHES_PER_GAP = 5
            _CLASH_A = 2.0
            _existing_heavy: list[tuple[str, str, int, str, tuple]] = []
            for _c in crystal[0]:
                for _r in _c:
                    if _r.id[0].strip():  # HETATM excluded
                        continue
                    for _a in _r:
                        if _a.element == "H":
                            continue
                        _existing_heavy.append((
                            _c.id, _r.resname, _r.id[1], _a.name,
                            (float(_a.coord[0]), float(_a.coord[1]), float(_a.coord[2])),
                        ))
            _new_clash_count = 0
            for _fr in fitted:
                for _fa in _fr:
                    if _fa.element == "H":
                        continue
                    _fx, _fy, _fz = float(_fa.coord[0]), float(_fa.coord[1]), float(_fa.coord[2])
                    for _ec_id, _er_name, _er_num, _ea_name, (_ex, _ey, _ez) in _existing_heavy:
                        # Skip same residue / adjacent-residue peptide neighbors
                        if _er_num == _fr.id[1] and _ec_id == cid:
                            continue
                        if _ec_id == cid and abs(_er_num - _fr.id[1]) <= 1:
                            continue
                        _dx = _fx - _ex; _dy = _fy - _ey; _dz = _fz - _ez
                        if _dx * _dx + _dy * _dy + _dz * _dz < _CLASH_A * _CLASH_A:
                            _new_clash_count += 1
                            if _new_clash_count > _MAX_NEW_CLASHES_PER_GAP:
                                break
                    if _new_clash_count > _MAX_NEW_CLASHES_PER_GAP:
                        break
                if _new_clash_count > _MAX_NEW_CLASHES_PER_GAP:
                    break
            if _new_clash_count > _MAX_NEW_CLASHES_PER_GAP and enable_rollback:
                continue  # rollback

            for r in fitted:
                if r.id in [c.id for c in crystal_chain]:
                    crystal_chain.detach_child(r.id)
                crystal_chain.add(r)

    # Bio.PDB Chain.add() appends without sorting; rebuild each chain so
    # residues are ordered (protein first by resi, then HETATM by resi).
    for model in crystal:
        for chain in list(model):
            sorted_residues = sorted(
                list(chain),
                key=lambda r: (r.id[0].strip() != "", r.id[1], r.id[2] or " "),
            )
            new_chain = BioChain(chain.id)
            for res in sorted_residues:
                res.detach_parent()
                new_chain.add(res)
            model.detach_child(chain.id)
            model.add(new_chain)

    out = Path(output_pdb_path)
    out.parent.mkdir(parents=True, exist_ok=True)
    io = PDBIO()
    io.set_structure(crystal)
    io.save(str(out))
    return out


def rollback_bad_gap_fills(
    input_pdb_path: str | Path,
    output_pdb_path: str | Path,
    gap_ranges_by_chain: list[tuple[str, int, int]],
    max_broken_bonds_per_gap: int = 0,
    max_clashes_per_gap: int = 5,
    peptide_bond_min_angstrom: float = 1.28,
    peptide_bond_max_angstrom: float = 1.40,
    clash_distance_angstrom: float = 2.0,
) -> tuple[Path, list[tuple[str, int, int]]]:
    """Post-refinement rollback: remove residues in gap ranges that STILL
    fail the peptide-bond or clash gate after a downstream refiner (e.g.
    MODELLER LoopModel) had a chance to fix them.

    ``gap_ranges_by_chain``: list of (chain_id, first_resnum, last_resnum)
    identifying regions that were AF-inserted.  For each range, we compute
    the same broken-bond + heavy-atom-clash metrics used at splice-time,
    and if the region still fails we detach those residues from the PDB
    (they revert to REMARK 465 missing residues, which the FRUTON reporter
    already flags for the reviewer).

    Returns ``(output_pdb_path, list_of_rolled_back_gaps)``.
    """
    import math
    input_path = Path(input_pdb_path)
    output_path = Path(output_pdb_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    structure = _PARSER.get_structure("m", str(input_path))
    # Build lookup: (chain_id, resnum) -> residue
    resmap: dict[tuple[str, int], object] = {}
    for chain in structure[0]:
        for r in chain:
            if r.id[0].strip():
                continue
            resmap[(chain.id, r.id[1])] = r

    # For clash detection across all heavy atoms outside the gap:
    all_heavy: list[tuple[str, int, str, tuple[float, float, float]]] = []
    for chain in structure[0]:
        for r in chain:
            if r.id[0].strip():
                continue
            for a in r:
                if a.element == "H":
                    continue
                all_heavy.append((
                    chain.id, r.id[1], a.name,
                    (float(a.coord[0]), float(a.coord[1]), float(a.coord[2])),
                ))

    rolled_back: list[tuple[str, int, int]] = []
    residues_to_drop: list[tuple[str, tuple]] = []

    for chain_id, lo, hi in gap_ranges_by_chain:
        # Count broken bonds inside the gap and at its boundaries
        broken_count = 0
        for resnum in range(lo, hi + 1):
            r_curr = resmap.get((chain_id, resnum))
            if r_curr is None:
                continue
            for neigh_num in (resnum - 1, resnum + 1):
                r_other = resmap.get((chain_id, neigh_num))
                if r_other is None:
                    continue
                # Direction: prev -> curr uses C(prev) - N(curr)
                if neigh_num == resnum - 1 and "C" in r_other and "N" in r_curr:
                    d = r_other["C"].coord - r_curr["N"].coord
                    dist = math.sqrt(float(d[0]) ** 2 + float(d[1]) ** 2 + float(d[2]) ** 2)
                    if dist < peptide_bond_min_angstrom or dist > peptide_bond_max_angstrom:
                        broken_count += 1
                elif neigh_num == resnum + 1 and "C" in r_curr and "N" in r_other:
                    d = r_curr["C"].coord - r_other["N"].coord
                    dist = math.sqrt(float(d[0]) ** 2 + float(d[1]) ** 2 + float(d[2]) ** 2)
                    if dist < peptide_bond_min_angstrom or dist > peptide_bond_max_angstrom:
                        broken_count += 1

        # Count clashes: any atom in gap residues within clash_distance of
        # a heavy atom in a non-adjacent residue.
        clash_count = 0
        for resnum in range(lo, hi + 1):
            r_curr = resmap.get((chain_id, resnum))
            if r_curr is None:
                continue
            for a in r_curr:
                if a.element == "H":
                    continue
                ax, ay, az = float(a.coord[0]), float(a.coord[1]), float(a.coord[2])
                for ec_id, er_num, _ea_name, (ex, ey, ez) in all_heavy:
                    if ec_id == chain_id and er_num == resnum:
                        continue
                    if ec_id == chain_id and abs(er_num - resnum) <= 2:
                        continue
                    dx = ax - ex; dy = ay - ey; dz = az - ez
                    if dx * dx + dy * dy + dz * dz < clash_distance_angstrom ** 2:
                        clash_count += 1
                        if clash_count > max_clashes_per_gap:
                            break
                if clash_count > max_clashes_per_gap:
                    break
            if clash_count > max_clashes_per_gap:
                break

        if broken_count > max_broken_bonds_per_gap or clash_count > max_clashes_per_gap:
            rolled_back.append((chain_id, lo, hi))
            for resnum in range(lo, hi + 1):
                r = resmap.get((chain_id, resnum))
                if r is not None:
                    residues_to_drop.append((chain_id, r.id))

    # 9E3M artefact fix (2026-08-22): when there is nothing to detach we
    # pass the input through byte-for-byte instead of rewriting via Bio.PDB
    # PDBIO.  The rewrite path shifts atom serials, TER records, and the
    # element column enough that downstream NeighborSearch flags 33 phantom
    # clashes on unchanged residues (documented in
    # docs/FRUTON_hardening_2026-08.md, section 4).  Bypass avoids the
    # artefact entirely and preserves the caller's file formatting.
    if not residues_to_drop:
        if str(output_path) != str(input_path):
            output_path.write_bytes(input_path.read_bytes())
        return output_path, rolled_back

    for chain_id, rid in residues_to_drop:
        try:
            structure[0][chain_id].detach_child(rid)
        except KeyError:
            pass

    io = PDBIO()
    io.set_structure(structure)
    io.save(str(output_path))
    return output_path, rolled_back
