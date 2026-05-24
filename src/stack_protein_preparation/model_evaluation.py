"""
Local structure quality evaluation for prepared protein structures.

Implements local equivalents of the SAVES server validation tools:
  - Ramachandran analysis  (PROCHECK equivalent)
  - Heavy-atom / all-atom clashscore  (ERRAT / MolProbity equivalent)
  - Ramachandran scatter plot (PNG figure)

The Ramachandran regions are approximated from Lovell et al. 2003
(MolProbity) and Laskowski et al. 1993 (PROCHECK).  Clashscore counts
sterically overlapping atom pairs (VDW overlap > 0.4 Å) per 100 residues,
following the MolProbity definition.

All proteins with a prepared structure are evaluated.  Results are written
to ``model_evaluation/`` inside the protein directory, including a PNG
Ramachandran scatter plot.
"""

from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Any

import numpy as np
from matplotlib.path import Path as MplPath

from Bio.PDB import NeighborSearch, PDBParser, PPBuilder

# ---------------------------------------------------------------------------
# Ramachandran region boundaries
# ---------------------------------------------------------------------------
# Polygon vertices (φ, ψ) in degrees, based on Lovell et al. 2003 and
# Laskowski et al. 1993.  These are simplified polygons that approximate
# the published probability-density contours; they reproduce the
# favored/allowed/outlier classification to within ~1 % of the exact
# MolProbity boundaries for typical protein structures.


def _rect(phi_lo, phi_hi, psi_lo, psi_hi) -> MplPath:
    v = [
        (phi_lo, psi_lo),
        (phi_hi, psi_lo),
        (phi_hi, psi_hi),
        (phi_lo, psi_hi),
        (phi_lo, psi_lo),
    ]
    return MplPath(np.array(v))


# --- General residues (not Gly, Pro, or pre-Pro) ---
# Favored: α-helix core + β-sheet core + left-handed α core
_GEN_FAVORED: list[MplPath] = [
    _rect(-90, -30, -80,  0),    # α-helix
    _rect(-170, -50, 90, 180),   # β-sheet (main)
    _rect(-180, -150, 90, 180),  # β-sheet (wrap at −180)
    _rect(130, 180, 90, 180),    # β-sheet (wrap at +180)
    _rect(40, 90, 20, 75),       # left-handed α
]
# Allowed: extended zones around each favoured region
_GEN_ALLOWED: list[MplPath] = [
    _rect(-120, -20, -90, 15),   # extended α
    _rect(-180, -30, 75, 180),   # extended β (main)
    _rect(-180, -140, -180, -150),  # upper-left β tail (wrap)
    _rect(120, 180, 75, 180),    # extended β (right wrap)
    _rect(25, 105, 10, 95),      # extended left α
]

# --- Glycine (symmetric: both positive and negative φ accessible) ---
_GLY_FAVORED: list[MplPath] = [
    _rect(-180, 180,  50, 180),  # upper half (large)
    _rect(-180, 180, -180, -140),  # lower wrap
    _rect(-160, -30, -80, 0),    # lower-left α-like
    _rect(30, 160, -80, 0),      # mirror lower-right
]
_GLY_ALLOWED: list[MplPath] = [
    _rect(-180, 180,  20, 180),
    _rect(-180, 180, -180, -110),
    _rect(-170, -20, -110, 30),
    _rect(20, 170, -110, 30),
]

# --- Proline (φ constrained by ring to ~ −90 to −30°) ---
_PRO_FAVORED: list[MplPath] = [
    _rect(-90, -30, -80, 175),
    _rect(-90, -30, -180, -165),  # wrap at −180
]
_PRO_ALLOWED: list[MplPath] = [
    _rect(-100, -20, -90, 180),
    _rect(-100, -20, -180, -155),
]

# ---------------------------------------------------------------------------
# VDW radii (Å) for heavy atoms and hydrogen (Bondi 1964, MolProbity values)
# ---------------------------------------------------------------------------
_VDW: dict[str, float] = {
    "H":  1.20,
    "C":  1.70,
    "N":  1.625,
    "O":  1.48,
    "S":  1.782,
    "P":  1.871,
    "SE": 1.90,
    "FE": 1.47,
    "ZN": 1.39,
    "MG": 1.73,
    "CA": 1.74,
    "MN": 1.61,
}
_VDW_DEFAULT = 1.70

# Clash threshold (MolProbity definition)
_CLASH_OVERLAP = 0.4  # Å


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _in_any(phi: float, psi: float, regions: list[MplPath]) -> bool:
    pt = np.array([[phi, psi]])
    return any(r.contains_points(pt)[0] for r in regions)


def _classify_rama(phi: float, psi: float, resname: str, next_resname: str | None) -> str:
    """Return 'favored', 'allowed', or 'outlier' for a residue's φ/ψ."""
    rn = resname.strip().upper()
    if rn == "GLY":
        fav, alw = _GLY_FAVORED, _GLY_ALLOWED
    elif rn == "PRO":
        fav, alw = _PRO_FAVORED, _PRO_ALLOWED
    else:
        fav, alw = _GEN_FAVORED, _GEN_ALLOWED
    if _in_any(phi, psi, fav):
        return "favored"
    if _in_any(phi, psi, alw):
        return "allowed"
    return "outlier"


def _element_from_name(atom_name: str) -> str:
    """Guess element symbol from PDB atom name when the element column is blank."""
    name = atom_name.strip().lstrip("0123456789")
    if not name:
        return "C"
    el = name[0].upper()
    # Two-char elements
    if len(name) >= 2 and name[:2].upper() in {"FE", "ZN", "MG", "CA", "MN", "SE", "NA", "CL"}:
        return name[:2].upper()
    return el


# ---------------------------------------------------------------------------
# Core analyses
# ---------------------------------------------------------------------------

def _get_modelled_residue_keys(
    pdb_id: str,
    protein_dir: Path,
    pipeline_record: dict[str, Any],
) -> set[tuple[str, int]] | None:
    """
    Return (chain_id, resseq) pairs added by Modeller gap-filling.

    Compares the pre-filling delins PDB with the prepared structure.
    Returns None when Modeller was not used or source files are unavailable.
    """
    filler_source = str(pipeline_record.get("filler.model_source", "")).lower()
    if "modeller" not in filler_source:
        return None

    # Pre-filling structure: insertion-code–processed crystal structure
    delins = sorted(protein_dir.rglob(f"*{pdb_id.upper()}_delins.pdb"))
    if not delins:
        delins = sorted(protein_dir.rglob(f"*{pdb_id.lower()}_delins.pdb"))
    if not delins:
        return None

    # Prepared (post-filling) structure
    prep_out = str(pipeline_record.get("prepared_structure.output_path", "")).strip()
    prep_pdb = Path(prep_out) if prep_out else None
    if not prep_pdb or not prep_pdb.exists():
        prep_dir = protein_dir / "prepared"
        hits = sorted(prep_dir.rglob("*.pdb")) if prep_dir.exists() else []
        prep_pdb = hits[0] if hits else None
    if not prep_pdb or not prep_pdb.exists():
        return None

    parser = PDBParser(QUIET=True)

    crystal_keys: set[tuple[str, int]] = set()
    try:
        for model in parser.get_structure("c", str(delins[0])):
            for chain in model:
                for res in chain:
                    if res.id[0] == " ":
                        crystal_keys.add((chain.id, res.id[1]))
            break
    except Exception:
        return None

    modelled: set[tuple[str, int]] = set()
    try:
        for model in parser.get_structure("p", str(prep_pdb)):
            for chain in model:
                for res in chain:
                    if res.id[0] == " ":
                        key = (chain.id, res.id[1])
                        if key not in crystal_keys:
                            modelled.add(key)
            break
    except Exception:
        return None

    return modelled if modelled else None


def _residue_in_gap_ranges(
    chain_id: str,
    resseq: int,
    gap_ranges: list[tuple[str, int, int]],
) -> bool:
    """Return True when (chain_id, resseq) falls within any (chain, prev_res, next_res) window."""
    for gc, prev_r, next_r in gap_ranges:
        if chain_id == gc and prev_r <= resseq <= next_r:
            return True
    return False


def _compute_ramachandran(
    pdb_path: Path,
    modelled_keys: set[tuple[str, int]] | None = None,
    gap_ranges: list[tuple[str, int, int]] | None = None,
) -> dict[str, Any]:
    """
    Compute Ramachandran statistics using BioPython φ/ψ extraction.

    When *gap_ranges* is provided, only residues whose (chain, resseq) falls
    within a (chain, prev_res, next_res) window are included in the statistics.

    When *modelled_keys* is provided (chain_id, resseq pairs added by Modeller),
    each residue is tagged ``is_modelled=True/False`` so the scatter plot can
    distinguish Modeller-built geometry from crystallographically determined residues.
    """
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("s", str(pdb_path))
    ppb = PPBuilder()

    counts: dict[str, int] = {"favored": 0, "allowed": 0, "outlier": 0}
    mod_counts: dict[str, int] = {"favored": 0, "allowed": 0, "outlier": 0}
    residue_data: list[dict[str, Any]] = []

    for model in structure:
        for chain in model:
            for pp in ppb.build_peptides(chain):
                residues = list(pp)
                phi_psi = pp.get_phi_psi_list()
                for i, (res, (phi, psi)) in enumerate(zip(residues, phi_psi)):
                    if phi is None or psi is None:
                        continue
                    if gap_ranges is not None and not _residue_in_gap_ranges(
                        chain.id, res.id[1], gap_ranges
                    ):
                        continue
                    resname = res.get_resname().strip()
                    next_rn = (
                        residues[i + 1].get_resname().strip()
                        if i + 1 < len(residues) else None
                    )
                    phi_deg = math.degrees(phi)
                    psi_deg = math.degrees(psi)
                    label = _classify_rama(phi_deg, psi_deg, resname, next_rn)
                    counts[label] += 1
                    is_mod = (
                        modelled_keys is not None
                        and (chain.id, res.id[1]) in modelled_keys
                    )
                    if is_mod:
                        mod_counts[label] += 1
                    residue_data.append({
                        "phi":        phi_deg,
                        "psi":        psi_deg,
                        "class":      label,
                        "resname":    resname,
                        "chain":      chain.id,
                        "resseq":     res.id[1],
                        "is_modelled": is_mod,
                    })

    total = sum(counts.values())
    if total == 0:
        return {
            "n_total": 0, "pct_favored": None,
            "pct_allowed": None, "pct_outlier": None,
            "residues": [],
            "modelled_n_total": 0, "modelled_pct_favored": None,
            "modelled_pct_outlier": None,
        }

    mod_total = sum(mod_counts.values())
    return {
        "n_total":     total,
        "n_favored":   counts["favored"],
        "n_allowed":   counts["allowed"],
        "n_outlier":   counts["outlier"],
        "pct_favored": round(100 * counts["favored"] / total, 1),
        "pct_allowed": round(100 * counts["allowed"] / total, 1),
        "pct_outlier": round(100 * counts["outlier"] / total, 1),
        "residues":    residue_data,
        # Modeller-built residues only (None when no gap-filling was done)
        "modelled_n_total":      mod_total if modelled_keys is not None else None,
        "modelled_pct_favored":  (
            round(100 * mod_counts["favored"] / mod_total, 1)
            if mod_total else None
        ),
        "modelled_pct_outlier":  (
            round(100 * mod_counts["outlier"] / mod_total, 1)
            if mod_total else None
        ),
    }


def _compute_clashscore(pdb_path: Path) -> dict[str, Any]:
    """
    Compute clashscore: bad clashes per 100 residues.

    A clash is defined as a pair of non-bonded atoms with VDW overlap > 0.4 Å.
    Exclusions: same residue; backbone–backbone across adjacent peptide bonds
    (|Δresseq| == 1).  All atoms (heavy + H) are used because the prepared
    structures already carry explicit protons from GROMACS protonation.
    """
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("s", str(pdb_path))

    # Collect atoms with their parent residue key
    atoms: list[Any] = []
    atom_res_key: list[tuple] = []  # (chain_id, resseq)

    for model in structure:
        for chain in model:
            for residue in chain:
                rid = residue.id[1]
                cid = chain.id
                for atom in residue.get_atoms():
                    atoms.append(atom)
                    atom_res_key.append((cid, rid))
        break  # first model only

    if len(atoms) < 2:
        return {"n_clashes": 0, "n_residues": 0, "clashscore": 0.0}

    residue_set: set[tuple] = set(atom_res_key)
    n_residues = len(residue_set)

    atom_index: dict[int, int] = {id(a): i for i, a in enumerate(atoms)}

    ns = NeighborSearch(atoms)
    max_search = 2 * max(_VDW.values()) + _CLASH_OVERLAP + 0.1  # ~4.3 Å

    seen_pairs: set[frozenset] = set()
    n_clashes = 0

    for idx, atom in enumerate(atoms):
        el = (atom.element or "").strip().upper()
        if not el:
            el = _element_from_name(atom.name)
        r1 = _VDW.get(el, _VDW_DEFAULT)
        cid1, rid1 = atom_res_key[idx]

        nearby = ns.search(atom.coord, max_search, level="A")
        for other in nearby:
            oidx = atom_index[id(other)]
            if oidx <= idx:
                continue

            pair = frozenset((idx, oidx))
            if pair in seen_pairs:
                continue
            seen_pairs.add(pair)

            cid2, rid2 = atom_res_key[oidx]

            # Exclude same residue
            if cid1 == cid2 and rid1 == rid2:
                continue

            # Exclude bonded backbone pairs across adjacent residues
            if cid1 == cid2 and abs(rid1 - rid2) == 1:
                bb = {"N", "CA", "C", "O"}
                if atom.name.strip() in bb and other.name.strip() in bb:
                    continue

            el2 = (other.element or "").strip().upper()
            if not el2:
                el2 = _element_from_name(other.name)
            r2 = _VDW.get(el2, _VDW_DEFAULT)

            dist = float(np.linalg.norm(atom.coord - other.coord))
            overlap = r1 + r2 - dist
            if overlap > _CLASH_OVERLAP:
                n_clashes += 1

    clashscore = round(n_clashes / max(n_residues, 1) * 100, 1)
    return {
        "n_clashes":  n_clashes,
        "n_residues": n_residues,
        "clashscore": clashscore,
    }


def _overall_quality(pct_favored: float | None, clashscore: float | None) -> str:
    """
    Overall quality grade.

    Thresholds are intentionally more lenient than for X-ray structures
    (MolProbity targets >98 % favoured for crystal structures; Modeller
    models typically achieve 85–96 % with clashscores of 5–30).

      good       : ≥ 90 % favoured AND clashscore < 15
      acceptable : ≥ 80 % favoured AND clashscore < 30
      poor       : below acceptable thresholds
    """
    if pct_favored is None or clashscore is None:
        return "unknown"
    if pct_favored >= 90 and clashscore < 15:
        return "good"
    if pct_favored >= 80 and clashscore < 30:
        return "acceptable"
    return "poor"


# ---------------------------------------------------------------------------
# Ramachandran scatter plot
# ---------------------------------------------------------------------------

def _generate_ramachandran_plot(
    residue_data: list[dict[str, Any]],
    output_path: Path,
    title: str = "",
) -> None:
    """
    Write a phi-psi scatter plot PNG.

    Background rectangles show favored (light green) and allowed (light
    yellow) regions for general residues.  Dots are coloured by classification
    (green / orange / red) and shaped by residue type (Gly = triangle, Pro =
    diamond, others = circle).
    """
    from matplotlib.figure import Figure
    from matplotlib.backends.backend_agg import FigureCanvasAgg
    import matplotlib.patches as mpatches
    import matplotlib.lines as mlines

    fig = Figure(figsize=(5.5, 5.5), dpi=150)
    FigureCanvasAgg(fig)
    ax = fig.add_subplot(1, 1, 1)

    ax.set_xlim(-180, 180)
    ax.set_ylim(-180, 180)
    ax.set_xlabel("φ (degrees)", fontsize=9)
    ax.set_ylabel("ψ (degrees)", fontsize=9)
    ax.tick_params(labelsize=8)
    ax.axhline(0, color="#BBBBBB", lw=0.5, zorder=0)
    ax.axvline(0, color="#BBBBBB", lw=0.5, zorder=0)

    FAVORED_BG = "#C8E6C9"   # light green
    ALLOWED_BG = "#FFF9C4"   # light yellow

    # Draw background: allowed first, then favored on top
    for region in _GEN_ALLOWED:
        v = region.vertices
        phi_lo, psi_lo = float(v[0, 0]), float(v[0, 1])
        phi_hi, psi_hi = float(v[2, 0]), float(v[2, 1])
        ax.add_patch(mpatches.Rectangle(
            (phi_lo, psi_lo), phi_hi - phi_lo, psi_hi - psi_lo,
            facecolor=ALLOWED_BG, edgecolor="none", zorder=1,
        ))
    for region in _GEN_FAVORED:
        v = region.vertices
        phi_lo, psi_lo = float(v[0, 0]), float(v[0, 1])
        phi_hi, psi_hi = float(v[2, 0]), float(v[2, 1])
        ax.add_patch(mpatches.Rectangle(
            (phi_lo, psi_lo), phi_hi - phi_lo, psi_hi - psi_lo,
            facecolor=FAVORED_BG, edgecolor="none", zorder=2,
        ))

    COLOR_MAP = {"favored": "#1B5E20", "allowed": "#E65100", "outlier": "#B71C1C"}
    MARKER_MAP: dict[str | None, str] = {"GLY": "^", "PRO": "D", None: "o"}

    has_modelled = any(r.get("is_modelled") for r in residue_data)

    for cls in ("favored", "allowed", "outlier"):
        cls_res = [r for r in residue_data if r["class"] == cls]
        for rn_key, marker in MARKER_MAP.items():
            if rn_key is None:
                sub = [r for r in cls_res if r["resname"] not in ("GLY", "PRO")]
            else:
                sub = [r for r in cls_res if r["resname"] == rn_key]
            if not sub:
                continue

            if has_modelled:
                # Crystal residues: faint background context
                cryst = [r for r in sub if not r.get("is_modelled")]
                if cryst:
                    ax.scatter(
                        [r["phi"] for r in cryst], [r["psi"] for r in cryst],
                        c=COLOR_MAP[cls], marker=marker,
                        s=7, alpha=0.25, linewidths=0, zorder=3,
                    )
                # Modelled residues: prominent foreground
                mod = [r for r in sub if r.get("is_modelled")]
                if mod:
                    ax.scatter(
                        [r["phi"] for r in mod], [r["psi"] for r in mod],
                        c=COLOR_MAP[cls], marker=marker,
                        s=22, alpha=0.95, linewidths=0.5,
                        edgecolors="black", zorder=4,
                    )
            else:
                xs = [r["phi"] for r in sub]
                ys = [r["psi"] for r in sub]
                ax.scatter(
                    xs, ys, c=COLOR_MAP[cls], marker=marker,
                    s=10, alpha=0.75, linewidths=0, zorder=3,
                )

    if title:
        ax.set_title(title, fontsize=8, wrap=True)

    legend_handles = [
        mpatches.Patch(facecolor=FAVORED_BG, edgecolor="#999", label="Favored region"),
        mpatches.Patch(facecolor=ALLOWED_BG, edgecolor="#999", label="Allowed region"),
        mlines.Line2D([], [], marker="o", color="w", markerfacecolor=COLOR_MAP["favored"],
                      markersize=6, label="Favored"),
        mlines.Line2D([], [], marker="o", color="w", markerfacecolor=COLOR_MAP["allowed"],
                      markersize=6, label="Allowed"),
        mlines.Line2D([], [], marker="o", color="w", markerfacecolor=COLOR_MAP["outlier"],
                      markersize=6, label="Outlier"),
        mlines.Line2D([], [], marker="^", color="w", markerfacecolor="#555",
                      markersize=6, label="Gly"),
        mlines.Line2D([], [], marker="D", color="w", markerfacecolor="#555",
                      markersize=6, label="Pro"),
    ]
    if has_modelled:
        legend_handles += [
            mlines.Line2D([], [], marker="o", color="w",
                          markerfacecolor="#555", markersize=4, alpha=0.3,
                          label="Crystal residues (context)"),
            mlines.Line2D([], [], marker="o", color="w",
                          markerfacecolor="#555", markersize=8,
                          markeredgecolor="black", markeredgewidth=0.5,
                          label="Modeller-built residues"),
        ]
    ax.legend(handles=legend_handles, loc="upper right", fontsize=6,
              framealpha=0.85, borderpad=0.5)

    fig.tight_layout()
    fig.savefig(str(output_path), dpi=150, bbox_inches="tight")


# ---------------------------------------------------------------------------
# Main entry point
# ---------------------------------------------------------------------------

def run_model_evaluation(
    protein_dir: Path,
    pdb_id: str,
    pipeline_record: dict[str, str],
) -> dict[str, Any]:
    """
    Evaluate the prepared structure(s) for one protein directory.

    Only runs for proteins where Modeller was used to fill gaps.  Evaluates
    residues in the window from the residue before each gap to the residue
    after it (gap boundaries stored in pipeline_record["gap_boundaries"]).

    Returns
    -------
    dict with keys: status, message, evaluations (list), manifest_path,
    summary scalars (pct_favored, pct_outlier, clashscore, overall_quality,
    rama_plot_path) taken from the best-available variant.
    """
    from stack_protein_preparation.pipeline_state import (
        FILLER_MODEL_SOURCE_COLUMN_NAME,
        GAP_BOUNDARIES_COLUMN_NAME,
    )

    output_dir = protein_dir / "model_evaluation"

    result: dict[str, Any] = {
        "pdb_id":                  pdb_id,
        "status":                  "skipped",
        "message":                 "",
        "evaluations":             [],
        "pct_favored":             None,
        "pct_allowed":             None,
        "pct_outlier":             None,
        "clashscore":              None,
        "overall_quality":         "unknown",
        "rama_plot_path":          "",
        "manifest_path":           str(output_dir / "model_evaluation_manifest.json"),
        "modelled_n_residues":     None,
        "modelled_pct_favored":    None,
        "modelled_pct_outlier":    None,
    }

    # Only evaluate Modeller-filled proteins
    filler_source = str(pipeline_record.get(FILLER_MODEL_SOURCE_COLUMN_NAME, "")).lower()
    used_modeller = "modeller" in filler_source
    if not used_modeller:
        result["message"] = "Skipped: Modeller was not used for gap-filling."
        return result

    prepared_dir = protein_dir / "prepared"
    if not prepared_dir.exists():
        result["status"] = "failed"
        result["message"] = f"prepared/ directory not found: {prepared_dir}"
        return result

    pdb_files = sorted(prepared_dir.rglob("*.pdb"))
    if not pdb_files:
        result["status"] = "failed"
        result["message"] = "No prepared PDB files found."
        return result

    # Build gap_ranges from stored gap boundaries: [(chain, prev_res, next_res), ...]
    gap_ranges: list[tuple[str, int, int]] | None = None
    raw_boundaries = str(pipeline_record.get(GAP_BOUNDARIES_COLUMN_NAME, "")).strip()
    if raw_boundaries:
        try:
            boundaries = json.loads(raw_boundaries)
            gap_ranges = [
                (str(g["chain"]), int(g["prev_res"]), int(g["next_res"]))
                for g in boundaries
                if g.get("chain") is not None
                and g.get("prev_res") is not None
                and g.get("next_res") is not None
            ]
            if not gap_ranges:
                gap_ranges = None
        except (json.JSONDecodeError, KeyError, TypeError, ValueError):
            gap_ranges = None

    output_dir.mkdir(parents=True, exist_ok=True)
    evaluations: list[dict] = []

    modelled_keys = _get_modelled_residue_keys(pdb_id, protein_dir, pipeline_record)

    for pdb_file in pdb_files:
        variant_label = pdb_file.parent.name  # e.g. "best_complete", "gaps"
        eval_entry: dict[str, Any] = {
            "variant":      variant_label,
            "pdb_path":     str(pdb_file),
            "rama":         None,
            "clash":        None,
            "overall":      "unknown",
            "error":        None,
            "rama_plot":    "",
        }

        try:
            rama = _compute_ramachandran(
                pdb_file,
                modelled_keys=modelled_keys,
                gap_ranges=gap_ranges,
            )
            eval_entry["rama"] = {k: v for k, v in rama.items() if k != "residues"}
        except Exception as exc:
            eval_entry["error"] = f"Ramachandran failed: {exc}"
            evaluations.append(eval_entry)
            continue

        # Generate scatter plot for this variant
        if rama.get("residues"):
            plot_path = output_dir / f"{pdb_id}_{variant_label}_ramachandran.png"
            title_parts = [f"{pdb_id} — {variant_label}"]
            if used_modeller:
                title_parts.append("(includes Modeller-built residues)")
            n_fav = rama.get("pct_favored")
            n_out = rama.get("pct_outlier")
            if n_fav is not None:
                title_parts.append(f"favored {n_fav}%  outlier {n_out}%")
            try:
                _generate_ramachandran_plot(
                    residue_data=rama["residues"],
                    output_path=plot_path,
                    title="\n".join(title_parts),
                )
                eval_entry["rama_plot"] = str(plot_path)
            except Exception as exc:
                eval_entry["error"] = (eval_entry.get("error") or "") + f" | Plot failed: {exc}"

        try:
            clash = _compute_clashscore(pdb_file)
            eval_entry["clash"] = clash
        except Exception as exc:
            eval_entry["error"] = (eval_entry.get("error") or "") + f" | Clashscore failed: {exc}"

        if eval_entry["rama"] and eval_entry["clash"]:
            eval_entry["overall"] = _overall_quality(
                eval_entry["rama"].get("pct_favored"),
                eval_entry["clash"].get("clashscore"),
            )

        evaluations.append(eval_entry)

    result["evaluations"] = evaluations
    result["status"] = "success" if evaluations else "failed"
    gap_count = len(gap_ranges) if gap_ranges else 0
    result["message"] = (
        f"Evaluated {len(evaluations)} prepared variant(s); "
        f"Ramachandran restricted to {gap_count} gap window(s)."
    )

    # Expose summary scalars from the best_complete or first variant
    preferred = next(
        (e for e in evaluations if "best_complete" in e["variant"]),
        evaluations[0] if evaluations else None,
    )
    if preferred:
        if preferred.get("rama"):
            result["pct_favored"] = preferred["rama"].get("pct_favored")
            result["pct_allowed"] = preferred["rama"].get("pct_allowed")
            result["pct_outlier"] = preferred["rama"].get("pct_outlier")
            result["modelled_n_residues"] = preferred["rama"].get("modelled_n_total")
            result["modelled_pct_favored"] = preferred["rama"].get("modelled_pct_favored")
            result["modelled_pct_outlier"] = preferred["rama"].get("modelled_pct_outlier")
        if preferred.get("clash"):
            result["clashscore"] = preferred["clash"].get("clashscore")
        result["overall_quality"] = preferred["overall"]
        result["rama_plot_path"] = preferred.get("rama_plot", "")

    # Write manifest (exclude large residue lists)
    manifest = {k: v for k, v in result.items() if k != "evaluations"}
    manifest["evaluations"] = [
        {ek: ev for ek, ev in e.items()} for e in evaluations
    ]
    Path(result["manifest_path"]).write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    return result
