"""Per-loop plDDT extraction for AlphaFold-filled models.

FRUTON's Step 9 filler can fill missing residues via MODELLER or AlphaFold.
For AlphaFold fills, the AF model stores per-residue plDDT confidence in the
B-factor column (values in [0, 100]).  This module extracts the mean plDDT per
filled gap segment so downstream stages (docking scoring, MD prep, reviewer
dashboards) can weight confidence.

For MODELLER fills there is no plDDT in the output (MODELLER writes its own
scoring in a separate log).  This module returns `None` for such loops.

Public entry point:
    compute_per_loop_plddt(af_model_path, gap_ranges) -> list[dict]

Where each result dict has:
    chain          : chain letter
    start          : first residue in the loop (int)
    end            : last residue in the loop (int)
    method         : "alphafold" | "modeller" | "unknown"
    mean_plddt     : float in [0, 100] or None
    mean_bfactor   : float (always populated when atoms present)
    n_atoms        : int (number of atoms averaged)
    plddt_note     : optional warning when mean_bfactor is outside [0, 100]
                     (would indicate the B-factor column carries crystal
                     temperature factors, not plDDT — happens when the AF
                     model was merged into a crystal frame)
"""
from __future__ import annotations

from pathlib import Path
from typing import Any

from Bio.PDB import PDBParser


_PARSER = PDBParser(QUIET=True)


def _extract_bfactors_in_range(
    pdb_path: Path,
    chain: str,
    start: int,
    end: int,
) -> list[float]:
    """Return CA B-factors for each residue in [start,end] on `chain`.
    Falls back to all-atom B-factors if no CA present."""
    try:
        struct = _PARSER.get_structure("s", str(pdb_path))
    except Exception:
        return []
    ca, any_atom = [], []
    for model in struct:
        if chain not in model:
            continue
        for res in model[chain]:
            hetflag, resi, _icode = res.id
            if hetflag.strip():
                continue
            if not (start <= resi <= end):
                continue
            for a in res:
                b = float(a.get_bfactor() or 0.0)
                any_atom.append(b)
                if a.get_name() == "CA":
                    ca.append(b)
    return ca if ca else any_atom


def compute_per_loop_plddt(
    af_model_path: str | Path,
    gap_ranges: list[dict[str, Any]],
    method: str = "alphafold",
) -> list[dict[str, Any]]:
    """For each gap range, extract per-residue B-factors from ``af_model_path``
    and return mean plDDT (if valid), mean B-factor, and a warning note when
    values are outside the physical [0, 100] plDDT range.

    Parameters
    ----------
    af_model_path : path to AlphaFold-filled model PDB (B-factor column = plDDT).
    gap_ranges    : list of {chain, start, end} dicts (integer residue numbers).
    method        : "alphafold" | "modeller" | "unknown"; MODELLER fills have
                    no plDDT so mean_plddt is always None even if extraction
                    succeeds.

    Returns
    -------
    List of {chain, start, end, method, mean_plddt, mean_bfactor, n_atoms,
    plddt_note?} — one entry per input gap range.
    """
    af_path = Path(af_model_path) if af_model_path else None
    out: list[dict[str, Any]] = []
    for g in gap_ranges:
        rec: dict[str, Any] = {
            "chain": g.get("chain"),
            "start": g.get("start"),
            "end": g.get("end"),
            "method": method,
            "mean_plddt": None,
            "mean_bfactor": None,
            "n_atoms": 0,
        }
        if not af_path or not af_path.is_file():
            out.append(rec)
            continue
        bfacs = _extract_bfactors_in_range(
            pdb_path=af_path,
            chain=rec["chain"],
            start=rec["start"],
            end=rec["end"],
        )
        if not bfacs:
            out.append(rec)
            continue
        mean_b = round(sum(bfacs) / len(bfacs), 1)
        rec["mean_bfactor"] = mean_b
        rec["n_atoms"] = len(bfacs)
        # AlphaFold plDDT is bounded [0, 100]. Values outside imply the
        # B-factor column carries crystallographic temperature factors
        # (typical of hybrid AF+crystal outputs, e.g. FRUTON merged fills).
        # Report a warning; do not claim plDDT in that case.
        if method == "alphafold" and 0 <= mean_b <= 100:
            rec["mean_plddt"] = mean_b
        elif method == "alphafold":
            rec["plddt_note"] = (
                "mean B-factor outside [0, 100] — value likely reflects "
                "crystallographic B-factors carried through the fill/merge; "
                "not a valid plDDT. Manual review recommended."
            )
        # For MODELLER (or method="unknown"), we deliberately leave mean_plddt
        # as None even if mean_bfactor is present — those B-factors do not
        # represent AlphaFold confidence.
        out.append(rec)
    return out


def parse_gap_ranges_from_remark_465(
    raw_pdb_path: str | Path,
) -> list[dict[str, Any]]:
    """Parse REMARK 465 (author-declared missing residues) from a raw PDB
    file and group consecutive residues into gap segments per chain.

    Returns list of {chain, start, end, size} — the same format
    ``compute_per_loop_plddt`` expects for ``gap_ranges``.

    Callers who need geometric breaks (undeclared physical gaps with
    contiguous residue numbering) should use
    ``gaps.geometric_breaks_by_chain_for_pdb`` instead.
    """
    try:
        struct = _PARSER.get_structure("s", str(raw_pdb_path))
    except Exception:
        return []
    missing = struct.header.get("missing_residues") or []
    # Bio.PDB missing_residues: list of {"model", "res_name", "chain", "ssseq",
    # "insertion"} dicts (older 'ssseq' spelling in Bio.PDB).
    parsed = []
    for r in missing:
        try:
            resi = int(r.get("ssseq") if r.get("ssseq") is not None else r.get("res_seq"))
        except (ValueError, TypeError):
            continue
        parsed.append({
            "chain": r.get("chain", "?"),
            "resi": resi,
        })
    if not parsed:
        return []
    parsed.sort(key=lambda x: (x["chain"], x["resi"]))
    gaps: list[dict[str, Any]] = []
    cur: dict[str, Any] | None = None
    for r in parsed:
        if cur is None:
            cur = {"chain": r["chain"], "start": r["resi"], "end": r["resi"]}
        elif r["chain"] == cur["chain"] and r["resi"] == cur["end"] + 1:
            cur["end"] = r["resi"]
        else:
            cur["size"] = cur["end"] - cur["start"] + 1
            gaps.append(cur)
            cur = {"chain": r["chain"], "start": r["resi"], "end": r["resi"]}
    if cur is not None:
        cur["size"] = cur["end"] - cur["start"] + 1
        gaps.append(cur)
    return gaps
