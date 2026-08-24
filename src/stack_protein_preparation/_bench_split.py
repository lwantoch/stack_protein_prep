"""Bench train / validate / holdout split registry.

USER MANDATE 2026-08-23: proper stat-hygiene split so FRUTON cannot be
tuned to its evaluation set.

Registered splits:

  train   — MMBSA_75 (75 PDBs).  Policy calibration + threshold tuning
            happens here.  If we EVER need a hyper-parameter, this is
            the only set it may be fit on.

  val     — MMBSA_125 (124 PDBs).  Disjoint from train.  Used to check
            that the policy transfers — same tier-distribution should
            appear on val as on train.  No tuning permitted here.

  holdout — 30 externally-picked PDBs, NOT in the MMBSA_200 CSV.
            Final blind test; only touched once per release.

The MMBSA_75 and MMBSA_125 directory-listings on LUSTRE
(``$LUSTRE/MMBSA_200/MMBSA_75/`` and ``.../MMBSA_125/``) are the
authoritative sources: this module reads the ``pdb_ids.csv`` inside
each (or falls back to directory glob).

The 30-PDB holdout is stored in the packaged file
``data/holdout30_pdb_ids.csv`` and can be regenerated via
``scripts/build_holdout30.py`` (uses PDB search with a diverse-family
seed; runs offline against a curated wishlist when the search fails).
"""
from __future__ import annotations

import csv
import re
from importlib import resources
from pathlib import Path
from typing import Iterable


PDB_ID_RE = re.compile(r"^[0-9][A-Za-z0-9]{3}$")


def _load_ids_from_dir(root: Path) -> list[str]:
    """Read a bench-set directory.  Prefers pdb_ids.csv, else globs subdirs."""
    csv_path = root / "pdb_ids.csv"
    if csv_path.is_file():
        with csv_path.open() as fh:
            reader = csv.reader(fh)
            first = next(reader, None)
            body = [row[0] for row in reader if row and PDB_ID_RE.match(row[0].strip())]
            # First row might be a header — include it if it matches the regex
            if first and first[0] and PDB_ID_RE.match(first[0].strip()):
                body.insert(0, first[0].strip())
            return sorted(set(body))
    return sorted(
        d.name for d in root.iterdir()
        if d.is_dir() and PDB_ID_RE.match(d.name)
    )


def load_holdout_ids() -> list[str]:
    """Load the packaged 30-PDB holdout list."""
    text = (
        resources.files("stack_protein_preparation")
        .joinpath("data/holdout30_pdb_ids.csv")
        .read_text()
    )
    lines = [ln.strip() for ln in text.splitlines() if ln.strip() and not ln.startswith("#")]
    # Skip header if first line isn't a valid PDB id
    if lines and not PDB_ID_RE.match(lines[0]):
        lines = lines[1:]
    return [ln for ln in lines if PDB_ID_RE.match(ln)]


def _load_packaged_csv_ids(filename: str) -> list[str]:
    text = (
        resources.files("stack_protein_preparation")
        .joinpath(f"data/{filename}").read_text()
    )
    lines = [ln.strip() for ln in text.splitlines() if ln.strip() and not ln.startswith("#")]
    if lines and not PDB_ID_RE.match(lines[0]):
        lines = lines[1:]
    return [ln for ln in lines if PDB_ID_RE.match(ln)]


def load_stresstest_30() -> list[str]:
    """USER MANDATE 2026-08-24: 30 random PDBs drawn from MMBSA_200
    (seed 20260824) used as STRESS TEST (BL-Pose fallback when no AF)."""
    return _load_packaged_csv_ids("stresstest_30_pdb_ids.csv")


def load_affinity_bench_27() -> list[str]:
    """USER MANDATE 2026-08-24: newbench_27 = the AFFINITY BENCHMARK
    (27 PDBs on /mnt/netapp1/Store_othcxlwa/newbench_27/)."""
    return _load_packaged_csv_ids("affinity_bench_27_pdb_ids.csv")


def load_split(
    split: str,
    lustre_root: Path | str | None = None,
) -> list[str]:
    """Return the list of PDB IDs for the named split.

    Args:
        split: 'train', 'val', or 'holdout'.
        lustre_root: override; defaults to ``$LUSTRE/MMBSA_200``.
    """
    import os
    if split == "holdout":
        return load_holdout_ids()
    if split == "stresstest_30":
        return load_stresstest_30()
    if split == "affinity_bench_27":
        return load_affinity_bench_27()

    lustre = Path(lustre_root or os.environ.get("LUSTRE", "")) / "MMBSA_200"
    if split == "train":
        return _load_ids_from_dir(lustre / "MMBSA_75")
    if split in ("val", "test"):
        return _load_ids_from_dir(lustre / "MMBSA_125")
    raise ValueError(
        f"unknown split {split!r}; "
        f"expect train|val|test|holdout|stresstest_30|affinity_bench_27"
    )


def filter_to_af_ready(pdb_ids: Iterable[str], fruton_root: Path | str) -> list[str]:
    """Restrict to PDBs that have an AF alignment on disk."""
    root = Path(fruton_root)
    ready = []
    for pid in pdb_ids:
        globbed = list(root.glob(f"{pid}/fasta/alignments/filler/*/alphafold/alphafold_aligned_model.pdb"))
        if globbed:
            ready.append(pid)
    return ready


def assert_splits_disjoint(train: Iterable[str], val: Iterable[str], holdout: Iterable[str]) -> None:
    """Raise if any two splits overlap — reviewer-critical invariant."""
    t, v, h = set(train), set(val), set(holdout)
    overlap_tv = t & v
    overlap_th = t & h
    overlap_vh = v & h
    problems = []
    if overlap_tv:
        problems.append(f"train ∩ val = {sorted(overlap_tv)}")
    if overlap_th:
        problems.append(f"train ∩ holdout = {sorted(overlap_th)}")
    if overlap_vh:
        problems.append(f"val ∩ holdout = {sorted(overlap_vh)}")
    if problems:
        raise ValueError("split contamination: " + "; ".join(problems))
