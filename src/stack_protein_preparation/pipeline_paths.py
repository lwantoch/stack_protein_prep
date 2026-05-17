"""Pure path-builder helpers for the FRUTON pipeline.

All functions in this module are pure: they only depend on ``re`` and
``pathlib.Path``.  No pipeline globals are referenced here.
"""

from __future__ import annotations

import re
from pathlib import Path


def _build_raw_protein_path(pdb_id: str, pdb_dir: Path) -> Path:
    return pdb_dir / "components" / f"{pdb_id}_protein.pdb"


def _build_monomer_protein_path(pdb_id: str, pdb_dir: Path) -> Path:
    return pdb_dir / "components" / f"{pdb_id}_protein_monomer.pdb"


def _build_representative_unit_path(pdb_id: str, pdb_dir: Path) -> Path:
    return pdb_dir / "components" / f"{pdb_id}_representative_unit.pdb"


def _build_sanitized_protonation_path(
    *,
    pdb_id: str,
    pdb_dir: Path,
    variant_label: str,
) -> Path:
    """Return a stable temporary sanitized input path for protonation.

    Sanitized structures are implementation artefacts of the GROMACS route.
    They are kept under ``components/.tmp_sanitized_protonation`` so they do not
    look like final FRUTON products, but they remain inspectable when a route
    fails. The variant label is normalized before it is used as a filename so
    fragment labels and future model-source labels cannot create nested paths.
    """

    safe_variant_label = re.sub(r"[^A-Za-z0-9_.-]+", "_", variant_label).strip("_")
    if not safe_variant_label:
        safe_variant_label = "variant"

    sanitize_dir = pdb_dir / "components" / ".tmp_sanitized_protonation"
    sanitize_dir.mkdir(parents=True, exist_ok=True)
    return sanitize_dir / f"{pdb_id}_{safe_variant_label}_sanitized.pdb"


def _build_sanitize_log_path(
    *,
    pdb_id: str,
    pdb_dir: Path,
    variant_label: str,
) -> Path:
    """Return the per-protein sanitizer log path for one protonation input."""

    safe_variant_label = re.sub(r"[^A-Za-z0-9_.-]+", "_", variant_label).strip("_")
    if not safe_variant_label:
        safe_variant_label = "variant"

    return pdb_dir.parent / "logs" / pdb_id / f"sanitize_{safe_variant_label}.log"


def _build_standard_residue_repair_path(
    *,
    pdb_id: str,
    pdb_dir: Path,
    variant_label: str,
) -> Path:
    """Return the temporary MODELLER repair path for one protonation input."""

    safe_variant_label = "".join(
        character if character.isalnum() or character in {"_", "-"} else "_"
        for character in variant_label.strip()
    ).strip("_")
    repair_dir = pdb_dir / "components" / ".tmp_standard_residue_repair"
    repair_dir.mkdir(parents=True, exist_ok=True)
    return repair_dir / f"{pdb_id}_{safe_variant_label}_standard_repaired.pdb"


def _build_standard_residue_repair_log_path(
    *,
    pdb_id: str,
    pdb_dir: Path,
    variant_label: str,
) -> Path:
    """Return the repair log path for one protonation input."""

    safe_variant_label = "".join(
        character if character.isalnum() or character in {"_", "-"} else "_"
        for character in variant_label.strip()
    ).strip("_")
    return pdb_dir.parent / "logs" / pdb_id / f"standard_residue_repair_{safe_variant_label}.log"


def _build_sanitized_monomer_protein_path(
    pdb_id: str,
    pdb_dir: Path,
) -> Path:
    """Return the stable sanitized representative protein path.

    This file is created after representative-unit selection and component
    splitting. It is the earliest safe place to sanitize the protein-only input
    because doing it before representative-unit/component splitting would either
    drop ligands/metals too early or sanitize the wrong multimeric structure.
    Downstream default protein selection prefers this file when it exists.
    """

    return pdb_dir / "components" / f"{pdb_id}_protein_monomer_sanitized.pdb"


def _build_representative_sanitize_log_path(
    pdb_id: str,
    pdb_dir: Path,
) -> Path:
    """Return the sanitizer log path for the representative protein input."""

    return pdb_dir.parent / "logs" / pdb_id / "sanitize_representative_protein.log"


def _build_sanitized_prepared_variant_path(
    *,
    pdb_id: str,
    pdb_dir: Path,
    variant_label: str,
) -> Path:
    """Return the sanitized protein-only audit path for one prepared variant.

    Step 12 sanitizes the produced prepared-structure model rather than the
    starting representative protein. The output is kept in an implementation
    directory because it is an audit artefact: it removes ligands, waters,
    metals, and other non-polymer heterogens so ``parameter_audit.py`` can test
    the actual produced protein model in a protein-only GROMACS probe.
    """

    safe_variant_label = re.sub(r"[^A-Za-z0-9_.-]+", "_", variant_label).strip("_")
    if not safe_variant_label:
        safe_variant_label = "variant"

    audit_dir = pdb_dir / "components" / ".tmp_sanitized_variant_audit"
    audit_dir.mkdir(parents=True, exist_ok=True)
    return audit_dir / f"{pdb_id}_{safe_variant_label}_prepared_sanitized.pdb"


def _build_prepared_variant_sanitize_log_path(
    *,
    pdb_id: str,
    pdb_dir: Path,
    variant_label: str,
) -> Path:
    """Return the sanitizer log path for a produced prepared variant."""

    safe_variant_label = re.sub(r"[^A-Za-z0-9_.-]+", "_", variant_label).strip("_")
    if not safe_variant_label:
        safe_variant_label = "variant"

    return pdb_dir.parent / "logs" / pdb_id / f"sanitize_prepared_{safe_variant_label}.log"


def _build_variant_parameter_audit_manifest_path(
    *,
    pdb_id: str,
    pdb_dir: Path,
    variant_label: str,
) -> Path:
    """Return the JSON manifest path for one prepared-variant audit."""

    safe_variant_label = re.sub(r"[^A-Za-z0-9_.-]+", "_", variant_label).strip("_")
    if not safe_variant_label:
        safe_variant_label = "variant"

    return pdb_dir / f"parameter_audit_{safe_variant_label}.json"


def _build_variant_parameter_audit_log_path(
    *,
    pdb_id: str,
    pdb_dir: Path,
    variant_label: str,
) -> Path:
    """Return the text log path for one prepared-variant audit."""

    safe_variant_label = re.sub(r"[^A-Za-z0-9_.-]+", "_", variant_label).strip("_")
    if not safe_variant_label:
        safe_variant_label = "variant"

    return pdb_dir.parent / "logs" / pdb_id / f"parameter_audit_{safe_variant_label}.log"


def _build_fragment_variant_label(
    base_variant_label: str,
    fragment_index: int,
) -> str:
    return f"{base_variant_label}_fragment_{fragment_index:02d}"
