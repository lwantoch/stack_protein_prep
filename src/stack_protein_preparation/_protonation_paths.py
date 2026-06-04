"""Path-builder helpers for the protonation module.

This private sub-module contains only path-building and labelling utilities
extracted from protonation.py to keep that module within a manageable size.
All public names are re-exported from protonation.py via ``from ._protonation_paths import *``.
"""
from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    pass

# ---------------------------------------------------------------------------
# Variant label sanitisation
# ---------------------------------------------------------------------------

def sanitize_variant_label(variant_label: str) -> str:
    """Convert a variant label into a filesystem-safe token.

    Variant labels are used in output filenames and log filenames. This helper
    keeps letters, digits, underscores, and hyphens unchanged while replacing
    everything else by underscores. Empty labels are rejected so accidental blank
    variant names do not overwrite standard protonation outputs.
    """

    cleaned = "".join(
        character if character.isalnum() or character in {"_", "-"} else "_"
        for character in variant_label.strip()
    )
    cleaned = cleaned.strip("_")

    if not cleaned:
        raise ValueError(f"Invalid empty variant label derived from {variant_label!r}")

    return cleaned


# ---------------------------------------------------------------------------
# Coordinate output paths
# ---------------------------------------------------------------------------

def build_standard_protonation_output_path(
    pdb_id: str,
    protein_dir: str | Path,
) -> Path:
    """Return the standard single-path protonated coordinate output.

    Example
    -------
    ``data/proteins/1ABC/components/1ABC_proteinH.pdb``
    """

    protein_dir = Path(protein_dir)
    components_dir = protein_dir / "components"
    return components_dir / f"{pdb_id}_proteinH.pdb"


def build_variant_protonation_output_path(
    pdb_id: str,
    protein_dir: str | Path,
    variant_label: str,
) -> Path:
    """Return the variant-specific protonated coordinate output.

    Example
    -------
    ``components/1ABC_proteinH_gaps.pdb``
    ``components/1ABC_proteinH_best_complete.pdb``
    """

    protein_dir = Path(protein_dir)
    components_dir = protein_dir / "components"
    safe_variant_label = sanitize_variant_label(variant_label)
    return components_dir / f"{pdb_id}_proteinH_{safe_variant_label}.pdb"


# ---------------------------------------------------------------------------
# GROMACS topology / restraints paths
# ---------------------------------------------------------------------------

def build_gromacs_topology_output_path(output_pdb: str | Path) -> Path:
    """Return the GROMACS topology path paired with a protonated coordinate file.

    ``pdb2gmx`` always generates a topology in addition to coordinates. Keeping
    the topology beside the protonated PDB makes the artefacts inspectable and
    prevents topology includes from being scattered across unrelated runtime
    directories. The topology is not used as the authoritative downstream
    topology unless later FRUTON stages explicitly choose to consume it.
    """

    output_pdb = Path(output_pdb)
    return output_pdb.with_suffix(".top")


def build_gromacs_position_restraints_output_path(output_pdb: str | Path) -> Path:
    """Return the position-restraints include path paired with the output PDB."""

    output_pdb = Path(output_pdb)
    return output_pdb.with_name(f"{output_pdb.stem}_posre.itp")


# ---------------------------------------------------------------------------
# Log paths
# ---------------------------------------------------------------------------

def build_protonation_stdout_log_path(
    pdb_id: str,
    protein_dir: str | Path,
    variant_label: str | None = None,
) -> Path:
    protein_dir = Path(protein_dir)
    components_dir = protein_dir / "components"

    if variant_label is None:
        return components_dir / f"{pdb_id}_protonation_stdout.log"

    safe_variant_label = sanitize_variant_label(variant_label)
    return components_dir / f"{pdb_id}_protonation_stdout_{safe_variant_label}.log"


def build_protonation_stderr_log_path(
    pdb_id: str,
    protein_dir: str | Path,
    variant_label: str | None = None,
) -> Path:
    protein_dir = Path(protein_dir)
    components_dir = protein_dir / "components"

    if variant_label is None:
        return components_dir / f"{pdb_id}_protonation_stderr.log"

    safe_variant_label = sanitize_variant_label(variant_label)
    return components_dir / f"{pdb_id}_protonation_stderr_{safe_variant_label}.log"


# ---------------------------------------------------------------------------
# Filler integration
# ---------------------------------------------------------------------------

DEFAULT_FILLER_FINAL_FILENAME = "final_filled_model.pdb"


def find_default_filler_final_model_path(
    filler_output_dir: str | Path,
) -> Path:
    """Return the default unified filler final model path.

    Expected filler convention
    --------------------------
    ``<filler_output_dir>/final_filled_model.pdb``
    """

    filler_output_dir = Path(filler_output_dir)
    return filler_output_dir / DEFAULT_FILLER_FINAL_FILENAME


__all__ = [
    "sanitize_variant_label",
    "build_standard_protonation_output_path",
    "build_variant_protonation_output_path",
    "build_gromacs_topology_output_path",
    "build_gromacs_position_restraints_output_path",
    "build_protonation_stdout_log_path",
    "build_protonation_stderr_log_path",
    "find_default_filler_final_model_path",
    "DEFAULT_FILLER_FINAL_FILENAME",
]
