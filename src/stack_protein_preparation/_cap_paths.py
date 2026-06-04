"""Path-builder helpers for the cap (capping) module.

These functions are extracted from cap.py to keep that module within a
manageable size. All public names are re-exported from cap.py.
"""
from __future__ import annotations

from pathlib import Path


def sanitize_variant_label(variant_label: str) -> str:
    """Convert a variant label into a filesystem-safe token."""
    cleaned = "".join(
        character if character.isalnum() or character in {"_", "-"} else "_"
        for character in variant_label.strip()
    )
    cleaned = cleaned.strip("_")
    if not cleaned:
        raise ValueError(f"Invalid empty variant label derived from {variant_label!r}")
    return cleaned


def build_standard_internal_capped_output_path(
    pdb_id: str,
    protein_dir: str | Path,
) -> Path:
    """Return the default internal-capped protein coordinate output path.

    Example: ``{protein_dir}/components/{pdb_id}_protein_internal_capped.pdb``
    """
    protein_dir = Path(protein_dir)
    return protein_dir / "components" / f"{pdb_id}_protein_internal_capped.pdb"


def build_variant_internal_capped_output_path(
    pdb_id: str,
    protein_dir: str | Path,
    variant_label: str,
) -> Path:
    """Return the variant-specific internal-capped protein coordinate output path.

    Example: ``{protein_dir}/components/{pdb_id}_protein_internal_capped_gaps.pdb``
    """
    protein_dir = Path(protein_dir)
    safe_variant_label = sanitize_variant_label(variant_label)
    return (
        protein_dir
        / "components"
        / f"{pdb_id}_protein_internal_capped_{safe_variant_label}.pdb"
    )


def build_nonstandard_capping_output_dir(
    pdb_id: str,
    protein_dir: str | Path,
) -> Path:
    """Return the output directory for nonstandard residue capping models."""
    protein_dir = Path(protein_dir)
    return protein_dir / "nonstandard_params"


def get_default_internal_capped_output_path(
    pdb_directory: str | Path,
    pdb_id: str,
) -> Path:
    """Return the default internal-capped output path for a PDB directory."""
    return build_standard_internal_capped_output_path(
        pdb_id=pdb_id,
        protein_dir=pdb_directory,
    )


def get_variant_internal_capped_output_path(
    pdb_directory: str | Path,
    pdb_id: str,
    variant_label: str,
) -> Path:
    """Return the variant-specific internal-capped output path."""
    return build_variant_internal_capped_output_path(
        pdb_id=pdb_id,
        protein_dir=pdb_directory,
        variant_label=variant_label,
    )


__all__ = [
    "sanitize_variant_label",
    "build_standard_internal_capped_output_path",
    "build_variant_internal_capped_output_path",
    "build_nonstandard_capping_output_dir",
    "get_default_internal_capped_output_path",
    "get_variant_internal_capped_output_path",
]
