"""Variant-planning helpers for the FRUTON pipeline.

These functions handle gap detection, model-source detection, available-model
field computation, and variant-plan construction.  They depend on
pipeline_state constants and pipeline_paths but have no fruton globals.
"""

from __future__ import annotations

from pathlib import Path

from stack_protein_preparation.pipeline_state import (
    AVAILABLE_MODELS_COLUMN_NAME,
    GAP_SIZES_COLUMN_NAME,
    HAS_METALS_COLUMN_NAME,
    PREPARED_ALPHAFOLD_OUTPUT_PATH_COLUMN_NAME,
    PREPARED_GAPS_OUTPUT_PATH_COLUMN_NAME,
    PREPARED_MODELLER_OUTPUT_PATH_COLUMN_NAME,
    PREPARED_STRUCTURE_VARIANT_COLUMN_NAME,
    RETAIN_ALPHAFOLD_VARIANT_COLUMN_NAME,
    RETAIN_GAPS_VARIANT_COLUMN_NAME,
    RETAIN_MODELLER_VARIANT_COLUMN_NAME,
    VARIANT_POLICY_COLUMN_NAME,
)
from stack_protein_preparation.pipeline_paths import (
    _build_monomer_protein_path,
    _build_raw_protein_path,
)


def _detect_model_source_from_result(
    *,
    modeller_model_path: Path | None,
    alphafold_model_path: Path | None,
) -> str:
    if alphafold_model_path is not None:
        return "alphafold"
    if modeller_model_path is not None:
        return "modeller"
    return ""


def _parse_gap_sizes(gap_sizes_text: str) -> list[int]:
    text = str(gap_sizes_text).strip().lower()
    if not text or text == "none":
        return []

    result: list[int] = []
    for token in text.split("|"):
        token = token.strip()
        if not token:
            continue
        try:
            result.append(int(token))
        except ValueError:
            continue
    return result


def _get_max_gap_size(pipeline_record: dict[str, str]) -> int:
    gap_sizes = _parse_gap_sizes(pipeline_record.get(GAP_SIZES_COLUMN_NAME, ""))
    return max(gap_sizes) if gap_sizes else 0


def _build_available_models_value(
    *,
    has_initial: bool,
    has_gaps: bool,
    has_modeller: bool,
    has_alphafold: bool,
) -> str:
    model_list: list[str] = []

    if has_initial:
        model_list.append("initial")
    if has_gaps:
        model_list.append("gaps")
    if has_modeller:
        model_list.append("modeller")
    if has_alphafold:
        model_list.append("alphafold")

    return " | ".join(model_list)


def _update_available_models_field(
    pipeline_record: dict[str, str],
) -> None:
    prepared_variant_text = str(
        pipeline_record.get(PREPARED_STRUCTURE_VARIANT_COLUMN_NAME, "")
    ).strip()
    prepared_variant_set = {
        token.strip() for token in prepared_variant_text.split(",") if token.strip()
    }

    has_initial = "single" in prepared_variant_set
    has_gaps = bool(
        str(pipeline_record.get(PREPARED_GAPS_OUTPUT_PATH_COLUMN_NAME, "")).strip()
    ) or ("gaps" in prepared_variant_set)

    has_modeller = bool(
        str(pipeline_record.get(PREPARED_MODELLER_OUTPUT_PATH_COLUMN_NAME, "")).strip()
    )
    has_alphafold = bool(
        str(pipeline_record.get(PREPARED_ALPHAFOLD_OUTPUT_PATH_COLUMN_NAME, "")).strip()
    )

    pipeline_record[AVAILABLE_MODELS_COLUMN_NAME] = _build_available_models_value(
        has_initial=has_initial,
        has_gaps=has_gaps,
        has_modeller=has_modeller,
        has_alphafold=has_alphafold,
    )


def _set_variant_policy_fields(
    pipeline_record: dict[str, str],
    n_gaps: int,
    max_gap_size: int,
    complete_model_source: str,
) -> None:
    if n_gaps == 0:
        pipeline_record[VARIANT_POLICY_COLUMN_NAME] = "single_best_available"
        pipeline_record[RETAIN_GAPS_VARIANT_COLUMN_NAME] = "no"
        pipeline_record[RETAIN_MODELLER_VARIANT_COLUMN_NAME] = "no"
        pipeline_record[RETAIN_ALPHAFOLD_VARIANT_COLUMN_NAME] = "no"
        return

    has_complete = bool(complete_model_source)

    if has_complete:
        if max_gap_size > 8:
            pipeline_record[VARIANT_POLICY_COLUMN_NAME] = "gaps_plus_large_gap_complete"
        else:
            pipeline_record[VARIANT_POLICY_COLUMN_NAME] = "gaps_plus_best_complete"
    else:
        pipeline_record[VARIANT_POLICY_COLUMN_NAME] = "gaps_only_no_complete_model"

    pipeline_record[RETAIN_GAPS_VARIANT_COLUMN_NAME] = "yes"
    pipeline_record[RETAIN_MODELLER_VARIANT_COLUMN_NAME] = (
        "yes" if complete_model_source == "modeller" else "no"
    )
    pipeline_record[RETAIN_ALPHAFOLD_VARIANT_COLUMN_NAME] = (
        "yes" if complete_model_source == "alphafold" else "no"
    )


def _choose_variant_structure_input_path(
    pdb_id: str,
    pdb_dir: Path,
    model_source: str,
    model_path: Path | None,
) -> tuple[Path, str]:
    if model_source == "default":
        return _select_default_protein_input_path(pdb_id=pdb_id, pdb_dir=pdb_dir), "protein"

    if model_path is None or not model_path.exists():
        raise FileNotFoundError(
            f"Complete-model input missing for {pdb_id}: {model_source=} {model_path=}"
        )

    if model_source == "modeller":
        return model_path, "modeller"

    if model_source == "alphafold":
        return model_path, "alphafold"

    raise ValueError(f"Unsupported model source: {model_source!r}")


def _select_default_protein_input_path(
    pdb_id: str,
    pdb_dir: Path,
) -> Path:
    """Return the original representative monomer protein input.

    Gap detection and fragment splitting must see the original residue numbering
    from ``components/<PDBID>_protein_monomer.pdb``. Sanitized copies are used
    later for parameter audit and GROMACS-facing protonation inputs, but they do
    not replace the canonical monomer path used for structural decisions.
    """

    monomer_protein_path = _build_monomer_protein_path(pdb_id=pdb_id, pdb_dir=pdb_dir)

    if monomer_protein_path.exists() and monomer_protein_path.stat().st_size > 0:
        return monomer_protein_path

    raw_protein_path = _build_raw_protein_path(pdb_id=pdb_id, pdb_dir=pdb_dir)

    raise FileNotFoundError(
        f"Representative monomer protein input not found for {pdb_id}: "
        f"{monomer_protein_path}. Raw protein path is not used as a fallback: "
        f"{raw_protein_path}."
    )


def _get_retained_variant_plan(
    n_gaps: int,
    max_gap_size: int,
    filler_model_source: str,
    filler_model_path: Path | None,
) -> list[dict[str, str | Path | None]]:
    variant_plan_list: list[dict[str, str | Path | None]] = []

    if n_gaps == 0:
        variant_plan_list.append(
            {
                "variant_label": "single",
                "model_source": "default",
                "model_path": None,
            }
        )
        return variant_plan_list

    variant_plan_list.append(
        {
            "variant_label": "gaps",
            "model_source": "default",
            "model_path": None,
        }
    )

    if filler_model_source and filler_model_path is not None:
        if max_gap_size > 8:
            variant_plan_list.append(
                {
                    "variant_label": "large_gap_complete",
                    "model_source": filler_model_source,
                    "model_path": filler_model_path,
                }
            )
        else:
            variant_plan_list.append(
                {
                    "variant_label": "best_complete",
                    "model_source": filler_model_source,
                    "model_path": filler_model_path,
                }
            )

    return variant_plan_list


def _choose_representative_successful_variant(
    successful_variant_result_list: list[dict[str, str]],
) -> dict[str, str]:
    """Choose which successful variant becomes the representative output.

    Gap-containing proteins can produce both a fragmented gaps variant and a
    complete filled variant. The fragmented gaps variant is useful to retain, but
    it should not become the main prepared-structure output when a complete model
    exists. This selector keeps the single variant as representative for
    proteins without gaps, prefers complete filled models for proteins with gaps,
    and only falls back to the gaps variant when no complete model was prepared.
    """

    priority = [
        "single",
        "best_complete",
        "large_gap_complete",
        "gaps",
    ]

    by_label = {
        result["variant_label"]: result
        for result in successful_variant_result_list
    }

    for variant_label in priority:
        if variant_label in by_label:
            return by_label[variant_label]

    raise RuntimeError("No successful prepared variant available.")


def _blocks_complete_model_generation(
    pipeline_record: dict[str, str],
) -> tuple[bool, str]:
    """Return whether this protein should skip complete-model generation.

    Complete-model generation means MODELLER or AlphaFold reconstruction. The
    runner only blocks this step for metal-containing structures because the
    current reconstruction route cannot preserve metal coordination chemistry
    reliably. Non-standard residues are audited separately by
    ``parameter_audit.py`` and must not suppress MODELLER/AlphaFold logs by
    themselves. Unsupported residue chemistry should be handled later through
    the residue-parameter workflow rather than by hiding the reconstruction
    attempt.
    """

    has_metals = str(
        pipeline_record.get(HAS_METALS_COLUMN_NAME, "")
    ).strip().lower()

    return False, ""
