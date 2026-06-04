"""Attach ACE/NME capping residues to artificial chain termini produced by internal gaps."""
from __future__ import annotations

import shutil
from pathlib import Path

from stack_protein_preparation.pipeline_logging import (
    append_key_value_block,
    build_module_log_paths,
)
from stack_protein_preparation._cap_paths import (  # noqa: F401  (re-exported)
    build_nonstandard_capping_output_dir,
    build_standard_internal_capped_output_path,
    build_variant_internal_capped_output_path,
    get_default_internal_capped_output_path,
    get_variant_internal_capped_output_path,
    sanitize_variant_label,
)
from stack_protein_preparation._cap_helpers import (  # noqa: F401  (re-exported)
    ChainCappingSummary,
    Fragment,
    InternalGapBoundary,
    MODULE_NAME,
    NonstandardCappingSummary,
    NonstandardResidueRecord,
    _ParsedFirstModel,
    _append_run_log,
    _build_log_path,
    _build_nme_residue,
    _build_ace_residue,
    _cap_with_ptmpsi,
    _has_peptide_bond,
    _is_polymer_residue,
    _is_protein_like_residue,
    _normalize_chain_id,
    _print_chain_summary_table,
    _print_end_box,
    _print_start_box,
    build_capped_nonstandard_residue_models,
    cap_internal_gap_boundaries,
    cap_internal_gaps_in_structure,
    cap_nonstandard_residues_for_pdb_directory,
    find_internal_fragments,
    find_internal_gap_boundaries,
    find_nonstandard_residues,
    summarize_capping_results,
)


def internal_cap_selected_structure(
    *,
    input_pdb_path: str | Path,
    output_pdb_path: str | Path,
    peptide_bond_max_distance: float = 1.8,
    variant_label: str | None = None,
) -> dict[str, str | bool | int | float]:
    input_pdb_path = Path(input_pdb_path).resolve()
    output_pdb_path = Path(output_pdb_path).resolve()
    log_path = _build_log_path(
        input_pdb_path=input_pdb_path,
        variant_label=variant_label,
    )

    _print_start_box(
        input_pdb_path=input_pdb_path,
        output_pdb_path=output_pdb_path,
        log_path=log_path,
        peptide_bond_max_distance=peptide_bond_max_distance,
        variant_label=variant_label,
    )

    summary_list = cap_internal_gaps_in_structure(
        input_pdb_path=input_pdb_path,
        output_pdb_path=output_pdb_path,
        peptide_bond_max_distance=peptide_bond_max_distance,
    )

    _append_run_log(
        log_path=log_path,
        input_pdb_path=input_pdb_path,
        output_pdb_path=output_pdb_path,
        peptide_bond_max_distance=peptide_bond_max_distance,
        variant_label=variant_label,
        summary_list=summary_list,
    )

    _print_chain_summary_table(summary_list)

    total_fragments = sum(x.n_fragments for x in summary_list)
    total_boundaries = sum(x.n_internal_boundaries for x in summary_list)
    total_capped = sum(x.n_boundaries_capped for x in summary_list)

    output_dict: dict[str, str | bool | int | float] = {
        "internal_capping_success": output_pdb_path.is_file()
        and output_pdb_path.stat().st_size > 0,
        "internal_capping_input_path": str(input_pdb_path),
        "internal_capping_output_path": str(output_pdb_path),
        "internal_capping_backend": "BioPython",
        "internal_capping_chain_count": len(summary_list),
        "internal_capping_total_fragments": total_fragments,
        "internal_capping_total_boundaries": total_boundaries,
        "internal_capping_total_capped_boundaries": total_capped,
        "internal_capping_peptide_bond_max_distance": peptide_bond_max_distance,
        "internal_capping_summary_text": summarize_capping_results(summary_list),
    }

    if variant_label is not None:
        output_dict["internal_capping_variant"] = variant_label

    _print_end_box(
        status="success" if bool(output_dict["internal_capping_success"]) else "failed",
        output_pdb_path=output_pdb_path,
        total_fragments=total_fragments,
        total_boundaries=total_boundaries,
        total_capped=total_capped,
    )

    return output_dict


def cap_fragment_termini_for_variant(
    pdb_directory: str | Path,
    pdb_id: str,
    *,
    variant_label: str,
    input_pdb_path: str | Path,
    output_pdb_path: str | Path,
    chain_id: str,
    add_ace: bool,
    add_nme: bool,
) -> dict[str, str | bool | int]:
    """Cap the artificial termini of one gap fragment with PTM-Psi.

    This helper is intentionally narrower than ``internal_cap_selected_structure``.
    The gaps route first splits a disconnected protein into physically connected
    fragments, so each fragment no longer contains an internal break by itself.
    The caller therefore decides from the fragment position whether the fragment
    needs an ACE cap on its N side, an NME cap on its C side, both caps, or no
    caps at all.

    If no cap is requested, the function copies the fragment unchanged to the
    requested output path. This keeps the runner simple because every fragment
    still has one normalized post-capping path before protonation. The function
    does not decide whether a fragment is chemically complete; it only applies
    the explicit cap flags supplied by the gap-fragment route.
    """
    protein_dir = Path(pdb_directory).resolve()
    normalized_pdb_id = str(pdb_id).strip().upper()
    normalized_variant_label = sanitize_variant_label(variant_label)
    input_pdb_path = Path(input_pdb_path).resolve()
    output_pdb_path = Path(output_pdb_path).resolve()
    normalized_chain_id = _normalize_chain_id(chain_id)

    if not input_pdb_path.exists():
        raise FileNotFoundError(f"Fragment PDB not found: {input_pdb_path}")

    log_paths = build_module_log_paths(
        protein_data_dir=protein_dir.parent,
        pdb_id=normalized_pdb_id,
        module_name=MODULE_NAME,
        variant_label=normalized_variant_label,
    )
    log_path = log_paths.module_log_path

    output_pdb_path.parent.mkdir(parents=True, exist_ok=True)

    if add_ace or add_nme:
        _cap_with_ptmpsi(
            input_pdb_path=input_pdb_path,
            output_pdb_path=output_pdb_path,
            chain_id=normalized_chain_id,
            add_ace=add_ace,
            add_nme=add_nme,
        )
    else:
        if input_pdb_path != output_pdb_path:
            shutil.copyfile(input_pdb_path, output_pdb_path)

    success = output_pdb_path.is_file() and output_pdb_path.stat().st_size > 0

    append_key_value_block(
        log_path=log_path,
        title="FRAGMENT TERMINI CAPPING",
        key_value_pairs=[
            ("module", MODULE_NAME),
            ("backend", "PTM-Psi" if (add_ace or add_nme) else "copy"),
            ("pdb_id", normalized_pdb_id),
            ("variant", normalized_variant_label),
            ("input_pdb_path", str(input_pdb_path)),
            ("output_pdb_path", str(output_pdb_path)),
            ("chain_id", normalized_chain_id),
            ("add_ace", str(add_ace)),
            ("add_nme", str(add_nme)),
            ("success", str(success)),
        ],
    )

    return {
        "fragment_capping_success": success,
        "fragment_capping_input_path": str(input_pdb_path),
        "fragment_capping_output_path": str(output_pdb_path),
        "fragment_capping_chain_id": normalized_chain_id,
        "fragment_capping_add_ace": add_ace,
        "fragment_capping_add_nme": add_nme,
        "fragment_capping_n_caps_requested": int(add_ace) + int(add_nme),
        "fragment_capping_variant": normalized_variant_label,
    }


def internal_cap_variant_structure(
    pdb_id: str,
    protein_dir: str | Path,
    *,
    variant_label: str,
    input_pdb_path: str | Path,
    peptide_bond_max_distance: float = 1.8,
) -> dict[str, str | bool | int | float]:
    output_pdb_path = build_variant_internal_capped_output_path(
        pdb_id=pdb_id,
        protein_dir=protein_dir,
        variant_label=variant_label,
    )

    return internal_cap_selected_structure(
        input_pdb_path=input_pdb_path,
        output_pdb_path=output_pdb_path,
        peptide_bond_max_distance=peptide_bond_max_distance,
        variant_label=variant_label,
    )


def internal_cap_protein_structure(
    pdb_id: str,
    protein_dir: str | Path,
    *,
    input_pdb_path: str | Path,
    peptide_bond_max_distance: float = 1.8,
) -> dict[str, str | bool | int | float]:
    output_pdb_path = build_standard_internal_capped_output_path(
        pdb_id=pdb_id,
        protein_dir=protein_dir,
    )

    return internal_cap_selected_structure(
        input_pdb_path=input_pdb_path,
        output_pdb_path=output_pdb_path,
        peptide_bond_max_distance=peptide_bond_max_distance,
        variant_label=None,
    )


def cap_internal_gaps_for_pdb_directory(
    pdb_directory: str | Path,
    pdb_id: str,
    input_pdb_path: str | Path,
    peptide_bond_max_distance: float = 1.8,
) -> tuple[Path, list[ChainCappingSummary]]:
    output_pdb_path = build_standard_internal_capped_output_path(
        protein_dir=pdb_directory,
        pdb_id=pdb_id,
    )

    summary_list = cap_internal_gaps_in_structure(
        input_pdb_path=input_pdb_path,
        output_pdb_path=output_pdb_path,
        peptide_bond_max_distance=peptide_bond_max_distance,
    )

    return output_pdb_path, summary_list


def cap_internal_gaps_for_variant(
    pdb_directory: str | Path,
    pdb_id: str,
    *,
    variant_label: str,
    input_pdb_path: str | Path,
    peptide_bond_max_distance: float = 1.8,
) -> dict[str, str | bool | int | float]:
    return internal_cap_variant_structure(
        pdb_id=pdb_id,
        protein_dir=pdb_directory,
        variant_label=variant_label,
        input_pdb_path=input_pdb_path,
        peptide_bond_max_distance=peptide_bond_max_distance,
    )


if __name__ == "__main__":
    import argparse

    argument_parser = argparse.ArgumentParser(
        description=(
            "Extract protein-like non-standard residues and build ACE/NME capped "
            "model compounds with PTM-Psi."
        )
    )
    argument_parser.add_argument(
        "pdb_directory",
        type=Path,
        help="Protein directory, e.g. data/proteins/1ABC",
    )
    argument_parser.add_argument(
        "pdb_id",
        type=str,
        help="PDB ID, e.g. 1ABC",
    )
    argument_parser.add_argument(
        "input_pdb",
        type=Path,
        help="Input representative-unit or protein PDB path.",
    )
    argument_parser.add_argument(
        "--internal-gaps",
        action="store_true",
        help="Also run PTM-Psi-based internal fragment capping.",
    )
    argument_parser.add_argument(
        "--peptide-bond-max-distance",
        type=float,
        default=1.8,
        help="Maximum C-N distance treated as peptide-bonded. Default: 1.8 Å",
    )
    arguments = argument_parser.parse_args()

    nonstandard_summary = cap_nonstandard_residues_for_pdb_directory(
        pdb_directory=arguments.pdb_directory,
        pdb_id=arguments.pdb_id,
        input_pdb_path=arguments.input_pdb,
    )

    print(f"nonstandard_manifest={nonstandard_summary.manifest_path}")
    print(f"n_nonstandard_residues={nonstandard_summary.n_nonstandard_residues}")
    print(f"n_capped_models_written={nonstandard_summary.n_capped_models_written}")

    if arguments.internal_gaps:
        result = internal_cap_protein_structure(
            pdb_id=arguments.pdb_id,
            protein_dir=arguments.pdb_directory,
            input_pdb_path=arguments.input_pdb,
            peptide_bond_max_distance=arguments.peptide_bond_max_distance,
        )
        print(f"internal_capping_output={result['internal_capping_output_path']}")
