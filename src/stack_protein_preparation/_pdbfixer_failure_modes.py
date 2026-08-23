"""Structured catalogue of the five documented pdbfixer failure modes.

Companion to ``docs/methods_no_pdbfixer.md``.  A reviewer, a downstream
consumer, or a CI check can import this module to get the machine-
readable form of the no-pdbfixer rationale.

Not a wrapper: does NOT run pdbfixer.  The hard-rule (established
2026-08-22, see [[feedback_no_pdbfixer]]) is that FRUTON never invokes
pdbfixer.  This module *describes* what pdbfixer would do wrong and
which FRUTON module supersedes it.
"""
from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class PdbfixerFailureMode:
    """One documented failure mode: chemistry, consequence, FRUTON substitute."""
    mode_id: str            # short kebab-case slug, e.g. "protonation-generic"
    title: str              # reviewer-facing one-liner
    chemistry: str          # what pdbfixer does wrong at the atomic level
    consequence_for_mmbsa: str
    fruton_alternative: str
    fruton_module: str      # dotted-path to the FRUTON module that supersedes it


FAILURE_MODES: tuple[PdbfixerFailureMode, ...] = (
    PdbfixerFailureMode(
        mode_id="protonation-generic",
        title="Active-site protonation is generic, not catalytic",
        chemistry=(
            "pdbfixer.addMissingHydrogens(pH=7.0) uses residue-level PROPKA-"
            "style pKa lookup with no active-site awareness. A catalytic-"
            "triad Asp that should be ASH is silently deprotonated because "
            "isolated-residue pKa is 3.7."
        ),
        consequence_for_mmbsa=(
            "Wrong charge on catalytic residue perturbs pocket electrostatics "
            "by 5-15 kcal/mol; reported ΔG is meaningless."
        ),
        fruton_alternative=(
            "_paper_override_suggest extracts protonation hints from paper "
            "evidence, emits a .suggested.json, and — after user approval — "
            "writes HID/HIE/HIP/ASH/GLH/CYM overrides into the AMBER prep chain."
        ),
        fruton_module="stack_protein_preparation._paper_override_suggest",
    ),
    PdbfixerFailureMode(
        mode_id="metal-his-tautomer",
        title="Metal-coordinating His tautomer is mis-assigned",
        chemistry=(
            "A His coordinating Zn²⁺ typically donates one imidazole N (Nδ or "
            "Nε, whichever is closest) and remains deprotonated on that atom. "
            "pdbfixer assigns the fully-protonated HIS at pH 7.0 because it "
            "has no metal-aware protonation logic."
        ),
        consequence_for_mmbsa=(
            "Metal coordination sphere destroyed during minimisation (extra H "
            "pushed away with 20+ kcal/mol of strain); metal-donor distances "
            "in trajectory are noise around a physically-impossible start."
        ),
        fruton_alternative=(
            "metal_hydrogen_cleanup.remove_metal_coordinated_hydrogens strips "
            "hydrogens from any N/O/S donor within Harding-2001 cutoff of a "
            "transition metal; tleap then re-assigns the correct HID/HIE "
            "variant by geometry."
        ),
        fruton_module="stack_protein_preparation.metal_hydrogen_cleanup",
    ),
    PdbfixerFailureMode(
        mode_id="cofactor-deletion",
        title="Cofactor cleanup deletes catalytic cofactors",
        chemistry=(
            "pdbfixer.removeHeterogens partitions on the HETATM record; every "
            "non-standard residue is a heterogen, including NAD/NAP/FAD/HEM/"
            "SAM/ATP/ADP/PLP and metalloporphyrins. No catalytic-cofactor "
            "allow-list."
        ),
        consequence_for_mmbsa=(
            "If the ligand binds in a cofactor pocket, silently removing the "
            "cofactor produces a completely different binding pocket; the ΔG "
            "answers a different question than the one being asked."
        ),
        fruton_alternative=(
            "_cofactor_params.parametrize_cofactors detects every non-standard "
            "residue that is not water/ion/protein-AA; routes known cofactors "
            "to the AMBER library and unknown ones to antechamber + parmchk2 "
            "+ GAFF2. Cofactor is never deleted; manifest is emitted."
        ),
        fruton_module="stack_protein_preparation._cofactor_params",
    ),
    PdbfixerFailureMode(
        mode_id="topology-bond-guessing",
        title="Distance-based bond guessing spawns spurious topology",
        chemistry=(
            "pdbfixer.findMissingAtoms + addMissingAtoms guess bond connectivity "
            "from atomic distances after inserting missing residues. Spurious "
            "bonds between residue i and residue i+3, or to adjacent cofactor "
            "atoms, are baked into the OpenMM Topology object."
        ),
        consequence_for_mmbsa=(
            "Either the pipeline fails opaquely (best case, tleap errors) or "
            "the wrong topology reaches pmemd (worst case: silent nonsense in "
            "trajectory)."
        ),
        fruton_alternative=(
            "_filler_af_splice uses template-driven splice with per-residue "
            "pLDDT gate + post-splice clash gate; connectivity is defined by "
            "AlphaFold template residue numbering. Junction validated by "
            "_filler_quality_check ω-planarity and improper-dihedral chirality."
        ),
        fruton_module="stack_protein_preparation._filler_af_splice",
    ),
    PdbfixerFailureMode(
        mode_id="gap-fill-no-template",
        title="Gap-fill has no template awareness and no quality gate",
        chemistry=(
            "pdbfixer.addMissingResidues uses a random-walk placement + short "
            "minimisation with no AlphaFold, no crystal-neighbour anchors, no "
            "MODELLER LoopModel refinement, no chirality guard."
        ),
        consequence_for_mmbsa=(
            "FRUTON pre-migration bench (27 proteins, iterations 4-6): ~30% "
            "of proteins with gaps ≥ 6 residues had catastrophic clashes, "
            "~15% had broken peptide bonds at splice junctions, ~2% had "
            "silently inverted D-amino-acid chirality."
        ),
        fruton_alternative=(
            "_filler_alphafold.fill_missing_residues: AF splice + pLDDT ≥ 50 "
            "gate + MobiDB IDR reject + post-splice clash reject + adaptive "
            "LoopModel refinement (fast → slow escalation) + chirality guard "
            "on conformer and final model + ω-planarity gate + un-fittable-"
            "gap rollback to REMARK 465."
        ),
        fruton_module="stack_protein_preparation._filler_alphafold",
    ),
)


def get_failure_mode(mode_id: str) -> PdbfixerFailureMode:
    """Look up one failure mode by slug or raise KeyError."""
    for m in FAILURE_MODES:
        if m.mode_id == mode_id:
            return m
    raise KeyError(f"unknown pdbfixer failure mode: {mode_id!r}")


def all_mode_ids() -> tuple[str, ...]:
    return tuple(m.mode_id for m in FAILURE_MODES)


def format_methods_paragraph() -> str:
    """One-paragraph methods-section blurb a reviewer would expect."""
    return (
        "Structural preparation was performed with FRUTON, which replaces the "
        "commonly-used pdbfixer for reasons of active-site protonation control, "
        "metal-coordination geometry preservation, cofactor retention, "
        "junction-topology fidelity, and template-driven gap-fill. pdbfixer's "
        "residue-level PROPKA-style protonation defaults, metal-blind hydrogen "
        "placement, and heterogen-partition-based cofactor deletion "
        "(removeHeterogens) are incompatible with catalytic-site MM-GBSA "
        "binding-free-energy pipelines. See docs/methods_no_pdbfixer.md in "
        "the FRUTON source for the full chemical rationale."
    )
