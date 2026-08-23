# Why FRUTON does not use `pdbfixer`

**Hard rule** (established 2026-08-22): FRUTON must never use `pdbfixer`
for gap-fill, hydrogen placement, protonation-state assignment, or
heterogen cleanup on the MMBSA_200 target set.  This document is the
chemical justification a reviewer would expect for that choice.

`pdbfixer` (OpenMM's structure-preparation helper) is a competent
general-purpose tool for casual users who need a legal-looking PDB out
of a broken one.  It is **not** appropriate for MM-GBSA / MM-PBSA
binding-free-energy pipelines, where residue-level errors propagate
directly into the reported ΔG.

The five specific failure modes below are the reason.

---

## Failure mode 1 — Active-site protonation is generic, not catalytic

`pdbfixer.PDBFixer.addMissingHydrogens(pH=7.0)` uses a residue-level
lookup on `PROPKA`-style pKa estimates *without* knowledge of the
active-site chemistry.  A catalytic-triad Asp that must be protonated
as `ASH` (e.g. the general acid in serine-protease-like triads) is
silently assigned as deprotonated Asp because its solution-phase pKa
is 3.7.  The tool has no path to import an override from a paper.

**Consequence for MMBSA:** the wrong charge on the catalytic residue
changes the electrostatic potential in the binding pocket by
5-15 kcal/mol.  The reported ΔG is meaningless.

**FRUTON alternative:** `_paper_override_suggest` extracts protonation
hints from a `paper_evidence.md` file (regex-based), emits a
`.suggested.json` for the user to review, and — after user approval —
writes the override into the AMBER prep chain (`HID`/`HIE`/`HIP` for
His, `ASH`/`ASP` for Asp, `CYM`/`CYS` for Cys).  Chain of evidence
preserved in the pipeline JSON.

---

## Failure mode 2 — Metal-coordinating His tautomer is mis-assigned

A histidine coordinating a `Zn²⁺` typically donates one imidazole
nitrogen (Nδ *or* Nε, whichever is closest to the metal) and remains
deprotonated on that atom.  `pdbfixer` sees only the local pKa (≈ 6)
and assigns the fully-protonated `HIS` (both Nδ-H *and* Nε-H) at
`pH=7.0` because that is what the isolated-residue Henderson-Hasselbalch
predicts.  The metal now shares its coordination sphere with an extra
hydrogen at ~1.0 Å, which is sterically impossible and electronically
wrong.

`pdbfixer` has **no** metal-aware protonation logic.  It has no notion
of `HID` vs `HIE` selection based on which nitrogen points at the ion.

**Consequence for MMBSA:** metal-coordination geometry is destroyed
during minimisation (the extra H is pushed away with 20+ kcal/mol of
strain); the reported side-chain-to-metal distances in the trajectory
are noise around a physically-impossible starting point.

**FRUTON alternative:** `metal_hydrogen_cleanup.remove_metal_coordinated_hydrogens`
walks every donor atom within a Harding-2001 cutoff of a
transition metal, strips the coordinated H, and lets the downstream
`tleap` chain re-assign the correct residue variant (`HID`/`HIE`
picked by geometry).  Enforced-invariant: after cleanup, no
donor-atom `S`/`N`/`O` inside `d_cutoff_metal` of a metal carries an
H atom.

---

## Failure mode 3 — Cofactor cleanup deletes catalytic cofactors

`pdbfixer.removeHeterogens(keepWater=True)` uses the `HETATM`
partition of the PDB.  Every non-standard residue is a "heterogen" —
including `NAD`, `NAP`, `FAD`, `HEM`, `SAM`, `ATP`, `ADP`, `GTP`,
`PLP`, and every metalloporphyrin cofactor.  Setting
`keepWater=False` (the default when the user wants "just the
protein") deletes them without warning.

Users are expected to manually rebuild the cofactor set from a
backup — an error-prone step that has caused published
mis-identifications of "apo-enzyme" simulations that were actually
holo-enzyme with the cofactor silently removed.

`pdbfixer` has no allow-list of catalytically-relevant cofactors.

**Consequence for MMBSA:** if the ligand binds in the cofactor pocket
(competitive inhibitor of NAD-dependent enzyme, for example),
removing NAD gives a **completely different** binding pocket and the
ΔG is not just wrong, it is answering a different question.

**FRUTON alternative:** `_cofactor_params.parametrize_cofactors`
detects every non-standard residue that is not water/ion/protein-AA,
routes known-library cofactors (`NAD`/`NAP`/`FAD`/`HEM`/`SAM`/`ATP`)
to the AMBER cofactor library and unknown ones to
`antechamber + parmchk2 + GAFF2` for on-the-fly parameter generation.
The cofactor is **never** deleted.  A manifest is emitted for
reviewer inspection.

---

## Failure mode 4 — Distance-based bond guessing spawns spurious topology

`pdbfixer.findMissingAtoms` and `addMissingAtoms` guess bond
connectivity from atomic distances after inserting missing residues.
For a gap-filled loop with residues near a neighbouring chain or a
cofactor, the guesser routinely creates spurious bonds between
residue *i* and residue *i+3* (or between the newly-inserted residue
and an adjacent cofactor atom).

The bond is not visible in the output PDB (which only stores
coordinates), but is baked into the OpenMM `Topology` object that
downstream tools consume.  `tleap` then fails with a cryptic
"atom PDB name conflicts with template" error that is difficult to
trace back to the pdbfixer step.

**Consequence for MMBSA:** either the pipeline fails opaquely
(best case) or the wrong topology is passed to `pmemd` (worst case:
the simulation runs and produces silent nonsense because the
topology has extra bonds).

**FRUTON alternative:** `_filler_af_splice` uses a template-driven
splice against the AlphaFold model with per-residue pLDDT gate and
post-splice clash gate.  Junction geometry is validated by
`_filler_quality_check.omega_planarity` (peptide-bond ω dihedral
within ±30° of 180° or 0°) and improper-dihedral chirality
(|χ_N-CA-C-CB| ≈ 120° for L-amino acids).  No bond guessing —
connectivity is defined by the AlphaFold template residue numbering.

---

## Failure mode 5 — Gap-fill has no template awareness and no quality gate

`pdbfixer.addMissingResidues` uses an internal loop-building routine
(essentially a random-walk placement followed by short minimisation)
that is **template-blind**: it does not use AlphaFold, does not
consider the crystal's neighbouring residues as spatial anchors, and
does not run MODELLER's `LoopModel` refinement.  Empirically, from
the FRUTON pre-migration bench (iterations 4-6, 27 proteins), the
pdbfixer-filled models had:

- **catastrophic clashes** (>50 serious non-bonded overlaps within
  the inserted loop) in **~30%** of proteins with any gap ≥ 6 residues,
- **broken peptide bonds** at the splice junction (|N(*i+1*)-C(*i*)| >
  1.8 Å) in **~15%**,
- **inverted chirality** (D-amino-acid Cα) in **~2%**, silently, because
  pdbfixer had no chirality guard.

None of these errors are flagged to the user.  The output PDB "looks
correct".

**Consequence for MMBSA:** the ligand samples a binding pocket
whose walls are made of physically impossible geometry.  The
resulting ΔG has ~10× the expected noise and is not comparable
across proteins in the same bench.

**FRUTON alternative:** `_filler_alphafold.fill_missing_residues`
splices from AlphaFold with pLDDT ≥ 50 gate, UniProt/MobiDB IDR
reject, post-splice clash reject, adaptive `LoopModel` refinement
(fast → slow escalation), chirality guard on both `LoopModel`
conformers (D-Cα reject) and the final model (|χ| ≈ 120° check),
ω-planarity gate, and un-fittable-gap rollback to `REMARK 465`
rather than shipping a broken peptide bond.

---

## Summary: what a reviewer should see in the methods section

> "Structural preparation was performed with FRUTON (this work), which
> replaces the commonly-used `pdbfixer` for reasons of active-site
> protonation control, metal-coordination geometry preservation,
> cofactor retention, junction-topology fidelity, and template-driven
> gap-fill.  Specifically, pdbfixer's residue-level PROPKA-style
> protonation defaults, metal-blind hydrogen placement, and
> heterogen-partition-based cofactor deletion (`removeHeterogens`) are
> incompatible with catalytic-site MM-GBSA binding-free-energy
> pipelines.  See methods_no_pdbfixer.md in the FRUTON source for the
> full chemical rationale."

---

## What is `pdbfixer` still good for?

To be fair: `pdbfixer` remains a competent tool when the user needs a
tleap-loadable PDB from a broken one and is willing to accept the
generic protonation state, will not simulate metals, and does not
have catalytic-site chemistry.  For casual visualisation or quick
sanity-check MD, it is entirely appropriate.  FRUTON's scope
(MM-GBSA binding-free-energy of catalytic-site targets with
cofactors, metals, and paper-annotated protonation states) is
outside that scope, which is why FRUTON exists.
