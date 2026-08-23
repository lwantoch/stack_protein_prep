# Why FRUTON does not use pdbfixer

**Status:** hard rule (see `feedback_no_pdbfixer` in memory). Applies to
gap-fill, missing-residue reconstruction, hydrogen assignment, and
side-chain heavy-atom repair.

**Historical note:** pdbfixer was removed from FRUTON on 2026-08-22 and
replaced by Bio.PDB + MODELLER + `OpenMM.app.Modeller.addHydrogens`
(the last of which is *not* pdbfixer — it is OpenMM's own topology-aware
protonation helper).

This document is the chemistry-defensible justification a reviewer would
demand. It is not an anti-pdbfixer polemic — pdbfixer is a fine tool
for *quickly getting a PDB into OpenMM*. The point is that FRUTON aims
for publication-grade model preparation for downstream MM/MMBSA / MD
work, and there the failure modes below become disqualifying.

## Failure modes

### 1. pH-based histidine tautomer assignment ignores experimental evidence

pdbfixer's `addMissingHydrogens(pH=7.0)` assigns every histidine to
HIE (Nε-H) by default. In a metalloenzyme where the paper explicitly
identifies a residue as the catalytic base or proton donor, this is
wrong.

Concrete cases in the MMBSA_200 bench:

- **Carbonic anhydrase II** — H94, H96, H119 coordinate Zn²⁺ via Nε;
  the fourth ligand rotates and gets a δ-protonation (HID) so that
  Nε is deprotonated and available to Zn. pdbfixer would give all
  three HIE → wrong charge on the active-site.
- **Serine proteases** — the "catalytic triad" His is HIP (both Nε
  and Nδ protonated) transiently during the acyl-enzyme step per
  Craik et al. 2011.
- **Cytochrome c oxidase** — coordinating His is HID (Nδ-H); the
  δ-proton hydrogen-bonds with a nearby carboxylate.

FRUTON's approach: `_paper_override_suggest.py` scans
`paper_evidence.md` per protein for regex phrases like "catalytic base",
"proton donor", "nucleophile", "thiolate", and emits a *suggested*
override JSON. The pipeline's owner reviews and promotes it to a real
`protonation_overrides.json`. See `[[reference_fruton_pipeline]]` and
memory `[[project_fruton_lit_review]]`.

### 2. Silent heteroatom removal breaks cofactor + metal coordination

`pdbfixer.removeHeterogens(keepWater=False)` strips *all* HETATM
records except water. In the MMBSA_200 audit
(`[[project_mmbsa200_nonstd_audit]]`), 40+ cofactor residue types were
detected across 200 proteins (NAD, NAP, FAD, HEM, SAM, ATP, ADP, GDP,
PLP, TPP, COA, biotin, …). Silent removal:

- Loses hydrogen-bond partners that the reviewer sees in Fig. 2 of
  the paper (e.g. NAD-nicotinamide stacking against a Tyr).
- Removes the metal ion itself (Zn²⁺, Mg²⁺, Fe²⁺/³⁺, Ni²⁺, Cu²⁺,
  Mn²⁺, Ca²⁺) → the downstream MM run runs "apo" without noticing.
- Removes the seleno-Cys sulfur, arsenic-substituted cysteines
  (CAS), and other exotic residues we deliberately parametrise.

FRUTON's approach: cofactors are *detected + preserved* (see
`_metalls_check_core.py`, `_metal_hydrogen_cleanup.py`); metal ions
are parametrised via MCPB workflow (see
`[[project_gpubench]]` / `_metall_params_mcpb.py`); user is told what
extra parameter sets need to be supplied (see
`feedback_fruton_cofactors_not_parametrized`).

### 3. Non-standard residue conversion loses chemistry

`pdbfixer.replaceNonstandardResidues()` maps modified residues to
their canonical parent by name-only (MSE → MET, SEP → SER, TPO → THR,
PTR → TYR, etc.). This is a chemistry-loss:

- **SEP / TPO / PTR** (phospho-Ser/Thr/Tyr) → the phosphate group
  disappears. In signalling-protein bench items (kinases + their
  substrates), that is the site of biological action.
- **MSE** (seleno-Met) → the selenium becomes sulfur; that changes
  electron density and, more importantly, the crystallographer
  chose MSE precisely because it was needed for SAD/MAD phasing —
  the biological equivalent is Met, but for MM-level accuracy
  the LJ params differ.
- **CSO / CSD / CAS** (oxidised Cys or arsenic-Cys) → collapses to
  CYS; loses the redox-active oxidation state that the paper
  documents.
- **KCX, LLP, PYL, ALY** and other post-translationally modified
  residues become plain LYS.

FRUTON's approach: `_nonstd_residue_params_core.py` catalogues 60+
modified residues with a curated action per residue (parametrise via
antechamber + parmchk2, or preserve as-is if AMBER already has params).
The 20-phospho + 5-Cys-ox + 6-exotic detected in the MMBSA200 audit
(`[[project_mmbsa200_nonstd_audit]]`) get correct parametrisation.

### 4. Backbone reconstruction has no chirality check

`pdbfixer.findMissingResidues` + `addMissingAtoms` uses template
geometry from an internal library. Two failure modes we observed on
5M7U during the removal (see also memory `[[project_fruton_af_hybrid]]`):

- New residues can be inserted with Cα in the D-chirality (signed
  volume of N-CA-C-CB tetrahedron flips sign). AMBER ff14SB assumes
  L-throughout; a D residue silently produces bogus geometry after
  minimisation.
- Backbone dihedrals of new residues are not checked against the
  Ramachandran plot — they land in `disallowed` regions ~10% of the
  time.

FRUTON's approach: `_filler_loop_refine.py` uses MODELLER LoopModel
with `reject_new_chirality_d=True` (signed-volume + improper-dihedral
check per commit `[JCTC R3]`); the post-fill quality gate
(`_filler_quality_check.py`) enforces Ramachandran + peptide-bond
planarity + chirality.

### 5. Terminal atom rebuild races the protonation state

pdbfixer's `_adjustBonds` step rebuilds terminal OXT / N-terminal
H₃⁺ atoms *before* the protonation-state assignment. On residues
whose terminal state was intentionally set (e.g. an aspartate at a
free C-terminus that should be neutral ASH because the crystallographer
observed the H at pH 5 in the deposition), the rebuild overwrites
the intent.

FRUTON's approach: OpenMM `Modeller.addHydrogens` receives the
override dict *before* terminal-atom rebuild; the terminal state
is honoured.

## What we use instead

| Task                                          | Replacement                                                     |
|-----------------------------------------------|-----------------------------------------------------------------|
| Missing-residue gap fill                      | MODELLER LoopModel + AF splice with rollback                    |
| Side-chain heavy-atom repair                  | Bio.PDB + MODELLER residue-topology library                     |
| Hydrogen assignment                           | OpenMM `Modeller.addHydrogens` with per-residue override dict   |
| Non-standard residue chemistry preservation   | `_nonstd_residue_params_core.py` + antechamber/parmchk2/gaff2   |
| Cofactor parametrisation                      | `_cofactor_params.py`                                           |
| Metal-ion parametrisation                     | 12-6-4 LJ (`_ion_params.py`) + MCPB (`_metall_params_mcpb.py`)  |

All replacements are free-license (BSD-3 / MIT / academic). MODELLER
requires a free academic license (Sali lab). OpenMM is MIT.
AmberTools (antechamber, parmchk2, tleap, sander) is GPL-compatible.

## References

- Craik et al. (2011) *Biochem. J.* — serine protease catalytic-triad
  protonation states.
- MacArthur & Thornton (1991) *J. Mol. Biol.* — cis-Pro frequency and
  the reviewer expectation of ≤ 1% cis-nonPro.
- Harding (2001) *Acta Cryst.* — metal-donor distance tolerances.
- Zheng et al. (2008) *CheckMyMetal* — coordination-number and
  geometry validation.
- Li & Merz (2013–2020) — 12-6-4 LJ parameter series for divalent
  and trivalent metal ions.
- Roe & Brooks (2020) *J. Chem. Theory Comput.* — Monte-Carlo
  barostat recommendation for AMBER production MD.
- Frishman & Argos (1995) *Proteins* — DSSP 8-state → 3-state
  mapping used by FRUTON's SS3 IDR cross-check.
- Touw et al. (2015) — DSSP modern maintenance / mkdssp.
