# FRUTON

<p align="center">

```text
███████╗██████╗ ██╗   ██╗████████╗ ██████╗ ███╗   ██╗
██╔════╝██╔══██╗██║   ██║╚══██╔══╝██╔═══██╗████╗  ██║
█████╗  ██████╔╝██║   ██║   ██║   ██║   ██║██╔██╗ ██║
██╔══╝  ██╔══██╗██║   ██║   ██║   ██║   ██║██║╚██╗██║
██║     ██║  ██║╚██████╔╝   ██║   ╚██████╔╝██║ ╚████║
╚═╝     ╚═╝  ╚═╝ ╚═════╝    ╚═╝    ╚═════╝ ╚═╝  ╚═══╝
```

</p>

<p align="center">
  <b>Framework for Reconstruction, UniProt alignment, and Topology-Oriented protein Normalization</b>
</p>

<p align="center">
  A modular protein-preparation pipeline for alignment-aware cleanup, representative-unit selection, gap handling, GROMACS protonation, final prepared-structure assembly, and later metal-site preparation.
</p>

<p align="center">
  <img alt="status" src="https://img.shields.io/badge/status-active%20development-2563eb">
  <img alt="python" src="https://img.shields.io/badge/python-3.12-0ea5e9">
  <img alt="environment" src="https://img.shields.io/badge/environment-pixi-7c3aed">
  <img alt="domain" src="https://img.shields.io/badge/domain-protein%20preparation-059669">
</p>

---

## Overview

**FRUTON** is a stepwise, state-driven protein-preparation framework for turning cropped or otherwise non-trivial structural inputs into reproducible, downstream-ready outputs.

It is designed for protein systems where preparation is not just “clean the PDB and continue”, but a sequence of structural, sequence-based, and chemistry-aware transformations. In practice, this includes cases with:

- insertion codes
- repeated coordinate chains
- heteromeric representative units
- missing internal residues
- UniProt mismatches
- ligands, water, and metals
- protonation ambiguities
- cropped residue ranges that create artificial termini
- disconnected internal fragments that may need later explicit chemical treatment

The project is built around a simple principle:

> **Every important transformation should leave visible evidence on disk.**

That means:

- intermediate files are intentionally written
- representative-unit files are preserved
- alignments and mapping files are preserved
- state is stored explicitly in JSON and XLSX
- downstream steps consume well-defined upstream outputs
- failures are meant to be inspectable, not mysterious

---

## Why FRUTON exists

Protein preparation for modelling or MD is rarely a single command. It is usually a chain of assumptions:

```text
structure → sequence interpretation → cleanup → representative-unit selection → reconstruction → chemistry → final prepared assembly
```

If those assumptions stay implicit, debugging becomes guesswork.

FRUTON exists to make that chain:

- **modular**
- **traceable**
- **reproducible**
- **easy to inspect**
- **safe to extend**

Rather than burying everything in a monolithic script, FRUTON treats preparation as a sequence of explicit modules with explicit files and explicit state updates.

---

## Representative-unit policy

FRUTON no longer treats “monomer” as “always one chain”.

The current rule is:

> **Write one representative protein unit by keeping one representative chain from each non-equivalent chain group.**

Examples:

```text
A/B equivalent
→ representative unit: A
```

```text
A/B equivalent and C/D equivalent
→ representative unit: A + C
```

This matters for heteromeric structures. A duplicated heteromer such as:

```text
A + C
B + D
```

should not become:

```text
A only
```

and should not stay:

```text
A + B + C + D
```

It should become:

```text
A + C
```

### Ligand handling during representative-unit selection

Representative-unit selection is performed on the full insertion-cleaned PDB, not on an already protein-only component file.

That means FRUTON can keep non-water HETATM records assigned to the selected representative chains.

Example:

```text
Detected chain groups:
A/B and C/D

Selected representative chains:
A and C

Ligands:
BBA on A
BBA on B

Representative unit:
protein A + protein C + BBA on A
```

The ligand on discarded duplicate chain `B` is not retained.

This is currently a chain-ID based rule. FRUTON does not yet perform distance-based ligand ownership assignment at this stage.

---

## Pipeline at a glance

```text
[ RAW PDB ]
     │
     ▼
[ FASTA FILES ]
     │
     ▼
[ UNIPROT ALIGNMENT ]
     │
     ▼
[ INSERTION CLEANUP ]
     │
     ▼
[ REPRESENTATIVE UNIT SELECTION ]
     │
     ├── remove duplicated equivalent chains
     ├── keep one representative chain per chain group
     └── keep ligands / cofactors assigned to selected representative chains
     │
     ▼
[ COMPONENT SPLIT FROM REPRESENTATIVE UNIT ]
     │
     ├── representative protein component
     ├── representative ligand component
     ├── representative water component
     └── representative metal component
     │
     ▼
[ GAP DETECTION ON REPRESENTATIVE PROTEIN ]
     │
     ├── no gaps ───────────────────► [ SINGLE VARIANT ] ───────────────────────────────┐
     │                                                                                   │
     └── gaps detected ─────────────┬──► [ GAPS VARIANT ] ───────────────────────────────┤
                                    │                                                    │
                                    └──► [ COMPLETE MODEL VARIANT IF FILLER SUCCEEDS ] ──┘
                                                                                   │
                                                                                   ▼
                                                                    [ GROMACS PROTONATION ]
                                                                    gmx pdb2gmx -ignh
                                                                                   │
                                                                                   ▼
                                                                         [ PREPARED STRUCTURE(S) ]
                                                                                   │
                                                                                   ▼
                                                                       [ LATER METAL PREPARATION ]
```

---

## Current gap-handling policy

FRUTON currently distinguishes retained variants using gap presence, maximum gap size, and whether the filler produced a usable complete model.

### No gaps

If no gaps are detected, FRUTON writes one best-available variant:

```text
single
```

### Gaps detected

If gaps are detected, FRUTON preserves a conservative gapped variant:

```text
gaps
```

If the filler stage produces a complete model, FRUTON also retains a complete variant:

```text
best_complete
```

or, for larger-gap cases:

```text
large_gap_complete
```

The runner currently tracks larger-gap cases using maximum gap size, not only raw gap count:

```text
max_gap_size > 5
max_gap_size >= 8
```

The representative final output prefers complete variants when they exist. The conservative `gaps` variant can still be retained, but it should not become the main representative prepared output when a complete model was successfully prepared.

---

## Current capabilities

| Stage | Purpose | Current output |
|---|---|---|
| **PDB sync** | synchronize dataset scope from CSV | protein directory tree |
| **FASTA generation** | generate structure- and reference-derived sequences | `SEQRES`, `ATOM`, `UniProt` FASTA files |
| **Sequence alignment** | compare PDB chains against UniProt | aligned FASTA, TSV mappings, PNGs |
| **Insertion cleanup** | remove insertion-code ambiguity | `<PDB_ID>_delins.pdb` |
| **Representative-unit selection** | remove duplicated equivalent chains while preserving one representative chain per chain group | `components/<PDB_ID>_representative_unit.pdb` |
| **Component split** | split the representative unit, not the full multimer | protein / water / ligand / metal component files |
| **Gap detection** | quantify missing internal regions in the representative protein | `n_gaps`, `gap_sizes`, `has_gaps` |
| **Filler** | generate reconstruction candidates when justified | MODELLER and/or AlphaFold-derived complete models |
| **Protonation** | add hydrogens / assign GROMACS-compatible protonation state | protonated protein PDB |
| **Prepared structure** | assemble final system for downstream use | `prepared/.../<PDB_ID>.pdb` |
| **Metal preparation** | planned / evolving branch for metal-site parametrization support | combined PDB, contacts output, future parameter-preparation inputs |

---

## Important current chemistry status

The current runner uses this downstream chemistry route:

```text
representative protein/model → gmx pdb2gmx -ignh → prepared_structure
```

Internal capping, AMBER renaming, and AMBER termini modules may exist in the repository, but the current runner path keeps them disabled while the representative-unit, gap, filler, and protonation behavior is being stabilized.

This is intentional for now. It keeps the current debugging surface smaller and avoids mixing representative-unit selection problems with later force-field naming or capping problems.

---

## Project architecture

FRUTON is easiest to understand as a stack of cooperating layers:

| Layer | Role |
|---|---|
| **Sequence layer** | FASTA generation, UniProt matching, alignment TSVs |
| **Structure layer** | insertion cleanup, representative-unit selection, component split, gap detection |
| **Reconstruction layer** | MODELLER / AlphaFold-based reconstruction candidates |
| **Chemistry layer** | current GROMACS protonation route; later AMBER normalization and capping |
| **Prepared-assembly layer** | final merge of prepared protein with waters / ligands / metals |
| **Metal branch** | future preparation for Gaussian / MCPB-style workflows |

This separation is intentional. It prevents alignment logic, structural cleanup, chemistry logic, prepared-assembly logic, and downstream parametrization logic from collapsing into one opaque block.

---

## Repository layout

```text
stack_protein_prep/
├── pixi.toml
├── README.md
├── scripts/
│   └── fruton.py
├── src/
│   └── stack_protein_preparation/
│       ├── fasta_files.py
│       ├── sequence_alignment.py
│       ├── insertion_codes.py
│       ├── monomer.py
│       ├── pdb_components.py
│       ├── gaps.py
│       ├── filler.py
│       ├── fragment_split.py
│       ├── protonation.py
│       ├── prepared_structure.py
│       ├── pipeline_state.py
│       ├── pipeline_table.py
│       └── pipeline_xlsx.py
└── data/
    └── proteins/
        ├── pdb_ids.csv
        ├── pipeline.json
        ├── pipeline.xlsx
        └── <PDB_ID>/
```

Some modules may exist for later or alternate chemistry paths, such as AMBER renaming, termini preparation, internal capping, or metal preparation. They are not necessarily active in the current main runner.

---

## Per-protein layout

A typical protein directory currently looks like this:

```text
data/proteins/<PDB_ID>/
├── <PDB_ID>.pdb
├── <PDB_ID>_delins.pdb
├── fasta/
│   ├── PDB-<PDB_ID>-SEQRES.fasta
│   ├── PDB-<PDB_ID>-ATOM.fasta
│   ├── UniProt_<UNIPROT_ID>.fasta
│   └── alignments/
│       ├── SEQRES_chain_A_vs_UniProt.input.fasta
│       ├── SEQRES_chain_A_vs_UniProt.aln.fasta
│       ├── SEQRES_chain_A_vs_UniProt.mapping.tsv
│       ├── SEQRES_chain_A_vs_UniProt.png
│       ├── ATOM_chain_A_vs_UniProt.input.fasta
│       ├── ATOM_chain_A_vs_UniProt.aln.fasta
│       ├── ATOM_chain_A_vs_UniProt.mapping.tsv
│       └── ATOM_chain_A_vs_UniProt.png
├── components/
│   ├── <PDB_ID>_representative_unit.pdb
│   ├── <PDB_ID>_protein.pdb
│   ├── <PDB_ID>_protein_monomer.pdb
│   ├── <PDB_ID>_water.pdb
│   ├── <PDB_ID>_ligand.pdb
│   └── <PDB_ID>_metals.pdb
├── protonation/
│   └── ...
├── prepared/
│   ├── <PDB_ID>.pdb
│   ├── gaps/
│   │   └── <PDB_ID>.pdb
│   └── complete/
│       └── <PDB_ID>.pdb
└── logs/
    └── ...
```

Not every protein will contain every file. The `prepared/` layout depends on gap detection, filler success, and retained variant logic.

---

## Core outputs and what they mean

| File | Interpretation |
|---|---|
| `<PDB_ID>_delins.pdb` | insertion-cleaned working PDB |
| `components/<PDB_ID>_representative_unit.pdb` | selected representative unit from the full insertion-cleaned structure |
| `components/<PDB_ID>_protein.pdb` | protein component split from the representative unit |
| `components/<PDB_ID>_protein_monomer.pdb` | strict downstream protein input; currently mirrors the representative protein component |
| `components/<PDB_ID>_ligand.pdb` | ligand component split from the representative unit |
| `components/<PDB_ID>_water.pdb` | water component split from the representative unit |
| `components/<PDB_ID>_metals.pdb` | metal component split from the representative unit |
| `PDB-<PDB_ID>-SEQRES.fasta` | sequence derived from SEQRES records |
| `PDB-<PDB_ID>-ATOM.fasta` | sequence derived from observed ATOM records |
| `UniProt_<UNIPROT_ID>.fasta` | external UniProt reference |
| `*_vs_UniProt.aln.fasta` | chain-specific alignment output |
| `*_vs_UniProt.mapping.tsv` | alignment mapping table |
| `protonation/...` | GROMACS protonation outputs |
| `prepared/<PDB_ID>.pdb` | final prepared structure for the single/default representative variant |
| `prepared/gaps/<PDB_ID>.pdb` | final prepared structure for the retained gapped variant |
| `prepared/complete/<PDB_ID>.pdb` | final prepared structure for a reconstructed complete variant |

---

## Why the individual stages matter

### Sequence-aware preparation

FRUTON treats sequence logic as core infrastructure, not as optional annotation. UniProt alignment is used to:

- compare observed structure against reference biology
- support chain-specific reasoning
- identify truncation and missing segments
- provide a stable reference layer for reconstruction

### Representative-unit selection

Many PDB files contain repeated coordinate copies. FRUTON removes duplicated equivalent chains before downstream preparation, but it does not assume that one biological unit is always one chain.

For homomers, this usually produces one representative chain.

For duplicated heteromers, this produces one representative chain from each non-equivalent chain group.

Ligands and cofactors are retained only when assigned to selected representative chains.

### Structural cleanup

Insertion-code cleanup and representative-unit component splitting simplify later reasoning. Once the representative structure is decomposed into explicit classes, downstream modules no longer need to guess whether a record belongs to protein, solvent, ligand, or metal.

### Reconstruction

The filler stage is where the pipeline stops being pure cleanup and becomes reconstruction-aware. Missing internal regions are classified and then handled explicitly through reconstruction candidates such as MODELLER and, in larger-gap cases, AlphaFold-derived models.

### Current chemistry-aware normalization

The current main runner uses GROMACS `pdb2gmx` for protonation and force-field-compatible protein preparation.

This step is currently deliberately simpler than the long-term chemistry plan. It allows the project to stabilize representative-unit selection, gap behavior, filler routing, and prepared-structure assembly before reintroducing more detailed AMBER naming, true termini handling, and internal fragment capping.

### Prepared structure assembly

The final prepared output is not just “the latest available protein PDB”. It is a deliberate system handoff file that combines the prepared protein with waters, ligands, and metals from the representative unit’s component split.

This avoids reintroducing ligands, cofactors, or metals from discarded duplicate chains.

### Metal branch

The metal branch is planned as a later parametrization layer. The intended direction is to detect metal identity, local coordination environment, geometry class, and quality/caution flags, then prepare downstream Gaussian / MCPB-style parameter derivation inputs.

---

## Prepared-structure output logic

For each protein directory:

```text
data/proteins/<PDB_ID>/
```

FRUTON writes final prepared structures according to retained variant logic.

### Case 1: one best-available path

```text
prepared/<PDB_ID>.pdb
```

### Case 2: retained gaps variant

```text
prepared/gaps/<PDB_ID>.pdb
```

### Case 3: retained reconstructed complete variant

```text
prepared/complete/<PDB_ID>.pdb
```

In all cases, the final prepared structure contains components derived from the representative unit:

- protonated protein variant
- representative-unit waters if present
- representative-unit ligands if present
- representative-unit metals if present

---

## Pipeline state

FRUTON writes state to:

```text
data/proteins/pipeline.json
data/proteins/pipeline.xlsx
```

### `pipeline.json`

The primary machine-readable state store.

### `pipeline.xlsx`

A human-readable overview for rapid inspection.

The current state model keeps fields useful for:

- orchestration
- JSON persistence
- XLSX export
- debugging failures
- tracking final prepared outputs

### Representative state fields

#### Identity and path fields

- `pdb_id`
- `range`
- `pdb_directory`
- `fasta_directory`
- `alignment_directory`
- `components_directory`
- `prepared_directory`
- `uniprot_id`

#### Structural summary

- `n_gaps`
- `gap_sizes`
- `has_gaps`
- `has_metals`
- `has_ligands`
- `has_nonstandard_residues`

#### Filler and prepared-variant fields

- `filler.directory`
- `filler.model_path`
- `filler.model_source`
- `filler.status`
- `variant_policy`
- `retain_gaps_variant`
- `retain_modeller_variant`
- `retain_alphafold_variant`
- `available_models`

#### Protonation and prepared-structure fields

- `protonation.status`
- `protonation.input_source`
- `protonation.input_path`
- `protonation.output_path`
- `prepared_structure.status`
- `prepared_structure.variant`
- `prepared_structure.protein_input_path`
- `prepared_structure.output_path`
- `prepared_structure.gaps_output_path`
- `prepared_structure.modeller_output_path`
- `prepared_structure.alphafold_output_path`

#### Step statuses

- `pdb_sync_done`
- `fasta_files_done`
- `sequence_alignment_done`
- `insertion_codes_done`
- `filler.status`
- `protonation.status`
- `prepared_structure.status`

---

## XLSX export

The Excel export is handled by:

```text
src/stack_protein_preparation/pipeline_xlsx.py
```

Current export behavior:

- writes a compact overview first
- keeps `pdb_id`, UniProt, range, and key decision/status fields easy to inspect
- exports detailed per-module information in follow-up worksheets
- colors status and key structural flags for rapid inspection

Current structural flag coloring:

- `has_ligands = no` -> normal green
- `has_ligands = yes` -> light green
- `has_nonstandard_residues = no` -> green
- `has_nonstandard_residues = yes` -> red
- `has_metals = no` -> green
- `has_metals = yes` -> red

Current status palette:

- `success` -> green
- `warning` -> yellow
- `required` -> red
- `skipped` -> grey
- `failed` -> orange-red

---

## Current strengths

FRUTON already provides a solid base for:

- alignment-aware protein preparation
- reproducible intermediate-file generation
- explicit step-by-step state tracking
- representative-unit selection
- chain-specific structural reasoning
- ligand retention per representative unit
- gap-aware variant retention
- MODELLER / AlphaFold-oriented reconstruction paths
- GROMACS-based protonation
- final prepared-structure assembly
- future branching toward metal-site preparation

---

## Current limitations

FRUTON is still actively evolving, especially in the later stages.

### Areas still being refined

- robust multi-chain filler handling
- filler branch ordering for MODELLER versus AlphaFold candidates
- terminal placeholder `X` handling in filler validation
- variant-specific downstream handling when several gap-related variants are retained
- reintroduction of AMBER renaming, AMBER termini, and internal capping after the current representative-unit logic is stable
- metal-selection and metal-contact parsing
- atom serial renumbering in merged metal-preparation PDBs
- automatic MCPB / Gaussian preparation

### Known practical issue

Old generated files can mislead interpretation when the runner logic changes. After changing representative-unit or component-split behavior, clean stale per-protein outputs before rerunning.

Example:

```bash
rm -rf data/proteins/<PDB_ID>/components
rm -rf data/proteins/<PDB_ID>/prepared
rm -rf data/proteins/<PDB_ID>/protonation
rm -f data/proteins/pipeline.json data/proteins/pipeline.xlsx
```

---

## Recommended next steps

1. **Validate representative-unit ligand handling across more PDB examples**
2. **Improve multi-chain filler handling**
3. **Repair filler branch ordering so AlphaFold candidates do not fail through MODELLER-only validation**
4. **Relax filler validation for terminal placeholder `X` residues**
5. **Make retained variant outputs fully variant-specific**
6. **Reintroduce AMBER renaming / termini / internal capping once the current route is stable**
7. **Improve metal branch logic and structured contact parsing**
8. **Build MCPB / Gaussian preparation automatically**

---

## Running FRUTON

Typical full pipeline run:

```bash
pixi run python scripts/fruton.py
```

For debugging a changed representative-unit or component-split rule, clean the relevant protein outputs first:

```bash
rm -rf data/proteins/<PDB_ID>/components
rm -rf data/proteins/<PDB_ID>/prepared
rm -rf data/proteins/<PDB_ID>/protonation
rm -f data/proteins/pipeline.json data/proteins/pipeline.xlsx
pixi run python scripts/fruton.py
```

---

## Quick representative-unit checks

Check which ligands survived representative-unit selection:

```bash
grep '^HETATM' data/proteins/<PDB_ID>/components/<PDB_ID>_representative_unit.pdb \
  | awk '{print substr($0,18,3), substr($0,22,1)}' \
  | sort \
  | uniq -c
```

Check the final split ligand component:

```bash
grep '^HETATM' data/proteins/<PDB_ID>/components/<PDB_ID>_ligand.pdb \
  | awk '{print substr($0,18,3), substr($0,22,1)}' \
  | sort \
  | uniq -c
```

For a case like duplicated ligand-bearing chains A/B where only A is selected, the split ligand component should show only the ligand on A.

---

## Main source files

### Orchestration

```text
scripts/fruton.py
```

### State and export

```text
src/stack_protein_preparation/pipeline_state.py
src/stack_protein_preparation/pipeline_table.py
src/stack_protein_preparation/pipeline_xlsx.py
```

### Sequence and alignment

```text
src/stack_protein_preparation/fasta_files.py
src/stack_protein_preparation/sequence_alignment.py
```

### Structure processing

```text
src/stack_protein_preparation/insertion_codes.py
src/stack_protein_preparation/monomer.py
src/stack_protein_preparation/pdb_components.py
src/stack_protein_preparation/gaps.py
src/stack_protein_preparation/fragment_split.py
src/stack_protein_preparation/filler.py
```

### Current chemistry and assembly

```text
src/stack_protein_preparation/protonation.py
src/stack_protein_preparation/prepared_structure.py
```

### Planned / inactive chemistry extensions

```text
src/stack_protein_preparation/amber_renaming.py
src/stack_protein_preparation/terminus.py
src/stack_protein_preparation/internal_capping.py
```

### Metal branch

```text
src/stack_protein_preparation/metall_params.py
```

---

## One-line summary

**FRUTON is a filesystem-explicit, state-driven protein-preparation framework that connects sequence-aware cleanup, representative-unit selection, gap-aware reconstruction, GROMACS protonation, prepared-structure assembly, and future metal-site preparation into one coherent workflow.**

