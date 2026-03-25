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
  A modular protein-preparation pipeline for alignment-aware cleanup, gap handling, protonation, AMBER normalization, termini preparation, internal fragment capping, final prepared-structure assembly, and metal-site preparation.
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
- missing internal residues
- UniProt mismatches
- ligands, water, and metals
- protonation ambiguities
- AMBER naming requirements
- cropped residue ranges that create artificial termini
- disconnected internal fragments that need explicit chemical treatment

The project is built around a simple principle:

> **Every important transformation should leave visible evidence on disk.**

That means:

- intermediate files are intentionally written
- alignments and mapping files are preserved
- state is stored explicitly in JSON and XLSX
- downstream steps consume well-defined upstream outputs
- failures are meant to be inspectable, not mysterious

---

## Why FRUTON exists

Protein preparation for modelling or MD is rarely a single command. It is usually a chain of assumptions:

```text
structure → sequence interpretation → cleanup → reconstruction → chemistry → final prepared assembly
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
[ COMPONENT SPLIT ]
     │
     ▼
[ GAP DETECTION ]
     │
     ▼
[ FILLER: MODELLER / ALPHAFOLD ]
     │
     ▼
[ PROTONATION ]
     │
     ▼
[ AMBER RENAMING ]
     │
     ▼
[ AMBER TERMINI ]
     │
     ▼
[ INTERNAL ACE/NME CAPPING ]
     │
     ▼
[ PREPARED STRUCTURE ]
     │
     └──────────────────────────────► [ METAL PREPARATION ]
```

---

## Current capabilities

| Stage | Purpose | Current output |
|---|---|---|
| **PDB sync** | synchronize dataset scope from CSV | protein directory tree |
| **FASTA generation** | generate structure- and reference-derived sequences | `SEQRES`, `ATOM`, `UniProt` FASTA files |
| **Sequence alignment** | compare PDB chains against UniProt | aligned FASTA, TSV mappings, PNGs |
| **Insertion cleanup** | remove insertion-code ambiguity | cleaned PDB |
| **Component split** | separate structural classes | protein / water / ligand / metal files |
| **Gap detection** | quantify missing internal regions | `n_gaps`, `gap_sizes`, `has_gaps` |
| **Filler** | rebuild missing regions when possible | MODELLER or AlphaFold-derived model |
| **Protonation** | add hydrogens / assign protonation | protonated protein PDB |
| **AMBER renaming** | assign AMBER-compatible residue naming | AMBER-style protein PDB |
| **AMBER termini** | convert true chain ends to AMBER terminal residues | `*_protein_amber_termini.pdb` |
| **Internal capping** | cap internal disconnected fragments with ACE/NME | `*_protein_internal_capped.pdb` |
| **Prepared structure** | assemble final system for downstream use | `prepared/.../<PDB_ID>.pdb` |
| **Metal preparation** | prepare metal systems for later parametrization | combined PDB, Chimera script, contacts output |

---

## Important chemical distinction

FRUTON now treats two related but different situations explicitly.

### 1. True chain termini

Real chain starts and chain ends are converted into **AMBER terminal residues**.

Examples:

- `ALA -> NALA`
- `GLY -> CGLY`

This is handled by:

```text
src/stack_protein_preparation/terminus.py
```

### 2. Internal missing segments / disconnected internal fragments

When a cropped or incomplete structure produces disconnected internal fragments inside one chain, these are **not** treated as true free chain termini.

Instead, the internal break boundaries are capped with:

- `NME` on the **left fragment end**
- `ACE` on the **right fragment start**

This is handled by:

```text
src/stack_protein_preparation/cap.py
```

So the rule is:

- **true chain ends** -> AMBER terminal residues
- **internal breaks** -> ACE/NME capping

---

## Project architecture

FRUTON is easiest to understand as a stack of cooperating layers:

| Layer | Role |
|---|---|
| **Sequence layer** | FASTA generation, UniProt matching, alignment TSVs |
| **Structure layer** | insertion cleanup, component split, gap detection |
| **Reconstruction layer** | MODELLER / AlphaFold-based filling |
| **Chemistry layer** | protonation, AMBER renaming, terminal normalization, internal capping |
| **Prepared-assembly layer** | final merge of prepared protein with waters / ligands / metals |
| **Metal branch** | preparation for later MCPB / Gaussian-style workflows |

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
│       ├── pdb_components.py
│       ├── gaps.py
│       ├── filler.py
│       ├── protonation.py
│       ├── amber_renaming.py
│       ├── terminus.py
│       ├── cap.py
│       ├── prepared_structure.py
│       ├── metall_params.py
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
│   ├── <PDB_ID>_protein.pdb
│   ├── <PDB_ID>_water.pdb
│   ├── <PDB_ID>_ligand.pdb
│   ├── <PDB_ID>_metal.pdb
│   ├── <PDB_ID>_proteinH.pdb
│   ├── <PDB_ID>_protein_as_Amber.pdb
│   ├── <PDB_ID>_protein_amber_termini.pdb
│   └── <PDB_ID>_protein_internal_capped.pdb
├── prepared/
│   ├── <PDB_ID>.pdb
│   ├── gaps/
│   │   └── <PDB_ID>.pdb
│   └── complete/
│       └── <PDB_ID>.pdb
└── metall_params/
    ├── tmp_param.pdb
    ├── metal_only.pdb
    ├── chimera_contacts.py
    ├── contacts.data
    └── chimera_run.log
```

Not every protein will contain every file. The `prepared/` layout depends on whether the system has gaps and whether a completed model variant is available.

---

## Core outputs and what they mean

| File | Interpretation |
|---|---|
| `<PDB_ID>_delins.pdb` | insertion-cleaned working PDB |
| `PDB-<PDB_ID>-SEQRES.fasta` | sequence derived from SEQRES records |
| `PDB-<PDB_ID>-ATOM.fasta` | sequence derived from observed ATOM records |
| `UniProt_<UNIPROT_ID>.fasta` | external UniProt reference |
| `*_vs_UniProt.aln.fasta` | chain-specific alignment output |
| `*_vs_UniProt.mapping.tsv` | alignment mapping table |
| `<PDB_ID>_protein.pdb` | structural protein component |
| `<PDB_ID>_proteinH.pdb` | protonated protein |
| `<PDB_ID>_protein_as_Amber.pdb` | AMBER-renamed protein |
| `<PDB_ID>_protein_amber_termini.pdb` | AMBER-compatible true terminal residues |
| `<PDB_ID>_protein_internal_capped.pdb` | internally ACE/NME-capped prepared protein |
| `prepared/<PDB_ID>.pdb` | final prepared structure without gap subdirectory |
| `prepared/gaps/<PDB_ID>.pdb` | final prepared structure for the gapped variant |
| `prepared/complete/<PDB_ID>.pdb` | final prepared structure for the completed variant |
| `<PDB_ID>_metal.pdb` | isolated metal component |
| `metall_params/tmp_param.pdb` | merged structural input for metal analysis |
| `metall_params/contacts.data` | Chimera contact/clash report |

---

## Why the individual stages matter

### Sequence-aware preparation

FRUTON treats sequence logic as core infrastructure, not as optional annotation. UniProt alignment is used to:

- compare observed structure against reference biology
- support chain-specific reasoning
- identify truncation and missing segments
- provide a stable reference layer for reconstruction

### Structural cleanup

Insertion-code cleanup and component splitting simplify later reasoning. Once the structure is decomposed into explicit classes, downstream modules no longer need to guess whether a record belongs to protein, solvent, ligand, or metal.

### Reconstruction

The filler stage is where the pipeline stops being pure cleanup and becomes reconstruction-aware. Missing internal regions are classified and then handled explicitly with MODELLER or AlphaFold fallback.

### Chemistry-aware normalization

Protonation and AMBER renaming remain separate steps because they solve related, but distinct, problems:

- protonation adds chemically meaningful hydrogens
- AMBER renaming encodes residue state in a force-field-compatible way

FRUTON now extends this chemistry layer with two additional explicit transformations:

- **AMBER termini** for true chain ends
- **internal ACE/NME capping** for disconnected internal fragments

### Prepared structure assembly

The final prepared output is not just “the latest available protein PDB”. It is a deliberate system handoff file that combines the chemically prepared protein with crystal waters, ligands, and metals in a predictable directory structure.

### Metal branch

The metal branch prepares later parametrization work by creating a combined system, isolating the metal component, and running first-pass local contact analysis in Chimera.

---

## Prepared-structure output logic

For each protein directory:

```text
data/proteins/<PDB_ID>/
```

FRUTON writes the final prepared structure as follows.

### Case 1: no gaps

```text
prepared/<PDB_ID>.pdb
```

### Case 2: gaps remain

```text
prepared/gaps/<PDB_ID>.pdb
```

### Case 3: completed variant available

```text
prepared/complete/<PDB_ID>.pdb
```

In all cases, the final prepared structure contains:

- the protonated, AMBER-normalized, terminally prepared protein
- crystal waters if present
- ligands if present
- metals if present

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

The current state model was reduced to keep only fields that are actually useful for:

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

#### Step outputs

- `filler.model_path`
- `protonation.output_path`
- `amber_renaming.output_path`
- `amber_termini.output_path`
- `internal_capping.output_path`
- `prepared_structure.output_path`

#### Step statuses

- `pdb_sync_done`
- `fasta_files_done`
- `sequence_alignment_done`
- `insertion_codes_done`
- `filler.status`
- `protonation.status`
- `amber_renaming.status`
- `amber_termini.status`
- `internal_capping.status`
- `prepared_structure.status`

---

## XLSX export

The Excel export is handled by:

```text
src/stack_protein_preparation/pipeline_xlsx.py
```

Current export behavior:

- exports only non-empty columns
- keeps `pdb_id` always
- colors only real status columns
- uses the current status palette:
  - `success` -> green
  - `warning` -> yellow
  - `required` -> red
  - `skipped` -> grey
  - `failed` -> orange-red

This keeps the sheet readable while still preserving the most important pipeline information.

---

## Current strengths

FRUTON already provides a solid base for:

- alignment-aware protein preparation
- reproducible intermediate-file generation
- explicit step-by-step state tracking
- modular orchestration
- chain-specific structural reasoning
- chemistry-aware normalization
- explicit handling of true termini versus internal fragment breaks
- final prepared-structure assembly
- branching toward metal-site preparation

---

## Current limitations

FRUTON is still actively evolving, especially in the later stages.

### Areas still being refined

- more robust multi-chain filler handling
- cleaner dual-variant support for `gaps` and `complete` all the way through chemistry steps
- Chimera metal-selection logic
- atom serial renumbering in metal-preparation merged PDBs
- parsing of metal-contact output into structured data

### Known practical issue

Merged files such as `metall_params/tmp_param.pdb` may still contain duplicate atom serial numbers because concatenation is currently used instead of full renumbering. This usually produces warnings rather than immediate failure, but it should eventually be cleaned up.

---

## Recommended next steps

1. **Make late-stage outputs fully variant-specific for `gaps` and `complete`**
2. **Improve Chimera metal-selection logic**
3. **Renumber merged atom serials in metal-preparation files**
4. **Parse `contacts.data` into structured output**
5. **Build MCPB input automatically**
6. **Extend the metal branch toward Gaussian preparation**
7. **Improve rerun granularity for late-stage chemistry steps**

---

## Running FRUTON

Typical full pipeline run:

```bash
pixi run python scripts/fruton.py
```

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
src/stack_protein_preparation/pdb_components.py
src/stack_protein_preparation/gaps.py
src/stack_protein_preparation/filler.py
```

### Chemistry and assembly

```text
src/stack_protein_preparation/protonation.py
src/stack_protein_preparation/amber_renaming.py
src/stack_protein_preparation/terminus.py
src/stack_protein_preparation/cap.py
src/stack_protein_preparation/prepared_structure.py
```

### Metal branch

```text
src/stack_protein_preparation/metall_params.py
```

---

## One-line summary

**FRUTON is a filesystem-explicit, state-driven protein-preparation framework that connects sequence-aware cleanup, structure reconstruction, chemistry-aware normalization, terminal handling, internal fragment capping, prepared-structure assembly, and metal-site preparation into one coherent workflow.**
