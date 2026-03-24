# `stack_protein_prep` — Pipeline Overview

## Goal

`stack_protein_prep` is a modular protein-preparation pipeline designed to take raw protein structure inputs and transform them into progressively more useful, reproducible, and downstream-ready files.

The long-term goal is to support a full preparation workflow for protein systems that may contain:

- missing residues / gaps
- insertion codes
- ligands
- water
- metals
- protonation-state ambiguities
- AMBER naming issues
- final numbering restoration requirements
- later metal parametrization and MD preparation

This repository is built as a **stepwise pipeline with explicit state tracking**, so that each protein can be processed reproducibly and partially rerun if needed.

---

# High-Level Design

The pipeline operates on a directory structure like:

```text
data/proteins/<PDB_ID>/
```

and stores global state in:

```text
data/proteins/pipeline.json
data/proteins/pipeline.xlsx
```

Each protein gets its own record in the pipeline state.

The pipeline is designed around these principles:

- each step produces explicit files
- each step updates explicit state fields
- failures should be local to a protein where possible
- downstream steps rely on well-defined upstream outputs
- files are preferred over implicit in-memory state
- the final pipeline should be easy to debug from filesystem + logs alone

---

# Core Directory Structure

## Global level

```text
data/proteins/
├── pdb_ids.csv
├── pipeline.json
├── pipeline.xlsx
├── 1IOO/
├── 2AFX/
├── 2P1T/
└── ...
```

---

## Per-protein level

Typical structure:

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
│       ├── ATOM_chain_A_vs_UniProt.png
│       └── <PDB_ID>_finalize_numbering.tsv
├── components/
│   ├── <PDB_ID>_protein.pdb
│   ├── <PDB_ID>_water.pdb
│   ├── <PDB_ID>_ligand.pdb
│   ├── <PDB_ID>_metal.pdb
│   ├── <PDB_ID>_proteinH.pdb
│   ├── <PDB_ID>_protein_as_Amber.pdb
│   └── <PDB_ID>_protein_final.pdb
├── MODELLER/
├── metall_params/
│   ├── tmp_param.pdb
│   ├── metal_only.pdb
│   ├── chimera_contacts.py
│   ├── contacts.data
│   └── chimera_run.log
└── ...
```

---

# Main Pipeline Script

Main entry point:

```text
/home/grheco/repositorios/stack_protein_prep/scripts/run_pipeline.py
```

This script orchestrates the currently implemented steps.

---

# Implemented Pipeline Steps

## Step 1 — `pdb_sync`

### Purpose

Synchronize the input CSV with the actual protein directories.

### Inputs

- `data/proteins/pdb_ids.csv`

### Outputs

- creates / updates protein directories
- ensures the working data tree is aligned with CSV contents

### Main idea

This step defines **which proteins are in scope** for the current run.

---

## Step 2 — Read input CSV

### Purpose

Load the CSV records that define the active protein set.

### Typical fields

At least:

- `pdb_id`
- `range` (optional residue-range restriction)

### Effect

These values are used to construct the initial pipeline records.

---

## Step 3 — Load existing pipeline state

### Purpose

Load previously saved state from:

```text
data/proteins/pipeline.json
```

### Effect

Allows continuation and partial reruns without rebuilding everything from scratch.

---

## Step 4 — Build and merge pipeline state

### Purpose

Construct fresh pipeline records from CSV and merge them with existing records.

### Output fields initialized / maintained

Examples:

- `pdb_id`
- `range`
- `pdb_directory`
- `fasta_directory`
- `alignment_directory`
- `components_directory`

---

## Step 5 — `fasta_files`

### Purpose

Generate the FASTA inputs required for alignment.

### Inputs

Typically based on raw or cleaned PDB structure.

### Outputs

```text
fasta/PDB-<PDB_ID>-SEQRES.fasta
fasta/PDB-<PDB_ID>-ATOM.fasta
fasta/UniProt_<UNIPROT_ID>.fasta
```

### Additional behavior

If no UniProt FASTA is available locally, the module can try to:

1. resolve a UniProt accession via RCSB
2. download UniProt FASTA

### Stored state

- `fasta_files_done`
- `uniprot_id`

---

## Step 6 — `sequence_alignment`

### Purpose

Align chain-specific PDB-derived sequences against UniProt.

### Comparisons currently created

- each `SEQRES` chain vs UniProt
- each `ATOM` chain vs UniProt

### Outputs per chain

```text
SEQRES_chain_A_vs_UniProt.input.fasta
SEQRES_chain_A_vs_UniProt.aln.fasta
SEQRES_chain_A_vs_UniProt.mapping.tsv
SEQRES_chain_A_vs_UniProt.png

ATOM_chain_A_vs_UniProt.input.fasta
ATOM_chain_A_vs_UniProt.aln.fasta
ATOM_chain_A_vs_UniProt.mapping.tsv
ATOM_chain_A_vs_UniProt.png
```

### Important note

The mapping TSV written here is an **alignment mapping**, not automatically the final numbering source for the final prepared model.

### Stored state

- `sequence_alignment_done`

---

## Step 7 — `insertion_codes`

### Purpose

Handle insertion codes in PDB residue numbering and create a cleaned structure.

### Inputs

- raw input PDB

### Output

```text
<PDB_ID>_delins.pdb
```

### Why this matters

Insertion codes often cause downstream mismatch problems in:

- sequence extraction
- residue mapping
- gap detection
- MODELLER alignment interpretation

### Stored state

- `insertion_codes_done`

---

## Step 8 — `component_split`

### Purpose

Split the cleaned PDB into structural components.

### Inputs

```text
<PDB_ID>_delins.pdb
```

### Outputs

```text
components/<PDB_ID>_protein.pdb
components/<PDB_ID>_water.pdb
components/<PDB_ID>_ligand.pdb
components/<PDB_ID>_metal.pdb
```

### Current component summary logic

The pipeline tracks booleans / statuses for:

- `has_metals`
- `has_ligands`
- `has_nonstandard_residues`

### Stored state

- `components_directory`
- `has_metals`
- `has_ligands`
- `has_nonstandard_residues`

---

## Step 9 — `gap_detection`

### Purpose

Detect structural gaps in the protein component.

### Input

```text
components/<PDB_ID>_protein.pdb
```

### Outputs in state

- `n_gaps`
- `gap_sizes`

### Example semantics

- `n_gaps = 0`
- `gap_sizes = "none"`

or

- `n_gaps = 2`
- `gap_sizes = "6|6"`

### Why this matters

Gap size classification affects whether:

- filler is skipped
- MODELLER is run
- AlphaFold fallback is preferred

---

## Step 10 — `filler`

### Purpose

Fill internal missing regions using either:

- MODELLER
- AlphaFold fallback

depending on the gap classification.

### Current gap logic in `filler.py`

- internal gaps only
- `1–5` residues → green
- `6–8` residues → yellow
- `>8` residues → AlphaFold candidate

### Inputs

- alignment files
- protein template PDB
- UniProt information
- optional residue range

### Outputs

Typical:

```text
fasta/alignments/filler/A/<PDB_ID>_protein_mod.pdb
```

### Current behavior

- no internal gaps → skip
- large internal gap → AlphaFold fallback if possible
- moderate / small internal gaps → MODELLER

### Stored state

- `filler_directory`
- `filler_model_path`
- `filler_status`

---

## Step 11 — `protonation`

### Purpose

Create a protonated protein structure.

### Input preference

- filler model if present
- otherwise protein component

### Output

```text
components/<PDB_ID>_proteinH.pdb
```

### Stored state

- `protonation.status`
- `protonation.input_source`
- `protonation.input_path`
- `protonation.output_path`
- `protonation.ph`
- `protonation.input_atom_count`
- `protonation.output_atom_count`
- `protonation.atom_count_increased`

---

## Step 12 — `amber_renaming`

### Purpose

Rename protonated residues into AMBER-compatible residue names.

### Input

```text
components/<PDB_ID>_proteinH.pdb
```

### Output

```text
components/<PDB_ID>_protein_as_Amber.pdb
```

### Examples of tracked renames

- HIS → HID / HIE / HIP
- ASP → ASH
- GLU → GLH
- CYS → CYM / CYX

### Stored state

- `amber_renaming.status`
- `amber_renaming.input_path`
- `amber_renaming.output_path`
- `amber_renaming.his_to_hid`
- `amber_renaming.his_to_hie`
- `amber_renaming.his_to_hip`
- `amber_renaming.asp_to_ash`
- `amber_renaming.glu_to_glh`
- `amber_renaming.cys_to_cym`
- `amber_renaming.cys_to_cyx`

---

## Step 13 — `finalize_protein`

### Purpose

Create the final prepared protein structure with final numbering.

### Current output

```text
components/<PDB_ID>_protein_final.pdb
```

### Current numbering logic

This step uses a **dedicated finalize-numbering TSV**, not directly the older alignment mapping file.

Typical file:

```text
fasta/alignments/<PDB_ID>_finalize_numbering.tsv
```

### Why this step is separate

The alignment mapping TSV from `sequence_alignment.py` describes alignment columns, but the final model may include:

- original residues
- modeled residues
- gap-filled residues
- AlphaFold/MODELLER-generated residues

Therefore a **final explicit numbering table** is more reliable.

### State fields currently reused

The pipeline currently stores this under the older `numbering_restore.*` fields:

- `numbering_restore.status`
- `numbering_restore.input_path`
- `numbering_restore.output_path`
- `numbering_restore.mapping_path`
- `numbering_restore.source`
- `numbering_restore.renumbered_atoms`
- `numbering_restore.renumbered_residues`
- `numbering_restore.message`

---

## Step 14 — `metall_params`

### Purpose

Prepare metal-containing systems for downstream metal parametrization.

### Inputs

Usually:

- `components/<PDB_ID>_protein_final.pdb`
- `components/<PDB_ID>_water.pdb`
- `components/<PDB_ID>_ligand.pdb`
- `components/<PDB_ID>_metal.pdb`

### Outputs

```text
metall_params/tmp_param.pdb
metall_params/metal_only.pdb
metall_params/chimera_contacts.py
metall_params/contacts.data
metall_params/chimera_run.log
```

### Current status

This step is still in active debugging.

### Current known issues addressed recently

- file naming mismatch: `metal.pdb` vs `metals.pdb`
- Chimera executable resolution
- incorrect Chimera selection logic (`#0 test other`)
- need for more explicit debug logging

### Current intended Chimera logic

- open full system as one model
- open metal-only file as another model
- run contact search from metal to system

---

## Step 15 — Save pipeline JSON

### Purpose

Persist the entire pipeline state to:

```text
data/proteins/pipeline.json
```

This is the main machine-readable state file.

---

## Step 16 — Write pipeline XLSX

### Purpose

Write a human-readable Excel summary to:

```text
data/proteins/pipeline.xlsx
```

Useful for quick review of:

- status overview
- paths
- gaps
- protonation
- AMBER renaming
- finalize status
- metal-prep status

---

# State Model

The pipeline uses a **flat record per protein**.

## Example record sections

### Identity and directory fields

- `pdb_id`
- `range`
- `pdb_directory`
- `fasta_directory`
- `alignment_directory`
- `components_directory`
- `uniprot_id`

### Structural summary

- `n_gaps`
- `gap_sizes`
- `has_metals`
- `has_ligands`
- `has_nonstandard_residues`

### Step status fields

- `pdb_sync_done`
- `fasta_files_done`
- `sequence_alignment_done`
- `insertion_codes_done`
- `filler_status`
- `protonation.status`
- `amber_renaming.status`
- `numbering_restore.status`
- `metall_params.status` (currently stored as flat string keys if present)

---

# File Semantics by Preparation Level

## `*_protein.pdb`

Basic protein-only structural component after component split.

## `*_proteinH.pdb`

Protonated protein.

## `*_protein_as_Amber.pdb`

Protonated protein with AMBER naming.

## `*_protein_final.pdb`

Final prepared protein after final numbering logic.

This is currently the preferred protein input for metal preparation.

---

# Current Metal Preparation Logic

## Why a separate `metall_params` directory?

Because metal parametrization is not just another cleanup step. It is the beginning of a separate downstream branch that may later create:

- MCPB input
- Gaussian input
- force-field files
- topology fragments

So it should stay isolated.

---

## Current assembled file

`tmp_param.pdb` is currently built by concatenating:

1. final protein
2. water
3. ligand
4. metal

This is intentionally simple and debuggable.

---

## Known caveat

The merged PDB may contain duplicate atom serial numbers because it is created by concatenation rather than full renumbering.

This can trigger Chimera warnings but is not necessarily fatal.

---

# Finalize / Numbering Philosophy

The project originally used alignment mapping files directly, but this proved insufficient once modelled residues were introduced.

The current direction is:

- keep alignment TSVs for sequence logic and diagnostics
- create a separate finalize-numbering TSV for final structural numbering
- let `finalize_protein.py` apply only that explicit numbering table

This separation is important because:

- alignment logic and numbering logic are not identical
- modeled residues need explicit policy
- final output should be deterministic

---

# Current Debugging Strategy

The project increasingly favors **explicit, verbose, step-local debugging**.

Especially for difficult steps like:

- filler
- finalize_protein
- metall_params

the preferred approach is:

- print paths
- print selected inputs
- print booleans for each decision step
- write logs to files
- store exact paths in pipeline state

This is consistent with the repository’s overall design philosophy:
**debug from filesystem + state, not from guesswork.**

---

# Current Status of the Pipeline

## Stable / mostly established

- PDB sync
- FASTA generation
- sequence alignment
- insertion-code cleanup
- component split
- gap detection
- protonation
- AMBER renaming

## Working but still evolving

- filler
- finalize_protein / finalize TSV logic
- metall_params / Chimera logic

---

# Typical End-to-End Flow for One Protein

A simplified progression for one protein looks like this:

```text
raw PDB
→ FASTA generation
→ PDB vs UniProt alignment
→ insertion-code cleanup
→ component split
→ gap detection
→ filler (MODELLER / AlphaFold if needed)
→ protonation
→ AMBER residue renaming
→ final numbering / final protein output
→ metal preparation (if metal present)
```

---

# Example Final Important Files

For a metal-free system:

```text
components/<PDB_ID>_protein_final.pdb
```

For a metal-containing system:

```text
components/<PDB_ID>_protein_final.pdb
metall_params/tmp_param.pdb
metall_params/contacts.data
```

---

# Naming Conventions Summary

| Kind | Example |
|---|---|
| cleaned PDB | `2AFX_delins.pdb` |
| protein component | `2AFX_protein.pdb` |
| protonated protein | `2AFX_proteinH.pdb` |
| AMBER-renamed protein | `2AFX_protein_as_Amber.pdb` |
| final protein | `2AFX_protein_final.pdb` |
| ligand component | `2AFX_ligand.pdb` |
| water component | `2AFX_water.pdb` |
| metal component | `2AFX_metal.pdb` |
| finalize numbering TSV | `2AFX_finalize_numbering.tsv` |

---

# Recommended Next Improvements

## 1. Formalize finalize TSV generation
The builder should be fully explicit and chain-aware.

## 2. Improve metal selection logic
Use metal-only model or explicit atom selection instead of whole-model `#0 test other`.

## 3. Renumber merged PDB atom serials
This would reduce Chimera warnings.

## 4. Parse `contacts.data`
Store structured contacts in JSON/TSV rather than raw text only.

## 5. Add MCPB input builder
Once metal contacts are robust, build the MCPB input step automatically.

## 6. Improve rerun granularity
Allow easy reruns from:
- filler onward
- finalize onward
- metall_params onward

---

# Repository Philosophy

This repository is not built as one monolithic black box.

It is built as a **traceable scientific workflow** where:

- every important step writes files
- every decision leaves traces
- intermediate files are intentionally inspectable
- pipeline state is explicit
- debugging is expected and supported

That is especially important for structural biology / MD preparation, where silent assumptions create bad downstream systems.

---

# Main Files of Interest

## Orchestration

```text
scripts/run_pipeline.py
```

## State handling

```text
src/stack_protein_preparation/pipeline_state.py
src/stack_protein_preparation/pipeline_table.py
src/stack_protein_preparation/pipeline_xlsx.py
```

## FASTA / alignment

```text
src/stack_protein_preparation/fasta_files.py
src/stack_protein_preparation/sequence_alignment.py
```

## Structural cleanup

```text
src/stack_protein_preparation/insertion_codes.py
src/stack_protein_preparation/pdb_components.py
src/stack_protein_preparation/gaps.py
```

## Model completion

```text
src/stack_protein_preparation/filler.py
```

## Chemistry-related preparation

```text
src/stack_protein_preparation/protonation.py
src/stack_protein_preparation/amber_renaming.py
src/stack_protein_preparation/finalize_protein.py
src/stack_protein_preparation/metall_params.py
```

---

# Summary

`stack_protein_prep` is a modular structural-preparation pipeline that currently:

- prepares protein systems step by step
- keeps explicit state in JSON/XLSX
- handles alignment, cleanup, gap filling, protonation, renaming, finalization
- is expanding toward metal parametrization support

The pipeline is already useful as a structured preparation framework, but its later steps — especially:

- `filler`
- `finalize_protein`
- `metall_params`

are still actively being refined.

The current architecture is solid because it already separates:

- structural cleanup
- sequence logic
- chemical naming
- final output generation
- metal-preparation branching

That separation will make future automation much easier.
