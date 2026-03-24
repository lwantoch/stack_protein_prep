# FRUTON

```text
███████╗██████╗ ██╗   ██╗████████╗ ██████╗ ███╗   ██╗
██╔════╝██╔══██╗██║   ██║╚══██╔══╝██╔═══██╗████╗  ██║
█████╗  ██████╔╝██║   ██║   ██║   ██║   ██║██╔██╗ ██║
██╔══╝  ██╔══██╗██║   ██║   ██║   ██║   ██║██║╚██╗██║
██║     ██║  ██║╚██████╔╝   ██║   ╚██████╔╝██║ ╚████║
╚═╝     ╚═╝  ╚═╝ ╚═════╝    ╚═╝    ╚═════╝ ╚═╝  ╚═══╝
```

> **Framework for Reconstruction, UniProt alignment, and Topology-Oriented protein Normalization**

A modular protein-preparation pipeline for **alignment-aware cleanup**, **gap handling**, **protonation**, **AMBER renaming**, **final numbering**, and **metal-site preparation**.

---

## Why FRUTON exists

Preparing a protein structure for downstream work is rarely a single cleanup step.  
Real systems are messy. They contain:

- insertion codes
- missing internal residues
- UniProt mismatches
- ligands
- waters
- metals
- protonation ambiguities
- force-field naming issues
- numbering inconsistencies after modelling

FRUTON exists to turn that mess into a **traceable**, **reproducible**, and **inspectable** workflow.

Instead of hiding everything inside one opaque script, FRUTON treats protein preparation as a chain of explicit transformations:

```text
raw PDB
→ FASTA generation
→ UniProt alignment
→ insertion cleanup
→ component split
→ gap detection
→ MODELLER / AlphaFold filler
→ protonation
→ AMBER renaming
→ final numbering
→ metal preparation
```

Each important step produces files, updates state, and leaves enough evidence on disk to debug what happened.

---

## What makes FRUTON different

FRUTON is built around a few strong design ideas:

| Principle | Meaning in practice |
|---|---|
| **Filesystem-first** | Intermediate files are intentionally written and kept. |
| **State-driven** | Every protein has an explicit record in `pipeline.json` and `pipeline.xlsx`. |
| **Alignment-aware** | Sequence logic is not an afterthought; UniProt alignment is central. |
| **Modular** | Each preparation stage lives in its own module. |
| **Debuggable** | Paths, outputs, statuses, and logs are visible and inspectable. |
| **Extensible** | Metal parametrization and downstream MD preparation can grow naturally from the current architecture. |

This makes FRUTON useful not only as a pipeline, but as a **scientific working framework**.

---

## Pipeline at a glance

### Structural flow

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
     ├─────────────── no internal gaps ───────────────┐
     │                                                │
     ▼                                                │
[ FILLER: MODELLER / ALPHAFOLD ]                      │
     │                                                │
     ▼                                                │
[ PROTONATION ]                                       │
     │                                                │
     ▼                                                │
[ AMBER RENAMING ]                                    │
     │                                                │
     ▼                                                │
[ FINAL NUMBERING ]                                   │
     │                                                │
     ▼                                                ▼
[ FINAL PROTEIN ] ----------------------------> [ METAL PREPARATION ]
```

### Conceptual layers

| Layer | Purpose |
|---|---|
| **Sequence layer** | FASTA generation, UniProt matching, mapping TSVs |
| **Structure layer** | insertion cleanup, component split, gap detection |
| **Reconstruction layer** | MODELLER / AlphaFold-based filling |
| **Chemistry layer** | protonation, AMBER-compatible residue naming |
| **Finalization layer** | final output numbering and normalization |
| **Metal branch** | preparation for later MCPB / Gaussian workflows |

---

## What FRUTON currently does

### Implemented steps

| Step | Module / Logic | Current role |
|---|---|---|
| 1 | `pdb_sync` | synchronize protein set from CSV |
| 2 | `fasta_files` | generate PDB- and UniProt-derived FASTA files |
| 3 | `sequence_alignment` | align chain-specific PDB sequences to UniProt |
| 4 | `insertion_codes` | remove insertion-code complications |
| 5 | `pdb_components` | split protein / water / ligand / metal |
| 6 | `gaps` | detect and summarize structural gaps |
| 7 | `filler` | use MODELLER or AlphaFold fallback |
| 8 | `protonation` | protonate structural protein |
| 9 | `amber_renaming` | assign AMBER-compatible residue names |
| 10 | `finalize_protein` | create final numbered output |
| 11 | `metall_params` | prepare metal-containing systems for later parametrization |

---

## Repository layout

```text
stack_protein_prep/
├── pixi.toml
├── scripts/
│   └── run_pipeline.py
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
│       ├── finalize_protein.py
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

## Per-protein data layout

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
└── metall_params/
    ├── tmp_param.pdb
    ├── metal_only.pdb
    ├── chimera_contacts.py
    ├── contacts.data
    └── chimera_run.log
```

---

## Key outputs and their meaning

| File | Meaning |
|---|---|
| `<PDB_ID>_delins.pdb` | insertion-cleaned PDB |
| `PDB-<PDB_ID>-SEQRES.fasta` | sequence extracted from SEQRES |
| `PDB-<PDB_ID>-ATOM.fasta` | sequence extracted from observed ATOM records |
| `UniProt_<UNIPROT_ID>.fasta` | UniProt reference sequence |
| `*_vs_UniProt.aln.fasta` | aligned chain-specific comparison |
| `*_vs_UniProt.mapping.tsv` | alignment mapping table |
| `<PDB_ID>_protein.pdb` | protein component only |
| `<PDB_ID>_proteinH.pdb` | protonated protein |
| `<PDB_ID>_protein_as_Amber.pdb` | AMBER-renamed protein |
| `<PDB_ID>_protein_final.pdb` | final normalized protein |
| `<PDB_ID>_metal.pdb` | metal component |
| `<PDB_ID>_finalize_numbering.tsv` | explicit final numbering table |
| `metall_params/tmp_param.pdb` | combined metal-analysis input |
| `metall_params/contacts.data` | Chimera clash/contact output |

---

## Why each stage matters

### FASTA generation and UniProt alignment

FRUTON does not treat sequence logic as optional decoration.  
Sequence alignment is one of the structural backbone layers of the whole project.

It is used to:

- compare observed structure to reference biology
- map chains against UniProt
- detect terminal truncation or missing regions
- support gap classification
- support later numbering logic

Without that layer, later “fixes” would be much less reliable.

---

### Insertion cleanup

Insertion codes are one of those structural annoyances that seem small until they break everything:

- residue mapping
- indexing
- chain sequence extraction
- modelling assumptions
- alignment interpretation

Cleaning them early reduces downstream ambiguity.

---

### Component split

Once a structure is split into:

- protein
- water
- ligand
- metal

the rest of the pipeline becomes far more manageable.

This step makes later logic cleaner because downstream modules no longer need to guess whether a residue is protein, solvent, ligand, or ion.

---

### Gap detection and filler

This is where FRUTON stops being a simple cleanup tool and becomes a reconstruction workflow.

Internal gaps are detected, classified, and used to decide whether to:

- skip filling
- run MODELLER
- use AlphaFold fallback

This makes the workflow explicit and inspectable instead of manual and ad hoc.

---

### Protonation and AMBER renaming

FRUTON separates:

1. **adding hydrogens / assigning protonation**
2. **renaming residues for AMBER-style semantics**

This is important because those are related, but not identical, tasks.

The result is a chemically more meaningful and downstream-ready structure.

---

### Finalize protein

The final structural output is not just “the latest PDB”.  
It is the result of a deliberate normalization step.

This includes:

- selecting the final structural source
- applying a dedicated numbering policy
- producing a canonical final protein file

That file is the current preferred handoff to later stages.

---

### Metal preparation

The metal branch is where FRUTON begins to connect structural preparation to later **parametrization workflows**.

Current goal:

- combine final protein + optional water + optional ligand + metal
- analyze local contacts around the metal with Chimera
- generate starting input for later MCPB / Gaussian preparation

This is still evolving, but the branch architecture is already in place.

---

## Pipeline state

FRUTON tracks pipeline state in:

```text
data/proteins/pipeline.json
data/proteins/pipeline.xlsx
```

### JSON
Machine-readable, complete, debug-friendly.

### XLSX
Human-readable overview for quick inspection.

### Why this matters

The flat state model makes it easy to answer:

- Which proteins still have gaps?
- Which proteins already have final outputs?
- Which proteins contain metals?
- Which ones failed protonation?
- Which ones are ready for metal preparation?

This is one of the most useful parts of the project.

---

## Example state fields

### Identity and directories

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

### Chemistry / model outputs

- `filler_directory`
- `filler_model_path`
- `protonation.output_path`
- `amber_renaming.output_path`
- `numbering_restore.output_path`

### Status fields

- `pdb_sync_done`
- `fasta_files_done`
- `sequence_alignment_done`
- `insertion_codes_done`
- `filler_status`
- `protonation.status`
- `amber_renaming.status`
- `numbering_restore.status`
- `metall_params.status`

---

## Current strengths

FRUTON already has a strong foundation for:

- alignment-aware structure preparation
- explicit intermediate files
- reproducible state tracking
- modular stepwise orchestration
- chain-specific sequence logic
- chemistry-aware structural normalization
- branching toward metal-site workflows

---

## Current limitations

FRUTON is still actively evolving, especially in the later stages.

### Areas still being refined

- finalize-numbering logic
- MODELLER / AlphaFold edge cases
- metal-site selection logic in Chimera
- atom serial renumbering in merged PDBs
- automatic contact parsing for metal workflows

### Known practical issue

Merged files such as `tmp_param.pdb` may contain duplicate atom serial numbers because concatenation is currently used instead of full renumbering.

That is usually survivable, but it should eventually be cleaned up.

---

## Recommended next steps

1. **formalize finalize numbering generation**
2. **improve Chimera metal-selection logic**
3. **renumber merged atom serials**
4. **parse `contacts.data` into structured machine-readable output**
5. **build MCPB input automatically**
6. **extend the metal branch toward Gaussian preparation**
7. **improve rerun granularity for later pipeline stages**

---

## Running FRUTON

Typical full pipeline run:

```bash
pixi run python scripts/run_pipeline.py
```

---

## Main source files

### Orchestration

```text
scripts/run_pipeline.py
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

### Chemistry and finalization

```text
src/stack_protein_preparation/protonation.py
src/stack_protein_preparation/amber_renaming.py
src/stack_protein_preparation/finalize_protein.py
src/stack_protein_preparation/metall_params.py
```

---

## In one sentence

**FRUTON is a filesystem-explicit, state-driven protein-preparation framework that connects sequence-aware cleanup, structure reconstruction, chemistry-aware normalization, and metal-site preparation into one coherent workflow.**

---
