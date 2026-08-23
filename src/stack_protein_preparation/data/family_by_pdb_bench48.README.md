# Bench48 family-assignment template

`family_by_pdb_bench48.template.json` holds all 48 PDB IDs from
`scripts/CESGA_SLURM/fruton_bench48_full_results.json` mapped to
`"__unassigned__"`.  Populate it with real family labels (from the
canonical vocabulary in `_family_stratification.FAMILY_LABELS`) before
running `scripts/plot_family_stratification.py` on this bench.

## Why a template instead of pre-filled labels?

Assigning families to modern PDB entries requires querying UniProt
per PDB and reading the enzyme-classification (EC) or the primary
paper.  Inventing labels without that lookup is worse than honest
`__unassigned__` because a wrong family assignment turns into a
wrong reviewer figure.  The template is the vehicle for the correct
assignments — it does not pretend to already know them.

## Populating

Two paths:

### A) UniProt SPARQL / REST lookup (recommended)

For each PDB in the template:

1. Query `https://data.rcsb.org/graphql` with
   `entries(entry_ids: "<PDB>") { rcsb_entry_container_identifiers { uniprot_ids } }`
   to get the UniProt accession(s).
2. Query `https://rest.uniprot.org/uniprotkb/<ACC>.json` for
   `proteinDescription` + `comments` + `keywords` fields.
3. Map to a canonical `FAMILY_LABELS` entry using keyword rules:
   - keyword `"Kinase"` → `kinase`
   - keyword `"G-protein coupled receptor"` → `gpcr`
   - keyword `"Serine protease"` / `"Aspartyl protease"` / `"Metalloprotease"` → `protease`
   - keyword `"Metal-binding"` + EC 3.x/4.x on a small metal-active domain → `metalloenzyme`
   - keyword `"Nuclear receptor"` → `nuclear_receptor`
   - keyword `"Protein phosphatase"` → `phosphatase`
   - EC 2.x (transferase) without kinase → `transferase`
   - EC 3.x (hydrolase) without protease → `hydrolase`
   - EC 1.x → `oxidoreductase`
   - `cofactor_dependent` when the paper describes NAD/FAD/HEM/SAM/PLP
   - `multi_domain` when the UniProt entry reports > 3 domains
     spanning distinct folds
   - otherwise `other`

### B) Manual per-PDB review

Open each PDB entry on `https://www.rcsb.org/structure/<PDB>` and
pick the most-specific canonical label.  Take ~30 seconds per PDB.

## Then

Copy the populated file to `family_by_pdb_bench48.json` (drop
`.template`) and run:

```bash
python scripts/plot_family_stratification.py \
    --bench-json scripts/CESGA_SLURM/fruton_bench48_full_results.json \
    --family-map src/stack_protein_preparation/data/family_by_pdb_bench48.json \
    --outdir bench_output/nature_r1_family_bench48/
```

## Existing well-known seed vs bench48

`family_by_pdb_seed.json` (24 entries) is the canonical
didactic-example map used by tests and demos.  It intentionally does
not overlap `bench48` because bench48 was assembled to stress-test
FRUTON's gap-filling on modern crystal structures with unresolved
loops — different selection criterion.
