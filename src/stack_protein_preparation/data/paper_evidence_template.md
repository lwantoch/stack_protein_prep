# Paper evidence — <PDB_ID> <SHORT NAME>

Use this template to record literature-derived active-site residue
information that FRUTON's `_paper_override_suggest.py` can convert
into a `.suggested.json` override manifest.  The extractor scans for
`### <RESNAME><resnum>` headings followed by `> quoted phrase` lines,
so the exact structure matters.

## Provenance

- **PDB entry**: <PDB_ID>
- **Primary paper**: <first-author> et al. (<year>) *<journal>*
  DOI: <10.xxxx/…>
- **Extracted by**: <your name / date>

## Active-site protonation evidence

Emit one `###` heading per residue.  Below each heading, quote one or
more sentences from the paper using `>`.  The extractor searches for
these keywords in the quote (case-insensitive):

- **ASP / GLU** → `ASH` / `GLH` on: `catalytic acid`,
  `proton donor`, `protonated asp/glu`, `neutral asp/glu`
- **HIS** → `HIP` on: `catalytic base`, `proton acceptor`,
  `protonated his`, `positively charged his`, `hip`
- **HIS** → `HID` on: `hid`, `nε2 protonated`, `ne2 protonated`
- **HIS** → `HIE` on: `hie`, `nδ1 protonated`, `nd1 protonated`
- **CYS** → `CYM` on: `nucleophile`, `thiolate`, `catalytic cys`,
  `deprotonated cys`, `cys thiolate`

Vague quotes (no matching keyword) go to `needs_review` and produce no
override suggestion; the pipeline owner should edit them manually.

### HIS<N>

> Direct quote from the paper describing the mechanistic role of this
> histidine.  Use terms from the keyword list above so the extractor
> emits a concrete override.

### ASP<N>

> Direct quote describing whether this aspartate is a catalytic acid
> (should be protonated ASH) or a general base (usually stays ASP).

### CYS<N>

> Direct quote about redox / nucleophile / thiolate character.

## Non-protonation notes (advisory, not extracted)

Free-text notes go here.  They are ignored by the extractor but useful
for reviewers and downstream analysis.  Metal ligands, cofactor
identity, disulfide connectivity, and any oxidation-state hints all
belong in this section.

## Suggested workflow

1. Fill this file in for the protein under study, save as
   `paper_evidence.md` next to the crystal PDB in the bench tree.
2. Run:
   ```bash
   python -c "
   from stack_protein_preparation._paper_override_suggest import write_suggested_overrides
   write_suggested_overrides('paper_evidence.md', 'active_site_overrides.suggested.json')
   "
   ```
3. Review the `.suggested.json` file, promote entries to
   `active_site_overrides.json` when correct.  FRUTON's protonation
   step reads only the real (`.json`) file, never `.suggested.json`.
