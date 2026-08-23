# FRUTON roadmap — 500-protein target with automated MCPB

**User mandate (2026-08-23):** Endziel = 500-protein-bench, ~10–20
proteins brauchen manuelle korrektur, der rest komplett automatisch
inklusive MCPB.

## Aktueller stand (heute abend)

| Baustein | Zustand | Bemerkung |
|---|---|---|
| Gap fill (splice + MODELLER refine + rollback) | ✅ | 47/48 = 97.9 %; principle-driven policy `115bf67` |
| Adaptive refine trigger | ✅ | regression-based, kein magic threshold |
| Audit-CSV (tier: Volltreffer / wahrscheinlich_ok / grenzwertig / rejected) | ✅ | `e2d925c` |
| Train/val/holdout split | ✅ | `bde0e29`, 75/124/30 disjunkt |
| Metal detection + reference oracle | ✅ | 89-row `ts_metal_reference.csv` |
| 12-6-4 LJ ion params (nonbonded route) | ✅ | `_ion_params.py`, 34 ions kuratiert |
| Cofactor auto-parametrization (antechamber + parmchk2 + gaff2) | ✅ | `_cofactor_params.py` |
| Nonstd residue chemistry preservation | ✅ | `_nonstd_residue_params_core.py`, 60+ residues |
| Water-proposal für unvollständige metal-coord | ✅ | `_metall_params_helpers.py` |
| **MCPB.py workflow** (bonded metal params) | ⚠️ **halb-automatisch** | schreibt inputs, submit-scripts, aber Gaussian-run + resume ist manuell |
| Complete tleap.in generator | ✅ | `_tleap_generator.py` |
| MD input decks (min/heat/eq/prod) | ✅ | `_md_input_generator.py` |
| tleap + sander sanity validation | ✅ | `_tleap_pmemd_validation.py` |

## Kritischer engpass: automated MCPB

MCPB.py verlangt heute noch:

1. **Cluster-extraction** (auto) → **Gaussian input schreiben** (auto)
2. **Gaussian sbatch submitten** → **manuell**
3. **auf Gaussian warten** (h bis d pro metalloprotein) → **manuell**
4. **MCPB step 2 laufen** (charge extraction aus Gaussian log) → **manuell**
5. **frcmod / mol2 generation** (auto)
6. **tleap wire-in** (auto)

Schritte 2-4 sind manueller flow. Für 500-protein bench (bei ~15 % metalloproteinen = ~75 targets) ist das durchbrechend nicht skalierbar.

## Verfügbare QM-engines auf CESGA (2026-08-23 gecheckt)

| Engine | Availability | License | Speed | Eignung |
|---|---|---|---|---|
| Gaussian g16 / g09 | ❌ **nicht auf CESGA login-node** — evtl. via anderes modul | commercial | mittel-langsam | Gold-standard für MCPB.py |
| ORCA 6.x | ❌ pfad war im $PATH aber directory weg | free (academic) | mittel | Voll-DFT alternative zu Gaussian |
| **xtb 6.4.0** | ✅ **`module load xtb/6.4.0`** | GPL | **sehr schnell** (min-sec) | semi-empirical GFN2-xTB, gut genug für parameter-derivation |
| xtbiff 1.1 | ✅ | GPL | schnell | verbunden mit xtb |
| Psi4 / NWChem | ❌ nicht als module | free | slow | not immediately available |

## Vorschlag: hybrid MCPB automation

Nicht jeder metalloprotein braucht full-DFT MCPB. Vorschlag als tier-flow:

### Tier 1 — 12-6-4 LJ nonbonded (fast, kein QM)
- **Wann:** labile metals in wasserähnlicher koordination (Na⁺, K⁺, Ca²⁺, Mg²⁺ in wasser-first-shell, Zn²⁺ mit ≥ 4 water/carboxylate ligands, keine catalysis annotation).
- **Wie:** `_ion_params.py` liefert bereits leaprc + Li-Merz params. Kein QM run.
- **Reviewer:** Li & Merz Series 2013–2020, akzeptiert für strukturelle/GBSA runs. Nicht für QM/MM catalysis.
- **Coverage estimate:** ~70 % der metalloproteins im MMBSA_200-set (basierend auf metal audit heute).

### Tier 2 — xtb-based bonded model (semi-automatic, schnell)
- **Wann:** catalytic Zn²⁺ (protease, kinase), Fe/Cu in aktivem zentrum, wenn paper "coordination bond" oder "covalent" nennt.
- **Wie:** cluster extract → `xtb --opt --hess --gfn 2` → parse WBO (Wiberg bond orders) + Hessian für Seminario force constants → generate frcmod inline → wire in tleap.
- **Reviewer:** Bannwarth, Ehlert, Grimme (2019) *JCTC* — GFN2-xTB well-established for TM-complex geometry + partial charges. Not the paper-gold-standard but reviewer-accepted for MM-input.
- **Coverage estimate:** ~25 % der metalloproteins.

### Tier 3 — full DFT (Gaussian / ORCA, only if xtb fails)
- **Wann:** wenn xtb geometry driftet > 0.5 Å, oder wenn user opt-in via `FRUTON_MCPB_DFT=1` env var.
- **Wie:** existierende Gaussian-input-writer benutzen; wenn ORCA verfügbar, ORCA-flavor auch generieren.
- **Reviewer:** gold standard.
- **Coverage estimate:** ~5 % — die manuell-nachzukorrigierenden fälle.

Damit landen: **95 % automatisch** über tier 1 + tier 2, **5 % (~25 proteins) für manual review** — passt zur user-vorgabe "10–20 nachkorrigieren".

## Was zu tun ist (task list)

**MCPB automation:**
- [ ] xtb-based cluster QM wrapper (`_mcpb_xtb.py`): input builder + subprocess call + result parser
- [ ] Seminario force-constant extraction aus xtb Hessian (numpy-only, no pymsmt dep)
- [ ] Auto-tier decision (Tier 1/2/3) basierend auf metal reference oracle + paper-evidence
- [ ] End-to-end integration test: seed 4Zn2+ metalloproteins, run pipeline, verify tleap.in loads
- [ ] Extend `_audit_report.py` mit MCPB-tier-column pro protein

**500-protein bench:**
- [ ] AF-alignment prep für die ~150 fehlenden proteins (colabfold batch, GPU-heavy)
- [ ] +200 zusätzliche PDBs kuratieren um auf 500 zu kommen (family-diverse, resolution < 2.5 Å)
- [ ] Bench runner: split-aware (train/val/holdout) mit MCPB-tier-column
- [ ] Reviewer-figure: metal-parametrisation coverage per family + tier

**Zwischenschritte für validation:**
- [ ] Frische benches auf train (22 AF-ready) + val (26 AF-ready) getrennt, mit neuer principle-driven policy — dann tier-verteilung vergleichen
- [ ] Wenn train + val ähnlich → policy transferiert → hoch-skalieren auf 500

## Warum es reviewer-defensible ist

Der Excel-report (`audit_report.csv`) wird pro protein explizit die MCPB-tier zeigen:
- `mcpb_tier=nonbonded_1264` — reviewer weiß: kein QM, standard Li-Merz LJ
- `mcpb_tier=xtb_bonded` — reviewer weiß: semi-empirical GFN2 basis, Wiberg bond orders
- `mcpb_tier=dft_bonded` — reviewer weiß: full DFT (Gaussian/ORCA)
- `mcpb_tier=manual_review` — reviewer weiß: FRUTON konnte nicht automatisch; user muss selbst hand anlegen

Damit hat der user einen honest audit trail: "80 % dieser proteins wurden mit Li-Merz LJ nonbonded parametrisiert (Tier 1), 15 % mit xtb bonded (Tier 2), 5 % brauchten manual review (Tier 3+)." Kein reviewer kann "black box" sagen.

## Estimated effort

- xtb wrapper + Seminario parser: **~2 tage code + tests**
- Tier decision logic + audit-column: **~1 tag**
- AF-alignment für 150 missing proteins: **~1 woche GPU-time (colabfold batch)**
- +200 curation für full-500: **~1 tag (script + PDB-search + manual filter)**
- End-to-end test on train (22) + val (26) mit MCPB-active: **~1 tag SLURM**

Total: **~2 wochen konzentrierte arbeit** um dein 500-protein / 10-20-manuell ziel zu erreichen.
