#!/usr/bin/env python3
"""Generate the 3 juplaunch-ready review notebooks.

Run once (or after editing this file) to (re)build:
    notebooks/1_review_newbench27_fills.ipynb
    notebooks/2_review_mmbsa200_fills.ipynb
    notebooks/3_review_protonation_side_by_side.ipynb

Kept as a script (not hand-edited notebooks) so cell content is
version-controlled + trivially regeneratable.
"""
from __future__ import annotations

import json
from pathlib import Path

import nbformat as nbf


HERE = Path(__file__).resolve().parent

# ---------------------------------------------------------------------------
# Shared setup cell used by all 3 notebooks
# ---------------------------------------------------------------------------

SETUP_MD = """\
# {title}

**Runs on CESGA via `juplaunch`**  (memory `[[reference_juplaunch]]`).

## Connection setup

On your **laptop** before launching this notebook:

```bash
# 1. Start PyMOL with the RPC server on the laptop
export PYMOL_PATH=/usr/lib/python3/dist-packages/pymol   # Debian-fix per feedback_pymol_remote_debian
pymol -R
# → PyMOL RPC listening on :9123

# 2. Set up the reverse tunnel to CESGA
ssh -R 9123:localhost:9123 ft3.cesga.es
```

On CESGA:

```bash
~/bin/juplaunch     # prints the JupyterLab URL to paste into your laptop browser
```

Then open this notebook in the browser and run the cells top-down.
The `PymolSession(hostname="localhost", port=9123)` call reaches through
the reverse tunnel to the laptop's PyMOL GUI.
"""

CONNECT_CELL = """\
# Connect to laptop PyMOL via reverse tunnel (localhost:9123)
from pymol_remote.client import PymolSession
import sys, pathlib

pm = PymolSession(hostname="localhost", port=9123, timeout=10.0)
print("Connected. Objects currently loaded:", pm.get_names())
"""

DISCOVER_NEWBENCH_CELL = """\
# Discover delivered newbench_27 proteins
import pathlib, json

DELIVERED_ROOT = pathlib.Path("/mnt/netapp1/Store_othcxlwa/newbench_27/delivered")
BENCH_ROOT     = pathlib.Path("/mnt/netapp1/Store_othcxlwa/newbench_27")

pdb_ids = sorted(p.name for p in DELIVERED_ROOT.iterdir() if p.is_dir())
print(f"{len(pdb_ids)} delivered newbench_27 proteins:")
print(" ".join(pdb_ids))
"""

DISCOVER_MMBSA200_CELL = """\
# Discover MMBSA_200 proteins with available filled models
import pathlib, glob

CRYSTAL_ROOT = pathlib.Path("/mnt/netapp1/Store_othcxlwa/FRUTON-NEW")
FILLED_GLOBS = [
    pathlib.Path("/mnt/lustre/scratch/nlsas/home/otras/hcx/lwa/MMBSA_200/bench_20260823_1620/artefacts"),
    pathlib.Path("/mnt/lustre/scratch/nlsas/home/otras/hcx/lwa/MMBSA_200/bench_20260823_1513/artefacts"),
]

crystal_pdbs = sorted(p.name for p in CRYSTAL_ROOT.iterdir() if p.is_dir())
filled_pdbs = set()
for fg in FILLED_GLOBS:
    if fg.is_dir():
        for f in fg.glob("*_final.pdb"):
            filled_pdbs.add(f.stem.replace("_final", ""))

pairs = sorted(p for p in crystal_pdbs if p in filled_pdbs)
print(f"{len(crystal_pdbs)} crystals total; {len(filled_pdbs)} have filled models; "
      f"{len(pairs)} pairs available for review.")
print("Pairs:", " ".join(pairs))
"""

# ---------------------------------------------------------------------------
# Notebook 1 — newbench_27 fills
# ---------------------------------------------------------------------------

NB1_LOADER_CELL = """\
# Loader: for a given PDB id, show crystal vs graft in PyMOL,
# highlighting the gap-region residues that FRUTON filled.
import pathlib

def load_newbench_fill(pdb_id: str) -> dict:
    entry = DELIVERED_ROOT / pdb_id
    crystal = BENCH_ROOT / pdb_id / f"{pdb_id}.pdb"
    graft = entry / "graft" / f"{pdb_id}_graft.pdb"
    receptor = entry / "receptor.pdb"
    manifest_path = entry / "graft" / f"{pdb_id}_graft_manifest.json"
    manifest = json.loads(manifest_path.read_text()) if manifest_path.is_file() else {}
    return {
        "pdb_id": pdb_id,
        "crystal": crystal if crystal.is_file() else None,
        "graft": graft if graft.is_file() else None,
        "receptor": receptor if receptor.is_file() else None,
        "manifest": manifest,
    }


def pymol_show_fill(pdb_id: str, gray_crystal: bool = True) -> dict:
    r = load_newbench_fill(pdb_id)
    pm.do("reinitialize")
    pm.do("bg_color white")
    if r["crystal"]:
        pm.do(f"load {r['crystal']}, crystal")
        if gray_crystal:
            pm.do("color grey70, crystal")
        pm.do("show cartoon, crystal")
    if r["receptor"]:
        pm.do(f"load {r['receptor']}, filled")
        pm.do("color skyblue, filled")
        pm.do("show cartoon, filled")
    # Highlight graft-region residues in the filled model, if manifest lists them.
    graft_ranges = r["manifest"].get("graft_ranges") or r["manifest"].get("ranges") or []
    for rng in graft_ranges:
        chain = rng.get("chain", "A")
        lo = rng.get("first_resnum") or rng.get("start")
        hi = rng.get("last_resnum") or rng.get("end")
        if lo is None or hi is None:
            continue
        sel = f"filled and chain {chain} and resi {lo}-{hi}"
        pm.do(f"color hotpink, {sel}")
        pm.do(f"show sticks, {sel} and (name CA+CB+N+C+O)")
    if r["crystal"] and r["receptor"]:
        pm.do("align filled and name CA, crystal and name CA")
    pm.do("orient")
    pm.do("zoom")
    print(f"{pdb_id}: crystal={bool(r['crystal'])}, filled={bool(r['receptor'])}, "
          f"n_graft_ranges={len(graft_ranges)}")
    return r
"""

NB1_WIDGET_CELL = """\
# Interactive selector — pick a PDB and PyMOL updates on the laptop
import ipywidgets as W
from IPython.display import display

dropdown = W.Dropdown(options=pdb_ids, description="PDB:", layout={"width": "260px"})
gray_toggle = W.Checkbox(value=True, description="Grey crystal")
btn = W.Button(description="Show in PyMOL", button_style="primary")
out = W.Output()

def _on_click(_):
    with out:
        out.clear_output()
        pymol_show_fill(dropdown.value, gray_crystal=gray_toggle.value)

btn.on_click(_on_click)
display(W.HBox([dropdown, gray_toggle, btn]), out)
"""

NB1_BATCH_CELL = """\
# Optional: walk every delivered newbench_27 protein one at a time,
# saving a session file per PDB so you can re-open later without re-running.
import pathlib
SESSION_DIR = pathlib.Path("~/pymol_sessions_newbench27").expanduser()
SESSION_DIR.mkdir(parents=True, exist_ok=True)

for pid in pdb_ids:
    pymol_show_fill(pid, gray_crystal=True)
    pm.do(f"save {SESSION_DIR / (pid + '.pse')}")
    print("  session ->", SESSION_DIR / (pid + ".pse"))
"""


# --- KI / Masse / Actives overview cells ---------------------------------

NB1_OVERVIEW_PARSE_CELL = """\
# Parse every delivered protein's manifest + graft_manifest and build a
# per-PDB summary dataframe.
import json, pathlib, re
import pandas as pd

def _residue_count(pdb_path: pathlib.Path) -> int:
    seen = set()
    if not pdb_path.is_file():
        return 0
    for line in pdb_path.read_text(errors='replace').splitlines():
        if line.startswith('ATOM') and len(line) > 26:
            key = (line[21], line[22:27].strip())  # chain + resnum+icode
            seen.add(key)
    return len(seen)

# Approximate mean AA molecular weight (Da) for quick MW estimate
_AVG_AA_MW = 110.0

def summarise(pdb_id: str) -> dict:
    entry = DELIVERED_ROOT / pdb_id
    manifest = {}
    if (entry / 'MANIFEST.json').is_file():
        manifest = json.loads((entry / 'MANIFEST.json').read_text())
    graft = {}
    gp = entry / 'graft' / f'{pdb_id}_graft_manifest.json'
    if gp.is_file():
        graft = json.loads(gp.read_text())
    receptor = entry / 'receptor.pdb'
    n_res = _residue_count(receptor) if receptor.is_file() else 0
    mw_kda = round(n_res * _AVG_AA_MW / 1000.0, 1)
    variant = manifest.get('variant', '—')
    method = graft.get('method', '') or ''
    # AI split heuristic — primary key = MANIFEST.json 'variant' field:
    #   'single'              → no fill needed (crystal is complete)
    #   'large_gap_graft'     → BL-Pose Kabsch graft (no AF)
    #   'large_gap_complete'  → AF-splice + MODELLER LoopModel
    #   'best_complete'       → AF-splice (best conformer, no explicit MODELLER)
    if variant == 'single':
        ai_source = 'no fill (crystal complete)'
    elif variant == 'large_gap_graft' or 'Kabsch' in method or 'Superimposer' in method:
        ai_source = 'BL-Pose (Kabsch graft, no AF)'
    elif variant == 'large_gap_complete' or 'MODELLER' in method:
        ai_source = 'AF-splice + MODELLER'
    elif variant == 'best_complete':
        ai_source = 'AF-splice (best conformer)'
    else:
        ai_source = 'unknown/other'
    n_breaks = manifest.get('backbone_continuity', {}).get('n_breaks', 0)
    n_warnings = len(manifest.get('warnings', []))
    return {
        'pdb_id': pdb_id,
        'n_residues': n_res,
        'MW_kDa_estimate': mw_kda,
        'variant': variant,
        'ai_source': ai_source,
        'n_gaps_grafted': graft.get('n_gaps_grafted', 0),
        'n_gaps_skipped': len(graft.get('skipped', []) or []),
        'n_backbone_breaks': n_breaks,
        'n_warnings': n_warnings,
    }

pdb_ids_all = sorted(p.name for p in DELIVERED_ROOT.iterdir() if p.is_dir())
summary_df = pd.DataFrame([summarise(p) for p in pdb_ids_all])
summary_df
"""

NB1_AI_SPLIT_CELL = """\
# KI-Aufteilung: how many proteins used AF-splice vs BL-Pose fallback?
counts = summary_df['ai_source'].value_counts()
print(counts.to_string())

import matplotlib.pyplot as plt
fig, ax = plt.subplots(figsize=(7, 3.5))
counts.plot.barh(ax=ax, color=['#2b6cb0', '#c05621', '#a0aec0'][:len(counts)])
ax.set_xlabel('n proteins')
ax.set_title('newbench_27 — AI source split (AF-splice vs BL-Pose fallback)')
for i, (label, v) in enumerate(counts.items()):
    ax.text(v + 0.2, i, f'  {v}', va='center')
plt.tight_layout()
plt.show()
"""

NB1_MASS_HISTOGRAM_CELL = """\
# Massenaufteilung: molecular-weight distribution across delivered receptors
import matplotlib.pyplot as plt
mw = summary_df['MW_kDa_estimate'].dropna()

fig, ax = plt.subplots(figsize=(9, 4))
ax.hist(mw, bins=15, color='#2b6cb0', edgecolor='none')
ax.set_xlabel('Molecular weight (kDa, estimate = 110 Da × n_residues)')
ax.set_ylabel('n proteins')
ax.set_title(f'newbench_27 — mass distribution '
             f'(n={len(mw)}, median={mw.median():.1f} kDa, '
             f'range {mw.min():.1f}–{mw.max():.1f} kDa)')
ax.axvline(mw.median(), color='grey', linestyle=':', linewidth=0.8,
           label=f'median {mw.median():.1f} kDa')
ax.legend()
plt.tight_layout()
plt.show()

# Also show the top-5 largest + smallest
print('Largest by MW:')
print(summary_df.nlargest(5, 'MW_kDa_estimate')[['pdb_id','MW_kDa_estimate','n_residues']].to_string(index=False))
print()
print('Smallest by MW:')
print(summary_df.nsmallest(5, 'MW_kDa_estimate')[['pdb_id','MW_kDa_estimate','n_residues']].to_string(index=False))
"""

NB1_ACTIVES_INACTIVES_CELL = """\
# Actives / inactives per target — scan for the standard file names.
# Newbench_27 does not (yet) ship ligand sets under delivered/;
# this cell looks in a few known locations and reports what's found.
import pathlib

ACTIVES_CANDIDATE_ROOTS = [
    pathlib.Path('/mnt/netapp1/Store_othcxlwa/DEKOIS-EXP'),
    pathlib.Path('/mnt/netapp1/Store_othcxlwa/newbench_27/actives'),
    pathlib.Path('/mnt/netapp1/Store_othcxlwa/newbench_27/ligands'),
]

def _count_lines(p: pathlib.Path) -> int:
    try:
        return sum(1 for _ in p.open())
    except Exception:
        return 0

def _find_ligand_files(pdb_id: str) -> dict:
    hits = {'actives_smi': None, 'inactives_smi': None,
            'n_actives': 0, 'n_inactives': 0}
    for root in ACTIVES_CANDIDATE_ROOTS:
        if not root.is_dir():
            continue
        # Try lowercase, uppercase, exact-case
        for variant in (pdb_id, pdb_id.lower(), pdb_id.upper()):
            for pattern in ('actives.smi', f'{variant}_actives.smi', 'actives.ism', f'{variant}.ism'):
                p = root / variant / pattern
                if p.is_file():
                    hits['actives_smi'] = str(p); hits['n_actives'] = _count_lines(p); break
            for pattern in ('inactives.smi', f'{variant}_inactives.smi', 'decoys.smi', f'{variant}_decoys.smi'):
                p = root / variant / pattern
                if p.is_file():
                    hits['inactives_smi'] = str(p); hits['n_inactives'] = _count_lines(p); break
    return hits

ligand_rows = []
for pid in pdb_ids_all:
    h = _find_ligand_files(pid)
    h['pdb_id'] = pid
    ligand_rows.append(h)
ligand_df = pd.DataFrame(ligand_rows)[['pdb_id', 'n_actives', 'n_inactives', 'actives_smi', 'inactives_smi']]

n_with = sum(1 for r in ligand_rows if r['n_actives'] > 0 or r['n_inactives'] > 0)
print(f'{n_with} / {len(ligand_rows)} proteins have ligand set data on disk.')
if n_with == 0:
    print(
        'No actives/inactives sets found in the standard locations.\\n'
        'If your ligands live elsewhere, extend ACTIVES_CANDIDATE_ROOTS in this cell.'
    )
ligand_df
"""

NB1_MERGED_TABLE_CELL = """\
# Merged summary table (KI + Masse + actives/inactives) per protein
merged = summary_df.merge(ligand_df, on='pdb_id', how='left')
merged = merged[['pdb_id', 'ai_source', 'MW_kDa_estimate', 'n_residues',
                 'n_gaps_grafted', 'n_gaps_skipped', 'n_backbone_breaks',
                 'n_warnings', 'n_actives', 'n_inactives', 'variant']]
merged
"""

# ---------------------------------------------------------------------------
# Notebook 2 — MMBSA_200 fills
# ---------------------------------------------------------------------------

NB2_LOADER_CELL = """\
# Loader for MMBSA_200: crystal from FRUTON-NEW/<PDB>/<PDB>.pdb,
# filled from SLURM bench_20260823_1620 (iter-4) or _1513 (iter-3) if newer.
import pathlib, json

def _pick_filled(pdb_id: str) -> pathlib.Path | None:
    for fg in FILLED_GLOBS:
        p = fg / f"{pdb_id}_final.pdb"
        if p.is_file():
            return p
    return None


def load_mmbsa200_fill(pdb_id: str) -> dict:
    return {
        "pdb_id": pdb_id,
        "crystal": (CRYSTAL_ROOT / pdb_id / f"{pdb_id}.pdb")
                   if (CRYSTAL_ROOT / pdb_id / f"{pdb_id}.pdb").is_file() else None,
        "filled": _pick_filled(pdb_id),
    }


def pymol_show_mmbsa_fill(pdb_id: str, gray_crystal: bool = True) -> dict:
    r = load_mmbsa200_fill(pdb_id)
    pm.do("reinitialize")
    pm.do("bg_color white")
    if r["crystal"]:
        pm.do(f"load {r['crystal']}, crystal")
        pm.do(("color grey70, crystal" if gray_crystal else "color yellow, crystal"))
        pm.do("show cartoon, crystal")
    if r["filled"]:
        pm.do(f"load {r['filled']}, filled")
        pm.do("color skyblue, filled")
        pm.do("show cartoon, filled")
        # Highlight residues present in filled but absent in crystal
        pm.do("select gap_fill, filled and not (byres filled within 0.5 of crystal)")
        pm.do("color hotpink, gap_fill")
        pm.do("show sticks, gap_fill and (name CA+CB+N+C+O)")
    if r["crystal"] and r["filled"]:
        pm.do("align filled and name CA, crystal and name CA")
    pm.do("orient")
    pm.do("zoom")
    print(f"{pdb_id}: crystal={bool(r['crystal'])}, filled={bool(r['filled'])}")
    return r
"""

NB2_WIDGET_CELL = """\
import ipywidgets as W
from IPython.display import display

if not pairs:
    print("No filled models found — check FILLED_GLOBS paths.")
else:
    dropdown = W.Dropdown(options=pairs, description="PDB:", layout={"width": "260px"})
    gray_toggle = W.Checkbox(value=True, description="Grey crystal")
    btn = W.Button(description="Show in PyMOL", button_style="primary")
    out = W.Output()

    def _on_click(_):
        with out:
            out.clear_output()
            pymol_show_mmbsa_fill(dropdown.value, gray_crystal=gray_toggle.value)

    btn.on_click(_on_click)
    display(W.HBox([dropdown, gray_toggle, btn]), out)
"""

# ---------------------------------------------------------------------------
# Notebook 3 — Protonation PROPKA vs Literature
# ---------------------------------------------------------------------------

NB3_PARSE_CELL = """\
# Parsers: PROPKA .pka files + paper_evidence.md
import pathlib, re
from dataclasses import dataclass, field

DELIVERED_ROOT = pathlib.Path("/mnt/netapp1/Store_othcxlwa/newbench_27/delivered")

_PKA_HEADER_RE = re.compile(r"^\\s*Group\\s+", re.IGNORECASE)
_PKA_LINE_RE = re.compile(
    r"^\\s*"
    r"(?P<resn>[A-Z]{3})\\s+"          # residue name
    r"(?P<resi>\\-?\\d+)\\s+"           # residue number
    r"(?P<chain>[A-Za-z0-9])\\s+"       # chain id
    r"(?P<pka>\\-?\\d+\\.\\d+)"          # pKa
)

@dataclass
class PkaEntry:
    chain: str
    resi: int
    resn: str
    pka: float

def parse_pka_file(path: pathlib.Path) -> list[PkaEntry]:
    if not path.is_file():
        return []
    out = []
    for line in path.read_text().splitlines():
        m = _PKA_LINE_RE.match(line)
        if m:
            try:
                out.append(PkaEntry(
                    chain=m.group("chain"), resi=int(m.group("resi")),
                    resn=m.group("resn"),   pka=float(m.group("pka")),
                ))
            except (ValueError, KeyError):
                continue
    return out

def merge_propka_dir(protonation_dir: pathlib.Path) -> dict[tuple[str,int], list[PkaEntry]]:
    out: dict[tuple[str,int], list[PkaEntry]] = {}
    if not protonation_dir.is_dir():
        return out
    for pka_file in sorted(protonation_dir.glob("*.pka")):
        for e in parse_pka_file(pka_file):
            out.setdefault((e.chain, e.resi), []).append(e)
    return out


_RESI_HEADING_RE = re.compile(r"^###\\s+([A-Z]{3})(\\d+)\\s*$", re.IGNORECASE)
_QUOTE_LINE_RE = re.compile(r"^>\\s*(.*)$")

@dataclass
class PaperNote:
    resn: str
    resi: int
    quote: str

def parse_paper_evidence(md_path: pathlib.Path) -> list[PaperNote]:
    if not md_path.is_file():
        return []
    text = md_path.read_text()
    notes: list[PaperNote] = []
    lines = text.splitlines()
    i = 0
    while i < len(lines):
        m = _RESI_HEADING_RE.match(lines[i].strip())
        if m:
            resn, resi = m.group(1).upper(), int(m.group(2))
            j = i + 1
            quotes = []
            while j < len(lines) and not _RESI_HEADING_RE.match(lines[j].strip()) \\
                                  and not lines[j].startswith("## "):
                q = _QUOTE_LINE_RE.match(lines[j])
                if q:
                    quotes.append(q.group(1).strip())
                j += 1
            i = j
            if quotes:
                notes.append(PaperNote(resn=resn, resi=resi,
                                       quote=" ".join(quotes)[:400]))
            continue
        i += 1
    return notes
"""

NB3_LOADER_CELL = """\
# Discover proteins that have both PROPKA output and paper_evidence.md
pdb_ids = sorted(p.name for p in DELIVERED_ROOT.iterdir() if p.is_dir())
usable = []
for pid in pdb_ids:
    entry = DELIVERED_ROOT / pid
    has_propka = (entry / "protonation").is_dir() \\
                 and any((entry / "protonation").glob("*.pka"))
    has_paper = (entry / "paper" / "paper_evidence.md").is_file()
    if has_propka and has_paper:
        usable.append(pid)
print(f"{len(usable)} of {len(pdb_ids)} delivered proteins have both PROPKA + paper_evidence:")
print(" ".join(usable))
"""

NB3_RENDER_CELL = """\
# Side-by-side renderer: for one PDB, PROPKA table (top) + paper snippets (bottom)
from IPython.display import Markdown, display

def render_protonation_side_by_side(pdb_id: str):
    entry = DELIVERED_ROOT / pdb_id
    propka = merge_propka_dir(entry / "protonation")
    paper = parse_paper_evidence(entry / "paper" / "paper_evidence.md")

    md_lines = [f"## {pdb_id} — PROPKA vs Literature", ""]

    # Union of residues from both sources
    propka_keys = set(propka.keys())
    paper_keys = {(("A"), p.resi) for p in paper}  # paper_evidence has no chain, assume A
    all_keys = sorted(propka_keys | paper_keys, key=lambda k: (k[0], k[1]))

    md_lines.append("| Chain | Resi | Resn | PROPKA pKa (min–max) | Paper snippet |")
    md_lines.append("|---|---:|---|---|---|")
    paper_by_resi = {p.resi: p for p in paper}
    for (chain, resi) in all_keys:
        pkas = propka.get((chain, resi), [])
        resn = pkas[0].resn if pkas else (paper_by_resi[resi].resn if resi in paper_by_resi else "—")
        if pkas:
            pks = [f"{p.pka:.2f}" for p in pkas]
            pka_str = pks[0] if len(pks) == 1 else f"{min(p.pka for p in pkas):.2f}–{max(p.pka for p in pkas):.2f}"
        else:
            pka_str = "—"
        snippet = paper_by_resi[resi].quote if resi in paper_by_resi else "—"
        snippet = snippet.replace("|", "\\\\|")[:120] + ("…" if resi in paper_by_resi and len(paper_by_resi[resi].quote) > 120 else "")
        md_lines.append(f"| {chain} | {resi} | {resn} | {pka_str} | {snippet} |")

    display(Markdown("\\n".join(md_lines)))

# Try the first available one:
if usable:
    render_protonation_side_by_side(usable[0])
"""

NB3_WIDGET_CELL = """\
import ipywidgets as W
from IPython.display import display, clear_output

if usable:
    dropdown = W.Dropdown(options=usable, description="PDB:", layout={"width": "260px"})
    btn = W.Button(description="Render", button_style="primary")
    out = W.Output()

    def _on_click(_):
        with out:
            clear_output()
            render_protonation_side_by_side(dropdown.value)

    btn.on_click(_on_click)
    display(W.HBox([dropdown, btn]), out)
else:
    print("No proteins have both PROPKA + paper_evidence yet.")
"""

NB3_ACTIVE_SITE_MD_CELL = """\
# Optional: display the pre-generated active_site_protonation.md for one protein
def show_active_site_md(pdb_id: str):
    p = DELIVERED_ROOT / pdb_id / "active_site_protonation.md"
    if not p.is_file():
        print(f"No active_site_protonation.md for {pdb_id}")
        return
    display(Markdown(p.read_text()))

if usable:
    show_active_site_md(usable[0])
"""

# ---------------------------------------------------------------------------
# Notebook builder
# ---------------------------------------------------------------------------


def _md(txt: str) -> nbf.NotebookNode:
    return nbf.v4.new_markdown_cell(txt)

def _code(txt: str) -> nbf.NotebookNode:
    return nbf.v4.new_code_cell(txt)


def build_nb_1() -> nbf.NotebookNode:
    nb = nbf.v4.new_notebook()
    nb.cells = [
        _md(SETUP_MD.format(title="newbench_27 — Residue-Fill Review + Overview via PyMOL")),
        _md("## Step 1 — Connect to laptop PyMOL"),
        _code(CONNECT_CELL),
        _md("## Step 2 — Discover delivered newbench_27 proteins"),
        _code(DISCOVER_NEWBENCH_CELL),
        _md("## Step 3 — Per-protein summary (KI-Aufteilung, Masse, Actives)"),
        _md("Parses every `MANIFEST.json` + `graft/*_graft_manifest.json` under "
            "`delivered/` and builds one row per protein with the AI source "
            "(AF-splice vs BL-Pose fallback), estimated molecular weight, gap "
            "counts, and any warnings the pipeline logged."),
        _code(NB1_OVERVIEW_PARSE_CELL),
        _md("### 3a — KI-Aufteilung (AF-splice vs BL-Pose fallback)"),
        _code(NB1_AI_SPLIT_CELL),
        _md("### 3b — Massenaufteilung (molecular-weight histogram)"),
        _code(NB1_MASS_HISTOGRAM_CELL),
        _md("### 3c — Actives / inactives ligand sets per target\n\n"
            "Scans the usual on-disk locations. If your ligand sets live "
            "elsewhere, extend `ACTIVES_CANDIDATE_ROOTS`."),
        _code(NB1_ACTIVES_INACTIVES_CELL),
        _md("### 3d — Merged summary table"),
        _code(NB1_MERGED_TABLE_CELL),
        _md("## Step 4 — Loader function\n\n"
            "Loads the crystal PDB + FRUTON receptor (with graft) into PyMOL, "
            "greys the crystal, colours the filled model sky-blue, and highlights "
            "the graft-region residues in hot pink with side-chain sticks."),
        _code(NB1_LOADER_CELL),
        _md("## Step 5 — Interactive selector"),
        _code(NB1_WIDGET_CELL),
        _md("## Optional — Batch walk (saves one `.pse` session per PDB)"),
        _code(NB1_BATCH_CELL),
    ]
    return nb


def build_nb_2() -> nbf.NotebookNode:
    nb = nbf.v4.new_notebook()
    nb.cells = [
        _md(SETUP_MD.format(title="MMBSA_200 — Residue-Fill Review via PyMOL")),
        _md("## Step 1 — Connect to laptop PyMOL"),
        _code(CONNECT_CELL),
        _md("## Step 2 — Discover MMBSA_200 pairs (crystal + filled)"),
        _code(DISCOVER_MMBSA200_CELL),
        _md("## Step 3 — Loader\n\n"
            "MMBSA_200 does not (yet) ship a graft-manifest listing gap ranges, "
            "so filled-region residues are inferred as anything in the filled "
            "model that isn't within 0.5 Å of a crystal atom."),
        _code(NB2_LOADER_CELL),
        _md("## Step 4 — Interactive selector"),
        _code(NB2_WIDGET_CELL),
    ]
    return nb


def build_nb_3() -> nbf.NotebookNode:
    nb = nbf.v4.new_notebook()
    nb.cells = [
        _md(SETUP_MD.format(title="Protonation Review — PROPKA vs Literature")),
        _md("*(PyMOL is not required for this notebook — the connection cell "
            "is still included so you can quickly load a residue afterwards.)*"),
        _md("## Step 1 — (Optional) connect to PyMOL for cross-checking"),
        _code(CONNECT_CELL),
        _md("## Step 2 — Parsers"),
        _code(NB3_PARSE_CELL),
        _md("## Step 3 — Discover usable proteins (need both PROPKA output + `paper_evidence.md`)"),
        _code(NB3_LOADER_CELL),
        _md("## Step 4 — Render side-by-side"),
        _code(NB3_RENDER_CELL),
        _md("## Step 5 — Interactive selector"),
        _code(NB3_WIDGET_CELL),
        _md("## Step 6 — Show the pre-generated `active_site_protonation.md`\n\n"
            "This is the shipped summary FRUTON already wrote per protein."),
        _code(NB3_ACTIVE_SITE_MD_CELL),
    ]
    return nb


def main() -> None:
    HERE.mkdir(parents=True, exist_ok=True)
    for name, nb in (
        ("1_review_newbench27_fills.ipynb", build_nb_1()),
        ("2_review_mmbsa200_fills.ipynb",   build_nb_2()),
        ("3_review_protonation_side_by_side.ipynb", build_nb_3()),
    ):
        path = HERE / name
        with path.open("w") as fh:
            nbf.write(nb, fh)
        print(f"wrote {path}")


if __name__ == "__main__":
    main()
