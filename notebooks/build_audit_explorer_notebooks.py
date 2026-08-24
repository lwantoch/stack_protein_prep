#!/usr/bin/env python3
"""Generate 3 juplaunch-ready audit-explorer notebooks — one per bench set.

    notebooks/4_mmbsa200_audit_explorer.ipynb         (train + test)
    notebooks/5_new_benchmark_audit_explorer.ipynb    (affinity_bench_27)
    notebooks/6_random30_audit_explorer.ipynb         (stresstest_30)

Each loads its audit_report.csv from $LUSTRE, shows tier distribution,
per-protein action items, component-type breakdown, and offers an
optional PyMOL-remote spot-check for any protein by id.

Run once (or re-run after editing this script) to (re)generate the
notebook files.  Kept as a script so cell content stays version-
controlled — hand-edited .ipynb files diff badly.
"""
from __future__ import annotations

from pathlib import Path

import nbformat as nbf

HERE = Path(__file__).resolve().parent


# ---------------------------------------------------------------------------
# Per-notebook config
# ---------------------------------------------------------------------------

NOTEBOOKS = [
    {
        "filename": "4_mmbsa200_audit_explorer.ipynb",
        "title": "MMBSA_200 audit explorer — train + test",
        "subtitle": (
            "22 AF-ready train (MMBSA_75) + 26 AF-ready test (MMBSA_125). "
            "Bench iter 7, 2026-08-24."
        ),
        "csv_paths": [
            "$LUSTRE/MMBSA_200/generic_train_20260824_0848/artefacts/audit_report.csv",
            "$LUSTRE/MMBSA_200/generic_test_20260824_0848/artefacts/audit_report.csv",
        ],
        "split_labels": ["train", "test"],
        "crystal_root": "/mnt/netapp1/Store_othcxlwa/FRUTON-NEW",
    },
    {
        "filename": "5_new_benchmark_audit_explorer.ipynb",
        "title": "New affinity benchmark audit explorer — newbench_27",
        "subtitle": (
            "27 PDBs (1 AF-ready + 26 BL-Pose fallback). "
            "Bench iter 7, 2026-08-24."
        ),
        "csv_paths": [
            "$LUSTRE/MMBSA_200/generic_affinity_bench_27_20260824_0918/artefacts/audit_report.csv",
        ],
        "split_labels": ["affinity_bench_27"],
        "crystal_root": "/mnt/netapp1/Store_othcxlwa/newbench_27",
    },
    {
        "filename": "6_random30_audit_explorer.ipynb",
        "title": "Stress test audit explorer — 30 random from MMBSA_200",
        "subtitle": (
            "30 randomly-picked (seed 20260824) PDBs: 5 AF-ready + 25 BL-Pose. "
            "Bench iter 7, 2026-08-24."
        ),
        "csv_paths": [
            "$LUSTRE/MMBSA_200/generic_stresstest_30_20260824_0903/artefacts/audit_report.csv",
        ],
        "split_labels": ["stresstest_30"],
        "crystal_root": "/mnt/netapp1/Store_othcxlwa/FRUTON-NEW",
    },
]


# ---------------------------------------------------------------------------
# Cell templates
# ---------------------------------------------------------------------------

def header_md(title: str, subtitle: str) -> str:
    return (
        f"# {title}\n\n"
        f"{subtitle}\n\n"
        "Runs on CESGA via `juplaunch` (memory `[[reference_juplaunch]]`).\n\n"
        "## Optional: connect to laptop PyMOL for spot-checks\n\n"
        "On your **laptop** before running this notebook:\n"
        "```bash\n"
        "export PYMOL_PATH=/usr/lib/python3/dist-packages/pymol\n"
        "pymol -R                            # RPC on :9123\n"
        "ssh -R 9123:localhost:9123 ft3.cesga.es\n"
        "```\n\n"
        "On CESGA: `~/bin/juplaunch` → paste the URL into your laptop browser.\n"
        "PyMOL connection is optional — the notebook works without it.\n"
    )


IMPORTS_CELL = """\
import os, json, glob, pathlib
import pandas as pd
import matplotlib.pyplot as plt

pd.set_option("display.max_rows", 200)
pd.set_option("display.max_columns", 60)
pd.set_option("display.width", 200)
"""


def load_cell(csv_paths: list[str], split_labels: list[str]) -> str:
    """Emit a cell that reads all listed CSVs + tags each row with split."""
    parts = ["frames = []"]
    for path, label in zip(csv_paths, split_labels):
        parts.append(
            f'p = os.path.expandvars("{path}")\n'
            f'if pathlib.Path(p).is_file():\n'
            f'    df = pd.read_csv(p)\n'
            f'    df["split"] = "{label}"\n'
            f'    frames.append(df)\n'
            f'    print(f"loaded {{len(df)}} rows from {label}")\n'
            f'else:\n'
            f'    print(f"MISSING {{p}}  (bench dir absent — did the SLURM job finish?)")'
        )
    parts.append(
        "\n"
        "df = pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()\n"
        'print(f"\\nTOTAL rows: {len(df)}")\n'
        "df.head()"
    )
    return "\n".join(parts)


TIER_CELL = """\
# Tier distribution per split
if len(df):
    tiers = df.groupby(["split", "overall_status"]).size().unstack(fill_value=0)
    tiers["total"] = tiers.sum(axis=1)
    for col in ("delivered_full_confidence", "delivered_with_notes",
                "delivered_needs_review", "not_delivered"):
        if col not in tiers.columns:
            tiers[col] = 0
    tiers["delivered_ok"] = tiers["delivered_full_confidence"] + tiers["delivered_with_notes"]
    tiers["ok_pct"] = (100 * tiers["delivered_ok"] / tiers["total"]).round(1)
    display(tiers[[
        "delivered_full_confidence", "delivered_with_notes",
        "delivered_needs_review", "not_delivered",
        "total", "delivered_ok", "ok_pct",
    ]])
"""

BARPLOT_CELL = """\
# Tier stacked bar
if len(df):
    fig, ax = plt.subplots(figsize=(7, 3.5))
    tier_order = ["delivered_full_confidence", "delivered_with_notes",
                  "delivered_needs_review", "not_delivered"]
    colours = ["#2b6cb0", "#ecc94b", "#dd6b20", "#c53030"]
    grouped = df.groupby(["split", "overall_status"]).size().unstack(fill_value=0)
    for col in tier_order:
        if col not in grouped.columns:
            grouped[col] = 0
    grouped = grouped[tier_order]
    grouped.plot(kind="barh", stacked=True, color=colours, ax=ax, edgecolor="none")
    ax.set_xlabel("# proteins")
    ax.set_title("Tier distribution per split")
    ax.legend(loc="lower right", fontsize=8)
    fig.tight_layout()
    plt.show()
"""

NEEDS_REVIEW_CELL = """\
# Proteins flagged for manual review (with action items)
if len(df):
    nr = df[df["overall_status"] == "delivered_needs_review"]
    if len(nr):
        print(f"{len(nr)} protein(s) need manual review:")
        display(nr[["pdb", "split", "n_gaps", "delta_n", "clash_pairs",
                    "omega_np", "omega_cis_nonpro", "rama_outliers",
                    "gate_pass", "action_items", "notes"]])
    else:
        print("Zero proteins flagged for manual review in this set.")
"""

WITH_NOTES_CELL = """\
# Proteins delivered_with_notes — action items reviewers should glance at
if len(df):
    wn = df[df["overall_status"] == "delivered_with_notes"]
    if len(wn):
        print(f"{len(wn)} protein(s) delivered_with_notes:")
        # Show shortened action_items column so the row fits
        display(
            wn[["pdb", "split", "n_high", "n_medium", "n_low", "n_failed", "action_items"]]
              .rename(columns={"action_items": "action_items (glance)"})
        )
    else:
        print("Zero proteins with notes.")
"""

COMPONENT_MIX_CELL = """\
# Component-type × confidence heatmap
# Requires re-reading the per-PDB JSONs (audit CSV aggregates but drops component detail).
def _load_components(csv_path: str) -> pd.DataFrame:
    art_dir = pathlib.Path(os.path.expandvars(csv_path)).parent
    rows = []
    for f in sorted(art_dir.glob("*.json")):
        if f.name.startswith("audit_") or f.name == "combined_results.json":
            continue
        try:
            d = json.loads(f.read_text())
        except Exception:
            continue
        pid = d.get("pdb", f.stem)
        for c in d.get("delivery", {}).get("components", []):
            rows.append({
                "pdb": pid,
                "component_type": c.get("component_type"),
                "confidence": c.get("confidence"),
                "method": c.get("method", ""),
                "reason": (c.get("reason") or "")[:120],
                "suggested_action": (c.get("suggested_action") or "")[:120],
            })
    return pd.DataFrame(rows)


comps_all = []
for p in CSV_PATHS:
    p = os.path.expandvars(p)
    cdf = _load_components(p)
    if len(cdf):
        cdf["source"] = pathlib.Path(p).parent.parent.name
        comps_all.append(cdf)

comps = pd.concat(comps_all, ignore_index=True) if comps_all else pd.DataFrame()
if len(comps):
    ct = comps.pivot_table(
        index="component_type", columns="confidence", aggfunc="size", fill_value=0,
    )
    for col in ("high", "medium", "low", "failed"):
        if col not in ct.columns:
            ct[col] = 0
    ct = ct[["high", "medium", "low", "failed"]]
    ct["total"] = ct.sum(axis=1)
    print("Component × confidence breakdown:")
    display(ct)
"""

METHOD_MIX_CELL = """\
# Most common per-component method labels
if len(comps):
    method_counts = (
        comps.groupby(["component_type", "method"]).size()
             .reset_index(name="n")
             .sort_values("n", ascending=False)
             .head(30)
    )
    display(method_counts)
"""

DRILLDOWN_CELL = """\
# Per-protein drill-down: pick any PDB id from the audit table
PDB_ID = "8Q68"   # ← change me; e.g. one of the needs_review pdbs shown above

if len(df):
    row = df[df["pdb"] == PDB_ID]
    if not len(row):
        print(f"{PDB_ID} not in this bench set")
    else:
        print(f"=== {PDB_ID} ({row.iloc[0]['split']}) ===")
        for col in ("overall_status", "n_gaps", "delta_n", "rolled_gaps",
                    "clash_pairs", "omega_np", "omega_cis_nonpro", "rama_outliers",
                    "gate_pass", "action_items", "notes"):
            if col in row.columns:
                print(f"  {col:<20}: {row.iloc[0][col]}")
        # Component detail from JSON side-load
        if len(comps):
            pc = comps[comps["pdb"] == PDB_ID]
            if len(pc):
                print("\\n  Components:")
                display(pc[["component_type", "confidence", "method", "reason", "suggested_action"]])
"""

PYMOL_CONNECT_CELL = """\
# Connect to laptop PyMOL via reverse tunnel (localhost:9123)
# Requires: pymol -R running on laptop + ssh -R 9123:localhost:9123 to ft3
from pymol_remote.client import PymolSession
pm = PymolSession(hostname="localhost", port=9123, timeout=10.0)
print("PyMOL connected. Loaded objects:", pm.get_names())
"""


PYMOL_HELPERS_CELL = """\
# --- Visual review helpers -------------------------------------------------
# All functions render into the laptop PyMOL GUI via the reverse tunnel.
# Each helper reads the per-PDB JSON emitted by FRUTON's generic bench
# driver + the actual crystal/final PDB on disk.
import pathlib, json, os

def _pdb_json(pdb_id: str) -> dict:
    for csv_path in CSV_PATHS:
        art_dir = pathlib.Path(os.path.expandvars(csv_path)).parent
        f = art_dir / f"{pdb_id}.json"
        if f.is_file():
            return json.loads(f.read_text())
    return {}


def _crystal_path(pdb_id: str) -> pathlib.Path | None:
    p = pathlib.Path(CRYSTAL_ROOT) / pdb_id / f"{pdb_id}.pdb"
    return p if p.is_file() else None


def _final_pdb_path(pdb_id: str) -> pathlib.Path | None:
    # Bench driver stores per-protein artefact PDBs under artefacts dir
    for csv_path in CSV_PATHS:
        art_dir = pathlib.Path(os.path.expandvars(csv_path)).parent
        # Match either <pdb>_final.pdb (splice path) or the crystal itself
        # (BL-Pose fallback ships crystal-as-is).
        cand = list(art_dir.parent.rglob(f"{pdb_id}_final.pdb"))
        if cand: return cand[0]
    return _crystal_path(pdb_id)


def _gap_ranges_from_json(d: dict) -> list[tuple[str, int, int]]:
    # Try to recover per-gap boundaries from the per-PDB record.  The
    # current bench driver emits n_gaps + rolled_gaps at protein level
    # but not per-gap ranges; fall back to '(chain, resnum-1, resnum+1)'
    # around REMARK 465 residues in the crystal.
    ranges = []
    crystal = _crystal_path(d.get("pdb", ""))
    if crystal:
        current: tuple[str, int, int] | None = None
        for line in crystal.read_text(errors="replace").splitlines():
            if not line.startswith("REMARK 465"):
                continue
            parts = line.split()
            if len(parts) < 5: continue
            try:
                resnum = int(parts[-1])
                chain = parts[-2]
            except ValueError:
                continue
            if current is None:
                current = (chain, resnum, resnum)
            elif chain == current[0] and resnum == current[2] + 1:
                current = (chain, current[1], resnum)
            else:
                ranges.append(current)
                current = (chain, resnum, resnum)
        if current: ranges.append(current)
    return ranges


def show_gaps(pdb_id: str):
    \"\"\"Load crystal (grey) + FRUTON final (blue), highlight gap residues red.\"\"\"
    d = _pdb_json(pdb_id)
    crystal = _crystal_path(pdb_id)
    final = _final_pdb_path(pdb_id)
    pm.do("reinitialize")
    pm.do("bg_color white")
    if crystal:
        pm.do(f"load {crystal}, crystal")
        pm.do("color grey70, crystal"); pm.do("show cartoon, crystal")
    if final and final != crystal:
        pm.do(f"load {final}, final")
        pm.do("color skyblue, final"); pm.do("show cartoon, final")
    # Highlight gap-region residues in the final model
    for chain, lo, hi in _gap_ranges_from_json(d):
        sel = f"final and chain {chain} and resi {lo}-{hi}"
        pm.do(f"color red, {sel}")
        pm.do(f"show sticks, {sel} and not name C+N+O+CA")
    pm.do("zoom polymer")
    print(f"{pdb_id}: crystal={crystal.name if crystal else 'None'}, "
          f"final={final.name if final else 'None'}, "
          f"gap_ranges={_gap_ranges_from_json(d)}")


def show_active_site_protonation(pdb_id: str, cutoff: float = 5.0):
    \"\"\"Colour residues around metals / active-site by their protonation state.

    HID orange, HIE yellow, HIP red, CYM green, CYS grey,
    ASH cyan (protonated Asp), GLH cyan (protonated Glu),
    TYM salmon (tyrosinate), LYN magenta (neutral Lys).
    \"\"\"
    final = _final_pdb_path(pdb_id)
    if not final:
        print(f"no final pdb for {pdb_id}"); return
    pm.do("reinitialize"); pm.do("bg_color white")
    pm.do(f"load {final}, mdl")
    pm.do("color grey80, mdl"); pm.do("show cartoon, mdl")
    # Detect metal-containing HETATMs so we can zoom on the active site
    pm.do("show spheres, mdl and hetatm and elem Zn+Fe+Cu+Mg+Ca+Mn+Ni+Co+Mo")
    pm.do("color orange, mdl and resn HID")
    pm.do("color yellow, mdl and resn HIE")
    pm.do("color red,    mdl and resn HIP")
    pm.do("color green,  mdl and resn CYM")
    pm.do("color grey40, mdl and resn CYS")
    pm.do("color cyan,   mdl and resn ASH+GLH")
    pm.do("color salmon, mdl and resn TYM")
    pm.do("color magenta,mdl and resn LYN")
    # Sticks on any residue whose Cα is within cutoff of a metal
    pm.do(f"show sticks, mdl and byres (mdl and elem Zn+Fe+Cu+Mg+Ca+Mn+Ni+Co+Mo) around {cutoff}")
    pm.do("zoom mdl and elem Zn+Fe+Cu+Mg+Ca+Mn+Ni+Co+Mo, 10")
    print(f"{pdb_id}: shown active site with protonation colour code")
    print("   HID=orange, HIE=yellow, HIP=red, CYM=green, ASH/GLH=cyan, "
          "TYM=salmon, LYN=magenta")


def show_literature_reference(pdb_id: str):
    \"\"\"If a paper_evidence.md or lit_image.png exists next to the crystal,
    display it inline in the notebook.\"\"\"
    from IPython.display import display, Markdown, Image
    crystal = _crystal_path(pdb_id)
    if not crystal:
        print(f"no crystal for {pdb_id}"); return
    pdir = crystal.parent
    md = pdir / "paper_evidence.md"
    img = None
    for candidate in (pdir / "lit_image.png", pdir / "figure.png",
                      pdir / "paper_figure.png", pdir / "site.png"):
        if candidate.is_file():
            img = candidate; break
    if md.is_file():
        display(Markdown(f"### Paper evidence for {pdb_id}\\n"))
        display(Markdown(md.read_text()))
    else:
        print(f"no paper_evidence.md for {pdb_id} (expected at {md})")
    if img:
        display(Markdown(f"### Literature image for {pdb_id}"))
        display(Image(str(img)))
    else:
        print(f"no lit image found for {pdb_id} in {pdir}")


print("Helpers loaded:")
print("  show_gaps(pdb_id)                       — crystal grey + final blue, gaps red")
print("  show_active_site_protonation(pdb_id)    — protonation colour code around metals")
print("  show_literature_reference(pdb_id)       — inline paper_evidence + lit image")
"""


PYMOL_DEMO_CELL = """\
# Demo — change PDB_ID to any protein from the audit table.
# Requires PyMOL connection above.
DEMO_PDB = "8Q68"    # ← change me
try:
    show_gaps(DEMO_PDB)
except Exception as e:
    print(f"show_gaps failed: {e!r}")
"""


PYMOL_PROTONATION_DEMO_CELL = """\
try:
    show_active_site_protonation(DEMO_PDB)
except Exception as e:
    print(f"show_active_site_protonation failed: {e!r}")
"""


PYMOL_LIT_DEMO_CELL = """\
try:
    show_literature_reference(DEMO_PDB)
except Exception as e:
    print(f"show_literature_reference failed: {e!r}")
"""


def build_notebook(cfg: dict) -> nbf.NotebookNode:
    nb = nbf.v4.new_notebook()
    cells = [
        nbf.v4.new_markdown_cell(header_md(cfg["title"], cfg["subtitle"])),
        nbf.v4.new_code_cell(IMPORTS_CELL),
        # Set CSV_PATHS + CRYSTAL_ROOT as module-level vars for later cells
        nbf.v4.new_code_cell(
            "CSV_PATHS = " + repr(cfg["csv_paths"]) + "\n"
            "CRYSTAL_ROOT = " + repr(cfg["crystal_root"]) + "\n"
            "print('CSV_PATHS:', CSV_PATHS)\n"
            "print('CRYSTAL_ROOT:', CRYSTAL_ROOT)"
        ),
        nbf.v4.new_markdown_cell("## Load audit CSV(s)"),
        nbf.v4.new_code_cell(load_cell(cfg["csv_paths"], cfg["split_labels"])),
        nbf.v4.new_markdown_cell("## Tier distribution"),
        nbf.v4.new_code_cell(TIER_CELL),
        nbf.v4.new_code_cell(BARPLOT_CELL),
        nbf.v4.new_markdown_cell("## Needs-review triage"),
        nbf.v4.new_code_cell(NEEDS_REVIEW_CELL),
        nbf.v4.new_markdown_cell("## Delivered-with-notes glance-list"),
        nbf.v4.new_code_cell(WITH_NOTES_CELL),
        nbf.v4.new_markdown_cell("## Component-type × confidence breakdown"),
        nbf.v4.new_code_cell(COMPONENT_MIX_CELL),
        nbf.v4.new_code_cell(METHOD_MIX_CELL),
        nbf.v4.new_markdown_cell("## Per-protein drill-down"),
        nbf.v4.new_code_cell(DRILLDOWN_CELL),
        nbf.v4.new_markdown_cell(
            "## Visual review — PyMOL windows\n\n"
            "The following cells render into the **laptop PyMOL GUI** via the "
            "reverse tunnel (`-R 9123:localhost:9123`).  Skip if you only want "
            "the audit CSV analysis.\n\n"
            "Provides three helpers:\n\n"
            "- `show_gaps(pdb_id)` — crystal (grey) + FRUTON final (blue), gap-fill residues in red\n"
            "- `show_active_site_protonation(pdb_id)` — final model, protonation colour-coded around metals\n"
            "- `show_literature_reference(pdb_id)` — inline paper_evidence.md + lit image if present\n"
        ),
        nbf.v4.new_code_cell(PYMOL_CONNECT_CELL),
        nbf.v4.new_code_cell(PYMOL_HELPERS_CELL),
        nbf.v4.new_markdown_cell("### Demo: reproduced gaps in PyMOL"),
        nbf.v4.new_code_cell(PYMOL_DEMO_CELL),
        nbf.v4.new_markdown_cell("### Demo: active-site protonation colour-code"),
        nbf.v4.new_code_cell(PYMOL_PROTONATION_DEMO_CELL),
        nbf.v4.new_markdown_cell("### Demo: paper evidence + literature image (inline)"),
        nbf.v4.new_code_cell(PYMOL_LIT_DEMO_CELL),
    ]
    nb["cells"] = cells
    return nb


def main() -> int:
    for cfg in NOTEBOOKS:
        nb = build_notebook(cfg)
        out = HERE / cfg["filename"]
        with out.open("w") as fh:
            nbf.write(nb, fh)
        print(f"wrote {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
