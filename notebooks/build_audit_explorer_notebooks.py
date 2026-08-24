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

PYMOL_CELL = """\
# OPTIONAL — connect to laptop PyMOL to spot-check any protein
# Requires the reverse tunnel (see top-of-notebook instructions).
try:
    from pymol_remote.client import PymolSession
    pm = PymolSession(hostname="localhost", port=9123, timeout=5.0)
    print("PyMOL connected. Loaded objects:", pm.get_names())

    def show_crystal(pdb_id: str):
        path = pathlib.Path(CRYSTAL_ROOT) / pdb_id / f"{pdb_id}.pdb"
        if not path.is_file():
            print(f"crystal not found: {path}")
            return
        pm.do("reinitialize")
        pm.do("bg_color white")
        pm.do(f"load {path}, crystal")
        pm.do("show cartoon")
        pm.do("color skyblue, crystal")
        print(f"loaded {path.name}")

    print("Usage: show_crystal('8Q68')  # to display any PDB in laptop PyMOL")
except Exception as e:
    print(f"PyMOL not available (reverse tunnel not set up?): {e!r}")
    print("Skip this cell if you don't need visual review.")
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
        nbf.v4.new_markdown_cell("## Optional — PyMOL spot-check"),
        nbf.v4.new_code_cell(PYMOL_CELL),
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
