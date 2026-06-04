"""Chart rendering and data extraction helpers for the FRUTON pipeline summary report."""
from __future__ import annotations

import io

from ._report_palette import (
    _IDIS_NAVY,
    _IDIS_TEXT,
    _IDIS_MUTED,
    _IDIS_LINE,
    _IDIS_PANEL,
    _IDIS_GOLD,
    _FRUTON_RED,
    _FRUTON_GREEN,
    _FRUTON_GREEN_SOFT,
    _AUDIT_WARN,
    _AUDIT_FAIL,
    _AUDIT_OK,
    _TRANSITION_METAL_ELEMENTS,
)


def _parse_residue_count(range_str: str) -> int:
    """Parse a residue range string like '16-243' or 'A16-B50' to residue count.

    The summary chart uses this number only as a scale cue, not as a strict
    biological residue count. FRUTON range fields may contain plain numeric
    ranges, chain-prefixed ranges, or empty values depending on how far the
    input record travelled through the pipeline. Invalid or absent values are
    converted to zero so the report can still be written for failed or partial
    runs. This keeps report generation independent from upstream validation.
    """
    if not range_str or "-" not in range_str:
        return 0
    parts = range_str.split("-", 1)
    if len(parts) != 2:
        return 0
    try:
        start_s = parts[0].lstrip("ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz ")
        end_s = parts[1].lstrip("ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz ")
        if start_s and end_s:
            return max(0, int(end_s) - int(start_s) + 1)
    except (ValueError, TypeError):
        pass
    return 0


def _extract_chart_data(pipeline_record_list: list[dict]) -> list[dict]:
    rows: list[dict] = []
    for rec in pipeline_record_list:
        pdb_id = str(rec.get("pdb_id", "")).strip()
        if not pdb_id:
            continue

        has_insertions = str(rec.get("insertion_codes_done", "")).strip().lower() == "success"

        try:
            n_gaps = int(str(rec.get("n_gaps", "0") or "0"))
        except (ValueError, TypeError):
            n_gaps = 0

        gap_sizes_raw = str(rec.get("gap_sizes", "") or "")
        gap_sizes_list: list[int] = []
        if gap_sizes_raw and gap_sizes_raw.lower() != "none":
            for part in gap_sizes_raw.split("|"):
                part = part.strip()
                if part and part.lower() != "none":
                    try:
                        gap_sizes_list.append(int(part))
                    except (ValueError, TypeError):
                        pass
        total_gap_length = sum(gap_sizes_list)

        has_cofactors = str(rec.get("has_cofactors", "") or "").strip().lower() == "yes"
        has_nonstd = str(rec.get("has_nonstandard_residues", "") or "").strip().lower() == "yes"

        has_metals_flag = str(rec.get("has_metals", "") or "").strip().lower() == "yes"
        has_transition_metal = False
        if has_metals_flag:
            ion_type_raw = str(rec.get("metals.ion_type", "") or "").strip()
            if ion_type_raw:
                elements = {
                    e.strip().upper()
                    for e in ion_type_raw.replace(",", "|").split("|")
                    if e.strip()
                }
                has_transition_metal = bool(elements & _TRANSITION_METAL_ELEMENTS)
            else:
                has_transition_metal = True

        n_residues = _parse_residue_count(str(rec.get("range", "") or ""))

        rows.append(
            {
                "pdb_id": pdb_id,
                "has_insertions": has_insertions,
                "n_gaps": n_gaps,
                "gap_sizes_list": gap_sizes_list,
                "total_gap_length": total_gap_length,
                "has_cofactors": has_cofactors,
                "has_nonstd": has_nonstd,
                "has_transition_metal": has_transition_metal,
                "n_residues": n_residues,
            }
        )
    return rows


def _render_chart_pngs(
    chart_data: list[dict],
    *,
    landscape: bool = False,
    chunk_size: int = 20,
) -> list[bytes]:
    """Render FRUTON-style per-protein feature charts as PNG byte strings.

    The chart intentionally avoids grid lines so it visually matches the cleaner
    FRUTON report pages. Residue count is shown as a pale navy background bar on
    a secondary axis, while audit features use the same navy, red, gold, green,
    and muted tones used by the per-protein reports and manual. The landscape
    chunking behavior is kept unchanged: large datasets are split into rotated
    pages so labels remain readable. The returned PNGs are temporary rendering
    artefacts and are not part of the pipeline state.
    """
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.patches as mpatches
    import matplotlib.pyplot as plt
    import numpy as np

    if not chart_data:
        fig, ax = plt.subplots(figsize=(8, 4))
        fig.patch.set_facecolor("white")
        ax.set_facecolor("white")
        ax.text(
            0.5,
            0.5,
            "No data",
            ha="center",
            va="center",
            color=_IDIS_MUTED,
            transform=ax.transAxes,
        )
        ax.set_xticks([])
        ax.set_yticks([])
        for spine in ax.spines.values():
            spine.set_visible(False)
        buf = io.BytesIO()
        fig.savefig(buf, format="png", dpi=180, bbox_inches="tight")
        plt.close(fig)
        return [buf.getvalue()]

    c_ins = _IDIS_GOLD
    c_gapn = _IDIS_NAVY
    c_gapr = "#B8A85A"
    c_met = _FRUTON_RED
    c_cof = _FRUTON_GREEN
    c_nst = "#A33A46"
    c_bg = _IDIS_NAVY

    bar_w = 0.13
    bar_gap = 0.015
    offsets = np.array([-2, -1, 0, 1, 2]) * (bar_w + bar_gap)
    bg_bar_w = 5 * bar_w + 4 * bar_gap + 0.06

    n_total = len(chart_data)
    chunks = (
        [chart_data[i : i + chunk_size] for i in range(0, n_total, chunk_size)]
        if landscape
        else [chart_data]
    )

    results: list[bytes] = []
    for chunk_idx, chunk in enumerate(chunks):
        nc = len(chunk)
        pdb_ids = [d["pdb_id"] for d in chunk]
        x = np.arange(nc, dtype=float)

        ins = np.array([int(d["has_insertions"]) for d in chunk], dtype=float)
        met = np.array([int(d["has_transition_metal"]) for d in chunk], dtype=float)
        cof = np.array([int(d["has_cofactors"]) for d in chunk], dtype=float)
        nst = np.array([int(d["has_nonstd"]) for d in chunk], dtype=float)
        nres = np.array([d["n_residues"] for d in chunk], dtype=float)
        gap_layers = [d["gap_sizes_list"] for d in chunk]
        max_gap_layers = max((len(g) for g in gap_layers), default=0)
        # Colour ramp: dark navy -> lighter steel blue, one shade per gap layer
        _navy = np.array([0x28 / 255, 0x32 / 255, 0x5A / 255])
        _light = np.array([0x7A / 255, 0xA0 / 255, 0xD0 / 255])

        def _gap_color(layer_idx: int) -> str:
            t = layer_idx / max(max_gap_layers - 1, 1)
            rgb = _navy + t * (_light - _navy)
            return "#{:02X}{:02X}{:02X}".format(
                int(rgb[0] * 255), int(rgb[1] * 255), int(rgb[2] * 255)
            )

        total_gap = np.array([d["total_gap_length"] for d in chunk], dtype=float)

        fig_w = max(14, nc * 0.72) if landscape else max(9, nc * 1.3)
        fig_h = 5.5 if landscape else 5.2

        fig, ax = plt.subplots(figsize=(fig_w, fig_h))
        ax2 = ax.twinx()
        fig.patch.set_facecolor("white")
        ax.set_facecolor("white")

        ax2.set_zorder(1)
        ax.set_zorder(2)
        ax.patch.set_visible(False)

        ax2.bar(x, nres, bg_bar_w, alpha=0.12, color=c_bg, zorder=0)
        ax2.set_ylabel("Starting range residues", fontsize=8, color=_IDIS_MUTED, labelpad=6)
        ax2.tick_params(axis="y", labelcolor=_IDIS_MUTED, labelsize=7, length=0)
        ax2.spines["top"].set_visible(False)
        ax2.spines["left"].set_visible(False)
        ax2.spines["right"].set_color(_IDIS_LINE)
        ax2.spines["right"].set_linewidth(0.6)
        ax2.yaxis.grid(False)
        nres_max = max(float(nres.max()) if nc else 1.0, 1.0)
        ax2.set_ylim(0, nres_max * 1.55)

        ax.bar(x + offsets[0], ins, bar_w, color=c_ins, zorder=2)
        # Gaps: one stacked segment per individual gap, darkest at bottom
        gap_bottoms = np.zeros(nc)
        for layer_idx in range(max_gap_layers):
            heights = np.array(
                [g[layer_idx] if layer_idx < len(g) else 0 for g in gap_layers],
                dtype=float,
            )
            ax.bar(x + offsets[1], heights, bar_w, bottom=gap_bottoms,
                   color=_gap_color(layer_idx), zorder=2)
            gap_bottoms += heights
        ax.bar(x + offsets[2], met, bar_w, color=c_met, zorder=2)
        ax.bar(x + offsets[3], cof, bar_w, color=c_cof, zorder=2)
        ax.bar(x + offsets[4], nst, bar_w, color=c_nst, zorder=2)

        fs_tick = 7 if nc > 15 else 9
        ax.set_xticks(x)
        ax.set_xticklabels(pdb_ids, fontsize=fs_tick, rotation=45, ha="right")
        ax.set_ylabel("Count / presence", fontsize=9, color=_IDIS_TEXT, labelpad=6)
        ax.tick_params(axis="y", labelcolor=_IDIS_TEXT, labelsize=8, length=0)
        ax.tick_params(axis="x", colors=_IDIS_TEXT, length=0)

        feat_max = max(float(total_gap.max()) if nc else 1.0, 1.0)
        ax.set_ylim(0, feat_max * 1.55)

        suffix = f" - page {chunk_idx + 1}" if len(chunks) > 1 else ""
        ax.set_title(
            f"Per-protein feature summary - {n_total} protein(s){suffix}",
            fontsize=10,
            fontweight="bold",
            color=_IDIS_NAVY,
            pad=8,
        )

        ax.axhline(0, color=_IDIS_LINE, linewidth=0.8, zorder=1)
        ax.grid(False)
        ax.yaxis.grid(False)
        ax.xaxis.grid(False)
        ax.set_axisbelow(False)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.spines["left"].set_visible(False)
        ax.spines["bottom"].set_color(_IDIS_LINE)
        ax.spines["bottom"].set_linewidth(0.8)

        legend_handles = [
            mpatches.Patch(color=c_ins, label="Insertions"),
            mpatches.Patch(color=_gap_color(0), label="Gap 1 (length, darkest)"),
            mpatches.Patch(color=_gap_color(max(max_gap_layers - 1, 0)), label="Gap N (length, lightest)"),
            mpatches.Patch(color=c_met, label="Transition metals"),
            mpatches.Patch(color=c_cof, label="Cofactors"),
            mpatches.Patch(color=c_nst, label="Non-standard residues"),
            mpatches.Patch(color=c_bg, alpha=0.30, label="Starting range residues"),
        ]
        legend = ax.legend(
            handles=legend_handles,
            fontsize=7,
            loc="upper left",
            framealpha=0.96,
            ncol=2,
            borderpad=0.6,
        )
        legend.get_frame().set_edgecolor(_IDIS_LINE)
        legend.get_frame().set_linewidth(0.5)
        legend.get_frame().set_facecolor("white")

        fig.tight_layout()

        buf = io.BytesIO()
        fig.savefig(buf, format="png", dpi=180, bbox_inches="tight")
        plt.close(fig)
        results.append(buf.getvalue())

    return results


def _compute_stats(
    pipeline_record_list: list[dict],
    chart_data: list[dict],
    gaussian_pending_ids: list[str] | None,
) -> dict:
    total = len(pipeline_record_list)
    prepared = sum(
        1
        for r in pipeline_record_list
        if str(r.get("prepared_structure.status", "") or "").strip().lower() == "success"
    )
    failed = sum(
        1
        for r in pipeline_record_list
        if str(r.get("prepared_structure.status", "") or "").strip().lower()
        in {"failed", "error"}
    )
    with_gaps = sum(1 for d in chart_data if d["n_gaps"] > 0)
    with_metals = sum(1 for d in chart_data if d["has_transition_metal"])
    with_nonstd = sum(1 for d in chart_data if d["has_nonstd"])
    with_cofactors = sum(1 for d in chart_data if d["has_cofactors"])
    gaussian_pending = len(gaussian_pending_ids) if gaussian_pending_ids else 0
    return {
        "total": total,
        "prepared": prepared,
        "failed": failed,
        "with_gaps": with_gaps,
        "with_metals": with_metals,
        "with_nonstd": with_nonstd,
        "with_cofactors": with_cofactors,
        "gaussian_pending": gaussian_pending,
    }


def _status_background(metric: str, value: int) -> str:
    if metric in {"prepared"}:
        return _AUDIT_OK if value else _IDIS_PANEL
    if metric in {"failed", "with_nonstd", "gaussian_pending"}:
        return _AUDIT_FAIL if value else _IDIS_PANEL
    if metric in {"with_gaps", "with_metals"}:
        return _AUDIT_WARN if value else _IDIS_PANEL
    if metric in {"with_cofactors"}:
        return _FRUTON_GREEN_SOFT if value else _IDIS_PANEL
    return _IDIS_PANEL
