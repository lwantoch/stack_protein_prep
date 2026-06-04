from __future__ import annotations

import io
import os
import tempfile
from datetime import datetime
from pathlib import Path
from typing import Any

_IDIS_NAVY = "#28325A"
_IDIS_NAVY_2 = "#323264"
_IDIS_TEXT = "#252A3F"
_IDIS_MUTED = "#697086"
_IDIS_LINE = "#D8DCE7"
_IDIS_PANEL = "#F5F7FB"
_IDIS_PANEL_2 = "#EEF1F7"
_IDIS_GOLD = "#8A730F"
_FRUTON_RED = "#E60014"
_FRUTON_RED_SOFT = "#FDECEE"
_FRUTON_YELLOW = "#FAE646"
_FRUTON_YELLOW_SOFT = "#FFF8CF"
_FRUTON_GREEN = "#50961E"
_FRUTON_GREEN_SOFT = "#EAF4E3"
_AUDIT_WARN = "#FFF3C4"
_AUDIT_FAIL = "#F9D7DB"
_AUDIT_OK = "#E7F3E1"

_FRUTON_EXPANSION = (
    "Framework for Reconstruction, UniProt alignment, and "
    "Topology-Oriented protein Normalization"
)

_TRANSITION_METAL_ELEMENTS = frozenset(
    {"ZN", "CU", "FE", "MN", "CO", "NI", "MO", "CD", "HG", "PD", "PT"}
)

_IDIS_LOGO_NAMES = (
    "logo_idis.png",
    "logo_IDIS_2020-1.png",
    "logo_IDIS_2020-1(1).png",
    "idis_logo.png",
)
_FRUTON_LOGO_NAMES = (
    "logo.png",
    "logo(1).png",
    "fruton_logo.png",
    "FRUTON_logo.png",
)


def _timestamp() -> str:
    return datetime.now().strftime("%Y-%m-%d %H:%M:%S")


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
        total_gap_length = 0
        if gap_sizes_raw and gap_sizes_raw.lower() != "none":
            for part in gap_sizes_raw.split("|"):
                part = part.strip()
                if part and part.lower() != "none":
                    try:
                        total_gap_length += int(part)
                    except (ValueError, TypeError):
                        pass

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
        gn = np.array([d["n_gaps"] for d in chunk], dtype=float)
        gr = np.array([d["total_gap_length"] for d in chunk], dtype=float)
        met = np.array([int(d["has_transition_metal"]) for d in chunk], dtype=float)
        cof = np.array([int(d["has_cofactors"]) for d in chunk], dtype=float)
        nst = np.array([int(d["has_nonstd"]) for d in chunk], dtype=float)
        nres = np.array([d["n_residues"] for d in chunk], dtype=float)

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
        ax.bar(x + offsets[1], gn, bar_w, color=c_gapn, zorder=2)
        ax.bar(x + offsets[1], gr, bar_w, bottom=gn, color=c_gapr, zorder=2)
        ax.bar(x + offsets[2], met, bar_w, color=c_met, zorder=2)
        ax.bar(x + offsets[3], cof, bar_w, color=c_cof, zorder=2)
        ax.bar(x + offsets[4], nst, bar_w, color=c_nst, zorder=2)

        fs_tick = 7 if nc > 15 else 9
        ax.set_xticks(x)
        ax.set_xticklabels(pdb_ids, fontsize=fs_tick, rotation=45, ha="right")
        ax.set_ylabel("Count / presence", fontsize=9, color=_IDIS_TEXT, labelpad=6)
        ax.tick_params(axis="y", labelcolor=_IDIS_TEXT, labelsize=8, length=0)
        ax.tick_params(axis="x", colors=_IDIS_TEXT, length=0)

        feat_max = max(float((gn + gr).max()) if nc else 1.0, 1.0)
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
            mpatches.Patch(color=c_gapn, label="Gaps (count)"),
            mpatches.Patch(color=c_gapr, label="Gap residues (stacked)"),
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


def _logo_search_roots(protein_data_dir: Path | None) -> list[Path]:
    module_dir = Path(__file__).resolve().parent
    repo_root = module_dir.parent.parent
    roots: list[Path] = []
    if protein_data_dir is not None:
        roots += [protein_data_dir, protein_data_dir.parent]
    roots += [
        repo_root / "data" / "assets",
        repo_root / "data",
        repo_root / "assets",
        repo_root,
        Path.cwd() / "data" / "assets",
        Path.cwd() / "data",
        Path.cwd(),
        module_dir / "assets",
        module_dir.parent / "assets",
    ]
    seen: set[Path] = set()
    unique: list[Path] = []
    for root in roots:
        try:
            resolved = root.resolve()
        except OSError:
            resolved = root
        if resolved not in seen:
            seen.add(resolved)
            unique.append(root)
    return unique


def _find_logo(
    env_var: str,
    names: tuple[str, ...],
    protein_data_dir: Path | None,
) -> Path | None:
    env_value = os.environ.get(env_var, "").strip()
    if env_value:
        p = Path(env_value).expanduser()
        if p.exists():
            return p
    for root in _logo_search_roots(protein_data_dir):
        for name in names:
            candidate = root / name
            if candidate.exists():
                return candidate
    return None


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


def generate_pipeline_summary_report(
    pipeline_record_list: list[dict],
    output_path: Path,
    *,
    gaussian_pending_ids: list[str] | None = None,
    protein_data_dir: Path | None = None,
    run_timestamp: str | None = None,
) -> dict:
    """Write a FRUTON-style PDF summary report for a full pipeline run.

    The function is intentionally importable and does not configure logging,
    parse command-line arguments, or inspect global pipeline state. It receives
    already-flattened pipeline records, computes only the small set of summary
    values needed for the PDF, and writes one ReportLab document. Page rotation
    logic is preserved: compact datasets keep the chart on the portrait summary
    page, while larger datasets are split into landscape chart pages. The return
    value follows the existing lightweight status dictionary used by the caller.
    """
    try:
        from reportlab.lib import colors
        from reportlab.lib.enums import TA_CENTER, TA_LEFT, TA_RIGHT
        from reportlab.lib.pagesizes import A4, landscape as rl_landscape
        from reportlab.lib.styles import ParagraphStyle
        from reportlab.lib.units import cm
        from reportlab.lib.utils import ImageReader
        from reportlab.platypus import (
            BaseDocTemplate,
            Frame,
            Image,
            NextPageTemplate,
            PageBreak,
            PageTemplate,
            Paragraph,
            Spacer,
            Table,
            TableStyle,
        )
    except ImportError as exc:
        return {"status": "failed", "message": f"reportlab not available: {exc}"}

    ts = run_timestamp or _timestamp()
    page_w, page_h = A4
    land_w, land_h = rl_landscape(A4)
    margin_x = 1.45 * cm
    bottom_margin = 1.20 * cm
    top_margin = 2.25 * cm

    idis_logo_path = _find_logo("IDIS_LOGO_PATH", _IDIS_LOGO_NAMES, protein_data_dir)
    fruton_logo_path = _find_logo("FRUTON_LOGO_PATH", _FRUTON_LOGO_NAMES, protein_data_dir)

    chart_data = _extract_chart_data(pipeline_record_list)
    stats = _compute_stats(pipeline_record_list, chart_data, gaussian_pending_ids)

    n_proteins = len(chart_data)
    use_landscape_charts = n_proteins > 15
    chart_chunk_size = 20

    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    usable_w = page_w - 2 * margin_x
    land_usable_w = land_w - 2 * margin_x

    style_kicker = ParagraphStyle(
        "FrutonKicker",
        fontName="Helvetica-Bold",
        fontSize=8,
        leading=10,
        textColor=colors.HexColor(_IDIS_GOLD),
        alignment=TA_LEFT,
    )
    style_title = ParagraphStyle(
        "FrutonTitle",
        fontName="Helvetica-Bold",
        fontSize=20,
        leading=23,
        textColor=colors.HexColor(_IDIS_NAVY),
        alignment=TA_LEFT,
        spaceAfter=2,
    )
    style_subtitle = ParagraphStyle(
        "FrutonSubtitle",
        fontName="Helvetica",
        fontSize=9,
        leading=11,
        textColor=colors.HexColor(_IDIS_MUTED),
        alignment=TA_LEFT,
    )
    style_normal = ParagraphStyle(
        "FrutonNormal",
        fontName="Helvetica",
        fontSize=8.5,
        leading=10.5,
        textColor=colors.HexColor(_IDIS_TEXT),
        alignment=TA_LEFT,
    )
    style_section = ParagraphStyle(
        "FrutonSection",
        fontName="Helvetica-Bold",
        fontSize=13,
        leading=15,
        textColor=colors.HexColor(_IDIS_NAVY),
        alignment=TA_LEFT,
    )
    style_warning_title = ParagraphStyle(
        "FrutonWarnTitle",
        fontName="Helvetica-Bold",
        fontSize=10,
        leading=12,
        textColor=colors.HexColor(_IDIS_NAVY),
        alignment=TA_LEFT,
    )
    style_warning_body = ParagraphStyle(
        "FrutonWarnBody",
        fontName="Helvetica",
        fontSize=8,
        leading=10,
        textColor=colors.HexColor(_IDIS_TEXT),
        alignment=TA_LEFT,
    )

    def _draw_logo(canvas: Any, logo_path: Path | None, x: float, y: float, max_w: float, max_h: float) -> None:
        if logo_path is None:
            return
        try:
            reader = ImageReader(str(logo_path))
            img_w, img_h = reader.getSize()
            scale = min(max_w / img_w, max_h / img_h)
            draw_w = img_w * scale
            draw_h = img_h * scale
            canvas.drawImage(
                reader,
                x,
                y + (max_h - draw_h) / 2,
                width=draw_w,
                height=draw_h,
                preserveAspectRatio=True,
                mask="auto",
            )
        except Exception:
            return

    def _draw_page_chrome(canvas: Any, width: float, height: float, subtitle: str) -> None:
        canvas.saveState()
        header_y = height - 1.10 * cm
        logo_h = 0.62 * cm
        _draw_logo(canvas, idis_logo_path, margin_x, header_y, 2.4 * cm, logo_h)
        _draw_logo(canvas, fruton_logo_path, width - margin_x - 2.4 * cm, header_y, 2.4 * cm, logo_h)

        canvas.setFillColor(colors.HexColor(_IDIS_NAVY))
        canvas.setFont("Helvetica-Bold", 8)
        canvas.drawCentredString(width / 2, height - 0.58 * cm, "FRUTON pipeline summary audit")
        canvas.setFillColor(colors.HexColor(_IDIS_MUTED))
        canvas.setFont("Helvetica", 6.5)
        canvas.drawCentredString(width / 2, height - 0.82 * cm, subtitle)

        line_y = height - 1.33 * cm
        line_w = width - 2 * margin_x
        canvas.setStrokeColor(colors.HexColor(_FRUTON_RED))
        canvas.setLineWidth(0.75)
        canvas.line(margin_x, line_y, margin_x + line_w * 0.32, line_y)
        canvas.setStrokeColor(colors.HexColor(_IDIS_GOLD))
        canvas.line(margin_x + line_w * 0.32, line_y, margin_x + line_w * 0.64, line_y)
        canvas.setStrokeColor(colors.HexColor(_IDIS_NAVY))
        canvas.line(margin_x + line_w * 0.64, line_y, margin_x + line_w, line_y)

        canvas.setFillColor(colors.HexColor(_IDIS_MUTED))
        canvas.setFont("Helvetica", 6.5)
        canvas.drawString(margin_x, 0.48 * cm, f"FRUTON - {_FRUTON_EXPANSION}")
        canvas.drawRightString(width - margin_x, 0.48 * cm, f"page {canvas.getPageNumber()}")
        canvas.restoreState()

    def _footer_portrait(canvas: Any, doc: Any) -> None:
        _draw_page_chrome(canvas, page_w, page_h, f"run summary - {ts}")

    def _footer_landscape(canvas: Any, doc: Any) -> None:
        _draw_page_chrome(canvas, land_w, land_h, f"feature chart - {ts}")

    def _section_heading(text: str, width: float) -> Table:
        tbl = Table([["", Paragraph(text, style_section)]], colWidths=[0.10 * cm, width - 0.10 * cm])
        tbl.setStyle(
            TableStyle(
                [
                    ("BACKGROUND", (0, 0), (0, 0), colors.HexColor(_FRUTON_RED)),
                    ("VALIGN", (0, 0), (-1, -1), "MIDDLE"),
                    ("LEFTPADDING", (0, 0), (-1, -1), 0),
                    ("RIGHTPADDING", (0, 0), (-1, -1), 0),
                    ("TOPPADDING", (0, 0), (-1, -1), 0),
                    ("BOTTOMPADDING", (0, 0), (-1, -1), 2),
                ]
            )
        )
        return tbl

    story: list[Any] = []
    story.append(Paragraph("FRUTON - IDIS visual report", style_kicker))
    story.append(Paragraph("Pipeline preparation summary audit", style_title))
    story.append(
        Paragraph(
            f"Generated: {ts} &nbsp;&middot;&nbsp; proteins: {stats['total']} &nbsp;&middot;&nbsp; "
            f"chart pages: {'landscape chunks' if use_landscape_charts else 'portrait summary'}",
            style_subtitle,
        )
    )
    story.append(Spacer(1, 0.28 * cm))

    dashboard_rows = [
        ["PROTEINS", str(stats["total"]), "PREPARED", str(stats["prepared"])],
        ["FAILED", str(stats["failed"]), "GAPS", str(stats["with_gaps"])],
        ["METALS", str(stats["with_metals"]), "NONSTANDARD", str(stats["with_nonstd"])],
        ["COFACTORS", str(stats["with_cofactors"]), "GAUSSIAN PENDING", str(stats["gaussian_pending"])],
    ]
    dashboard_tbl = Table(
        dashboard_rows,
        colWidths=[usable_w * 0.22, usable_w * 0.28, usable_w * 0.22, usable_w * 0.28],
    )
    dashboard_style: list[Any] = [
        ("FONTNAME", (0, 0), (-1, -1), "Helvetica-Bold"),
        ("FONTSIZE", (0, 0), (-1, -1), 8),
        ("TEXTCOLOR", (0, 0), (-1, -1), colors.HexColor(_IDIS_TEXT)),
        ("BACKGROUND", (0, 0), (0, -1), colors.HexColor(_IDIS_PANEL_2)),
        ("BACKGROUND", (2, 0), (2, -1), colors.HexColor(_IDIS_PANEL_2)),
        ("TEXTCOLOR", (0, 0), (0, -1), colors.HexColor(_IDIS_MUTED)),
        ("TEXTCOLOR", (2, 0), (2, -1), colors.HexColor(_IDIS_MUTED)),
        ("GRID", (0, 0), (-1, -1), 0.35, colors.HexColor(_IDIS_LINE)),
        ("VALIGN", (0, 0), (-1, -1), "MIDDLE"),
        ("TOPPADDING", (0, 0), (-1, -1), 6),
        ("BOTTOMPADDING", (0, 0), (-1, -1), 6),
        ("LEFTPADDING", (0, 0), (-1, -1), 7),
        ("RIGHTPADDING", (0, 0), (-1, -1), 7),
    ]
    value_cell_map = {
        (1, 0): ("total", stats["total"]),
        (3, 0): ("prepared", stats["prepared"]),
        (1, 1): ("failed", stats["failed"]),
        (3, 1): ("with_gaps", stats["with_gaps"]),
        (1, 2): ("with_metals", stats["with_metals"]),
        (3, 2): ("with_nonstd", stats["with_nonstd"]),
        (1, 3): ("with_cofactors", stats["with_cofactors"]),
        (3, 3): ("gaussian_pending", stats["gaussian_pending"]),
    }
    for (col, row), (metric, value) in value_cell_map.items():
        dashboard_style.append(
            ("BACKGROUND", (col, row), (col, row), colors.HexColor(_status_background(metric, value)))
        )
    dashboard_tbl.setStyle(TableStyle(dashboard_style))
    story.append(dashboard_tbl)
    story.append(Spacer(1, 0.42 * cm))

    story.append(_section_heading("Run interpretation", usable_w))
    story.append(
        Paragraph(
            "This summary report is the dataset-level audit surface. It does not replace "
            "the per-protein reports; it gives a compact cross-protein view of prepared "
            "status, gaps, transition metals, cofactors, non-standard residues, and pending "
            "Gaussian work before the detailed evidence is inspected.",
            style_normal,
        )
    )
    story.append(Spacer(1, 0.38 * cm))

    if gaussian_pending_ids:
        story.append(_section_heading("Gaussian continuation required", usable_w))
        warning_rows = [
            [Paragraph("Pending proteins", style_warning_title), Paragraph(", ".join(gaussian_pending_ids), style_warning_body)],
            [Paragraph("Submit", style_warning_title), Paragraph("Run submit_gaussian.sh (RESP) or MCPB_submit.sh (metals) for each affected protein on the HPC cluster.", style_warning_body)],
            [Paragraph("Resume", style_warning_title), Paragraph("After Gaussian finishes, run run_after_gaussian.sh in each nonstandard_params/<label>/resp/ directory and resume with fruton.py --from-step 15.", style_warning_body)],
        ]
        warning_tbl = Table(warning_rows, colWidths=[usable_w * 0.24, usable_w * 0.76])
        warning_tbl.setStyle(
            TableStyle(
                [
                    ("BACKGROUND", (0, 0), (-1, -1), colors.HexColor(_FRUTON_RED_SOFT)),
                    ("BOX", (0, 0), (-1, -1), 0.45, colors.HexColor(_FRUTON_RED)),
                    ("INNERGRID", (0, 0), (-1, -1), 0.25, colors.HexColor(_IDIS_LINE)),
                    ("VALIGN", (0, 0), (-1, -1), "TOP"),
                    ("TOPPADDING", (0, 0), (-1, -1), 5),
                    ("BOTTOMPADDING", (0, 0), (-1, -1), 5),
                    ("LEFTPADDING", (0, 0), (-1, -1), 7),
                    ("RIGHTPADDING", (0, 0), (-1, -1), 7),
                ]
            )
        )
        story.append(warning_tbl)
        story.append(Spacer(1, 0.35 * cm))

    tmp_paths: list[str] = []
    try:
        chart_pngs = _render_chart_pngs(
            chart_data,
            landscape=use_landscape_charts,
            chunk_size=chart_chunk_size,
        )
        for png_bytes in chart_pngs:
            tmp = tempfile.NamedTemporaryFile(suffix=".png", delete=False)
            tmp.write(png_bytes)
            tmp.close()
            tmp_paths.append(tmp.name)

        if not use_landscape_charts:
            story.append(_section_heading("Per-protein feature chart", usable_w))
            chart_img = Image(tmp_paths[0])
            aspect = chart_img.imageWidth / chart_img.imageHeight
            chart_img.drawWidth = usable_w * 0.96
            chart_img.drawHeight = chart_img.drawWidth / aspect
            chart_img.hAlign = "CENTER"
            story.append(chart_img)
        else:
            for tmp_path in tmp_paths:
                story.append(NextPageTemplate("landscape"))
                story.append(PageBreak())
                chart_img = Image(tmp_path)
                aspect = chart_img.imageWidth / chart_img.imageHeight
                chart_img.drawWidth = land_usable_w * 0.96
                chart_img.drawHeight = chart_img.drawWidth / aspect
                chart_img.hAlign = "CENTER"
                story.append(chart_img)
    except Exception as exc:
        story.append(Paragraph(f"Chart could not be rendered: {exc}", style_normal))

    portrait_frame = Frame(
        margin_x,
        bottom_margin,
        page_w - 2 * margin_x,
        page_h - top_margin - bottom_margin,
        id="portrait_frame",
    )
    portrait_template = PageTemplate(
        id="portrait",
        frames=[portrait_frame],
        onPage=_footer_portrait,
        pagesize=A4,
    )

    land_frame = Frame(
        margin_x,
        bottom_margin,
        land_w - 2 * margin_x,
        land_h - top_margin - bottom_margin,
        id="landscape_frame",
    )
    landscape_template = PageTemplate(
        id="landscape",
        frames=[land_frame],
        onPage=_footer_landscape,
        pagesize=rl_landscape(A4),
    )

    doc = BaseDocTemplate(
        str(output_path),
        pagesize=A4,
        leftMargin=margin_x,
        rightMargin=margin_x,
        topMargin=top_margin,
        bottomMargin=bottom_margin,
    )
    doc.addPageTemplates([portrait_template, landscape_template])
    doc.build(story)

    for tmp_path in tmp_paths:
        try:
            os.unlink(tmp_path)
        except Exception:
            pass

    return {"status": "success", "report_path": str(output_path)}
