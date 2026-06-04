from __future__ import annotations

import io
import tempfile
from datetime import datetime
from pathlib import Path
from typing import Any

_IDIS_NAVY = "#28325A"
_IDIS_TEXT = "#252A3F"
_IDIS_MUTED = "#697086"
_IDIS_LINE = "#D8DCE7"
_IDIS_PANEL = "#F5F7FB"
_IDIS_GOLD = "#8A730F"
_FRUTON_RED = "#E60014"
_FRUTON_RED_SOFT = "#FDECEE"
_FRUTON_GREEN = "#50961E"

_TRANSITION_METAL_ELEMENTS = frozenset(
    {"ZN", "CU", "FE", "MN", "CO", "NI", "MO", "CD", "HG", "PD", "PT"}
)


def _timestamp() -> str:
    return datetime.now().strftime("%Y-%m-%d %H:%M:%S")


def _hex_to_rgb_float(hex_color: str) -> tuple[float, float, float]:
    h = hex_color.lstrip("#")
    r, g, b = int(h[0:2], 16), int(h[2:4], 16), int(h[4:6], 16)
    return r / 255.0, g / 255.0, b / 255.0


def _extract_chart_data(pipeline_record_list: list[dict]) -> list[dict]:
    rows: list[dict] = []
    for rec in pipeline_record_list:
        pdb_id = str(rec.get("pdb_id", "")).strip()
        if not pdb_id:
            continue

        has_insertions = str(rec.get("insertion_codes.done", "")).strip().lower() == "success"

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

        rows.append(
            {
                "pdb_id": pdb_id,
                "has_insertions": has_insertions,
                "n_gaps": n_gaps,
                "total_gap_length": total_gap_length,
                "has_cofactors": has_cofactors,
                "has_nonstd": has_nonstd,
                "has_transition_metal": has_transition_metal,
            }
        )
    return rows


def _render_chart_png(chart_data: list[dict]) -> bytes:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import matplotlib.patches as mpatches
    import numpy as np

    if not chart_data:
        fig, ax = plt.subplots(figsize=(8, 4))
        ax.text(0.5, 0.5, "No data", ha="center", va="center", transform=ax.transAxes)
        buf = io.BytesIO()
        fig.savefig(buf, format="png", dpi=150, bbox_inches="tight")
        plt.close(fig)
        return buf.getvalue()

    pdb_ids = [d["pdb_id"] for d in chart_data]
    n = len(pdb_ids)
    x = np.arange(n)

    has_insertions = np.array([1.0 if d["has_insertions"] else 0.0 for d in chart_data])
    n_gaps = np.array([d["n_gaps"] for d in chart_data], dtype=float)
    total_gap_length = np.array([d["total_gap_length"] for d in chart_data], dtype=float)
    has_cofactors = np.array([1.0 if d["has_cofactors"] else 0.0 for d in chart_data])
    has_nonstd = np.array([1.0 if d["has_nonstd"] else 0.0 for d in chart_data])
    has_transition_metal = np.array([1.0 if d["has_transition_metal"] else 0.0 for d in chart_data])

    fig, axes = plt.subplots(
        5, 1,
        figsize=(max(24, n * 0.45), 11),
        sharex=True,
        gridspec_kw={"height_ratios": [1, 2.5, 1, 1, 1]},
    )
    fig.subplots_adjust(hspace=0.08)

    def _style_binary_ax(ax: Any, values: np.ndarray, color: str, ylabel: str) -> None:
        ax.bar(x, values, color=color, width=0.7, zorder=2)
        ax.set_ylim(-0.05, 1.3)
        ax.set_yticks([0, 1])
        ax.set_yticklabels(["no", "yes"], fontsize=7)
        ax.set_ylabel(ylabel, fontsize=8, labelpad=4)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.yaxis.grid(True, color="#E0E0E0", linewidth=0.5, zorder=1)
        ax.set_axisbelow(True)

    _style_binary_ax(axes[0], has_insertions, "#E07B28", "insertions")

    # Stacked gap bar
    gap_ax = axes[1]
    gap_ax.bar(x, n_gaps, color="#2C5BAA", width=0.7, label="n_gaps", zorder=2)
    gap_ax.bar(x, total_gap_length, bottom=n_gaps, color="#7BA8D6", width=0.7, label="total_gap_length", zorder=2)
    gap_ax.set_ylabel("gaps", fontsize=8, labelpad=4)
    gap_ax.spines["top"].set_visible(False)
    gap_ax.spines["right"].set_visible(False)
    gap_ax.yaxis.grid(True, color="#E0E0E0", linewidth=0.5, zorder=1)
    gap_ax.set_axisbelow(True)
    legend_patches = [
        mpatches.Patch(color="#2C5BAA", label="n_gaps"),
        mpatches.Patch(color="#7BA8D6", label="total_gap_length"),
    ]
    gap_ax.legend(handles=legend_patches, fontsize=7, loc="upper right", framealpha=0.7)

    _style_binary_ax(axes[2], has_cofactors, "#50961E", "cofactors")
    _style_binary_ax(axes[3], has_nonstd, "#7B3FA0", "non-std res")
    _style_binary_ax(axes[4], has_transition_metal, "#E60014", "trans. metals")

    axes[4].set_xticks(x)
    axes[4].set_xticklabels(pdb_ids, rotation=45, ha="right", fontsize=max(5, min(8, 200 // n)))

    fig.patch.set_facecolor("white")

    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=150, bbox_inches="tight")
    plt.close(fig)
    return buf.getvalue()


def _compute_stats(
    pipeline_record_list: list[dict],
    chart_data: list[dict],
    gaussian_pending_ids: list[str] | None,
) -> dict:
    total = len(pipeline_record_list)
    prepared = sum(
        1 for r in pipeline_record_list
        if str(r.get("prepared_structure.status", "") or "").strip().lower() == "success"
    )
    failed = sum(
        1 for r in pipeline_record_list
        if str(r.get("prepared_structure.status", "") or "").strip().lower() in {"failed", "error"}
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


def generate_pipeline_summary_report(
    pipeline_record_list: list[dict],
    output_path: Path,
    *,
    gaussian_pending_ids: list[str] | None = None,
    protein_data_dir: Path | None = None,
    run_timestamp: str | None = None,
) -> dict:
    """Return {"status": "success", "report_path": str} or {"status": "failed", "message": str}"""
    try:
        from reportlab.lib import colors
        from reportlab.lib.pagesizes import A4, landscape
        from reportlab.lib.units import cm
        from reportlab.platypus import (
            SimpleDocTemplate,
            Spacer,
            Table,
            TableStyle,
            Image,
            Paragraph,
        )
        from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
        from reportlab.lib.enums import TA_CENTER, TA_LEFT, TA_RIGHT
    except ImportError as exc:
        return {"status": "failed", "message": f"reportlab not available: {exc}"}

    ts = run_timestamp or _timestamp()
    page_w, page_h = landscape(A4)
    margin = 1.5 * cm

    chart_data = _extract_chart_data(pipeline_record_list)
    stats = _compute_stats(pipeline_record_list, chart_data, gaussian_pending_ids)

    try:
        chart_png = _render_chart_png(chart_data)
    except Exception as exc:
        chart_png = None

    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    doc = SimpleDocTemplate(
        str(output_path),
        pagesize=landscape(A4),
        leftMargin=margin,
        rightMargin=margin,
        topMargin=margin,
        bottomMargin=margin,
    )

    usable_w = page_w - 2 * margin

    styles = getSampleStyleSheet()
    style_normal = ParagraphStyle(
        "FrutonNormal",
        fontName="Helvetica",
        fontSize=9,
        textColor=colors.HexColor(_IDIS_TEXT),
    )
    style_header_center = ParagraphStyle(
        "FrutonHeaderCenter",
        fontName="Helvetica-Bold",
        fontSize=13,
        textColor=colors.white,
        alignment=TA_CENTER,
    )
    style_header_left = ParagraphStyle(
        "FrutonHeaderLeft",
        fontName="Helvetica-Bold",
        fontSize=10,
        textColor=colors.white,
        alignment=TA_LEFT,
    )
    style_header_right = ParagraphStyle(
        "FrutonHeaderRight",
        fontName="Helvetica",
        fontSize=8,
        textColor=colors.HexColor(_IDIS_LINE),
        alignment=TA_RIGHT,
    )
    style_warning_title = ParagraphStyle(
        "FrutonWarnTitle",
        fontName="Helvetica-Bold",
        fontSize=10,
        textColor=colors.white,
    )
    style_warning_body = ParagraphStyle(
        "FrutonWarnBody",
        fontName="Helvetica",
        fontSize=8,
        textColor=colors.HexColor(_IDIS_TEXT),
        leftIndent=4,
    )

    story: list[Any] = []

    # Header row
    header_col_w = usable_w / 3
    header_data = [[
        Paragraph("FRUTON", style_header_left),
        Paragraph("Pipeline Summary Report", style_header_center),
        Paragraph(ts, style_header_right),
    ]]
    header_table = Table(header_data, colWidths=[header_col_w, header_col_w, header_col_w])
    header_table.setStyle(TableStyle([
        ("BACKGROUND", (0, 0), (-1, -1), colors.HexColor(_IDIS_NAVY)),
        ("VALIGN", (0, 0), (-1, -1), "MIDDLE"),
        ("TOPPADDING", (0, 0), (-1, -1), 8),
        ("BOTTOMPADDING", (0, 0), (-1, -1), 8),
        ("LEFTPADDING", (0, 0), (0, -1), 10),
        ("RIGHTPADDING", (-1, 0), (-1, -1), 10),
        ("LINEBELOW", (0, 0), (-1, -1), 2, colors.HexColor(_FRUTON_RED)),
    ]))
    story.append(header_table)
    story.append(Spacer(1, 0.4 * cm))

    # Stats table (2-column layout)
    stats_rows = [
        ["Metric", "Value"],
        ["Proteins processed", str(stats["total"])],
        ["Proteins prepared (success)", str(stats["prepared"])],
        ["Proteins failed", str(stats["failed"])],
        ["Proteins with gaps", str(stats["with_gaps"])],
        ["Proteins with transition metals", str(stats["with_metals"])],
        ["Proteins with non-standard residues", str(stats["with_nonstd"])],
        ["Proteins with cofactors", str(stats["with_cofactors"])],
        ["Gaussian pending (HPC required)", str(stats["gaussian_pending"])],
    ]
    col_l = usable_w * 0.55
    col_r = usable_w * 0.45
    stats_table = Table(stats_rows, colWidths=[col_l, col_r])
    stats_style = [
        ("BACKGROUND", (0, 0), (-1, 0), colors.HexColor(_IDIS_NAVY)),
        ("TEXTCOLOR", (0, 0), (-1, 0), colors.white),
        ("FONTNAME", (0, 0), (-1, 0), "Helvetica-Bold"),
        ("FONTSIZE", (0, 0), (-1, 0), 9),
        ("BACKGROUND", (0, 1), (-1, -1), colors.HexColor(_IDIS_PANEL)),
        ("ROWBACKGROUNDS", (0, 1), (-1, -1), [colors.white, colors.HexColor(_IDIS_PANEL)]),
        ("FONTNAME", (0, 1), (-1, -1), "Helvetica"),
        ("FONTSIZE", (0, 1), (-1, -1), 8),
        ("TEXTCOLOR", (0, 1), (-1, -1), colors.HexColor(_IDIS_TEXT)),
        ("ALIGN", (1, 0), (1, -1), "RIGHT"),
        ("FONTNAME", (0, 1), (0, -1), "Helvetica"),
        ("GRID", (0, 0), (-1, -1), 0.5, colors.HexColor(_IDIS_LINE)),
        ("TOPPADDING", (0, 0), (-1, -1), 4),
        ("BOTTOMPADDING", (0, 0), (-1, -1), 4),
        ("LEFTPADDING", (0, 0), (-1, -1), 6),
        ("RIGHTPADDING", (0, 0), (-1, -1), 6),
    ]
    if stats["gaussian_pending"] > 0:
        stats_style.append(
            ("BACKGROUND", (0, -1), (-1, -1), colors.HexColor(_FRUTON_RED_SOFT))
        )
    stats_table.setStyle(TableStyle(stats_style))
    story.append(stats_table)
    story.append(Spacer(1, 0.5 * cm))

    # Chart
    if chart_png is not None:
        try:
            with tempfile.NamedTemporaryFile(suffix=".png", delete=False) as tmp:
                tmp.write(chart_png)
                tmp_path = tmp.name
            chart_img = Image(tmp_path)
            chart_img_w = usable_w * 0.98
            aspect = chart_img.imageWidth / chart_img.imageHeight
            chart_img.drawWidth = chart_img_w
            chart_img.drawHeight = chart_img_w / aspect
            story.append(chart_img)
        except Exception:
            story.append(Paragraph("Chart could not be embedded.", style_normal))
    else:
        story.append(Paragraph("Chart could not be rendered (matplotlib unavailable).", style_normal))

    # Gaussian pending section
    if gaussian_pending_ids:
        story.append(Spacer(1, 0.5 * cm))
        warn_title_data = [[Paragraph("Gaussian HPC runs required before pipeline is complete", style_warning_title)]]
        warn_title_table = Table(warn_title_data, colWidths=[usable_w])
        warn_title_table.setStyle(TableStyle([
            ("BACKGROUND", (0, 0), (-1, -1), colors.HexColor(_FRUTON_RED)),
            ("TOPPADDING", (0, 0), (-1, -1), 6),
            ("BOTTOMPADDING", (0, 0), (-1, -1), 6),
            ("LEFTPADDING", (0, 0), (-1, -1), 8),
        ]))
        story.append(warn_title_table)
        pending_lines = [
            "The following proteins have Gaussian input files prepared but not yet computed on HPC:",
            "  " + ", ".join(gaussian_pending_ids),
            "",
            "To complete parametrization:",
            "  1. Submit: run submit_gaussian.sh (RESP) or MCPB_submit.sh (metals) for each protein on the HPC cluster",
            "  2. After Gaussian jobs finish: run run_after_gaussian.sh in each nonstandard_params/<label>/resp/ directory",
            "  3. Resume the pipeline: fruton.py --from-step 15",
        ]
        warn_body_data = [[Paragraph(line if line else " ", style_warning_body)] for line in pending_lines]
        warn_body_table = Table(warn_body_data, colWidths=[usable_w])
        warn_body_table.setStyle(TableStyle([
            ("BACKGROUND", (0, 0), (-1, -1), colors.HexColor(_FRUTON_RED_SOFT)),
            ("TOPPADDING", (0, 0), (-1, -1), 2),
            ("BOTTOMPADDING", (0, 0), (-1, -1), 2),
            ("LEFTPADDING", (0, 0), (-1, -1), 8),
        ]))
        story.append(warn_body_table)

    def _page_footer(canvas: Any, doc: Any) -> None:
        canvas.saveState()
        canvas.setFont("Helvetica", 7)
        canvas.setFillColor(colors.HexColor(_IDIS_MUTED))
        footer_text = f"Generated by FRUTON  ·  {ts}"
        canvas.drawCentredString(page_w / 2, 0.6 * cm, footer_text)
        canvas.restoreState()

    doc.build(story, onFirstPage=_page_footer, onLaterPages=_page_footer)

    return {"status": "success", "report_path": str(output_path)}
