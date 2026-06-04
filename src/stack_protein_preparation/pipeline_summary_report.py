"""Pipeline-level PDF summary report for the FRUTON protein preparation pipeline.

Generates one PDF covering all proteins in a pipeline run, with a dashboard
of preparation status counts and a per-protein feature chart.

Entry point: generate_pipeline_summary_report()
"""
from __future__ import annotations

import os
import tempfile
from datetime import datetime
from pathlib import Path
from typing import Any

from ._report_palette import (
    _IDIS_NAVY,
    _IDIS_TEXT,
    _IDIS_MUTED,
    _IDIS_LINE,
    _IDIS_PANEL,
    _IDIS_PANEL_2,
    _IDIS_GOLD,
    _FRUTON_RED,
    _FRUTON_RED_SOFT,
    _FRUTON_GREEN_SOFT,
    _FRUTON_EXPANSION,
)
from ._report_assets import find_idis_logo, find_fruton_logo
from ._report_chart import (
    _extract_chart_data,
    _render_chart_pngs,
    _compute_stats,
    _status_background,
)


def _timestamp() -> str:
    return datetime.now().strftime("%Y-%m-%d %H:%M:%S")


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

    idis_logo_path = find_idis_logo(protein_data_dir)
    fruton_logo_path = find_fruton_logo(protein_data_dir)

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
