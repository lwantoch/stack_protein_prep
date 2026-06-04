from __future__ import annotations

import io
import os
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


_IDIS_LOGO_NAMES = ("logo_idis.png", "logo_IDIS_2020-1.png", "idis_logo.png")
_FRUTON_LOGO_NAMES = ("logo.png", "fruton_logo.png", "FRUTON_logo.png", "logo(1).png")


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
            BaseDocTemplate,
            Frame,
            PageTemplate,
            NextPageTemplate,
            PageBreak,
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

    # Thresholds: portrait fits up to 60 proteins in the chart; landscape fits up to 90.
    _PORTRAIT_CHART_MAX = 60
    _PORTRAIT_CHUNK = 55
    _LANDSCAPE_CHUNK = 90

    ts = run_timestamp or _timestamp()
    margin = 1.5 * cm
    footer_h = 1.2 * cm

    port_w, port_h = A4
    land_w, land_h = landscape(A4)

    idis_logo_path = _find_logo("IDIS_LOGO_PATH", _IDIS_LOGO_NAMES, protein_data_dir)
    fruton_logo_path = _find_logo("FRUTON_LOGO_PATH", _FRUTON_LOGO_NAMES, protein_data_dir)

    chart_data = _extract_chart_data(pipeline_record_list)
    stats = _compute_stats(pipeline_record_list, chart_data, gaussian_pending_ids)

    n_bars = len(chart_data)
    use_landscape_chart = n_bars > _PORTRAIT_CHART_MAX
    chunk_size = _LANDSCAPE_CHUNK if use_landscape_chart else _PORTRAIT_CHUNK
    chart_chunks: list[list[dict]] = (
        [chart_data[i: i + chunk_size] for i in range(0, n_bars, chunk_size)]
        if n_bars > 0
        else [chart_data]
    )

    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    # Page templates
    port_usable_w = port_w - 2 * margin
    port_usable_h = port_h - 2 * margin - footer_h
    land_usable_w = land_w - 2 * margin
    land_usable_h = land_h - 2 * margin - footer_h

    def _make_footer(page_width: float) -> Any:
        def _footer(canvas: Any, doc: Any) -> None:
            canvas.saveState()
            canvas.setFont("Helvetica", 7)
            canvas.setFillColor(colors.HexColor(_IDIS_MUTED))
            canvas.drawCentredString(page_width / 2, 0.6 * cm, f"Generated by FRUTON  ·  {ts}")
            canvas.restoreState()
        return _footer

    port_frame = Frame(margin, margin + footer_h, port_usable_w, port_usable_h, id="portrait_frame")
    port_tmpl = PageTemplate(id="portrait", frames=[port_frame], pagesize=A4,
                              onPage=_make_footer(port_w))
    land_frame = Frame(margin, margin + footer_h, land_usable_w, land_usable_h, id="landscape_frame")
    land_tmpl = PageTemplate(id="landscape", frames=[land_frame], pagesize=landscape(A4),
                              onPage=_make_footer(land_w))

    doc = BaseDocTemplate(
        str(output_path),
        pageTemplates=[port_tmpl, land_tmpl],
        leftMargin=margin,
        rightMargin=margin,
        topMargin=margin,
        bottomMargin=margin + footer_h,
    )

    # Paragraph styles
    style_normal = ParagraphStyle(
        "FrutonNormal", fontName="Helvetica", fontSize=9,
        textColor=colors.HexColor(_IDIS_TEXT),
    )
    style_header_center = ParagraphStyle(
        "FrutonHeaderCenter", fontName="Helvetica-Bold", fontSize=13,
        textColor=colors.white, alignment=TA_CENTER,
    )
    style_header_left = ParagraphStyle(
        "FrutonHeaderLeft", fontName="Helvetica-Bold", fontSize=10,
        textColor=colors.white, alignment=TA_LEFT,
    )
    style_header_right = ParagraphStyle(
        "FrutonHeaderRight", fontName="Helvetica", fontSize=8,
        textColor=colors.HexColor(_IDIS_LINE), alignment=TA_RIGHT,
    )
    style_warning_title = ParagraphStyle(
        "FrutonWarnTitle", fontName="Helvetica-Bold", fontSize=10, textColor=colors.white,
    )
    style_warning_body = ParagraphStyle(
        "FrutonWarnBody", fontName="Helvetica", fontSize=8,
        textColor=colors.HexColor(_IDIS_TEXT), leftIndent=4,
    )

    # Helper: scale logo to target height
    _LOGO_H = 32

    def _header_logo(logo_path: Path | None, align: str = "LEFT") -> Any:
        if logo_path is None:
            return None
        try:
            img = Image(str(logo_path))
            scale = _LOGO_H / img.imageHeight
            img.drawWidth = img.imageWidth * scale
            img.drawHeight = _LOGO_H
            img.hAlign = align
            return img
        except Exception:
            return None

    def _build_header(usable_w: float) -> list[Any]:
        idis_img = _header_logo(idis_logo_path, "LEFT")
        fruton_img = _header_logo(fruton_logo_path, "RIGHT")
        left_cell: Any = idis_img if idis_img is not None else Paragraph("FRUTON", style_header_left)
        right_cell: Any = fruton_img if fruton_img is not None else Paragraph(ts, style_header_right)
        col_w = usable_w / 3
        tbl = Table(
            [[left_cell, Paragraph("Pipeline Summary Report", style_header_center), right_cell]],
            colWidths=[col_w, col_w, col_w],
        )
        tbl.setStyle(TableStyle([
            ("BACKGROUND", (0, 0), (-1, -1), colors.HexColor(_IDIS_NAVY)),
            ("VALIGN", (0, 0), (-1, -1), "MIDDLE"),
            ("TOPPADDING", (0, 0), (-1, -1), 6),
            ("BOTTOMPADDING", (0, 0), (-1, -1), 6),
            ("LEFTPADDING", (0, 0), (0, -1), 10),
            ("RIGHTPADDING", (-1, 0), (-1, -1), 10),
            ("LINEBELOW", (0, 0), (-1, -1), 2, colors.HexColor(_FRUTON_RED)),
        ]))
        elems: list[Any] = [tbl]
        if fruton_img is not None:
            elems.append(Paragraph(
                f'<font size="7" color="{_IDIS_MUTED}">{ts}</font>',
                ParagraphStyle("FrutonTS", fontName="Helvetica", fontSize=7,
                               textColor=colors.HexColor(_IDIS_MUTED), alignment=TA_RIGHT),
            ))
        elems.append(Spacer(1, 0.3 * cm))
        return elems

    def _embed_chart(chunk: list[dict], usable_w: float) -> Any:
        try:
            png = _render_chart_png(chunk)
            with tempfile.NamedTemporaryFile(suffix=".png", delete=False) as tmp:
                tmp.write(png)
                tmp_path = tmp.name
            img = Image(tmp_path)
            aspect = img.imageWidth / img.imageHeight
            img.drawWidth = usable_w * 0.99
            img.drawHeight = img.drawWidth / aspect
            img.hAlign = "CENTER"
            return img
        except Exception:
            return Paragraph("Chart could not be rendered.", style_normal)

    story: list[Any] = []

    # --- Portrait page: header + stats (+ gaussian warning if any) ---
    story.append(NextPageTemplate("portrait"))
    story += _build_header(port_usable_w)

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
    col_l = port_usable_w * 0.65
    col_r = port_usable_w * 0.35
    stats_table = Table(stats_rows, colWidths=[col_l, col_r])
    stats_style = [
        ("BACKGROUND", (0, 0), (-1, 0), colors.HexColor(_IDIS_NAVY)),
        ("TEXTCOLOR", (0, 0), (-1, 0), colors.white),
        ("FONTNAME", (0, 0), (-1, 0), "Helvetica-Bold"),
        ("FONTSIZE", (0, 0), (-1, 0), 9),
        ("ROWBACKGROUNDS", (0, 1), (-1, -1), [colors.white, colors.HexColor(_IDIS_PANEL)]),
        ("FONTNAME", (0, 1), (-1, -1), "Helvetica"),
        ("FONTSIZE", (0, 1), (-1, -1), 8),
        ("TEXTCOLOR", (0, 1), (-1, -1), colors.HexColor(_IDIS_TEXT)),
        ("ALIGN", (1, 0), (1, -1), "RIGHT"),
        ("GRID", (0, 0), (-1, -1), 0.5, colors.HexColor(_IDIS_LINE)),
        ("TOPPADDING", (0, 0), (-1, -1), 4),
        ("BOTTOMPADDING", (0, 0), (-1, -1), 4),
        ("LEFTPADDING", (0, 0), (-1, -1), 6),
        ("RIGHTPADDING", (0, 0), (-1, -1), 6),
    ]
    if stats["gaussian_pending"] > 0:
        stats_style.append(("BACKGROUND", (0, -1), (-1, -1), colors.HexColor(_FRUTON_RED_SOFT)))
    stats_table.setStyle(TableStyle(stats_style))
    story.append(stats_table)

    if gaussian_pending_ids:
        story.append(Spacer(1, 0.5 * cm))
        warn_title_data = [[Paragraph(
            "Gaussian HPC runs required before pipeline is complete", style_warning_title
        )]]
        warn_title_tbl = Table(warn_title_data, colWidths=[port_usable_w])
        warn_title_tbl.setStyle(TableStyle([
            ("BACKGROUND", (0, 0), (-1, -1), colors.HexColor(_FRUTON_RED)),
            ("TOPPADDING", (0, 0), (-1, -1), 6),
            ("BOTTOMPADDING", (0, 0), (-1, -1), 6),
            ("LEFTPADDING", (0, 0), (-1, -1), 8),
        ]))
        story.append(warn_title_tbl)
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
        warn_body_tbl = Table(warn_body_data, colWidths=[port_usable_w])
        warn_body_tbl.setStyle(TableStyle([
            ("BACKGROUND", (0, 0), (-1, -1), colors.HexColor(_FRUTON_RED_SOFT)),
            ("TOPPADDING", (0, 0), (-1, -1), 2),
            ("BOTTOMPADDING", (0, 0), (-1, -1), 2),
            ("LEFTPADDING", (0, 0), (-1, -1), 8),
        ]))
        story.append(warn_body_tbl)

    # --- Chart pages ---
    chart_template = "landscape" if use_landscape_chart else "portrait"
    chart_usable_w = land_usable_w if use_landscape_chart else port_usable_w
    story.append(NextPageTemplate(chart_template))
    story.append(PageBreak())

    for i, chunk in enumerate(chart_chunks):
        story.append(_embed_chart(chunk, chart_usable_w))
        if i < len(chart_chunks) - 1:
            story.append(PageBreak())

    doc.build(story)

    return {"status": "success", "report_path": str(output_path)}
