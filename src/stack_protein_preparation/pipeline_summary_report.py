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


def _parse_residue_count(range_str: str) -> int:
    """Parse a residue range string like '16-243' or 'A16-B50' to residue count."""
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
    """Per-protein grouped bar chart with twin y-axis (residue count background bars).

    Returns one PNG bytes object per page chunk.
    """
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
        return [buf.getvalue()]

    _C_INS  = "#E07B28"   # orange — insertions
    _C_GAPN = "#2C5BAA"   # blue   — gap count (bottom segment)
    _C_GAPR = "#7BA8D6"   # light blue — gap residues (stacked)
    _C_MET  = "#E60014"   # red    — transition metals
    _C_COF  = "#50961E"   # green  — cofactors
    _C_NST  = "#7B3FA0"   # purple — non-standard residues
    _C_BG   = _IDIS_NAVY  # navy   — residue count background bar

    bar_w = 0.13
    bar_gap = 0.015
    offsets = np.array([-2, -1, 0, 1, 2]) * (bar_w + bar_gap)
    bg_bar_w = 5 * bar_w + 4 * bar_gap + 0.06   # full group width + small padding

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

        ins  = np.array([int(d["has_insertions"])      for d in chunk], dtype=float)
        gn   = np.array([d["n_gaps"]                   for d in chunk], dtype=float)
        gr   = np.array([d["total_gap_length"]          for d in chunk], dtype=float)
        met  = np.array([int(d["has_transition_metal"]) for d in chunk], dtype=float)
        cof  = np.array([int(d["has_cofactors"])        for d in chunk], dtype=float)
        nst  = np.array([int(d["has_nonstd"])           for d in chunk], dtype=float)
        nres = np.array([d["n_residues"]                for d in chunk], dtype=float)

        fig_w = max(14, nc * 0.72) if landscape else max(9, nc * 1.3)
        fig_h = 5.5 if landscape else 5.2

        fig, ax = plt.subplots(figsize=(fig_w, fig_h))
        ax2 = ax.twinx()

        # Draw background (residue count) on ax2 behind the feature bars
        ax2.set_zorder(1)
        ax.set_zorder(2)
        ax.patch.set_visible(False)

        ax2.bar(x, nres, bg_bar_w, alpha=0.18, color=_C_BG, zorder=0)
        ax2.set_ylabel("Starting range residues", fontsize=8, color=_IDIS_MUTED, labelpad=6)
        ax2.tick_params(axis="y", labelcolor=_IDIS_MUTED, labelsize=7)
        ax2.spines["top"].set_visible(False)
        ax2.spines["right"].set_color(_IDIS_LINE)
        ax2.yaxis.grid(False)
        nres_max = max(float(nres.max()) if nc else 1.0, 1.0)
        ax2.set_ylim(0, nres_max * 1.55)

        # Feature bars on primary axis
        ax.bar(x + offsets[0], ins,  bar_w, color=_C_INS,  zorder=2)
        ax.bar(x + offsets[1], gn,   bar_w, color=_C_GAPN, zorder=2)
        ax.bar(x + offsets[1], gr,   bar_w, bottom=gn, color=_C_GAPR, zorder=2)
        ax.bar(x + offsets[2], met,  bar_w, color=_C_MET,  zorder=2)
        ax.bar(x + offsets[3], cof,  bar_w, color=_C_COF,  zorder=2)
        ax.bar(x + offsets[4], nst,  bar_w, color=_C_NST,  zorder=2)

        fs_tick = 7 if nc > 15 else 9
        ax.set_xticks(x)
        ax.set_xticklabels(pdb_ids, fontsize=fs_tick, rotation=45, ha="right")
        ax.set_ylabel("Count / presence (0 or 1)", fontsize=9, color=_IDIS_TEXT, labelpad=6)
        ax.tick_params(axis="y", labelcolor=_IDIS_TEXT, labelsize=8)

        feat_max = max(float((gn + gr).max()) if nc else 1.0, 1.0)
        ax.set_ylim(0, feat_max * 1.55)

        suffix = f" — page {chunk_idx + 1}" if len(chunks) > 1 else ""
        ax.set_title(
            f"Per-protein feature summary  ·  {n_total} protein(s){suffix}",
            fontsize=10, fontweight="bold", color=_IDIS_NAVY, pad=8,
        )

        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.spines["left"].set_color(_IDIS_LINE)
        ax.spines["bottom"].set_color(_IDIS_LINE)
        ax.yaxis.grid(True, color="#E0E0E0", linewidth=0.5, zorder=1)
        ax.set_axisbelow(True)

        legend_handles = [
            mpatches.Patch(color=_C_INS,  label="Insertions"),
            mpatches.Patch(color=_C_GAPN, label="Gaps (count)"),
            mpatches.Patch(color=_C_GAPR, label="Gap residues (stacked)"),
            mpatches.Patch(color=_C_MET,  label="Transition metals"),
            mpatches.Patch(color=_C_COF,  label="Cofactors"),
            mpatches.Patch(color=_C_NST,  label="Non-standard residues"),
            mpatches.Patch(color=_C_BG, alpha=0.35, label="Starting range residues (right axis)"),
        ]
        ax.legend(handles=legend_handles, fontsize=7, loc="upper left",
                  framealpha=0.88, ncol=2, borderpad=0.6)

        fig.patch.set_facecolor("white")
        fig.tight_layout()

        buf = io.BytesIO()
        fig.savefig(buf, format="png", dpi=150, bbox_inches="tight")
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
        from reportlab.lib.pagesizes import A4, landscape as rl_landscape
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
        from reportlab.lib.styles import ParagraphStyle
        from reportlab.lib.enums import TA_CENTER, TA_LEFT, TA_RIGHT
    except ImportError as exc:
        return {"status": "failed", "message": f"reportlab not available: {exc}"}

    ts = run_timestamp or _timestamp()
    margin = 1.5 * cm

    page_w, page_h = A4
    land_w, land_h = rl_landscape(A4)

    idis_logo_path = _find_logo("IDIS_LOGO_PATH", _IDIS_LOGO_NAMES, protein_data_dir)
    fruton_logo_path = _find_logo("FRUTON_LOGO_PATH", _FRUTON_LOGO_NAMES, protein_data_dir)

    chart_data = _extract_chart_data(pipeline_record_list)
    stats = _compute_stats(pipeline_record_list, chart_data, gaussian_pending_ids)

    n_proteins = len(chart_data)
    use_landscape_charts = n_proteins > 15
    _CHART_CHUNK = 20

    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    usable_w = page_w - 2 * margin
    land_usable_w = land_w - 2 * margin

    # --- Paragraph styles ---
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

    idis_img = _header_logo(idis_logo_path, "LEFT")
    fruton_img = _header_logo(fruton_logo_path, "RIGHT")
    left_cell: Any = idis_img if idis_img is not None else Paragraph("FRUTON", style_header_left)
    right_cell: Any = fruton_img if fruton_img is not None else Paragraph(ts, style_header_right)
    col_w = usable_w / 3
    header_tbl = Table(
        [[left_cell, Paragraph("Pipeline Summary Report", style_header_center), right_cell]],
        colWidths=[col_w, col_w, col_w],
    )
    header_tbl.setStyle(TableStyle([
        ("BACKGROUND", (0, 0), (-1, -1), colors.HexColor(_IDIS_NAVY)),
        ("VALIGN", (0, 0), (-1, -1), "MIDDLE"),
        ("TOPPADDING", (0, 0), (-1, -1), 6),
        ("BOTTOMPADDING", (0, 0), (-1, -1), 6),
        ("LEFTPADDING", (0, 0), (0, -1), 10),
        ("RIGHTPADDING", (-1, 0), (-1, -1), 10),
        ("LINEBELOW", (0, 0), (-1, -1), 2, colors.HexColor(_FRUTON_RED)),
    ]))

    story: list[Any] = [header_tbl]
    if fruton_img is not None:
        story.append(Paragraph(
            f'<font size="7" color="{_IDIS_MUTED}">{ts}</font>',
            ParagraphStyle("FrutonTS", fontName="Helvetica", fontSize=7,
                           textColor=colors.HexColor(_IDIS_MUTED), alignment=TA_RIGHT),
        ))
    story.append(Spacer(1, 0.35 * cm))

    # Stats table
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
    col_l = usable_w * 0.65
    col_r = usable_w * 0.35
    stats_tbl = Table(stats_rows, colWidths=[col_l, col_r])
    stats_style: list[Any] = [
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
    stats_tbl.setStyle(TableStyle(stats_style))
    story.append(stats_tbl)
    story.append(Spacer(1, 0.5 * cm))

    # Gaussian warning (on stats page before chart pages)
    if gaussian_pending_ids:
        warn_title_tbl = Table(
            [[Paragraph("Gaussian HPC runs required before pipeline is complete", style_warning_title)]],
            colWidths=[usable_w],
        )
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
        warn_body_tbl = Table(
            [[Paragraph(line if line else " ", style_warning_body)] for line in pending_lines],
            colWidths=[usable_w],
        )
        warn_body_tbl.setStyle(TableStyle([
            ("BACKGROUND", (0, 0), (-1, -1), colors.HexColor(_FRUTON_RED_SOFT)),
            ("TOPPADDING", (0, 0), (-1, -1), 2),
            ("BOTTOMPADDING", (0, 0), (-1, -1), 2),
            ("LEFTPADDING", (0, 0), (-1, -1), 8),
        ]))
        story.append(warn_body_tbl)
        story.append(Spacer(1, 0.3 * cm))

    # --- Chart pages ---
    tmp_paths: list[str] = []
    try:
        chart_pngs = _render_chart_pngs(
            chart_data, landscape=use_landscape_charts, chunk_size=_CHART_CHUNK
        )
        for png_bytes in chart_pngs:
            tmp = tempfile.NamedTemporaryFile(suffix=".png", delete=False)
            tmp.write(png_bytes)
            tmp.close()
            tmp_paths.append(tmp.name)

        if not use_landscape_charts:
            # Single portrait chart embedded on stats page
            chart_img = Image(tmp_paths[0])
            aspect = chart_img.imageWidth / chart_img.imageHeight
            chart_img.drawWidth = usable_w * 0.94
            chart_img.drawHeight = chart_img.drawWidth / aspect
            chart_img.hAlign = "CENTER"
            story.append(chart_img)
        else:
            # One landscape page per chart chunk
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

    # --- Page templates ---
    def _footer_portrait(canvas: Any, doc: Any) -> None:
        canvas.saveState()
        canvas.setFont("Helvetica", 7)
        canvas.setFillColor(colors.HexColor(_IDIS_MUTED))
        canvas.drawCentredString(page_w / 2, 0.5 * cm, f"Generated by FRUTON  ·  {ts}")
        canvas.restoreState()

    def _footer_landscape(canvas: Any, doc: Any) -> None:
        canvas.saveState()
        canvas.setFont("Helvetica", 7)
        canvas.setFillColor(colors.HexColor(_IDIS_MUTED))
        canvas.drawCentredString(land_w / 2, 0.5 * cm, f"Generated by FRUTON  ·  {ts}")
        canvas.restoreState()

    portrait_frame = Frame(
        margin, margin, page_w - 2 * margin, page_h - 2 * margin, id="portrait_frame"
    )
    portrait_template = PageTemplate(
        id="portrait", frames=[portrait_frame], onPage=_footer_portrait
    )

    land_frame = Frame(
        margin, margin, land_w - 2 * margin, land_h - 2 * margin, id="landscape_frame"
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
        leftMargin=margin,
        rightMargin=margin,
        topMargin=margin,
        bottomMargin=margin,
    )
    doc.addPageTemplates([portrait_template, landscape_template])
    doc.build(story)

    for p in tmp_paths:
        try:
            os.unlink(p)
        except Exception:
            pass

    return {"status": "success", "report_path": str(output_path)}
