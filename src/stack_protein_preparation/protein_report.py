"""
Per-protein PDF report for the FRUTON protein preparation pipeline.

Generates one PDF per protein containing:
  - A PyMOL-rendered figure of the starting structure (centered on ligand
    when one is present)
  - A prose narrative for every pipeline step that ran, with inline citations
  - A bibliography section with full references
  - A shared BibTeX file written to the proteins directory

Entry point: generate_protein_report()
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

from ._report_assets import find_idis_logo, find_fruton_logo
from ._report_bibliography import (
    _INLINE_REFS,
    _fetch_pdb_citation,
    _format_pdb_citation_inline,
    _collect_ref_keys,
    _write_bib,
)
from ._report_narratives import _build_narratives, _cite
from ._report_pymol import (
    _render_pymol_figure,
    _render_metal_pocket_figures,
    _render_nonstd_mol2_figure,
    _parse_resp_mol2,
    _find_metal_h_contacts,
)
from ._report_pdf_builders import (
    _safe_pdf_text,
    _record_text,
    _boolish_yes,
    _discover_fasta_records,
    _discover_alignment_assets,
    _short_path_for_report,
    _collect_path_evidence,
    _collect_prefixed_evidence,
    _preparation_decision,
    _status_palette_for_value,
    _quality_summary,
    _make_page_chrome,
    _section_block,
    _build_summary_table,
    _build_narrative_table,
    _build_simple_evidence_table,
    _collect_residue_renames,
    _build_residue_rename_table,
    _build_sequence_evidence_flowables,
)
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
    _FRUTON_YELLOW_SOFT,
    _FRUTON_GREEN_SOFT,
    _AUDIT_GREY,
)


# ---------------------------------------------------------------------------
# PDF builder (lazy reportlab import)
# ---------------------------------------------------------------------------

def _build_pdf(
    pdb_id: str,
    pipeline_record: dict[str, str],
    output_pdf: Path,
    figure_png: Path | None,
    ref_keys: list[str],
    metal_png_h: Path | None = None,
    metal_png_dist: Path | None = None,
    ion_type: str = "",
    pdb_citation_inline: str = "",
    nonstd_residue_data: list[dict] | None = None,
    rerun_markers: list[dict] | None = None,
) -> None:
    from reportlab.lib import colors
    from reportlab.lib.pagesizes import A4
    from reportlab.lib.styles import ParagraphStyle, getSampleStyleSheet
    from reportlab.lib.units import cm
    from reportlab.platypus import (
        HRFlowable,
        Image,
        KeepTogether,
        Paragraph,
        SimpleDocTemplate,
        Spacer,
        Table,
        TableStyle,
    )

    report_dir = output_pdf.parent
    protein_dir = report_dir.parent
    idis_logo = find_idis_logo(protein_dir)
    fruton_logo = find_fruton_logo(protein_dir)

    doc = SimpleDocTemplate(
        str(output_pdf),
        pagesize=A4,
        leftMargin=1.75 * cm,
        rightMargin=1.75 * cm,
        topMargin=2.25 * cm,
        bottomMargin=1.70 * cm,
        title=f"FRUTON protein preparation audit: {pdb_id}",
        author="FRUTON pipeline",
        subject="Per-protein preparation report",
    )
    usable_width = A4[0] - doc.leftMargin - doc.rightMargin
    styles = getSampleStyleSheet()

    styles.add(ParagraphStyle(
        "Eyebrow",
        parent=styles["Normal"],
        fontName="Helvetica-Bold",
        fontSize=7.5,
        leading=9,
        textColor=colors.HexColor(_IDIS_GOLD),
        spaceAfter=2,
        uppercase=True,
    ))
    styles.add(ParagraphStyle(
        "ReportTitle",
        parent=styles["Title"],
        fontName="Helvetica-Bold",
        fontSize=21,
        leading=24,
        alignment=0,
        spaceBefore=0,
        spaceAfter=4,
        textColor=colors.HexColor(_IDIS_NAVY),
    ))
    styles.add(ParagraphStyle(
        "ReportSubtitle",
        parent=styles["Normal"],
        fontSize=8.8,
        leading=12,
        textColor=colors.HexColor(_IDIS_MUTED),
        spaceAfter=8,
    ))
    styles.add(ParagraphStyle(
        "SectionHeading",
        parent=styles["Heading2"],
        fontName="Helvetica-Bold",
        fontSize=12.4,
        leading=15,
        spaceBefore=0,
        spaceAfter=1,
        textColor=colors.HexColor(_IDIS_NAVY),
    ))
    styles.add(ParagraphStyle(
        "SectionSubtitle",
        parent=styles["Normal"],
        fontSize=7.8,
        leading=10,
        textColor=colors.HexColor(_IDIS_MUTED),
        spaceAfter=0,
    ))
    styles.add(ParagraphStyle(
        "BodyAudit",
        parent=styles["Normal"],
        fontSize=9.2,
        leading=13.4,
        spaceAfter=6,
        textColor=colors.HexColor(_IDIS_TEXT),
    ))
    styles.add(ParagraphStyle(
        "CaptionAudit",
        parent=styles["Normal"],
        fontSize=7.5,
        leading=9.5,
        textColor=colors.HexColor(_IDIS_MUTED),
        spaceAfter=2,
    ))
    styles.add(ParagraphStyle(
        "ReferenceAudit",
        parent=styles["Normal"],
        fontSize=7.1,
        leading=9.2,
        leftIndent=8,
        firstLineIndent=-8,
        spaceAfter=2.4,
        textColor=colors.HexColor(_IDIS_TEXT),
    ))
    styles.add(ParagraphStyle(
        "MetricLabel",
        parent=styles["Normal"],
        fontName="Helvetica-Bold",
        fontSize=6.5,
        leading=8,
        textColor=colors.HexColor(_IDIS_MUTED),
    ))
    styles.add(ParagraphStyle(
        "MetricValue",
        parent=styles["Normal"],
        fontName="Helvetica-Bold",
        fontSize=8.3,
        leading=10.5,
        textColor=colors.HexColor(_IDIS_TEXT),
    ))
    styles.add(ParagraphStyle(
        "SmallNote",
        parent=styles["Normal"],
        fontSize=7.4,
        leading=9.5,
        textColor=colors.HexColor(_IDIS_MUTED),
    ))

    styles.add(ParagraphStyle(
        "MonoSequence",
        parent=styles["Normal"],
        fontName="Courier",
        fontSize=6.7,
        leading=8.2,
        textColor=colors.HexColor(_IDIS_TEXT),
        spaceAfter=3,
    ))
    styles.add(ParagraphStyle(
        "EvidenceHeader",
        parent=styles["Normal"],
        fontName="Helvetica-Bold",
        fontSize=7.0,
        leading=8.5,
        textColor=colors.HexColor(_IDIS_NAVY),
    ))
    styles.add(ParagraphStyle(
        "EvidenceCell",
        parent=styles["Normal"],
        fontSize=7.0,
        leading=9.0,
        textColor=colors.HexColor(_IDIS_TEXT),
    ))

    thin_rule = HRFlowable(
        width="100%",
        thickness=0.6,
        color=colors.HexColor(_IDIS_LINE),
        spaceBefore=8,
        spaceAfter=6,
    )

    story: list[Any] = []

    # ------------------------------------------------------------------
    # Title and summary block
    # ------------------------------------------------------------------
    story.append(Paragraph("FRUTON · IDIS visual report", styles["Eyebrow"]))
    story.append(Paragraph(f"Protein preparation audit: {_safe_pdf_text(pdb_id)}", styles["ReportTitle"]))

    uniprot = _record_text(pipeline_record, "uniprot_id")
    residue_range = _record_text(pipeline_record, "range")
    actual_range = _record_text(pipeline_record, "prepared_structure.actual_range")
    subtitle_parts = []
    if uniprot:
        subtitle_parts.append(f"UniProt: {_safe_pdf_text(uniprot)}")
    if residue_range and actual_range and actual_range != residue_range:
        subtitle_parts.append(
            f"Range: {_safe_pdf_text(residue_range)} "
            f"(actual: {_safe_pdf_text(actual_range)})"
        )
    elif residue_range:
        subtitle_parts.append(f"Range: {_safe_pdf_text(residue_range)}")
    subtitle_parts.append(_safe_pdf_text(_quality_summary(pipeline_record)))
    story.append(Paragraph(" &nbsp;·&nbsp; ".join(subtitle_parts), styles["ReportSubtitle"]))

    story.append(_build_summary_table(pipeline_record, styles, usable_width))
    story.append(Spacer(1, 0.32 * cm))

    # ------------------------------------------------------------------
    # Structure figure
    # ------------------------------------------------------------------
    story.append(_section_block(
        "Starting structure",
        "Initial coordinate evidence rendered for visual inspection before chemistry-specific interpretation.",
        styles,
        usable_width,
    ))
    story.append(Spacer(1, 0.08 * cm))
    if figure_png and figure_png.exists():
        ligand_note = (
            "centered on ligand"
            if _boolish_yes(_record_text(pipeline_record, "has_ligands"))
            else "default orientation"
        )
        img = Image(str(figure_png), width=13.7 * cm, height=10.25 * cm)
        fig_table = Table([[img]], colWidths=[usable_width])
        fig_table.setStyle(TableStyle([
            ("ALIGN", (0, 0), (-1, -1), "CENTER"),
            ("BACKGROUND", (0, 0), (-1, -1), colors.HexColor(_IDIS_PANEL)),
            ("BOX", (0, 0), (-1, -1), 0.6, colors.HexColor(_IDIS_LINE)),
            ("TOPPADDING", (0, 0), (-1, -1), 8),
            ("BOTTOMPADDING", (0, 0), (-1, -1), 8),
        ]))
        story.append(fig_table)
        story.append(Paragraph(
            f"<i>Figure 1. PyMOL{_cite('schrodinger2015pymol')} rendering of the starting "
            f"structure ({_safe_pdf_text(pdb_id)}, {_safe_pdf_text(ligand_note)}). "
            "Helices use the IDIS navy; strands use the IDIS gold accent; "
            "loops are muted grey; non-water heteroatoms are shown as sticks.</i>",
            styles["CaptionAudit"],
        ))
    else:
        story.append(Paragraph(
            "Structure figure not available. PyMOL rendering was skipped or failed, but the textual audit remains usable.",
            styles["BodyAudit"],
        ))
    story.append(Spacer(1, 0.28 * cm))

    # ------------------------------------------------------------------
    # Sequence and alignment evidence
    # ------------------------------------------------------------------
    figure_num = 2
    fasta_records = _discover_fasta_records(protein_dir, pdb_id)
    story.append(_section_block(
        "Sequence evidence: FASTA records",
        "Canonical, deposited, and coordinate-resolved sequences used by FRUTON before alignment and gap interpretation.",
        styles,
        usable_width,
    ))
    story.append(Spacer(1, 0.08 * cm))
    if fasta_records:
        story.extend(_build_sequence_evidence_flowables(fasta_records, styles, usable_width, protein_dir))
    else:
        story.append(Paragraph(
            "No FASTA files were found under the protein folder. This usually means the report was generated without the fasta_files.py outputs or from a reduced artifact bundle.",
            styles["BodyAudit"],
        ))
    story.append(Spacer(1, 0.22 * cm))

    alignment_assets = _discover_alignment_assets(protein_dir)
    if alignment_assets["images"] or alignment_assets["mappings"] or alignment_assets["alignments"]:
        story.append(_section_block(
            "Alignment evidence",
            "Files produced by the PDB-vs-UniProt alignment layer; figures expose matches, gaps, and unresolved regions.",
            styles,
            usable_width,
        ))
        story.append(Spacer(1, 0.08 * cm))
        alignment_rows = [
            ("alignment PNG files", ", ".join(_short_path_for_report(path, protein_dir) for path in alignment_assets["images"]) or "none"),
            ("mapping TSV files", ", ".join(_short_path_for_report(path, protein_dir) for path in alignment_assets["mappings"][:6]) or "none"),
            ("alignment FASTA files", ", ".join(_short_path_for_report(path, protein_dir) for path in alignment_assets["alignments"][:6]) or "none"),
        ]
        story.append(_build_simple_evidence_table(alignment_rows, styles, usable_width, left_header="Evidence", right_header="Files"))
        story.append(Spacer(1, 0.12 * cm))
        for alignment_png in alignment_assets["images"]:
            story.append(KeepTogether([
                Image(str(alignment_png), width=13.7 * cm, height=5.4 * cm),
                Paragraph(
                    f"<i>Figure {figure_num}. Sequence-alignment visualization generated by FRUTON. "
                    "The figure should be used to verify sequence coverage, terminal padding, internal gaps, and chain-specific PDB-to-UniProt consistency before interpreting gap-filling decisions.</i>",
                    styles["CaptionAudit"],
                ),
            ]))
            figure_num += 1
            story.append(Spacer(1, 0.14 * cm))
        story.append(Spacer(1, 0.18 * cm))

    # ------------------------------------------------------------------
    # Pipeline narrative
    # ------------------------------------------------------------------
    story.append(_section_block(
        "Preparation pipeline",
        "Step narrative follows the FRUTON runner order and keeps citations local to the action being reported.",
        styles,
        usable_width,
    ))
    story.append(Spacer(1, 0.08 * cm))

    # Build sorted list of (step, timestamp) rerun markers -- each becomes a red box.
    _markers: list[tuple[int, str]] = []
    for m in (rerun_markers or []):
        try:
            _markers.append((int(m["step"]), str(m.get("timestamp", ""))))
        except (KeyError, TypeError, ValueError):
            pass
    _markers.sort(key=lambda x: x[0])
    _pending_markers = list(_markers)  # consumed as we walk through paragraphs

    def _make_marker_box(step: int, ts: str) -> list[Any]:
        ts_note = f" ({_safe_pdf_text(ts)})" if ts else ""
        t = Table(
            [[Paragraph(
                f"FILES CHANGED BY USER. — pipeline re-run from step {step}{ts_note}",
                styles["MetricValue"],
            )]],
            colWidths=[usable_width - 0.9 * cm],
        )
        t.setStyle(TableStyle([
            ("BACKGROUND", (0, 0), (-1, -1), colors.HexColor(_FRUTON_RED_SOFT)),
            ("TEXTCOLOR", (0, 0), (-1, -1), colors.HexColor(_FRUTON_RED)),
            ("BOX", (0, 0), (-1, -1), 1.2, colors.HexColor(_FRUTON_RED)),
            ("LEFTPADDING", (0, 0), (-1, -1), 8),
            ("RIGHTPADDING", (0, 0), (-1, -1), 8),
            ("TOPPADDING", (0, 0), (-1, -1), 6),
            ("BOTTOMPADDING", (0, 0), (-1, -1), 6),
        ]))
        return [t, Spacer(1, 0.08 * cm)]

    narrative_flowables: list[Any] = []
    tagged = _build_narratives(pdb_id, pipeline_record, pdb_citation_inline=pdb_citation_inline)
    if tagged:
        for step_num, para_text in tagged:
            # Insert all pending markers whose step <= current narrative step
            while _pending_markers and _pending_markers[0][0] <= step_num:
                m_step, m_ts = _pending_markers.pop(0)
                narrative_flowables.extend(_make_marker_box(m_step, m_ts))
            narrative_flowables.append(Paragraph(para_text, styles["BodyAudit"]))
        # Any markers beyond the last narrative step go at the end
        for m_step, m_ts in _pending_markers:
            narrative_flowables.extend(_make_marker_box(m_step, m_ts))
    else:
        for m_step, m_ts in _pending_markers:
            narrative_flowables.extend(_make_marker_box(m_step, m_ts))
        narrative_flowables.append(Paragraph("No completed pipeline steps were recorded.", styles["BodyAudit"]))
    story.append(_build_narrative_table(narrative_flowables, usable_width))
    story.append(Spacer(1, 0.30 * cm))

    # ------------------------------------------------------------------
    # Protonation residue rename table
    # ------------------------------------------------------------------
    prot_in_path_str = _record_text(pipeline_record, "protonation.input_path")
    prot_out_path_str = _record_text(pipeline_record, "protonation.output_path")
    if prot_in_path_str and prot_out_path_str:
        prot_in = Path(prot_in_path_str)
        prot_out = Path(prot_out_path_str)
        residue_renames = _collect_residue_renames(prot_in, prot_out)
        if residue_renames:
            story.append(_section_block(
                "Protonation: residue assignments",
                "All residue name changes introduced by the PROPKA + gmx pdb2gmx protonation step. "
                "PROPKA-driven rows (highlighted) reflect pKa-based histidine state selection. "
                "GROMACS-driven rows reflect force-field template substitutions (e.g. disulfide "
                "bridges, selenomethionine conversion).",
                styles,
                usable_width,
            ))
            story.append(Spacer(1, 0.08 * cm))
            story.append(_build_residue_rename_table(residue_renames, styles, usable_width))
            story.append(Spacer(1, 0.22 * cm))

    # ------------------------------------------------------------------
    # Explicit audit decision and file evidence
    # ------------------------------------------------------------------
    decision, reason = _preparation_decision(pipeline_record)
    decision_rows = [
        ("preparation decision", decision),
        ("reason", reason),
        ("prepared status", _record_text(pipeline_record, "prepared_structure.status") or "not recorded"),
        ("prepared variant", _record_text(pipeline_record, "prepared_structure.variant") or "not recorded"),
        ("Ramachandran favored", (_record_text(pipeline_record, "model_eval.rama_pct_favored") + "%") if _record_text(pipeline_record, "model_eval.rama_pct_favored") else "not recorded"),
        ("Ramachandran outlier", (_record_text(pipeline_record, "model_eval.rama_pct_outlier") + "%") if _record_text(pipeline_record, "model_eval.rama_pct_outlier") else "not recorded"),
        ("clashscore", _record_text(pipeline_record, "model_eval.clashscore") or "not recorded"),
    ]
    story.append(_section_block(
        "Preparation decision",
        "Conservative audit status derived from prepared-structure, quality, metal, and non-standard-residue evidence.",
        styles,
        usable_width,
    ))
    story.append(Spacer(1, 0.08 * cm))
    story.append(_build_simple_evidence_table(decision_rows, styles, usable_width, left_header="Decision field", right_header="Recorded value"))
    story.append(Spacer(1, 0.22 * cm))

    path_rows = _collect_path_evidence(pipeline_record, protein_dir)
    if path_rows:
        story.append(_section_block(
            "Output and file evidence",
            "Important paths recorded by the flattened FRUTON state; these are the files a reviewer should open before accepting a structure.",
            styles,
            usable_width,
        ))
        story.append(Spacer(1, 0.08 * cm))
        story.append(_build_simple_evidence_table(path_rows, styles, usable_width, left_header="Record key", right_header="Path"))
        story.append(Spacer(1, 0.22 * cm))

    metal_evidence_rows = _collect_prefixed_evidence(
        pipeline_record,
        ("metals.", "metall_params.", "nonstd_residue_params.", "gaussian_params."),
    )
    if _boolish_yes(_record_text(pipeline_record, "has_metals")) or _boolish_yes(_record_text(pipeline_record, "has_nonstandard_residues")) or metal_evidence_rows:
        story.append(_section_block(
            "Parameterization evidence (metals, non-standard residues, Gaussian)",
            "Metal detection, non-standard residue identification, Gaussian QM calculations, and generated force-field parameters.",
            styles,
            usable_width,
        ))
        story.append(Spacer(1, 0.08 * cm))
        if not metal_evidence_rows:
            metal_evidence_rows = [("has_metals", _record_text(pipeline_record, "has_metals") or "yes")]
        story.append(_build_simple_evidence_table(metal_evidence_rows[:18], styles, usable_width, left_header="Record key", right_header="Value"))
        story.append(Paragraph(
            "Interpretation rule: a rendered metal pocket is evidence of detection and inspection, not proof that bonded metal parameters were accepted. MCPB.py/Gaussian acceptance should be based on explicit status fields and generated force-field files.",
            styles["SmallNote"],
        ))
        story.append(Spacer(1, 0.22 * cm))

    # ------------------------------------------------------------------
    # Metal coordination pocket figures
    # ------------------------------------------------------------------
    if metal_png_h or metal_png_dist:
        ion_label = f" ({_safe_pdf_text(ion_type)})" if ion_type else ""
        story.append(_section_block(
            f"Metal coordination pocket{ion_label}",
            "Local metal-site view for protonation-state and donor-geometry inspection. "
            "Coordinating water molecules are shown in blue; protein residues in gold.",
            styles,
            usable_width,
        ))
        story.append(Spacer(1, 0.08 * cm))

        if metal_png_h and metal_png_h.exists():
            story.append(KeepTogether([
                Image(str(metal_png_h), width=9.6 * cm, height=9.6 * cm),
                Paragraph(
                    f"<i>Figure {figure_num}. Metal coordination pocket{ion_label} with all hydrogen atoms visible. "
                    "The protein backbone is a transparent grey cartoon, the metal is a FRUTON-red sphere, "
                    "coordinating protein residues use the gold carbon accent, and water molecules are shown in "
                    "light blue (oxygen) and white (hydrogen). "
                    "Hydrogen positions expose donor protonation states and water orientation relative to the metal. "
                    "Water H atoms pointing toward the metal require reorientation before MCPB.py/Gaussian.</i>",
                    styles["CaptionAudit"],
                ),
            ]))
            figure_num += 1
            story.append(Spacer(1, 0.20 * cm))

        if metal_png_dist and metal_png_dist.exists():
            story.append(KeepTogether([
                Image(str(metal_png_dist), width=9.6 * cm, height=9.6 * cm),
                Paragraph(
                    f"<i>Figure {figure_num}. Metal coordination geometry{ion_label}. "
                    "Dashed FRUTON-yellow lines show distances (A) between the metal ion and donor heavy atoms "
                    "(N, O, S) within the direct-donor cutoff. Coordinating water oxygens are shown as blue sticks. "
                    "Hydrogen atoms are shown; water H atoms pointing toward the metal (M-H &lt; M-O) "
                    "are candidates for reorientation (see table below).</i>",
                    styles["CaptionAudit"],
                ),
            ]))
            figure_num += 1
            story.append(Spacer(1, 0.20 * cm))

        # ------ H-contact table (water H atoms toward metal) ------
        prepared_out_for_h = str(pipeline_record.get("prepared_structure.output_path", "")).strip()
        h_contacts: list[dict] = []
        if prepared_out_for_h:
            _prep_pdb_for_h = Path(prepared_out_for_h)
            if _prep_pdb_for_h.exists():
                h_contacts = _find_metal_h_contacts(_prep_pdb_for_h)

        if h_contacts:
            from reportlab.platypus import Paragraph as _Para, Table as _Tbl, TableStyle as _TS
            from reportlab.lib import colors as _colors

            story.append(Paragraph(
                "Water H atoms requiring reorientation",
                styles["SectionSubtitle"],
            ))
            story.append(Paragraph(
                "The table below lists hydrogen atoms whose metal-H distance (M-H) is shorter than "
                "the corresponding metal-oxygen distance (M-O).  These H atoms point toward the "
                "metal rather than away from it, which is physically incorrect for a coordinated "
                "water molecule.  They must be reoriented (rotated away from the metal) before "
                "running MCPB.py or Gaussian — they should never be removed, as the water "
                "molecule remains a metal ligand.",
                styles["SmallNote"],
            ))
            story.append(Spacer(1, 0.08 * cm))

            hdr = [
                _Para("Metal", styles["EvidenceHeader"]),
                _Para("Residue", styles["EvidenceHeader"]),
                _Para("H atom", styles["EvidenceHeader"]),
                _Para("Parent (O/N)", styles["EvidenceHeader"]),
                _Para("M-H (A)", styles["EvidenceHeader"]),
                _Para("M-O (A)", styles["EvidenceHeader"]),
                _Para("Status", styles["EvidenceHeader"]),
            ]
            rows: list[list] = [hdr]
            _toward = [r for r in h_contacts if r["status"] == "toward_metal"]
            _near   = [r for r in h_contacts if r["status"] != "toward_metal"]
            for r in _toward + _near:
                mo_str = f"{r['mo_dist']:.3f}" if r["mo_dist"] is not None else "n/a"
                chain_str = r["h_chain"] or "—"
                res_str = f"{r['h_res']} {chain_str}{r['h_resnum']}"
                status_str = "toward metal" if r["status"] == "toward_metal" else "near metal"
                rows.append([
                    _Para(_safe_pdf_text(r["metal_label"]), styles["EvidenceCell"]),
                    _Para(_safe_pdf_text(res_str), styles["EvidenceCell"]),
                    _Para(_safe_pdf_text(r["h_name"]), styles["EvidenceCell"]),
                    _Para(_safe_pdf_text(r["parent_o_name"]), styles["EvidenceCell"]),
                    _Para(_safe_pdf_text(f"{r['mh_dist']:.3f}"), styles["EvidenceCell"]),
                    _Para(_safe_pdf_text(mo_str), styles["EvidenceCell"]),
                    _Para(_safe_pdf_text(status_str), styles["EvidenceCell"]),
                ])

            col_w = [usable_width * f for f in (0.22, 0.18, 0.10, 0.12, 0.12, 0.12, 0.14)]
            ts = _TS([
                ("BACKGROUND",  (0, 0), (-1, 0), _colors.HexColor("#2D3259")),
                ("TEXTCOLOR",   (0, 0), (-1, 0), _colors.white),
                ("FONTSIZE",    (0, 0), (-1, -1), 7),
                ("ROWBACKGROUNDS", (0, 1), (-1, -1),
                 [_colors.HexColor("#F5F5F5"), _colors.white]),
                ("GRID", (0, 0), (-1, -1), 0.3, _colors.HexColor("#CCCCCC")),
                ("VALIGN", (0, 0), (-1, -1), "MIDDLE"),
                ("TOPPADDING",  (0, 0), (-1, -1), 2),
                ("BOTTOMPADDING", (0, 0), (-1, -1), 2),
            ])
            # Highlight "toward_metal" rows in light orange
            for i, r in enumerate(_toward, start=1):
                ts.add("BACKGROUND", (0, i), (-1, i), _colors.HexColor("#FFF0D0"))
            tbl = _Tbl(rows, colWidths=col_w)
            tbl.setStyle(ts)
            story.append(tbl)
            story.append(Paragraph(
                "Highlighted rows (toward metal) require water-molecule rotation before "
                "MCPB.py/Gaussian parametrization.",
                styles["SmallNote"],
            ))
            story.append(Spacer(1, 0.20 * cm))

        story.append(Spacer(1, 0.10 * cm))

    # ------------------------------------------------------------------
    # Ramachandran scatter plot
    # ------------------------------------------------------------------
    rama_plot_path_str = _record_text(pipeline_record, "model_eval.rama_plot_path")
    rama_plot_file = Path(rama_plot_path_str) if rama_plot_path_str else None
    if rama_plot_file and rama_plot_file.exists():
        story.append(_section_block(
            "Structure quality: Ramachandran plot",
            "Stereochemical quality evidence for the prepared or model-completed structure.",
            styles,
            usable_width,
        ))
        story.append(Spacer(1, 0.08 * cm))
        story.append(Image(str(rama_plot_file), width=9.8 * cm, height=9.8 * cm))
        pct_fav_val = _record_text(pipeline_record, "model_eval.rama_pct_favored")
        pct_out_val = _record_text(pipeline_record, "model_eval.rama_pct_outlier")
        cs_val = _record_text(pipeline_record, "model_eval.clashscore")
        quality_val = _record_text(pipeline_record, "model_eval.overall_quality")
        stats_line = ""
        if pct_fav_val:
            stats_line = (
                f" Ramachandran: {_safe_pdf_text(pct_fav_val)}% favored, "
                f"{_safe_pdf_text(pct_out_val)}% outlier. Clashscore: {_safe_pdf_text(cs_val)}. "
                f"Overall quality: {_safe_pdf_text(quality_val)}."
            )
        story.append(Paragraph(
            f"<i>Figure {figure_num}. Ramachandran phi/psi scatter plot"
            f"{_cite('lovell2003rama', 'williams2018molprobity', 'laskowski1993procheck')}. "
            "Background: green = favored region, yellow = allowed region. "
            "Residue markers: circle = general, triangle = Gly, diamond = Pro. "
            "Dot colour: dark green = favored, orange = allowed, red = outlier."
            + stats_line + "</i>",
            styles["CaptionAudit"],
        ))
        figure_num += 1
        story.append(Spacer(1, 0.30 * cm))

    # ------------------------------------------------------------------
    # Non-standard residue RESP charge figures
    # ------------------------------------------------------------------
    if nonstd_residue_data:
        from reportlab.platypus import Paragraph as _Para, Table as _Tbl, TableStyle as _TS
        from reportlab.lib import colors as _colors

        story.append(_section_block(
            "Non-standard residue parameters",
            "Atom types and RESP partial charges for each non-standard residue. "
            "Charges derive from two-stage RESP fitting to HF/6-31G* electrostatic potentials "
            "or from a curated force-field database.",
            styles,
            usable_width,
        ))
        story.append(Spacer(1, 0.08 * cm))

        for _nrd in nonstd_residue_data:
            _resname      = _nrd["resname"]
            _label        = _nrd["label"]
            _atoms        = _nrd["atoms"]
            _total_charge = _nrd["total_charge"]
            _formal_q     = _nrd.get("formal_charge")
            _src          = _nrd.get("charge_source", "")
            _in_db        = _nrd.get("in_database", False)
            _fig_png      = _nrd.get("figure_png")

            _src_note = ""
            if _in_db:
                _src_note = "Charges from curated force-field database."
            elif _src:
                _src_note = f"Charge source: {_src}."

            _formal_note = ""
            if _formal_q is not None:
                _formal_note = f" Formal charge: {_formal_q:+d}."

            _sub_header = (
                f"<b>{_safe_pdf_text(_resname)}</b> ({_safe_pdf_text(_label)}) — "
                f"RESP total charge: {_total_charge:+.4f}.{_formal_note}"
            )
            if _src_note:
                _sub_header += f" {_safe_pdf_text(_src_note)}"

            story.append(Paragraph(_sub_header, styles["SectionSubtitle"]))
            story.append(Spacer(1, 0.05 * cm))

            # Build atom charge table (always shown alongside the figure)
            _tbl_hdr = [
                _Para("Atom", styles["EvidenceHeader"]),
                _Para("SYBYL type", styles["EvidenceHeader"]),
                _Para("RESP charge (e)", styles["EvidenceHeader"]),
            ]
            _tbl_rows: list[list] = [_tbl_hdr]
            for _at in _atoms:
                _q = _at["charge"]
                _q_str = f"{_q:+.6f}"
                _tbl_rows.append([
                    _Para(_safe_pdf_text(_at["name"]), styles["EvidenceCell"]),
                    _Para(_safe_pdf_text(_at["sybyl_type"]), styles["EvidenceCell"]),
                    _Para(_safe_pdf_text(_q_str), styles["EvidenceCell"]),
                ])
            _tbl_rows.append([
                _Para("<b>Total</b>", styles["EvidenceCell"]),
                _Para("", styles["EvidenceCell"]),
                _Para(f"<b>{_total_charge:+.6f}</b>", styles["EvidenceCell"]),
            ])

            _charge_ts = _TS([
                ("BACKGROUND",    (0, 0), (-1, 0), _colors.HexColor("#2D3259")),
                ("TEXTCOLOR",     (0, 0), (-1, 0), _colors.white),
                ("FONTSIZE",      (0, 0), (-1, -1), 7),
                ("ROWBACKGROUNDS",(0, 1), (-1, -2),
                 [_colors.HexColor("#F5F5F5"), _colors.white]),
                ("BACKGROUND",    (0, -1), (-1, -1), _colors.HexColor(_IDIS_PANEL_2)),
                ("GRID",          (0, 0), (-1, -1), 0.3, _colors.HexColor("#CCCCCC")),
                ("VALIGN",        (0, 0), (-1, -1), "MIDDLE"),
                ("TOPPADDING",    (0, 0), (-1, -1), 2),
                ("BOTTOMPADDING", (0, 0), (-1, -1), 2),
            ])

            _fig_png_path = Path(_fig_png) if _fig_png else None
            if _fig_png_path and _fig_png_path.exists():
                # Figure is 1200x900 (4:3); allocate ~60% of page width for it
                _img_w = usable_width * 0.60
                _img_h = _img_w * (900 / 1200)
                _img = Image(str(_fig_png_path), width=_img_w, height=_img_h)
                _tbl_w = usable_width - _img_w - 0.3 * cm
                _col_w = [_tbl_w * 0.28, _tbl_w * 0.32, _tbl_w * 0.40]
                _charge_tbl = _Tbl(_tbl_rows, colWidths=_col_w)
                _charge_tbl.setStyle(_charge_ts)
                _side_tbl = _Tbl(
                    [[_img, _charge_tbl]],
                    colWidths=[_img_w, usable_width - _img_w],
                    style=_TS([("VALIGN", (0, 0), (-1, -1), "TOP")]),
                )
                story.append(_side_tbl)
            else:
                _col_w = [usable_width * 0.28, usable_width * 0.30, usable_width * 0.30]
                _charge_tbl = _Tbl(_tbl_rows, colWidths=_col_w)
                _charge_tbl.setStyle(_charge_ts)
                story.append(_charge_tbl)

            story.append(Paragraph(
                f"<i>Figure {figure_num}. {_safe_pdf_text(_resname)} ({_safe_pdf_text(_label)}): "
                "atoms coloured by RESP partial charge (red = negative, white = ~0, blue = positive); "
                "heavy atoms labelled with atom name and charge value (white label background). "
                f"RESP total charge = {_total_charge:+.4f} e.</i>",
                styles["CaptionAudit"],
            ))
            figure_num += 1
            story.append(Spacer(1, 0.25 * cm))

    # ------------------------------------------------------------------
    # References
    # ------------------------------------------------------------------
    story.append(thin_rule)
    story.append(_section_block(
        "References",
        "Only methods triggered by the recorded pipeline state are included here.",
        styles,
        usable_width,
    ))
    seen: set[str] = set()
    for key in ref_keys:
        if key in _INLINE_REFS and key not in seen:
            seen.add(key)
            story.append(Paragraph(
                f"[{_safe_pdf_text(key)}]&nbsp;&nbsp;{_safe_pdf_text(_INLINE_REFS[key])}",
                styles["ReferenceAudit"],
            ))
    if not seen:
        story.append(Paragraph("No references collected.", styles["ReferenceAudit"]))

    chrome = _make_page_chrome(pdb_id=pdb_id, idis_logo=idis_logo, fruton_logo=fruton_logo)
    doc.build(story, onFirstPage=chrome, onLaterPages=chrome)


# ---------------------------------------------------------------------------
# Public entry point
# ---------------------------------------------------------------------------

def generate_protein_report(
    protein_dir: Path,
    pdb_id: str,
    pipeline_record: dict[str, str],
    *,
    global_bib_path: Path | None = None,
) -> dict[str, Any]:
    """
    Generate a per-protein PDF report with a PyMOL figure, step narratives,
    and bibliography.

    Returns
    -------
    dict with keys: status, report_pdf_path, figure_png_path, bib_path, message.
    """
    report_dir = protein_dir / "report"
    report_dir.mkdir(parents=True, exist_ok=True)

    output_pdf = report_dir / f"{pdb_id}_report.pdf"
    output_png = report_dir / f"{pdb_id}_structure.png"

    messages: list[str] = []
    result: dict[str, Any] = {
        "pdb_id":           pdb_id,
        "status":           "failed",
        "message":          "",
        "report_pdf_path":  str(output_pdf),
        "figure_png_path":  str(output_png),
        "bib_path":         str(global_bib_path) if global_bib_path else "",
    }

    # 0. Fetch PDB crystal paper citation from RCSB (best-effort, non-blocking)
    pdb_citation_inline = ""
    try:
        pdb_cit = _fetch_pdb_citation(pdb_id)
        if pdb_cit:
            pdb_citation_inline = _format_pdb_citation_inline(pdb_cit, pdb_id)
    except Exception as exc:
        messages.append(f"PDB citation fetch: {exc}")

    # 1. Render PyMOL overview figure (best-effort)
    fig_error = _render_pymol_figure(protein_dir, pdb_id, pipeline_record, output_png)
    if fig_error:
        messages.append(f"Figure: {fig_error}")

    # 2. Render metal pocket figures when metals were parametrised (best-effort)
    metal_png_h: Path | None = None
    metal_png_dist: Path | None = None
    ion_type_str = ""
    has_metals = str(pipeline_record.get("has_metals", "")).strip().lower()
    prepared_out = str(pipeline_record.get("prepared_structure.output_path", "")).strip()
    if has_metals in ("yes", "true") and prepared_out:
        prepared_pdb = Path(prepared_out)
        if not prepared_pdb.exists():
            # Fall back: look for any PDB under prepared/
            prep_dir = protein_dir / "prepared"
            hits = sorted(prep_dir.rglob("*.pdb")) if prep_dir.exists() else []
            prepared_pdb = hits[0] if hits else prepared_pdb

        if prepared_pdb.exists():
            ion_type_str = str(pipeline_record.get("metals.ion_type", "")).strip().upper()
            if not ion_type_str:
                ion_type_str = "METAL"
            metal_out_dir = report_dir / "metal_figures"
            try:
                pocket = _render_metal_pocket_figures(
                    pdb_path=prepared_pdb,
                    ion_type=ion_type_str,
                    output_dir=metal_out_dir,
                    pdb_id=pdb_id,
                )
                metal_png_h    = pocket.get("with_h")
                metal_png_dist = pocket.get("distances")
            except Exception as exc:
                messages.append(f"Metal figures: {exc}")

    # 3. Render non-standard residue mol2 figures and collect RESP charge data
    nonstd_residue_data: list[dict] = []
    nonstd_manifest_str = str(pipeline_record.get("nonstd_residue_params.manifest_path", "")).strip()
    if nonstd_manifest_str:
        nonstd_manifest_path = Path(nonstd_manifest_str)
        if nonstd_manifest_path.exists():
            try:
                import json as _json_nonstd
                _nonstd_manifest = _json_nonstd.loads(
                    nonstd_manifest_path.read_text(encoding="utf-8")
                )
                nonstd_fig_dir = report_dir / "nonstd_figures"
                nonstd_fig_dir.mkdir(parents=True, exist_ok=True)
                for _res in _nonstd_manifest.get("residues", []):
                    _label   = _res.get("label", "")
                    _resname = _res.get("resname", "")
                    _mol2 = (
                        nonstd_manifest_path.parent
                        / _label / "resp" / "step03_resp_params"
                        / f"{_resname}.mol2"
                    )
                    if not _mol2.exists():
                        continue
                    _atoms = _parse_resp_mol2(_mol2)
                    if not _atoms:
                        continue
                    _total_charge = sum(a["charge"] for a in _atoms)
                    _fig_png = nonstd_fig_dir / f"{_resname}_{_label}_resp.png"
                    _fig_err = _render_nonstd_mol2_figure(_mol2, _fig_png)
                    if _fig_err:
                        messages.append(f"Nonstd figure {_resname}: {_fig_err}")
                    nonstd_residue_data.append({
                        "resname":       _resname,
                        "label":         _label,
                        "formal_charge": _res.get("formal_charge"),
                        "charge_source": _res.get("charge_source", ""),
                        "in_database":   _res.get("in_database", False),
                        "figure_png":    _fig_png if _fig_png.exists() else None,
                        "atoms":         _atoms,
                        "total_charge":  _total_charge,
                    })
            except Exception as exc:
                messages.append(f"Nonstd residue data: {exc}")

    # 4. Collect citation keys for the bibliography section
    ref_keys = _collect_ref_keys(pipeline_record)

    # 5. Load per-protein rerun markers (accumulated across runs by --from-step)
    rerun_markers: list[dict] = []
    markers_path = protein_dir / "rerun_markers.json"
    try:
        if markers_path.is_file():
            import json as _json
            rerun_markers = _json.loads(markers_path.read_text())
    except Exception as exc:
        messages.append(f"rerun_markers.json read warning: {exc}")

    # 6. Build the PDF
    try:
        _build_pdf(
            pdb_id=pdb_id,
            pipeline_record=pipeline_record,
            output_pdf=output_pdf,
            figure_png=output_png if output_png.exists() else None,
            ref_keys=ref_keys,
            metal_png_h=metal_png_h,
            metal_png_dist=metal_png_dist,
            ion_type=ion_type_str,
            pdb_citation_inline=pdb_citation_inline,
            nonstd_residue_data=nonstd_residue_data,
            rerun_markers=rerun_markers,
        )
    except Exception as exc:
        messages.append(f"PDF generation failed: {exc}")
        result["message"] = " | ".join(messages)
        return result

    # 7. Write/update shared BibTeX file
    if global_bib_path:
        try:
            _write_bib(global_bib_path)
        except Exception as exc:
            messages.append(f"BibTeX write warning: {exc}")

    result["status"] = "success"
    result["message"] = " | ".join(messages) if messages else f"Report: {output_pdf.name}"
    return result
