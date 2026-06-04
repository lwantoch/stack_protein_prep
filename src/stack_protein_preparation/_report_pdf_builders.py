"""ReportLab helper functions and PDF table/flowable builders for FRUTON reports."""
from __future__ import annotations

from pathlib import Path
from typing import Any
from xml.sax.saxutils import escape

from ._report_palette import (
    _IDIS_NAVY,
    _IDIS_NAVY_2,
    _IDIS_TEXT,
    _IDIS_MUTED,
    _IDIS_LINE,
    _IDIS_PANEL,
    _IDIS_PANEL_2,
    _IDIS_GOLD,
    _FRUTON_RED,
    _FRUTON_RED_SOFT,
    _FRUTON_YELLOW,
    _FRUTON_YELLOW_SOFT,
    _FRUTON_GREEN,
    _FRUTON_GREEN_SOFT,
    _AUDIT_GREY,
    _AUDIT_WARN,
    _AUDIT_FAIL,
    _AUDIT_OK,
    _FRUTON_EXPANSION,
)
from ._report_assets import find_idis_logo, find_fruton_logo


# ---------------------------------------------------------------------------
# Text utilities
# ---------------------------------------------------------------------------

def _safe_pdf_text(value: Any) -> str:
    """Return text escaped for ReportLab Paragraph input."""
    return escape(str(value or "").strip())


def _record_text(record: dict[str, str], key: str, default: str = "") -> str:
    """Read one pipeline-record value as stripped text."""
    value = str(record.get(key, default) or default).strip()
    return value


def _boolish_yes(value: str) -> bool:
    """Return True for the yes/true/1 forms used in flattened pipeline state."""
    return value.strip().lower() in {"yes", "true", "1", "y"}


# ---------------------------------------------------------------------------
# FASTA / sequence helpers
# ---------------------------------------------------------------------------

def _read_fasta_records(fasta_path: Path) -> list[dict[str, str]]:
    """Read FASTA entries from *fasta_path* without adding a BioPython dependency here.

    The report module is intentionally a presentation layer.  It should consume the
    files produced by ``fasta_files.py`` and ``sequence_alignment.py`` but should not
    re-run sequence logic or assume one exact file naming convention.  A tiny FASTA
    parser is sufficient for report evidence because it only needs headers,
    sequences, lengths, and source paths.
    """
    records: list[dict[str, str]] = []
    header = ""
    seq_parts: list[str] = []

    try:
        lines = fasta_path.read_text(encoding="utf-8", errors="replace").splitlines()
    except OSError:
        return records

    def flush() -> None:
        nonlocal header, seq_parts
        if not header and not seq_parts:
            return
        sequence = "".join(part.strip() for part in seq_parts if part.strip()).replace(" ", "")
        if sequence:
            records.append({
                "header": header.strip() or fasta_path.stem,
                "sequence": sequence,
                "path": str(fasta_path),
            })
        header = ""
        seq_parts = []

    for line in lines:
        stripped = line.strip()
        if not stripped:
            continue
        if stripped.startswith(">"):
            flush()
            header = stripped[1:].strip()
        else:
            seq_parts.append(stripped)
    flush()
    return records


def _classify_fasta_source(path: Path, header: str) -> str:
    """Return a human-readable source label for a FASTA record."""
    text = f"{path.name} {header}".lower()
    if "uniprot" in text:
        return "UniProt canonical sequence"
    if "seqres" in text:
        return "PDB SEQRES deposited sequence"
    if "atom" in text:
        return "PDB ATOM resolved sequence"
    if "alignment" in [part.lower() for part in path.parts] or "aln" in text:
        return "Alignment FASTA"
    return "FASTA sequence"


def _wrap_sequence_for_pdf(sequence: str, line_width: int = 60, max_chars: int = 1400) -> tuple[str, bool]:
    """Wrap a sequence for a ReportLab paragraph and report whether it was truncated."""
    clean = "".join(sequence.split())
    truncated = len(clean) > max_chars
    if truncated:
        keep_tail = min(180, max_chars // 5)
        keep_head = max_chars - keep_tail
        clean = clean[:keep_head] + "..." + clean[-keep_tail:]
    lines = [clean[i:i + line_width] for i in range(0, len(clean), line_width)] or [""]
    return "<br/>".join(_safe_pdf_text(line) for line in lines), truncated


def _discover_fasta_records(protein_dir: Path, pdb_id: str, max_records: int = 8) -> list[dict[str, str]]:
    """Find the most useful FASTA records produced for one FRUTON protein folder."""
    fasta_roots = [protein_dir / "fasta", protein_dir]
    candidates: list[Path] = []
    seen: set[Path] = set()

    for root in fasta_roots:
        if not root.exists():
            continue
        for pattern in ("*.fasta", "*.fa", "*.faa"):
            for hit in root.rglob(pattern):
                lower_name = hit.name.lower()
                lower_parts = {part.lower() for part in hit.parts}
                if "alignments" in lower_parts:
                    continue
                if lower_name.endswith(".aln.fasta") or ".aln." in lower_name:
                    continue
                if "alignment" in lower_name or "_vs_" in lower_name:
                    continue
                try:
                    resolved = hit.resolve()
                except OSError:
                    resolved = hit
                if resolved not in seen:
                    seen.add(resolved)
                    candidates.append(hit)

    def sort_key(path: Path) -> tuple[int, str]:
        name = path.name.lower()
        if "uniprot" in name:
            priority = 0
        elif "atom" in name:
            priority = 1
        elif "seqres" in name:
            priority = 2
        else:
            priority = 3
        return (priority, str(path))

    records: list[dict[str, str]] = []
    for fasta_path in sorted(candidates, key=sort_key):
        for rec in _read_fasta_records(fasta_path):
            rec_path = Path(rec["path"])
            rec["source"] = _classify_fasta_source(rec_path, rec["header"])
            rec["length"] = str(len(rec["sequence"]))
            rec["unknown_x"] = str(rec["sequence"].upper().count("X"))
            records.append(rec)
            if len(records) >= max_records:
                return records
    return records


def _discover_alignment_assets(protein_dir: Path, max_images: int = 3) -> dict[str, list[Path]]:
    """Return alignment PNG and mapping TSV evidence generated by FRUTON."""
    alignment_roots = [protein_dir / "fasta" / "alignments", protein_dir / "alignments"]
    images: list[Path] = []
    mappings: list[Path] = []
    fastas: list[Path] = []

    for root in alignment_roots:
        if not root.exists():
            continue
        images.extend(sorted(root.glob("*.png")))
        mappings.extend(sorted(root.glob("*.mapping.tsv")))
        mappings.extend(sorted(root.glob("*mapping*.tsv")))
        fastas.extend(sorted(root.glob("*.aln.fasta")))

    # De-duplicate while preserving order.
    def unique(paths: list[Path]) -> list[Path]:
        out: list[Path] = []
        seen: set[Path] = set()
        for path in paths:
            try:
                resolved = path.resolve()
            except OSError:
                resolved = path
            if resolved not in seen:
                seen.add(resolved)
                out.append(path)
        return out

    return {
        "images": unique(images)[:max_images],
        "mappings": unique(mappings),
        "alignments": unique(fastas),
    }


# ---------------------------------------------------------------------------
# Path evidence helpers
# ---------------------------------------------------------------------------

def _short_path_for_report(path: str | Path, protein_dir: Path) -> str:
    """Return a compact path string for PDF tables."""
    p = Path(path)
    try:
        return str(p.relative_to(protein_dir))
    except ValueError:
        try:
            return str(p.relative_to(protein_dir.parent))
        except ValueError:
            return str(p)


def _collect_path_evidence(record: dict[str, str], protein_dir: Path, max_rows: int = 18) -> list[tuple[str, str]]:
    """Collect non-empty path-like fields from the flattened FRUTON record."""
    preferred = [
        "pdb_path",
        "processed_pdb_path",
        "fasta.seqres_path",
        "fasta.atom_path",
        "fasta.uniprot_path",
        "sequence_alignment.mapping_path",
        "components_dir",
        "components.protein_path",
        "components.water_path",
        "components.ligand_path",
        "components.metals_path",
        "filler.final_model_path",
        "protonation.output_path",
        "prepared_structure.output_path",
        "model_eval.rama_plot_path",
        "metall_params_directory",
        "metall_params.manifest_path",
        "nonstd_residue_params.manifest_path",
        "gaussian_params.manifest_path",
    ]
    rows: list[tuple[str, str]] = []
    seen_keys: set[str] = set()

    for key in preferred:
        value = _record_text(record, key)
        if value:
            rows.append((key, _short_path_for_report(value, protein_dir)))
            seen_keys.add(key)

    for key in sorted(record):
        if key in seen_keys:
            continue
        key_l = key.lower()
        if not any(token in key_l for token in ("path", "dir")):
            continue
        value = _record_text(record, key)
        if not value:
            continue
        rows.append((key, _short_path_for_report(value, protein_dir)))
        if len(rows) >= max_rows:
            break

    return rows[:max_rows]


def _collect_prefixed_evidence(record: dict[str, str], prefixes: tuple[str, ...]) -> list[tuple[str, str]]:
    """Collect non-empty flattened record values matching any prefix."""
    rows: list[tuple[str, str]] = []
    for key in sorted(record):
        if not key.startswith(prefixes):
            continue
        value = _record_text(record, key)
        if value:
            rows.append((key, value))
    return rows


# ---------------------------------------------------------------------------
# Preparation decision audit
# ---------------------------------------------------------------------------

def _preparation_decision(record: dict[str, str]) -> tuple[str, str]:
    """Return a conservative audit decision and its reason."""
    prep_status = _record_text(record, "prepared_structure.status").lower()
    quality = _record_text(record, "model_eval.overall_quality").lower()
    metals = _record_text(record, "has_metals").lower()
    metal_status = _record_text(record, "metall_params.status").lower()
    nonstd = _record_text(record, "has_nonstandard_residues").lower()
    nonstd_status = _record_text(record, "nonstd_residue_params.status").lower()

    gauss_status = _record_text(record, "gaussian_params.status").lower()
    needs_gaussian = (
        (metals in {"yes", "true"} and metal_status in {"success", "ready", "passed"})
        or (nonstd in {"yes", "true"} and nonstd_status in {"success", "ready", "passed"})
    )

    if prep_status in {"failed", "failure", "error"}:
        return "blocked", "prepared structure generation failed"
    if metals in {"yes", "true"} and metal_status not in {"success", "ready", "passed", "warning"}:
        return "requires metal-parameter review", "metal detected without accepted metal-parameter evidence"
    if nonstd in {"yes", "true"} and nonstd_status not in {"success", "ready", "passed", "warning"}:
        return "requires non-standard-residue review", "non-standard residue detected without accepted parameter evidence"
    if needs_gaussian and gauss_status not in {"success", "ready", "passed", "partial"}:
        return "requires Gaussian review", "MCPB/RESP Gaussian parametrization not yet complete or accepted"
    if quality in {"poor", "bad", "failed", "fail"}:
        return "requires structure-quality review", "Ramachandran/clash metrics are outside the preferred range"
    if prep_status in {"success", "warning"}:
        return "prepared with caveats", "review the evidence tables before downstream MD/GBSA"
    return "not enough evidence", "prepared-structure status or validation metrics are missing"


# ---------------------------------------------------------------------------
# Status palette
# ---------------------------------------------------------------------------

def _status_palette_for_value(field: str, value: str) -> tuple[str, str]:
    """Return background/text colours for compact audit cells.

    Status fields are process outcomes; semantic fields are domain flags.  This
    distinction follows the FRUTON workbook/audit policy: ligands are mild,
    metals and non-standard residues require review, and gaps are structural
    caution rather than automatic failure.
    """
    raw = value.strip().lower()
    field_l = field.strip().lower()

    if field_l in {"has_ligands", "ligands"}:
        return (_FRUTON_GREEN_SOFT, _IDIS_TEXT) if raw in {"yes", "true"} else ("#FFFFFF", _IDIS_MUTED)
    if field_l in {"has_metals", "metals", "has_nonstandard_residues", "nonstandard"}:
        return (_AUDIT_FAIL, _IDIS_TEXT) if raw in {"yes", "true"} else (_AUDIT_OK, _IDIS_TEXT)
    if field_l in {"has_gaps", "gaps"}:
        return (_AUDIT_WARN, _IDIS_TEXT) if raw in {"yes", "true"} else (_AUDIT_OK, _IDIS_TEXT)
    if field_l in {"overall_quality", "quality"}:
        if raw in {"good", "excellent", "acceptable", "pass", "passed"}:
            return (_AUDIT_OK, _IDIS_TEXT)
        if raw in {"poor", "bad", "failed", "fail"}:
            return (_AUDIT_FAIL, _IDIS_TEXT)
        if raw:
            return (_AUDIT_WARN, _IDIS_TEXT)
        return (_AUDIT_GREY, _IDIS_MUTED)

    if raw in {"success", "passed", "pass", "ready", "yes", "true"}:
        return (_AUDIT_OK, _IDIS_TEXT)
    if raw in {"warning", "required", "review", "caution"}:
        return (_AUDIT_WARN, _IDIS_TEXT)
    if raw in {"failed", "failure", "error", "no", "false"}:
        return (_AUDIT_FAIL, _IDIS_TEXT)
    if raw in {"skipped", "skip", "none", ""}:
        return (_AUDIT_GREY, _IDIS_MUTED)
    return (_IDIS_PANEL, _IDIS_TEXT)


# ---------------------------------------------------------------------------
# Quality summary
# ---------------------------------------------------------------------------

def _quality_summary(record: dict[str, str]) -> str:
    """Build one compact model-quality sentence for the title card."""
    favored = _record_text(record, "model_eval.rama_pct_favored")
    outlier = _record_text(record, "model_eval.rama_pct_outlier")
    clash = _record_text(record, "model_eval.clashscore")
    quality = _record_text(record, "model_eval.overall_quality")
    parts: list[str] = []
    if quality:
        parts.append(f"quality: {quality}")
    if favored:
        parts.append(f"Rama favored: {favored}%")
    if outlier:
        parts.append(f"outlier: {outlier}%")
    if clash:
        parts.append(f"clashscore: {clash}")
    return " · ".join(parts) if parts else "quality metrics not recorded"


# ---------------------------------------------------------------------------
# Page chrome (header/footer callback)
# ---------------------------------------------------------------------------

def _make_page_chrome(
    *,
    pdb_id: str,
    idis_logo: Path | None,
    fruton_logo: Path | None,
):
    """Return a ReportLab page callback with IDIS/FRUTON header and footer."""
    from reportlab.lib import colors
    from reportlab.lib.pagesizes import A4
    from reportlab.lib.units import cm

    page_w, page_h = A4

    def _draw(canvas, doc) -> None:  # type: ignore[no-untyped-def]
        canvas.saveState()
        left = doc.leftMargin
        right = page_w - doc.rightMargin
        top = page_h - 0.75 * cm

        if idis_logo and idis_logo.exists():
            canvas.drawImage(
                str(idis_logo),
                left,
                page_h - 1.45 * cm,
                width=2.1 * cm,
                height=0.85 * cm,
                preserveAspectRatio=True,
                mask="auto",
            )
        else:
            canvas.setFillColor(colors.HexColor(_IDIS_NAVY))
            canvas.setFont("Helvetica-Bold", 18)
            canvas.drawString(left, page_h - 1.20 * cm, "idis")

        canvas.setFillColor(colors.HexColor(_IDIS_NAVY))
        canvas.setFont("Helvetica-Bold", 8.5)
        canvas.drawCentredString(page_w / 2, top - 0.05 * cm, "FRUTON protein preparation audit")
        canvas.setFillColor(colors.HexColor(_IDIS_MUTED))
        canvas.setFont("Helvetica", 6.6)
        canvas.drawCentredString(page_w / 2, top - 0.35 * cm, f"per-protein report · {pdb_id}")

        if fruton_logo and fruton_logo.exists():
            canvas.drawImage(
                str(fruton_logo),
                right - 1.45 * cm,
                page_h - 1.45 * cm,
                width=1.45 * cm,
                height=0.95 * cm,
                preserveAspectRatio=True,
                mask="auto",
            )
        else:
            canvas.setFillColor(colors.HexColor(_FRUTON_RED))
            canvas.setFont("Helvetica-Bold", 9)
            canvas.drawRightString(right, page_h - 1.05 * cm, "FRUTON")

        y = page_h - 1.63 * cm
        span = right - left
        canvas.setLineWidth(1.0)
        canvas.setStrokeColor(colors.HexColor(_FRUTON_RED))
        canvas.line(left, y, left + 0.25 * span, y)
        canvas.setStrokeColor(colors.HexColor(_IDIS_GOLD))
        canvas.line(left + 0.26 * span, y, left + 0.49 * span, y)
        canvas.setStrokeColor(colors.HexColor(_IDIS_NAVY))
        canvas.line(left + 0.50 * span, y, right, y)

        footer_y = 0.95 * cm
        canvas.setStrokeColor(colors.HexColor(_IDIS_LINE))
        canvas.setLineWidth(0.45)
        canvas.line(left, footer_y + 0.22 * cm, right, footer_y + 0.22 * cm)
        canvas.setFillColor(colors.HexColor(_IDIS_MUTED))
        canvas.setFont("Helvetica", 6.5)
        canvas.drawString(left, footer_y, f"FRUTON - {_FRUTON_EXPANSION}")
        canvas.drawRightString(right, footer_y, f"page {doc.page}")
        canvas.restoreState()

    return _draw


# ---------------------------------------------------------------------------
# Section block
# ---------------------------------------------------------------------------

def _section_block(title: str, subtitle: str | None, styles: dict[str, Any], width: float):
    """Build a compact section title with a FRUTON red side rule."""
    from reportlab.lib import colors
    from reportlab.platypus import Paragraph, Table, TableStyle

    content = [Paragraph(_safe_pdf_text(title), styles["SectionHeading"])]
    if subtitle:
        content.append(Paragraph(_safe_pdf_text(subtitle), styles["SectionSubtitle"]))
    block = Table([[[*content]]], colWidths=[width])
    block.setStyle(TableStyle([
        ("LEFTPADDING", (0, 0), (-1, -1), 8),
        ("RIGHTPADDING", (0, 0), (-1, -1), 0),
        ("TOPPADDING", (0, 0), (-1, -1), 4),
        ("BOTTOMPADDING", (0, 0), (-1, -1), 4),
        ("LINEBEFORE", (0, 0), (0, 0), 3.0, colors.HexColor(_FRUTON_RED)),
    ]))
    return block


# ---------------------------------------------------------------------------
# Summary table (triage card)
# ---------------------------------------------------------------------------

def _build_summary_table(record: dict[str, str], styles: dict[str, Any], width: float):
    """Build the per-protein triage card shown on the first page."""
    from reportlab.lib import colors
    from reportlab.platypus import Paragraph, Table, TableStyle

    variant = _record_text(record, "prepared_structure.variant") or _record_text(record, "prepared_structure.final_type") or "not recorded"
    prep_status = _record_text(record, "prepared_structure.status") or "not recorded"
    models = _record_text(record, "available_models") or _record_text(record, "final_model_type") or "not recorded"
    n_gaps = _record_text(record, "n_gaps")
    gap_sizes = _record_text(record, "gap_sizes")
    gaps_value = "yes" if (n_gaps and n_gaps != "0") or _boolish_yes(_record_text(record, "has_gaps")) else "no"
    gaps_display = f"{n_gaps or '0'}" + (f" ({gap_sizes})" if gap_sizes else "")
    metal_value = _record_text(record, "has_metals") or "no"
    ligand_value = _record_text(record, "has_ligands") or "no"
    nonstd_value = _record_text(record, "has_nonstandard_residues") or "no"
    quality = _record_text(record, "model_eval.overall_quality") or "not recorded"

    rows = [
        ["prepared", prep_status, "variant", variant],
        ["models", models, "gaps", gaps_display],
        ["ligands", ligand_value, "metals", metal_value],
        ["nonstandard", nonstd_value, "quality", quality],
    ]
    data = []
    for label_a, value_a, label_b, value_b in rows:
        data.append([
            Paragraph(_safe_pdf_text(label_a.upper()), styles["MetricLabel"]),
            Paragraph(_safe_pdf_text(value_a), styles["MetricValue"]),
            Paragraph(_safe_pdf_text(label_b.upper()), styles["MetricLabel"]),
            Paragraph(_safe_pdf_text(value_b), styles["MetricValue"]),
        ])

    table = Table(data, colWidths=[2.35 * width / 14, 4.65 * width / 14, 2.35 * width / 14, 4.65 * width / 14])
    style_cmds: list[tuple[Any, ...]] = [
        ("BOX", (0, 0), (-1, -1), 0.75, colors.HexColor(_IDIS_LINE)),
        ("INNERGRID", (0, 0), (-1, -1), 0.35, colors.HexColor(_IDIS_LINE)),
        ("BACKGROUND", (0, 0), (-1, -1), colors.HexColor("#FFFFFF")),
        ("BACKGROUND", (0, 0), (0, -1), colors.HexColor(_IDIS_PANEL_2)),
        ("BACKGROUND", (2, 0), (2, -1), colors.HexColor(_IDIS_PANEL_2)),
        ("VALIGN", (0, 0), (-1, -1), "MIDDLE"),
        ("LEFTPADDING", (0, 0), (-1, -1), 6),
        ("RIGHTPADDING", (0, 0), (-1, -1), 6),
        ("TOPPADDING", (0, 0), (-1, -1), 5),
        ("BOTTOMPADDING", (0, 0), (-1, -1), 5),
    ]
    semantic_cells = {
        (1, 0): ("prepared", prep_status),
        (3, 1): ("gaps", gaps_value),
        (1, 2): ("has_ligands", ligand_value),
        (3, 2): ("has_metals", metal_value),
        (1, 3): ("has_nonstandard_residues", nonstd_value),
        (3, 3): ("quality", quality),
    }
    for cell, (field, value) in semantic_cells.items():
        bg, fg = _status_palette_for_value(field, value)
        style_cmds.append(("BACKGROUND", cell, cell, colors.HexColor(bg)))
        style_cmds.append(("TEXTCOLOR", cell, cell, colors.HexColor(fg)))
    table.setStyle(TableStyle(style_cmds))
    return table


# ---------------------------------------------------------------------------
# Narrative table
# ---------------------------------------------------------------------------

def _build_narrative_table(paragraphs: list[Any], width: float):
    """Wrap pipeline narrative paragraphs in a light audit panel."""
    from reportlab.lib import colors
    from reportlab.platypus import Table, TableStyle

    table = Table([[paragraphs]], colWidths=[width])
    table.setStyle(TableStyle([
        ("BACKGROUND", (0, 0), (-1, -1), colors.HexColor("#FFFFFF")),
        ("BOX", (0, 0), (-1, -1), 0.75, colors.HexColor(_IDIS_LINE)),
        ("LINEBEFORE", (0, 0), (-1, -1), 2.0, colors.HexColor(_IDIS_GOLD)),
        ("LEFTPADDING", (0, 0), (-1, -1), 10),
        ("RIGHTPADDING", (0, 0), (-1, -1), 10),
        ("TOPPADDING", (0, 0), (-1, -1), 8),
        ("BOTTOMPADDING", (0, 0), (-1, -1), 6),
    ]))
    return table


# ---------------------------------------------------------------------------
# Simple evidence table (two-column)
# ---------------------------------------------------------------------------

def _build_simple_evidence_table(
    rows: list[tuple[str, str]],
    styles: dict[str, Any],
    width: float,
    *,
    left_header: str = "Field",
    right_header: str = "Value",
):
    """Build a two-column evidence table for paths, parameters, and decisions."""
    from reportlab.lib import colors
    from reportlab.platypus import Paragraph, Table, TableStyle

    if not rows:
        rows = [("status", "no evidence recorded")]

    data = [[
        Paragraph(_safe_pdf_text(left_header), styles["EvidenceHeader"]),
        Paragraph(_safe_pdf_text(right_header), styles["EvidenceHeader"]),
    ]]
    for key, value in rows:
        data.append([
            Paragraph(_safe_pdf_text(key), styles["EvidenceCell"]),
            Paragraph(_safe_pdf_text(value), styles["EvidenceCell"]),
        ])

    table = Table(data, colWidths=[0.32 * width, 0.68 * width], repeatRows=1)
    table.setStyle(TableStyle([
        ("BACKGROUND", (0, 0), (-1, 0), colors.HexColor(_IDIS_PANEL_2)),
        ("TEXTCOLOR", (0, 0), (-1, 0), colors.HexColor(_IDIS_NAVY)),
        ("BOX", (0, 0), (-1, -1), 0.6, colors.HexColor(_IDIS_LINE)),
        ("INNERGRID", (0, 0), (-1, -1), 0.3, colors.HexColor(_IDIS_LINE)),
        ("VALIGN", (0, 0), (-1, -1), "TOP"),
        ("LEFTPADDING", (0, 0), (-1, -1), 5),
        ("RIGHTPADDING", (0, 0), (-1, -1), 5),
        ("TOPPADDING", (0, 0), (-1, -1), 4),
        ("BOTTOMPADDING", (0, 0), (-1, -1), 4),
    ]))
    return table


# ---------------------------------------------------------------------------
# Residue rename table (protonation step)
# ---------------------------------------------------------------------------

def _collect_residue_renames(
    input_pdb: Path,
    output_pdb: Path,
) -> list[dict[str, str]]:
    """Return per-residue rename records by diffing protonation input vs output.

    Both files are read and residues are matched by (chain, resseq, icode).
    A record is emitted only when the residue name changes. Returns an empty
    list when either file is missing or unreadable.

    Each record has keys: chain, res_num, icode, from_name, to_name, source.
    source is 'PROPKA' for HIP assignments, 'GROMACS' for others (e.g.
    disulfide CYS2, MSE->MET).
    """
    def _parse_residues(pdb_path: Path) -> dict[tuple[str, int, str], str]:
        res: dict[tuple[str, int, str], str] = {}
        try:
            with pdb_path.open("r", encoding="utf-8", errors="replace") as fh:
                for line in fh:
                    if not (line.startswith("ATOM  ") or line.startswith("HETATM")):
                        continue
                    try:
                        res_name = line[17:20].strip()
                        chain = line[21]
                        res_num = int(line[22:26])
                        icode = line[26].strip() if len(line) > 26 else ""
                        key = (chain, res_num, icode)
                        if key not in res:
                            res[key] = res_name
                    except (ValueError, IndexError):
                        pass
        except OSError:
            pass
        return res

    if not input_pdb.is_file() or not output_pdb.is_file():
        return []

    before = _parse_residues(input_pdb)
    after = _parse_residues(output_pdb)

    _PROPKA_TARGETS = {"HIP", "HIE", "HID"}
    renames: list[dict[str, str]] = []
    for key, from_name in sorted(before.items(), key=lambda kv: kv[0]):
        to_name = after.get(key)
        if to_name is None or to_name == from_name:
            continue
        chain, res_num, icode = key
        source = "PROPKA" if from_name == "HIS" and to_name in _PROPKA_TARGETS else "GROMACS"
        renames.append({
            "chain": chain,
            "res_num": str(res_num),
            "icode": icode or "",
            "from_name": from_name,
            "to_name": to_name,
            "source": source,
        })
    return renames


def _build_residue_rename_table(
    renames: list[dict[str, str]],
    styles: dict[str, Any],
    width: float,
) -> Any:
    """Build a table of protonation-step residue renamings."""
    from reportlab.lib import colors
    from reportlab.platypus import Paragraph, Table, TableStyle

    data: list[list[Any]] = [[
        Paragraph("Chain", styles["EvidenceHeader"]),
        Paragraph("Res.", styles["EvidenceHeader"]),
        Paragraph("Ins.", styles["EvidenceHeader"]),
        Paragraph("From", styles["EvidenceHeader"]),
        Paragraph("To", styles["EvidenceHeader"]),
        Paragraph("Driver", styles["EvidenceHeader"]),
        Paragraph("Note", styles["EvidenceHeader"]),
    ]]

    _notes = {
        "HIP": "doubly protonated His; net charge +1",
        "HIE": "neutral His, epsilon-tautomer (N-epsilon protonated); net charge 0",
        "HID": "neutral His, delta-tautomer (N-delta protonated); net charge 0",
        "CYS2": "disulfide-bridged Cys",
        "MET": "selenomethionine converted to methionine",
    }

    propka_rows: set[int] = set()
    for r in renames:
        row_idx = len(data)
        note = _notes.get(r["to_name"], "")
        data.append([
            Paragraph(_safe_pdf_text(r["chain"]), styles["EvidenceCell"]),
            Paragraph(_safe_pdf_text(r["res_num"]), styles["EvidenceCell"]),
            Paragraph(_safe_pdf_text(r["icode"] or "-"), styles["EvidenceCell"]),
            Paragraph(_safe_pdf_text(r["from_name"]), styles["EvidenceCell"]),
            Paragraph(_safe_pdf_text(r["to_name"]), styles["EvidenceCell"]),
            Paragraph(_safe_pdf_text(r["source"]), styles["EvidenceCell"]),
            Paragraph(_safe_pdf_text(note), styles["EvidenceCell"]),
        ])
        if r["source"] == "PROPKA":
            propka_rows.add(row_idx)

    col_w = width
    col_widths = [
        0.07 * col_w,   # Chain
        0.07 * col_w,   # Res.
        0.06 * col_w,   # Ins.
        0.10 * col_w,   # From
        0.10 * col_w,   # To
        0.12 * col_w,   # Driver
        0.48 * col_w,   # Note
    ]

    table = Table(data, colWidths=col_widths, repeatRows=1)
    style_cmds = [
        ("BACKGROUND", (0, 0), (-1, 0), colors.HexColor(_IDIS_PANEL_2)),
        ("TEXTCOLOR", (0, 0), (-1, 0), colors.HexColor(_IDIS_NAVY)),
        ("BOX", (0, 0), (-1, -1), 0.6, colors.HexColor(_IDIS_LINE)),
        ("INNERGRID", (0, 0), (-1, -1), 0.3, colors.HexColor(_IDIS_LINE)),
        ("VALIGN", (0, 0), (-1, -1), "TOP"),
        ("LEFTPADDING", (0, 0), (-1, -1), 5),
        ("RIGHTPADDING", (0, 0), (-1, -1), 5),
        ("TOPPADDING", (0, 0), (-1, -1), 4),
        ("BOTTOMPADDING", (0, 0), (-1, -1), 4),
    ]
    for row_idx in propka_rows:
        style_cmds.append(
            ("BACKGROUND", (0, row_idx), (-1, row_idx), colors.HexColor(_FRUTON_YELLOW_SOFT))
        )
    table.setStyle(TableStyle(style_cmds))
    return table


# ---------------------------------------------------------------------------
# Sequence evidence flowables
# ---------------------------------------------------------------------------

def _build_sequence_evidence_flowables(
    fasta_records: list[dict[str, str]],
    styles: dict[str, Any],
    width: float,
    protein_dir: Path,
) -> list[Any]:
    """Build flowables showing FASTA sources and wrapped sequence content."""
    from reportlab.lib import colors
    from reportlab.platypus import Paragraph, Spacer, Table, TableStyle
    from reportlab.lib.units import cm

    flowables: list[Any] = []
    summary_rows: list[list[Any]] = [[
        Paragraph("Source", styles["EvidenceHeader"]),
        Paragraph("Record", styles["EvidenceHeader"]),
        Paragraph("Length", styles["EvidenceHeader"]),
        Paragraph("X", styles["EvidenceHeader"]),
        Paragraph("File", styles["EvidenceHeader"]),
    ]]

    for idx, rec in enumerate(fasta_records, start=1):
        header = rec.get("header", "")
        short_header = header if len(header) <= 70 else header[:67] + "..."
        summary_rows.append([
            Paragraph(_safe_pdf_text(rec.get("source", "FASTA")), styles["EvidenceCell"]),
            Paragraph(_safe_pdf_text(short_header), styles["EvidenceCell"]),
            Paragraph(_safe_pdf_text(rec.get("length", "")), styles["EvidenceCell"]),
            Paragraph(_safe_pdf_text(rec.get("unknown_x", "0")), styles["EvidenceCell"]),
            Paragraph(_safe_pdf_text(_short_path_for_report(rec.get("path", ""), protein_dir)), styles["EvidenceCell"]),
        ])

    summary = Table(summary_rows, colWidths=[0.23 * width, 0.30 * width, 0.08 * width, 0.05 * width, 0.34 * width], repeatRows=1)
    summary.setStyle(TableStyle([
        ("BACKGROUND", (0, 0), (-1, 0), colors.HexColor(_IDIS_PANEL_2)),
        ("TEXTCOLOR", (0, 0), (-1, 0), colors.HexColor(_IDIS_NAVY)),
        ("BOX", (0, 0), (-1, -1), 0.6, colors.HexColor(_IDIS_LINE)),
        ("INNERGRID", (0, 0), (-1, -1), 0.25, colors.HexColor(_IDIS_LINE)),
        ("VALIGN", (0, 0), (-1, -1), "TOP"),
        ("LEFTPADDING", (0, 0), (-1, -1), 4),
        ("RIGHTPADDING", (0, 0), (-1, -1), 4),
        ("TOPPADDING", (0, 0), (-1, -1), 3.5),
        ("BOTTOMPADDING", (0, 0), (-1, -1), 3.5),
    ]))
    flowables.extend([summary, Spacer(1, 0.13 * cm)])

    for idx, rec in enumerate(fasta_records, start=1):
        wrapped, truncated = _wrap_sequence_for_pdf(rec.get("sequence", ""))
        trunc_note = " — sequence truncated for report length" if truncated else ""
        title = (
            f"FASTA {idx}: {rec.get('source', 'FASTA')} · "
            f"length {rec.get('length', '')} · X={rec.get('unknown_x', '0')}{trunc_note}"
        )
        seq_block = Table([[
            [
                Paragraph(_safe_pdf_text(title), styles["EvidenceHeader"]),
                Paragraph(wrapped, styles["MonoSequence"]),
            ]
        ]], colWidths=[width])
        seq_block.setStyle(TableStyle([
            ("BACKGROUND", (0, 0), (-1, -1), colors.HexColor(_IDIS_PANEL)),
            ("BOX", (0, 0), (-1, -1), 0.45, colors.HexColor(_IDIS_LINE)),
            ("LEFTPADDING", (0, 0), (-1, -1), 7),
            ("RIGHTPADDING", (0, 0), (-1, -1), 7),
            ("TOPPADDING", (0, 0), (-1, -1), 5),
            ("BOTTOMPADDING", (0, 0), (-1, -1), 5),
        ]))
        flowables.extend([seq_block, Spacer(1, 0.10 * cm)])

    return flowables
