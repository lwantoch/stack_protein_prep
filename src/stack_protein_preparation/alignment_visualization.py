# /home/grheco/repositorios/stack_protein_prep/src/stack_protein_preparation/alignment_visualization.py

"""
Render pairwise sequence alignments as publication-style PNG images.

Purpose
-------
- read a two-sequence FASTA alignment produced by MAFFT
- render the alignment as a multi-line block image
- show residue letters directly in the figure
- color exact matches, mismatches, and gaps distinctly

Default visual style
--------------------
- blue letters
- yellow background for exact matches
- grey background for residue mismatches
- red background for gaps / missing residues
- light page background

Expected input
--------------
A FASTA alignment with exactly two records, for example:
    >PDB|1ABC|ATOM|chain_A
    MKT--AI...
    >sp|P12345|PROT_HUMAN
    MKTEQAI...

Public API
----------
Main function expected by older pipeline code:
    alignment_to_image(alignment_fasta_path, output_png_path, ...)

Convenience alias:
    render_alignment_png(...)

Notes
-----
- this module does NOT compute alignments
- this module only visualizes already aligned sequences
"""

from __future__ import annotations

from dataclasses import dataclass
from math import ceil
from pathlib import Path
from typing import Iterable

from PIL import Image, ImageDraw, ImageFont

from stack_protein_preparation.pipeline_logging import (
    append_exception_traceback,
    append_key_value_block,
    append_log_text,
    capture_python_output_to_log,
    environment_snapshot_text,
    log_and_print_start,
    log_failure,
    log_success,
    print_light_table_box,
    print_mixed_box,
    print_step_end_box,
)

# ---------------------------------------------------------------------------
# Default visual style
# ---------------------------------------------------------------------------

PAGE_BACKGROUND = (245, 245, 245)
TEXT_BLUE = (48, 55, 97)
MATCH_BACKGROUND = (138, 115, 15)
MISMATCH_BACKGROUND = (136, 136, 136)
GAP_BACKGROUND = (178, 58, 58)

HEADER_TEXT = (40, 40, 40)
GUIDE_TEXT = (70, 70, 70)
BORDER_COLOR = (220, 220, 220)

DEFAULT_FONT_SIZE = 18
DEFAULT_HEADER_FONT_SIZE = 20
DEFAULT_LINE_SPACING = 16
DEFAULT_BLOCK_SPACING = 28
DEFAULT_CELL_PADDING_X = 6
DEFAULT_CELL_PADDING_Y = 4
DEFAULT_RESIDUES_PER_LINE = 60
DEFAULT_LABEL_WIDTH = 24
DEFAULT_LEFT_MARGIN = 40
DEFAULT_TOP_MARGIN = 36
DEFAULT_RIGHT_MARGIN = 40
DEFAULT_BOTTOM_MARGIN = 36

MODULE_NAME = "alignment_visualization"

# ---------------------------------------------------------------------------
# Small containers
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class FastaRecord:
    header: str
    sequence: str


@dataclass(frozen=True)
class AlignmentRenderConfig:
    residues_per_line: int = DEFAULT_RESIDUES_PER_LINE
    font_size: int = DEFAULT_FONT_SIZE
    header_font_size: int = DEFAULT_HEADER_FONT_SIZE
    line_spacing: int = DEFAULT_LINE_SPACING
    block_spacing: int = DEFAULT_BLOCK_SPACING
    cell_padding_x: int = DEFAULT_CELL_PADDING_X
    cell_padding_y: int = DEFAULT_CELL_PADDING_Y
    label_width: int = DEFAULT_LABEL_WIDTH
    left_margin: int = DEFAULT_LEFT_MARGIN
    top_margin: int = DEFAULT_TOP_MARGIN
    right_margin: int = DEFAULT_RIGHT_MARGIN
    bottom_margin: int = DEFAULT_BOTTOM_MARGIN


# ---------------------------------------------------------------------------
# Internal logging/path helpers
# ---------------------------------------------------------------------------


def _infer_pipeline_context(
    alignment_fasta_path: str | Path,
) -> tuple[Path | None, str | None]:
    """
    Try to infer:
    - protein_data_dir -> data/proteins
    - pdb_id          -> e.g. 1UOU

    Expected pipeline path shape:
        .../data/proteins/<PDBID>/fasta/alignments/<file>.fasta

    Returns
    -------
    tuple[Path | None, str | None]
        (protein_data_dir, pdb_id)
    """
    path = Path(alignment_fasta_path).resolve()

    parts = path.parts
    try:
        proteins_index = parts.index("proteins")
    except ValueError:
        return None, None

    if proteins_index + 1 >= len(parts):
        return None, None

    protein_data_dir = Path(*parts[: proteins_index + 1])
    pdb_id = parts[proteins_index + 1].strip().upper()

    if not pdb_id:
        return None, None

    return protein_data_dir, pdb_id


def _build_log_context_lines(
    *,
    alignment_fasta_path: Path,
    output_png_path: Path,
    residues_per_line: int,
    font_size: int,
    header_font_size: int,
    label_width: int,
) -> list[str]:
    return [
        f"alignment_fasta_path={alignment_fasta_path}",
        f"output_png_path={output_png_path}",
        f"residues_per_line={residues_per_line}",
        f"font_size={font_size}",
        f"header_font_size={header_font_size}",
        f"label_width={label_width}",
        *environment_snapshot_text(),
    ]


def _print_render_overview_box(
    *,
    alignment_fasta_path: Path,
    output_png_path: Path,
    residues_per_line: int,
    font_size: int,
    header_font_size: int,
    label_width: int,
) -> None:
    print_light_table_box(
        title="alignment_visualization · render plan",
        row_list=[
            ("Input FASTA", str(alignment_fasta_path)),
            ("Output PNG", str(output_png_path)),
            ("Residues/line", str(residues_per_line)),
            ("Font size", str(font_size)),
            ("Header font", str(header_font_size)),
            ("Label width", str(label_width)),
        ],
        left_header="Field",
        right_header="Value",
    )


# ---------------------------------------------------------------------------
# FASTA parsing
# ---------------------------------------------------------------------------


def read_fasta_records(fasta_path: str | Path) -> list[FastaRecord]:
    fasta_path = Path(fasta_path)

    if not fasta_path.is_file():
        raise FileNotFoundError(f"Alignment FASTA not found: {fasta_path}")

    text = fasta_path.read_text(encoding="utf-8").strip()
    if not text:
        raise ValueError(f"Alignment FASTA is empty: {fasta_path}")

    records: list[FastaRecord] = []
    current_header: str | None = None
    current_sequence_lines: list[str] = []

    for raw_line in text.splitlines():
        line = raw_line.strip()
        if not line:
            continue

        if line.startswith(">"):
            if current_header is not None:
                records.append(
                    FastaRecord(
                        header=current_header,
                        sequence="".join(current_sequence_lines).replace(" ", ""),
                    )
                )
            current_header = line[1:].strip()
            current_sequence_lines = []
            continue

        current_sequence_lines.append(line)

    if current_header is not None:
        records.append(
            FastaRecord(
                header=current_header,
                sequence="".join(current_sequence_lines).replace(" ", ""),
            )
        )

    if not records:
        raise ValueError(f"No FASTA records found in: {fasta_path}")

    return records


def read_pairwise_alignment(
    alignment_fasta_path: str | Path,
) -> tuple[FastaRecord, FastaRecord]:
    records = read_fasta_records(alignment_fasta_path)

    if len(records) != 2:
        raise ValueError(
            f"Expected exactly 2 FASTA records in pairwise alignment, found {len(records)}"
        )

    first_record, second_record = records

    if len(first_record.sequence) != len(second_record.sequence):
        raise ValueError(
            "Aligned sequences do not have the same length: "
            f"{len(first_record.sequence)} vs {len(second_record.sequence)}"
        )

    return first_record, second_record


# ---------------------------------------------------------------------------
# Text helpers
# ---------------------------------------------------------------------------


def _truncate_label(header: str, label_width: int) -> str:
    if len(header) <= label_width:
        return header
    if label_width <= 3:
        return header[:label_width]
    return header[: label_width - 3] + "..."


def _load_font(font_size: int) -> ImageFont.FreeTypeFont | ImageFont.ImageFont:
    preferred_paths = [
        "/usr/share/fonts/truetype/dejavu/DejaVuSansMono.ttf",
        "/usr/share/fonts/truetype/liberation2/LiberationMono-Regular.ttf",
    ]

    for font_path in preferred_paths:
        path = Path(font_path)
        if path.exists():
            try:
                return ImageFont.truetype(str(path), font_size)
            except Exception:
                pass

    return ImageFont.load_default()


def _measure_text(
    draw: ImageDraw.ImageDraw,
    text: str,
    font: ImageFont.FreeTypeFont | ImageFont.ImageFont,
) -> tuple[int, int]:
    bbox = draw.textbbox((0, 0), text, font=font)
    return bbox[2] - bbox[0], bbox[3] - bbox[1]


# ---------------------------------------------------------------------------
# Coloring logic
# ---------------------------------------------------------------------------


def _cell_background_color(char_a: str, char_b: str) -> tuple[int, int, int]:
    if char_a == "-" or char_b == "-":
        return GAP_BACKGROUND
    if char_a == char_b:
        return MATCH_BACKGROUND
    return MISMATCH_BACKGROUND


def _iter_alignment_blocks(
    seq_a: str,
    seq_b: str,
    residues_per_line: int,
) -> Iterable[tuple[int, int, str, str]]:
    total_length = len(seq_a)
    for start in range(0, total_length, residues_per_line):
        end = min(start + residues_per_line, total_length)
        yield start, end, seq_a[start:end], seq_b[start:end]


# ---------------------------------------------------------------------------
# Rendering
# ---------------------------------------------------------------------------


def render_alignment_png(
    alignment_fasta_path: str | Path,
    output_png_path: str | Path,
    *,
    residues_per_line: int = DEFAULT_RESIDUES_PER_LINE,
    font_size: int = DEFAULT_FONT_SIZE,
    header_font_size: int = DEFAULT_HEADER_FONT_SIZE,
    label_width: int = DEFAULT_LABEL_WIDTH,
) -> Path:
    """
    Render a two-sequence FASTA alignment as a PNG.

    Parameters
    ----------
    alignment_fasta_path
        Path to the aligned FASTA file with exactly two records.
    output_png_path
        Destination PNG path.
    residues_per_line
        Number of alignment columns shown per horizontal block.
    font_size
        Monospace font size for sequence letters.
    header_font_size
        Font size for title/header text.
    label_width
        Maximum visible width for row labels.
    """
    alignment_fasta_path = Path(alignment_fasta_path)
    output_png_path = Path(output_png_path)

    protein_data_dir, pdb_id = _infer_pipeline_context(alignment_fasta_path)
    log_path: Path | None = None

    if protein_data_dir is not None and pdb_id is not None:
        log_path = log_and_print_start(
            protein_data_dir=protein_data_dir,
            pdb_id=pdb_id,
            module_name=MODULE_NAME,
            extra_log_lines=_build_log_context_lines(
                alignment_fasta_path=alignment_fasta_path,
                output_png_path=output_png_path,
                residues_per_line=residues_per_line,
                font_size=font_size,
                header_font_size=header_font_size,
                label_width=label_width,
            ),
        )
    else:
        print_mixed_box(
            title=f"{MODULE_NAME} · START",
            line_list=[
                f"Input    : {alignment_fasta_path}",
                f"Output   : {output_png_path}",
            ],
        )

    _print_render_overview_box(
        alignment_fasta_path=alignment_fasta_path,
        output_png_path=output_png_path,
        residues_per_line=residues_per_line,
        font_size=font_size,
        header_font_size=header_font_size,
        label_width=label_width,
    )

    try:
        if log_path is not None:
            with capture_python_output_to_log(log_path):
                record_a, record_b = read_pairwise_alignment(alignment_fasta_path)
        else:
            record_a, record_b = read_pairwise_alignment(alignment_fasta_path)

        config = AlignmentRenderConfig(
            residues_per_line=residues_per_line,
            font_size=font_size,
            header_font_size=header_font_size,
            label_width=label_width,
        )

        mono_font = _load_font(config.font_size)
        header_font = _load_font(config.header_font_size)

        scratch_image = Image.new("RGB", (10, 10), PAGE_BACKGROUND)
        scratch_draw = ImageDraw.Draw(scratch_image)

        char_width, char_height = _measure_text(scratch_draw, "M", mono_font)
        title_height = _measure_text(scratch_draw, "Ag", header_font)[1]

        cell_width = char_width + 2 * config.cell_padding_x
        cell_height = char_height + 2 * config.cell_padding_y

        label_a = _truncate_label(record_a.header, config.label_width)
        label_b = _truncate_label(record_b.header, config.label_width)
        label_pixel_width = max(
            _measure_text(scratch_draw, label_a, mono_font)[0],
            _measure_text(scratch_draw, label_b, mono_font)[0],
            _measure_text(scratch_draw, "Pos", mono_font)[0],
        )

        n_blocks = ceil(len(record_a.sequence) / config.residues_per_line)

        block_height = (
            cell_height
            + config.line_spacing
            + cell_height
            + config.line_spacing
            + cell_height
        )

        image_width = (
            config.left_margin
            + label_pixel_width
            + 24
            + config.residues_per_line * cell_width
            + config.right_margin
        )

        image_height = (
            config.top_margin
            + title_height
            + 20
            + n_blocks * block_height
            + max(0, n_blocks - 1) * config.block_spacing
            + config.bottom_margin
        )

        if log_path is not None:
            append_key_value_block(
                log_path=log_path,
                title="RENDER METRICS",
                key_value_pairs=[
                    ("record_a_header", record_a.header),
                    ("record_b_header", record_b.header),
                    ("alignment_length", str(len(record_a.sequence))),
                    ("n_blocks", str(n_blocks)),
                    ("cell_width", str(cell_width)),
                    ("cell_height", str(cell_height)),
                    ("image_width", str(image_width)),
                    ("image_height", str(image_height)),
                ],
            )

        print_light_table_box(
            title="alignment_visualization · parsed alignment",
            row_list=[
                ("Sequence A header", record_a.header),
                ("Sequence B header", record_b.header),
                ("Alignment length", str(len(record_a.sequence))),
                ("Blocks", str(n_blocks)),
                ("Image size", f"{image_width} x {image_height}"),
            ],
            left_header="Item",
            right_header="Value",
        )

        image = Image.new("RGB", (image_width, image_height), PAGE_BACKGROUND)
        draw = ImageDraw.Draw(image)

        title = alignment_fasta_path.name
        draw.text(
            (config.left_margin, config.top_margin),
            title,
            fill=HEADER_TEXT,
            font=header_font,
        )

        current_y = config.top_margin + title_height + 20
        seq_x0 = config.left_margin + label_pixel_width + 24

        for start, end, block_a, block_b in _iter_alignment_blocks(
            record_a.sequence,
            record_b.sequence,
            config.residues_per_line,
        ):
            draw.text(
                (config.left_margin, current_y + config.cell_padding_y),
                "Pos",
                fill=GUIDE_TEXT,
                font=mono_font,
            )

            for i in range(len(block_a)):
                absolute_index = start + i + 1
                x = seq_x0 + i * cell_width

                if i == 0 or absolute_index % 10 == 0:
                    draw.text(
                        (x + 1, current_y + config.cell_padding_y),
                        str(absolute_index),
                        fill=GUIDE_TEXT,
                        font=mono_font,
                    )

            row_a_y = current_y + cell_height + config.line_spacing
            draw.text(
                (config.left_margin, row_a_y + config.cell_padding_y),
                label_a,
                fill=HEADER_TEXT,
                font=mono_font,
            )

            row_b_y = row_a_y + cell_height + config.line_spacing
            draw.text(
                (config.left_margin, row_b_y + config.cell_padding_y),
                label_b,
                fill=HEADER_TEXT,
                font=mono_font,
            )

            for i, (char_a, char_b) in enumerate(zip(block_a, block_b)):
                x = seq_x0 + i * cell_width
                cell_color = _cell_background_color(char_a, char_b)

                draw.rectangle(
                    [x, row_a_y, x + cell_width, row_a_y + cell_height],
                    fill=cell_color,
                    outline=BORDER_COLOR,
                )
                draw.text(
                    (x + config.cell_padding_x, row_a_y + config.cell_padding_y),
                    char_a,
                    fill=TEXT_BLUE,
                    font=mono_font,
                )

                draw.rectangle(
                    [x, row_b_y, x + cell_width, row_b_y + cell_height],
                    fill=cell_color,
                    outline=BORDER_COLOR,
                )
                draw.text(
                    (x + config.cell_padding_x, row_b_y + config.cell_padding_y),
                    char_b,
                    fill=TEXT_BLUE,
                    font=mono_font,
                )

            current_y += block_height + config.block_spacing

        output_png_path.parent.mkdir(parents=True, exist_ok=True)
        image.save(output_png_path)

        if log_path is not None:
            append_log_text(log_path, f"[OUTPUT PNG] {output_png_path}")
            log_success(log_path, f"rendered alignment PNG -> {output_png_path}")

        if pdb_id is not None:
            print_step_end_box(
                module_name=MODULE_NAME,
                pdb_id=pdb_id,
                status="success",
                detail_lines=[
                    f"Alignment length : {len(record_a.sequence)}",
                    f"Blocks           : {n_blocks}",
                    f"PNG              : {output_png_path}",
                ],
            )
        else:
            print_mixed_box(
                title=f"{MODULE_NAME} · END",
                line_list=[
                    "Status           : success",
                    f"Alignment length : {len(record_a.sequence)}",
                    f"Blocks           : {n_blocks}",
                    f"PNG              : {output_png_path}",
                ],
            )

        return output_png_path

    except Exception as error:
        if log_path is not None:
            log_failure(log_path, "alignment rendering failed", error)
            append_exception_traceback(log_path, error)

        if pdb_id is not None:
            print_step_end_box(
                module_name=MODULE_NAME,
                pdb_id=pdb_id,
                status="failed",
                detail_lines=[f"Error   : {error!r}"],
            )
        else:
            print_mixed_box(
                title=f"{MODULE_NAME} · END",
                line_list=[
                    "Status  : failed",
                    f"Error   : {error!r}",
                ],
            )

        raise


# ---------------------------------------------------------------------------
# Backward-compatible public name expected by older code
# ---------------------------------------------------------------------------


def alignment_to_image(
    alignment_fasta_path: str | Path,
    output_png_path: str | Path,
    *,
    residues_per_line: int = DEFAULT_RESIDUES_PER_LINE,
    font_size: int = DEFAULT_FONT_SIZE,
    header_font_size: int = DEFAULT_HEADER_FONT_SIZE,
    label_width: int = DEFAULT_LABEL_WIDTH,
) -> Path:
    """
    Backward-compatible wrapper used by older pipeline code.
    """
    return render_alignment_png(
        alignment_fasta_path=alignment_fasta_path,
        output_png_path=output_png_path,
        residues_per_line=residues_per_line,
        font_size=font_size,
        header_font_size=header_font_size,
        label_width=label_width,
    )
