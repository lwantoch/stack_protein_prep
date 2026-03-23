from __future__ import annotations

import math
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Sequence

import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle


@dataclass
class AlignmentEntry:
    header: str
    sequence: str


def _read_fasta_alignment(fasta_path: str | Path) -> list[AlignmentEntry]:
    """
    Read an aligned FASTA file and return all entries in order.

    Parameters
    ----------
    fasta_path
        Path to aligned FASTA file.

    Returns
    -------
    list[AlignmentEntry]
        Parsed FASTA entries.

    Raises
    ------
    FileNotFoundError
        If the FASTA file does not exist.
    ValueError
        If the FASTA file is empty or contains no entries.
    """
    fasta_path = Path(fasta_path)

    if not fasta_path.exists():
        raise FileNotFoundError(f"Alignment FASTA not found: {fasta_path}")

    entries: list[AlignmentEntry] = []
    current_header: str | None = None
    current_seq: list[str] = []

    with fasta_path.open("r", encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue

            if line.startswith(">"):
                if current_header is not None:
                    entries.append(
                        AlignmentEntry(
                            header=current_header,
                            sequence="".join(current_seq),
                        )
                    )
                current_header = line[1:].strip()
                current_seq = []
            else:
                current_seq.append(line)

    if current_header is not None:
        entries.append(
            AlignmentEntry(
                header=current_header,
                sequence="".join(current_seq),
            )
        )

    if not entries:
        raise ValueError(f"No FASTA entries found in {fasta_path}")

    return entries


def _slugify_header(header: str) -> str:
    slug = re.sub(r"[^A-Za-z0-9._-]+", "_", header.strip())
    slug = slug.strip("_.")
    return slug or "sequence"


def _looks_like_uniprot_header(header: str) -> bool:
    header_lower = header.lower()
    return (
        "uniprot" in header_lower
        or header_lower.startswith("sp|")
        or header_lower.startswith("tr|")
        or "|uniprot|" in header_lower
    )


def _extract_chain_label(header: str) -> str:
    """
    Try to infer a compact chain label from FASTA header.

    Examples
    --------
    '1A5H_CHAIN_A' -> 'A'
    'chain A' -> 'A'
    """
    patterns = [
        r"\bCHAIN[_\s-]?([A-Za-z0-9])\b",
        r"\bchain[_\s-]?([A-Za-z0-9])\b",
        r"_([A-Za-z0-9])$",
    ]

    for pattern in patterns:
        match = re.search(pattern, header)
        if match:
            return match.group(1)

    return _slugify_header(header)


def _pair_alignment_entries(
    entries: list[AlignmentEntry],
) -> list[tuple[AlignmentEntry, AlignmentEntry, str]]:
    """
    Split FASTA entries into pairwise plots.

    Returns
    -------
    list[tuple[AlignmentEntry, AlignmentEntry, str]]
        Tuples of:
        (entry_a, entry_b, output_suffix)

    Rules
    -----
    - exactly 2 entries:
        render one pair
    - exactly 1 UniProt entry + n non-UniProt entries:
        render n pairs (chain vs UniProt)
    - even number of entries without UniProt detection:
        render in 0/1, 2/3, 4/5 order
    """
    if len(entries) == 2:
        suffix = _slugify_header(entries[0].header)
        return [(entries[0], entries[1], suffix)]

    uniprot_entries = [
        entry for entry in entries if _looks_like_uniprot_header(entry.header)
    ]

    if len(uniprot_entries) == 1:
        uniprot_entry = uniprot_entries[0]
        other_entries = [entry for entry in entries if entry is not uniprot_entry]

        if not other_entries:
            raise ValueError(
                "Found UniProt entry, but no matching chain entries to render."
            )

        pairs: list[tuple[AlignmentEntry, AlignmentEntry, str]] = []
        for other in other_entries:
            suffix = f"chain_{_extract_chain_label(other.header)}"
            pairs.append((other, uniprot_entry, suffix))
        return pairs

    if len(entries) % 2 == 0:
        pairs: list[tuple[AlignmentEntry, AlignmentEntry, str]] = []
        for i in range(0, len(entries), 2):
            entry_a = entries[i]
            entry_b = entries[i + 1]
            suffix = _slugify_header(entry_a.header)
            pairs.append((entry_a, entry_b, suffix))
        return pairs

    headers = [entry.header for entry in entries]
    raise ValueError(
        "Could not determine how to split alignment entries into pairwise plots. "
        f"Found {len(entries)} entries: {headers}"
    )


def _truncate_header(header: str, max_len: int = 60) -> str:
    if len(header) <= max_len:
        return header
    return f"{header[: max_len - 3]}..."


def _residue_bg_color(res_a: str, res_b: str) -> str:
    """
    Background coloring rule.

    - exact match (not gap): yellow
    - any gap: red
    - mismatch: light grey
    """
    if res_a == "-" or res_b == "-":
        return "#f4a6a6"  # soft red
    if res_a == res_b:
        return "#fff3a3"  # soft yellow
    return "#d9d9d9"  # light grey


def _draw_sequence_block(
    ax,
    seq_a: str,
    seq_b: str,
    header_a: str,
    header_b: str,
    start_idx: int,
    end_idx: int,
    y_top: float,
    row_height: float,
    left_margin: float,
    x_offset: float = 0.0,
) -> None:
    """
    Draw one alignment block covering residues [start_idx:end_idx].
    """
    block_a = seq_a[start_idx:end_idx]
    block_b = seq_b[start_idx:end_idx]

    label_x = left_margin - 1.0
    seq_start_x = left_margin + x_offset

    # Left labels
    ax.text(
        label_x,
        y_top,
        _truncate_header(header_a, 28),
        ha="right",
        va="center",
        fontsize=8,
        fontfamily="monospace",
        color="black",
    )
    ax.text(
        label_x,
        y_top - row_height,
        _truncate_header(header_b, 28),
        ha="right",
        va="center",
        fontsize=8,
        fontfamily="monospace",
        color="black",
    )

    # Residue positions
    ax.text(
        label_x - 0.25,
        y_top,
        str(start_idx + 1),
        ha="right",
        va="center",
        fontsize=7,
        fontfamily="monospace",
        color="black",
    )
    ax.text(
        label_x - 0.25,
        y_top - row_height,
        str(start_idx + 1),
        ha="right",
        va="center",
        fontsize=7,
        fontfamily="monospace",
        color="black",
    )

    for i, (res_a, res_b) in enumerate(zip(block_a, block_b)):
        x = seq_start_x + i
        bg = _residue_bg_color(res_a, res_b)

        # top row cell
        ax.add_patch(
            Rectangle(
                (x - 0.5, y_top - 0.42),
                width=1.0,
                height=0.84,
                facecolor=bg,
                edgecolor="none",
            )
        )

        # bottom row cell
        ax.add_patch(
            Rectangle(
                (x - 0.5, y_top - row_height - 0.42),
                width=1.0,
                height=0.84,
                facecolor=bg,
                edgecolor="none",
            )
        )

        # residues
        ax.text(
            x,
            y_top,
            res_a,
            ha="center",
            va="center",
            fontsize=8,
            fontfamily="monospace",
            color="#1f4ed8",  # blue letters
        )
        ax.text(
            x,
            y_top - row_height,
            res_b,
            ha="center",
            va="center",
            fontsize=8,
            fontfamily="monospace",
            color="#1f4ed8",  # blue letters
        )

        # match marker row
        if res_a == res_b and res_a != "-":
            marker = "|"
        elif res_a == "-" or res_b == "-":
            marker = " "
        else:
            marker = "."

        ax.text(
            x,
            y_top - (row_height / 2.0),
            marker,
            ha="center",
            va="center",
            fontsize=7,
            fontfamily="monospace",
            color="black",
        )

    # end positions
    ax.text(
        seq_start_x + len(block_a) + 0.25,
        y_top,
        str(end_idx),
        ha="left",
        va="center",
        fontsize=7,
        fontfamily="monospace",
        color="black",
    )
    ax.text(
        seq_start_x + len(block_b) + 0.25,
        y_top - row_height,
        str(end_idx),
        ha="left",
        va="center",
        fontsize=7,
        fontfamily="monospace",
        color="black",
    )


def _render_alignment_core(
    header_a: str,
    seq_a: str,
    header_b: str,
    seq_b: str,
    output_png_path: str | Path,
    block_size: int = 80,
    dpi: int = 300,
) -> Path:
    """
    Render a pairwise alignment as a PNG.

    Visual style:
    - blue letters
    - yellow background for exact matches
    - red background for gaps
    - grey background for mismatches
    """
    output_png_path = Path(output_png_path)

    if len(seq_a) != len(seq_b):
        raise ValueError(
            f"Aligned sequences must have the same length, got {len(seq_a)} and {len(seq_b)} "
            f"for '{header_a}' vs '{header_b}'."
        )

    if block_size <= 0:
        raise ValueError(f"block_size must be > 0, got {block_size}")

    n_cols = len(seq_a)
    n_blocks = math.ceil(n_cols / block_size)

    # layout parameters in data coordinates
    row_height = 1.2
    block_vertical_gap = 1.4
    top_margin = 1.5
    bottom_margin = 1.0
    left_margin = 8.5
    right_margin = 2.5

    total_height = (
        top_margin + bottom_margin + n_blocks * (2 * row_height + block_vertical_gap)
    )
    total_width = left_margin + right_margin + block_size + 3

    # convert rough data-space to figure inches
    fig_width = max(12, total_width * 0.16)
    fig_height = max(4, total_height * 0.28)

    fig, ax = plt.subplots(figsize=(fig_width, fig_height), dpi=dpi)
    ax.set_xlim(0, total_width)
    ax.set_ylim(0, total_height)
    ax.axis("off")

    title = f"{_truncate_header(header_a, 50)} vs {_truncate_header(header_b, 50)}"
    ax.text(
        left_margin,
        total_height - 0.7,
        title,
        ha="left",
        va="center",
        fontsize=11,
        fontweight="bold",
        color="black",
    )

    for block_idx in range(n_blocks):
        start_idx = block_idx * block_size
        end_idx = min((block_idx + 1) * block_size, n_cols)

        y_top = (
            total_height
            - top_margin
            - block_idx * (2 * row_height + block_vertical_gap)
        )

        _draw_sequence_block(
            ax=ax,
            seq_a=seq_a,
            seq_b=seq_b,
            header_a=header_a,
            header_b=header_b,
            start_idx=start_idx,
            end_idx=end_idx,
            y_top=y_top,
            row_height=row_height,
            left_margin=left_margin,
        )

    output_png_path.parent.mkdir(parents=True, exist_ok=True)
    plt.tight_layout()
    fig.savefig(output_png_path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)

    return output_png_path


def render_pairwise_alignment_entries_png(
    entry_a: AlignmentEntry,
    entry_b: AlignmentEntry,
    output_png_path: str | Path,
    block_size: int = 80,
    dpi: int = 300,
) -> Path:
    """
    Render two already-parsed alignment entries into a PNG.
    """
    return _render_alignment_core(
        header_a=entry_a.header,
        seq_a=entry_a.sequence,
        header_b=entry_b.header,
        seq_b=entry_b.sequence,
        output_png_path=output_png_path,
        block_size=block_size,
        dpi=dpi,
    )


def render_pairwise_alignment_png(
    fasta_path: str | Path,
    output_png_path: str | Path,
    block_size: int = 80,
    dpi: int = 300,
) -> Path:
    """
    Render a FASTA alignment containing exactly 2 entries into a PNG.
    """
    entries = _read_fasta_alignment(fasta_path)

    if len(entries) != 2:
        raise ValueError(
            f"Expected exactly 2 aligned FASTA entries in {fasta_path}, got {len(entries)}."
        )

    return render_pairwise_alignment_entries_png(
        entry_a=entries[0],
        entry_b=entries[1],
        output_png_path=output_png_path,
        block_size=block_size,
        dpi=dpi,
    )


def alignment_to_image(
    fasta_path: str | Path,
    output_png_path: str | Path,
    block_size: int = 80,
    dpi: int = 300,
) -> list[Path]:
    """
    Convert an aligned FASTA file into one or more PNG images.

    Behavior
    --------
    - 2 entries:
        exactly 1 PNG
    - 1 UniProt + n chains:
        n PNGs, one chain vs UniProt
    - even number of entries without UniProt header:
        pairwise PNGs in order

    Returns
    -------
    list[Path]
        Paths to generated PNG files.
    """
    fasta_path = Path(fasta_path)
    output_png_path = Path(output_png_path)

    entries = _read_fasta_alignment(fasta_path)
    pairs = _pair_alignment_entries(entries)

    generated_paths: list[Path] = []

    if len(pairs) == 1:
        entry_a, entry_b, _suffix = pairs[0]
        png_path = render_pairwise_alignment_entries_png(
            entry_a=entry_a,
            entry_b=entry_b,
            output_png_path=output_png_path,
            block_size=block_size,
            dpi=dpi,
        )
        generated_paths.append(png_path)
        return generated_paths

    stem = output_png_path.stem
    suffix = output_png_path.suffix or ".png"
    parent = output_png_path.parent

    for entry_a, entry_b, pair_suffix in pairs:
        pair_output = parent / f"{stem}_{pair_suffix}{suffix}"
        png_path = render_pairwise_alignment_entries_png(
            entry_a=entry_a,
            entry_b=entry_b,
            output_png_path=pair_output,
            block_size=block_size,
            dpi=dpi,
        )
        generated_paths.append(png_path)

    return generated_paths


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description="Render aligned FASTA file into one or more alignment PNGs."
    )
    parser.add_argument("fasta_path", type=Path, help="Path to aligned FASTA file")
    parser.add_argument("output_png_path", type=Path, help="Base output PNG path")
    parser.add_argument(
        "--block-size",
        type=int,
        default=80,
        help="Number of aligned residues per line block",
    )
    parser.add_argument(
        "--dpi",
        type=int,
        default=300,
        help="Output DPI",
    )

    args = parser.parse_args()

    generated = alignment_to_image(
        fasta_path=args.fasta_path,
        output_png_path=args.output_png_path,
        block_size=args.block_size,
        dpi=args.dpi,
    )

    for path in generated:
        print(path)
