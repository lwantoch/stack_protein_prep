# /home/grheco/repositorios/stack_protein_prep/src/stack_protein_preparation/fragment_split.py

"""
Split protein-only PDB files into physically disconnected fragments for the
gaps-protonation route.

Purpose
-------
The gaps route should NOT send one artificially capped multi-fragment structure
directly into pdb2pqr. Instead, the protein should first be split into true
connected backbone fragments, each fragment protonated separately, and only
after that merged again.

This module is responsible for the first step only:
- read a protein-only PDB
- detect backbone disconnections inside each chain
- write one fragment PDB per connected fragment
- write a manifest JSON describing the split
- write a per-run log file

Important design choice
-----------------------
This module does NOT rely on UniProt alignments or filler output.
It uses the actual input PDB geometry and chain continuity.

A boundary between two consecutive residues in the same chain is introduced when
the residues are NOT peptide-connected, i.e.:
- left residue has no backbone C atom, or
- right residue has no backbone N atom, or
- C(left) ... N(right) distance is greater than the cutoff

This means:
- numbering artifacts with a still-connected backbone do NOT force a split
- true structural breaks with bad/missing peptide geometry DO force a split

Recommended place in pipeline
-----------------------------
For the "gaps" variant, use this BEFORE protonation:

    original protein-only PDB
        -> split_gap_variant_structure(...)
        -> protonate each fragment separately
        -> merge protonated fragments
        -> continue with amber_renaming etc.

Expected input
--------------
- protein-only PDB
- typically:
    data/proteins/<PDB_ID>/components/<PDB_ID>_protein.pdb

Default output layout
---------------------
Given:
    data/proteins/<PDB_ID>/components/<PDB_ID>_protein.pdb

This module writes by default:
    data/proteins/<PDB_ID>/components/fragments/gaps/
        <PDB_ID>_chain_A_fragment_01.pdb
        <PDB_ID>_chain_A_fragment_02.pdb
        ...
        fragment_manifest.json

Logging
-------
Default log path:
    data/proteins/logs/<PDB_ID>/fragment_split.gaps.log
"""

from __future__ import annotations

import argparse
import json
import math
from dataclasses import dataclass, field
from datetime import datetime
from pathlib import Path
from typing import Any

MODULE_NAME = "fragment_split"
DEFAULT_VARIANT_LABEL = "gaps"
DEFAULT_CN_DISTANCE_CUTOFF = 1.8

__all__ = [
    "split_protein_for_gap_protonation",
    "split_gap_variant_structure",
]


# ---------------------------------------------------------------------------
# Terminal helpers
# ---------------------------------------------------------------------------


def _screen_header(title: str, lines: list[str]) -> None:
    width = max([len(title), *[len(line) for line in lines]]) if lines else len(title)
    print(f"╔{'═' * (width + 2)}╗")
    print(f"║ {title.ljust(width)} ║")
    if lines:
        print(f"╠{'═' * (width + 2)}╣")
        for line in lines:
            print(f"║ {line.ljust(width)} ║")
    print(f"╚{'═' * (width + 2)}╝")


def _screen_sub(message: str) -> None:
    print(f"╟─ {message}")


def _screen_result(message: str) -> None:
    print(f"╚═ {message}")


# ---------------------------------------------------------------------------
# Logging helpers
# ---------------------------------------------------------------------------


def _timestamp() -> str:
    return datetime.now().strftime("%Y-%m-%d %H:%M:%S")


def _append_log(log_path: Path, title: str, lines: list[str]) -> None:
    log_path.parent.mkdir(parents=True, exist_ok=True)
    with log_path.open("a", encoding="utf-8") as handle:
        handle.write("═" * 100 + "\n")
        handle.write(f"[{_timestamp()}] {title}\n")
        handle.write("─" * 100 + "\n")
        for line in lines:
            handle.write(f"{line}\n")
        handle.write("\n")


# ---------------------------------------------------------------------------
# Data structures
# ---------------------------------------------------------------------------


@dataclass
class AtomRecord:
    atom_name: str
    altloc: str
    x: float
    y: float
    z: float
    line: str


@dataclass
class ResidueRecord:
    chain_id: str
    resseq: int
    icode: str
    resname: str
    lines: list[str] = field(default_factory=list)
    atoms_by_name: dict[str, list[AtomRecord]] = field(default_factory=dict)

    @property
    def chain_label(self) -> str:
        return self.chain_id if self.chain_id else "(blank)"

    @property
    def residue_label(self) -> str:
        insertion = self.icode if self.icode else ""
        return f"{self.resname} {self.chain_label} {self.resseq}{insertion}"

    def add_atom(self, atom: AtomRecord) -> None:
        self.lines.append(atom.line)
        self.atoms_by_name.setdefault(atom.atom_name, []).append(atom)

    def get_preferred_atom(self, atom_name: str) -> AtomRecord | None:
        atom_list = self.atoms_by_name.get(atom_name, [])
        if not atom_list:
            return None

        altloc_priority = {"": 0, "A": 1, "1": 2}

        def sort_key(atom: AtomRecord) -> tuple[int, str]:
            return (altloc_priority.get(atom.altloc, 99), atom.altloc)

        return sorted(atom_list, key=sort_key)[0]

    def atom_count(self) -> int:
        return len(self.lines)


@dataclass
class BoundaryRecord:
    chain_id: str
    left_residue_label: str
    right_residue_label: str
    left_resseq: int
    right_resseq: int
    residue_jump: int
    cn_distance: float | None
    reason: str


@dataclass
class FragmentRecord:
    chain_id: str
    fragment_index_global: int
    fragment_index_in_chain: int
    residues: list[ResidueRecord]
    output_pdb_path: Path

    @property
    def chain_label(self) -> str:
        return self.chain_id if self.chain_id else "blank"

    @property
    def start_resseq(self) -> int:
        return self.residues[0].resseq

    @property
    def end_resseq(self) -> int:
        return self.residues[-1].resseq

    @property
    def start_resname(self) -> str:
        return self.residues[0].resname

    @property
    def end_resname(self) -> str:
        return self.residues[-1].resname

    @property
    def residue_count(self) -> int:
        return len(self.residues)

    @property
    def atom_count(self) -> int:
        return sum(residue.atom_count() for residue in self.residues)


# ---------------------------------------------------------------------------
# Path helpers
# ---------------------------------------------------------------------------


def _infer_pdb_id_from_input_path(input_pdb_path: Path) -> str:
    stem = input_pdb_path.stem
    if stem.endswith("_protein"):
        return stem.removesuffix("_protein").upper()
    if "_protein_" in stem:
        return stem.split("_protein_", 1)[0].upper()
    return stem.split("_", 1)[0].upper()


def _infer_protein_dir(input_pdb_path: Path) -> Path:
    if input_pdb_path.parent.name == "components":
        return input_pdb_path.parent.parent
    if input_pdb_path.parent.parent.name == "components":
        return input_pdb_path.parent.parent.parent
    return input_pdb_path.parent


def _default_output_dir(
    input_pdb_path: Path,
    variant_label: str,
) -> Path:
    return input_pdb_path.parent / "fragments" / variant_label


def _default_log_path(
    input_pdb_path: Path,
    pdb_id: str,
    variant_label: str,
) -> Path:
    protein_dir = _infer_protein_dir(input_pdb_path)
    logs_dir = protein_dir.parent / "logs" / pdb_id
    return logs_dir / f"{MODULE_NAME}.{variant_label}.log"


def _default_input_pdb_path_for_variant(
    protein_dir: Path,
    pdb_id: str,
    variant_label: str,
) -> Path:
    components_dir = protein_dir / "components"

    if variant_label == "gaps":
        return components_dir / f"{pdb_id}_protein.pdb"

    return components_dir / f"{pdb_id}_protein_{variant_label}.pdb"


# ---------------------------------------------------------------------------
# PDB parsing helpers
# ---------------------------------------------------------------------------


def _parse_float(field: str) -> float:
    return float(field.strip())


def _parse_atom_line(
    line: str,
) -> tuple[str, int, str, str, str, str, float, float, float]:
    chain_id = line[21].strip()
    resseq = int(line[22:26].strip())
    icode = line[26].strip()
    resname = line[17:20].strip().upper()
    atom_name = line[12:16].strip().upper()
    altloc = line[16].strip()
    x = _parse_float(line[30:38])
    y = _parse_float(line[38:46])
    z = _parse_float(line[46:54])
    return chain_id, resseq, icode, resname, atom_name, altloc, x, y, z


def read_residue_records_from_protein_pdb(input_pdb_path: Path) -> list[ResidueRecord]:
    residue_list: list[ResidueRecord] = []

    current_key: tuple[str, int, str] | None = None
    current_residue: ResidueRecord | None = None

    with input_pdb_path.open("r", encoding="utf-8") as handle:
        for line in handle:
            if not line.startswith("ATOM"):
                continue

            (
                chain_id,
                resseq,
                icode,
                resname,
                atom_name,
                altloc,
                x,
                y,
                z,
            ) = _parse_atom_line(line)

            residue_key = (chain_id, resseq, icode)

            if current_key != residue_key:
                if current_residue is not None:
                    residue_list.append(current_residue)

                current_residue = ResidueRecord(
                    chain_id=chain_id,
                    resseq=resseq,
                    icode=icode,
                    resname=resname,
                )
                current_key = residue_key

            assert current_residue is not None

            current_residue.add_atom(
                AtomRecord(
                    atom_name=atom_name,
                    altloc=altloc,
                    x=x,
                    y=y,
                    z=z,
                    line=line,
                )
            )

    if current_residue is not None:
        residue_list.append(current_residue)

    return residue_list


# ---------------------------------------------------------------------------
# Geometry helpers
# ---------------------------------------------------------------------------


def _distance(atom_a: AtomRecord, atom_b: AtomRecord) -> float:
    return math.sqrt(
        (atom_a.x - atom_b.x) ** 2
        + (atom_a.y - atom_b.y) ** 2
        + (atom_a.z - atom_b.z) ** 2
    )


def assess_peptide_connection(
    left_residue: ResidueRecord,
    right_residue: ResidueRecord,
    cn_distance_cutoff: float,
) -> tuple[bool, BoundaryRecord | None]:
    left_c = left_residue.get_preferred_atom("C")
    right_n = right_residue.get_preferred_atom("N")
    residue_jump = right_residue.resseq - left_residue.resseq

    if left_c is None:
        return (
            False,
            BoundaryRecord(
                chain_id=left_residue.chain_id,
                left_residue_label=left_residue.residue_label,
                right_residue_label=right_residue.residue_label,
                left_resseq=left_residue.resseq,
                right_resseq=right_residue.resseq,
                residue_jump=residue_jump,
                cn_distance=None,
                reason="missing_left_backbone_C",
            ),
        )

    if right_n is None:
        return (
            False,
            BoundaryRecord(
                chain_id=left_residue.chain_id,
                left_residue_label=left_residue.residue_label,
                right_residue_label=right_residue.residue_label,
                left_resseq=left_residue.resseq,
                right_resseq=right_residue.resseq,
                residue_jump=residue_jump,
                cn_distance=None,
                reason="missing_right_backbone_N",
            ),
        )

    cn_distance = _distance(left_c, right_n)

    if cn_distance <= cn_distance_cutoff:
        return (True, None)

    return (
        False,
        BoundaryRecord(
            chain_id=left_residue.chain_id,
            left_residue_label=left_residue.residue_label,
            right_residue_label=right_residue.residue_label,
            left_resseq=left_residue.resseq,
            right_resseq=right_residue.resseq,
            residue_jump=residue_jump,
            cn_distance=cn_distance,
            reason=f"cn_distance_exceeds_cutoff({cn_distance:.3f} > {cn_distance_cutoff:.3f})",
        ),
    )


# ---------------------------------------------------------------------------
# Fragment building
# ---------------------------------------------------------------------------


def split_residues_into_fragments(
    residue_list: list[ResidueRecord],
    cn_distance_cutoff: float,
) -> tuple[list[list[ResidueRecord]], list[BoundaryRecord]]:
    chain_to_residue_list: dict[str, list[ResidueRecord]] = {}
    chain_order: list[str] = []

    for residue in residue_list:
        if residue.chain_id not in chain_to_residue_list:
            chain_to_residue_list[residue.chain_id] = []
            chain_order.append(residue.chain_id)
        chain_to_residue_list[residue.chain_id].append(residue)

    all_fragment_residue_lists: list[list[ResidueRecord]] = []
    boundary_record_list: list[BoundaryRecord] = []

    for chain_id in chain_order:
        residues = chain_to_residue_list[chain_id]
        if not residues:
            continue

        current_fragment: list[ResidueRecord] = [residues[0]]

        for left_residue, right_residue in zip(residues, residues[1:]):
            is_connected, boundary_record = assess_peptide_connection(
                left_residue=left_residue,
                right_residue=right_residue,
                cn_distance_cutoff=cn_distance_cutoff,
            )

            if is_connected:
                current_fragment.append(right_residue)
                continue

            all_fragment_residue_lists.append(current_fragment)
            current_fragment = [right_residue]

            assert boundary_record is not None
            boundary_record_list.append(boundary_record)

        all_fragment_residue_lists.append(current_fragment)

    return all_fragment_residue_lists, boundary_record_list


def write_fragment_pdb(
    fragment_residues: list[ResidueRecord], output_pdb_path: Path
) -> None:
    output_pdb_path.parent.mkdir(parents=True, exist_ok=True)

    with output_pdb_path.open("w", encoding="utf-8") as handle:
        for residue in fragment_residues:
            for line in residue.lines:
                handle.write(line)
        handle.write("TER\n")
        handle.write("END\n")


def _add_public_summary_aliases(summary: dict[str, Any]) -> dict[str, Any]:
    fragment_paths = [
        record["output_pdb_path"] for record in summary.get("fragment_records", [])
    ]

    aliased = dict(summary)
    aliased["success"] = summary.get("status") == "success"
    aliased["split_success"] = summary.get("status") == "success"
    aliased["fragment_output_dir"] = summary.get("output_dir", "")
    aliased["fragment_manifest_path"] = summary.get("manifest_path", "")
    aliased["fragment_log_path"] = summary.get("log_path", "")
    aliased["fragment_paths"] = fragment_paths
    aliased["fragment_pdb_paths"] = fragment_paths
    return aliased


# ---------------------------------------------------------------------------
# Main public function
# ---------------------------------------------------------------------------


def split_protein_for_gap_protonation(
    input_pdb_path: str | Path,
    output_dir: str | Path | None = None,
    pdb_id: str | None = None,
    variant_label: str = DEFAULT_VARIANT_LABEL,
    cn_distance_cutoff: float = DEFAULT_CN_DISTANCE_CUTOFF,
    log_path: str | Path | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """
    Split a protein-only PDB into physically connected fragments for the gaps route.

    Parameters
    ----------
    input_pdb_path
        Protein-only input PDB.
    output_dir
        Directory where fragment PDB files and manifest JSON will be written.
        If omitted, defaults to:
            components/fragments/<variant_label>/
    pdb_id
        Optional explicit PDB ID. If omitted, inferred from input filename.
    variant_label
        Usually "gaps".
    cn_distance_cutoff
        Maximum C(left)-N(right) distance still considered peptide-connected.
    log_path
        Optional explicit log path.
    verbose
        If True, print compact terminal output.

    Returns
    -------
    dict[str, Any]
        Structured split summary.
    """
    input_pdb_path = Path(input_pdb_path).resolve()

    if not input_pdb_path.exists():
        raise FileNotFoundError(f"Input PDB does not exist: {input_pdb_path}")

    resolved_pdb_id = (
        pdb_id.upper() if pdb_id else _infer_pdb_id_from_input_path(input_pdb_path)
    )
    resolved_output_dir = (
        Path(output_dir).resolve()
        if output_dir is not None
        else _default_output_dir(input_pdb_path, variant_label).resolve()
    )
    resolved_log_path = (
        Path(log_path).resolve()
        if log_path is not None
        else _default_log_path(input_pdb_path, resolved_pdb_id, variant_label).resolve()
    )

    if verbose:
        _screen_header(
            "fragment_split · START",
            [
                f"Input PDB    : {input_pdb_path}",
                f"Output dir   : {resolved_output_dir}",
                f"Log file     : {resolved_log_path}",
                f"PDB ID       : {resolved_pdb_id}",
                f"Variant      : {variant_label}",
                f"C-N cutoff   : {cn_distance_cutoff}",
            ],
        )

    _append_log(
        resolved_log_path,
        "split_protein_for_gap_protonation:start",
        [
            f"input_pdb_path      : {input_pdb_path}",
            f"output_dir          : {resolved_output_dir}",
            f"pdb_id              : {resolved_pdb_id}",
            f"variant_label       : {variant_label}",
            f"cn_distance_cutoff  : {cn_distance_cutoff}",
        ],
    )

    residue_list = read_residue_records_from_protein_pdb(input_pdb_path)

    if not residue_list:
        raise ValueError(f"No ATOM residues found in input PDB: {input_pdb_path}")

    fragment_residue_lists, boundary_record_list = split_residues_into_fragments(
        residue_list=residue_list,
        cn_distance_cutoff=cn_distance_cutoff,
    )

    resolved_output_dir.mkdir(parents=True, exist_ok=True)

    fragment_record_list: list[FragmentRecord] = []
    chain_fragment_counter: dict[str, int] = {}
    fragment_index_global = 0

    for fragment_residues in fragment_residue_lists:
        fragment_index_global += 1

        chain_id = fragment_residues[0].chain_id
        chain_fragment_counter[chain_id] = chain_fragment_counter.get(chain_id, 0) + 1
        fragment_index_in_chain = chain_fragment_counter[chain_id]

        chain_label_for_filename = chain_id if chain_id else "blank"

        output_pdb_path = (
            resolved_output_dir
            / f"{resolved_pdb_id}_chain_{chain_label_for_filename}_fragment_{fragment_index_in_chain:02d}.pdb"
        )

        write_fragment_pdb(
            fragment_residues=fragment_residues,
            output_pdb_path=output_pdb_path,
        )

        fragment_record = FragmentRecord(
            chain_id=chain_id,
            fragment_index_global=fragment_index_global,
            fragment_index_in_chain=fragment_index_in_chain,
            residues=fragment_residues,
            output_pdb_path=output_pdb_path,
        )
        fragment_record_list.append(fragment_record)

        if verbose:
            _screen_sub(
                "fragment "
                f"{fragment_index_global:02d} | "
                f"chain={fragment_record.chain_label} | "
                f"{fragment_record.start_resname} {fragment_record.start_resseq}"
                f" -> {fragment_record.end_resname} {fragment_record.end_resseq} | "
                f"residues={fragment_record.residue_count} | "
                f"atoms={fragment_record.atom_count}"
            )

    manifest_path = resolved_output_dir / "fragment_manifest.json"

    summary: dict[str, Any] = {
        "status": "success",
        "module": MODULE_NAME,
        "pdb_id": resolved_pdb_id,
        "variant_label": variant_label,
        "input_pdb_path": str(input_pdb_path),
        "output_dir": str(resolved_output_dir),
        "manifest_path": str(manifest_path),
        "log_path": str(resolved_log_path),
        "cn_distance_cutoff": cn_distance_cutoff,
        "n_input_residues": len(residue_list),
        "n_fragments": len(fragment_record_list),
        "n_boundaries": len(boundary_record_list),
        "fragment_records": [
            {
                "fragment_index_global": fragment.fragment_index_global,
                "fragment_index_in_chain": fragment.fragment_index_in_chain,
                "chain_id": fragment.chain_id,
                "chain_label": fragment.chain_label,
                "start_resseq": fragment.start_resseq,
                "end_resseq": fragment.end_resseq,
                "start_resname": fragment.start_resname,
                "end_resname": fragment.end_resname,
                "residue_count": fragment.residue_count,
                "atom_count": fragment.atom_count,
                "output_pdb_path": str(fragment.output_pdb_path),
            }
            for fragment in fragment_record_list
        ],
        "boundary_records": [
            {
                "chain_id": boundary.chain_id,
                "left_residue_label": boundary.left_residue_label,
                "right_residue_label": boundary.right_residue_label,
                "left_resseq": boundary.left_resseq,
                "right_resseq": boundary.right_resseq,
                "residue_jump": boundary.residue_jump,
                "cn_distance": boundary.cn_distance,
                "reason": boundary.reason,
            }
            for boundary in boundary_record_list
        ],
    }

    with manifest_path.open("w", encoding="utf-8") as handle:
        json.dump(summary, handle, indent=2)

    _append_log(
        resolved_log_path,
        "split_protein_for_gap_protonation:summary",
        [
            f"n_input_residues : {summary['n_input_residues']}",
            f"n_fragments      : {summary['n_fragments']}",
            f"n_boundaries     : {summary['n_boundaries']}",
            f"manifest_path    : {manifest_path}",
        ],
    )

    for boundary in boundary_record_list:
        _append_log(
            resolved_log_path,
            "split_protein_for_gap_protonation:boundary",
            [
                f"chain_id           : {boundary.chain_id if boundary.chain_id else '(blank)'}",
                f"left_residue       : {boundary.left_residue_label}",
                f"right_residue      : {boundary.right_residue_label}",
                f"left_resseq        : {boundary.left_resseq}",
                f"right_resseq       : {boundary.right_resseq}",
                f"residue_jump       : {boundary.residue_jump}",
                f"cn_distance        : {boundary.cn_distance}",
                f"reason             : {boundary.reason}",
            ],
        )

    if verbose:
        _screen_header(
            "fragment_split · END",
            [
                "Status        : success",
                f"Fragments     : {summary['n_fragments']}",
                f"Boundaries    : {summary['n_boundaries']}",
                f"Manifest JSON : {manifest_path}",
            ],
        )

    return _add_public_summary_aliases(summary)


def split_gap_variant_structure(
    pdb_id: str,
    protein_dir: str | Path,
    variant_label: str = DEFAULT_VARIANT_LABEL,
    input_pdb_path: str | Path | None = None,
    output_dir: str | Path | None = None,
    cn_distance_cutoff: float = DEFAULT_CN_DISTANCE_CUTOFF,
    log_path: str | Path | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """
    Backward-/pipeline-compatible wrapper expected by fruton.py.

    Expected call style
    -------------------
    Typical FRUTON usage:

        split_gap_variant_structure(
            pdb_id="1UOU",
            protein_dir=Path("data/proteins/1UOU"),
            variant_label="gaps",
            input_pdb_path=Path(".../components/1UOU_protein.pdb"),
        )

    If input_pdb_path is omitted, the default is:
        <protein_dir>/components/<PDB_ID>_protein.pdb
    """
    protein_dir = Path(protein_dir).resolve()
    resolved_pdb_id = pdb_id.strip().upper()

    resolved_input_pdb_path = (
        Path(input_pdb_path).resolve()
        if input_pdb_path is not None
        else _default_input_pdb_path_for_variant(
            protein_dir=protein_dir,
            pdb_id=resolved_pdb_id,
            variant_label=variant_label,
        ).resolve()
    )

    resolved_output_dir = (
        Path(output_dir).resolve()
        if output_dir is not None
        else protein_dir / "components" / "fragments" / variant_label
    )

    resolved_log_path = (
        Path(log_path).resolve()
        if log_path is not None
        else protein_dir.parent
        / "logs"
        / resolved_pdb_id
        / f"{MODULE_NAME}.{variant_label}.log"
    )

    summary = split_protein_for_gap_protonation(
        input_pdb_path=resolved_input_pdb_path,
        output_dir=resolved_output_dir,
        pdb_id=resolved_pdb_id,
        variant_label=variant_label,
        cn_distance_cutoff=cn_distance_cutoff,
        log_path=resolved_log_path,
        verbose=verbose,
    )

    summary["protein_dir"] = str(protein_dir)
    summary["input_path"] = summary["input_pdb_path"]
    summary["fragment_directory"] = summary["output_dir"]

    return summary


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------


def build_argparser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Split a protein-only PDB into connected fragments for the gaps route."
    )
    parser.add_argument("input_pdb", type=Path, help="Protein-only input PDB path.")
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=None,
        help="Optional output directory for fragment PDB files.",
    )
    parser.add_argument(
        "--pdb-id",
        type=str,
        default=None,
        help="Optional explicit PDB ID.",
    )
    parser.add_argument(
        "--variant-label",
        type=str,
        default=DEFAULT_VARIANT_LABEL,
        help=f"Variant label. Default: {DEFAULT_VARIANT_LABEL}",
    )
    parser.add_argument(
        "--cn-distance-cutoff",
        type=float,
        default=DEFAULT_CN_DISTANCE_CUTOFF,
        help=f"Backbone C-N cutoff in Å. Default: {DEFAULT_CN_DISTANCE_CUTOFF}",
    )
    parser.add_argument(
        "--log-path",
        type=Path,
        default=None,
        help="Optional explicit log path.",
    )
    parser.add_argument(
        "--quiet",
        action="store_true",
        help="Suppress terminal output.",
    )
    return parser


def main() -> int:
    parser = build_argparser()
    args = parser.parse_args()

    split_protein_for_gap_protonation(
        input_pdb_path=args.input_pdb,
        output_dir=args.output_dir,
        pdb_id=args.pdb_id,
        variant_label=args.variant_label,
        cn_distance_cutoff=args.cn_distance_cutoff,
        log_path=args.log_path,
        verbose=not args.quiet,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
