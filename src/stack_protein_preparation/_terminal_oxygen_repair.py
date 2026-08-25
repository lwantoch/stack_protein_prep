"""Rebuild C-terminal backbone O + OXT atoms if missing.

BUG CONTEXT
-----------
FRUTON's Stage 10 protonation pipeline (gmx pdb2gmx or OpenMM Modeller
addHydrogens) silently drops the terminal backbone O and OXT atoms on
some C-terminal residues.  Downstream tools that expect standard PDB
terminals (Meeko for docking prep, RDKit-based tools, etc.) then fail
with 'No template matched residue X' errors.

FRUTON's own tleap-based sanity check does not catch this because tleap
auto-adds terminal oxygens during topology building. GROMACS `pdb2gmx`
also auto-adds. Meeko does NOT.

This module reconstructs the missing terminal oxygens geometrically from
the backbone N-CA-C plane using idealized bond geometry:
- C=O bond length 1.24 A
- OXT-C bond length 1.24 A
- O-C-OXT angle 122.4° (standard carboxylate)
- Both oxygens in the N-CA-C plane

Reference: standard peptide C-terminus geometry, Engh & Huber
Acta Cryst A47:392 (1991) DOI 10.1107/S0108767391001071.
"""
from __future__ import annotations

import math
from pathlib import Path
from typing import Iterable


def _dot(a, b):
    return sum(x * y for x, y in zip(a, b))


def _sub(a, b):
    return [ax - bx for ax, bx in zip(a, b)]


def _add(a, b):
    return [ax + bx for ax, bx in zip(a, b)]


def _scale(a, s):
    return [x * s for x in a]


def _cross(a, b):
    return [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]


def _norm(a):
    return math.sqrt(_dot(a, a))


def _unit(a):
    n = _norm(a)
    return [x / n for x in a] if n > 0 else a


def _rotate_around_axis(v, axis, angle_rad):
    """Rodrigues rotation: rotate v around unit axis by angle_rad."""
    c, s = math.cos(angle_rad), math.sin(angle_rad)
    dot = _dot(axis, v)
    cross = _cross(axis, v)
    return [
        v[i] * c + cross[i] * s + axis[i] * dot * (1 - c)
        for i in range(3)
    ]


def reconstruct_terminal_oxygens(
    N: list[float], CA: list[float], C: list[float]
) -> tuple[list[float], list[float]]:
    """Return (O, OXT) coords for a C-terminal residue given backbone N, CA, C.

    Places O and OXT symmetrically in the N-CA-C plane, at 122.4° apart
    from each other (standard carboxylate), each 1.24 A from C, symmetric
    around the extension of the CA-C bond.
    """
    # Normal to N-CA-C plane
    v_ca_n = _sub(N, CA)
    v_ca_c = _sub(C, CA)
    normal = _unit(_cross(v_ca_c, v_ca_n))

    # Direction pointing "away" from CA along CA-C extension
    dir_out = _unit(v_ca_c)

    # O and OXT at ±61.2° from the CA-C extension (122.4° total between them)
    half_angle = math.radians(122.4 / 2)
    dir_o = _rotate_around_axis(dir_out, normal, half_angle)
    dir_oxt = _rotate_around_axis(dir_out, normal, -half_angle)

    bond_len = 1.24
    O = _add(C, _scale(dir_o, bond_len))
    OXT = _add(C, _scale(dir_oxt, bond_len))
    return O, OXT


def _parse_pdb_atom(line: str) -> dict | None:
    if not (line.startswith("ATOM") or line.startswith("HETATM")):
        return None
    try:
        return {
            "record": line[0:6].strip(),
            "serial": int(line[6:11]),
            "atom": line[12:16].strip(),
            "altloc": line[16],
            "resname": line[17:20].strip(),
            "chain": line[21],
            "resnum": int(line[22:26]),
            "icode": line[26],
            "x": float(line[30:38]),
            "y": float(line[38:46]),
            "z": float(line[46:54]),
            "occ": line[54:60],
            "bfac": line[60:66],
            "element": line[76:78].strip() if len(line) >= 78 else "",
            "raw": line,
        }
    except (ValueError, IndexError):
        return None


def _format_atom_line(
    serial: int, atom: str, resname: str, chain: str, resnum: int,
    x: float, y: float, z: float, element: str = ""
) -> str:
    # Atom name field width 4, right-just for 4-char names, left-just for shorter
    if len(atom) < 4:
        atom_field = f" {atom:<3}"
    else:
        atom_field = atom
    return (
        f"ATOM  {serial:>5} {atom_field}{'':1}{resname:>3} {chain}{resnum:>4}{'':1}   "
        f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00          {element:>2}\n"
    )


def repair_pdb_terminal_oxygens(input_pdb: Path, output_pdb: Path) -> dict:
    """Read `input_pdb`, add missing terminal O + OXT to each C-terminal
    protein residue per chain, write `output_pdb`. Return summary dict.
    """
    lines = Path(input_pdb).read_text().splitlines(keepends=True)

    # Pass 1: per-chain, find last ATOM residue (protein, not HOH)
    per_chain: dict[str, dict] = {}
    for i, ln in enumerate(lines):
        a = _parse_pdb_atom(ln)
        if a is None or a["record"] != "ATOM":
            continue
        if a["resname"] in ("HOH", "WAT"):
            continue
        ch = a["chain"]
        key = (ch, a["resnum"])
        rec = per_chain.setdefault(ch, {"last_key": None, "atoms": {}, "lines": []})
        if rec["last_key"] is None or a["resnum"] > rec["last_key"][1]:
            rec["last_key"] = key
            rec["atoms"] = {a["atom"]: a}
            rec["lines"] = [i]
        elif key == rec["last_key"]:
            rec["atoms"][a["atom"]] = a
            rec["lines"].append(i)

    # Pass 2: for each C-terminal, check O + OXT presence
    to_add: list[tuple[int, dict]] = []  # (insert_after_line_idx, atom_dict)
    max_serial = 0
    for ln in lines:
        a = _parse_pdb_atom(ln)
        if a is not None and a["record"] in ("ATOM", "HETATM"):
            max_serial = max(max_serial, a["serial"])
    next_serial = max_serial + 1

    summary = {"chains_processed": 0, "chains_needing_repair": 0, "added": []}
    for ch, rec in per_chain.items():
        summary["chains_processed"] += 1
        atoms = rec["atoms"]
        needs_o = "O" not in atoms
        needs_oxt = "OXT" not in atoms
        if not (needs_o or needs_oxt):
            continue
        # Need N, CA, C
        if not all(k in atoms for k in ("N", "CA", "C")):
            summary["added"].append({
                "chain": ch, "resnum": rec["last_key"][1],
                "skipped_reason": f"missing backbone atoms {sorted(atoms.keys())}"
            })
            continue

        summary["chains_needing_repair"] += 1
        N = [atoms["N"]["x"], atoms["N"]["y"], atoms["N"]["z"]]
        CA = [atoms["CA"]["x"], atoms["CA"]["y"], atoms["CA"]["z"]]
        C = [atoms["C"]["x"], atoms["C"]["y"], atoms["C"]["z"]]
        O_xyz, OXT_xyz = reconstruct_terminal_oxygens(N, CA, C)

        template = atoms["C"]
        insert_after = max(rec["lines"])

        if needs_o:
            to_add.append((insert_after, {
                "serial": next_serial, "atom": "O", "resname": template["resname"],
                "chain": ch, "resnum": template["resnum"],
                "x": O_xyz[0], "y": O_xyz[1], "z": O_xyz[2], "element": "O",
            }))
            next_serial += 1
        if needs_oxt:
            to_add.append((insert_after, {
                "serial": next_serial, "atom": "OXT", "resname": template["resname"],
                "chain": ch, "resnum": template["resnum"],
                "x": OXT_xyz[0], "y": OXT_xyz[1], "z": OXT_xyz[2], "element": "O",
            }))
            next_serial += 1
        summary["added"].append({
            "chain": ch, "resnum": template["resnum"],
            "resname": template["resname"],
            "added_O": needs_o, "added_OXT": needs_oxt,
        })

    # Pass 3: write with new atoms inserted after their target lines
    to_add.sort(key=lambda p: p[0])
    result_lines: list[str] = []
    inserts_by_line: dict[int, list[dict]] = {}
    for insert_after, atom_dict in to_add:
        inserts_by_line.setdefault(insert_after, []).append(atom_dict)

    for i, ln in enumerate(lines):
        result_lines.append(ln)
        if i in inserts_by_line:
            for atom_dict in inserts_by_line[i]:
                result_lines.append(_format_atom_line(
                    atom_dict["serial"], atom_dict["atom"], atom_dict["resname"],
                    atom_dict["chain"], atom_dict["resnum"],
                    atom_dict["x"], atom_dict["y"], atom_dict["z"],
                    atom_dict.get("element", ""),
                ))

    Path(output_pdb).write_text("".join(result_lines))
    return summary


def repair_pdb_terminal_oxygens_in_place(pdb_path: Path) -> dict:
    return repair_pdb_terminal_oxygens(pdb_path, pdb_path)
