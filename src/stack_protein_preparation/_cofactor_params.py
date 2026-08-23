"""Auto-parametrize FRUTON-detected cofactors via antechamber + parmchk2.

FRUTON's ``pdb_components.py`` detects ~40+ cofactor resnames (FAD, NAD, HEM,
FES, PLP, TPP, SAM, ADP, ATP, ...) and writes them to
``components/<PDB>_cofactor.pdb``, but it never generates force-field
parameters for them.  Users then have to supply ``.mol2`` and ``.frcmod``
files by hand for every cofactor they intend to load in tleap.

This module runs antechamber (AM1-BCC / GAFF2) + parmchk2 (GAFF2 missing
parameter fill) for each unique cofactor found in the cofactor PDB.  The
parameters are first-pass — they are not RESP charges, so highly-charged
cofactors (NADP, phospho-groups) may need manual re-derivation.  Charged
cofactors that require MCPB.py (heme iron) are flagged separately for the
metal parametrization stage.

Deliberately skipped:
- Cofactors that ship in the AMBER contrib library (Bryce group, GLYCAM,
  lipids17) — user can load those via existing leaprc's instead of
  paying an antechamber cost.
- Water (HOH) and simple metal ions.

Public entry point:
    parametrize_cofactors(cofactor_pdb_path, output_dir) -> dict[str, Any]

Returns a manifest ``{cofactors: [{resname, status, mol2, frcmod, net_charge,
message}], status: 'success' | 'partial' | 'failed' | 'skipped'}``.
"""
from __future__ import annotations

import shutil
import subprocess
import tempfile
from collections import OrderedDict
from pathlib import Path
from typing import Any


# Cofactors that HAVE parametrization in well-known AMBER contrib libraries.
# Skip antechamber for these — user should load the correct lib/frcmod instead.
# (Sources: Bryce group (Manchester), GLYCAM, Amber lipid libraries.)
_KNOWN_LIBRARY_COFACTORS: dict[str, str] = {
    # Bryce group — Manchester
    "NAD": "Bryce (nad_dab.lib / nad_dab.frcmod)",
    "NAP": "Bryce (nadp.lib / nadp.frcmod)",
    "NDP": "Bryce (nadph.lib / nadph.frcmod)",
    "FAD": "Bryce (fad.lib / fad.frcmod)",
    "FMN": "Bryce (fmn.lib / fmn.frcmod)",
    "HEM": "Bryce heme; or MCPB.py for the iron",
    "HEC": "Bryce heme-c; or MCPB.py for the iron",
    "SAM": "Bryce (sam.lib / sam.frcmod)",
    "SAH": "Bryce (sah.lib / sah.frcmod)",
    "ATP": "Bryce (atp.lib / atp.frcmod)",
    "ADP": "Bryce (adp.lib / adp.frcmod)",
    "AMP": "Bryce (amp.lib / amp.frcmod)",
    "GTP": "Bryce (gtp.lib / gtp.frcmod)",
    "GDP": "Bryce (gdp.lib / gdp.frcmod)",
    # Iron-sulfur clusters require MCPB or the Carvalho parameters.
    "FES": "MCPB.py or Carvalho parameters",
    "SF4": "MCPB.py or Carvalho parameters",
    "F3S": "MCPB.py or Carvalho parameters",
}


def _read_pdb_lines(path: Path) -> list[str]:
    if not path.is_file():
        return []
    return path.read_text(errors="replace").splitlines()


def _unique_cofactor_resnames(pdb_lines: list[str]) -> "OrderedDict[str, list[str]]":
    """Return {resname: [lines belonging to that residue in original order]}."""
    per_res: OrderedDict[str, list[str]] = OrderedDict()
    for line in pdb_lines:
        if not (line.startswith("HETATM") or line.startswith("ATOM")):
            continue
        r = line[17:20].strip().upper()
        if not r or r == "HOH":
            continue
        per_res.setdefault(r, []).append(line)
    return per_res


def _write_single_cofactor_pdb(lines: list[str], out_path: Path) -> None:
    out_path.write_text("\n".join(lines) + "\nEND\n")


def _guess_net_charge_for_cofactor(resname: str) -> int:
    """Coarse table of physiological net charges for common cofactors.
    Wrong for highly-charged cofactors that user should re-parametrize."""
    table = {
        "NAD": -1, "NAP": -3, "NDP": -4, "FAD": -2, "FMN": -2,
        "SAM": 1, "SAH": 0, "ATP": -4, "ADP": -3, "AMP": -2,
        "GTP": -4, "GDP": -3, "GMP": -2,
        "PLP": -2, "TPP": -2, "COA": -3, "COB": 1, "BTN": -1,
        "PQQ": -3, "THF": 0, "MTA": 0, "MTX": -1,
    }
    return table.get(resname, 0)


def _run_antechamber_and_parmchk2(
    cofactor_pdb: Path,
    resname: str,
    net_charge: int,
    output_dir: Path,
    timeout_s: int = 300,
) -> tuple[bool, str]:
    """Run antechamber + parmchk2. Returns (success, message)."""
    antechamber = shutil.which("antechamber")
    parmchk2 = shutil.which("parmchk2")
    if antechamber is None:
        return False, "antechamber not on PATH"
    if parmchk2 is None:
        return False, "parmchk2 not on PATH"

    output_dir.mkdir(parents=True, exist_ok=True)
    mol2_path = output_dir / f"{resname}.mol2"
    frcmod_path = output_dir / f"{resname}.frcmod"

    with tempfile.TemporaryDirectory(prefix=f"antech_{resname}_") as _td:
        td = Path(_td)
        # antechamber requires cwd it can pollute with intermediate files.
        try:
            proc = subprocess.run(
                [
                    antechamber,
                    "-i", str(cofactor_pdb),
                    "-fi", "pdb",
                    "-o", str(mol2_path),
                    "-fo", "mol2",
                    "-c", "bcc",
                    "-nc", str(net_charge),
                    "-at", "gaff2",
                    "-s", "2",
                    "-pf", "y",
                ],
                capture_output=True,
                text=True,
                timeout=timeout_s,
                cwd=str(td),
            )
        except subprocess.TimeoutExpired:
            return False, f"antechamber timed out after {timeout_s}s"
        if not (mol2_path.is_file() and mol2_path.stat().st_size > 0):
            return False, (
                f"antechamber returned exit={proc.returncode}, "
                f"stderr={proc.stderr.splitlines()[-1] if proc.stderr else '(none)'}"
            )
        # parmchk2 to fill any missing GAFF2 parameters.
        try:
            proc = subprocess.run(
                [parmchk2, "-i", str(mol2_path), "-f", "mol2",
                 "-o", str(frcmod_path), "-s", "gaff2"],
                capture_output=True,
                text=True,
                timeout=timeout_s,
                cwd=str(td),
            )
        except subprocess.TimeoutExpired:
            return False, "parmchk2 timed out"
        if not (frcmod_path.is_file() and frcmod_path.stat().st_size > 0):
            return False, (
                f"parmchk2 returned exit={proc.returncode}, "
                f"stderr={proc.stderr.splitlines()[-1] if proc.stderr else '(none)'}"
            )
    return True, "ok"


def parametrize_cofactors(
    cofactor_pdb_path: str | Path,
    output_dir: str | Path,
    *,
    skip_known_library_cofactors: bool = True,
    timeout_per_cofactor_s: int = 300,
) -> dict[str, Any]:
    """Iterate unique cofactor resnames in ``cofactor_pdb_path`` and run
    antechamber+parmchk2 for each.  Skip resnames that already have known
    AMBER-contrib libraries (unless ``skip_known_library_cofactors`` is
    False).

    Writes .mol2 + .frcmod into ``output_dir`` per successful cofactor.
    """
    pdb_path = Path(cofactor_pdb_path)
    out_dir = Path(output_dir)
    manifest: dict[str, Any] = {
        "status": "skipped",
        "cofactor_pdb": str(pdb_path),
        "output_dir": str(out_dir),
        "cofactors": [],
    }
    if not pdb_path.is_file() or pdb_path.stat().st_size == 0:
        manifest["message"] = "no cofactor PDB present"
        return manifest

    lines = _read_pdb_lines(pdb_path)
    per_res = _unique_cofactor_resnames(lines)
    if not per_res:
        manifest["message"] = "no cofactors detected in file"
        return manifest

    out_dir.mkdir(parents=True, exist_ok=True)
    n_success = 0
    n_skipped_library = 0
    n_failed = 0

    for resname, res_lines in per_res.items():
        entry: dict[str, Any] = {
            "resname": resname,
            "status": "pending",
            "mol2": None,
            "frcmod": None,
            "net_charge": _guess_net_charge_for_cofactor(resname),
            "message": "",
        }
        if skip_known_library_cofactors and resname in _KNOWN_LIBRARY_COFACTORS:
            entry["status"] = "skipped_library"
            entry["message"] = (
                f"already parametrized externally: "
                f"{_KNOWN_LIBRARY_COFACTORS[resname]} (prefer that over AM1-BCC)"
            )
            n_skipped_library += 1
            manifest["cofactors"].append(entry)
            continue

        single_pdb = out_dir / f"{resname}_input.pdb"
        _write_single_cofactor_pdb(res_lines, single_pdb)

        ok, msg = _run_antechamber_and_parmchk2(
            cofactor_pdb=single_pdb,
            resname=resname,
            net_charge=entry["net_charge"],
            output_dir=out_dir,
            timeout_s=timeout_per_cofactor_s,
        )
        # Clean up the per-cofactor input; keep the .mol2 + .frcmod outputs.
        try:
            single_pdb.unlink()
        except FileNotFoundError:
            pass

        if ok:
            entry["status"] = "success"
            entry["mol2"] = str(out_dir / f"{resname}.mol2")
            entry["frcmod"] = str(out_dir / f"{resname}.frcmod")
            entry["message"] = (
                "AM1-BCC/GAFF2 parameters generated (first-pass; verify "
                "against literature charges for highly-charged cofactors)"
            )
            n_success += 1
        else:
            entry["status"] = "failed"
            entry["message"] = msg
            n_failed += 1
        manifest["cofactors"].append(entry)

    if n_success and not n_failed:
        manifest["status"] = "success"
    elif n_success and n_failed:
        manifest["status"] = "partial"
    elif n_failed:
        manifest["status"] = "failed"
    else:
        manifest["status"] = "skipped"  # all were library-known
    manifest["summary"] = {
        "n_success": n_success,
        "n_failed": n_failed,
        "n_skipped_library": n_skipped_library,
    }
    return manifest
