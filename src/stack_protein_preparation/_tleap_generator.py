"""Assemble a complete ``tleap.in`` script from FRUTON's prep outputs.

This module is the "glue" that consumes:
  * a final protonated protein PDB (``<PDB>_proteinH.pdb`` from the
    protonation stage);
  * a list of detected ion resnames (from ``pdb_components``) that get
    routed to the appropriate ``leaprc.water.tip3p_ion12-6-4*`` via
    :mod:`stack_protein_preparation._ion_params`;
  * an optional cofactor parameter directory (``<resname>.mol2`` +
    ``<resname>.frcmod`` files produced by
    :mod:`stack_protein_preparation._cofactor_params`);
  * an optional non-standard residue directory (same schema, produced
    by :mod:`stack_protein_preparation._nonstd_residue_params_core`);

...and writes a ``tleap.in`` that a user can run directly on a SLURM
node with AmberTools loaded to produce a ready-to-simulate topology
(``system.prmtop`` + ``system.rst7``).

Design decisions:

* **Emitted script is self-contained** -- one file the user runs
  verbatim, no post-hoc editing.  All ``source`` and ``loadamberparams``
  paths are absolute so tleap can be launched from any cwd.

* **Never invokes tleap.**  The pipeline stage that RUNS tleap belongs
  in a separate validation module (``COMPLETENESS 5 / task #91``).
  This module is deterministic string generation only -- no subprocess,
  no side effect beyond writing one text file.

* **Fail-open on missing components.**  If the cofactor / nonstandard
  directories don't exist, their sections are omitted with a comment
  explaining why.  User can rerun after supplying missing params.

* **Ion-source dedup** via :func:`_ion_params.tleap_source_lines`, so
  a protein with Zn + Mg + Ca gets one 12-6-4 source line, and adding
  a Gd adds the HFE-tuned lanthanide leaprc as a second line.

* **Solvent / neutralisation** parameters are user-settable but come
  with sensible defaults (TIP3P, 10 A truncated octahedron, Na+/Cl-
  neutralisation).

License-free: pure Python string generation, csv-driven data lookups.
No external tool invocation.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable


DEFAULT_PROTEIN_FF = "leaprc.protein.ff14SB"
DEFAULT_WATER_FF = "leaprc.water.tip3p"
DEFAULT_GAFF = "leaprc.gaff2"
DEFAULT_SOLVATE_MARGIN_ANGSTROM = 10.0
DEFAULT_TOPOLOGY_BASENAME = "system"


@dataclass
class TleapGenerationResult:
    tleap_in_path: Path
    lines: list[str]
    n_cofactors_loaded: int = 0
    n_nonstandard_loaded: int = 0
    n_ion_leaprc_lines: int = 0
    warnings: list[str] = field(default_factory=list)

    def to_dict(self) -> dict:
        return {
            "tleap_in_path": str(self.tleap_in_path),
            "n_cofactors_loaded": self.n_cofactors_loaded,
            "n_nonstandard_loaded": self.n_nonstandard_loaded,
            "n_ion_leaprc_lines": self.n_ion_leaprc_lines,
            "warnings": list(self.warnings),
        }


def _find_param_pairs(param_dir: Path | None) -> list[tuple[str, Path, Path]]:
    """Return ``[(resname, mol2, frcmod), ...]`` for every complete pair
    in ``param_dir``.  Silently skip orphan files (mol2 without frcmod
    or vice versa)."""
    if param_dir is None or not param_dir.is_dir():
        return []
    out: list[tuple[str, Path, Path]] = []
    for mol2 in sorted(param_dir.glob("*.mol2")):
        resname = mol2.stem.upper()
        frcmod = param_dir / f"{mol2.stem}.frcmod"
        if frcmod.is_file():
            out.append((resname, mol2.resolve(), frcmod.resolve()))
    return out


def write_tleap_script(
    final_model_pdb: str | Path,
    output_dir: str | Path,
    *,
    ion_resnames: Iterable[str] = (),
    cofactor_param_dir: str | Path | None = None,
    nonstandard_param_dir: str | Path | None = None,
    protein_ff: str = DEFAULT_PROTEIN_FF,
    water_ff: str = DEFAULT_WATER_FF,
    solvate_margin_angstrom: float = DEFAULT_SOLVATE_MARGIN_ANGSTROM,
    topology_basename: str = DEFAULT_TOPOLOGY_BASENAME,
    neutralise: bool = True,
) -> TleapGenerationResult:
    """Emit ``<output_dir>/tleap.in`` covering the full setup and return
    the generated line list for logging.  See module docstring for
    the semantics of each argument.
    """
    from stack_protein_preparation._ion_params import (
        tleap_source_lines as _ion_tleap_lines,
    )

    model_path = Path(final_model_pdb).resolve()
    out_dir = Path(output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    tleap_in = out_dir / "tleap.in"

    warnings: list[str] = []
    if not model_path.is_file():
        warnings.append(
            f"final_model_pdb {model_path!r} does not exist -- "
            f"tleap will fail at loadpdb"
        )

    cofactor_pairs = _find_param_pairs(
        Path(cofactor_param_dir) if cofactor_param_dir else None
    )
    nonstd_pairs = _find_param_pairs(
        Path(nonstandard_param_dir) if nonstandard_param_dir else None
    )
    ion_lines = _ion_tleap_lines(list(ion_resnames))

    L: list[str] = []
    L.append("# =========================================================")
    L.append(f"# tleap.in -- auto-generated by FRUTON _tleap_generator")
    L.append(f"# model:            {model_path}")
    L.append(f"# protein_ff:       {protein_ff}")
    L.append(f"# water_model:      {water_ff}")
    L.append(f"# solvate margin:   {solvate_margin_angstrom} A")
    L.append(f"# neutralise:       {neutralise}")
    L.append(f"# ion resnames:     {sorted(set(r.upper() for r in ion_resnames)) or '(none)'}")
    L.append(f"# cofactor params:  {len(cofactor_pairs)} pair(s)")
    L.append(f"# nonstd params:    {len(nonstd_pairs)} pair(s)")
    L.append("# =========================================================")
    L.append("")

    # Force fields
    L.append(f"source {protein_ff}")
    L.append(f"source {water_ff}")
    if cofactor_pairs or nonstd_pairs:
        L.append(f"source {DEFAULT_GAFF}")
    for ion_line in ion_lines:
        L.append(ion_line)
    L.append("")

    # Cofactor params
    if cofactor_pairs:
        L.append("# --- cofactor parameters (antechamber + parmchk2) ---")
        for resname, mol2, frcmod in cofactor_pairs:
            L.append(f"{resname} = loadmol2 {mol2}")
            L.append(f"loadamberparams {frcmod}")
        L.append("")
    else:
        L.append("# (no cofactor params loaded)")
        L.append("")

    # Non-standard residue params (phospho, ALY, etc.)
    if nonstd_pairs:
        L.append("# --- non-standard residue parameters ---")
        for resname, mol2, frcmod in nonstd_pairs:
            L.append(f"{resname} = loadmol2 {mol2}")
            L.append(f"loadamberparams {frcmod}")
        L.append("")
    else:
        L.append("# (no non-standard residue params loaded)")
        L.append("")

    # Model
    L.append(f"# --- load prepared model ---")
    L.append(f"mol = loadpdb {model_path}")
    L.append("")

    # Solvation + neutralisation
    water_box_name = "TIP3PBOX" if "tip3p" in water_ff.lower() else "TIP3PBOX"
    L.append("# --- solvate + neutralise ---")
    L.append(f"solvatebox mol {water_box_name} {solvate_margin_angstrom:.1f}")
    if neutralise:
        L.append("addions mol Na+ 0")
        L.append("addions mol Cl- 0")
    L.append("")

    # Sanity check
    L.append("# --- sanity check + save ---")
    L.append("check mol")
    L.append(f"saveamberparm mol {topology_basename}.prmtop {topology_basename}.rst7")
    L.append(f"savepdb mol {topology_basename}.pdb")
    L.append("quit")
    L.append("")

    tleap_in.write_text("\n".join(L), encoding="utf-8")

    return TleapGenerationResult(
        tleap_in_path=tleap_in,
        lines=L,
        n_cofactors_loaded=len(cofactor_pairs),
        n_nonstandard_loaded=len(nonstd_pairs),
        n_ion_leaprc_lines=len(ion_lines),
        warnings=warnings,
    )
