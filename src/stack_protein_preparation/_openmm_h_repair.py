"""OpenMM-driven hydrogen re-placement for pdb2gmx outputs that fail Meeko.

Called as a post-success retry engine from ``protonation.py`` when the
Meeko compatibility gate reports an ``AtomValenceException``. The retry
does not touch heavy atoms — only Hs are stripped and re-added by
``openmm.app.Modeller.addHydrogens`` using the amber14 protein forcefield
at the caller-supplied pH.

Why OpenMM Modeller specifically:
  - No external tool (no pdbfixer, no reduce, no tleap subprocess).
  - Deterministic and residue-aware — respects HIS/GLU/ASP/LYS/CYS ionizable
    states via the pH argument.
  - Uses amber14 heavy-atom templates that overlap AMBER's ff14SB/99SB-ILDN
    conventions, so downstream Meeko/AutoDock preparation stops seeing the
    valence conflicts from GROMACS's slightly different H-naming.
  - Already a hard-dependency of FRUTON (see pixi.toml ``openmm>=8.4.0``).
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class OpenMMHRepairResult:
    """Outcome of a Modeller.addHydrogens retry."""

    ok: bool
    """True when a repaired PDB was written."""

    output_path: Path | None
    """Path of the rewritten PDB (== input to this function on success)."""

    n_hydrogens_added: int
    """Number of H atoms in the repaired structure (heavy atoms preserved)."""

    error_message: str
    """Empty on success; otherwise the exception's string form."""


def repair_hydrogens_openmm(
    pdb_path: str | Path,
    *,
    ph: float = 7.4,
    forcefield_files: tuple[str, ...] = ("amber14-all.xml", "amber14/tip3pfb.xml"),
) -> OpenMMHRepairResult:
    """Strip and re-add hydrogens on ``pdb_path`` using OpenMM Modeller.

    Overwrites ``pdb_path`` in place on success. Never raises — every failure
    mode becomes an ``OpenMMHRepairResult`` with ``ok=False``.

    The ``forcefield_files`` default matches AMBER14 + TIP3P (Fox-Berthiaume);
    override if the caller has a specific forcefield already loaded elsewhere.
    """

    path = Path(pdb_path)
    if not path.is_file() or path.stat().st_size == 0:
        return OpenMMHRepairResult(
            ok=False,
            output_path=None,
            n_hydrogens_added=0,
            error_message=f"pdb file missing or empty: {path}",
        )

    try:
        from openmm.app import ForceField, Modeller, PDBFile
        from openmm.app.element import hydrogen
    except ImportError as exc:
        return OpenMMHRepairResult(
            ok=False,
            output_path=None,
            n_hydrogens_added=0,
            error_message=f"openmm not importable: {exc}",
        )

    try:
        pdb = PDBFile(str(path))
        modeller = Modeller(pdb.topology, pdb.positions)
        # Strip any hydrogens pdb2gmx placed. Modeller.addHydrogens rebuilds
        # them from the FF template so the H-naming matches AMBER exactly,
        # which is what Meeko's Polymer template matcher expects.
        modeller.delete([a for a in modeller.topology.atoms() if a.element == hydrogen])
        forcefield = ForceField(*forcefield_files)
        modeller.addHydrogens(forcefield, pH=ph)
        with path.open("w", encoding="utf-8") as fh:
            PDBFile.writeFile(modeller.topology, modeller.positions, fh, keepIds=True)
        n_h = sum(1 for a in modeller.topology.atoms() if a.element == hydrogen)
        return OpenMMHRepairResult(
            ok=True, output_path=path, n_hydrogens_added=n_h, error_message=""
        )
    except Exception as exc:  # noqa: BLE001 — validation shim, must swallow
        return OpenMMHRepairResult(
            ok=False,
            output_path=None,
            n_hydrogens_added=0,
            error_message=f"{type(exc).__name__}: {exc}",
        )
