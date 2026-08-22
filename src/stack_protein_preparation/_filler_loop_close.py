"""MODELLER-based loop closure for AF-spliced gap regions.

After ``splice_af_gaps_into_crystal`` inserts AF-modelled residues into the
crystal template, two failure modes remain:

* Mid-loop peptide bond breaks (when AF loop length in cartesian space does
  not match the crystal gap distance) — the splice anchors both ends
  cleanly but a discontinuity appears somewhere inside the loop.
* Steric clashes between AF-modelled side chains and the crystal
  environment (or between AF residues themselves).

An unrestrained MD/min would drag the crystal along.  A restrained OpenMM
min (see ``_filler_junction_relax``) works for isolated chains but has
robustness issues on multi-domain crystals with undeclared chain breaks.

MODELLER's own optimiser handles this cleanly: load the PDB, build
``stereo`` restraints for the whole model (bonds, angles, planarity,
Ramachandran), define a Selection covering ``gap ± flank`` residues, and
run ConjugateGradients on just that Selection.  Everything outside is
held rigid.  MODELLER also refines through 1-3 and 1-4 interactions with
its own soft-sphere potential so we don't need explicit clash handling.

HETATMs (ligands, cofactors, waters) are preserved by copying them back
from the input after MODELLER writes -- MODELLER's atom library only
knows standard amino acids.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable


DEFAULT_FLANK_RESIDUES = 2
DEFAULT_MAX_ITERATIONS = 500


@dataclass
class LoopCloseResult:
    output_pdb_path: Path
    ran: bool
    fallback_reason: str | None = None
    n_free_residues: int = 0
    n_iterations_run: int | None = None
    final_energy: float | None = None

    def to_dict(self) -> dict:
        return {
            "ran": self.ran,
            "fallback_reason": self.fallback_reason,
            "n_free_residues": self.n_free_residues,
            "n_iterations_run": self.n_iterations_run,
            "final_energy": self.final_energy,
        }


def _preserve_hetatms(input_pdb: Path, output_pdb: Path) -> None:
    """Append HETATM/CONECT records from ``input_pdb`` into ``output_pdb``
    before its END record.  MODELLER drops these because top_heav.lib only
    knows standard amino acids.
    """
    het_lines: list[str] = []
    for line in input_pdb.read_text(encoding="utf-8", errors="replace").splitlines():
        if line.startswith("HETATM") or line.startswith("CONECT"):
            het_lines.append(line)
    if not het_lines:
        return
    existing = output_pdb.read_text(encoding="utf-8")
    if "END\n" in existing:
        existing = existing.replace(
            "END\n", "\n".join(het_lines) + "\nEND\n", 1
        )
    else:
        existing = existing.rstrip() + "\n" + "\n".join(het_lines) + "\nEND\n"
    output_pdb.write_text(existing, encoding="utf-8")


def close_loops_via_modeller(
    input_pdb_path: str | Path,
    output_pdb_path: str | Path,
    gap_ranges: Iterable[tuple[int, int]],
    flank_residues: int = DEFAULT_FLANK_RESIDUES,
    max_iterations: int = DEFAULT_MAX_ITERATIONS,
) -> LoopCloseResult:
    """Refine only ``gap ± flank`` residues via MODELLER ConjugateGradients.

    ``gap_ranges`` is an iterable of ``(first_resnum, last_resnum)`` pairs
    identifying the AF-inserted residue ranges.  Flank residues on each
    side are also freed so peptide bonds at the crystal/AF junction can
    absorb strain.  Every atom outside the free set is held rigid.

    Failure modes fall back to identity (copy input to output) so the
    downstream pipeline always sees the expected artifact.
    """
    input_path = Path(input_pdb_path)
    output_path = Path(output_pdb_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    try:
        from modeller import Environ, Model, Selection
        from modeller.optimizers import ConjugateGradients
    except Exception as exc:  # noqa: BLE001
        return LoopCloseResult(
            output_pdb_path=input_path,
            ran=False,
            fallback_reason=f"modeller unavailable: {exc!r}",
        )

    env = Environ()
    env.io.hetatm = True  # keep HETATM references in topology (positional only)
    env.libs.topology.read(file="$(LIB)/top_heav.lib")
    env.libs.parameters.read(file="$(LIB)/par.lib")
    env.edat.dynamic_sphere = True
    env.edat.dynamic_lennard = True

    try:
        mdl = Model(env, file=str(input_path))
        mdl.generate_topology(sequence=None)
        mdl.transfer_xyz(mdl, cluster_cut=-1.0)
    except Exception as exc:  # noqa: BLE001
        # Some MODELLER versions require a Sequence argument for
        # generate_topology; fall back to the simpler complete_pdb path.
        try:
            from modeller.scripts import complete_pdb as _complete_pdb
            mdl = _complete_pdb(env, str(input_path))
        except Exception as exc2:  # noqa: BLE001
            output_path.write_bytes(input_path.read_bytes())
            return LoopCloseResult(
                output_pdb_path=output_path,
                ran=False,
                fallback_reason=f"MODELLER Model+topology load failed: {exc!r} / {exc2!r}",
            )

    try:
        mdl.restraints.make(mdl, restraint_type="stereo", spline_on_site=False)
    except Exception as exc:  # noqa: BLE001
        output_path.write_bytes(input_path.read_bytes())
        return LoopCloseResult(
            output_pdb_path=output_path,
            ran=False,
            fallback_reason=f"MODELLER restraints.make failed: {exc!r}",
        )

    # Build free-residue selection: gap ± flank on every chain.
    free_sel = Selection()
    n_added = 0
    for chain in mdl.chains:
        chain_name = chain.name
        for lo, hi in gap_ranges:
            free_lo = lo - flank_residues
            free_hi = hi + flank_residues
            try:
                res_range = mdl.residue_range(
                    f"{free_lo}:{chain_name}", f"{free_hi}:{chain_name}"
                )
                free_sel.add(res_range)
                n_added += 1
            except (KeyError, IndexError, Exception):  # noqa: BLE001
                # residue range not present in this chain -- skip
                continue

    if n_added == 0:
        output_path.write_bytes(input_path.read_bytes())
        return LoopCloseResult(
            output_pdb_path=output_path,
            ran=False,
            fallback_reason="no gap residues resolved in any chain",
        )

    n_free_residues = 0
    try:
        n_free_residues = len(list(free_sel.only_atom_types("CA")))
    except Exception:  # noqa: BLE001
        n_free_residues = 0

    try:
        # MODELLER's ConjugateGradients.optimize(atmsel, ...) moves only the
        # atoms in atmsel; all other atoms are treated as fixed for the
        # duration of the optimisation.  No explicit fix-all step needed.
        cg = ConjugateGradients(output="NO_REPORT")
        cg.optimize(free_sel, max_iterations=max_iterations, actions=None)
    except Exception as exc:  # noqa: BLE001
        output_path.write_bytes(input_path.read_bytes())
        return LoopCloseResult(
            output_pdb_path=output_path,
            ran=False,
            fallback_reason=f"MODELLER ConjugateGradients failed: {exc!r}",
            n_free_residues=n_free_residues,
        )

    final_energy = None
    try:
        final_energy = float(mdl.energy()[0])
    except Exception:  # noqa: BLE001
        pass

    try:
        mdl.write(file=str(output_path))
    except Exception as exc:  # noqa: BLE001
        output_path.write_bytes(input_path.read_bytes())
        return LoopCloseResult(
            output_pdb_path=output_path,
            ran=False,
            fallback_reason=f"MODELLER write failed: {exc!r}",
            n_free_residues=n_free_residues,
        )

    # Divergence check: local minimisation on badly-fit loops (mid-loop
    # peptide bond breaks > 10 A across a physical gap) can push atoms to
    # infinity.  PDB writer emits '****' when coords exceed the fixed-width
    # field.  If we see any, the optimisation blew up -- revert to input.
    _diverged = "****" in output_path.read_text(encoding="utf-8", errors="replace")
    if _diverged:
        output_path.write_bytes(input_path.read_bytes())
        return LoopCloseResult(
            output_pdb_path=output_path,
            ran=False,
            fallback_reason=(
                "MODELLER optimisation diverged (coords > 9999 A) -- reverted to input. "
                "Cause: gap has broken peptide bonds > 10 A that local minimisation "
                "cannot bridge; needs LoopModel or KIC."
            ),
            n_free_residues=n_free_residues,
        )

    _preserve_hetatms(input_path, output_path)

    return LoopCloseResult(
        output_pdb_path=output_path,
        ran=True,
        fallback_reason=None,
        n_free_residues=n_free_residues,
        n_iterations_run=max_iterations,
        final_energy=final_energy,
    )
