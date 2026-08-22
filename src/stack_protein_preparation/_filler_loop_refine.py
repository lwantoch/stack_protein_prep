"""MODELLER LoopModel-based loop refinement for AF-spliced FRUTON output.

Replaces the previous _filler_loop_close.py which used ConjugateGradients
directly and diverged on badly-fit loops.  MODELLER's LoopModel does
proper stochastic sampling with peptide-bond distance restraints built
into the internal coordinate system, so peptide bonds always end near
1.33 A even when the starting coordinates are pathological (broken
splice fits from AF-loop / crystal-gap length mismatches).

The refinement samples multiple loop conformers, scores each by DOPE,
and picks the best.  Trade-off: 3-5 conformers with ``refine.slow``
takes ~30-90 s per gap on modern CPUs, but produces reviewer-grade
loop geometry that no local minimisation can match.

Chirality guard: MODELLER can produce D-amino-acid outliers (~4/1085
observed in a 5M7U test).  We detect and reject conformers whose
Cα signed-volume chirality has more D-outliers than the input.

Interface mirrors ``_filler_junction_relax`` / ``_filler_loop_close``:
returns a dataclass with ``ran`` / ``fallback_reason`` and a numeric
diagnostic block, and always writes ``output_pdb`` so downstream stages
find their expected artefact.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable


DEFAULT_N_CONFORMERS = 3
DEFAULT_REFINE_LEVEL = "slow"  # very_fast | fast | slow | slow_large | very_slow


@dataclass
class LoopRefineResult:
    output_pdb_path: Path
    ran: bool
    fallback_reason: str | None = None
    n_conformers_built: int = 0
    n_conformers_kept: int = 0
    best_dope: float | None = None
    diagnostics: list[str] = field(default_factory=list)

    def to_dict(self) -> dict:
        return {
            "ran": self.ran,
            "fallback_reason": self.fallback_reason,
            "n_conformers_built": self.n_conformers_built,
            "n_conformers_kept": self.n_conformers_kept,
            "best_dope": self.best_dope,
            "diagnostics": self.diagnostics,
        }


def _chirality_d_count(pdb_path: Path) -> int:
    """Count Cα atoms whose signed tetrahedron volume is negative
    (D-amino-acid).  Uses the same convention as _filler_quality_check.
    """
    import math
    from Bio.PDB import PDBParser
    try:
        s = PDBParser(QUIET=True).get_structure("m", str(pdb_path))
    except Exception:
        return -1
    d = 0
    for chain in s[0]:
        for res in chain:
            if res.id[0].strip() or res.resname == "GLY":
                continue
            if not all(a in res for a in ("N", "CA", "C", "CB")):
                continue
            ca = res["CA"].coord
            vn = res["N"].coord - ca
            vc = res["C"].coord - ca
            vcb = res["CB"].coord - ca
            cross = (
                vn[1] * vc[2] - vn[2] * vc[1],
                vn[2] * vc[0] - vn[0] * vc[2],
                vn[0] * vc[1] - vn[1] * vc[0],
            )
            sv = float(cross[0] * vcb[0] + cross[1] * vcb[1] + cross[2] * vcb[2])
            if sv < 0:
                d += 1
    return d


def refine_loops_via_modeller(
    input_pdb_path: str | Path,
    output_pdb_path: str | Path,
    gap_ranges_by_chain: Iterable[tuple[str, int, int]],
    n_conformers: int = DEFAULT_N_CONFORMERS,
    refine_level: str = DEFAULT_REFINE_LEVEL,
    reject_new_chirality_d: bool = True,
) -> LoopRefineResult:
    """Refine loop regions via MODELLER LoopModel.

    ``gap_ranges_by_chain``: iterable of ``(chain_id, first_resnum, last_resnum)``
    tuples naming the AF-inserted regions to refine.

    Multiple conformers are built (``n_conformers``) with the chosen
    ``refine_level`` schedule; the lowest-DOPE model is written to
    ``output_pdb_path``.  If ``reject_new_chirality_d`` is True and every
    conformer introduces additional D-amino-acid outliers vs input,
    fall back to input (safer to ship un-refined than chirality-broken).

    Always writes output_pdb (falls back to input on failure) so downstream
    stages find the expected artefact.
    """
    import os
    import tempfile

    input_path = Path(input_pdb_path)
    output_path = Path(output_pdb_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    gap_list = list(gap_ranges_by_chain)
    if not gap_list:
        output_path.write_bytes(input_path.read_bytes())
        return LoopRefineResult(
            output_pdb_path=output_path,
            ran=False,
            fallback_reason="no gap ranges provided",
        )

    try:
        from modeller import Environ, Selection
        from modeller.automodel import LoopModel, refine as _refine_mod, assess
    except Exception as exc:  # noqa: BLE001
        output_path.write_bytes(input_path.read_bytes())
        return LoopRefineResult(
            output_pdb_path=output_path,
            ran=False,
            fallback_reason=f"modeller unavailable: {exc!r}",
        )

    _refine_map = {
        "very_fast": _refine_mod.very_fast,
        "fast": _refine_mod.fast,
        "slow": _refine_mod.slow,
        "slow_large": _refine_mod.slow_large,
        "very_slow": _refine_mod.very_slow,
    }
    refine_choice = _refine_map.get(refine_level, _refine_mod.slow)

    env = Environ()
    env.io.hetatm = True
    env.libs.topology.read(file="$(LIB)/top_heav.lib")
    env.libs.parameters.read(file="$(LIB)/par.lib")

    # LoopModel writes into cwd -- run in a temp dir to isolate artefacts.
    workdir = Path(tempfile.mkdtemp(prefix="fruton_loop_refine_"))
    orig_cwd = os.getcwd()

    # Copy input to workdir so MODELLER can find it via basename.
    local_input = workdir / input_path.name
    local_input.write_bytes(input_path.read_bytes())

    baseline_d = _chirality_d_count(input_path)

    _gaps_local = list(gap_list)
    _sequence_name = "target_" + input_path.stem.replace(".", "_")

    class _Refiner(LoopModel):
        def select_loop_atoms(self):
            sel = Selection()
            for chain_id, lo, hi in _gaps_local:
                try:
                    sel.add(self.residue_range(f"{lo}:{chain_id}", f"{hi}:{chain_id}"))
                except Exception:  # noqa: BLE001
                    continue
            return sel

    try:
        os.chdir(workdir)
        m = _Refiner(
            env,
            inimodel=local_input.name,
            sequence=_sequence_name,
            loop_assess_methods=(assess.DOPE,),
        )
        m.loop.starting_model = 1
        m.loop.ending_model = max(1, int(n_conformers))
        m.loop.md_level = refine_choice
        m.make()
    except Exception as exc:  # noqa: BLE001
        os.chdir(orig_cwd)
        output_path.write_bytes(input_path.read_bytes())
        return LoopRefineResult(
            output_pdb_path=output_path,
            ran=False,
            fallback_reason=f"MODELLER LoopModel.make failed: {exc!r}",
        )
    finally:
        try:
            os.chdir(orig_cwd)
        except Exception:
            pass

    # Collect conformers with their DOPE scores from m.loop.outputs.
    conformer_scores: list[tuple[float, Path]] = []
    for entry in (getattr(m.loop, "outputs", []) or []):
        try:
            p = workdir / entry["name"]
            score = float(entry.get("DOPE score", 0.0))
        except Exception:
            continue
        if not p.is_file():
            continue
        conformer_scores.append((score, p))

    n_built = len(conformer_scores)
    if n_built == 0:
        output_path.write_bytes(input_path.read_bytes())
        return LoopRefineResult(
            output_pdb_path=output_path,
            ran=False,
            fallback_reason="MODELLER LoopModel produced 0 conformers",
        )

    # Chirality filter: drop conformers with more D-outliers than input.
    kept: list[tuple[float, Path, int]] = []
    diagnostics: list[str] = []
    for score, p in conformer_scores:
        d = _chirality_d_count(p)
        diagnostics.append(f"{p.name}: DOPE={score:.1f}, chirality_D={d}")
        if reject_new_chirality_d and baseline_d >= 0 and d > baseline_d:
            continue
        kept.append((score, p, d))

    if not kept:
        output_path.write_bytes(input_path.read_bytes())
        return LoopRefineResult(
            output_pdb_path=output_path,
            ran=False,
            fallback_reason=(
                f"all {n_built} conformers introduced chirality D-outliers "
                f"(baseline {baseline_d}); reverted to input"
            ),
            n_conformers_built=n_built,
            n_conformers_kept=0,
            diagnostics=diagnostics,
        )

    # Pick lowest DOPE among survivors.
    kept.sort(key=lambda x: x[0])
    best_score, best_path, best_d = kept[0]
    output_path.write_bytes(best_path.read_bytes())

    return LoopRefineResult(
        output_pdb_path=output_path,
        ran=True,
        fallback_reason=None,
        n_conformers_built=n_built,
        n_conformers_kept=len(kept),
        best_dope=best_score,
        diagnostics=diagnostics,
    )
