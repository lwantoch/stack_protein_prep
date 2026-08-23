"""Emit an Amber-standard 5-stage MD input suite.

FRUTON is a model-prep tool; it does not run MD.  But every user of a
prepared model needs the same starting-point deck (min1, min2, heat,
eq, prod) with sensible defaults.  Writing those files by hand is
error-prone and produces inconsistent conventions across projects.
This module emits the canonical 5-stage suite so the user's next
command after ``tleap -f tleap.in`` is:

    pmemd.cuda -O -i min1.in -o min1.out -c system.rst7 -p system.prmtop \
               -r min1.rst7 -ref system.rst7

...and each subsequent stage chains from the previous ``.rst7``.

Stage defaults (derived from Amber Best Practices and Roe & Brooks
JCTC 2020 for equilibration protocol):

    min1 -- restrained heavy-atom SD+CG, 500+500 steps, wt=10.0
    min2 -- unrestrained SD+CG, 1000+1000 steps
    heat -- NVT, 100 K -> 300 K over 20 ps, wt-restrain Cα (wt=5.0)
    eq   -- NPT, 300 K, 100 ps, Berendsen barostat, tau_p=2.0
    prod -- NPT, 300 K, 10 ns default, Langevin gamma_ln=2.0

Rationale for each choice:

* **SHAKE on H (ntc=2, ntf=2)** everywhere -- allows 2 fs step.
* **Langevin thermostat gamma_ln=2.0 ps^-1** -- Loncharich et al.
  Biopolymers 1992; the value used in the AMBER tutorials.
* **Berendsen barostat during eq, Monte-Carlo barostat during prod**
  -- Roe & Brooks JCTC 2020 recommend MC-barostat for production to
  avoid the well-documented Berendsen artefacts.
* **cut=10.0 A, PME iwrap=1** -- standard for solvated globular
  proteins; user can tighten to 8.0 A for very large systems.
* **restraintmask uses ``!:WAT & !@H=`` for min1 and ``@CA`` for
  heat** -- protects the crystal / spliced backbone during initial
  water equilibration, then releases everything for production.

Public entry:
    write_md_input_decks(output_dir, prefix='system', ...) -> dict

License-free: pure Python string templates.  No external tool.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path


DEFAULT_TOPOLOGY_BASENAME = "system"
DEFAULT_PROD_NANOSECONDS = 10.0
DEFAULT_TIMESTEP_FS = 2.0
DEFAULT_TEMPERATURE_K = 300.0
DEFAULT_NONBONDED_CUTOFF_ANGSTROM = 10.0
DEFAULT_LANGEVIN_GAMMA_PS = 2.0


@dataclass
class MdInputDeckResult:
    output_dir: Path
    files: dict[str, Path] = field(default_factory=dict)
    warnings: list[str] = field(default_factory=list)

    def to_dict(self) -> dict:
        return {
            "output_dir": str(self.output_dir),
            "files": {stage: str(p) for stage, p in self.files.items()},
            "warnings": list(self.warnings),
        }


def _min1_template(prefix: str) -> str:
    return f"""\
{prefix}: min1 -- restrained heavy-atom minimization (500 SD + 500 CG)
 &cntrl
  imin = 1,                 ! minimization
  maxcyc = 1000,            ! total steps
  ncyc = 500,               ! SD steps before switching to CG
  ntb = 1,                  ! periodic boundary with constant V
  cut = {DEFAULT_NONBONDED_CUTOFF_ANGSTROM:.1f},
  ntr = 1,                  ! positional restraints on
  restraint_wt = 10.0,      ! kcal/mol/A^2
  restraintmask = '!:WAT & !@H=',   ! restrain all non-water heavy atoms
  ntpr = 100,               ! log every 100 steps
 /
"""


def _min2_template(prefix: str) -> str:
    return f"""\
{prefix}: min2 -- unrestrained minimization (1000 SD + 1000 CG)
 &cntrl
  imin = 1,
  maxcyc = 2000,
  ncyc = 1000,
  ntb = 1,
  cut = {DEFAULT_NONBONDED_CUTOFF_ANGSTROM:.1f},
  ntr = 0,                  ! release restraints
  ntpr = 100,
 /
"""


def _heat_template(
    prefix: str,
    duration_ps: float = 20.0,
    dt_fs: float = DEFAULT_TIMESTEP_FS,
    target_temperature_K: float = DEFAULT_TEMPERATURE_K,
) -> str:
    n_steps = int(duration_ps * 1000.0 / dt_fs)
    dt_ps = dt_fs / 1000.0
    return f"""\
{prefix}: heat -- NVT 100 K -> {target_temperature_K:.0f} K over {duration_ps:.0f} ps
 &cntrl
  imin = 0,                 ! MD
  irest = 0,                ! new run
  ntx = 1,                  ! read coordinates only
  nstlim = {n_steps},           ! {duration_ps:.0f} ps at {dt_fs:.1f} fs step
  dt = {dt_ps:.4f},
  ntb = 1,                  ! constant V, periodic
  cut = {DEFAULT_NONBONDED_CUTOFF_ANGSTROM:.1f},
  ntc = 2,                  ! SHAKE on H bonds
  ntf = 2,                  ! omit H-bond force calculations
  ntt = 3,                  ! Langevin thermostat
  gamma_ln = {DEFAULT_LANGEVIN_GAMMA_PS:.1f},
  tempi = 100.0,
  temp0 = {target_temperature_K:.1f},
  nmropt = 1,               ! use tempi/temp0 ramp defined below
  ntr = 1,                  ! weak restraint on Cα
  restraint_wt = 5.0,
  restraintmask = '@CA',
  ntpr = 500, ntwr = 5000, ntwx = 5000,
  ig = -1,                  ! random seed
 /
 &wt
  type = 'TEMP0', istep1 = 0, istep2 = {n_steps},
  value1 = 100.0, value2 = {target_temperature_K:.1f},
 /
 &wt
  type = 'END'
 /
"""


def _eq_template(
    prefix: str,
    duration_ps: float = 100.0,
    dt_fs: float = DEFAULT_TIMESTEP_FS,
    temperature_K: float = DEFAULT_TEMPERATURE_K,
) -> str:
    n_steps = int(duration_ps * 1000.0 / dt_fs)
    dt_ps = dt_fs / 1000.0
    return f"""\
{prefix}: eq -- NPT equilibration {duration_ps:.0f} ps at {temperature_K:.0f} K
 &cntrl
  imin = 0,
  irest = 1,                ! restart from heat
  ntx = 5,                  ! read v + x from rst7
  nstlim = {n_steps},
  dt = {dt_ps:.4f},
  ntb = 2,                  ! constant P
  cut = {DEFAULT_NONBONDED_CUTOFF_ANGSTROM:.1f},
  ntc = 2,
  ntf = 2,
  ntt = 3,
  gamma_ln = {DEFAULT_LANGEVIN_GAMMA_PS:.1f},
  temp0 = {temperature_K:.1f},
  ntp = 1,                  ! isotropic P scaling
  barostat = 1,             ! Berendsen for eq (fast relax)
  taup = 2.0,
  ntr = 0,
  ntpr = 500, ntwr = 5000, ntwx = 5000,
  ig = -1,
 /
"""


def _prod_template(
    prefix: str,
    duration_ns: float = DEFAULT_PROD_NANOSECONDS,
    dt_fs: float = DEFAULT_TIMESTEP_FS,
    temperature_K: float = DEFAULT_TEMPERATURE_K,
) -> str:
    n_steps = int(duration_ns * 1_000_000.0 / dt_fs)
    dt_ps = dt_fs / 1000.0
    return f"""\
{prefix}: prod -- NPT production {duration_ns:.1f} ns at {temperature_K:.0f} K
 &cntrl
  imin = 0,
  irest = 1,
  ntx = 5,
  nstlim = {n_steps},
  dt = {dt_ps:.4f},
  ntb = 2,
  cut = {DEFAULT_NONBONDED_CUTOFF_ANGSTROM:.1f},
  ntc = 2,
  ntf = 2,
  ntt = 3,
  gamma_ln = {DEFAULT_LANGEVIN_GAMMA_PS:.1f},
  temp0 = {temperature_K:.1f},
  ntp = 1,
  barostat = 2,             ! Monte-Carlo barostat (Roe & Brooks JCTC 2020)
  mcbarint = 100,
  taup = 2.0,
  ntpr = 5000, ntwr = 50000, ntwx = 5000,
  iwrap = 1,                ! wrap coords into primary unit cell
  ig = -1,
 /
"""


def _readme(prefix: str, prod_ns: float) -> str:
    return f"""\
# {prefix} — Amber production MD input suite (auto-generated by FRUTON)

Emit order (chain the ``-c`` / ``-r`` args):

    pmemd.cuda -O -i min1.in -o min1.out \\
        -p {prefix}.prmtop -c {prefix}.rst7 -r min1.rst7 -ref {prefix}.rst7

    pmemd.cuda -O -i min2.in -o min2.out \\
        -p {prefix}.prmtop -c min1.rst7   -r min2.rst7

    pmemd.cuda -O -i heat.in -o heat.out \\
        -p {prefix}.prmtop -c min2.rst7   -r heat.rst7 -x heat.nc

    pmemd.cuda -O -i eq.in   -o eq.out \\
        -p {prefix}.prmtop -c heat.rst7   -r eq.rst7   -x eq.nc

    pmemd.cuda -O -i prod.in -o prod.out \\
        -p {prefix}.prmtop -c eq.rst7     -r prod.rst7 -x prod.nc

Notes:
* Requires the ``{prefix}.prmtop`` + ``{prefix}.rst7`` produced by the
  ``tleap.in`` in this directory (FRUTON _tleap_generator).
* Prod length: {prod_ns:.1f} ns at 2 fs step ({int(prod_ns * 500_000)} steps).
  Adjust nstlim in prod.in if you need a different length.
* Barostat: Berendsen for eq (fast relax) → Monte-Carlo for prod
  (Roe & Brooks JCTC 2020 to avoid Berendsen-artefact bias).
"""


def write_md_input_decks(
    output_dir: str | Path,
    prefix: str = DEFAULT_TOPOLOGY_BASENAME,
    *,
    prod_duration_ns: float = DEFAULT_PROD_NANOSECONDS,
    heat_duration_ps: float = 20.0,
    eq_duration_ps: float = 100.0,
    dt_fs: float = DEFAULT_TIMESTEP_FS,
    target_temperature_K: float = DEFAULT_TEMPERATURE_K,
    write_readme: bool = True,
) -> MdInputDeckResult:
    """Emit the 5-file production suite (+ optional README) into ``output_dir``.

    Returns a result dataclass with the paths keyed by stage name.  Never
    raises; missing output_dir is created.  Fail-open: emits every file
    unconditionally so downstream automation can rely on their presence.
    """
    out_dir = Path(output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    warnings: list[str] = []

    # Sanity clamps
    if dt_fs > 2.5:
        warnings.append(
            f"dt_fs={dt_fs} > 2.5 fs is aggressive for SHAKE-on-H; "
            f"expect instability without HMR (hydrogen mass repartition)"
        )
    if prod_duration_ns > 1000:
        warnings.append(
            f"prod_duration_ns={prod_duration_ns} is > 1 us; verify disk / wallclock budget"
        )

    files: dict[str, Path] = {}
    files["min1"] = out_dir / "min1.in"
    files["min1"].write_text(_min1_template(prefix), encoding="utf-8")

    files["min2"] = out_dir / "min2.in"
    files["min2"].write_text(_min2_template(prefix), encoding="utf-8")

    files["heat"] = out_dir / "heat.in"
    files["heat"].write_text(
        _heat_template(
            prefix, duration_ps=heat_duration_ps,
            dt_fs=dt_fs, target_temperature_K=target_temperature_K,
        ),
        encoding="utf-8",
    )

    files["eq"] = out_dir / "eq.in"
    files["eq"].write_text(
        _eq_template(
            prefix, duration_ps=eq_duration_ps,
            dt_fs=dt_fs, temperature_K=target_temperature_K,
        ),
        encoding="utf-8",
    )

    files["prod"] = out_dir / "prod.in"
    files["prod"].write_text(
        _prod_template(
            prefix, duration_ns=prod_duration_ns,
            dt_fs=dt_fs, temperature_K=target_temperature_K,
        ),
        encoding="utf-8",
    )

    if write_readme:
        files["readme"] = out_dir / "MD_README.md"
        files["readme"].write_text(_readme(prefix, prod_duration_ns), encoding="utf-8")

    return MdInputDeckResult(
        output_dir=out_dir,
        files=files,
        warnings=warnings,
    )
