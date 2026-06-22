"""MCPB.py scaffold file writers for metal parametrization.

These functions generate the shell scripts and input files needed to run
MCPB.py → Gaussian → AMBER parametrization for a metal center. All symbols
are imported back into metall_params via
``from ._metall_params_mcpb import *``.
"""
from __future__ import annotations

from pathlib import Path

from stack_protein_preparation.metalls_check import (
    is_true_metal_atom,
    read_pdb_atom_records,
)

# Plausible spin multiplicities per (element_symbol_upper, formal_charge).
# Order matters: the first entry is the biological default that FRUTON writes
# into the MCPB.in scaffold (and submits as small_opt's charge_m_sm spin).
# Subsequent entries are documented as alternatives in the scaffold comment
# so the user can swap them when the geometry/coordination clearly demands it.
#
# Naming choice: HS (high-spin) is the biological-protein default for first-row
# open-shell metals (Mn(II), Fe(II/III), Co(II), Ni(II) octahedral) because most
# metalloenzymes operate weak-field ligand sets (His/Asp/Glu/H2O/Cys).
# Square-planar Ni(II) (LS) and porphyrinic LS Fe(II) are real exceptions —
# users should override after inspecting coordination.
#
# Closed-shell ions (d0, d10, d6-LS, etc.) have a single entry → unambiguous,
# behaves like the original _UNAMBIGUOUS_SPIN.
_PLAUSIBLE_SPINS: dict[tuple[str, int], tuple[int, ...]] = {
    # d10 closed-shell — unambiguous
    ("ZN", 2): (1,),
    ("CD", 2): (1,),
    ("HG", 2): (1,),
    ("CU", 1): (1,),
    # d9 — Cu(II): always one unpaired electron
    ("CU", 2): (2,),
    # d8 — Ni(II): octahedral HS triplet, square-planar LS singlet
    ("NI", 2): (3, 1),
    # d7 — Co(II): HS quartet dominant biologically, LS doublet rarer
    ("CO", 2): (4, 2),
    # d6 — Fe(II) HS quintet (deoxy-heme, non-heme), LS singlet (CO/CN/strong field)
    ("FE", 2): (5, 1),
    # d6 — Co(III) typically LS singlet (strong-field d6), HS quintet rare
    ("CO", 3): (1, 5),
    # d5 — Fe(III) HS sextet (most protein sites), LS doublet (cytochromes, low-spin heme)
    ("FE", 3): (6, 2),
    # d5 — Mn(II) HS sextet (almost always), LS doublet only in very strong fields
    ("MN", 2): (6, 2),
    # d4 — Mn(III) HS quintet, LS triplet possible
    ("MN", 3): (5, 3),
    # d4 — Fe(IV) HS quintet, LS triplet in compound I/II of peroxidases (LS triplet common)
    ("FE", 4): (3, 5),
    # d3 — Cr(III) quartet (essentially fixed)
    ("CR", 3): (4,),
    # d3 — Mn(IV) quartet
    ("MN", 4): (4,),
    # d2 — V(III) triplet
    ("V", 3): (3,),
    # d1 — Ti(III) doublet
    ("TI", 3): (2,),
    # d1 — V(IV) doublet
    ("V", 4): (2,),
    # d0 — closed shell
    ("TI", 4): (1,),
    ("V", 5): (1,),
    ("ZR", 4): (1,),
    ("MO", 6): (1,),
    # Group 1/2 ions — closed shell
    ("MG", 2): (1,),
    ("CA", 2): (1,),
    ("NA", 1): (1,),
    ("K", 1): (1,),
}


def _get_spin_candidates(element: str, formal_charge: int | None) -> tuple[int, ...]:
    """Return the ordered tuple of plausible multiplicities for the given metal.

    Empty tuple means FRUTON has no opinion — the .in scaffold will fall back
    to a TODO_spin placeholder requiring the user to set it by hand.
    """
    if formal_charge is None:
        return ()
    return _PLAUSIBLE_SPINS.get((element.upper(), formal_charge), ())


# Back-compat shim: external callers (metall_params.py, _metall_params_helpers.py)
# still query ``_UNAMBIGUOUS_SPIN.get(element)`` and expect either an int or None.
# We synthesise a dict-like view that returns the default spin (first candidate)
# only when the element-charge pair is unambiguous (single candidate); otherwise
# returns None so the caller writes a clearly-marked alternative-spin scaffold.
class _UnambiguousSpinView:
    """Element-only view: returns the default multiplicity only when the
    metal has exactly one plausible spin across its plausible oxidation states.

    This preserves the original ``_UNAMBIGUOUS_SPIN.get("ZN") -> 1`` semantics
    for d10 closed-shell ions while explicitly returning ``None`` for ions like
    Fe/Mn/Co/Ni where the chemistry is genuinely ambiguous and the caller must
    handle alternative-spin enumeration (see Task #7 phase 2).
    """

    def get(self, element: str, default: int | None = None) -> int | None:
        candidates_by_ox = [
            spins
            for (el, _ox), spins in _PLAUSIBLE_SPINS.items()
            if el == element.upper()
        ]
        if not candidates_by_ox:
            return default
        unique_spins = {s for tup in candidates_by_ox for s in tup}
        if len(unique_spins) == 1:
            return next(iter(unique_spins))
        return default


_UNAMBIGUOUS_SPIN = _UnambiguousSpinView()

# ---------------------------------------------------------------------------
# MCPB scaffold writers (03_mcpb)
# ---------------------------------------------------------------------------

def _write_mcpb_in(
    output_path: Path,
    pdb_id: str,
    element: str,
    chain: str,
    resseq: int,
    max_contact_dist: float,
    ion_serial: int | None,
    ion_mol2file: str,
    formal_charge: int | None,
    spin_guess: int | None,
    naa_mol2files: list[str] | None = None,
    spin_alternatives: tuple[int, ...] = (),
) -> None:
    """Write a first-pass MCPB.py ``.in`` scaffold.

    Fields that FRUTON can resolve automatically (serial number, mol2 filename,
    metal formal charge, spin for d10 ions) are filled in directly.  Fields
    that require manual chemistry input (spin for open-shell metals, ligand
    partial charges) are written as ``TODO`` with an explanatory comment.

    ``spin_alternatives`` is the ordered tuple of additional plausible
    multiplicities (excluding ``spin_guess``) for this (element, oxidation)
    pair.  When present, the scaffold comment lists them so the user can swap
    in HS/LS without recomputing the table by hand.
    """
    group_name = f"{pdb_id}_{element}_{chain}_{resseq}"
    cut_off = max_contact_dist + 0.5

    # ion_ids — resolved from actual serial in merged PDB
    if ion_serial is not None:
        ion_ids_line = f"ion_ids {ion_serial}"
        ion_ids_comment = f"# Serial of {element} in {pdb_id}_mcpb.pdb (auto-resolved by FRUTON)."
    else:
        ion_ids_line = "ion_ids TODO_metal_serial_in_mcpb_pdb"
        ion_ids_comment = f"# TODO: Serial number of {element} in {pdb_id}_mcpb.pdb."

    # ion_mol2files — generated by FRUTON when charge is known
    if not ion_mol2file.startswith("TODO"):
        mol2_comment = "# mol2 file for the metal ion (auto-generated by FRUTON from known formal charge)."
    else:
        mol2_comment = "# TODO: Provide mol2 file for the metal ion and any non-standard ligands."

    # charge_m_sm / charge_m_lm
    if formal_charge is not None and spin_guess is not None:
        charge_sm_line = f"charge_m_sm {formal_charge} {spin_guess}"
        charge_lm_line = f"charge_m_lm {formal_charge} {spin_guess}"
        if spin_alternatives:
            alt_pretty = ", ".join(str(s) for s in spin_alternatives)
            charge_comment = (
                f"# Metal formal charge {formal_charge}+ and spin {spin_guess} "
                f"auto-filled (biological default for {element}{formal_charge}+).\n"
                f"# Alternative plausible multiplicities for this ion: {alt_pretty}.\n"
                f"# If the coordination geometry / ligand-field clearly favours an\n"
                f"# alternative (e.g. square-planar Ni(II) → 1, low-spin heme Fe(II) → 1),\n"
                f"# substitute it in the charge_m_sm / charge_m_lm lines below.\n"
                "# NOTE: Add net charge of any ionised ligand residues in the QM region."
            )
        else:
            charge_comment = (
                f"# Metal formal charge {formal_charge}+ and spin {spin_guess} "
                f"(unambiguous for {element}, auto-filled).\n"
                "# NOTE: Add net charge of any ionised ligand residues in the QM region."
            )
    elif formal_charge is not None:
        charge_sm_line = f"charge_m_sm {formal_charge} TODO_spin"
        charge_lm_line = f"charge_m_lm {formal_charge} TODO_spin"
        charge_comment = (
            f"# Metal formal charge {formal_charge}+ auto-filled.\n"
            f"# TODO: Set spin multiplicity — depends on oxidation state and\n"
            f"#       ligand-field splitting of {element}. "
            "1=singlet, 2=doublet, 3=triplet, …\n"
            "# NOTE: Add net charge of any ionised ligand residues in the QM region."
        )
    else:
        charge_sm_line = "charge_m_sm TODO_charge TODO_spin"
        charge_lm_line = "charge_m_lm TODO_charge TODO_spin"
        charge_comment = (
            "# TODO: Set total charge and spin multiplicity for the QM region.\n"
            "# Spin: 1=singlet, 2=doublet, 3=triplet, …"
        )

    # naa_mol2files — non-standard organic ligands in the coordination shell
    naa_block = ""
    if naa_mol2files:
        naa_list = " ".join(naa_mol2files)
        naa_block = (
            "\n# mol2 file(s) for non-standard organic ligand(s) in the QM region.\n"
            f"naa_mol2files {naa_list}\n"
        )

    content = f"""\
# MCPB.py input file — first-pass scaffold generated by FRUTON
# Review carefully before running. Remaining TODOs are marked explicitly.

original_pdb {pdb_id}_mcpb.pdb
group_name {group_name}

# cut_off derived from max donor-metal distance + 0.5 Å buffer. Verify.
cut_off {cut_off:.1f}

{ion_ids_comment}
{ion_ids_line}

{mol2_comment}
ion_mol2files {ion_mol2file}
{naa_block}
# 1 = optimise only H positions in the large model (recommended).
large_opt 1

force_field ff99SB

# Verify Gaussian version on your system (g09 or g16).
software_version g16

water_model tip3p

{charge_comment}
{charge_sm_line}
{charge_lm_line}
"""
    output_path.write_text(content, encoding="utf-8")


# ---------------------------------------------------------------------------
# MCPB input/output consistency audit
# ---------------------------------------------------------------------------

def _parse_charge_mult_from_in_file(in_path: Path) -> tuple[int, int] | None:
    """Return (charge, multiplicity) from ``charge_m_sm <C> <M>`` in an MCPB .in
    file.  Returns None if the line is absent or marked TODO.
    """
    if not in_path.is_file():
        return None
    for line in in_path.read_text(errors="replace").splitlines():
        stripped = line.strip()
        if not stripped.startswith("charge_m_sm"):
            continue
        parts = stripped.split()
        if len(parts) < 3:
            return None
        try:
            return int(parts[1]), int(parts[2])
        except ValueError:
            return None
    return None


def _parse_charge_mult_from_com_file(com_path: Path) -> tuple[int, int] | None:
    """Return (charge, multiplicity) from a Gaussian ``.com`` input.

    The line follows the route + blank + title + blank block.  We scan the first
    ~14 lines and return the first two-integer line found.
    """
    if not com_path.is_file():
        return None
    import re as _re
    int_pair = _re.compile(r"^(-?\d+)\s+(\d+)\s*$")
    for idx, line in enumerate(com_path.read_text(errors="replace").splitlines()):
        if idx > 14:
            break
        m = int_pair.match(line.strip())
        if m:
            return int(m.group(1)), int(m.group(2))
    return None


def audit_mcpb_input_consistency(mcpb_dir: Path) -> list[dict[str, object]]:
    """Audit charge/multiplicity consistency between an MCPB ``.in`` and its
    Gaussian ``.com`` files.

    Returns a list of issue dicts.  An empty list means every .com's leading
    ``<charge> <mult>`` line matches what the .in's ``charge_m_sm`` prescribed.
    Detects:

    * .com hand-edited after MCPB.py emitted it (most common cause when a user
      tweaks a route to retry stuck convergence and accidentally changes the
      electronic state).
    * .com missing the charge/mult line entirely (Gaussian will crash on read).
    * .in still carrying a TODO_charge / TODO_spin placeholder.

    Re-running ``MCPB.py -s 1`` regenerates .com from .in and clears drift.
    """
    issues: list[dict[str, object]] = []
    in_files = sorted(mcpb_dir.glob("*.in"))
    if not in_files:
        return issues
    in_path = in_files[0]
    expected = _parse_charge_mult_from_in_file(in_path)
    if expected is None:
        issues.append(
            {
                "kind": "in_unparseable_or_todo",
                "in_path": str(in_path),
                "message": (
                    "Could not parse charge_m_sm from the .in (TODO "
                    "placeholder or missing line).  MCPB.py needs explicit "
                    "charge and multiplicity to emit Gaussian inputs."
                ),
            }
        )
        return issues
    expected_charge, expected_mult = expected
    gauss_dir = mcpb_dir / "step02_gaussian"
    if not gauss_dir.is_dir():
        gauss_dir = mcpb_dir / "step01_gen_inputs"
    for com_path in sorted(gauss_dir.glob("*.com")):
        # small_fc reads small_opt.chk via Geom=AllCheckpoint and legitimately
        # omits an explicit charge/mult line; skip.
        if com_path.name.endswith("_small_fc.com"):
            continue
        observed = _parse_charge_mult_from_com_file(com_path)
        if observed is None:
            issues.append(
                {
                    "kind": "com_missing_charge_mult",
                    "com_path": str(com_path),
                    "expected_charge": expected_charge,
                    "expected_mult": expected_mult,
                    "message": (
                        f"No charge/multiplicity line found in {com_path.name}; "
                        "Gaussian will fail to read the input."
                    ),
                }
            )
            continue
        observed_charge, observed_mult = observed
        if observed_charge != expected_charge or observed_mult != expected_mult:
            issues.append(
                {
                    "kind": "com_in_charge_mult_mismatch",
                    "com_path": str(com_path),
                    "in_path": str(in_path),
                    "expected_charge": expected_charge,
                    "expected_mult": expected_mult,
                    "observed_charge": observed_charge,
                    "observed_mult": observed_mult,
                    "message": (
                        f"{com_path.name} has charge={observed_charge}, "
                        f"mult={observed_mult} but the .in prescribes "
                        f"{expected_charge}/{expected_mult}.  Either the .com "
                        "was hand-edited after MCPB.py -s 1, or MCPB.py "
                        "auto-derived a different charge from the QM-region "
                        "residue composition.  Re-run MCPB.py -s 1 to "
                        "regenerate .com from .in, or update .in to match the "
                        "intended electronic state."
                    ),
                }
            )
    return issues


def _write_commands_sh(
    output_path: Path,
    pdb_id: str,
    group_name: str,
    mol2_filename: str,
    naa_mol2files: list[str] | None = None,
) -> None:
    """Write the step-1 script: creates step01_gen_inputs/, runs MCPB.py -s 1,
    then stages .com files into step02_gaussian/."""
    # Build antechamber generation + copy block for naa mol2 files (if any).
    # Each mol2 is generated on-the-fly if it doesn't already exist so that
    # commands.sh can be run directly in the AmberTools environment on HPC
    # without needing a separate preparation step.
    naa_gen_and_copy_lines = ""
    if naa_mol2files:
        lines: list[str] = [
            "\n# -- Generate mol2 for non-standard ligand(s) if not already present --------",
            "# IMPORTANT: crystal ligand PDB records almost never include hydrogens.",
            "# antechamber does NOT add hydrogens itself; running it on a heavy-atom-only",
            "# PDB produces a mol2 with the missing-H valence structure, and MCPB.py then",
            "# computes an incorrect electron count (typically odd, forcing a wrong-state",
            "# doublet Berny optimisation that does not converge in Gaussian).  We add",
            "# hydrogens with OpenBabel first when available, then run antechamber.",
        ]
        for mol2_f in naa_mol2files:
            resname = mol2_f.replace(".mol2", "")
            lines.append(
                f"if [ ! -f {mol2_f} ]; then\n"
                f"    # TODO: verify net charge (default 0) — update -nc if ligand is ionised\n"
                f"    if command -v obabel >/dev/null 2>&1; then\n"
                f"        obabel {resname}.pdb -O {resname}_H.pdb -h\n"
                f"        antechamber -i {resname}_H.pdb -fi pdb -o {mol2_f} -fo mol2"
                f" -c bcc -nc 0 -at gaff2 -s 2 -pf y\n"
                f"    else\n"
                f"        echo \"WARNING: obabel not in PATH; running antechamber on H-less PDB.\" >&2\n"
                f"        echo \"         Verify the resulting mol2 has the expected hydrogen count.\" >&2\n"
                f"        antechamber -i {resname}.pdb -fi pdb -o {mol2_f} -fo mol2"
                f" -c bcc -nc 0 -at gaff2 -s 2 -pf y\n"
                f"    fi\n"
                f"fi"
            )
        copy_str = " ".join(naa_mol2files)
        lines.append(f"cp {copy_str} step01_gen_inputs/")
        naa_gen_and_copy_lines = "\n".join(lines) + "\n"

    content = f"""\
#!/usr/bin/env bash
# =============================================================================
# MCPB.py — Step 1  |  {pdb_id}  |  site {group_name}
# =============================================================================
# Creates step01_gen_inputs/, copies input files, runs MCPB.py -s 1.
# Stages the three Gaussian .com files into step02_gaussian/.
#
# After this script:
#   1. Review .com files in step02_gaussian/.
#   2. Submit Gaussian jobs:  bash submit_gaussian.sh
#   3. Once all Gaussian jobs complete:  bash run_after_gaussian.sh
# =============================================================================
set -euo pipefail
{naa_gen_and_copy_lines}
# -- Create step01 directory and copy all required inputs --------------------
mkdir -p step01_gen_inputs
cp {pdb_id}.in {pdb_id}_mcpb.pdb {mol2_filename} step01_gen_inputs/
# -- Run MCPB.py step 1 (generates Gaussian input files and model PDBs) ------
cd step01_gen_inputs
MCPB.py -i {pdb_id}.in -s 1
cd ..

# -- Stage Gaussian input files into step02_gaussian/ ------------------------
mkdir -p step02_gaussian
cp step01_gen_inputs/{group_name}_small_opt.com step02_gaussian/
cp step01_gen_inputs/{group_name}_small_fc.com  step02_gaussian/
cp step01_gen_inputs/{group_name}_large_mk.com  step02_gaussian/

echo ""
echo "Step 1 done.  Gaussian inputs staged in step02_gaussian/:"
echo "  {group_name}_small_opt.com   geometry optimisation (small model)"
echo "  {group_name}_small_fc.com    force constants / Hessian (small model)"
echo "  {group_name}_large_mk.com    RESP charges (large model)"
echo ""
echo "Review each .com file, then run:  bash submit_gaussian.sh"
"""
    output_path.write_text(content, encoding="utf-8")
    output_path.chmod(0o755)


def _write_submit_gaussian_sh(output_path: Path, group_name: str) -> None:
    """Write a CESGA-compatible MCPB Gaussian submission script.

    MCPB requires three Gaussian jobs in a specific order:
      - small_opt  (geometry optimisation → writes .chk)
      - large_mk   (MK ESP, independent of small_opt)
      - small_fc   (frequency/force-constants → reads small_opt.chk)

    small_opt and large_mk are submitted immediately; small_fc is submitted
    with ``--dependency=afterok:<small_opt_job_id>`` so SLURM holds it until
    small_opt succeeds.  This avoids the "submit all, let one fail, resubmit"
    anti-pattern.
    """
    content = f"""\
#!/usr/bin/env bash
# =============================================================================
# Gaussian HPC submission — {group_name}
# =============================================================================
# Submits the three MCPB Gaussian jobs with the correct SLURM dependency:
#   small_opt  ──┐ (parallel)
#   large_mk     │
#   small_fc  ←──┘ afterok:small_opt  (needs small_opt.chk)
#
# Run ONLY after bash commands.sh has completed successfully.
# After all jobs complete:  bash run_after_gaussian.sh
# =============================================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${{BASH_SOURCE[0]}}")" && pwd)"
# Locate the central CESGA submission script.
# Priority: $FRUTON_SCRIPTS env var (set this when working from $LUSTRE
# where the git repo is not present), then git-based discovery, then
# directory walk-up.
if [[ -n "${{FRUTON_SCRIPTS:-}}" && -f "${{FRUTON_SCRIPTS}}/submit_gaussian.sh" ]]; then
    CENTRAL="${{FRUTON_SCRIPTS}}/submit_gaussian.sh"
elif REPO_ROOT="$(git -C "$SCRIPT_DIR" rev-parse --show-toplevel 2>/dev/null)"; then
    CENTRAL="${{REPO_ROOT}}/scripts/CESGA_SLURM/submit_gaussian.sh"
else
    CENTRAL=""
    DIR="$SCRIPT_DIR"
    for _ in 1 2 3 4 5 6 7 8; do
        if [[ -f "${{DIR}}/scripts/CESGA_SLURM/submit_gaussian.sh" ]]; then
            CENTRAL="${{DIR}}/scripts/CESGA_SLURM/submit_gaussian.sh"
            break
        fi
        DIR="$(dirname "$DIR")"
    done
fi
[[ -n "${{CENTRAL:-}}" && -f "$CENTRAL" ]] || {{ echo "ERROR: central submit_gaussian.sh not found. Set FRUTON_SCRIPTS=/path/to/scripts/CESGA_SLURM" >&2; exit 2; }}

GAUSS_DIR="$SCRIPT_DIR/step02_gaussian"
SLURM_OPTS=(--partition short --time 06:00:00 --mem 32G --cpus 16)

# 1. small_opt — geometry optimisation (produces .chk needed by small_fc)
OPT_OUT=$(bash "$CENTRAL" "$GAUSS_DIR/{group_name}_small_opt.com" "${{SLURM_OPTS[@]}}" "$@")
echo "$OPT_OUT"
OPT_JOB=$(echo "$OPT_OUT" | grep -oP '(?<=Submitted batch job )\\d+' || true)
if [[ -z "$OPT_JOB" ]]; then
    echo "WARNING: could not parse small_opt job ID — small_fc will be submitted without dependency" >&2
fi

# 2. large_mk — ESP calculation, independent of small_opt
bash "$CENTRAL" "$GAUSS_DIR/{group_name}_large_mk.com" "${{SLURM_OPTS[@]}}" "$@"

# 3. small_fc — force constants, depends on small_opt.chk
if [[ -n "$OPT_JOB" ]]; then
    bash "$CENTRAL" "$GAUSS_DIR/{group_name}_small_fc.com" "${{SLURM_OPTS[@]}}" \
        --dependency "afterok:$OPT_JOB" "$@"
else
    bash "$CENTRAL" "$GAUSS_DIR/{group_name}_small_fc.com" "${{SLURM_OPTS[@]}}" "$@"
fi

echo ""
echo "Submitted: small_opt (job $OPT_JOB), large_mk, small_fc (afterok:$OPT_JOB)"
echo "Once all complete:  bash run_after_gaussian.sh"
"""
    output_path.write_text(content, encoding="utf-8")
    output_path.chmod(0o755)


def _write_run_after_gaussian_sh(output_path: Path, pdb_id: str, group_name: str) -> None:
    """Write steps 2-4 script: creates step03_amber_params/, merges step01 +
    step02 outputs, runs MCPB.py -s 2,3,4 and tleap."""
    content = f"""\
#!/usr/bin/env bash
# =============================================================================
# MCPB.py — Steps 2-4  |  {pdb_id}  |  site {group_name}
# =============================================================================
# Creates step03_amber_params/, assembles all required inputs, then runs
# MCPB.py -s 2, -s 3, -s 4, and tleap.
#
# Run ONLY after all three Gaussian jobs in step02_gaussian/ have completed.
#
# Required files in step02_gaussian/ after Gaussian:
#   {group_name}_small_opt.fchk   (from: formchk _small_opt.chk)
#   {group_name}_small_fc.log
#   {group_name}_large_mk.log
#
# Final outputs in step03_amber_params/:
#   {group_name}.frcmod            AMBER bond/angle force constants
#   {group_name}_addmed.mol2       AMBER partial charges (RESP)
#   {group_name}_tleap.in / {group_name}_tleap.out
# =============================================================================
set -euo pipefail

# -- Create step03 directory and assemble all required inputs ----------------
mkdir -p step03_amber_params

# All step01 outputs (model PDBs, .in, mol2, MCPB intermediate files)
cp step01_gen_inputs/* step03_amber_params/

# Gaussian outputs from step02
cp step02_gaussian/{group_name}_small_opt.fchk step03_amber_params/
cp step02_gaussian/{group_name}_small_fc.log   step03_amber_params/
cp step02_gaussian/{group_name}_large_mk.log   step03_amber_params/

# -- Run MCPB.py steps 2-4 and tleap ----------------------------------------
cd step03_amber_params

# Step 2: Derive bond/angle force constants (Seminario method).
MCPB.py -i {pdb_id}.in -s 2

# Step 3: Fit RESP charges.
MCPB.py -i {pdb_id}.in -s 3

# Step 4: Generate LEaP topology input.
MCPB.py -i {pdb_id}.in -s 4
tleap -s -f {group_name}_tleap.in > {group_name}_tleap.out

echo ""
echo "Parametrization complete.  Key outputs in step03_amber_params/:"
echo "  {group_name}.frcmod"
echo "  {group_name}_addmed.mol2"
echo "  {group_name}_tleap.in / {group_name}_tleap.out"
"""
    output_path.write_text(content, encoding="utf-8")
    output_path.chmod(0o755)


def _write_metal_mol2(
    output_path: Path,
    element: str,
    formal_charge: int,
    x: float,
    y: float,
    z: float,
) -> None:
    """Write a minimal TRIPOS mol2 for a bare transition-metal ion.

    The atom type uses the AMBER/GAFF convention (element symbol in title case,
    e.g. ``Zn``, ``Fe``).  The formal charge is taken from
    :data:`METAL_RESIDUE_IDENTITY` and is written as a float in the charge
    column so that MCPB.py can read it directly.
    """
    atom_type = element.capitalize()
    charge_int = formal_charge
    # subst_name must match the PDB residue name (e.g. "ZN") so that MCPB.py's
    # chargedict lookup (keyed by residue name) finds the charge from this mol2.
    subst_name = element
    content = (
        "@<TRIPOS>MOLECULE\n"
        f"{element}\n"
        " 1 0 0 0 0\n"
        "SMALL\n"
        "NO_CHARGES\n"
        "\n"
        "\n"
        "@<TRIPOS>ATOM\n"
        f"      1 {element:<8s}{x:10.4f}{y:10.4f}{z:10.4f} "
        f"{atom_type:<8s}  1  {subst_name:<8s}{float(formal_charge):10.4f}\n"
        "@<TRIPOS>BOND\n"
    )
    output_path.write_text(content, encoding="utf-8")


def _find_metal_serial_in_pdb(pdb_path: Path, element: str) -> int | None:
    """Return the serial number of the first metal atom matching *element*."""
    try:
        for rec in read_pdb_atom_records(pdb_path, source_role="serial_lookup"):
            if is_true_metal_atom(rec) and rec.element == element:
                return rec.serial
    except Exception:
        pass
    return None


def _write_mcpb_readme(output_path: Path, pdb_id: str, group_name: str) -> None:
    """Write a README.md explaining the full MCPB multi-step workflow."""
    content = f"""\
# MCPB Parametrization — {pdb_id}  |  site {group_name}

## Overview

MCPB.py (Metal Center Parameter Builder) derives bonded-model AMBER parameters
for transition-metal centers using quantum-chemical calculations.

Reference: doi:10.1021/acs.jcim.5b00674

## Directory layout

```
03_mcpb/                         ← run all scripts from here
├── {pdb_id}.in
├── {pdb_id}_mcpb.pdb
├── *.mol2
├── commands.sh           Step 1
├── submit_gaussian.sh    Gaussian HPC submission
├── run_after_gaussian.sh Steps 2-4
│
├── step01_gen_inputs/    created by commands.sh
│   ├── {group_name}_small_opt.com
│   ├── {group_name}_small_fc.com
│   ├── {group_name}_large_mk.com
│   └── {group_name}_small.pdb / _large.pdb  (model PDBs)
│
├── step02_gaussian/      .com files staged here; Gaussian logs written here
│   ├── {group_name}_small_opt.log / .fchk
│   ├── {group_name}_small_fc.log
│   └── {group_name}_large_mk.log
│
└── step03_amber_params/  created by run_after_gaussian.sh; final outputs here
    ├── {group_name}.frcmod
    ├── {group_name}_addmed.mol2
    └── {group_name}_tleap.in / _tleap.out
```

## Workflow (3 scripts, 4 MCPB.py invocations)

```
bash commands.sh
  → mkdir step01_gen_inputs/, run MCPB.py -s 1
  → mkdir step02_gaussian/, copy .com files

bash submit_gaussian.sh
  → sbatch all three Gaussian jobs from step02_gaussian/

bash run_after_gaussian.sh
  → mkdir step03_amber_params/, merge step01 + step02 outputs
  → MCPB.py -s 2, -s 3, -s 4, tleap
```

## Before running Step 1

Review ``{pdb_id}.in``.  Fields auto-filled by FRUTON are marked in comments.
Fields that require manual chemistry input carry ``TODO`` placeholders:

- ``software_version`` — Gaussian version on your cluster (g09 or g16)
- ``charge_m_sm / charge_m_lm`` — total charge + spin multiplicity of the QM
  region; spin is unambiguous only for closed-shell d10 ions (Zn²⁺, Cu⁺)

## Gaussian

Gaussian is a licensed external program; FRUTON does not run it automatically.
``submit_gaussian.sh`` provides a SLURM template — adjust ``--cpus``,
``--gpus``, ``--mem``, and ``--time`` to match your cluster's policies.

All three Gaussian jobs are independent and can be submitted simultaneously.
Steps 2-4 require all three to finish before starting.
"""
    output_path.write_text(content, encoding="utf-8")
