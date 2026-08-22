"""Junction-only energy minimization after AF-splice.

After ``splice_af_gaps_into_crystal`` inserts AF-modelled residues into the
crystal template, the peptide bond at the crystal/AF boundary sometimes carries
residual strain (bad omega, sub-optimal C-N distance).  A restrained
minimization that lets only ``gap ± N`` residues move -- with every other
heavy atom held by a stiff harmonic restraint -- absorbs the strain into the
one region we intend to change while leaving the crystal untouched.

Rationale: Croll et al. IUCr Acta D 2025 (10.1107/S2059798325008496) shows that
an unrestrained minimization silently masks bad junctions -- the bond geometry
becomes ideal while Ramachandran remains outlier.  A boundary-only relaxation
gives the strain somewhere to go without disturbing MolProbity-favoured
regions of the crystal.

We use OpenMM with amber14/ff14SB + implicit solvent (GBn2) so this runs
inside the FRUTON pixi environment without requiring the CESGA AmberTools
module load.  Non-standard residues (phosphos, ALY, CAS, ...) are skipped
gracefully: if force-field parametrisation fails the input is returned
unchanged and the failure reason is persisted for downstream inspection.
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

from Bio.PDB import PDBIO, PDBParser
from Bio.PDB.Chain import Chain as _BioChain


DEFAULT_FLANK_RESIDUES = 2
# Peptide-bond C-N distance beyond which OpenMM's `createStandardBonds` will
# NOT infer a bond (its internal threshold is ~2.0 A).  We split the chain
# there so both residues become chain termini.
CHAIN_BREAK_C_N_ANGSTROM = 2.0
# CA-CA fallback threshold when the C or N atom is missing (very rare after
# splice+cleanup, but be defensive).  Peptide-bond CA-CA is 3.80 A nominal.
CHAIN_BREAK_CA_CA_ANGSTROM = 4.5
# 10 kcal/mol/A^2 in kJ/mol/nm^2 (positional restraint strength).  Matches the
# AlphaFold reference relax used by DeepMind's amber_minimize.py.
DEFAULT_RESTRAINT_K_KJ_PER_MOL_NM2 = 10.0 * 4.184 * 100.0
# LocalEnergyMinimizer tolerance in kJ/mol/nm.  AF2 uses ~2.39.
DEFAULT_MIN_TOLERANCE_KJ_PER_MOL_NM = 2.39
DEFAULT_MAX_ITERATIONS = 2500


@dataclass
class JunctionRelaxResult:
    output_pdb_path: Path
    ran: bool
    fallback_reason: str | None
    junction_rmsd_angstrom: float | None
    restrained_rmsd_angstrom: float | None
    final_energy_kj_per_mol: float | None
    n_free_atoms: int
    n_restrained_atoms: int

    def to_dict(self) -> dict:
        return {
            "ran": self.ran,
            "fallback_reason": self.fallback_reason,
            "junction_rmsd_angstrom": self.junction_rmsd_angstrom,
            "restrained_rmsd_angstrom": self.restrained_rmsd_angstrom,
            "final_energy_kj_per_mol": self.final_energy_kj_per_mol,
            "n_free_atoms": self.n_free_atoms,
            "n_restrained_atoms": self.n_restrained_atoms,
        }


def _free_residue_ids(
    gap_ranges: Iterable[tuple[int, int]], flank: int
) -> set[int]:
    free: set[int] = set()
    for lo, hi in gap_ranges:
        for r in range(lo - flank, hi + flank + 1):
            free.add(r)
    return free


def _split_chains_at_ca_breaks(
    input_pdb: Path, output_pdb: Path
) -> tuple[Path, dict]:
    """Re-write ``input_pdb`` so consecutive protein residues whose CA atoms
    are further apart than ``CHAIN_BREAK_CA_CA_ANGSTROM`` end up in separate
    Bio.PDB chains.  This gives OpenMM proper TER records so it can recognise
    the residues on either side as chain termini and place N-terminal /
    C-terminal capping atoms correctly.

    Returns (output_path, {original_chain_id: [split_chain_ids...]}).
    HETATM records and their chains are copied through unchanged.
    """
    import math as _math
    s = PDBParser(QUIET=True).get_structure("s", str(input_pdb))
    split_map: dict[str, list[str]] = {}
    for model in s:
        for chain in list(model):
            residues = list(chain)
            protein = [r for r in residues if not r.id[0].strip()]
            if len(protein) < 2:
                continue
            breaks: list[int] = []
            for i in range(1, len(protein)):
                r_prev = protein[i - 1]
                r_curr = protein[i]
                is_break = False
                if "C" in r_prev and "N" in r_curr:
                    d = r_prev["C"].coord - r_curr["N"].coord
                    if _math.sqrt(float(d[0]) ** 2 + float(d[1]) ** 2 + float(d[2]) ** 2) > CHAIN_BREAK_C_N_ANGSTROM:
                        is_break = True
                elif "CA" in r_prev and "CA" in r_curr:
                    d = r_prev["CA"].coord - r_curr["CA"].coord
                    if _math.sqrt(float(d[0]) ** 2 + float(d[1]) ** 2 + float(d[2]) ** 2) > CHAIN_BREAK_CA_CA_ANGSTROM:
                        is_break = True
                if is_break:
                    breaks.append(i)
            if not breaks:
                continue
            segments: list[list] = []
            prev = 0
            for b in breaks:
                segments.append(protein[prev:b])
                prev = b
            segments.append(protein[prev:])
            hetatm_list = [r for r in residues if r.id[0].strip()]
            model.detach_child(chain.id)
            # PDB format requires single-char chain IDs.  Use the original ID
            # for the first segment; then draw from lowercase letters + digits
            # for subsequent segments, skipping IDs already used by other
            # chains in the model.
            existing_ids = {c.id for c in model}
            single_char_pool = (
                [chr(c) for c in range(ord("a"), ord("z") + 1)]
                + list("0123456789")
                + [chr(c) for c in range(ord("A"), ord("Z") + 1)]
            )
            new_ids: list[str] = []
            def _pick_id() -> str:
                for cand in single_char_pool:
                    if cand not in existing_ids and cand not in new_ids:
                        return cand
                raise ValueError("no single-char chain id available")
            for idx, seg in enumerate(segments):
                new_cid = chain.id if idx == 0 else _pick_id()
                new_chain = _BioChain(new_cid)
                for r in seg:
                    r.detach_parent()
                    new_chain.add(r)
                # HETATMs are attached to the LAST segment so they keep their
                # relationship with the C-terminal environment.
                if idx == len(segments) - 1:
                    for r in hetatm_list:
                        r.detach_parent()
                        new_chain.add(r)
                model.add(new_chain)
                new_ids.append(new_cid)
            split_map[chain.id] = new_ids
    io = PDBIO()
    io.set_structure(s)
    io.save(str(output_pdb))
    return output_pdb, split_map


def _copy_original_numbering(
    template_pdb: Path, relaxed_pdb: Path, output_pdb: Path
) -> None:
    """Overwrite residue numbers in ``relaxed_pdb`` with those from
    ``template_pdb`` (in file order).  OpenMM's PDBFile writer preserves
    numbering from the input topology; this is a safety net for the rare case
    where a chain-ordering shuffle occurs.
    """
    tpl = PDBParser(QUIET=True).get_structure("t", str(template_pdb))
    rlx = PDBParser(QUIET=True).get_structure("r", str(relaxed_pdb))
    tpl_res = [r for r in tpl.get_residues()]
    rlx_res = [r for r in rlx.get_residues()]
    for t, r in zip(tpl_res, rlx_res):
        r.id = t.id
    io = PDBIO()
    io.set_structure(rlx)
    io.save(str(output_pdb))


def _rmsd_angstrom(coords_before, coords_after) -> float:
    import math
    n = len(coords_before)
    if n == 0:
        return 0.0
    ss = 0.0
    for (x0, y0, z0), (x1, y1, z1) in zip(coords_before, coords_after):
        ss += (x1 - x0) ** 2 + (y1 - y0) ** 2 + (z1 - z0) ** 2
    # coords in nm inside OpenMM; convert nm^2 -> A^2 by *100
    return math.sqrt(ss / n) * 10.0


def relax_junctions(
    input_pdb_path: str | Path,
    output_pdb_path: str | Path,
    gap_ranges: Iterable[tuple[int, int]],
    flank_residues: int = DEFAULT_FLANK_RESIDUES,
    restraint_k: float = DEFAULT_RESTRAINT_K_KJ_PER_MOL_NM2,
    tolerance: float = DEFAULT_MIN_TOLERANCE_KJ_PER_MOL_NM,
    max_iterations: int = DEFAULT_MAX_ITERATIONS,
) -> JunctionRelaxResult:
    """Restrained OpenMM minimization holding everything but ``gap ± N``
    heavy atoms rigid.

    ``gap_ranges`` is an iterable of ``(first_resnum, last_resnum)`` pairs
    naming residues that were AF-inserted.  Flanking crystal residues within
    ``flank_residues`` on each side are also allowed to move so the strained
    peptide bonds at the junction can relax.

    Chain identity is ignored when applying the free-set (residue numbers are
    globally unique within the FRUTON post-splice PDB per chain-wise numbering
    scheme).  If a resnum collides across chains, both are freed -- safer than
    freezing a real junction.
    """
    input_path = Path(input_pdb_path)
    output_path = Path(output_pdb_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    try:
        from openmm import (
            CustomExternalForce,
            LangevinIntegrator,
            Platform,
            unit,
        )
        from openmm.app import (
            ForceField,
            HBonds,
            Modeller,
            NoCutoff,
            PDBFile,
            Simulation,
        )
    except Exception as exc:  # noqa: BLE001
        return JunctionRelaxResult(
            output_pdb_path=input_path,
            ran=False,
            fallback_reason=f"openmm import failed: {exc!r}",
            junction_rmsd_angstrom=None,
            restrained_rmsd_angstrom=None,
            final_energy_kj_per_mol=None,
            n_free_atoms=0,
            n_restrained_atoms=0,
        )

    # Pre-process: normalise chain-breaks that lack TER records so OpenMM
    # can distinguish real chain termini from mid-chain residues.  Multi-
    # domain crystals often ship as one chain with 4-100 A gaps and no TER;
    # OpenMM's Modeller.addHydrogens then complains that the residue on the
    # far side of the gap "is missing a terminal capping group".
    import tempfile as _tempfile
    _norm_path = Path(_tempfile.mkstemp(suffix="_normchains.pdb")[1])
    try:
        _split_chains_at_ca_breaks(input_path, _norm_path)
        pdb_path_for_openmm = _norm_path
    except Exception:  # noqa: BLE001 -- fall back to input if pre-proc fails
        pdb_path_for_openmm = input_path

    pdb = PDBFile(str(pdb_path_for_openmm))
    forcefield = ForceField(
        "amber14/protein.ff14SB.xml",
        "implicit/gbn2.xml",
    )

    # Heavy-atom-only PDBs need H-addition before ff14SB can build a system.
    # OpenMM's Modeller.addHydrogens is NOT pdbfixer -- it is the app-level
    # hydrogen-placement utility that uses ff14SB templates directly.  The
    # added H atoms are stripped from the output before write so downstream
    # FRUTON stages see the same heavy-atom-only interface they expect.
    modeller = Modeller(pdb.topology, pdb.positions)
    n_heavy_input = sum(1 for a in pdb.topology.atoms() if a.element.symbol != "H")
    try:
        modeller.addHydrogens(forcefield, pH=7.0)
    except Exception as exc:  # noqa: BLE001
        output_path.write_bytes(input_path.read_bytes())
        return JunctionRelaxResult(
            output_pdb_path=output_path,
            ran=False,
            fallback_reason=f"Modeller.addHydrogens failed: {exc!r}",
            junction_rmsd_angstrom=None,
            restrained_rmsd_angstrom=None,
            final_energy_kj_per_mol=None,
            n_free_atoms=0,
            n_restrained_atoms=0,
        )

    try:
        system = forcefield.createSystem(
            modeller.topology,
            nonbondedMethod=NoCutoff,
            constraints=HBonds,
        )
    except Exception as exc:  # noqa: BLE001
        # Usually a non-standard residue (PTR/SEP/TPO/ALY/CAS/...) with no
        # ff14SB template.  Return input untouched.
        output_path.write_bytes(input_path.read_bytes())
        return JunctionRelaxResult(
            output_pdb_path=output_path,
            ran=False,
            fallback_reason=f"ff14SB createSystem failed: {exc!r}",
            junction_rmsd_angstrom=None,
            restrained_rmsd_angstrom=None,
            final_energy_kj_per_mol=None,
            n_free_atoms=0,
            n_restrained_atoms=0,
        )

    free_ids = _free_residue_ids(gap_ranges, flank_residues)

    restraint = CustomExternalForce(
        "k*periodicdistance(x,y,z,x0,y0,z0)^2"
    )
    restraint.addGlobalParameter("k", restraint_k * unit.kilojoule_per_mole / unit.nanometer**2)
    restraint.addPerParticleParameter("x0")
    restraint.addPerParticleParameter("y0")
    restraint.addPerParticleParameter("z0")

    positions = modeller.positions
    n_restrained = 0
    n_free = 0
    restrained_indices: list[int] = []
    free_indices: list[int] = []
    for atom in modeller.topology.atoms():
        resid_num = atom.residue.id
        try:
            resid = int(resid_num)
        except (TypeError, ValueError):
            resid = -1
        if resid in free_ids:
            free_indices.append(atom.index)
            n_free += 1
            continue
        pos = positions[atom.index]
        x0 = pos.x if hasattr(pos, "x") else pos[0].value_in_unit(unit.nanometer)
        y0 = pos.y if hasattr(pos, "y") else pos[1].value_in_unit(unit.nanometer)
        z0 = pos.z if hasattr(pos, "z") else pos[2].value_in_unit(unit.nanometer)
        restraint.addParticle(atom.index, [x0, y0, z0])
        restrained_indices.append(atom.index)
        n_restrained += 1

    system.addForce(restraint)

    integrator = LangevinIntegrator(
        300 * unit.kelvin, 1.0 / unit.picosecond, 0.002 * unit.picoseconds
    )
    try:
        platform = Platform.getPlatformByName("Reference")
    except Exception:
        platform = None
    simulation = (
        Simulation(modeller.topology, system, integrator, platform)
        if platform is not None
        else Simulation(modeller.topology, system, integrator)
    )
    simulation.context.setPositions(positions)

    pre_state = simulation.context.getState(getPositions=True)
    pre_positions = pre_state.getPositions(asNumpy=False)
    pre_free = [
        (pre_positions[i].x, pre_positions[i].y, pre_positions[i].z)
        for i in free_indices
    ]
    pre_restrained = [
        (pre_positions[i].x, pre_positions[i].y, pre_positions[i].z)
        for i in restrained_indices
    ]

    simulation.minimizeEnergy(
        tolerance=tolerance * unit.kilojoule_per_mole / unit.nanometer,
        maxIterations=max_iterations,
    )

    post_state = simulation.context.getState(getPositions=True, getEnergy=True)
    post_positions = post_state.getPositions(asNumpy=False)
    post_free = [
        (post_positions[i].x, post_positions[i].y, post_positions[i].z)
        for i in free_indices
    ]
    post_restrained = [
        (post_positions[i].x, post_positions[i].y, post_positions[i].z)
        for i in restrained_indices
    ]

    junction_rmsd = _rmsd_angstrom(pre_free, post_free)
    restrained_rmsd = _rmsd_angstrom(pre_restrained, post_restrained)
    final_energy = post_state.getPotentialEnergy().value_in_unit(
        unit.kilojoule_per_mole
    )

    # Write heavy atoms only (strip added Hs) so downstream FRUTON stages
    # see the same file shape they saw on input.  Filter by element symbol.
    from openmm.app import PDBFile as _PDBFile
    heavy_indices = [
        a.index for a in modeller.topology.atoms() if a.element.symbol != "H"
    ]
    if len(heavy_indices) == n_heavy_input:
        heavy_positions = [post_positions[i] for i in heavy_indices]
        # Rebuild a heavy-only Modeller for writeFile
        heavy_modeller = Modeller(modeller.topology, post_positions)
        h_atoms = [a for a in heavy_modeller.topology.atoms() if a.element.symbol == "H"]
        heavy_modeller.delete(h_atoms)
        with output_path.open("w") as fh:
            _PDBFile.writeFile(heavy_modeller.topology, heavy_modeller.positions, fh)
    else:
        # Unexpected mismatch: write the full protonated model as a safe fallback.
        with output_path.open("w") as fh:
            _PDBFile.writeFile(modeller.topology, post_positions, fh)
    _copy_original_numbering(input_path, output_path, output_path)

    return JunctionRelaxResult(
        output_pdb_path=output_path,
        ran=True,
        fallback_reason=None,
        junction_rmsd_angstrom=junction_rmsd,
        restrained_rmsd_angstrom=restrained_rmsd,
        final_energy_kj_per_mol=final_energy,
        n_free_atoms=n_free,
        n_restrained_atoms=n_restrained,
    )
