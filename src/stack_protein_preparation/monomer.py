"""Analyze and write representative monomer units from PDB structures.

This module identifies equivalent coordinate chains in a PDB file and writes
representative protein units for downstream preparation workflows. It is
intentionally conservative: it does not infer biological assemblies, rebuild
missing residues, renumber residues, assign chemistry, or decide which
ligand/cofactor is biologically essential. The module is meant as a preprocessing
helper for cases where a PDB file may contain repeated copies of the same
protein chain and downstream workflow steps should receive a representative
protein unit rather than the whole multichain structure.

Chain equivalence is based on coordinate-derived amino-acid sequences. Sequences
are extracted from the first model only, because most docking and parametrization
workflows expect a single coordinate model. Pairwise global alignments are used
to calculate sequence identity and bidirectional coverage, and equivalent chains
are grouped conservatively using complete-linkage logic. The write-out step uses
BioPython's PDBIO selection mechanism so coordinate output stays delegated to
BioPython instead of relying on local PDB line parsing.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any

from Bio import Align
from Bio.Data.PDBData import protein_letters_3to1_extended
from Bio.PDB import PDBIO, PDBParser
from Bio.PDB.PDBIO import Select
from Bio.PDB.Polypeptide import is_aa

DEFAULT_IDENTITY_THRESHOLD = 0.95
DEFAULT_COVERAGE_THRESHOLD = 0.90


@dataclass(frozen=True)
class ChainSequence:
    """Store the coordinate-derived sequence for one PDB chain.

    ``chain_id`` is the identifier used by BioPython for the chain in the first
    model. ``sequence`` is built from residues that BioPython recognizes as
    amino acids with ``standard=False``, which allows common modified residues to
    contribute to the chain sequence. Unknown recognized residue names are mapped
    to ``X`` instead of failing, because some structures contain unusual residue
    names that should be reported without crashing the full workflow. The
    counters are used for deterministic representative-chain selection.
    """

    chain_id: str
    sequence: str
    residue_count: int
    unknown_residue_count: int

    @property
    def unknown_fraction(self) -> float:
        """Return the fraction of residues represented as ``X``.

        This is mainly useful for diagnostics and ranking representative chains.
        A shorter unknown fraction usually means that BioPython could map more
        residue names to meaningful one-letter amino-acid codes. Empty sequences
        are not emitted by the extractor, but the property remains defensive.
        """

        if self.residue_count == 0:
            return 0.0
        return self.unknown_residue_count / self.residue_count


@dataclass(frozen=True)
class ChainPairComparison:
    """Store the result of a pairwise chain-sequence comparison.

    ``sequence_identity`` is calculated over aligned residue-residue pairs and
    does not count ``X`` versus ``X`` as an identity match. ``coverage_a`` and
    ``coverage_b`` measure how much of each original chain sequence is paired
    against a residue in the other chain instead of aligned to a gap. This avoids
    grouping a heavily truncated chain with a complete chain merely because the
    overlapping part is identical. The aligned and matching pair counters are
    included so callers can inspect borderline cases.
    """

    chain_id_a: str
    chain_id_b: str
    sequence_identity: float
    coverage_a: float
    coverage_b: float
    aligned_residue_pairs: int
    matching_residue_pairs: int


@dataclass(frozen=True)
class MonomerGroup:
    """Represent one detected chain-equivalence group.

    A group contains chains that passed the configured identity and coverage
    thresholds against every other chain in the same group. The representative is
    selected deterministically by preferring the longest chain, then the chain
    with fewer unknown residues, then the lexicographically first chain ID. This
    makes repeated runs stable while usually selecting the most complete chain as
    output representative. The representative sequence is stored for inspection
    and reporting.
    """

    representative_chain_id: str
    chain_ids: tuple[str, ...]
    representative_sequence: str


@dataclass(frozen=True)
class MonomerAnalysisResult:
    """Describe the complete monomer analysis for one PDB file.

    ``chains`` contains the coordinate-derived sequences extracted from the first
    model. ``comparisons`` contains all pairwise comparisons used during
    grouping. ``groups`` is the conservative final clustering result and is used
    by the representative write helpers. The thresholds are stored so JSON
    summaries or test assertions can report exactly which criteria produced the
    grouping.
    """

    input_pdb: Path
    model_id: Any
    chains: tuple[ChainSequence, ...]
    comparisons: tuple[ChainPairComparison, ...]
    groups: tuple[MonomerGroup, ...]
    identity_threshold: float
    coverage_threshold: float

    @property
    def is_single_monomer_type(self) -> bool:
        """Return whether all detected chains belong to one chain type."""

        return len(self.groups) == 1


@dataclass(frozen=True)
class MonomerWriteResult:
    """Describe one representative protein-unit PDB written to disk.

    ``representative_chain_ids`` are the actual chains copied to the output file.
    ``represented_chain_groups`` records which equivalent source chains each
    copied chain represents. For a homomeric structure this usually contains one
    representative chain. For a heteromeric repeated unit, such as A/B plus C/D,
    this contains one representative chain per non-equivalent group, for example
    A and C.

    The record counters are calculated from the written file so they reflect
    what BioPython actually emitted. Non-protein HETATM records are retained only
    when they belong to one of the selected representative chains and
    ``keep_non_protein_hetero`` is enabled. Compatibility properties are provided
    for older callers that expected a single representative chain.
    """

    output_pdb: Path
    source_pdb: Path
    model_id: Any
    representative_chain_ids: tuple[str, ...]
    represented_chain_groups: tuple[tuple[str, ...], ...]
    atom_records_written: int
    hetatm_records_written: int
    ter_records_written: int

    @property
    def coordinate_records_written(self) -> int:
        """Return the total number of written coordinate records."""

        return self.atom_records_written + self.hetatm_records_written

    @property
    def representative_chain_id(self) -> str:
        """Return a legacy string representation of representative chains."""

        return "/".join(self.representative_chain_ids)

    @property
    def represented_chain_ids(self) -> tuple[str, ...]:
        """Return all represented source chain IDs as a flattened tuple."""

        return tuple(
            chain_id
            for group in self.represented_chain_groups
            for chain_id in group
        )


def analyze_monomer_units(
    input_pdb: str | Path,
    *,
    identity_threshold: float = DEFAULT_IDENTITY_THRESHOLD,
    coverage_threshold: float = DEFAULT_COVERAGE_THRESHOLD,
) -> MonomerAnalysisResult:
    """Analyze a PDB file and group equivalent protein chains.

    The function parses the input PDB with BioPython, extracts one
    coordinate-derived amino-acid sequence per chain from the first model, and
    compares all chains pairwise. Chains are considered equivalent only when they
    pass both the sequence-identity threshold and the bidirectional coverage
    threshold. The grouping is deliberately conservative and uses complete
    linkage so a chain can join a group only if it is equivalent to every current
    group member.

    This function performs no coordinate write-out. It is the diagnostic entry
    point to check whether a structure behaves like a homomer, heteromer, or a
    mixed structure with truncated/non-equivalent chains. A ``ValueError`` is
    raised when no protein-like coordinate chains are found.
    """

    input_path = Path(input_pdb)
    _validate_input_pdb(input_path)
    _validate_threshold("identity_threshold", identity_threshold)
    _validate_threshold("coverage_threshold", coverage_threshold)

    parsed = _parse_first_model(input_path)
    chains = extract_chain_sequences(input_path)

    if not chains:
        raise ValueError(f"No protein-like coordinate chains were found in {input_path}.")

    comparisons = compare_chain_sequences(chains)
    groups = _group_equivalent_chains(
        chains=chains,
        comparisons=comparisons,
        identity_threshold=identity_threshold,
        coverage_threshold=coverage_threshold,
    )

    return MonomerAnalysisResult(
        input_pdb=input_path,
        model_id=parsed.model.id,
        chains=chains,
        comparisons=comparisons,
        groups=groups,
        identity_threshold=identity_threshold,
        coverage_threshold=coverage_threshold,
    )


def extract_chain_sequences(input_pdb: str | Path) -> tuple[ChainSequence, ...]:
    """Extract coordinate-derived amino-acid sequences from the first model.

    The extractor iterates through BioPython chain and residue objects rather
    than reading PDB columns manually. Residues are included when BioPython
    recognizes them as amino acids with ``standard=False``, which includes common
    modified amino acids such as selenomethionine. The sequence reflects only
    residues present in the coordinate model, not SEQRES records. Blank chain
    identifiers are skipped because output filenames and downstream decisions
    should remain explicit and reproducible.
    """

    input_path = Path(input_pdb)
    _validate_input_pdb(input_path)

    parsed = _parse_first_model(input_path)
    chain_sequences: list[ChainSequence] = []

    for chain in parsed.model:
        chain_id = str(chain.id)
        if not chain_id or chain_id.isspace():
            continue

        residues: list[str] = []
        unknown_residue_count = 0

        for residue in chain:
            if not is_aa(residue, standard=False):
                continue

            one_letter = _residue_to_one_letter(residue)
            residues.append(one_letter)

            if one_letter == "X":
                unknown_residue_count += 1

        if not residues:
            continue

        sequence = "".join(residues)
        chain_sequences.append(
            ChainSequence(
                chain_id=chain_id,
                sequence=sequence,
                residue_count=len(sequence),
                unknown_residue_count=unknown_residue_count,
            )
        )

    return tuple(sorted(chain_sequences, key=lambda item: item.chain_id))


def compare_chain_sequences(
    chains: tuple[ChainSequence, ...],
) -> tuple[ChainPairComparison, ...]:
    """Compare all chain sequences pairwise by global alignment.

    The comparison uses BioPython's ``PairwiseAligner`` instead of a local
    alignment implementation. The scoring is intentionally simple because the
    goal is to identify equivalent chain copies, not to publish a sequence
    alignment. The returned comparisons are sorted by the input chain order,
    which should already be deterministic because chain extraction sorts by
    chain ID.
    """

    aligner = _make_pairwise_aligner()
    comparisons: list[ChainPairComparison] = []

    for index_a, chain_a in enumerate(chains):
        for chain_b in chains[index_a + 1 :]:
            comparisons.append(
                _compare_two_chain_sequences(
                    chain_a=chain_a,
                    chain_b=chain_b,
                    aligner=aligner,
                )
            )

    return tuple(comparisons)


def write_representative_monomer_units(
    input_pdb: str | Path,
    output_dir: str | Path,
    *,
    identity_threshold: float = DEFAULT_IDENTITY_THRESHOLD,
    coverage_threshold: float = DEFAULT_COVERAGE_THRESHOLD,
    keep_modified_residues: bool = True,
    keep_non_protein_hetero: bool = True,
    keep_waters: bool = False,
    overwrite: bool = True,
) -> tuple[MonomerWriteResult, ...]:
    """Write one representative PDB file per detected chain group.

    The function first runs ``analyze_monomer_units`` and then writes only the
    representative chain for each detected group. For a homomeric structure this
    normally produces one output file. For a heteromeric structure this produces
    one representative file per non-equivalent chain group, which is useful for
    diagnostics and for workflows that intentionally want separate chain-type
    files.

    Standard amino-acid residues are always retained for the representative
    chain. Modified amino-acid residues can be retained separately from generic
    non-protein hetero residues so users do not need to choose between losing
    common modified residues and carrying every ligand or ion into the monomer
    output. Water retention is separate and remains disabled by default.
    """

    input_path = Path(input_pdb)
    output_path = Path(output_dir)

    analysis = analyze_monomer_units(
        input_path,
        identity_threshold=identity_threshold,
        coverage_threshold=coverage_threshold,
    )

    output_path.mkdir(parents=True, exist_ok=True)

    results: list[MonomerWriteResult] = []
    for group in analysis.groups:
        chain_label = _safe_chain_label(group.representative_chain_id)
        output_pdb = output_path / f"{input_path.stem}_monomer_chain_{chain_label}.pdb"

        results.append(
            write_chain_unit(
                input_pdb=input_path,
                output_pdb=output_pdb,
                chain_id=group.representative_chain_id,
                represented_chain_ids=group.chain_ids,
                keep_modified_residues=keep_modified_residues,
                keep_non_protein_hetero=keep_non_protein_hetero,
                keep_waters=keep_waters,
                overwrite=overwrite,
            )
        )

    return tuple(results)


def write_representative_monomer_unit(
    input_pdb: str | Path,
    output_pdb: str | Path,
    *,
    identity_threshold: float = DEFAULT_IDENTITY_THRESHOLD,
    coverage_threshold: float = DEFAULT_COVERAGE_THRESHOLD,
    keep_modified_residues: bool = True,
    keep_non_protein_hetero: bool = True,
    keep_waters: bool = False,
    overwrite: bool = True,
) -> MonomerWriteResult:
    """Write one representative protein unit to a single PDB file.

    The output contains one representative chain from each detected chain group.
    For a homomeric structure such as A/B, this normally writes only A. For a
    heteromeric repeated unit such as A/B plus C/D, this writes A and C into the
    same output file. Non-protein HETATM records are retained only when they are
    assigned to one of the selected representative chains, so ligands are kept
    per representative unit instead of copied from discarded duplicate chains.

    The function does not infer the biological assembly from REMARK records or
    external assembly metadata. It only removes duplicated equivalent coordinate
    chains according to sequence identity and coverage. Water retention remains
    explicit because crystallographic waters are often not intended to define the
    representative protein unit.
    """

    input_path = Path(input_pdb)
    output_path = Path(output_pdb)

    analysis = analyze_monomer_units(
        input_path,
        identity_threshold=identity_threshold,
        coverage_threshold=coverage_threshold,
    )

    representative_chain_ids = tuple(
        group.representative_chain_id for group in analysis.groups
    )
    represented_chain_groups = tuple(group.chain_ids for group in analysis.groups)

    return write_chain_unit_set(
        input_pdb=input_path,
        output_pdb=output_path,
        chain_ids=representative_chain_ids,
        represented_chain_groups=represented_chain_groups,
        keep_modified_residues=keep_modified_residues,
        keep_non_protein_hetero=keep_non_protein_hetero,
        keep_waters=keep_waters,
        overwrite=overwrite,
    )


def write_single_representative_monomer_unit(
    input_pdb: str | Path,
    output_pdb: str | Path,
    *,
    identity_threshold: float = DEFAULT_IDENTITY_THRESHOLD,
    coverage_threshold: float = DEFAULT_COVERAGE_THRESHOLD,
    keep_modified_residues: bool = True,
    keep_non_protein_hetero: bool = True,
    keep_waters: bool = False,
    overwrite: bool = True,
) -> MonomerWriteResult:
    """Write the single representative unit expected by downstream workflows.

    The historical name is kept for compatibility with the FRUTON runner and
    earlier tests. The function no longer means "write exactly one chain". It
    means "write exactly one representative unit file", and that unit can contain
    multiple chains when the structure has multiple non-equivalent chain groups.

    For example, duplicated heteromeric units grouped as A/A-B and C/C-D are
    written as one output containing chains A and C. Ligands assigned to A and C
    are retained by default, while ligands assigned to discarded duplicate chains
    such as B and D are not written. Use ``write_representative_monomer_units``
    when separate output files per chain group are wanted instead.
    """

    return write_representative_monomer_unit(
        input_pdb=input_pdb,
        output_pdb=output_pdb,
        identity_threshold=identity_threshold,
        coverage_threshold=coverage_threshold,
        keep_modified_residues=keep_modified_residues,
        keep_non_protein_hetero=keep_non_protein_hetero,
        keep_waters=keep_waters,
        overwrite=overwrite,
    )


def write_chain_unit(
    input_pdb: str | Path,
    output_pdb: str | Path,
    chain_id: str,
    *,
    represented_chain_ids: tuple[str, ...] | None = None,
    keep_modified_residues: bool = True,
    keep_non_protein_hetero: bool = False,
    keep_waters: bool = False,
    overwrite: bool = True,
) -> MonomerWriteResult:
    """Write one selected chain from the first model to a PDB file.

    This is the low-level write helper used after representative-chain selection
    has already happened. It uses BioPython ``PDBIO`` with a custom ``Select``
    object instead of manually copying PDB text lines. The selected model and
    chain are retained, while residue-level filtering controls modified amino
    acids, generic hetero residues, and waters.

    This function does not decide whether the requested chain is biologically
    representative. That decision belongs to ``analyze_monomer_units`` and the
    representative write wrappers. The low-level default keeps non-protein HETATM
    disabled so direct callers remain explicit about ligand/cofactor retention.
    """

    return write_chain_unit_set(
        input_pdb=input_pdb,
        output_pdb=output_pdb,
        chain_ids=(chain_id,),
        represented_chain_groups=(represented_chain_ids or (chain_id,),),
        keep_modified_residues=keep_modified_residues,
        keep_non_protein_hetero=keep_non_protein_hetero,
        keep_waters=keep_waters,
        overwrite=overwrite,
    )


def write_chain_unit_set(
    input_pdb: str | Path,
    output_pdb: str | Path,
    chain_ids: tuple[str, ...],
    *,
    represented_chain_groups: tuple[tuple[str, ...], ...] | None = None,
    keep_modified_residues: bool = True,
    keep_non_protein_hetero: bool = False,
    keep_waters: bool = False,
    overwrite: bool = True,
) -> MonomerWriteResult:
    """Write a selected set of chains from the first model to one PDB file.

    This helper is the shared output path for both single-chain and heteromeric
    representative-unit writing. It keeps all requested chains in one file while
    applying the same residue-level filtering rules to each chain. Non-protein
    HETATM records are kept only when their chain belongs to ``chain_ids`` and
    ``keep_non_protein_hetero`` is enabled.

    The function does not compare chains or choose representatives. Callers must
    pass already selected chain IDs, usually from ``analyze_monomer_units``. A
    ``ValueError`` is raised when the output contains no coordinate records,
    because that would indicate either an invalid selection or overly aggressive
    residue filtering.
    """

    input_path = Path(input_pdb)
    output_path = Path(output_pdb)

    _validate_input_pdb(input_path)
    _validate_output_pdb(output_path, overwrite=overwrite)
    _validate_chain_ids(chain_ids)

    parsed = _parse_first_model(input_path)
    for chain_id in chain_ids:
        _validate_chain_exists(parsed.model, chain_id)

    output_path.parent.mkdir(parents=True, exist_ok=True)

    io = PDBIO()
    io.set_structure(parsed.structure)
    io.save(
        str(output_path),
        select=_RepresentativeChainsSelect(
            model_id=parsed.model.id,
            chain_ids=chain_ids,
            keep_modified_residues=keep_modified_residues,
            keep_non_protein_hetero=keep_non_protein_hetero,
            keep_waters=keep_waters,
        ),
        preserve_atom_numbering=True,
    )

    atom_count, hetatm_count, ter_count = _count_written_pdb_records(output_path)

    if atom_count + hetatm_count == 0:
        output_path.unlink(missing_ok=True)
        chain_text = "/".join(chain_ids)
        raise ValueError(
            f"No coordinate records were written for chains {chain_text!r} "
            f"from {input_path}."
        )

    if represented_chain_groups is None:
        represented_chain_groups = tuple((chain_id,) for chain_id in chain_ids)

    return MonomerWriteResult(
        output_pdb=output_path,
        source_pdb=input_path,
        model_id=parsed.model.id,
        representative_chain_ids=chain_ids,
        represented_chain_groups=represented_chain_groups,
        atom_records_written=atom_count,
        hetatm_records_written=hetatm_count,
        ter_records_written=ter_count,
    )


@dataclass(frozen=True)
class _ParsedFirstModel:
    """Internal container for a parsed BioPython structure and its first model."""

    structure: Any
    model: Any


class _RepresentativeChainsSelect(Select):
    """BioPython PDBIO selector for one representative model and chain set."""

    def __init__(
        self,
        *,
        model_id: Any,
        chain_ids: tuple[str, ...],
        keep_modified_residues: bool,
        keep_non_protein_hetero: bool,
        keep_waters: bool,
    ) -> None:
        self._model_id = model_id
        self._chain_ids = set(chain_ids)
        self._keep_modified_residues = keep_modified_residues
        self._keep_non_protein_hetero = keep_non_protein_hetero
        self._keep_waters = keep_waters

    def accept_model(self, model: Any) -> bool:
        return model.id == self._model_id

    def accept_chain(self, chain: Any) -> bool:
        return str(chain.id) in self._chain_ids

    def accept_residue(self, residue: Any) -> bool:
        hetero_flag = residue.id[0]

        if hetero_flag == " ":
            if is_aa(residue, standard=True):
                return True
            if is_aa(residue, standard=False):
                return self._keep_modified_residues
            return False

        if hetero_flag == "W":
            return self._keep_waters

        if is_aa(residue, standard=False):
            return self._keep_modified_residues

        return self._keep_non_protein_hetero


def _parse_first_model(input_pdb: Path) -> _ParsedFirstModel:
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure(input_pdb.stem, str(input_pdb))

    models = list(structure.get_models())
    if not models:
        raise ValueError(f"No coordinate models were found in {input_pdb}.")

    return _ParsedFirstModel(structure=structure, model=models[0])


def _residue_to_one_letter(residue: Any) -> str:
    residue_name = residue.get_resname().strip().upper()
    return protein_letters_3to1_extended.get(residue_name, "X")


def _make_pairwise_aligner() -> Align.PairwiseAligner:
    aligner = Align.PairwiseAligner()
    aligner.mode = "global"
    aligner.match_score = 2.0
    aligner.mismatch_score = -1.0
    aligner.open_gap_score = -2.0
    aligner.extend_gap_score = -0.5
    return aligner


def _compare_two_chain_sequences(
    *,
    chain_a: ChainSequence,
    chain_b: ChainSequence,
    aligner: Align.PairwiseAligner,
) -> ChainPairComparison:
    alignment = aligner.align(chain_a.sequence, chain_b.sequence)[0]

    aligned_residue_pairs = 0
    matching_residue_pairs = 0
    paired_residues_a = 0
    paired_residues_b = 0

    aligned_blocks_a, aligned_blocks_b = alignment.aligned

    for block_a, block_b in zip(aligned_blocks_a, aligned_blocks_b):
        start_a, end_a = int(block_a[0]), int(block_a[1])
        start_b, end_b = int(block_b[0]), int(block_b[1])

        block_length = min(end_a - start_a, end_b - start_b)

        for offset in range(block_length):
            residue_a = chain_a.sequence[start_a + offset]
            residue_b = chain_b.sequence[start_b + offset]

            aligned_residue_pairs += 1
            paired_residues_a += 1
            paired_residues_b += 1

            if residue_a == residue_b and residue_a != "X":
                matching_residue_pairs += 1

    sequence_identity = (
        matching_residue_pairs / aligned_residue_pairs
        if aligned_residue_pairs
        else 0.0
    )
    coverage_a = paired_residues_a / chain_a.residue_count
    coverage_b = paired_residues_b / chain_b.residue_count

    return ChainPairComparison(
        chain_id_a=chain_a.chain_id,
        chain_id_b=chain_b.chain_id,
        sequence_identity=sequence_identity,
        coverage_a=coverage_a,
        coverage_b=coverage_b,
        aligned_residue_pairs=aligned_residue_pairs,
        matching_residue_pairs=matching_residue_pairs,
    )


def _group_equivalent_chains(
    *,
    chains: tuple[ChainSequence, ...],
    comparisons: tuple[ChainPairComparison, ...],
    identity_threshold: float,
    coverage_threshold: float,
) -> tuple[MonomerGroup, ...]:
    comparison_by_pair = {
        frozenset((comparison.chain_id_a, comparison.chain_id_b)): comparison
        for comparison in comparisons
    }

    raw_groups: list[list[ChainSequence]] = []

    for chain in sorted(chains, key=lambda item: item.chain_id):
        matching_group: list[ChainSequence] | None = None

        for group in raw_groups:
            if all(
                _chains_are_equivalent(
                    chain.chain_id,
                    member.chain_id,
                    comparison_by_pair=comparison_by_pair,
                    identity_threshold=identity_threshold,
                    coverage_threshold=coverage_threshold,
                )
                for member in group
            ):
                matching_group = group
                break

        if matching_group is None:
            raw_groups.append([chain])
        else:
            matching_group.append(chain)

    groups: list[MonomerGroup] = []
    for group in raw_groups:
        representative = sorted(
            group,
            key=lambda item: (
                -item.residue_count,
                item.unknown_residue_count,
                item.chain_id,
            ),
        )[0]

        groups.append(
            MonomerGroup(
                representative_chain_id=representative.chain_id,
                chain_ids=tuple(sorted(item.chain_id for item in group)),
                representative_sequence=representative.sequence,
            )
        )

    return tuple(sorted(groups, key=lambda item: item.representative_chain_id))


def _chains_are_equivalent(
    chain_id_a: str,
    chain_id_b: str,
    *,
    comparison_by_pair: dict[frozenset[str], ChainPairComparison],
    identity_threshold: float,
    coverage_threshold: float,
) -> bool:
    if chain_id_a == chain_id_b:
        return True

    comparison = comparison_by_pair[frozenset((chain_id_a, chain_id_b))]

    return (
        comparison.sequence_identity >= identity_threshold
        and comparison.coverage_a >= coverage_threshold
        and comparison.coverage_b >= coverage_threshold
    )


def _validate_input_pdb(input_pdb: Path) -> None:
    if not input_pdb.exists():
        raise FileNotFoundError(f"Input PDB does not exist: {input_pdb}")
    if not input_pdb.is_file():
        raise ValueError(f"Input PDB is not a file: {input_pdb}")
    if input_pdb.suffix.lower() != ".pdb":
        raise ValueError(f"Expected a .pdb input file, got: {input_pdb}")


def _validate_output_pdb(output_pdb: Path, *, overwrite: bool) -> None:
    if output_pdb.suffix.lower() != ".pdb":
        raise ValueError(f"Expected a .pdb output file, got: {output_pdb}")
    if output_pdb.exists() and not overwrite:
        raise FileExistsError(f"Output PDB already exists: {output_pdb}")


def _validate_chain_id(chain_id: str) -> None:
    if len(chain_id) != 1:
        raise ValueError(f"chain_id must be exactly one character, got {chain_id!r}.")
    if chain_id.isspace():
        raise ValueError("chain_id must not be blank.")


def _validate_chain_ids(chain_ids: tuple[str, ...]) -> None:
    if not chain_ids:
        raise ValueError("At least one chain ID must be selected.")

    seen_chain_ids: set[str] = set()
    for chain_id in chain_ids:
        _validate_chain_id(chain_id)
        if chain_id in seen_chain_ids:
            raise ValueError(f"Duplicate selected chain ID: {chain_id!r}.")
        seen_chain_ids.add(chain_id)


def _validate_chain_exists(model: Any, chain_id: str) -> None:
    if chain_id not in {str(chain.id) for chain in model}:
        available = ", ".join(str(chain.id) for chain in model)
        raise ValueError(
            f"Chain {chain_id!r} was not found in first model. "
            f"Available chains: {available or '<none>'}."
        )


def _validate_threshold(name: str, value: float) -> None:
    if not 0.0 <= value <= 1.0:
        raise ValueError(f"{name} must be between 0.0 and 1.0, got {value!r}.")


def _safe_chain_label(chain_id: str) -> str:
    if chain_id.isalnum():
        return chain_id
    return f"ord_{ord(chain_id):02x}"


def _count_written_pdb_records(output_pdb: Path) -> tuple[int, int, int]:
    atom_count = 0
    hetatm_count = 0
    ter_count = 0

    for line in output_pdb.read_text(encoding="utf-8").splitlines():
        if line.startswith("ATOM  "):
            atom_count += 1
        elif line.startswith("HETATM"):
            hetatm_count += 1
        elif line.startswith("TER"):
            ter_count += 1

    return atom_count, hetatm_count, ter_count