from __future__ import annotations

import json
import shutil
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any

from Bio.Data.PDBData import protein_letters_3to1, protein_letters_3to1_extended
from Bio.PDB import PDBIO, Chain, Model, PDBParser, Residue, Structure
from Bio.PDB.PDBIO import Select
from Bio.PDB.Polypeptide import is_aa

from stack_protein_preparation.pipeline_logging import (
    append_key_value_block,
    append_log_text,
    build_module_log_paths,
    print_double_header_box,
    print_light_table_box,
    print_mixed_box,
)

MODULE_NAME = "nonstandard_capping"

FORCE_FIELD_RESIDUE_ALIASES = {
    "ASH": "ASP",
    "GLH": "GLU",
    "HID": "HIS",
    "HIE": "HIS",
    "HIP": "HIS",
    "CYM": "CYS",
    "CYX": "CYS",
    "LYN": "LYS",
}

# Non-standard backbone residues that need RESP parametrization.
# BioPython's is_aa(standard=False) misses many of these (e.g. CSD, OCS, MLY),
# so we check this set explicitly in _is_protein_like_residue.
# Protonation-state variants (HID, HIE, HIP, CYM, CYX, ASH, GLH, LYN) are
# intentionally excluded — they are handled by standard force fields.
_BACKBONE_NONSTD_RESNAMES: frozenset[str] = frozenset({
    "MSE", "SEC", "PYL",          # selenium / pyrrolysine
    "SEP", "TPO", "PTR",          # phosphoserine/threonine/tyrosine
    "S1P", "T1P", "Y1P", "NEP",   # monoprotonated phospho variants
    "CSO", "CSD", "OCS", "CSS", "CME", "SCY", "SMC", "P1L", "CAF",  # Cys oxidations
    "MHO", "SME", "OMT",          # Met oxidations
    "TYS", "NIY", "TYI", "IYR", "TPQ", "TYQ",  # Tyr modifications
    "MLZ", "MLY", "M3L", "ALY", "KCR", "SLL", "KCX", "LYZ", "LLP",  # Lys modifications
    "AGM", "DA2", "MMO", "MMA",   # Arg methylations
    "CIR",                         # citrulline
    "CGU",                         # gamma-carboxyglutamate
    "HYP", "4FB",                  # Pro hydroxylation
    "4HT", "TRO", "HTR",           # Trp modifications
    "HIC",                         # 4-methyl-His
    "PCA",                         # pyroglutamate
    "SAR", "MAA", "MEA", "MEN",    # backbone N-methylation
    "SAC",                         # N-terminal acetylation
    "CRO", "DYG",                  # GFP chromophores
    "H2D", "H1D", "H2E", "H1E",   # phosphohistidine (Bryce group)
})


@dataclass(slots=True, frozen=True)
class Fragment:
    chain_id: str
    start_index: int
    end_index: int
    start_resseq: int
    end_resseq: int
    start_resname: str
    end_resname: str


@dataclass(slots=True, frozen=True)
class InternalGapBoundary:
    chain_id: str
    left_residue_index: int
    right_residue_index: int
    left_resseq: int
    right_resseq: int
    left_resname: str
    right_resname: str


@dataclass(slots=True, frozen=True)
class ChainCappingSummary:
    chain_id: str
    n_fragments: int
    n_internal_boundaries: int
    n_boundaries_capped: int


@dataclass(slots=True, frozen=True)
class NonstandardResidueRecord:
    pdb_id: str
    chain_id: str
    residue_name: str
    residue_number: int
    insertion_code: str
    residue_label: str
    source_pdb_path: Path
    extracted_pdb_path: Path
    capped_pdb_path: Path
    manifest_path: Path
    backend: str
    success: bool
    error_message: str = ""


@dataclass(slots=True, frozen=True)
class NonstandardCappingSummary:
    pdb_id: str
    input_pdb_path: Path
    output_dir: Path
    n_nonstandard_residues: int
    n_capped_models_written: int
    records: tuple[NonstandardResidueRecord, ...]
    manifest_path: Path


@dataclass(slots=True, frozen=True)
class _ParsedFirstModel:
    structure: Any
    model: Any


class _SingleResidueSelect(Select):
    """BioPython selector for one residue in one chain from the first model.

    This selector is used to write the isolated non-standard residue before the
    PTM-Psi capping backend is called. It avoids manual PDB text slicing and
    keeps residue selection tied to BioPython's parsed chain and residue
    objects. The selected residue may be ATOM or HETATM as long as it was
    classified as protein-like by the caller.
    """

    def __init__(
        self,
        *,
        model_id: Any,
        chain_id: str,
        residue_id: tuple[str, int, str],
    ) -> None:
        self._model_id = model_id
        self._chain_id = _normalize_chain_id(chain_id)
        self._residue_id = residue_id

    def accept_model(self, model: Any) -> bool:
        return model.id == self._model_id

    def accept_chain(self, chain: Any) -> bool:
        return _normalize_chain_id(str(chain.id)) == self._chain_id

    def accept_residue(self, residue: Any) -> bool:
        return residue.id == self._residue_id


class _FragmentSelect(Select):
    """BioPython selector for one connected residue fragment.

    Internal capping is implemented by splitting a chain into connected
    fragments, writing each fragment, and asking PTM-Psi to add ACE/NME to the
    fragment termini that correspond to internal breaks. The selector writes the
    requested fragment without direct PDB column parsing. Only protein-like
    residues are emitted.
    """

    def __init__(
        self,
        *,
        model_id: Any,
        chain_id: str,
        residue_ids: set[tuple[str, int, str]],
    ) -> None:
        self._model_id = model_id
        self._chain_id = _normalize_chain_id(chain_id)
        self._residue_ids = residue_ids

    def accept_model(self, model: Any) -> bool:
        return model.id == self._model_id

    def accept_chain(self, chain: Any) -> bool:
        return _normalize_chain_id(str(chain.id)) == self._chain_id

    def accept_residue(self, residue: Any) -> bool:
        return residue.id in self._residue_ids and _is_protein_like_residue(residue)


def sanitize_variant_label(variant_label: str) -> str:
    cleaned = "".join(
        character if character.isalnum() or character in {"_", "-"} else "_"
        for character in variant_label.strip()
    )
    cleaned = cleaned.strip("_")
    if not cleaned:
        raise ValueError(f"Invalid empty variant label derived from {variant_label!r}")
    return cleaned


def build_standard_internal_capped_output_path(
    pdb_id: str,
    protein_dir: str | Path,
) -> Path:
    protein_dir = Path(protein_dir)
    return protein_dir / "components" / f"{pdb_id}_protein_internal_capped.pdb"


def build_variant_internal_capped_output_path(
    pdb_id: str,
    protein_dir: str | Path,
    variant_label: str,
) -> Path:
    protein_dir = Path(protein_dir)
    safe_variant_label = sanitize_variant_label(variant_label)
    return (
        protein_dir
        / "components"
        / f"{pdb_id}_protein_internal_capped_{safe_variant_label}.pdb"
    )


def build_nonstandard_capping_output_dir(
    pdb_id: str,
    protein_dir: str | Path,
) -> Path:
    protein_dir = Path(protein_dir)
    return protein_dir / "nonstandard_params"


def _build_log_path(
    *,
    input_pdb_path: Path,
    variant_label: str | None,
) -> Path:
    protein_dir = input_pdb_path.resolve().parents[1]
    protein_data_dir = protein_dir.parent
    pdb_id = protein_dir.name.upper()

    log_paths = build_module_log_paths(
        protein_data_dir=protein_data_dir,
        pdb_id=pdb_id,
        module_name=MODULE_NAME,
        variant_label=variant_label,
    )
    return log_paths.module_log_path


def _normalize_chain_id(chain_id: str | None) -> str:
    return (chain_id or "").strip() or "_"


def _normalize_residue_name_for_sequence(residue_name: str) -> str:
    normalized_name = residue_name.strip().upper()
    if normalized_name in FORCE_FIELD_RESIDUE_ALIASES:
        return FORCE_FIELD_RESIDUE_ALIASES[normalized_name]

    if len(normalized_name) == 4 and normalized_name[0] in {"N", "C"}:
        candidate = normalized_name[1:]
        if candidate in protein_letters_3to1_extended:
            return candidate

    return normalized_name


def _is_standard_or_forcefield_variant(residue: Residue.Residue) -> bool:
    if is_aa(residue, standard=True):
        return True

    residue_name = residue.get_resname().strip().upper()
    normalized_name = _normalize_residue_name_for_sequence(residue_name)

    return normalized_name in protein_letters_3to1


def _is_protein_like_residue(residue: Residue.Residue) -> bool:
    residue_name = residue.get_resname().strip().upper()
    normalized_name = _normalize_residue_name_for_sequence(residue_name)

    return (
        is_aa(residue, standard=False)
        or normalized_name in protein_letters_3to1_extended
        or residue_name in _BACKBONE_NONSTD_RESNAMES
    )


def _is_nonstandard_protein_residue(residue: Residue.Residue) -> bool:
    if not _is_protein_like_residue(residue):
        return False

    return not _is_standard_or_forcefield_variant(residue)


def _residue_has_backbone_for_capping(residue: Residue.Residue) -> bool:
    atom_names = {str(atom.id).strip().upper() for atom in residue.get_atoms()}
    return {"N", "CA", "C"}.issubset(atom_names)


def _parse_first_model(input_pdb_path: Path) -> _ParsedFirstModel:
    if not input_pdb_path.exists():
        raise FileNotFoundError(f"Input PDB not found: {input_pdb_path}")

    parser = PDBParser(QUIET=True)
    structure = parser.get_structure(input_pdb_path.stem, str(input_pdb_path))
    models = list(structure.get_models())

    if not models:
        raise ValueError(f"No coordinate model found in input PDB: {input_pdb_path}")

    return _ParsedFirstModel(structure=structure, model=models[0])


def _residue_label(
    *,
    chain_id: str,
    residue_name: str,
    residue_number: int,
    insertion_code: str,
) -> str:
    icode = insertion_code.strip()
    suffix = f"{residue_number}{icode}" if icode else str(residue_number)
    return f"{chain_id}_{residue_name}_{suffix}"


def _safe_residue_label(label: str) -> str:
    return "".join(
        character if character.isalnum() or character in {"_", "-"} else "_"
        for character in label
    )


def _load_ptmpsi_protein_class() -> Any:
    try:
        from ptmpsi.protein import Protein
    except ImportError as exc:
        raise ImportError(
            "PTM-Psi is required for capping in this module. Install it in the "
            "active environment, for example from the PTMPSI source repository."
        ) from exc

    return Protein


def _cap_with_ptmpsi(
    *,
    input_pdb_path: Path,
    output_pdb_path: Path,
    chain_id: str,
    add_ace: bool,
    add_nme: bool,
) -> None:
    Protein = _load_ptmpsi_protein_class()

    protein = Protein(
        filename=str(input_pdb_path),
        interactive=False,
        delwat=True,
        delhet=False,
    )

    if add_ace:
        protein.prepend(_normalize_chain_id(chain_id), "ACE")

    if add_nme:
        protein.append(_normalize_chain_id(chain_id), "NME")

    output_pdb_path.parent.mkdir(parents=True, exist_ok=True)
    protein.write_pdb(str(output_pdb_path))

    if not output_pdb_path.exists() or output_pdb_path.stat().st_size == 0:
        raise RuntimeError(f"PTM-Psi did not write capped PDB: {output_pdb_path}")


def _write_single_residue_pdb(
    *,
    input_pdb_path: Path,
    output_pdb_path: Path,
    chain_id: str,
    residue_id: tuple[str, int, str],
) -> Path:
    parsed = _parse_first_model(input_pdb_path)
    output_pdb_path.parent.mkdir(parents=True, exist_ok=True)

    io = PDBIO()
    io.set_structure(parsed.structure)
    io.save(
        str(output_pdb_path),
        select=_SingleResidueSelect(
            model_id=parsed.model.id,
            chain_id=chain_id,
            residue_id=residue_id,
        ),
        preserve_atom_numbering=True,
    )

    if not output_pdb_path.exists() or output_pdb_path.stat().st_size == 0:
        raise RuntimeError(f"Could not write extracted residue PDB: {output_pdb_path}")

    return output_pdb_path


def _write_fragment_pdb(
    *,
    input_pdb_path: Path,
    output_pdb_path: Path,
    chain_id: str,
    residue_ids: set[tuple[str, int, str]],
) -> Path:
    parsed = _parse_first_model(input_pdb_path)
    output_pdb_path.parent.mkdir(parents=True, exist_ok=True)

    io = PDBIO()
    io.set_structure(parsed.structure)
    io.save(
        str(output_pdb_path),
        select=_FragmentSelect(
            model_id=parsed.model.id,
            chain_id=chain_id,
            residue_ids=residue_ids,
        ),
        preserve_atom_numbering=True,
    )

    if not output_pdb_path.exists() or output_pdb_path.stat().st_size == 0:
        raise RuntimeError(f"Could not write fragment PDB: {output_pdb_path}")

    return output_pdb_path


def find_nonstandard_residues(
    input_pdb_path: str | Path,
) -> list[tuple[str, Residue.Residue]]:
    input_pdb_path = Path(input_pdb_path)
    parsed = _parse_first_model(input_pdb_path)

    result: list[tuple[str, Residue.Residue]] = []

    for chain in parsed.model.get_chains():
        chain_id = _normalize_chain_id(str(chain.id))

        for residue in chain.get_residues():
            if not _is_nonstandard_protein_residue(residue):
                continue

            if not _residue_has_backbone_for_capping(residue):
                continue

            result.append((chain_id, residue))

    return result


def build_capped_nonstandard_residue_models(
    *,
    pdb_id: str,
    input_pdb_path: str | Path,
    output_dir: str | Path,
    overwrite: bool = True,
) -> NonstandardCappingSummary:
    input_pdb_path = Path(input_pdb_path).resolve()
    output_dir = Path(output_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    if overwrite and output_dir.exists():
        for child in output_dir.iterdir():
            if child.is_dir():
                shutil.rmtree(child)
            else:
                child.unlink()

    log_path = _build_log_path(input_pdb_path=input_pdb_path, variant_label="nonstandard")

    _print_nonstandard_start_box(
        input_pdb_path=input_pdb_path,
        output_dir=output_dir,
        log_path=log_path,
    )

    append_key_value_block(
        log_path=log_path,
        title="NONSTANDARD CAPPING START",
        key_value_pairs=[
            ("pdb_id", pdb_id),
            ("input_pdb_path", str(input_pdb_path)),
            ("output_dir", str(output_dir)),
            ("backend", "PTM-Psi"),
        ],
    )

    records: list[NonstandardResidueRecord] = []

    for chain_id, residue in find_nonstandard_residues(input_pdb_path):
        residue_name = residue.get_resname().strip().upper()
        residue_number = int(residue.id[1])
        insertion_code = str(residue.id[2]).strip()

        label = _safe_residue_label(
            _residue_label(
                chain_id=chain_id,
                residue_name=residue_name,
                residue_number=residue_number,
                insertion_code=insertion_code,
            )
        )

        residue_dir = output_dir / label
        extracted_pdb_path = residue_dir / "extracted_residue.pdb"
        capped_pdb_path = residue_dir / "capped_model_ace_nme.pdb"
        manifest_path = residue_dir / "manifest.json"

        success = True
        error_message = ""

        try:
            _write_single_residue_pdb(
                input_pdb_path=input_pdb_path,
                output_pdb_path=extracted_pdb_path,
                chain_id=chain_id,
                residue_id=residue.id,
            )

            _cap_with_ptmpsi(
                input_pdb_path=extracted_pdb_path,
                output_pdb_path=capped_pdb_path,
                chain_id=chain_id,
                add_ace=True,
                add_nme=True,
            )

        except Exception as exc:
            success = False
            error_message = f"{type(exc).__name__}: {exc}"

        record = NonstandardResidueRecord(
            pdb_id=pdb_id,
            chain_id=chain_id,
            residue_name=residue_name,
            residue_number=residue_number,
            insertion_code=insertion_code,
            residue_label=label,
            source_pdb_path=input_pdb_path,
            extracted_pdb_path=extracted_pdb_path,
            capped_pdb_path=capped_pdb_path,
            manifest_path=manifest_path,
            backend="PTM-Psi",
            success=success,
            error_message=error_message,
        )

        residue_dir.mkdir(parents=True, exist_ok=True)
        manifest_path.write_text(
            json.dumps(_nonstandard_record_to_json(record), indent=2),
            encoding="utf-8",
        )

        append_key_value_block(
            log_path=log_path,
            title=f"NONSTANDARD RESIDUE {label}",
            key_value_pairs=[
                ("chain_id", chain_id),
                ("residue_name", residue_name),
                ("residue_number", str(residue_number)),
                ("insertion_code", insertion_code or "(none)"),
                ("extracted_pdb_path", str(extracted_pdb_path)),
                ("capped_pdb_path", str(capped_pdb_path)),
                ("success", str(success)),
                ("error_message", error_message or "(none)"),
            ],
        )

        records.append(record)

    summary_manifest_path = output_dir / "nonstandard_capping_manifest.json"

    summary = NonstandardCappingSummary(
        pdb_id=pdb_id,
        input_pdb_path=input_pdb_path,
        output_dir=output_dir,
        n_nonstandard_residues=len(records),
        n_capped_models_written=sum(1 for record in records if record.success),
        records=tuple(records),
        manifest_path=summary_manifest_path,
    )

    summary_manifest_path.write_text(
        json.dumps(_nonstandard_summary_to_json(summary), indent=2),
        encoding="utf-8",
    )

    _append_nonstandard_summary_log(log_path=log_path, summary=summary)
    _print_nonstandard_summary_table(summary)

    return summary


def _nonstandard_record_to_json(record: NonstandardResidueRecord) -> dict[str, Any]:
    payload = asdict(record)
    for key in (
        "source_pdb_path",
        "extracted_pdb_path",
        "capped_pdb_path",
        "manifest_path",
    ):
        payload[key] = str(payload[key])
    return payload


def _nonstandard_summary_to_json(summary: NonstandardCappingSummary) -> dict[str, Any]:
    return {
        "pdb_id": summary.pdb_id,
        "input_pdb_path": str(summary.input_pdb_path),
        "output_dir": str(summary.output_dir),
        "n_nonstandard_residues": summary.n_nonstandard_residues,
        "n_capped_models_written": summary.n_capped_models_written,
        "records": [
            _nonstandard_record_to_json(record) for record in summary.records
        ],
        "manifest_path": str(summary.manifest_path),
    }


def _append_nonstandard_summary_log(
    *,
    log_path: Path,
    summary: NonstandardCappingSummary,
) -> None:
    append_key_value_block(
        log_path=log_path,
        title="NONSTANDARD CAPPING SUMMARY",
        key_value_pairs=[
            ("pdb_id", summary.pdb_id),
            ("input_pdb_path", str(summary.input_pdb_path)),
            ("output_dir", str(summary.output_dir)),
            ("n_nonstandard_residues", str(summary.n_nonstandard_residues)),
            ("n_capped_models_written", str(summary.n_capped_models_written)),
            ("manifest_path", str(summary.manifest_path)),
        ],
    )


def _print_nonstandard_start_box(
    *,
    input_pdb_path: Path,
    output_dir: Path,
    log_path: Path,
) -> None:
    print_mixed_box(
        title=f"{MODULE_NAME} · START",
        line_list=[
            f"Input       : {input_pdb_path}",
            f"Output dir  : {output_dir}",
            f"Backend     : PTM-Psi",
            f"Log file    : {log_path}",
        ],
    )


def _print_nonstandard_summary_table(summary: NonstandardCappingSummary) -> None:
    row_list = [
        ("input", str(summary.input_pdb_path)),
        ("output_dir", str(summary.output_dir)),
        ("nonstandard residues", str(summary.n_nonstandard_residues)),
        ("capped models", str(summary.n_capped_models_written)),
        ("manifest", str(summary.manifest_path)),
    ]

    print_light_table_box(
        title=f"{MODULE_NAME} · nonstandard summary",
        row_list=row_list,
        left_header="Field",
        right_header="Value",
    )


def _has_peptide_bond(
    left_residue: Residue.Residue,
    right_residue: Residue.Residue,
    peptide_bond_max_distance: float,
) -> bool:
    if "C" not in left_residue or "N" not in right_residue:
        return False

    distance = float(left_residue["C"] - right_residue["N"])
    return distance <= peptide_bond_max_distance


def find_internal_fragments(
    chain: Chain.Chain,
    peptide_bond_max_distance: float = 1.8,
) -> list[Fragment]:
    residue_list = list(chain.get_residues())

    polymer_positions = [
        index
        for index, residue in enumerate(residue_list)
        if _is_protein_like_residue(residue)
    ]

    if not polymer_positions:
        return []

    fragment_list: list[Fragment] = []
    current_start_index = polymer_positions[0]
    current_end_index = polymer_positions[0]

    for previous_index, current_index in zip(polymer_positions[:-1], polymer_positions[1:]):
        previous_residue = residue_list[previous_index]
        current_residue = residue_list[current_index]

        if _has_peptide_bond(
            left_residue=previous_residue,
            right_residue=current_residue,
            peptide_bond_max_distance=peptide_bond_max_distance,
        ):
            current_end_index = current_index
            continue

        start_residue = residue_list[current_start_index]
        end_residue = residue_list[current_end_index]

        fragment_list.append(
            Fragment(
                chain_id=_normalize_chain_id(str(chain.id)),
                start_index=current_start_index,
                end_index=current_end_index,
                start_resseq=int(start_residue.id[1]),
                end_resseq=int(end_residue.id[1]),
                start_resname=start_residue.get_resname().strip(),
                end_resname=end_residue.get_resname().strip(),
            )
        )

        current_start_index = current_index
        current_end_index = current_index

    start_residue = residue_list[current_start_index]
    end_residue = residue_list[current_end_index]

    fragment_list.append(
        Fragment(
            chain_id=_normalize_chain_id(str(chain.id)),
            start_index=current_start_index,
            end_index=current_end_index,
            start_resseq=int(start_residue.id[1]),
            end_resseq=int(end_residue.id[1]),
            start_resname=start_residue.get_resname().strip(),
            end_resname=end_residue.get_resname().strip(),
        )
    )

    return fragment_list


def find_internal_gap_boundaries(
    chain: Chain.Chain,
    peptide_bond_max_distance: float = 1.8,
) -> list[InternalGapBoundary]:
    residue_list = list(chain.get_residues())
    fragment_list = find_internal_fragments(
        chain=chain,
        peptide_bond_max_distance=peptide_bond_max_distance,
    )

    if len(fragment_list) < 2:
        return []

    boundary_list: list[InternalGapBoundary] = []

    for left_fragment, right_fragment in zip(fragment_list[:-1], fragment_list[1:]):
        left_residue = residue_list[left_fragment.end_index]
        right_residue = residue_list[right_fragment.start_index]

        boundary_list.append(
            InternalGapBoundary(
                chain_id=_normalize_chain_id(str(chain.id)),
                left_residue_index=left_fragment.end_index,
                right_residue_index=right_fragment.start_index,
                left_resseq=int(left_residue.id[1]),
                right_resseq=int(right_residue.id[1]),
                left_resname=left_residue.get_resname().strip(),
                right_resname=right_residue.get_resname().strip(),
            )
        )

    return boundary_list


def cap_internal_gaps_in_structure(
    input_pdb_path: str | Path,
    output_pdb_path: str | Path,
    peptide_bond_max_distance: float = 1.8,
) -> list[ChainCappingSummary]:
    input_pdb_path = Path(input_pdb_path).resolve()
    output_pdb_path = Path(output_pdb_path).resolve()
    parsed = _parse_first_model(input_pdb_path)

    tmp_dir = output_pdb_path.parent / f".{output_pdb_path.stem}_ptmpsi_fragments"
    if tmp_dir.exists():
        shutil.rmtree(tmp_dir)
    tmp_dir.mkdir(parents=True, exist_ok=True)

    capped_fragment_paths: list[Path] = []
    summary_list: list[ChainCappingSummary] = []

    for chain in parsed.model.get_chains():
        chain_id = _normalize_chain_id(str(chain.id))
        residue_list = list(chain.get_residues())

        fragment_list = find_internal_fragments(
            chain=chain,
            peptide_bond_max_distance=peptide_bond_max_distance,
        )
        boundary_list = find_internal_gap_boundaries(
            chain=chain,
            peptide_bond_max_distance=peptide_bond_max_distance,
        )

        for fragment_index, fragment in enumerate(fragment_list):
            residue_ids = {
                residue_list[index].id
                for index in range(fragment.start_index, fragment.end_index + 1)
            }

            fragment_input = tmp_dir / f"{chain_id}_fragment_{fragment_index + 1:02d}.pdb"
            fragment_output = tmp_dir / f"{chain_id}_fragment_{fragment_index + 1:02d}_capped.pdb"

            _write_fragment_pdb(
                input_pdb_path=input_pdb_path,
                output_pdb_path=fragment_input,
                chain_id=chain_id,
                residue_ids=residue_ids,
            )

            add_ace = fragment_index > 0
            add_nme = fragment_index < len(fragment_list) - 1

            _cap_with_ptmpsi(
                input_pdb_path=fragment_input,
                output_pdb_path=fragment_output,
                chain_id=chain_id,
                add_ace=add_ace,
                add_nme=add_nme,
            )

            capped_fragment_paths.append(fragment_output)

        summary_list.append(
            ChainCappingSummary(
                chain_id=chain_id,
                n_fragments=len(fragment_list),
                n_internal_boundaries=len(boundary_list),
                n_boundaries_capped=len(boundary_list),
            )
        )

    _merge_pdb_fragments_with_biopython(
        input_pdb_paths=capped_fragment_paths,
        output_pdb_path=output_pdb_path,
    )

    return summary_list


def _merge_pdb_fragments_with_biopython(
    *,
    input_pdb_paths: list[Path],
    output_pdb_path: Path,
) -> None:
    if not input_pdb_paths:
        raise ValueError("No capped fragment PDB files were produced.")

    output_structure = Structure.Structure("ptmpsi_capped_structure")
    output_model = Model.Model(0)
    output_structure.add(output_model)

    parser = PDBParser(QUIET=True)
    serial_counter = 0
    chain_counter = 0

    for input_path in input_pdb_paths:
        fragment_structure = parser.get_structure(input_path.stem, str(input_path))
        fragment_model = next(fragment_structure.get_models())

        for source_chain in fragment_model.get_chains():
            output_chain = Chain.Chain(source_chain.id)
            if source_chain.id in {chain.id for chain in output_model.get_chains()}:
                chain_counter += 1
                output_chain = Chain.Chain(chr(ord("A") + (chain_counter % 26)))

            for source_residue in source_chain.get_residues():
                copied_residue = Residue.Residue(
                    id=source_residue.id,
                    resname=source_residue.get_resname(),
                    segid=source_residue.get_segid(),
                )

                for source_atom in source_residue.get_atoms():
                    copied_atom = source_atom.copy()
                    serial_counter += 1
                    copied_atom.serial_number = serial_counter
                    copied_residue.add(copied_atom)

                output_chain.add(copied_residue)

            output_model.add(output_chain)

    output_pdb_path.parent.mkdir(parents=True, exist_ok=True)

    io = PDBIO()
    io.set_structure(output_structure)
    io.save(str(output_pdb_path))


def summarize_capping_results(summary_list: list[ChainCappingSummary]) -> str:
    return "\n".join(
        (
            f"chain={summary.chain_id} "
            f"n_fragments={summary.n_fragments} "
            f"n_internal_boundaries={summary.n_internal_boundaries} "
            f"n_boundaries_capped={summary.n_boundaries_capped}"
        )
        for summary in summary_list
    )


def _print_start_box(
    *,
    input_pdb_path: Path,
    output_pdb_path: Path,
    log_path: Path,
    peptide_bond_max_distance: float,
    variant_label: str | None,
) -> None:
    print_mixed_box(
        title=f"{MODULE_NAME} · internal START",
        line_list=[
            f"Input       : {input_pdb_path}",
            f"Output      : {output_pdb_path}",
            f"Backend     : PTM-Psi",
            f"Log file    : {log_path}",
            f"Variant     : {variant_label or 'standard'}",
            f"C-N max     : {peptide_bond_max_distance}",
        ],
    )


def _print_chain_summary_table(summary_list: list[ChainCappingSummary]) -> None:
    row_list = [
        (
            f"chain {summary.chain_id}",
            (
                f"fragments={summary.n_fragments}, "
                f"boundaries={summary.n_internal_boundaries}, "
                f"capped={summary.n_boundaries_capped}"
            ),
        )
        for summary in summary_list
    ]

    if not row_list:
        row_list = [("none", "no polymer chains found")]

    print_light_table_box(
        title=f"{MODULE_NAME} · chain summary",
        row_list=row_list,
        left_header="Chain",
        right_header="Summary",
    )


def _print_end_box(
    *,
    status: str,
    output_pdb_path: Path,
    total_fragments: int,
    total_boundaries: int,
    total_capped: int,
) -> None:
    print_double_header_box(
        title=f"{MODULE_NAME} · internal END",
        line_list=[
            f"Status      : {status}",
            f"Output PDB  : {output_pdb_path}",
            f"Fragments   : {total_fragments}",
            f"Boundaries  : {total_boundaries}",
            f"Capped      : {total_capped}",
        ],
    )


def _append_run_log(
    *,
    log_path: Path,
    input_pdb_path: Path,
    output_pdb_path: Path,
    peptide_bond_max_distance: float,
    variant_label: str | None,
    summary_list: list[ChainCappingSummary],
) -> None:
    append_key_value_block(
        log_path=log_path,
        title="RUN INFO",
        key_value_pairs=[
            ("module", MODULE_NAME),
            ("backend", "PTM-Psi"),
            ("variant", variant_label or "standard"),
            ("input_pdb_path", str(input_pdb_path)),
            ("output_pdb_path", str(output_pdb_path)),
            ("peptide_bond_max_distance", str(peptide_bond_max_distance)),
            ("chain_count", str(len(summary_list))),
            ("total_fragments", str(sum(x.n_fragments for x in summary_list))),
            (
                "total_internal_boundaries",
                str(sum(x.n_internal_boundaries for x in summary_list)),
            ),
            (
                "total_capped_boundaries",
                str(sum(x.n_boundaries_capped for x in summary_list)),
            ),
        ],
    )

    append_log_text(log_path, "[CHAIN SUMMARIES]")
    if summary_list:
        for summary in summary_list:
            append_log_text(
                log_path,
                (
                    f"chain={summary.chain_id} "
                    f"n_fragments={summary.n_fragments} "
                    f"n_internal_boundaries={summary.n_internal_boundaries} "
                    f"n_boundaries_capped={summary.n_boundaries_capped}"
                ),
            )
    else:
        append_log_text(log_path, "none")
    append_log_text(log_path, "")


def internal_cap_selected_structure(
    *,
    input_pdb_path: str | Path,
    output_pdb_path: str | Path,
    peptide_bond_max_distance: float = 1.8,
    variant_label: str | None = None,
) -> dict[str, str | bool | int | float]:
    input_pdb_path = Path(input_pdb_path).resolve()
    output_pdb_path = Path(output_pdb_path).resolve()
    log_path = _build_log_path(
        input_pdb_path=input_pdb_path,
        variant_label=variant_label,
    )

    _print_start_box(
        input_pdb_path=input_pdb_path,
        output_pdb_path=output_pdb_path,
        log_path=log_path,
        peptide_bond_max_distance=peptide_bond_max_distance,
        variant_label=variant_label,
    )

    summary_list = cap_internal_gaps_in_structure(
        input_pdb_path=input_pdb_path,
        output_pdb_path=output_pdb_path,
        peptide_bond_max_distance=peptide_bond_max_distance,
    )

    _append_run_log(
        log_path=log_path,
        input_pdb_path=input_pdb_path,
        output_pdb_path=output_pdb_path,
        peptide_bond_max_distance=peptide_bond_max_distance,
        variant_label=variant_label,
        summary_list=summary_list,
    )

    _print_chain_summary_table(summary_list)

    total_fragments = sum(x.n_fragments for x in summary_list)
    total_boundaries = sum(x.n_internal_boundaries for x in summary_list)
    total_capped = sum(x.n_boundaries_capped for x in summary_list)

    output_dict: dict[str, str | bool | int | float] = {
        "internal_capping_success": output_pdb_path.is_file()
        and output_pdb_path.stat().st_size > 0,
        "internal_capping_input_path": str(input_pdb_path),
        "internal_capping_output_path": str(output_pdb_path),
        "internal_capping_backend": "PTM-Psi",
        "internal_capping_chain_count": len(summary_list),
        "internal_capping_total_fragments": total_fragments,
        "internal_capping_total_boundaries": total_boundaries,
        "internal_capping_total_capped_boundaries": total_capped,
        "internal_capping_peptide_bond_max_distance": peptide_bond_max_distance,
        "internal_capping_summary_text": summarize_capping_results(summary_list),
    }

    if variant_label is not None:
        output_dict["internal_capping_variant"] = variant_label

    _print_end_box(
        status="success" if bool(output_dict["internal_capping_success"]) else "failed",
        output_pdb_path=output_pdb_path,
        total_fragments=total_fragments,
        total_boundaries=total_boundaries,
        total_capped=total_capped,
    )

    return output_dict


def cap_fragment_termini_for_variant(
    pdb_directory: str | Path,
    pdb_id: str,
    *,
    variant_label: str,
    input_pdb_path: str | Path,
    output_pdb_path: str | Path,
    chain_id: str,
    add_ace: bool,
    add_nme: bool,
) -> dict[str, str | bool | int]:
    """Cap the artificial termini of one gap fragment with PTM-Psi.

    This helper is intentionally narrower than ``internal_cap_selected_structure``.
    The gaps route first splits a disconnected protein into physically connected
    fragments, so each fragment no longer contains an internal break by itself.
    The caller therefore decides from the fragment position whether the fragment
    needs an ACE cap on its N side, an NME cap on its C side, both caps, or no
    caps at all.

    If no cap is requested, the function copies the fragment unchanged to the
    requested output path. This keeps the runner simple because every fragment
    still has one normalized post-capping path before protonation. The function
    does not decide whether a fragment is chemically complete; it only applies
    the explicit cap flags supplied by the gap-fragment route.
    """
    protein_dir = Path(pdb_directory).resolve()
    normalized_pdb_id = str(pdb_id).strip().upper()
    normalized_variant_label = sanitize_variant_label(variant_label)
    input_pdb_path = Path(input_pdb_path).resolve()
    output_pdb_path = Path(output_pdb_path).resolve()
    normalized_chain_id = _normalize_chain_id(chain_id)

    if not input_pdb_path.exists():
        raise FileNotFoundError(f"Fragment PDB not found: {input_pdb_path}")

    log_paths = build_module_log_paths(
        protein_data_dir=protein_dir.parent,
        pdb_id=normalized_pdb_id,
        module_name=MODULE_NAME,
        variant_label=normalized_variant_label,
    )
    log_path = log_paths.module_log_path

    output_pdb_path.parent.mkdir(parents=True, exist_ok=True)

    if add_ace or add_nme:
        _cap_with_ptmpsi(
            input_pdb_path=input_pdb_path,
            output_pdb_path=output_pdb_path,
            chain_id=normalized_chain_id,
            add_ace=add_ace,
            add_nme=add_nme,
        )
    else:
        if input_pdb_path != output_pdb_path:
            shutil.copyfile(input_pdb_path, output_pdb_path)

    success = output_pdb_path.is_file() and output_pdb_path.stat().st_size > 0

    append_key_value_block(
        log_path=log_path,
        title="FRAGMENT TERMINI CAPPING",
        key_value_pairs=[
            ("module", MODULE_NAME),
            ("backend", "PTM-Psi" if (add_ace or add_nme) else "copy"),
            ("pdb_id", normalized_pdb_id),
            ("variant", normalized_variant_label),
            ("input_pdb_path", str(input_pdb_path)),
            ("output_pdb_path", str(output_pdb_path)),
            ("chain_id", normalized_chain_id),
            ("add_ace", str(add_ace)),
            ("add_nme", str(add_nme)),
            ("success", str(success)),
        ],
    )

    return {
        "fragment_capping_success": success,
        "fragment_capping_input_path": str(input_pdb_path),
        "fragment_capping_output_path": str(output_pdb_path),
        "fragment_capping_chain_id": normalized_chain_id,
        "fragment_capping_add_ace": add_ace,
        "fragment_capping_add_nme": add_nme,
        "fragment_capping_n_caps_requested": int(add_ace) + int(add_nme),
        "fragment_capping_variant": normalized_variant_label,
    }


def internal_cap_variant_structure(
    pdb_id: str,
    protein_dir: str | Path,
    *,
    variant_label: str,
    input_pdb_path: str | Path,
    peptide_bond_max_distance: float = 1.8,
) -> dict[str, str | bool | int | float]:
    output_pdb_path = build_variant_internal_capped_output_path(
        pdb_id=pdb_id,
        protein_dir=protein_dir,
        variant_label=variant_label,
    )

    return internal_cap_selected_structure(
        input_pdb_path=input_pdb_path,
        output_pdb_path=output_pdb_path,
        peptide_bond_max_distance=peptide_bond_max_distance,
        variant_label=variant_label,
    )


def internal_cap_protein_structure(
    pdb_id: str,
    protein_dir: str | Path,
    *,
    input_pdb_path: str | Path,
    peptide_bond_max_distance: float = 1.8,
) -> dict[str, str | bool | int | float]:
    output_pdb_path = build_standard_internal_capped_output_path(
        pdb_id=pdb_id,
        protein_dir=protein_dir,
    )

    return internal_cap_selected_structure(
        input_pdb_path=input_pdb_path,
        output_pdb_path=output_pdb_path,
        peptide_bond_max_distance=peptide_bond_max_distance,
        variant_label=None,
    )


def cap_internal_gaps_for_pdb_directory(
    pdb_directory: str | Path,
    pdb_id: str,
    input_pdb_path: str | Path,
    peptide_bond_max_distance: float = 1.8,
) -> tuple[Path, list[ChainCappingSummary]]:
    output_pdb_path = build_standard_internal_capped_output_path(
        protein_dir=pdb_directory,
        pdb_id=pdb_id,
    )

    summary_list = cap_internal_gaps_in_structure(
        input_pdb_path=input_pdb_path,
        output_pdb_path=output_pdb_path,
        peptide_bond_max_distance=peptide_bond_max_distance,
    )

    return output_pdb_path, summary_list


def cap_internal_gaps_for_variant(
    pdb_directory: str | Path,
    pdb_id: str,
    *,
    variant_label: str,
    input_pdb_path: str | Path,
    peptide_bond_max_distance: float = 1.8,
) -> dict[str, str | bool | int | float]:
    return internal_cap_variant_structure(
        pdb_id=pdb_id,
        protein_dir=pdb_directory,
        variant_label=variant_label,
        input_pdb_path=input_pdb_path,
        peptide_bond_max_distance=peptide_bond_max_distance,
    )


def cap_nonstandard_residues_for_pdb_directory(
    pdb_directory: str | Path,
    pdb_id: str,
    input_pdb_path: str | Path,
    overwrite: bool = True,
) -> NonstandardCappingSummary:
    output_dir = build_nonstandard_capping_output_dir(
        pdb_id=pdb_id,
        protein_dir=pdb_directory,
    )

    return build_capped_nonstandard_residue_models(
        pdb_id=pdb_id,
        input_pdb_path=input_pdb_path,
        output_dir=output_dir,
        overwrite=overwrite,
    )


if __name__ == "__main__":
    import argparse

    argument_parser = argparse.ArgumentParser(
        description=(
            "Extract protein-like non-standard residues and build ACE/NME capped "
            "model compounds with PTM-Psi."
        )
    )
    argument_parser.add_argument(
        "pdb_directory",
        type=Path,
        help="Protein directory, e.g. data/proteins/1ABC",
    )
    argument_parser.add_argument(
        "pdb_id",
        type=str,
        help="PDB ID, e.g. 1ABC",
    )
    argument_parser.add_argument(
        "input_pdb",
        type=Path,
        help="Input representative-unit or protein PDB path.",
    )
    argument_parser.add_argument(
        "--internal-gaps",
        action="store_true",
        help="Also run PTM-Psi-based internal fragment capping.",
    )
    argument_parser.add_argument(
        "--peptide-bond-max-distance",
        type=float,
        default=1.8,
        help="Maximum C-N distance treated as peptide-bonded. Default: 1.8 Å",
    )
    arguments = argument_parser.parse_args()

    nonstandard_summary = cap_nonstandard_residues_for_pdb_directory(
        pdb_directory=arguments.pdb_directory,
        pdb_id=arguments.pdb_id,
        input_pdb_path=arguments.input_pdb,
    )

    print(f"nonstandard_manifest={nonstandard_summary.manifest_path}")
    print(f"n_nonstandard_residues={nonstandard_summary.n_nonstandard_residues}")
    print(f"n_capped_models_written={nonstandard_summary.n_capped_models_written}")

    if arguments.internal_gaps:
        result = internal_cap_protein_structure(
            pdb_id=arguments.pdb_id,
            protein_dir=arguments.pdb_directory,
            input_pdb_path=arguments.input_pdb,
            peptide_bond_max_distance=arguments.peptide_bond_max_distance,
        )
        print(f"internal_capping_output={result['internal_capping_output_path']}")
        print(result["internal_capping_summary_text"])