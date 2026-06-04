"""Internal helpers for cap.py — dataclasses, BioPython utilities, and capping logic."""
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
from stack_protein_preparation._cap_paths import (
    build_nonstandard_capping_output_dir,
    sanitize_variant_label,
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

_BACKBONE_NONSTD_RESNAMES: frozenset[str] = frozenset({
    "MSE", "SEC", "PYL",
    "SEP", "TPO", "PTR",
    "S1P", "T1P", "Y1P", "NEP",
    "CSO", "CSD", "OCS", "CSS", "CME", "SCY", "SMC", "P1L", "CAF",
    "MHO", "SME", "OMT",
    "TYS", "NIY", "TYI", "IYR", "TPQ", "TYQ",
    "MLZ", "MLY", "M3L", "ALY", "KCR", "SLL", "KCX", "LYZ", "LLP",
    "AGM", "DA2", "MMO", "MMA",
    "CIR",
    "CGU",
    "HYP", "4FB",
    "4HT", "TRO", "HTR",
    "HIC",
    "PCA",
    "SAR", "MAA", "MEA", "MEN",
    "SAC",
    "CRO", "DYG",
    "H2D", "H1D", "H2E", "H1E",
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
    """BioPython selector for one residue in one chain from the first model."""

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
    """BioPython selector for one connected residue fragment."""

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


def _is_polymer_residue(residue: Residue.Residue) -> bool:
    """Return True if the residue is a standard polymer (ATOM record, not water/HETATM)."""
    hetflag = residue.get_id()[0]
    if hetflag not in {" ", ""}:
        return False
    resname = residue.get_resname().strip().upper()
    if resname in {"HOH", "WAT", "SOL"}:
        return False
    return True


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

    output_structure = Structure.Structure("capped_structure")
    output_model = Model.Model(0)
    output_structure.add(output_model)

    serial_counter = [0]
    summary_list: list[ChainCappingSummary] = []

    for chain in parsed.model.get_chains():
        new_chain, summary = cap_internal_gap_boundaries(
            chain=chain,
            serial_counter=serial_counter,
            peptide_bond_max_distance=peptide_bond_max_distance,
        )
        output_model.add(new_chain)
        summary_list.append(summary)

    output_pdb_path.parent.mkdir(parents=True, exist_ok=True)
    io = PDBIO()
    io.set_structure(output_structure)
    io.save(str(output_pdb_path))

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


def _build_nme_residue(
    chain_id: str,
    resseq: int,
    serial_counter: list[int],
    reference_residue: "Residue.Residue",
) -> "Residue.Residue":
    """Build a minimal NME capping residue placed near the C-terminus of the previous residue."""
    import numpy as np

    nme = Residue.Residue((" ", resseq, " "), "NME", "")

    c_pos = None
    if "C" in reference_residue:
        c_atom = reference_residue["C"]
        c_pos = c_atom.get_vector().get_array()

    if c_pos is None:
        for atom in reference_residue.get_atoms():
            c_pos = atom.get_vector().get_array()
            break

    if c_pos is None:
        c_pos = np.array([0.0, 0.0, 0.0])

    offset = np.array([1.5, 0.0, 0.0])

    from Bio.PDB.Atom import Atom as _Atom

    for atom_name, delta in [("N", offset), ("CH3", offset + np.array([1.5, 0.0, 0.0]))]:
        serial_counter[0] += 1
        pos = c_pos + delta
        atom = _Atom(
            name=atom_name,
            coord=pos,
            bfactor=0.0,
            occupancy=1.0,
            altloc=" ",
            fullname=f" {atom_name:<3}",
            serial_number=serial_counter[0],
            element="N" if atom_name == "N" else "C",
        )
        nme.add(atom)

    return nme


def _build_ace_residue(
    chain_id: str,
    resseq: int,
    serial_counter: list[int],
    reference_residue: "Residue.Residue",
) -> "Residue.Residue":
    """Build a minimal ACE capping residue placed near the N-terminus of the next residue."""
    import numpy as np

    ace = Residue.Residue((" ", resseq, " "), "ACE", "")

    n_pos = None
    if "N" in reference_residue:
        n_atom = reference_residue["N"]
        n_pos = n_atom.get_vector().get_array()

    if n_pos is None:
        for atom in reference_residue.get_atoms():
            n_pos = atom.get_vector().get_array()
            break

    if n_pos is None:
        n_pos = np.array([0.0, 0.0, 0.0])

    offset = np.array([-1.5, 0.0, 0.0])

    from Bio.PDB.Atom import Atom as _Atom

    for atom_name, delta in [
        ("CH3", offset + np.array([-1.5, 0.0, 0.0])),
        ("C", offset),
        ("O", offset + np.array([0.0, 1.2, 0.0])),
    ]:
        serial_counter[0] += 1
        pos = n_pos + delta
        atom = _Atom(
            name=atom_name,
            coord=pos,
            bfactor=0.0,
            occupancy=1.0,
            altloc=" ",
            fullname=f" {atom_name:<3}",
            serial_number=serial_counter[0],
            element="C" if atom_name in {"CH3", "C"} else "O",
        )
        ace.add(atom)

    return ace


def cap_internal_gap_boundaries(
    chain: "Chain.Chain",
    serial_counter: list[int],
    peptide_bond_max_distance: float = 1.8,
) -> tuple["Chain.Chain", ChainCappingSummary]:
    """Cap internal gap boundaries in a BioPython chain with ACE/NME residues.

    This is a pure BioPython implementation that does not require PTM-Psi.
    For each internal gap boundary the C-terminal fragment receives an NME cap
    and the N-terminal fragment of the next block receives an ACE cap.

    Parameters
    ----------
    chain
        BioPython Chain to process (modified in-place via a new Chain).
    serial_counter
        Single-element list used as a mutable integer counter for atom serials.
    peptide_bond_max_distance
        Maximum C-N distance (Angstrom) considered a peptide bond.

    Returns
    -------
    tuple[Chain, ChainCappingSummary]
        Modified chain and capping summary.
    """
    chain_id = _normalize_chain_id(str(chain.id))
    residue_list = [r for r in chain.get_residues() if _is_polymer_residue(r)]

    fragments = find_internal_fragments(chain, peptide_bond_max_distance=peptide_bond_max_distance)
    boundaries = find_internal_gap_boundaries(chain, peptide_bond_max_distance=peptide_bond_max_distance)

    new_chain = Chain.Chain(chain.id)

    all_residues = list(chain.get_residues())
    frag_residue_sets: list[list[Residue.Residue]] = []
    for frag in fragments:
        frag_residues = [
            all_residues[i]
            for i in range(frag.start_index, frag.end_index + 1)
            if _is_polymer_residue(all_residues[i])
        ]
        frag_residue_sets.append(frag_residues)

    n_boundaries_capped = 0
    resseq_counter = 1

    for frag_idx, frag_residues in enumerate(frag_residue_sets):
        is_first_frag = frag_idx == 0
        is_last_frag = frag_idx == len(frag_residue_sets) - 1

        if not is_first_frag and frag_residues:
            ace = _build_ace_residue(
                chain_id=chain_id,
                resseq=resseq_counter,
                serial_counter=serial_counter,
                reference_residue=frag_residues[0],
            )
            resseq_counter += 1
            new_chain.add(ace)
            n_boundaries_capped += 1

        for res_idx, res in enumerate(frag_residues):
            is_first_in_frag = res_idx == 0
            is_last_in_frag = res_idx == len(frag_residues) - 1

            new_res = Residue.Residue((" ", resseq_counter, " "), res.get_resname(), res.get_segid())
            resseq_counter += 1

            for atom in res.get_atoms():
                atom_name = atom.get_name().strip()
                if is_last_in_frag and not is_last_frag and atom_name == "OXT":
                    continue
                if is_first_in_frag and not is_first_frag and atom_name == "H1":
                    continue
                serial_counter[0] += 1
                copied = atom.copy()
                copied.serial_number = serial_counter[0]
                new_res.add(copied)

            new_chain.add(new_res)

        if not is_last_frag and frag_residues:
            nme = _build_nme_residue(
                chain_id=chain_id,
                resseq=resseq_counter,
                serial_counter=serial_counter,
                reference_residue=frag_residues[-1],
            )
            resseq_counter += 1
            new_chain.add(nme)

    summary = ChainCappingSummary(
        chain_id=chain_id,
        n_fragments=len(fragments),
        n_internal_boundaries=len(boundaries),
        n_boundaries_capped=n_boundaries_capped,
    )

    return new_chain, summary


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
