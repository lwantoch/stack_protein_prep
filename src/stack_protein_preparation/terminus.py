# /home/grheco/repositorios/stack_protein_prep/src/stack_protein_preparation/terminus.py

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Tuple

STANDARD_AMINO_ACID_RESNAMES = {
    "ALA",
    "ARG",
    "ASN",
    "ASP",
    "ASH",
    "CYS",
    "CYM",
    "CYX",
    "GLN",
    "GLU",
    "GLH",
    "GLY",
    "HIS",
    "HID",
    "HIE",
    "HIP",
    "ILE",
    "LEU",
    "LYS",
    "LYN",
    "MET",
    "PHE",
    "PRO",
    "SER",
    "THR",
    "TRP",
    "TYR",
    "VAL",
}

PROTEIN_LIKE_NONSTANDARD_RESNAMES = {
    "MSE",
    "SEC",
    "PYL",
    "SEP",
    "TPO",
    "PTR",
    "CSO",
    "CSD",
    "OCS",
    "KCX",
    "LLP",
    "MLY",
    "HYP",
}

POLYMER_PROTEIN_RESNAMES = (
    STANDARD_AMINO_ACID_RESNAMES | PROTEIN_LIKE_NONSTANDARD_RESNAMES
)

INTERNAL_CAP_RESNAMES = {"ACE", "NME"}

N_TERMINAL_AMBER_NAME_MAP = {
    "ALA": "NALA",
    "ARG": "NARG",
    "ASN": "NASN",
    "ASP": "NASP",
    "ASH": "NASH",
    "CYS": "NCYS",
    "CYM": "NCYM",
    "CYX": "NCYX",
    "GLN": "NGLN",
    "GLU": "NGLU",
    "GLH": "NGLH",
    "GLY": "NGLY",
    "HIS": "NHIS",
    "HID": "NHID",
    "HIE": "NHIE",
    "HIP": "NHIP",
    "ILE": "NILE",
    "LEU": "NLEU",
    "LYS": "NLYS",
    "LYN": "NLYN",
    "MET": "NMET",
    "PHE": "NPHE",
    "PRO": "NPRO",
    "SER": "NSER",
    "THR": "NTHR",
    "TRP": "NTRP",
    "TYR": "NTYR",
    "VAL": "NVAL",
}

C_TERMINAL_AMBER_NAME_MAP = {
    "ALA": "CALA",
    "ARG": "CARG",
    "ASN": "CASN",
    "ASP": "CASP",
    "ASH": "CASH",
    "CYS": "CCYS",
    "CYM": "CCYM",
    "CYX": "CCYX",
    "GLN": "CGLN",
    "GLU": "CGLU",
    "GLH": "CGLH",
    "GLY": "CGLY",
    "HIS": "CHIS",
    "HID": "CHID",
    "HIE": "CHIE",
    "HIP": "CHIP",
    "ILE": "CILE",
    "LEU": "CLEU",
    "LYS": "CLYS",
    "LYN": "CLYN",
    "MET": "CMET",
    "PHE": "CPHE",
    "PRO": "CPRO",
    "SER": "CSER",
    "THR": "CTHR",
    "TRP": "CTRP",
    "TYR": "CTYR",
    "VAL": "CVAL",
}


@dataclass(frozen=True)
class TerminusConversionResult:
    chain_id: str
    first_resseq: int
    last_resseq: int
    first_old_resname: str
    first_new_resname: str
    last_old_resname: str
    last_new_resname: str
    single_residue_chain: bool


@dataclass
class AtomRecord:
    record_name: str
    serial: int
    atom_name: str
    altloc: str
    resname: str
    chain_id: str
    resseq: int
    icode: str
    x: float
    y: float
    z: float
    occupancy: float
    tempfactor: float
    element: str
    charge: str
    original_line: str

    @property
    def residue_key(self) -> Tuple[str, int, str]:
        return (self.chain_id, self.resseq, self.icode)

    def format_pdb_line(self) -> str:
        record_name = f"{self.record_name:<6}"
        atom_name = _format_atom_name(self.atom_name, self.element)
        altloc = (self.altloc or " ")[:1]
        resname = f"{self.resname:>4}"[:4]
        chain_id = (self.chain_id or " ")[:1]
        icode = (self.icode or " ")[:1]
        element = f"{(self.element or '').strip():>2}"[:2]
        charge = f"{(self.charge or '').strip():>2}"[:2]

        return (
            f"{record_name}"
            f"{self.serial:5d} "
            f"{atom_name}"
            f"{altloc}"
            f"{resname}"
            f"{chain_id}"
            f"{self.resseq:4d}"
            f"{icode}   "
            f"{self.x:8.3f}"
            f"{self.y:8.3f}"
            f"{self.z:8.3f}"
            f"{self.occupancy:6.2f}"
            f"{self.tempfactor:6.2f}"
            f"          "
            f"{element}"
            f"{charge}"
        )


@dataclass
class RawLine:
    line: str


ParsedEntry = AtomRecord | RawLine


def sanitize_variant_label(variant_label: str) -> str:
    cleaned = "".join(
        character if character.isalnum() or character in {"_", "-"} else "_"
        for character in variant_label.strip()
    )
    cleaned = cleaned.strip("_")
    if not cleaned:
        raise ValueError(f"Invalid empty variant label derived from {variant_label!r}")
    return cleaned


def build_default_terminus_output_path(
    pdb_directory: Path | str,
    pdb_id: str,
) -> Path:
    pdb_directory = Path(pdb_directory)
    return pdb_directory / "components" / f"{pdb_id}_protein_amber_termini.pdb"


def build_variant_terminus_output_path(
    pdb_directory: Path | str,
    pdb_id: str,
    variant_label: str,
) -> Path:
    pdb_directory = Path(pdb_directory)
    safe_variant_label = sanitize_variant_label(variant_label)
    return (
        pdb_directory
        / "components"
        / f"{pdb_id}_protein_amber_termini_{safe_variant_label}.pdb"
    )


def convert_protein_termini_to_amber(
    input_pdb_path: Path | str,
    output_pdb_path: Path | str,
) -> List[TerminusConversionResult]:
    input_pdb_path = Path(input_pdb_path)
    output_pdb_path = Path(output_pdb_path)

    parsed_entries = _read_pdb_entries(input_pdb_path)
    atom_records = [entry for entry in parsed_entries if isinstance(entry, AtomRecord)]

    residue_map = _collect_residue_records(atom_records)
    chain_polymer_residue_keys = _collect_polymer_residue_keys_by_chain(residue_map)

    terminus_result_list: List[TerminusConversionResult] = []
    n_terminal_keys: set[Tuple[str, int, str]] = set()
    c_terminal_keys: set[Tuple[str, int, str]] = set()

    for chain_id, residue_key_list in chain_polymer_residue_keys.items():
        if not residue_key_list:
            continue

        first_key = residue_key_list[0]
        last_key = residue_key_list[-1]

        first_record = residue_map[first_key][0]
        last_record = residue_map[last_key][0]

        first_old_resname = first_record.resname
        last_old_resname = last_record.resname
        single_residue_chain = first_key == last_key

        if single_residue_chain:
            first_new_resname = first_old_resname
            last_new_resname = last_old_resname
        else:
            first_new_resname = N_TERMINAL_AMBER_NAME_MAP.get(
                first_old_resname,
                first_old_resname,
            )
            last_new_resname = C_TERMINAL_AMBER_NAME_MAP.get(
                last_old_resname,
                last_old_resname,
            )

        _apply_residue_name_to_key(
            residue_map=residue_map,
            residue_key=first_key,
            new_resname=first_new_resname,
        )
        if last_key != first_key:
            _apply_residue_name_to_key(
                residue_map=residue_map,
                residue_key=last_key,
                new_resname=last_new_resname,
            )

        terminus_result_list.append(
            TerminusConversionResult(
                chain_id=chain_id,
                first_resseq=first_record.resseq,
                last_resseq=last_record.resseq,
                first_old_resname=first_old_resname,
                first_new_resname=first_new_resname,
                last_old_resname=last_old_resname,
                last_new_resname=last_new_resname,
                single_residue_chain=single_residue_chain,
            )
        )

        n_terminal_keys.add(first_key)
        c_terminal_keys.add(last_key)

    new_atom_records = _normalize_terminal_atom_names(
        atom_records=atom_records,
        n_terminal_keys=n_terminal_keys,
        c_terminal_keys=c_terminal_keys,
        residue_map=residue_map,
    )

    parsed_entries.extend(new_atom_records)

    output_pdb_path.parent.mkdir(parents=True, exist_ok=True)
    _write_pdb_entries(parsed_entries, output_pdb_path)

    return terminus_result_list


def convert_variant_protein_termini_to_amber(
    pdb_id: str,
    pdb_directory: Path | str,
    *,
    variant_label: str,
    input_pdb_path: Path | str,
) -> dict[str, str | bool | int]:
    pdb_directory = Path(pdb_directory)
    input_pdb_path = Path(input_pdb_path)
    output_pdb_path = build_variant_terminus_output_path(
        pdb_directory=pdb_directory,
        pdb_id=pdb_id,
        variant_label=variant_label,
    )

    result_list = convert_protein_termini_to_amber(
        input_pdb_path=input_pdb_path,
        output_pdb_path=output_pdb_path,
    )

    return {
        "amber_termini_success": output_pdb_path.is_file()
        and output_pdb_path.stat().st_size > 0,
        "amber_termini_input_path": str(input_pdb_path.resolve()),
        "amber_termini_output_path": str(output_pdb_path.resolve()),
        "amber_termini_variant": variant_label,
        "amber_termini_chain_count": len(result_list),
    }


def _read_pdb_entries(pdb_path: Path) -> List[ParsedEntry]:
    parsed_entries: List[ParsedEntry] = []

    for line in pdb_path.read_text(encoding="utf-8").splitlines():
        if line.startswith(("ATOM  ", "HETATM")):
            parsed_entries.append(_parse_atom_line(line))
        else:
            parsed_entries.append(RawLine(line=line))

    return parsed_entries


def _write_pdb_entries(
    parsed_entries: Iterable[ParsedEntry],
    output_pdb_path: Path,
) -> None:
    output_line_list: List[str] = []

    for entry in parsed_entries:
        if isinstance(entry, AtomRecord):
            output_line_list.append(entry.format_pdb_line())
        else:
            output_line_list.append(entry.line)

    output_pdb_path.write_text(
        "\n".join(output_line_list) + "\n",
        encoding="utf-8",
    )


def _parse_atom_line(line: str) -> AtomRecord:
    padded = line.rstrip("\n")
    if len(padded) < 80:
        padded = padded.ljust(80)

    record_name = padded[0:6].strip()
    serial_text = padded[6:11].strip()
    atom_name = padded[12:16]
    altloc = padded[16]
    resname = padded[17:21].strip()
    chain_id = padded[21]
    resseq_text = padded[22:26].strip()
    icode = padded[26]
    x_text = padded[30:38].strip()
    y_text = padded[38:46].strip()
    z_text = padded[46:54].strip()
    occupancy_text = padded[54:60].strip()
    tempfactor_text = padded[60:66].strip()
    element = padded[76:78].strip()
    charge = padded[78:80].strip()

    return AtomRecord(
        record_name=record_name,
        serial=int(serial_text) if serial_text else 0,
        atom_name=atom_name.strip(),
        altloc=altloc if altloc != " " else "",
        resname=resname,
        chain_id=chain_id if chain_id != " " else "",
        resseq=int(resseq_text) if resseq_text else 0,
        icode=icode if icode != " " else "",
        x=float(x_text) if x_text else 0.0,
        y=float(y_text) if y_text else 0.0,
        z=float(z_text) if z_text else 0.0,
        occupancy=float(occupancy_text) if occupancy_text else 1.0,
        tempfactor=float(tempfactor_text) if tempfactor_text else 0.0,
        element=element or _infer_element_from_atom_name(atom_name.strip()),
        charge=charge,
        original_line=line,
    )


def _collect_residue_records(
    atom_records: Iterable[AtomRecord],
) -> Dict[Tuple[str, int, str], List[AtomRecord]]:
    residue_map: Dict[Tuple[str, int, str], List[AtomRecord]] = {}

    for atom_record in atom_records:
        residue_map.setdefault(atom_record.residue_key, []).append(atom_record)

    return residue_map


def _collect_polymer_residue_keys_by_chain(
    residue_map: Dict[Tuple[str, int, str], List[AtomRecord]],
) -> Dict[str, List[Tuple[str, int, str]]]:
    chain_polymer_residue_keys: Dict[str, List[Tuple[str, int, str]]] = {}

    for residue_key, residue_atom_records in residue_map.items():
        representative = residue_atom_records[0]
        resname = representative.resname
        chain_id = representative.chain_id

        if resname in INTERNAL_CAP_RESNAMES:
            continue

        if resname not in POLYMER_PROTEIN_RESNAMES:
            continue

        chain_polymer_residue_keys.setdefault(chain_id, []).append(residue_key)

    for chain_id in chain_polymer_residue_keys:
        chain_polymer_residue_keys[chain_id].sort(key=lambda key: (key[1], key[2]))

    return chain_polymer_residue_keys


def _apply_residue_name_to_key(
    residue_map: Dict[Tuple[str, int, str], List[AtomRecord]],
    residue_key: Tuple[str, int, str],
    new_resname: str,
) -> None:
    for atom_record in residue_map.get(residue_key, []):
        atom_record.resname = new_resname


def _normalize_terminal_atom_names(
    atom_records: Iterable[AtomRecord],
    n_terminal_keys: set[Tuple[str, int, str]],
    c_terminal_keys: set[Tuple[str, int, str]],
    residue_map: Dict[Tuple[str, int, str], List[AtomRecord]],
) -> List[AtomRecord]:
    residue_to_atoms: Dict[Tuple[str, int, str], List[AtomRecord]] = {}
    for atom_record in atom_records:
        residue_to_atoms.setdefault(atom_record.residue_key, []).append(atom_record)

    serial_counter = [max((atom.serial for atom in atom_records), default=0)]
    single_residue_keys = n_terminal_keys & c_terminal_keys
    new_atom_records: List[AtomRecord] = []

    for residue_key in n_terminal_keys:
        residue_atom_records = residue_to_atoms.get(residue_key, [])
        _normalize_n_terminal_hydrogens(residue_atom_records)

    for residue_key in c_terminal_keys:
        residue_atom_records = residue_to_atoms.get(residue_key, [])
        new_oxt_atom = _normalize_c_terminal_oxygens(
            residue_atom_records,
            serial_counter=serial_counter,
            add_oxt_if_missing=residue_key in single_residue_keys,
        )
        if new_oxt_atom is not None:
            residue_map.setdefault(residue_key, []).append(new_oxt_atom)
            new_atom_records.append(new_oxt_atom)

    for residue_key in single_residue_keys:
        residue_atom_records = residue_to_atoms.get(residue_key, [])
        _ensure_single_residue_terminal_hydrogen(residue_atom_records)

    return new_atom_records


def _normalize_n_terminal_hydrogens(
    residue_atom_records: List[AtomRecord],
) -> None:
    n_hydrogen_records: List[AtomRecord] = []

    for atom_record in residue_atom_records:
        stripped_name = atom_record.atom_name.strip().upper()

        if stripped_name in {
            "H",
            "H1",
            "H2",
            "H3",
            "1H",
            "2H",
            "3H",
            "HT1",
            "HT2",
            "HT3",
        }:
            n_hydrogen_records.append(atom_record)

    if not n_hydrogen_records:
        return

    canonical_name_order = ["H1", "H2", "H3"]
    normalized_names: List[str] = []

    for atom_record in n_hydrogen_records:
        current = atom_record.atom_name.strip().upper()
        if current in {"H1", "H2", "H3"}:
            normalized_names.append(current)
        elif current in {"1H", "HT1"}:
            normalized_names.append("H1")
        elif current in {"2H", "HT2"}:
            normalized_names.append("H2")
        elif current in {"3H", "HT3"}:
            normalized_names.append("H3")
        elif current == "H":
            normalized_names.append("H")
        else:
            normalized_names.append(current)

    unique_non_plain_h = [name for name in normalized_names if name != "H"]
    if len(n_hydrogen_records) == 1:
        if normalized_names[0] != "H":
            n_hydrogen_records[0].atom_name = "H"
        return

    assigned_name_list: List[str] = []

    for preferred_name in canonical_name_order:
        if (
            preferred_name in unique_non_plain_h
            and preferred_name not in assigned_name_list
        ):
            assigned_name_list.append(preferred_name)

    for preferred_name in canonical_name_order:
        if len(assigned_name_list) >= len(n_hydrogen_records):
            break
        if preferred_name not in assigned_name_list:
            assigned_name_list.append(preferred_name)

    for atom_record, assigned_name in zip(n_hydrogen_records, assigned_name_list):
        atom_record.atom_name = assigned_name


def _normalize_c_terminal_oxygens(
    residue_atom_records: List[AtomRecord],
    *,
    serial_counter: List[int],
    add_oxt_if_missing: bool = False,
) -> AtomRecord | None:
    oxygen_like_records: List[AtomRecord] = []

    for atom_record in residue_atom_records:
        if _infer_element_from_atom_name(atom_record.atom_name) != "O":
            continue

        stripped_name = atom_record.atom_name.strip().upper()
        if stripped_name in {"O", "O1", "OT1", "OXT", "O2", "OT2"}:
            oxygen_like_records.append(atom_record)

    if not oxygen_like_records:
        return None

    if len(oxygen_like_records) == 1:
        oxygen_like_records[0].atom_name = "O"
        if add_oxt_if_missing:
            return _duplicate_atom_as_oxt(
                source_atom=oxygen_like_records[0],
                serial_counter=serial_counter,
            )
        return None

    primary_assigned = False
    secondary_assigned = False

    for atom_record in oxygen_like_records:
        stripped_name = atom_record.atom_name.strip().upper()

        if stripped_name in {"O", "O1", "OT1"} and not primary_assigned:
            atom_record.atom_name = "O"
            primary_assigned = True
            continue

        if stripped_name in {"OXT", "O2", "OT2"} and not secondary_assigned:
            atom_record.atom_name = "OXT"
            secondary_assigned = True
            continue

    for atom_record in oxygen_like_records:
        stripped_name = atom_record.atom_name.strip().upper()

        if not primary_assigned:
            atom_record.atom_name = "O"
            primary_assigned = True
            continue

        if not secondary_assigned and stripped_name != "O":
            atom_record.atom_name = "OXT"
            secondary_assigned = True
            continue

    seen_o = False
    seen_oxt = False
    for atom_record in oxygen_like_records:
        stripped_name = atom_record.atom_name.strip().upper()
        if stripped_name == "O":
            if not seen_o:
                seen_o = True
            elif not seen_oxt:
                atom_record.atom_name = "OXT"
                seen_oxt = True
        elif stripped_name == "OXT":
            seen_oxt = True

    if add_oxt_if_missing and not any(
        atom.atom_name.strip().upper() == "OXT" for atom in residue_atom_records
    ):
        first_o = next(
            (
                atom
                for atom in residue_atom_records
                if atom.atom_name.strip().upper() == "O"
            ),
            None,
        )
        if first_o is not None:
            return _duplicate_atom_as_oxt(
                source_atom=first_o,
                serial_counter=serial_counter,
            )

    return None


def _ensure_single_residue_terminal_hydrogen(
    residue_atom_records: List[AtomRecord],
) -> None:
    hydrogen_records = [
        atom
        for atom in residue_atom_records
        if atom.atom_name.strip().upper() in {"H", "H1", "1H", "HT1"}
    ]

    if not hydrogen_records:
        return

    hydrogen_records[0].atom_name = "H1"


def _duplicate_atom_as_oxt(
    source_atom: AtomRecord,
    serial_counter: List[int],
) -> AtomRecord:
    serial_counter[0] += 1
    duplicated = AtomRecord(
        record_name=source_atom.record_name,
        serial=serial_counter[0],
        atom_name="OXT",
        altloc=source_atom.altloc,
        resname=source_atom.resname,
        chain_id=source_atom.chain_id,
        resseq=source_atom.resseq,
        icode=source_atom.icode,
        x=source_atom.x,
        y=source_atom.y,
        z=source_atom.z,
        occupancy=source_atom.occupancy,
        tempfactor=source_atom.tempfactor,
        element="O",
        charge=source_atom.charge,
        original_line=source_atom.original_line,
    )
    return duplicated


def _format_atom_name(atom_name: str, element: str) -> str:
    atom_name = atom_name.strip()
    if not atom_name:
        return "    "

    inferred_element = (
        (element or _infer_element_from_atom_name(atom_name)).strip().upper()
    )

    if len(atom_name) >= 4:
        return atom_name[:4]

    if len(inferred_element) == 1 and len(atom_name) < 4:
        return f"{atom_name:>4}"

    return f"{atom_name:<4}"


def _infer_element_from_atom_name(atom_name: str) -> str:
    stripped = atom_name.strip()
    if not stripped:
        return ""

    letters_only = "".join(character for character in stripped if character.isalpha())
    if not letters_only:
        return ""

    letters_only = letters_only.upper()

    two_letter_elements = {
        "ZN",
        "FE",
        "MG",
        "MN",
        "CL",
        "BR",
        "NA",
        "CA",
        "CU",
        "CO",
        "NI",
        "CD",
        "HG",
    }

    if len(letters_only) >= 2 and letters_only[:2] in two_letter_elements:
        return letters_only[:2]

    return letters_only[0]
