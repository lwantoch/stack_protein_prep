#!/usr/bin/env python3
"""Stage FRUTON prepared structures for input to gbsa-grub / gbsa-pipeline.

For each protein with a completed prepared structure, this script:
  1. Writes {PDBID}.pdb — the prepared structure with the cocrystal ligand
     and crystal waters removed (waters go to {PDBID}_WAT.pdb).
  2. Writes {PDBID}_WAT.pdb — HOH/WAT/SOL residues extracted from the
     prepared structure for pre-insertion into the solvation box.
  3. Converts each biochemical cofactor (e.g. FAD, NAD, ADP) from
     components/{PDBID}_cofactor.pdb to one SDF per molecule in cofactors/.
     FRUTON's _cofactor.pdb is treated as authoritative — all molecules in it
     are staged as cofactors.
  4. Stages phosaa14SB.xml for proteins with phosphorylated residues.
  5. Writes a per-variant manifest.json.

Output layout (matches generate_grid.sh):
  {output_dir}/{dir_id}/{dir_id}.pdb          # protein without waters
  {output_dir}/{dir_id}/{dir_id}_WAT.pdb      # crystal waters (if any)
  {output_dir}/{dir_id}/cofactors/*.sdf       # cofactors (from _cofactor.pdb, if non-empty)
  {output_dir}/{dir_id}/ligands/              # populated separately with DEKOIS SDFs
  {output_dir}/{dir_id}/extra_ff/             # phosaa14SB.xml if needed

Where dir_id is {PDBID} for default/single variants, {PDBID}_{variant}
for named variants (gap, alpha, mod).

Usage
-----
    pixi run python scripts/stage_for_gbsa.py \\
        --proteins-dir data/proteins \\
        --output-dir /path/to/proteins

    # Restrict to specific PDB IDs:
    pixi run python scripts/stage_for_gbsa.py \\
        --proteins-dir data/proteins \\
        --output-dir /path/to/proteins \\
        --pdb-ids 1A5H 3OLL
"""
from __future__ import annotations

import argparse
import json
import logging
import shutil
import sys
from pathlib import Path

logger = logging.getLogger(__name__)

_WATER_RESNAMES: frozenset[str] = frozenset({"HOH", "WAT", "SOL", "DOD"})


# ---------------------------------------------------------------------------
# PDB helpers
# ---------------------------------------------------------------------------

def _resnames_in_pdb(path: Path) -> set[str]:
    names: set[str] = set()
    if not path.is_file():
        return names
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        if line[:6] in ("ATOM  ", "HETATM") and len(line) >= 20:
            names.add(line[17:20].strip().upper())
    return names


def _partition_pdb_lines(
    prepared_pdb: Path,
    ligand_resnames: set[str],
) -> tuple[list[str], list[str]]:
    """Split prepared PDB into (protein_lines, water_lines).

    Cocrystal-ligand HETATM records and crystal waters are removed from the
    protein output; waters are collected separately. Crystal waters,
    metals, and backbone non-standard residues (PTR, SEP …) that are not
    cocrystal ligands remain in the protein output.
    """
    protein_lines: list[str] = []
    water_lines: list[str] = []

    for raw in prepared_pdb.read_text(encoding="utf-8", errors="replace").splitlines():
        line = raw.rstrip()
        if not line:
            continue
        rec = line[:6]
        if rec in ("ATOM  ", "HETATM") and len(line) >= 20:
            resname = line[17:20].strip().upper()
            if rec == "HETATM" and resname in ligand_resnames:
                continue
            if resname in _WATER_RESNAMES:
                water_lines.append(line)
                continue
        protein_lines.append(line)

    if protein_lines and protein_lines[-1] != "END":
        protein_lines.append("END")

    return protein_lines, water_lines


def _write_water_pdb(water_lines: list[str], path: Path) -> None:
    lines = water_lines[:]
    if lines and lines[-1] != "END":
        lines.append("END")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


# ---------------------------------------------------------------------------
# Cofactor PDB → SDF (one SDF per distinct molecule instance)
# ---------------------------------------------------------------------------

def _group_cofactor_molecules(cofactor_pdb: Path) -> dict[str, list[str]]:
    """Group cofactor PDB lines by (resname, chain, resseq) molecule key."""
    molecules: dict[str, list[str]] = {}
    for raw in cofactor_pdb.read_text(encoding="utf-8", errors="replace").splitlines():
        line = raw.rstrip()
        if not line or line[:6] not in ("ATOM  ", "HETATM") or len(line) < 26:
            continue
        resname = line[17:20].strip().upper()
        chain = line[21:22].strip() or "_"
        try:
            resseq = int(line[22:26].strip())
        except ValueError:
            resseq = 0
        key = f"{resname}_{chain}_{resseq}"
        molecules.setdefault(key, []).append(line)
    return molecules


def _pdb_lines_to_sdf(pdb_lines: list[str], sdf_path: Path) -> bool:
    """Convert PDB atom lines to SDF using RDKit. Returns True on success.

    RDKit infers bond orders from atom connectivity and valence rules.
    sanitize=False + explicit SanitizeMol lets partial failures surface as
    warnings rather than silently returning None.
    """
    try:
        from rdkit import Chem
    except ImportError:
        logger.warning("rdkit not available — cannot convert %s", sdf_path.name)
        return False

    pdb_str = "\n".join(pdb_lines) + "\nEND\n"

    mol = Chem.MolFromPDBBlock(pdb_str, removeHs=False, sanitize=True)
    if mol is None:
        mol = Chem.MolFromPDBBlock(pdb_str, removeHs=False, sanitize=False)
        if mol is None:
            logger.warning("RDKit could not parse PDB block for %s", sdf_path.name)
            return False
        Chem.SanitizeMol(mol, catchErrors=True)

    try:
        writer = Chem.SDWriter(str(sdf_path))
        writer.write(mol)
        writer.close()
    except Exception as exc:
        logger.warning("RDKit SDF write failed for %s: %s", sdf_path.name, exc)
        return False

    return sdf_path.is_file() and sdf_path.stat().st_size > 0


# ---------------------------------------------------------------------------
# Gaussian / MCPB / RESP force-field file collection
# ---------------------------------------------------------------------------

def _collect_gaussian_ff_files(pdb_dir: Path, gaussian_manifest_path: str | None) -> list[Path]:
    """Return resolved paths to frcmod and mol2 files from a successful gaussian_params_manifest."""
    if not gaussian_manifest_path:
        return []
    mp = Path(gaussian_manifest_path)
    if not mp.is_file():
        return []
    try:
        manifest = json.loads(mp.read_text(encoding="utf-8"))
    except Exception:
        return []
    if manifest.get("status") != "success":
        return []

    files: list[Path] = []

    def _resolve(raw: str) -> Path | None:
        p = Path(raw)
        resolved = p if p.is_absolute() else pdb_dir / p
        return resolved if resolved.is_file() else None

    for site in manifest.get("mcpb_result", {}).get("site_results", []):
        if site.get("status") != "success":
            continue
        if site.get("frcmod"):
            p = _resolve(site["frcmod"])
            if p:
                files.append(p)
        for mol2 in site.get("mol2_files", []):
            p = _resolve(mol2)
            if p:
                files.append(p)

    for res in manifest.get("resp_result", {}).get("residue_results", []):
        if res.get("status") != "success":
            continue
        for key in ("frcmod_path", "mol2_path"):
            if res.get(key):
                p = _resolve(res[key])
                if p:
                    files.append(p)

    return files


# ---------------------------------------------------------------------------
# phosaa14SB detection and XML staging
# ---------------------------------------------------------------------------

def _has_phosaa14sb(nonstd_manifest_path: str | None) -> bool:
    if not nonstd_manifest_path:
        return False
    p = Path(nonstd_manifest_path)
    if not p.is_file():
        return False
    try:
        manifest = json.loads(p.read_text(encoding="utf-8"))
    except Exception:
        return False
    return any(
        str(res.get("database", "")).strip().lower() == "phosaa14sb"
        for res in manifest.get("residues", [])
    )


def _find_phosaa14sb_xml() -> Path | None:
    try:
        import openmmforcefields as _off
        p = Path(_off.__file__).parent / "ffxml" / "amber" / "phosaa14SB.xml"
        if p.is_file():
            return p
    except ImportError:
        pass
    # Fallback: known gbsa-grub pixi environment path on this cluster
    fallback = (
        Path("/mnt/netapp2/Store_uni/home/otras/hcx/lwa/repos/gbsa-grub")
        / ".pixi/envs/default/lib/python3.12/site-packages"
        / "openmmforcefields/ffxml/amber/phosaa14SB.xml"
    )
    return fallback if fallback.is_file() else None


# ---------------------------------------------------------------------------
# Directory ID helper
# ---------------------------------------------------------------------------

def _dir_id(pdb_id: str, variant: str) -> str:
    """Return the output directory base name for a given protein variant.

    Single/default variants use just the PDB ID; named variants (gap, alpha,
    mod) append the variant as a suffix so generate_grid.sh can treat each
    model as an independent protein directory.
    """
    if not variant or variant.lower() in ("default", "single"):
        return pdb_id
    return f"{pdb_id}_{variant}"


# ---------------------------------------------------------------------------
# Stage one variant
# ---------------------------------------------------------------------------

def stage_variant(
    *,
    pdb_id: str,
    variant: str,
    prepared_pdb: Path,
    components_dir: Path,
    nonstd_manifest_path: str | None,
    gaussian_ff_files: list[Path],
    out_dir: Path,
    dir_id: str,
) -> dict:
    """Stage one prepared variant into out_dir. Returns the manifest dict."""
    out_dir.mkdir(parents=True, exist_ok=True)

    # --- {dir_id}.pdb and {dir_id}_WAT.pdb --------------------------------
    ligand_resnames = _resnames_in_pdb(components_dir / f"{pdb_id}_ligand.pdb")
    protein_lines, water_lines = _partition_pdb_lines(prepared_pdb, ligand_resnames)

    protein_pdb_out = out_dir / f"{dir_id}.pdb"
    protein_pdb_out.write_text("\n".join(protein_lines) + "\n", encoding="utf-8")
    logger.info("[%s] %s written (%d lines)", dir_id, protein_pdb_out.name, len(protein_lines))

    crystal_waters_pdb_out: Path | None = None
    if water_lines:
        crystal_waters_pdb_out = out_dir / f"{dir_id}_WAT.pdb"
        _write_water_pdb(water_lines, crystal_waters_pdb_out)
        logger.info("[%s] %s written (%d water records)", dir_id, crystal_waters_pdb_out.name, len(water_lines))

    # --- cofactor SDFs in cofactors/ ---------------------------------------
    # Stage whatever FRUTON put in _cofactor.pdb — that file is FRUTON's
    # authoritative classification. If it has atoms, they are cofactors.
    cofactor_sdf_paths: list[str] = []
    cofactor_component = components_dir / f"{pdb_id}_cofactor.pdb"
    if cofactor_component.is_file() and cofactor_component.stat().st_size > 0:
        molecules = _group_cofactor_molecules(cofactor_component)
        if molecules:
            cofactors_dir = out_dir / "cofactors"
            cofactors_dir.mkdir(exist_ok=True)
            for mol_key, mol_lines in molecules.items():
                sdf_out = cofactors_dir / f"{mol_key}.sdf"
                if _pdb_lines_to_sdf(mol_lines, sdf_out):
                    cofactor_sdf_paths.append(str(sdf_out.resolve()))
                    logger.info("[%s] cofactor SDF: %s", dir_id, sdf_out.name)
                else:
                    fallback = cofactors_dir / f"{mol_key}.pdb"
                    fallback.write_text("\n".join(mol_lines) + "\nEND\n", encoding="utf-8")
                    logger.warning(
                        "[%s] SDF conversion failed for %s — PDB saved: %s",
                        dir_id, mol_key, fallback.name,
                    )

    # --- extra_ff/: phosaa14SB.xml + MCPB/RESP frcmod + mol2 files -------
    extra_ff_files: list[str] = []

    def _stage_extra_ff(src: Path, label: str) -> None:
        nonlocal extra_ff_files
        extra_ff_dir = out_dir / "extra_ff"
        extra_ff_dir.mkdir(exist_ok=True)
        dst = extra_ff_dir / src.name
        if not dst.exists():
            shutil.copy2(src, dst)
        dst_str = str(dst.resolve())
        if dst_str not in extra_ff_files:
            extra_ff_files.append(dst_str)
        logger.info("[%s] extra_ff: %s", dir_id, label)

    if _has_phosaa14sb(nonstd_manifest_path):
        xml_src = _find_phosaa14sb_xml()
        if xml_src is not None:
            _stage_extra_ff(xml_src, "phosaa14SB.xml")
        else:
            logger.warning(
                "[%s] phosaa14SB detected but XML not found — check openmmforcefields install",
                dir_id,
            )

    for ff_file in gaussian_ff_files:
        _stage_extra_ff(ff_file, ff_file.name)

    # --- manifest.json ----------------------------------------------------
    manifest: dict = {
        "pdb_id": pdb_id,
        "variant": variant,
        "dir_id": dir_id,
        "protein_pdb": str(protein_pdb_out.resolve()),
        "crystal_waters_pdb": str(crystal_waters_pdb_out.resolve()) if crystal_waters_pdb_out else None,
        "cofactor_sdfs": cofactor_sdf_paths,
        "extra_ff_files": extra_ff_files,
        "has_cofactors": len(cofactor_sdf_paths) > 0,
        "has_phosaa14sb": len(extra_ff_files) > 0,
    }
    (out_dir / "manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
    )
    return manifest


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def run(proteins_dir: Path, output_dir: Path, pdb_ids: list[str] | None, force: bool = False) -> None:
    pipeline_json = proteins_dir / "pipeline.json"
    if not pipeline_json.is_file():
        sys.exit(f"ERROR: pipeline.json not found at {pipeline_json}")

    records: list[dict] = json.loads(pipeline_json.read_text(encoding="utf-8"))

    if pdb_ids:
        requested = {p.upper() for p in pdb_ids}
        records = [r for r in records if r.get("pdb_id", "").upper() in requested]

    all_manifests: list[dict] = []
    n_staged = n_skipped = 0

    for rec in records:
        pdb_id = rec.get("pdb_id", "").strip().upper()
        if not pdb_id:
            continue

        status = rec.get("prepared_structure.status")
        if status != "success" and not force:
            logger.debug("Skipping %s — prepared_structure.status=%s", pdb_id, status)
            n_skipped += 1
            continue

        pdb_dir = Path(rec["pdb_directory"])
        components_dir = pdb_dir / "components"
        nonstd_manifest_path = rec.get("nonstd_residue_params.manifest_path") or None
        gauss_manifest_path = rec.get("gaussian_params.manifest_path") or None
        gaussian_ff_files = _collect_gaussian_ff_files(pdb_dir, gauss_manifest_path)

        variant_str = rec.get("prepared_structure.variant", "") or ""
        variants = [v.strip() for v in variant_str.split(",") if v.strip()]

        if not variants:
            prepared_dir = pdb_dir / "prepared"
            variants = [
                d.name
                for d in prepared_dir.iterdir()
                if d.is_dir() and (d / f"{pdb_id}.pdb").is_file()
            ] if prepared_dir.is_dir() else []
            if variants:
                logger.debug("[%s] variant(s) discovered from filesystem: %s", pdb_id, variants)

        for variant in variants:
            prepared_pdb = pdb_dir / "prepared" / variant / f"{pdb_id}.pdb"
            if not prepared_pdb.is_file():
                logger.warning("Missing prepared PDB: %s", prepared_pdb)
                n_skipped += 1
                continue

            did = _dir_id(pdb_id, variant)
            out_dir = output_dir / did
            try:
                manifest = stage_variant(
                    pdb_id=pdb_id,
                    variant=variant,
                    prepared_pdb=prepared_pdb,
                    components_dir=components_dir,
                    nonstd_manifest_path=nonstd_manifest_path,
                    gaussian_ff_files=gaussian_ff_files,
                    out_dir=out_dir,
                    dir_id=did,
                )
                all_manifests.append(manifest)
                n_staged += 1
            except Exception as exc:
                logger.error("Failed to stage %s (%s): %s", pdb_id, variant, exc, exc_info=True)
                n_skipped += 1

    output_dir.mkdir(parents=True, exist_ok=True)
    combined_path = output_dir / "staging_manifest.json"
    combined_path.write_text(json.dumps(all_manifests, indent=2) + "\n", encoding="utf-8")
    print(f"Staged {n_staged} variant(s), skipped {n_skipped}.")
    print(f"Combined manifest → {combined_path}")


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Stage FRUTON prepared structures for gbsa-grub input."
    )
    parser.add_argument(
        "--proteins-dir", type=Path, required=True,
        help="proteins data directory (contains pipeline.json)",
    )
    parser.add_argument(
        "--output-dir", type=Path, required=True,
        help="Directory where staged files are written (the proteins/ dir for generate_grid.sh)",
    )
    parser.add_argument(
        "--pdb-ids", nargs="*",
        help="Subset of PDB IDs to stage (default: all with success status)",
    )
    parser.add_argument(
        "--force", action="store_true",
        help="Stage proteins even if prepared_structure.status != success",
    )
    parser.add_argument("--verbose", action="store_true")
    args = parser.parse_args(argv)

    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(levelname)s: %(message)s",
    )

    if not args.proteins_dir.is_dir():
        print(f"ERROR: proteins-dir not found: {args.proteins_dir}", file=sys.stderr)
        return 1

    run(
        proteins_dir=args.proteins_dir,
        output_dir=args.output_dir,
        pdb_ids=args.pdb_ids,
        force=args.force,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
