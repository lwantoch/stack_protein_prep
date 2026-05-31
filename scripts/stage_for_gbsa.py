#!/usr/bin/env python3
"""Stage FRUTON prepared structures for input to gbsa-grub / gbsa-pipeline.

For each protein with a completed prepared structure, this script:
  1. Writes protein.pdb — the prepared structure with the cocrystal ligand
     removed (replaced by the docked pose in gbsa-grub).
  2. Converts each biochemical cofactor (e.g. FAD, NAD) from
     components/{PDBID}_cofactor.pdb to one SDF file per molecule.
  3. Stages phosaa14SB.xml for proteins with phosphorylated residues.
  4. Writes a per-variant manifest.json with absolute paths ready for
     constructing gbsa-grub statepoints.

Usage
-----
    pixi run python scripts/stage_for_gbsa.py \\
        --proteins-dir data/proteins \\
        --output-dir /path/to/gbsa_staging

    # Restrict to specific PDB IDs:
    pixi run python scripts/stage_for_gbsa.py \\
        --proteins-dir data/proteins \\
        --output-dir /path/to/gbsa_staging \\
        --pdb-ids 1A5H 3OLL
"""
from __future__ import annotations

import argparse
import json
import logging
import shutil
import sys
import tempfile
from pathlib import Path

logger = logging.getLogger(__name__)


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


def _strip_ligand_from_pdb(prepared_pdb: Path, ligand_resnames: set[str]) -> list[str]:
    """Return PDB lines with cocrystal-ligand HETATM records removed.

    Crystal waters, metals, and backbone non-standard residues (PTR, SEP …)
    are intentionally kept — gbsa-pipeline extracts waters separately and
    OpenMM's AMBER FF handles simple ions via amber14-all.xml.
    """
    out: list[str] = []
    for raw in prepared_pdb.read_text(encoding="utf-8", errors="replace").splitlines():
        line = raw.rstrip()
        if not line:
            continue
        rec = line[:6]
        if rec == "HETATM" and len(line) >= 20:
            resname = line[17:20].strip().upper()
            if resname in ligand_resnames:
                continue
        out.append(line)
    if out and out[-1] != "END":
        out.append("END")
    return out


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
    """Convert PDB atom lines to SDF via PyMOL. Returns True on success."""
    try:
        import pymol2
    except ImportError:
        logger.warning("pymol2 not available — cannot convert %s", sdf_path.name)
        return False

    pdb_str = "\n".join(pdb_lines) + "\nEND\n"
    with tempfile.NamedTemporaryFile(
        suffix=".pdb", mode="w", encoding="utf-8", delete=False
    ) as tmp:
        tmp.write(pdb_str)
        tmp_path = Path(tmp.name)

    try:
        with pymol2.PyMOL() as p:
            p.cmd.load(str(tmp_path), "cof")
            p.cmd.save(str(sdf_path), "cof")
        return sdf_path.is_file() and sdf_path.stat().st_size > 0
    except Exception as exc:
        logger.warning("PyMOL SDF conversion failed for %s: %s", sdf_path.name, exc)
        return False
    finally:
        tmp_path.unlink(missing_ok=True)


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
# Stage one variant
# ---------------------------------------------------------------------------

def stage_variant(
    *,
    pdb_id: str,
    variant: str,
    prepared_pdb: Path,
    components_dir: Path,
    nonstd_manifest_path: str | None,
    out_dir: Path,
) -> dict:
    """Stage one prepared variant into out_dir. Returns the manifest dict."""
    out_dir.mkdir(parents=True, exist_ok=True)

    # --- protein.pdb --------------------------------------------------
    ligand_resnames = _resnames_in_pdb(components_dir / f"{pdb_id}_ligand.pdb")
    protein_lines = _strip_ligand_from_pdb(prepared_pdb, ligand_resnames)
    protein_pdb_out = out_dir / "protein.pdb"
    protein_pdb_out.write_text("\n".join(protein_lines) + "\n", encoding="utf-8")
    logger.info("[%s/%s] protein.pdb written (%d lines)", pdb_id, variant, len(protein_lines))

    # --- cofactor SDFs ------------------------------------------------
    cofactor_sdf_paths: list[str] = []
    cofactor_component = components_dir / f"{pdb_id}_cofactor.pdb"
    if cofactor_component.is_file() and cofactor_component.stat().st_size > 0:
        molecules = _group_cofactor_molecules(cofactor_component)
        for mol_key, mol_lines in molecules.items():
            sdf_out = out_dir / f"cofactor_{mol_key}.sdf"
            if _pdb_lines_to_sdf(mol_lines, sdf_out):
                cofactor_sdf_paths.append(str(sdf_out.resolve()))
                logger.info("[%s/%s] cofactor SDF: %s", pdb_id, variant, sdf_out.name)
            else:
                # Write PDB fallback for manual conversion
                fallback = out_dir / f"cofactor_{mol_key}.pdb"
                fallback.write_text("\n".join(mol_lines) + "\nEND\n", encoding="utf-8")
                logger.warning(
                    "[%s/%s] SDF conversion failed for %s — PDB saved: %s",
                    pdb_id, variant, mol_key, fallback.name,
                )

    # --- phosaa14SB.xml -----------------------------------------------
    extra_ff_files: list[str] = []
    if _has_phosaa14sb(nonstd_manifest_path):
        xml_src = _find_phosaa14sb_xml()
        if xml_src is not None:
            extra_ff_dir = out_dir / "extra_ff"
            extra_ff_dir.mkdir(exist_ok=True)
            dst = extra_ff_dir / "phosaa14SB.xml"
            if not dst.exists():
                shutil.copy2(xml_src, dst)
            extra_ff_files.append(str(dst.resolve()))
            logger.info("[%s/%s] phosaa14SB.xml staged", pdb_id, variant)
        else:
            logger.warning(
                "[%s/%s] phosaa14SB detected but XML not found — check openmmforcefields install",
                pdb_id, variant,
            )

    # --- manifest.json ------------------------------------------------
    manifest: dict = {
        "pdb_id": pdb_id,
        "variant": variant,
        "protein_pdb": str(protein_pdb_out.resolve()),
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

        variant_str = rec.get("prepared_structure.variant", "") or ""
        variants = [v.strip() for v in variant_str.split(",") if v.strip()]

        # If the pipeline table has no variant info (e.g. status=failed but file
        # was written), discover variants from the prepared/ directory on disk.
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

            out_dir = output_dir / pdb_id / variant
            try:
                manifest = stage_variant(
                    pdb_id=pdb_id,
                    variant=variant,
                    prepared_pdb=prepared_pdb,
                    components_dir=components_dir,
                    nonstd_manifest_path=nonstd_manifest_path,
                    out_dir=out_dir,
                )
                all_manifests.append(manifest)
                n_staged += 1
            except Exception as exc:
                logger.error("Failed to stage %s/%s: %s", pdb_id, variant, exc, exc_info=True)
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
        help="Directory where staged files are written",
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
