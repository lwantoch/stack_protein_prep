"""PyMOL launcher discovery and rendering helpers for FRUTON reports."""
from __future__ import annotations

import os
import subprocess
import sys
import tempfile
from pathlib import Path

# ---------------------------------------------------------------------------
# PyMOL launcher discovery
# ---------------------------------------------------------------------------

# Env vars required by the Debian pymol package (wrapper script values).
_DEBIAN_PYMOL_ENV: dict[str, str] = {
    "PYMOL_PATH":    "/usr/share/pymol",
    "PYMOL_DATA":    "/usr/share/pymol/data",
    "PYMOL_SCRIPTS": "/usr/share/pymol/scripts",
    "CHEMPY_DATA":   "/usr/share/pymol/data/chempy",
}

_pymol_cmd_cache: list[str] | None = None   # [python3_path, "-m", "pymol.__init__"]
_pymol_env_cache: dict[str, str] | None = None


def _find_pymol_launcher() -> tuple[list[str], dict[str, str]] | None:
    """
    Return ([python3, '-m', 'pymol.__init__'], extra_env) for a working PyMOL,
    or None if no working interpreter is found.

    Tries the current interpreter, system /usr/bin/python3, and 'python3' on
    PATH; each with and without the Debian env-var block.  Result is cached.
    """
    global _pymol_cmd_cache, _pymol_env_cache
    if _pymol_cmd_cache is not None:
        return _pymol_cmd_cache, _pymol_env_cache  # type: ignore[return-value]

    candidates = [sys.executable, "/usr/bin/python3", "python3"]
    for py in candidates:
        for extra in ({}, _DEBIAN_PYMOL_ENV):
            env = {**os.environ, **extra}
            try:
                r = subprocess.run(
                    [py, "-c", "import pymol"],
                    env=env, capture_output=True, timeout=15,
                )
                if r.returncode == 0:
                    cmd = [py, "-m", "pymol.__init__"]
                    _pymol_cmd_cache = cmd
                    _pymol_env_cache = extra
                    return cmd, extra
            except (FileNotFoundError, subprocess.TimeoutExpired):
                pass

    _pymol_cmd_cache = []   # sentinel: searched, not found
    _pymol_env_cache = {}
    return None


def _run_pml(pml_lines: list[str], output_png: Path, timeout: int = 180) -> str | None:
    """
    Write *pml_lines* to a temp .pml file, run PyMOL headless, and check that
    *output_png* was created.  Returns an error string or None on success.
    """
    launcher = _find_pymol_launcher()
    if launcher is None:
        return "PyMOL not found in any Python environment."
    cmd, extra_env = launcher

    env = {**os.environ, **extra_env}

    with tempfile.NamedTemporaryFile(
        mode="w", suffix=".pml", delete=False, encoding="utf-8"
    ) as tmp:
        tmp.write("\n".join(pml_lines) + "\n")
        pml_path = Path(tmp.name)

    try:
        subprocess.run(
            cmd + ["-cq", str(pml_path)],
            env=env, capture_output=True, timeout=timeout,
        )
    except subprocess.TimeoutExpired:
        return f"PyMOL rendering timed out after {timeout} s."
    except Exception as exc:
        return f"PyMOL subprocess error: {exc}"
    finally:
        pml_path.unlink(missing_ok=True)

    if not output_png.exists():
        return "PyMOL ran but did not produce the expected PNG."
    return None


# ---------------------------------------------------------------------------
# Non-standard residue helpers
# ---------------------------------------------------------------------------

def _parse_resp_mol2(mol2_path: Path) -> list[dict]:
    """Parse @<TRIPOS>ATOM section of a mol2 file, returning per-atom dicts.

    Each dict has keys: name (str), sybyl_type (str), charge (float).
    """
    atoms: list[dict] = []
    in_atom = False
    try:
        text = mol2_path.read_text(encoding="utf-8", errors="replace")
    except OSError:
        return atoms
    for line in text.splitlines():
        stripped = line.strip()
        if stripped == "@<TRIPOS>ATOM":
            in_atom = True
            continue
        if stripped.startswith("@<TRIPOS>") and in_atom:
            break
        if in_atom and stripped:
            parts = stripped.split()
            if len(parts) >= 9:
                try:
                    atoms.append({
                        "name":       parts[1],
                        "sybyl_type": parts[5],
                        "charge":     float(parts[8]),
                    })
                except (ValueError, IndexError):
                    pass
    return atoms


def _render_nonstd_mol2_figure(mol2_path: Path, output_png: Path) -> str | None:
    """Render non-standard residue mol2 with RESP charge labels via headless PyMOL.

    Atoms are coloured by partial charge (red=negative -> white=neutral -> blue=positive).
    Heavy atoms are labelled with their atom name and RESP charge value.
    Labels have a white opaque background for readability.
    """
    pml = [
        f"load {mol2_path}, mol",
        "bg_color white",
        "hide everything, mol",
        "show sticks, mol",
        # Small spheres on heavy atoms so labels have a clear anchor point
        "show spheres, mol and not elem H",
        "set sphere_scale, 0.18, mol",
        "set stick_radius, 0.10",
        # Colour by RESP partial charge: red=negative, white=~0, blue=positive
        "spectrum partial_charge, red_white_blue, mol",
        # Label heavy atoms: atom name + charge on one line
        'label mol and not elem H, "%s %.3f" % (name, partial_charge)',
        # White opaque label background
        "set label_bg_color, white",
        "set label_bg_transparency, 0.0",
        "set label_color, black",
        "set label_size, 13",
        "set label_font_id, 7",
        "set label_position, (0, 0, 0.8)",
        "orient mol",
        "zoom mol, 2.5",
        "set ray_opaque_background, 1",
        "set antialias, 2",
        "set ray_shadows, 0",
        "set ray_trace_mode, 0",
        f"png {output_png}, width=1200, height=900, dpi=150, ray=1",
        "quit",
    ]
    return _run_pml(pml, output_png)


# ---------------------------------------------------------------------------
# PyMOL figure renderers
# ---------------------------------------------------------------------------

def _render_pymol_figure(
    protein_dir: Path,
    pdb_id: str,
    pipeline_record: dict[str, str],
    output_png: Path,
) -> str | None:
    """
    Render a PNG overview of the starting structure via headless PyMOL.

    Searches for {pdb_id}_delins.pdb first (post-insertion-code step), then
    the raw downloaded PDB.  Returns an error string on failure, None on success.
    """
    # Prefer the delins PDB (cleanest starting point post-renumbering)
    candidates = sorted(protein_dir.rglob(f"*{pdb_id.upper()}_delins.pdb"))
    if not candidates:
        candidates = sorted(protein_dir.rglob(f"*{pdb_id.lower()}_delins.pdb"))
    if not candidates:
        candidates = sorted(protein_dir.rglob(f"{pdb_id.upper()}.pdb"))
    if not candidates:
        candidates = sorted(protein_dir.rglob(f"{pdb_id.lower()}.pdb"))
    if not candidates:
        return f"No input PDB found for {pdb_id} in {protein_dir}."

    input_pdb = candidates[0]
    has_ligands = str(pipeline_record.get("has_ligands", "")).strip().lower()
    center_on_ligand = has_ligands in ("yes", "true")

    pml = [
        f"load {input_pdb}, s",
        "bg_color white",
        "hide everything, s",
        "set_color fruton_navy, [0.157, 0.196, 0.353]",
        "set_color fruton_gold, [0.541, 0.451, 0.059]",
        "set_color fruton_red, [0.902, 0.000, 0.078]",
        "set_color fruton_softgrey, [0.720, 0.735, 0.770]",
        # Cartoon style settings for a clean, illustrative look
        "set cartoon_fancy_helices, 1",
        "set cartoon_smooth_loops, 1",
        "set cartoon_tube_radius, 0.35",
        "set cartoon_loop_radius, 0.20",
        "set specular, 0",
        "show cartoon, s and polymer",
        "color fruton_navy, s and ss h",
        "color fruton_gold, s and ss s",
        "color fruton_softgrey, s and ss l+''",
    ]
    if center_on_ligand:
        pml += [
            "select lig, s and hetatm and not (resn HOH or resn WAT)",
            "show sticks, lig",
            "color atomic, lig and not elem C",
            "color fruton_red, lig and elem C",
            "zoom lig, 12",
        ]
    else:
        pml += ["orient s", "zoom s"]
    pml += [
        "set ray_opaque_background, 1",
        "set antialias, 2",
        "set ray_shadows, 1",
        "set ray_trace_mode, 3",
        "set ray_trace_color, black",
        "set ray_trace_gain, 1.5",
        f"png {output_png}, width=1200, height=900, dpi=150, ray=1",
        "quit",
    ]
    return _run_pml(pml, output_png)


def _ion_element(ion_type: str) -> str:
    """
    Extract a 1-2 character PDB element symbol from an ion type string.

    Examples: 'Ca2+' -> 'CA', 'Zn2+' -> 'ZN', 'Fe' -> 'FE', 'MG' -> 'MG'.
    Used to build a reliable PyMOL ``elem`` selection that avoids the
    CA (alpha-carbon vs calcium) ambiguity of atom-name-based selections.
    """
    import re
    clean = re.sub(r"[0-9+\-\s]", "", ion_type).upper()
    return clean[:2] if clean else ion_type.upper()[:2]


def _find_metal_h_contacts(
    pdb_path: Path,
    donor_cutoff: float = 2.8,
    h_cutoff: float = 2.5,
) -> list[dict]:
    """Return H atoms pointing toward each metal ion using BioPython NeighborSearch.

    For each metal ion found in *pdb_path*, locates all H atoms within
    *h_cutoff* A.  An H atom at distance M-H < M-O (its parent oxygen) is
    pointing toward the metal and should be flagged for removal or reorientation
    before MCPB.py / Gaussian runs.

    Returns a list of dicts with keys:
        metal_label, metal_elem, h_name, h_res, h_chain, h_resnum,
        mh_dist, mo_dist, parent_o_name, status
    """
    try:
        from Bio.PDB import PDBParser, NeighborSearch
        import numpy as np
    except ImportError:
        return []

    try:
        parser = PDBParser(QUIET=True)
        struct = parser.get_structure("s", str(pdb_path))
    except Exception:
        return []

    all_atoms = list(struct.get_atoms())
    ns = NeighborSearch(all_atoms)
    results: list[dict] = []

    metal_elements = {
        "MG", "CA", "ZN", "FE", "MN", "CU", "CO", "NI", "NA", "K",
        "SR", "BA", "CD", "HG", "CR", "V", "MO", "W", "AL",
    }

    for metal_atom in all_atoms:
        if metal_atom.element.upper() not in metal_elements:
            continue
        metal_pos = np.array(metal_atom.get_vector().get_array(), dtype="d")
        metal_res = metal_atom.get_parent()
        metal_chain = metal_res.get_parent().id
        metal_label = (
            f"{metal_atom.element} {metal_chain}{metal_res.id[1]} {metal_atom.name}"
        )

        near = ns.search(metal_pos, max(donor_cutoff, h_cutoff) + 0.5, "A")
        # Build quick lookup of O/N distances from this metal
        donor_dist: dict[tuple, float] = {}
        for a in near:
            if a.element.upper() in {"N", "O", "S"}:
                d = float(np.linalg.norm(
                    np.array(a.get_vector().get_array(), dtype="d") - metal_pos
                ))
                if d <= donor_cutoff:
                    key = (a.get_parent().get_parent().id,
                           a.get_parent().id[1], a.name)
                    donor_dist[key] = d

        for h_atom in near:
            if h_atom.element.upper() != "H":
                continue
            mh = float(np.linalg.norm(
                np.array(h_atom.get_vector().get_array(), dtype="d") - metal_pos
            ))
            if mh > h_cutoff:
                continue

            h_res = h_atom.get_parent()
            h_chain = h_res.get_parent().id
            # Find parent heavy atom (O or N within 1.2 A of this H)
            h_pos = np.array(h_atom.get_vector().get_array(), dtype="d")
            parent_candidates = [
                a for a in ns.search(h_pos, 1.2, "A")
                if a is not h_atom and a.element.upper() in {"N", "O", "S"}
            ]
            if not parent_candidates:
                parent_o_name = "?"
                mo_dist = float("nan")
            else:
                parent_o = min(
                    parent_candidates,
                    key=lambda a: np.linalg.norm(
                        np.array(a.get_vector().get_array(), dtype="d") - h_pos
                    ),
                )
                parent_o_name = parent_o.name
                mo_dist = float(np.linalg.norm(
                    np.array(parent_o.get_vector().get_array(), dtype="d") - metal_pos
                ))

            status = "toward_metal" if (not np.isnan(mo_dist) and mh < mo_dist) else "near_metal"
            results.append({
                "metal_label": metal_label,
                "metal_elem": metal_atom.element.upper(),
                "h_name": h_atom.name,
                "h_res": h_res.resname,
                "h_chain": h_chain,
                "h_resnum": h_res.id[1],
                "mh_dist": round(mh, 3),
                "mo_dist": round(mo_dist, 3) if not np.isnan(mo_dist) else None,
                "parent_o_name": parent_o_name,
                "status": status,
            })

    results.sort(key=lambda r: (r["metal_label"], r["mh_dist"]))
    return results


def _render_metal_pocket_figures(
    pdb_path: Path,
    ion_type: str,
    output_dir: Path,
    pdb_id: str,
) -> dict[str, Path | None]:
    """
    Render two metal-pocket PNGs for one ion type in *pdb_path*:

    ``with_h``    -- coordination shell with H visible; water molecules shown
    ``distances`` -- coordination shell without H, with labelled M-donor distances

    Water molecules within the coordination shell are always shown alongside
    coordinating protein residues.  Returns dict with keys 'with_h' and
    'distances'; value is Path or None.
    """
    output_dir.mkdir(parents=True, exist_ok=True)
    ion_elem = _ion_element(ion_type)   # reliable element symbol, e.g. 'CA', 'ZN'
    ion_tag  = ion_elem.replace(" ", "")

    png_h = output_dir / f"{pdb_id}_metal_{ion_tag}_with_h.png"
    png_d = output_dir / f"{pdb_id}_metal_{ion_tag}_distances.png"
    results: dict[str, Path | None] = {"with_h": None, "distances": None}

    # Coordination cutoffs
    COORD_SHELL_CUTOFF = 3.5   # residue/water selection: grab full sidechains
    DONOR_CUTOFF       = 2.8   # strict donor atoms (direct bonds only)
    # Water residue names used in GROMACS-prepared structures
    WATER_RESNAMES = "HOH+SOL+WAT+TIP3+TIP"

    # Common setup block -- elem selector avoids the CA (alpha-C) / Ca2+ collision
    setup = [
        f"load {pdb_path}, s",
        "bg_color white",
        "hide everything",
        "set_color fruton_navy,    [0.157, 0.196, 0.353]",
        "set_color fruton_gold,    [0.541, 0.451, 0.059]",
        "set_color fruton_red,     [0.902, 0.000, 0.078]",
        "set_color fruton_yellow,  [0.980, 0.902, 0.275]",
        "set_color fruton_softgrey,[0.760, 0.770, 0.790]",
        "set_color fruton_water,   [0.400, 0.780, 1.000]",
        # Thin transparent cartoon for context -- fancy cartoon style
        "set cartoon_fancy_helices, 1",
        "set cartoon_smooth_loops, 1",
        "set specular, 0",
        "show cartoon, s and polymer",
        "color fruton_softgrey, s and polymer",
        "set cartoon_transparency, 0.65",
        # Metal: HETATM record with matching element
        f"select metal_at, (s and hetatm and elem {ion_elem})",
        "show sphere, metal_at",
        "color fruton_red, metal_at",
        "set sphere_scale, 0.45, metal_at",
        # Coordinating protein residues: full sidechain within shell cutoff
        f"select coord_res, byres ((s and polymer) within {COORD_SHELL_CUTOFF} of metal_at)",
        "show sticks, coord_res",
        "color fruton_gold, coord_res and elem C",
        "color atomic, coord_res and not elem C",
        # Coordinating water molecules (always include; water often IS the coord sphere)
        f"select coord_water, byres ((s and resn {WATER_RESNAMES}) within {COORD_SHELL_CUTOFF} of metal_at)",
        "show sticks, coord_water",
        "color fruton_water, coord_water and elem O",
        "color white, coord_water and elem H",
        "center metal_at",
        # Zoom to the union of metal, protein shell, and water shell
        "zoom (metal_at or coord_res or coord_water), 3",
    ]

    # ---- Scene 1: with hydrogens ----
    pml_h = setup + [
        # Protein coordination-shell H atoms
        "show sticks, coord_res and elem H",
        "color white, coord_res and elem H",
        # Water H atoms already shown via coord_water above
        "set ray_opaque_background, 1",
        "set antialias, 2",
        "set ray_shadows, 1",
        "set ray_trace_mode, 3",
        "set ray_trace_color, black",
        "set ray_trace_gain, 1.5",
        f"png {png_h}, width=900, height=900, dpi=150, ray=1",
        "quit",
    ]
    err = _run_pml(pml_h, png_h)
    if err is None:
        results["with_h"] = png_h

    # ---- Scene 2: coordination distances with H visible ----
    # H atoms are kept visible so the viewer can see which water H atoms point
    # toward the metal (short M-H, lying between metal and OW) and need
    # reorientation.  Direct donor heavy atoms get distance labels; H atoms
    # that are wrongly oriented appear close to the metal and can be compared
    # directly with the H-reorientation table in the report.
    pml_d = setup + [
        # Protein coordination-shell H atoms (already shown for water via setup)
        "show sticks, coord_res and elem H",
        "color white, coord_res and elem H",
        # Direct donor heavy atoms (N/O/S/Se) from the full structure within cutoff.
        # NOTE: 'donors' is a reserved PyMOL keyword -- use 'metal_donors'.
        (f"select metal_donors, (s within {DONOR_CUTOFF} of metal_at) and not metal_at "
         "and not elem H and (elem N or elem O or elem S or elem SE)"),
        # Ensure donors are visible as sticks so distance labels have an anchor atom
        "show sticks, metal_donors",
        # distance name, sel1, sel2, cutoff -- cutoff prevents extra cross-pairs
        f"distance coord_dist, metal_at, metal_donors, {DONOR_CUTOFF}",
        "show dashes, coord_dist",
        "enable coord_dist",
        "set dash_gap, 0.15",
        "set dash_radius, 0.07",
        "set label_size, 16",
        "set label_color, black",
        "color fruton_yellow, coord_dist",
        "set ray_opaque_background, 1",
        "set antialias, 2",
        "set ray_shadows, 1",
        "set ray_trace_mode, 3",
        "set ray_trace_color, black",
        "set ray_trace_gain, 1.5",
        f"png {png_d}, width=900, height=900, dpi=150, ray=1",
        "quit",
    ]
    err = _run_pml(pml_d, png_d)
    if err is None:
        results["distances"] = png_d

    return results
