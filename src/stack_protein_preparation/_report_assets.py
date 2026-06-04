"""Logo and file-asset discovery for FRUTON report modules.

Both per-protein and summary reports need to locate IDIS/FRUTON logo images
at render time without knowing the caller's working directory or installation
layout. This module provides a single unified finder used by all report modules.
"""
from __future__ import annotations

import os
from pathlib import Path

from ._report_palette import _IDIS_LOGO_NAMES, _FRUTON_LOGO_NAMES


def _asset_search_roots(base_dir: Path | None) -> list[Path]:
    """Return candidate directories to search for report assets.

    Checked in order:
      1. base_dir and its parent (protein dir or run dir, when provided)
      2. Repo-root-anchored data/assets — canonical install location
      3. cwd fallbacks for notebook/interactive runs
      4. Module-adjacent paths
    """
    module_dir = Path(__file__).resolve().parent
    repo_root = module_dir.parent.parent
    roots: list[Path] = []

    if base_dir is not None:
        roots += [
            base_dir / "report",
            base_dir,
            base_dir.parent,
            base_dir.parent / "assets",
            base_dir.parent.parent / "assets",
            base_dir.parent.parent / "data",
        ]

    roots += [
        repo_root / "data" / "assets",
        repo_root / "data",
        repo_root / "assets",
        repo_root,
        Path.cwd() / "data" / "assets",
        Path.cwd() / "data",
        Path.cwd(),
        module_dir / "assets",
        module_dir.parent / "assets",
        Path("/mnt/data"),
    ]

    seen: set[Path] = set()
    unique: list[Path] = []
    for root in roots:
        try:
            resolved = root.resolve()
        except OSError:
            resolved = root
        if resolved not in seen:
            seen.add(resolved)
            unique.append(root)
    return unique


def find_asset(
    env_var: str,
    names: tuple[str, ...],
    base_dir: Path | None = None,
) -> Path | None:
    """Return the path to a named image asset, or None if not found.

    Checks *env_var* first so the caller can override via the environment,
    then searches standard roots derived from *base_dir*.
    """
    env_value = os.environ.get(env_var, "").strip()
    if env_value:
        p = Path(env_value).expanduser()
        if p.exists():
            return p
    for root in _asset_search_roots(base_dir):
        for name in names:
            candidate = root / name
            if candidate.exists():
                return candidate
    return None


def find_idis_logo(base_dir: Path | None = None) -> Path | None:
    return find_asset("IDIS_LOGO_PATH", _IDIS_LOGO_NAMES, base_dir)


def find_fruton_logo(base_dir: Path | None = None) -> Path | None:
    return find_asset("FRUTON_LOGO_PATH", _FRUTON_LOGO_NAMES, base_dir)
