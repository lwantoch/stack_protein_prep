"""Dependency checks for the FRUTON command-line runner.

Standard-library only — runs before the rest of FRUTON is imported.

Policy
------
All dependencies are required.  Python packages are installed automatically
with pip when missing.  External tools (GROMACS, AmberTools, Modeller,
Gaussian) are checked with shutil.which; when missing the user receives a
specific actionable message instead of a silent skip.
"""

from __future__ import annotations

import importlib.util
import shutil
import subprocess
import sys
from dataclasses import dataclass


@dataclass(frozen=True, slots=True)
class DependencyCheck:
    name: str
    kind: str          # "python" | "command"
    available: bool
    installed: bool    # True only when this run installed it via pip
    import_name: str = ""
    command_name: str = ""
    install_spec: str = ""
    hint: str = ""     # shown when not available


@dataclass(frozen=True, slots=True)
class DependencyReport:
    checks: tuple[DependencyCheck, ...]
    python_executable: str

    @property
    def ok(self) -> bool:
        return all(check.available for check in self.checks)

    @property
    def missing(self) -> tuple[DependencyCheck, ...]:
        return tuple(check for check in self.checks if not check.available)

    def to_log_lines(self) -> list[str]:
        lines = [
            f"python_executable : {self.python_executable}",
            f"overall_ok        : {self.ok}",
        ]
        for check in self.checks:
            lines.append(
                " | ".join([
                    f"name={check.name}",
                    f"kind={check.kind}",
                    f"available={check.available}",
                    f"installed={check.installed}",
                    f"hint={check.hint}",
                ])
            )
        return lines


@dataclass(frozen=True, slots=True)
class _PythonDependency:
    name: str
    import_name: str
    install_spec: str
    no_deps: bool = False


@dataclass(frozen=True, slots=True)
class _CommandDependency:
    name: str
    command_name: str
    hint: str   # shown when the command is missing


PYTHON_DEPENDENCIES: tuple[_PythonDependency, ...] = (
    _PythonDependency("BioPython",  "Bio",      "biopython"),
    _PythonDependency("openpyxl",   "openpyxl", "openpyxl"),
    _PythonDependency("OpenMM",     "openmm",   "openmm"),
)

COMMAND_DEPENDENCIES: tuple[_CommandDependency, ...] = (
    _CommandDependency(
        "GROMACS",
        "gmx",
        "GROMACS (gmx) is not in PATH.\n"
        "  On CESGA: module load gromacs/2025.3\n"
        "  Add to your fruton wrapper or ~/.bashrc.",
    ),
    _CommandDependency(
        "MODELLER",
        "mod10.8",
        "MODELLER is not in PATH.  Fix: pixi install\n"
        "  Make sure KEY_MODELLER is set: export KEY_MODELLER=<your_key>\n"
        "  Register free academic key at: https://salilab.org/modeller/registration.html",
    ),
    _CommandDependency(
        "MAFFT",
        "mafft",
        "mafft is not in PATH.\n"
        "  On CESGA: module load mafft/7.525-with-extensions\n"
        "  Add to your fruton wrapper or ~/.bashrc.",
    ),
    _CommandDependency(
        "AmberTools — antechamber",
        "antechamber",
        "antechamber is not in PATH.\n"
        "  On CESGA: module load amber/20.13-AmberTools-22.2\n"
        "  (may require prerequisites — check: module spider amber/20.13-AmberTools-22.2)",
    ),
    _CommandDependency(
        "AmberTools — prepgen",
        "prepgen",
        "prepgen is not in PATH.\n"
        "  On CESGA: module load amber/20.13-AmberTools-22.2",
    ),
    _CommandDependency(
        "AmberTools — parmchk2",
        "parmchk2",
        "parmchk2 is not in PATH.\n"
        "  On CESGA: module load amber/20.13-AmberTools-22.2",
    ),
    _CommandDependency(
        "AmberTools — tleap",
        "tleap",
        "tleap is not in PATH.\n"
        "  On CESGA: module load amber/20.13-AmberTools-22.2",
    ),
    _CommandDependency(
        "AmberTools — MCPB.py",
        "MCPB.py",
        "MCPB.py is not in PATH.\n"
        "  On CESGA: module load amber/20.13-AmberTools-22.2",
    ),
    _CommandDependency(
        "Gaussian 16",
        "g16",
        "g16 is not in PATH.\n"
        "  On CESGA: module load cesga/2020 g16/c1\n"
        "  Gaussian is licensed software — load on the HPC cluster before running.",
    ),
)


def ensure_fruton_dependencies() -> DependencyReport:
    """Check all FRUTON dependencies.  Auto-install missing Python packages."""
    checks: list[DependencyCheck] = []

    for dep in PYTHON_DEPENDENCIES:
        checks.append(_ensure_python_dependency(dep))

    for dep in COMMAND_DEPENDENCIES:
        checks.append(_check_command_dependency(dep))

    return DependencyReport(
        checks=tuple(checks),
        python_executable=sys.executable,
    )


def _ensure_python_dependency(dep: _PythonDependency) -> DependencyCheck:
    if _importable(dep.import_name):
        return DependencyCheck(
            name=dep.name, kind="python",
            available=True, installed=False,
            import_name=dep.import_name, install_spec=dep.install_spec,
        )

    # Attempt auto-install
    try:
        _ensure_pip()
        cmd = [sys.executable, "-m", "pip", "install"]
        if dep.no_deps:
            cmd.append("--no-deps")
        cmd.append(dep.install_spec)
        subprocess.run(cmd, check=True, capture_output=True)
    except Exception as exc:
        return DependencyCheck(
            name=dep.name, kind="python",
            available=False, installed=False,
            import_name=dep.import_name, install_spec=dep.install_spec,
            hint=(
                f"Automatic install failed ({type(exc).__name__}: {exc}).\n"
                f"  Fix: pixi run python -m pip install {dep.install_spec}"
            ),
        )

    ok = _importable(dep.import_name)
    return DependencyCheck(
        name=dep.name, kind="python",
        available=ok, installed=ok,
        import_name=dep.import_name, install_spec=dep.install_spec,
        hint="" if ok else f"pip finished but import of {dep.import_name!r} still fails.",
    )


def _check_command_dependency(dep: _CommandDependency) -> DependencyCheck:
    found = shutil.which(dep.command_name) is not None
    return DependencyCheck(
        name=dep.name, kind="command",
        available=found, installed=False,
        command_name=dep.command_name,
        hint="" if found else dep.hint,
    )


def _importable(import_name: str) -> bool:
    return importlib.util.find_spec(import_name) is not None


def _ensure_pip() -> None:
    if subprocess.run(
        [sys.executable, "-m", "pip", "--version"],
        stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
    ).returncode == 0:
        return
    if subprocess.run(
        [sys.executable, "-m", "ensurepip", "--upgrade"],
        stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
    ).returncode != 0:
        raise RuntimeError(
            "pip is not available and ensurepip failed.  Fix: pixi add pip"
        )


if __name__ == "__main__":
    report = ensure_fruton_dependencies()
    for line in report.to_log_lines():
        print(line)
