"""Dependency checks for the FRUTON command-line runner.

This module is deliberately small and standard-library only because it runs
before the rest of the FRUTON package is imported. The runner can therefore use
it to install or report missing Python packages before modules such as
``monomer.py`` or ``pipeline_xlsx.py`` import BioPython or openpyxl. The module
also checks external command-line tools, but it does not try to install tools
that normally come from conda/pixi or system software, such as GROMACS,
AmberTools, MODELLER, or Gaussian.

The intended policy is conservative: install only small Python dependencies
that are safe to install with pip in the active interpreter, and report all
non-Python tools with actionable hints. This keeps the runner useful for end
users while avoiding hidden changes to large scientific toolchains. The return
object is structured so ``fruton.py`` can print a compact screen summary and
write a fuller dependency report into the run log.
"""

from __future__ import annotations

import importlib.util
import shutil
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Literal

DependencyKind = Literal["python", "command"]
DependencySeverity = Literal["required", "optional"]


@dataclass(frozen=True, slots=True)
class DependencyCheck:
    """Single dependency check result.

    ``name`` is the user-facing dependency label. ``kind`` distinguishes Python
    imports from external executables. ``available`` reports the final state
    after any attempted installation. ``installed`` is true only when this run
    successfully installed a missing Python dependency. ``hint`` is written for
    missing external tools or failed installations so users can repair the
    environment without reading a traceback.
    """

    name: str
    kind: DependencyKind
    severity: DependencySeverity
    available: bool
    installed: bool
    import_name: str = ""
    command_name: str = ""
    install_spec: str = ""
    hint: str = ""


@dataclass(frozen=True, slots=True)
class DependencyReport:
    """Structured dependency report returned to the runner."""

    checks: tuple[DependencyCheck, ...]
    install_missing: bool
    python_executable: str

    @property
    def ok(self) -> bool:
        return all(
            check.available or check.severity != "required" for check in self.checks
        )

    @property
    def missing_required(self) -> tuple[DependencyCheck, ...]:
        return tuple(
            check
            for check in self.checks
            if check.severity == "required" and not check.available
        )

    @property
    def missing_optional(self) -> tuple[DependencyCheck, ...]:
        return tuple(
            check
            for check in self.checks
            if check.severity == "optional" and not check.available
        )

    def to_log_lines(self) -> list[str]:
        lines = [
            f"python_executable : {self.python_executable}",
            f"install_missing   : {self.install_missing}",
            f"overall_ok        : {self.ok}",
        ]
        for check in self.checks:
            lines.append(
                " | ".join(
                    [
                        f"name={check.name}",
                        f"kind={check.kind}",
                        f"severity={check.severity}",
                        f"available={check.available}",
                        f"installed={check.installed}",
                        f"import={check.import_name}",
                        f"command={check.command_name}",
                        f"hint={check.hint}",
                    ]
                )
            )
        return lines

    def to_screen_lines(self) -> list[str]:
        installed_count = sum(1 for check in self.checks if check.installed)
        missing_required_count = len(self.missing_required)
        missing_optional_count = len(self.missing_optional)
        return [
            f"required missing : {missing_required_count}",
            f"optional missing : {missing_optional_count}",
            f"installed now    : {installed_count}",
        ]


@dataclass(frozen=True, slots=True)
class _PythonDependency:
    name: str
    import_name: str
    install_spec: str
    severity: DependencySeverity = "required"
    no_deps: bool = False


@dataclass(frozen=True, slots=True)
class _CommandDependency:
    name: str
    command_name: str
    severity: DependencySeverity = "optional"
    hint: str = ""


CORE_PYTHON_DEPENDENCIES: tuple[_PythonDependency, ...] = (
    _PythonDependency("BioPython", "Bio", "biopython", "required"),
    _PythonDependency("openpyxl", "openpyxl", "openpyxl", "required"),
    _PythonDependency("OpenMM", "openmm", "openmm", "required"),
)

OPTIONAL_PYTHON_DEPENDENCIES: tuple[_PythonDependency, ...] = (
    _PythonDependency(
        "PTM-Psi",
        "ptmpsi",
        "git+https://github.com/pnnl/PTMPSI.git",
        "optional",
        no_deps=True,
    ),
)

COMMAND_DEPENDENCIES: tuple[_CommandDependency, ...] = (
    _CommandDependency(
        "GROMACS",
        "gmx",
        "required",
        "Install GROMACS in the pixi/conda environment or provide a gmx executable in PATH.",
    ),
    _CommandDependency(
        "MODELLER Python",
        "mod10.8",
        "optional",
        "MODELLER is optional; FRUTON can still use AlphaFold fallback for large gaps.",
    ),
    _CommandDependency(
        "AmberTools antechamber",
        "antechamber",
        "optional",
        "Install AmberTools in the environment before running custom residue parameterization.",
    ),
    _CommandDependency(
        "AmberTools prepgen",
        "prepgen",
        "optional",
        "Install AmberTools in the environment before building uncapped residue templates.",
    ),
    _CommandDependency(
        "AmberTools parmchk2",
        "parmchk2",
        "optional",
        "Install AmberTools in the environment before generating frcmod files.",
    ),
    _CommandDependency(
        "AmberTools tleap",
        "tleap",
        "optional",
        "Install AmberTools in the environment before validating custom residue templates.",
    ),
    _CommandDependency(
        "Gaussian",
        "g16",
        "optional",
        "Gaussian is normally licensed software; FRUTON writes input but does not install it.",
    ),
)


def ensure_fruton_dependencies(
    *,
    install_missing: bool = True,
    include_optional_python: bool = False,
    strict: bool = False,
) -> DependencyReport:
    """Check and optionally install FRUTON runtime dependencies.

    Required Python packages are checked before the rest of the pipeline imports
    BioPython-, OpenMM-, or openpyxl-dependent modules. If ``install_missing`` is
    true, the current interpreter is used to call pip. External scientific tools
    are checked with ``shutil.which`` and are reported with hints, but this
    function does not try to install them because they normally come from pixi,
    conda, system modules, or licensed software. If ``strict`` is true, missing
    required dependencies raise ``RuntimeError`` after the report is built.
    """

    checks: list[DependencyCheck] = []

    python_dependencies = list(CORE_PYTHON_DEPENDENCIES)
    if include_optional_python:
        python_dependencies.extend(OPTIONAL_PYTHON_DEPENDENCIES)

    for dependency in python_dependencies:
        checks.append(
            _ensure_python_dependency(
                dependency=dependency,
                install_missing=install_missing,
            )
        )

    for dependency in COMMAND_DEPENDENCIES:
        checks.append(_check_command_dependency(dependency))

    report = DependencyReport(
        checks=tuple(checks),
        install_missing=install_missing,
        python_executable=sys.executable,
    )

    if strict and report.missing_required:
        missing = ", ".join(check.name for check in report.missing_required)
        hints = "; ".join(
            check.hint for check in report.missing_required if check.hint
        )
        raise RuntimeError(f"Missing required FRUTON dependencies: {missing}. {hints}")

    return report


def _ensure_python_dependency(
    *,
    dependency: _PythonDependency,
    install_missing: bool,
) -> DependencyCheck:
    if _python_import_available(dependency.import_name):
        return DependencyCheck(
            name=dependency.name,
            kind="python",
            severity=dependency.severity,
            available=True,
            installed=False,
            import_name=dependency.import_name,
            install_spec=dependency.install_spec,
        )

    if not install_missing:
        return DependencyCheck(
            name=dependency.name,
            kind="python",
            severity=dependency.severity,
            available=False,
            installed=False,
            import_name=dependency.import_name,
            install_spec=dependency.install_spec,
            hint=f"Install with: {sys.executable} -m pip install {dependency.install_spec}",
        )

    try:
        _ensure_pip_available()
        command = [sys.executable, "-m", "pip", "install"]
        if dependency.no_deps:
            command.append("--no-deps")
        command.append(dependency.install_spec)
        subprocess.run(command, check=True)
    except Exception as exc:
        return DependencyCheck(
            name=dependency.name,
            kind="python",
            severity=dependency.severity,
            available=False,
            installed=False,
            import_name=dependency.import_name,
            install_spec=dependency.install_spec,
            hint=(
                f"Automatic install failed: {type(exc).__name__}: {exc}. "
                f"Try: pixi add pip && pixi run python -m pip install {dependency.install_spec}"
            ),
        )

    available_after_install = _python_import_available(dependency.import_name)
    return DependencyCheck(
        name=dependency.name,
        kind="python",
        severity=dependency.severity,
        available=available_after_install,
        installed=available_after_install,
        import_name=dependency.import_name,
        install_spec=dependency.install_spec,
        hint=(
            ""
            if available_after_install
            else f"pip completed but import still failed for {dependency.import_name!r}."
        ),
    )


def _python_import_available(import_name: str) -> bool:
    return importlib.util.find_spec(import_name) is not None


def _ensure_pip_available() -> None:
    pip_check = subprocess.run(
        [sys.executable, "-m", "pip", "--version"],
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
        check=False,
    )
    if pip_check.returncode == 0:
        return

    ensurepip_check = subprocess.run(
        [sys.executable, "-m", "ensurepip", "--upgrade"],
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
        check=False,
    )
    if ensurepip_check.returncode != 0:
        raise RuntimeError(
            "pip is not available in this Python environment and ensurepip failed. "
            "In pixi, run: pixi add pip"
        )


def _check_command_dependency(dependency: _CommandDependency) -> DependencyCheck:
    path = shutil.which(dependency.command_name)
    return DependencyCheck(
        name=dependency.name,
        kind="command",
        severity=dependency.severity,
        available=path is not None,
        installed=False,
        command_name=dependency.command_name,
        hint="" if path is not None else dependency.hint,
    )


if __name__ == "__main__":
    report = ensure_fruton_dependencies(
        install_missing=False,
        include_optional_python=True,
        strict=False,
    )
    for line in report.to_log_lines():
        print(line)
