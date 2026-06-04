"""Module-level injectable logger and settings for the pipeline runner."""
from __future__ import annotations

from typing import Callable

DEFAULT_FORCE_FIELD = "amber99sb-ildn"
DEFAULT_WATER_MODEL = "tip3p"
DEFAULT_PH = 7.0

_log: Callable[[str, list[str]], None] = lambda _tag, _lines: None
_force_field: str = DEFAULT_FORCE_FIELD
_water_model: str = DEFAULT_WATER_MODEL
_ph: float = DEFAULT_PH


def init_runner(
    log_fn: Callable[[str, list[str]], None],
    force_field: str = DEFAULT_FORCE_FIELD,
    water_model: str = DEFAULT_WATER_MODEL,
    ph: float = DEFAULT_PH,
) -> None:
    global _log, _force_field, _water_model, _ph
    _log = log_fn
    _force_field = force_field
    _water_model = water_model
    _ph = ph
