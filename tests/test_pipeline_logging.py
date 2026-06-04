"""Tests for pipeline_logging utilities."""
from __future__ import annotations

import re
from pathlib import Path

import pytest

from stack_protein_preparation.pipeline_logging import (
    ModuleLogPaths,
    append_key_value_block,
    append_log_text,
    build_module_log_paths,
    current_timestamp,
    current_timestamp_for_filename,
)


# ---------------------------------------------------------------------------
# Timestamp helpers
# ---------------------------------------------------------------------------

def test_current_timestamp_format() -> None:
    ts = current_timestamp()
    assert re.fullmatch(r"\d{4}-\d{2}-\d{2} \d{2}:\d{2}:\d{2}", ts)


def test_current_timestamp_for_filename_format() -> None:
    ts = current_timestamp_for_filename()
    assert re.fullmatch(r"\d{8}_\d{6}", ts)


# ---------------------------------------------------------------------------
# build_module_log_paths
# ---------------------------------------------------------------------------

def test_log_path_without_variant(tmp_path: Path) -> None:
    paths = build_module_log_paths(
        protein_data_dir=tmp_path,
        pdb_id="1ABC",
        module_name="protonation",
    )
    assert paths.module_log_path.name == "protonation.log"
    assert paths.protein_log_dir.name == "1ABC"
    assert paths.logs_root_dir == tmp_path / "logs"


def test_log_path_with_variant(tmp_path: Path) -> None:
    paths = build_module_log_paths(
        protein_data_dir=tmp_path,
        pdb_id="1ABC",
        module_name="protonation",
        variant_label="gaps",
    )
    assert paths.module_log_path.name == "protonation.gaps.log"


def test_pdb_id_is_uppercased(tmp_path: Path) -> None:
    paths = build_module_log_paths(
        protein_data_dir=tmp_path,
        pdb_id="1abc",
        module_name="filler",
    )
    assert paths.protein_log_dir.name == "1ABC"


def test_protein_log_dir_created(tmp_path: Path) -> None:
    paths = build_module_log_paths(
        protein_data_dir=tmp_path,
        pdb_id="2XYZ",
        module_name="gaps",
    )
    assert paths.protein_log_dir.is_dir()


def test_empty_pdb_id_raises(tmp_path: Path) -> None:
    with pytest.raises(ValueError, match="pdb_id"):
        build_module_log_paths(
            protein_data_dir=tmp_path,
            pdb_id="",
            module_name="gaps",
        )


def test_empty_module_name_raises(tmp_path: Path) -> None:
    with pytest.raises(ValueError, match="module_name"):
        build_module_log_paths(
            protein_data_dir=tmp_path,
            pdb_id="1ABC",
            module_name="",
        )


def test_returns_module_log_paths_dataclass(tmp_path: Path) -> None:
    result = build_module_log_paths(
        protein_data_dir=tmp_path,
        pdb_id="1ABC",
        module_name="cap",
    )
    assert isinstance(result, ModuleLogPaths)


# ---------------------------------------------------------------------------
# append_log_text
# ---------------------------------------------------------------------------

def test_append_creates_file(tmp_path: Path) -> None:
    log = tmp_path / "sub" / "test.log"
    append_log_text(log, "hello")
    assert log.is_file()
    assert "hello" in log.read_text()


def test_append_adds_newline_if_missing(tmp_path: Path) -> None:
    log = tmp_path / "test.log"
    append_log_text(log, "no newline")
    content = log.read_text()
    assert content.endswith("\n")


def test_append_multiple_calls_accumulate(tmp_path: Path) -> None:
    log = tmp_path / "test.log"
    append_log_text(log, "line1")
    append_log_text(log, "line2")
    content = log.read_text()
    assert "line1" in content
    assert "line2" in content


def test_append_empty_string_writes_nothing(tmp_path: Path) -> None:
    log = tmp_path / "test.log"
    append_log_text(log, "")
    assert not log.exists()


# ---------------------------------------------------------------------------
# append_key_value_block
# ---------------------------------------------------------------------------

def test_key_value_block_written(tmp_path: Path) -> None:
    log = tmp_path / "test.log"
    append_key_value_block(log, "STEP", [("key1", "val1"), ("key2", "val2")])
    content = log.read_text()
    assert "STEP" in content
    assert "key1: val1" in content
    assert "key2: val2" in content


def test_key_value_block_ends_with_blank_line(tmp_path: Path) -> None:
    log = tmp_path / "test.log"
    append_key_value_block(log, "TITLE", [("a", "b")])
    content = log.read_text()
    assert content.endswith("\n\n") or content.endswith("\n")
