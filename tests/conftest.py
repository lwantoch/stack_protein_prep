"""
conftest.py

Pytest configuration for this project.

Purpose
-------
Add the src directory to sys.path so that pytest can import the package
during test collection.
"""

import sys
from pathlib import Path

TESTS_DIR = Path(__file__).resolve().parent
PROJECT_ROOT_DIR = TESTS_DIR.parent
SRC_DIR = PROJECT_ROOT_DIR / "src"

if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))
