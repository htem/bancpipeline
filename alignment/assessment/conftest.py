"""Shared pytest fixtures for alignment/assessment/.

The main aligner is `banc-alignment-run.py` — a script filename with a
hyphen, which Python cannot import directly. Load it as a module once via
importlib and expose helper accessors.
"""
from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import pytest

_ALIGN_DIR = Path(__file__).resolve().parent.parent
_ALIGN_PY = _ALIGN_DIR / "banc-alignment-run.py"


def _load_align_module():
    if "_banc_align" in sys.modules:
        return sys.modules["_banc_align"]
    if str(_ALIGN_DIR) not in sys.path:
        sys.path.insert(0, str(_ALIGN_DIR))
    spec = importlib.util.spec_from_file_location("_banc_align", _ALIGN_PY)
    mod = importlib.util.module_from_spec(spec)
    sys.modules["_banc_align"] = mod
    spec.loader.exec_module(mod)
    return mod


@pytest.fixture(scope="session")
def align():
    """Aligner module loaded from banc-alignment-run.py."""
    return _load_align_module()


@pytest.fixture(scope="session")
def splits():
    """alignment_splits module (importable by name)."""
    if str(_ALIGN_DIR) not in sys.path:
        sys.path.insert(0, str(_ALIGN_DIR))
    import alignment_splits  # noqa: E402
    return alignment_splits
