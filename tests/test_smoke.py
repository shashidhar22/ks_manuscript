"""Smoke test: every figure module renders to a non-empty PDF."""
from __future__ import annotations

from pathlib import Path

import pytest

from ksplots import figures


@pytest.mark.parametrize("name", figures.available())
def test_render(tmp_path: Path, name: str) -> None:
    out = figures.render(name, tmp_path)
    assert out.exists(), f"{name} did not produce a file"
    assert out.stat().st_size > 0, f"{name} produced an empty file"
