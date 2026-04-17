"""Figure registry.

Each figure module exposes ``render(out_dir: Path) -> Path`` and is registered
here so the CLI can dispatch by name/number.
"""
from __future__ import annotations

import importlib
from pathlib import Path
from typing import Callable, Dict

# (number, kind) -> module path. kind is "main" or "supp".
_REGISTRY: Dict[str, str] = {
    "1": "ksplots.figures.figure_01",
    "2": "ksplots.figures.figure_02",
    "3": "ksplots.figures.figure_03",
    "s1": "ksplots.figures.supplementary.sfig_01",
    "s2": "ksplots.figures.supplementary.sfig_02",
    "4": "ksplots.figures.figure_04_wetlab",
    "5": "ksplots.figures.figure_05_wetlab",
    "6": "ksplots.figures.figure_06_pel",
    "s3": "ksplots.figures.supplementary.sfig_03_alluvial",
    "s4": "ksplots.figures.supplementary.sfig_04_cos7",
}


def available() -> list[str]:
    return list(_REGISTRY.keys())


def render(name: str, out_dir: Path) -> Path:
    name = str(name).lower()
    if name not in _REGISTRY:
        raise KeyError(f"Unknown figure {name!r}. Available: {available()}")
    mod = importlib.import_module(_REGISTRY[name])
    out_dir.mkdir(parents=True, exist_ok=True)
    fn: Callable[[Path], Path] = mod.render
    return fn(out_dir)
