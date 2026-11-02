"""Supplementary Figure 6 — column graphs (static image).

This figure is two column-graph screenshots stitched together; there is nothing
to compute. The image is committed as a packaged asset
(``src/ksplots/assets/sfig_06.png``) and displayed here edge-to-edge, so the
figure is reproducible from the repo via ``ksplots figure s6``.

The source screenshots are ~200 dpi; the figure is intended for a compact
(~3.3 in wide) supplementary panel, where it meets the 300 dpi print spec.
"""
from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt

_ASSET = Path(__file__).resolve().parents[2] / "assets" / "sfig_06.png"
_DPI = 300


def render(out_dir: Path) -> Path:
    img = plt.imread(str(_ASSET))
    h, w = img.shape[:2]

    # Figure sized so the raster maps 1:1 at _DPI; axes fill the whole canvas
    # (no margins/ticks) so the image is reproduced exactly.
    fig = plt.figure(figsize=(w / _DPI, h / _DPI))
    ax = fig.add_axes([0, 0, 1, 1])
    ax.imshow(img)
    ax.axis("off")

    out = Path(out_dir) / "Supplementary_Figure_06.pdf"
    fig.savefig(out, dpi=_DPI)
    plt.close(fig)
    return out
