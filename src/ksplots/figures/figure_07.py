"""Figure 7 — schematic (static image).

This figure is a single BioRender schematic with two sections (A and B); there
is nothing to compute. The image is committed as a packaged asset
(``src/ksplots/assets/fig_07.png``) and displayed here edge-to-edge, so the
figure is reproducible from the repo via ``ksplots figure 7``.

The source is a 600 dpi BioRender PNG export with a transparent background; it
is flattened onto white here so the output has no alpha (matching the PLOS
flattened-RGB print spec).
"""
from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt

_ASSET = Path(__file__).resolve().parents[1] / "assets" / "fig_07.png"
_DPI = 600


def render(out_dir: Path) -> Path:
    img = plt.imread(str(_ASSET))

    # Flatten any alpha channel onto a white background -> HxWx3 RGB.
    if img.ndim == 3 and img.shape[2] == 4:
        rgb, alpha = img[..., :3], img[..., 3:4]
        img = rgb * alpha + (1.0 - alpha)

    h, w = img.shape[:2]

    # Figure sized so the raster maps 1:1 at _DPI; axes fill the whole canvas
    # (no margins/ticks) so the image is reproduced exactly.
    fig = plt.figure(figsize=(w / _DPI, h / _DPI))
    ax = fig.add_axes([0, 0, 1, 1])
    ax.imshow(img)
    ax.axis("off")

    out = Path(out_dir) / "Figure_07.pdf"
    fig.savefig(out, dpi=_DPI)
    plt.close(fig)
    return out
