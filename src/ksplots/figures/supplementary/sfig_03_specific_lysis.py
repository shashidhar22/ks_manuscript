"""Supplementary Figure 3 — Specific-lysis titration curves (CRA).

Ported from ``specific_lysis_curve.R``. Two stacked panels (A over B), each a
chromium-release-style peptide titration:

* **Panel A** — HIV Pol₉₈₂₋₉₉₀ (IIKDYGKQM, HLA-B*42:01; TRB CASSVDKGGTDTQYF)
* **Panel B** — HIV Vpr₃₄₋₄₂ (FPRPWLHGL, HLA-B*42:01; TRB CASSLWAGGSNEQFF)

x = peptide concentration (nM, log10); y = % specific lysis. The 0 nM
(no-peptide) control is plotted off the log scale, past a ``//`` axis break.

The assay values are inline (hardcoded from the source titrations, matching the
R script), so this module has no external data dependencies.
"""
from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

# ---------------------------------------------------------------------------
# Styling
# ---------------------------------------------------------------------------

_CURVE_COLOR = "#2563eb"   # blue
_X0 = 0.0025               # off-scale x position for the 0 nM control point
_XBREAK = 0.005            # x location of the axis-break "//" glyph
_XLIM = (0.0018, 100.0)

_AXIS_LABEL_SIZE = 18
_TICK_LABEL_SIZE = 14
_ANNOT_SIZE = 16
_PANEL_LETTER_SIZE = 34
_MARKERSIZE = 11
_LINEWIDTH = 2.4

# ---------------------------------------------------------------------------
# Per-panel data / styling (mirrors panelA / panelB in the R script)
# ---------------------------------------------------------------------------

_PANELS = [
    dict(
        letter="A",
        conc=[100.0, 10.0, 1.0, 0.1, 0.01],
        lysis=[64.3, 62.4, 38.9, 16.1, 11.8],
        ctrl=11.2,
        ylim_top=70,
        labels=[
            r"$\mathbf{HIV\ Pol_{982\text{-}990}}$",
            "IIKDYGKQM (HLA-B*42:01)",
            "TRB: CASSVDKGGTDTQYF",
        ],
    ),
    dict(
        letter="B",
        conc=[100.0, 10.0, 1.0, 0.1, 0.01],
        lysis=[85.3, 76.4, 42.8, 24.0, 19.9],
        ctrl=11.2,
        ylim_top=90,
        labels=[
            r"$\mathbf{HIV\ Vpr_{34\text{-}42}}$",
            "FPRPWLHGL (HLA-B*42:01)",
            "TRB: CASSLWAGGSNEQFF",
        ],
    ),
]


def _draw_lysis(ax, conc, lysis, ctrl, ylim_top, letter, labels) -> None:
    conc = np.asarray(conc, dtype=float)
    lysis = np.asarray(lysis, dtype=float)

    ax.set_xscale("log")
    ax.set_xlim(*_XLIM)
    ax.set_ylim(0, ylim_top)

    # Main titration curve (points + connecting line).
    ax.plot(conc, lysis, "-o", color=_CURVE_COLOR,
            markersize=_MARKERSIZE, linewidth=_LINEWIDTH, zorder=3)
    # 0 nM (no-peptide) control, off the log scale, NOT joined to the curve.
    ax.plot([_X0], [ctrl], "o", color=_CURVE_COLOR,
            markersize=_MARKERSIZE, zorder=3)

    # X ticks: labelled decades + the off-scale "0" tick.
    major = [0.01, 0.1, 1, 10, 100]
    ax.set_xticks([_X0] + major)
    ax.set_xticklabels(["0", "0.01", "0.1", "1", "10", "100"])
    minor = [v for v in np.outer(np.arange(2, 10), 10.0 ** np.arange(-3, 2)).ravel()
             if 0.01 <= v <= 100]
    ax.set_xticks(minor, minor=True)

    ax.set_yticks(range(0, ylim_top + 1, 10))
    ax.tick_params(axis="both", which="major", direction="out",
                   length=6, labelsize=_TICK_LABEL_SIZE)
    ax.tick_params(axis="x", which="minor", direction="out", length=3)

    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)

    ax.set_xlabel("Peptide concentration (nM)",
                  fontsize=_AXIS_LABEL_SIZE, fontweight="bold")
    ax.set_ylabel("% Specific lysis",
                  fontsize=_AXIS_LABEL_SIZE, fontweight="bold")

    # Big bold panel letter, out in the left margin above the y-axis.
    ax.text(-0.17, 1.02, letter, transform=ax.transAxes,
            fontsize=_PANEL_LETTER_SIZE, fontweight="bold",
            ha="left", va="bottom")

    # Upper-left annotation lines (peptide identity / sequence / TRB).
    ytop = ylim_top - 0.03 * ylim_top
    ystep = 0.085 * ylim_top
    for i, lbl in enumerate(labels):
        ax.text(0.011, ytop - i * ystep, lbl, ha="left", va="top",
                fontsize=_ANNOT_SIZE, fontweight="bold")

    # --- Axis-break "//" glyph straddling the bottom axis, left of 0.01 ---
    w = 1.12                       # multiplicative half-width (log axis)
    dy = 0.023 * ylim_top          # vertical half-height of each slash
    # white gap in the axis line
    ax.plot([_XBREAK / 1.35, _XBREAK * 1.35], [0, 0],
            color="white", linewidth=4, zorder=4, clip_on=False)
    for xc in (_XBREAK / 1.15, _XBREAK * 1.15):
        ax.plot([xc / w, xc * w], [-dy, dy],
                color="black", linewidth=1.5, zorder=5, clip_on=False)


def render(out_dir: Path) -> Path:
    # 7.5 x 11 in -> 2250 x 3300 px at 300 dpi (matches the R combined figure).
    fig, (ax_a, ax_b) = plt.subplots(2, 1, figsize=(7.5, 11))
    fig.subplots_adjust(left=0.15, right=0.95, top=0.95, bottom=0.08, hspace=0.35)

    for ax, p in zip((ax_a, ax_b), _PANELS):
        _draw_lysis(ax, p["conc"], p["lysis"], p["ctrl"],
                    p["ylim_top"], p["letter"], p["labels"])

    out = Path(out_dir) / "Supplementary_Figure_03.pdf"
    fig.savefig(out)
    plt.close(fig)
    return out
