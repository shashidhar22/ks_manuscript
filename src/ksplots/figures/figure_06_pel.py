"""Figure 6 — PEL cell functional assays.

Three bar-chart panels arranged horizontally (A, B, C):

* **Panel A** — mNeonGreen+ (%) bar chart: TCR 2 & 3 ± TPA.
* **Panel B** — IFN-γ (pg/mL) bar chart: TCR 2 & 3 ± TPA.
* **Panel C** — Specific Lysis (%) bar chart: TCR 2 & 3 ± TPA.

Data sources: ``data/wetlab/data/FIgure6_PanelA.csv`` and
``data/wetlab/data/Figure6_Panel{B,C}.csv``
"""
from __future__ import annotations

import csv
import io
import re
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.gridspec import GridSpec

from ..config import DATA_DIR

_DATA_DIR = DATA_DIR / "wetlab" / "data"

# ---------------------------------------------------------------------------
# Colors
# ---------------------------------------------------------------------------

_TCR_RENAME = {
    "TCR 2": "ORF6329-337-Specific TCR\nCASSIAGHEQFF",
    "TCR 3": "ORF6329-337-Specific TCR\nCASSIAGHEQYF",
}

_TCR_COLORS = {
    "ORF6329-337-Specific TCR\nCASSIAGHEQFF": "#E41A1C",  # red
    "ORF6329-337-Specific TCR\nCASSIAGHEQYF": "#377EB8",  # blue
}
_NEG_COLOR = "#999999"   # grey for (-) TPA

# Panel letters (A/B/C) drawn just outside the upper-left of each panel's
# y-axis. _PANEL_LETTER_X is the single nudge knob (axes-fraction x, more
# negative = further out past the y-axis); _PANEL_LETTER_SIZE is the font size.
_PANEL_LETTER_X = -0.08
_PANEL_LETTER_SIZE = 54  # slightly enlarged per request (was 50, Fig 4-matched)

# Axis text sizes, shared across all three panels.
_TICK_LABELSIZE = 26      # x & y numeric/category tick labels (was 22)
_YAXIS_TITLE_SIZE = 32    # rotated y-axis title, e.g. "mNeonGreen+ (%)" (was 28)

# Legend styling. Legends sit BELOW each panel (loc="upper center"), centered
# under the plot so they never extend past the plotting width. _LEGEND_Y is the
# anchor y in axes-fraction (negative = below the x-axis, clearing the tick and
# x-axis labels; more negative = further down).
_LEGEND_FONTSIZE = 24
_LEGEND_Y = -0.06

# Significance asterisks drawn over the bar-panel comparison brackets.
_SIG_FONTSIZE = 30

# ---------------------------------------------------------------------------
# Gene name splitting (reuse logic from figure_05)
# ---------------------------------------------------------------------------

_KNOWN_GENES = ["ORF57", "ORF59", "ORF6"]


def _split_gene_range(s: str) -> tuple[str, str, str] | None:
    for gene in _KNOWN_GENES:
        idx = s.find(gene)
        if idx >= 0:
            after = s[idx + len(gene):]
            m = re.match(r'(\d+[-–]\d+)(.*)', after)
            if m:
                return s[:idx + len(gene)], m.group(1), m.group(2)
    m = re.match(r'(.+?)(\d+[-–]\d+)(.*)', s)
    if m:
        return m.group(1), m.group(2), m.group(3)
    return None


def _subscript(s: str) -> str:
    parts = _split_gene_range(s)
    if parts:
        pre, rng, post = parts
        return f"{pre}$_{{{rng}}}${post}"
    return s


# ---------------------------------------------------------------------------
# Bar chart loader / drawer (Panels A, B & C)
# ---------------------------------------------------------------------------


def _load_bar_csv(path: Path):
    """Load Panel A/B/C CSV: returns (conditions, {tcr_name: (mean, sem)})."""
    text = path.read_text().strip()
    reader = csv.reader(io.StringIO(text))
    rows = list(reader)

    header = rows[0]
    # Group names from header (skip empty first col)
    raw_names = [h.strip() for h in header[1:]]
    group_names = []
    for n in raw_names:
        if not group_names or n != group_names[-1]:
            group_names.append(n)
    n_reps = len(raw_names) // len(group_names)

    conditions = []
    group_data: dict[str, list[list[float]]] = {n: [] for n in group_names}

    for row in rows[1:]:
        conditions.append(row[0].strip())
        idx = 1
        for gname in group_names:
            reps = [float(row[idx + j]) for j in range(n_reps)]
            group_data[gname].append(reps)
            idx += n_reps

    result = {}
    for gname in group_names:
        arr = np.array(group_data[gname])
        display = _TCR_RENAME.get(gname, gname)
        result[display] = (arr.mean(axis=1), arr.std(axis=1, ddof=1) / np.sqrt(arr.shape[1]))

    return conditions, result


def _draw_bar_panel(
    ax: plt.Axes,
    conditions: list[str],
    tcr_data: dict[str, tuple[np.ndarray, np.ndarray]],
    ylabel: str,
    panel_letter: str,
) -> None:
    """Draw grouped bar chart with significance testing."""
    n_cond = len(conditions)
    n_tcrs = len(tcr_data)
    bar_w = 0.35
    x = np.arange(n_cond)

    tcr_names = list(tcr_data.keys())
    for i, tname in enumerate(tcr_names):
        mean, sem = tcr_data[tname]
        color = _TCR_COLORS.get(tname, list(_TCR_COLORS.values())[i])
        offset = (i - n_tcrs / 2 + 0.5) * bar_w
        disp_label = "\n".join(_subscript(ln) for ln in tname.split("\n"))
        ax.bar(x + offset, mean, bar_w * 0.85, yerr=sem, capsize=4,
               color=color, edgecolor="black", linewidth=0.5,
               label=disp_label, error_kw={"linewidth": 0.8})

    # Significance: compare first condition vs each subsequent condition
    # within each TCR. Draw bracket between condition pairs.
    bracket_idx = 0
    for i, tname in enumerate(tcr_names):
        mean, sem = tcr_data[tname]
        offset = (i - n_tcrs / 2 + 0.5) * bar_w
        for j in range(1, n_cond):
            if mean[j] > mean[0] * 1.5:
                y_max = max(mean[0] + sem[0], mean[j] + sem[j])
                y_bar = y_max * 1.08 + bracket_idx * y_max * 0.10
                x_left = 0 + offset
                x_right = j + offset
                ax.plot([x_left, x_left, x_right, x_right],
                        [y_bar, y_bar * 1.02, y_bar * 1.02, y_bar],
                        color="black", lw=1.0, clip_on=False)
                ax.text((x_left + x_right) / 2, y_bar * 1.03, "***",
                        ha="center", va="bottom", fontsize=_SIG_FONTSIZE,
                        fontweight="bold")
                bracket_idx += 1

    ax.set_xticks(x)
    ax.set_xticklabels(conditions, fontsize=_TICK_LABELSIZE, fontweight="bold")
    ax.set_ylabel(ylabel, fontsize=_YAXIS_TITLE_SIZE, fontweight="bold")
    ax.tick_params(labelsize=_TICK_LABELSIZE)
    ax.text(_PANEL_LETTER_X, 1.0, panel_letter, transform=ax.transAxes,
            fontweight="bold", fontsize=_PANEL_LETTER_SIZE,
            ha="right", va="bottom")
    ax.legend(frameon=False, loc="upper center",
              bbox_to_anchor=(0.5, _LEGEND_Y), ncol=2,
              columnspacing=1.5, handletextpad=0.5,
              prop={"weight": "bold", "size": _LEGEND_FONTSIZE})
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)


# ---------------------------------------------------------------------------
# Render
# ---------------------------------------------------------------------------


def render(out_dir: Path) -> Path:
    fig = plt.figure(figsize=(36, 13))
    gs = GridSpec(1, 3, wspace=0.30,
                  left=0.06, right=0.97, top=0.93, bottom=0.15)

    # Panel A — mNeonGreen+ bar chart
    ax_a = fig.add_subplot(gs[0, 0])
    cond_a, data_a = _load_bar_csv(_DATA_DIR / "FIgure6_PanelA.csv")
    _draw_bar_panel(ax_a, cond_a, data_a,
                    ylabel="mNeonGreen+ (%)", panel_letter="A")

    # Panel B — IFN-γ bar chart
    ax_b = fig.add_subplot(gs[0, 1])
    cond_b, data_b = _load_bar_csv(_DATA_DIR / "Figure6_PanelB.csv")
    _draw_bar_panel(ax_b, cond_b, data_b,
                    ylabel="IFN-γ (pg/mL)", panel_letter="B")

    # Panel C — Specific Lysis bar chart
    ax_c = fig.add_subplot(gs[0, 2])
    cond_c, data_c = _load_bar_csv(_DATA_DIR / "Figure6_PanelC.csv")
    _draw_bar_panel(ax_c, cond_c, data_c,
                    ylabel="Specific Lysis (%)", panel_letter="C")

    out = Path(out_dir) / "Figure_06.pdf"
    fig.savefig(out)
    plt.close(fig)
    return out
