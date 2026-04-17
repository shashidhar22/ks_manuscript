"""Figure 4 (wetlab) — HIV-specific TCR functional validation.

Reproduces FIGURES_IT.pptx Slide 1:

* **Panel A** — HLA restriction: mNeonGreen+ (%) for 5 TCRs with
  HLA-B*45:01 (off-target, near-zero) vs HLA-B*42:01 (on-target, ~60-75%)
  at (+) 100 nM peptide.
* **Panel C** — CRA dose-response: Specific Lysis (%) vs peptide
  concentration (nM) for 5 HIV-specific TCRs with sigmoidal fits and
  half-max EC50 annotations.

Data sources:
  - Panel A: Master Prism (off-target bars + on-target from titration @ 100 nM)
  - Panel C: CRA Prism (``PUBLICATION_CRA_All 5_FINAL_14Jun2024``)
"""
from __future__ import annotations

from pathlib import Path

import re

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.gridspec import GridSpec
from scipy.optimize import curve_fit

import pandas as pd

from ..config import DATA_DIR


def _subscript_gene_range(s: str) -> str:
    """Convert gene number ranges to subscript, e.g. Nef71-79 → Nef$_{71-79}$."""
    return re.sub(r'(\d+[\-–]\d+)', r'$_{\1}$', s)


# TCR display order (matching PPTX legend)
_TCR_ORDER = [
    "Nef-Specific TCR 1",
    "Nef-Specific TCR 2",
    "Vpr-Specific TCR 1",
    "Vpr-Specific TCR 2",
    "Pol-Specific TCR",
]

# CRA column name → display name (for dose-response Panel C)
_CRA_COLS = [
    ("Nef71-79-Specific TCR", "Nef-Specific TCR 1"),
    ("Nef71-79-Specific TCR (2)", "Nef-Specific TCR 2"),
    ("Vpr34-42-Specific TCR (2)", "Vpr-Specific TCR 1"),
    ("Vpr34-42-Specific TCR", "Vpr-Specific TCR 2"),
    ("Pol982-990-Specific TCR", "Pol-Specific TCR"),
]

_MARKERS = {
    "Nef-Specific TCR 1": "s",
    "Nef-Specific TCR 2": "o",
    "Vpr-Specific TCR 1": "^",
    "Vpr-Specific TCR 2": "v",
    "Pol-Specific TCR": "D",
}


# ---------------------------------------------------------------------------
# Panel A — HLA restriction bars
# ---------------------------------------------------------------------------

_PANEL_A_CSV = DATA_DIR / "wetlab" / "data" / "Figure4_PanelA.csv"

# TCR column headers in CSV → display names
_PANEL_A_TCR_MAP = {
    "Nef71-79": "Nef-Specific TCR",
    "Vpr34-42": "Vpr-Specific TCR",
}

# 8 distinct TCR colors for use across figures
TCR_PALETTE = {
    "TCR1": "#E41A1C",  # red
    "TCR2": "#377EB8",  # blue
    "TCR3": "#4DAF4A",  # green
    "TCR4": "#FF7F00",  # orange
    "TCR5": "#984EA3",  # purple
    "TCR6": "#A65628",  # brown
    "TCR7": "#F781BF",  # pink
    "TCR8": "#999999",  # grey
}

_PANEL_A_COLORS = [
    TCR_PALETTE["TCR1"],  # Nef71-79 TCR 1
    TCR_PALETTE["TCR2"],  # Nef71-79 TCR 2
    TCR_PALETTE["TCR3"],  # Vpr34-42 TCR
]


def _load_panel_a_data() -> tuple[list[str], list[str], np.ndarray]:
    """Parse the Panel A CSV.

    Returns (conditions, tcr_labels, values) where values is
    shape (n_conditions, n_tcrs).
    """
    lines = _PANEL_A_CSV.read_text().strip().replace("\r\n", "\n")
    # The CSV uses embedded newlines within quoted fields; split on unquoted
    # newlines by re-joining then parsing.
    import csv
    import io
    reader = csv.reader(io.StringIO(lines))
    rows = list(reader)

    # Row 0 is the header: first cell empty, then TCR column names
    tcr_labels = [c.strip().replace("\n", "\n") for c in rows[0][1:]]
    conditions = []
    values = []
    for row in rows[1:]:
        conditions.append(row[0].strip().replace("\n", "\n"))
        values.append([float(v) for v in row[1:]])

    return conditions, tcr_labels, np.array(values)


def _draw_panel_a(ax) -> None:
    """Horizontal grouped bar chart from Panel A CSV data."""
    conditions, tcr_labels, values = _load_panel_a_data()

    n_cond = len(conditions)
    n_tcrs = len(tcr_labels)
    bar_w = 0.22
    y = np.arange(n_cond)

    for i in range(n_tcrs):
        color = _PANEL_A_COLORS[i]
        offset = (i - n_tcrs / 2 + 0.5) * bar_w
        label = "\n".join(_subscript_gene_range(ln) for ln in tcr_labels[i].split("\n"))
        ax.barh(y + offset, values[:, i], bar_w * 0.85, color=color,
                edgecolor="black", linewidth=0.5, label=label)

    ax.set_yticks(y)
    ax.set_yticklabels(conditions, fontsize=16, fontweight="bold")
    ax.set_xlabel("mNeonGreen+ (%)", fontsize=20, fontweight="bold")
    ax.set_xlim(0, max(values.max() * 1.1, 80))
    ax.tick_params(axis="x", labelsize=18)
    leg = ax.legend(fontsize=15, loc="lower right", frameon=False, ncol=1,
                    prop={"weight": "bold"})
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)


# ---------------------------------------------------------------------------
# Panel C — CRA dose-response
# ---------------------------------------------------------------------------

_PANEL_B_CSV = DATA_DIR / "wetlab" / "data" / "Figure4_PanelB.csv"

# TCR grouping: (csv column header, display name from CSV, color, marker)
_PANEL_B_TCRS = [
    ("Nef71-79 -Specific TCR\nCASSQEFPGALYNEQFF\nCAHGSSNTGKLIF",
     "Nef71-79 -Specific TCR\nCASSQEFPGALYNEQFF\nCAHGSSNTGKLIF",
     TCR_PALETTE["TCR1"], "s"),
    ("Nef71-79 -Specific TCR\nCASSQEFPGAAYNEQFF\nCVPLSNTGKLIF",
     "Nef71-79 -Specific TCR\nCASSQEFPGAAYNEQFF\nCVPLSNTGKLIF",
     TCR_PALETTE["TCR2"], "o"),
    ("Vpr34-42 -Specific TCR \nCASSLWGGPSNEQFF\nCAYSGAGSYQLTF",
     "Vpr34-42 -Specific TCR\nCASSLWGGPSNEQFF\nCAYSGAGSYQLTF",
     TCR_PALETTE["TCR3"], "^"),
]


def _sigmoidal_4pl(x, bottom, top, ic50, hill):
    """4-parameter logistic: Y = Bottom + (Top-Bottom)/(1+(IC50/X)^Hill)."""
    return bottom + (top - bottom) / (1 + (ic50 / x) ** hill)


def _fit_ec50(conc: np.ndarray, mean: np.ndarray) -> float | None:
    """Fit 4PL to dose-response data, return IC50 (nM) or None on failure."""
    mask = conc > 0
    try:
        popt, _ = curve_fit(
            _sigmoidal_4pl, conc[mask], mean[mask],
            p0=[5, 75, 0.5, 1], maxfev=10000,
            bounds=([0, 0, 1e-4, 0.1], [50, 100, 100, 10]),
        )
        return round(popt[2], 2)  # IC50
    except Exception:
        return None


def _load_panel_b_data():
    """Parse Panel B CSV: returns (concentrations, {display_name: (mean, sem)})."""
    import csv
    import io

    text = _PANEL_B_CSV.read_text().strip().replace("\r\n", "\n")
    reader = csv.reader(io.StringIO(text))
    rows = list(reader)

    # Header row: first two cols empty, then 3 replicates per TCR (9 total)
    header = rows[0]

    # Data rows: col 0 = label (e.g. "102"), col 1 = concentration, cols 2-10 = values
    concentrations = []
    # Collect raw replicate columns per TCR
    tcr_data = {name: [] for _, name, _, _ in _PANEL_B_TCRS}

    for row in rows[1:]:
        concentrations.append(float(row[1]))
        idx = 2  # start of data columns
        for col_header, name, _, _ in _PANEL_B_TCRS:
            reps = [float(row[idx + j]) for j in range(3)]
            tcr_data[name].append(reps)
            idx += 3

    conc = np.array(concentrations)
    result = {}
    for _, name, _, _ in _PANEL_B_TCRS:
        arr = np.array(tcr_data[name])  # shape (n_conc, 3)
        mean = arr.mean(axis=1)
        sem = arr.std(axis=1, ddof=1) / np.sqrt(arr.shape[1])
        result[name] = (mean, sem)

    return conc, result


def _draw_panel_b(ax) -> None:
    """CRA dose-response from CSV with half-max marker lines and text box."""
    X_ZERO = 3e-3
    conc, tcr_data = _load_panel_b_data()

    curve_data = {}
    ec50_values = {}
    for _, name, color, marker in _PANEL_B_TCRS:
        mean, sem = tcr_data[name]
        x_plot = np.where(conc == 0, X_ZERO, conc)
        disp_name = "\n".join(_subscript_gene_range(ln) for ln in name.split("\n"))
        ax.errorbar(x_plot, mean, yerr=sem,
                    fmt=f"{marker}-", color=color, markersize=8,
                    capsize=3, linewidth=1.8, label=disp_name, zorder=3)
        curve_data[name] = (x_plot, mean)
        ec50_values[name] = _fit_ec50(conc, mean)

    # Draw 1/2 max marker lines
    for _, name, color, _ in _PANEL_B_TCRS:
        hm_x = ec50_values.get(name)
        if hm_x is None:
            continue
        x_pts, y_pts = curve_data[name]

        # Interpolate y at the half-max concentration
        pos = x_pts > 0
        x_log = np.log10(x_pts[pos])
        y_pos = y_pts[pos]
        order = np.argsort(x_log)
        half_y = float(np.interp(np.log10(hm_x), x_log[order], y_pos[order]))

        # Horizontal dotted line from y-axis to half-max x
        ax.plot([300, hm_x], [half_y, half_y],
                ls=":", color=color, lw=1.2, alpha=0.7, zorder=2)
        # Vertical dotted line from half_y down to x-axis
        ax.plot([hm_x, hm_x], [half_y, 0],
                ls=":", color=color, lw=1.2, alpha=0.7, zorder=2)

    # Text box with computed 1/2 max EC50 values
    hm_lines = [r"$\bf{1/2\ max\ EC_{50}:}$"]
    for _, name, _, _ in _PANEL_B_TCRS:
        hm_x = ec50_values.get(name)
        if hm_x is None:
            continue
        short_name = name.split("\n")[0].strip()
        # Build bold text with subscript gene range in mathtext
        def _bold_with_sub(s):
            m = re.match(r'(.+?)(\d+[-–]\d+)(.*)', s)
            if m:
                pre = m.group(1).replace(" ", r"\ ").replace("-", r"\text{-}")
                rng = m.group(2).replace("-", r"\text{-}")
                post = m.group(3).replace(" ", r"\ ").replace("-", r"\text{-}")
                return rf"$\bf{{{pre}}}$$\bf{{_{{{rng}}}}}$$\bf{{{post}}}$"
            safe = s.replace(" ", r"\ ").replace("-", r"\text{-}")
            return rf"$\bf{{{safe}}}$"
        hm_lines.append(f"  {_bold_with_sub(short_name)}: {hm_x} nM")
    ax.text(0.98, 0.98, "\n".join(hm_lines), transform=ax.transAxes,
            fontsize=14, verticalalignment="top", horizontalalignment="right",
            bbox=dict(boxstyle="round,pad=0.4", facecolor="white",
                      edgecolor="grey", alpha=0.9))

    ax.set_xscale("log")
    ax.set_xlabel("Peptide Concentration (nM)", fontsize=20)
    ax.set_ylabel("Specific Lysis (%)", fontsize=20, fontweight="bold")
    ax.set_xlim(1e-3, 5e2)
    ax.set_ylim(0, 100)
    ax.invert_xaxis()
    ax.tick_params(labelsize=18)

    ax.set_xticks([100, 10, 1, 0.1, 0.01, X_ZERO])
    ax.set_xticklabels(
        ["$10^{2}$", "$10^{1}$", "$10^{0}$", "$10^{-1}$", "$10^{-2}$", "0"]
    )

    ax.legend(fontsize=16, loc="center left", bbox_to_anchor=(1.02, 0.5),
              frameon=False, prop={"weight": "bold"})
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)


# ---------------------------------------------------------------------------
# Render
# ---------------------------------------------------------------------------

def render(out_dir: Path) -> Path:
    fig = plt.figure(figsize=(20, 8))
    gs = GridSpec(1, 2, width_ratios=[0.42, 0.58], wspace=0.30,
                  left=0.14, right=0.85, top=0.92, bottom=0.15)

    # Panel A: HLA restriction bars
    ax_a = fig.add_subplot(gs[0])
    _draw_panel_a(ax_a)
    ax_a.set_title("A", loc="left", fontweight="bold", fontsize=24)

    # Panel B: CRA dose-response
    ax_b = fig.add_subplot(gs[1])
    _draw_panel_b(ax_b)
    ax_b.set_title("B", loc="left", fontweight="bold", fontsize=24)

    out = Path(out_dir) / "Figure_04.pdf"
    fig.savefig(out)
    plt.close(fig)
    return out
