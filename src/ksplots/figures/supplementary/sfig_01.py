"""Supplementary Figure 1 — Targeted viral expression heatmap (Log2FC).

Reproduces ``reference/Supplementary_figure_one.pdf`` from NanoString data:

    data/nanostring/nanostring_counts.xlsx  (Sheet 4)

Layout
------
* Rows = patients, grouped by phenotype (Endemic KS top, Epidemic KS bottom).
* Columns = probes, grouped by target category with black-bordered banners:
  Endothelial, EBV, CMV, KSHV, HIV1, HIV2.
* Colorbar: viridis, label "Log2FC".
* Cohort labels ("Endemic KS", "Epidemic KS") on the right side.
"""
from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import Rectangle

from ...config import DATA_DIR

NANOSTRING = DATA_DIR / "nanostring" / "nanostring_counts.xlsx"

# Sample name normalisation (matching notebook logic).
SAMPLE_RENAMES = {
    "8020": "008_020",
    "030-B": "008_030",
    "028-B": "008_028",
    "024-B": "008_024",
    "013-B": "008_013",
    "008-C": "008_008",
    "001-B": "008_001",
}

# Display order: Endemic KS first (top), then Epidemic KS (bottom, reversed).
ENDEMIC = ["008_097", "008_061"]
EPIDEMIC = [
    "008_106", "008_101", "008_099", "008_093", "008_092", "008_088",
    "008_085", "008_075", "008_062", "008_059", "008_052", "008_037",
    "008_036", "008_030", "008_028", "008_024", "008_020", "008_013",
    "008_008", "008_001",
]

# Target display order and label mapping.
TARGET_ORDER = [
    ("Endothelial", "Endothelial"),
    ("HHV4", "EBV"),
    ("HHV5", "CMV"),
    ("HHV8", "KSHV"),
    ("HIV1", "HIV1"),
    ("HIV2", "HIV2"),
]


def _draw_boxed_label(ax, x, y, w, h, text, fontsize=14, fontweight="bold"):
    ax.add_patch(Rectangle((x, y), w, h, transform=ax.transAxes,
                            facecolor="white", edgecolor="black",
                            linewidth=1.2, clip_on=False, zorder=4))
    ax.text(x + w / 2, y + h / 2, text, transform=ax.transAxes,
            ha="center", va="center", fontsize=fontsize,
            fontweight=fontweight, zorder=5)


def render(out_dir: Path) -> Path:
    raw = pd.read_excel(NANOSTRING, sheet_name=3)

    # Normalise sample column names.
    sample_cols = [c for c in raw.columns if c not in ("Target", "Probe Name")]
    rename_map = {}
    for c in sample_cols:
        if c in SAMPLE_RENAMES:
            rename_map[c] = SAMPLE_RENAMES[c]
        else:
            rename_map[c] = c.replace("-", "_")
    raw = raw.rename(columns=rename_map)

    # Build ordered probe list and gene-group spans.
    probes: list[str] = []
    probe_spans: list[tuple[str, int, int]] = []
    for target_id, display_label in TARGET_ORDER:
        sub = raw[raw["Target"] == target_id]["Probe Name"].tolist()
        start = len(probes)
        probes.extend(sub)
        probe_spans.append((display_label, start, len(probes)))

    # Build matrix: rows = samples (ordered), cols = probes.
    sample_order = ENDEMIC + EPIDEMIC
    mat = pd.DataFrame(index=sample_order, columns=probes, dtype=float)
    for _, row in raw.iterrows():
        probe = row["Probe Name"]
        if probe not in probes:
            continue
        for s in sample_order:
            if s in row.index:
                mat.loc[s, probe] = row[s]

    arr = mat.to_numpy(dtype=float)

    fig, ax = plt.subplots(figsize=(20, 12))
    # Clip negative values to 0 (reference uses viridis with breaks 2–10).
    arr = np.clip(arr, 0, None)
    im = ax.imshow(arr, aspect="auto", cmap="viridis", interpolation="nearest",
                   vmin=0, vmax=10)

    # Gene names along the bottom.
    ax.set_xticks(np.arange(len(probes)))
    ax.set_xticklabels(probes, fontsize=13, rotation=90, fontweight="bold")
    ax.tick_params(axis="x", length=2, pad=2)

    # Sample names on the left.
    ax.set_yticks(np.arange(len(sample_order)))
    ax.set_yticklabels(sample_order, fontsize=14, fontweight="bold")
    ax.tick_params(axis="y", length=0, pad=4)

    # White divider between Endemic and Epidemic groups.
    ax.axhline(len(ENDEMIC) - 0.5, color="white", lw=2.5)

    # White dividers between probe groups.
    for _, start, end in probe_spans:
        if start > 0:
            ax.axvline(start - 0.5, color="white", lw=2.0)

    # Target group banners above columns.
    n_cols = len(probes)
    for label, start, end in probe_spans:
        x0 = start / n_cols
        x1 = end / n_cols
        _draw_boxed_label(ax, x=x0, y=1.005, w=x1 - x0, h=0.035,
                          text=label, fontsize=16)

    # Cohort label boxes on the right side.
    n_rows = len(sample_order)
    # Endemic box
    y0_end = 1.0 - len(ENDEMIC) / n_rows
    h_end = len(ENDEMIC) / n_rows
    _draw_boxed_label(ax, x=1.005, y=y0_end, w=0.06, h=h_end,
                      text="Endemic\nKS", fontsize=14)
    # Epidemic box
    y0_epi = 0.0
    h_epi = len(EPIDEMIC) / n_rows
    _draw_boxed_label(ax, x=1.005, y=y0_epi, w=0.06, h=h_epi,
                      text="Epidemic\nKS", fontsize=14)

    ax.set_title("", pad=45)  # space for banners

    # Colorbar.
    cb = fig.colorbar(im, ax=ax, shrink=0.4, pad=0.12)
    cb.set_label("Log2FC", fontsize=15)
    cb.set_ticks([2, 4, 6, 8, 10])
    cb.ax.tick_params(labelsize=13)

    fig.tight_layout()
    out = Path(out_dir) / "Supplementary_Figure_01.pdf"
    fig.savefig(out)
    plt.close(fig)
    return out
