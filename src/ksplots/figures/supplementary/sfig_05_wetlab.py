"""Supplementary Figure (wetlab) — ORF6 fine epitope mapping.

Reproduces FIGURES_IT.pptx Slide 4:

* **Panel A** — ORF6 fine-mapping for KS TCR 2 (ORF6-Specific TCR 1):
  ORF6 full, 1-454, 311-375, 311-339, 329-337, AEQALHIGA peptide, OKT3.
* **Panel B** — Same for KS TCR 3 (ORF6-Specific TCR 2).

Data source: ``data/wetlab/Experiments/.../Data Analysis.prism``
"""
from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.gridspec import GridSpec

from ...config import PRISM_MASTER, WETLAB_TCR_COLORS
from ...io.prism import extract_grouped_mean_sem

_TCR2_COLOR = WETLAB_TCR_COLORS["ORF6-Specific TCR 1"]
_TCR3_COLOR = WETLAB_TCR_COLORS["ORF6-Specific TCR 2"]
_NEG_COLOR = "#333333"


def _draw_fine_mapping(ax, df, pos_color, neg_color=_NEG_COLOR):
    """Draw ORF6 fine-mapping paired bars."""
    labels = df["label"].tolist()
    n = len(labels)
    x = np.arange(n)
    bar_w = 0.35

    mean_cols = [c for c in df.columns if c.endswith("_mean")]
    sem_cols = [c for c in df.columns if c.endswith("_sem")]

    if len(mean_cols) >= 2:
        neg_mean = df[mean_cols[0]].values
        neg_sem = df[sem_cols[0]].values
        pos_mean = df[mean_cols[1]].values
        pos_sem = df[sem_cols[1]].values

        neg_label = mean_cols[0].replace("_mean", "")
        pos_label = mean_cols[1].replace("_mean", "")

        ax.bar(x - bar_w / 2, neg_mean, bar_w, yerr=neg_sem, capsize=3,
               color=neg_color, edgecolor="black", linewidth=0.5,
               label=neg_label, error_kw={"linewidth": 0.8})
        ax.bar(x + bar_w / 2, pos_mean, bar_w, yerr=pos_sem, capsize=3,
               color=pos_color, edgecolor="black", linewidth=0.5,
               label=pos_label, error_kw={"linewidth": 0.8})

    # Add bracket for "ORF6 residues (aa)" under the truncation bars
    # Find indices for residue-labeled bars
    residue_indices = [
        i for i, label in enumerate(labels)
        if any(c.isdigit() for c in label) and "OKT3" not in label
        and "AEQ" not in label and "ORF6" not in label
    ]
    if len(residue_indices) >= 2:
        x_start = residue_indices[0] - 0.5
        x_end = residue_indices[-1] + 0.5
        y_bracket = -0.12 * ax.get_ylim()[1]
        ax.annotate("", xy=(x_start, y_bracket), xytext=(x_end, y_bracket),
                    arrowprops=dict(arrowstyle="-", lw=1.2),
                    annotation_clip=False)
        ax.text((x_start + x_end) / 2, y_bracket * 1.4,
                "ORF6 residues\n(aa)", ha="center", va="top",
                fontsize=11, fontweight="bold", clip_on=False)

    ax.set_xticks(x)
    ax.set_xticklabels(labels, fontsize=12, rotation=30, ha="right")
    ax.set_ylabel("mNeonGreen+ (%)", fontsize=14, fontweight="bold")
    ax.set_ylim(0, 80)
    ax.tick_params(labelsize=12)
    ax.legend(fontsize=10, frameon=False)
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)


def render(out_dir: Path) -> Path:
    df_tcr2 = extract_grouped_mean_sem(
        PRISM_MASTER, "ORF6 Segments to AEQ_KS TCR 2_FINAL"
    )
    df_tcr3 = extract_grouped_mean_sem(
        PRISM_MASTER, "ORF6 Segments to AEQ_KS TCR 3_FINAL"
    )

    fig, axes = plt.subplots(1, 2, figsize=(18, 7))

    _draw_fine_mapping(axes[0], df_tcr2, _TCR2_COLOR)
    axes[0].set_title("A", loc="left", fontweight="bold", fontsize=18)
    axes[0].text(0.5, 1.05, "ORF6-Specific TCR 1", transform=axes[0].transAxes,
                 ha="center", fontsize=14, fontweight="bold", color=_TCR2_COLOR)

    _draw_fine_mapping(axes[1], df_tcr3, _TCR3_COLOR)
    axes[1].set_title("B", loc="left", fontweight="bold", fontsize=18)
    axes[1].text(0.5, 1.05, "ORF6-Specific TCR 2", transform=axes[1].transAxes,
                 ha="center", fontsize=14, fontweight="bold", color=_TCR3_COLOR)

    fig.tight_layout(rect=[0, 0.05, 1, 0.95])
    out = Path(out_dir) / "Supplementary_Figure_ORF6_mapping.pdf"
    fig.savefig(out)
    plt.close(fig)
    return out
