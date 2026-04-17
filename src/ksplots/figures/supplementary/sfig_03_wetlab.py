"""Supplementary Figure (wetlab) — COS7 Jurkat reporter assays.

Reproduces FIGURES_IT.pptx Slide 3:

* **Panel A** — ORF57 + segments: mNeonGreen+ (%) ± HLA-A*66:01 (KS TCR 21).
* **Panel B** — ORF59 peptide titration (KS TCR 24).
* **Panels C–D** — ORF6 segment mapping: mNeonGreen+ (%) ± HLA-B*45:01
  for KS TCR 2 (C) and KS TCR 3 (D), showing ORF6 full, 1-454, 455-1132.

Data source: ``data/wetlab/Experiments/.../Data Analysis.prism``
"""
from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.gridspec import GridSpec
from scipy.stats import ttest_ind

from ...config import PRISM_MASTER, WETLAB_TCR_COLORS
from ...io.prism import extract_grouped_mean_sem, extract_table

# Panel colors
_TCR2_COLOR = WETLAB_TCR_COLORS["ORF6-Specific TCR 1"]   # blue
_TCR3_COLOR = WETLAB_TCR_COLORS["ORF6-Specific TCR 2"]   # light blue
_TCR21_COLOR = WETLAB_TCR_COLORS["ORF57-Specific TCR"]    # dark blue
_NEG_COLOR = "#333333"  # dark grey for (-) HLA


def _draw_grouped_bars(
    ax,
    df,
    pos_color: str,
    neg_color: str = _NEG_COLOR,
    ylabel: str = "mNeonGreen+ (%)",
    ylim: float = 80,
) -> None:
    """Draw paired bar chart: (-) HLA vs (+) HLA."""
    labels = df["label"].tolist()
    n = len(labels)
    x = np.arange(n)
    bar_w = 0.35

    # Find mean/sem columns
    mean_cols = [c for c in df.columns if c.endswith("_mean")]
    sem_cols = [c for c in df.columns if c.endswith("_sem")]

    if len(mean_cols) >= 2:
        # Two conditions: (-) HLA, (+) HLA
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

        # Significance testing between (-) and (+) HLA for each condition
        raw_df = _get_raw_for_significance(df, mean_cols)
        if raw_df is not None:
            for i in range(n):
                neg_vals = raw_df.iloc[i, :3].dropna().values
                pos_vals = raw_df.iloc[i, 3:6].dropna().values
                if len(neg_vals) >= 2 and len(pos_vals) >= 2:
                    _, p = ttest_ind(neg_vals, pos_vals)
                    if p < 0.001:
                        _add_significance_bracket(
                            ax, i - bar_w / 2, i + bar_w / 2,
                            max(pos_mean[i] + pos_sem[i],
                                neg_mean[i] + neg_sem[i]),
                            "***", ylim,
                        )

    elif len(mean_cols) == 1:
        # Single condition (e.g., peptide titration)
        vals = df[mean_cols[0]].values
        errs = df[sem_cols[0]].values if sem_cols else np.zeros(n)
        ax.bar(x, vals, bar_w * 1.5, yerr=errs, capsize=3,
               color=pos_color, edgecolor="black", linewidth=0.5,
               error_kw={"linewidth": 0.8})

    ax.set_xticks(x)
    ax.set_xticklabels(labels, fontsize=12, rotation=30, ha="right")
    ax.set_ylabel(ylabel, fontsize=14, fontweight="bold")
    ax.set_ylim(0, ylim)
    ax.tick_params(labelsize=12)
    handles, labels = ax.get_legend_handles_labels()
    if handles:
        ax.legend(fontsize=10, frameon=False, loc="upper left")
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)


def _get_raw_for_significance(df, mean_cols):
    """Get raw replicate values for significance testing."""
    # The mean/sem were computed from replicate columns
    # We need the original replicate data
    rep_cols = [c for c in df.columns
                if c != "label" and not c.endswith("_mean") and not c.endswith("_sem")]
    if rep_cols:
        return df[rep_cols]
    return None


def _add_significance_bracket(ax, x1, x2, y_max, text, ylim):
    """Draw a significance bracket above bars."""
    y_bar = min(y_max * 1.05, ylim * 0.92)
    ax.plot([x1, x1, x2, x2],
            [y_bar, y_bar * 1.03, y_bar * 1.03, y_bar],
            color="black", lw=1.0, clip_on=False)
    ax.text((x1 + x2) / 2, y_bar * 1.04, text,
            ha="center", va="bottom", fontsize=12, fontweight="bold")


def render(out_dir: Path) -> Path:
    # Panel A: ORF57 + segments (KS TCR 21)
    df_orf57 = extract_grouped_mean_sem(
        PRISM_MASTER,
        "FINAL_ORF57 and SEGMENTS and peptide_KS TCR 21",
    )

    # Panel B: ORF59 (KS TCR 24) — this is a titration, use as bars
    df_orf59 = extract_grouped_mean_sem(
        PRISM_MASTER,
        "KS TCR 24 ORF59 Peptide Titration_FINAL",
    )

    # Panel C: ORF6 segments (KS TCR 2)
    df_orf6_tcr2 = extract_grouped_mean_sem(
        PRISM_MASTER,
        "FINAL_KS TCR 2_ORF6_Segment A and B_Segment A5",
    )

    # Panel D: ORF6 segments (KS TCR 3)
    df_orf6_tcr3 = extract_grouped_mean_sem(
        PRISM_MASTER,
        "FINAL_KS TCR 3_ORF6_Segment A and B_Segment A5",
    )

    fig = plt.figure(figsize=(22, 14))
    gs = GridSpec(2, 3, hspace=0.45, wspace=0.35,
                  left=0.07, right=0.95, top=0.94, bottom=0.08)

    # Panel A: ORF57
    ax_a = fig.add_subplot(gs[0, 0])
    _draw_grouped_bars(ax_a, df_orf57, _TCR21_COLOR, ylim=60)
    ax_a.set_title("A", loc="left", fontweight="bold", fontsize=18)

    # Panel B: ORF59
    ax_b = fig.add_subplot(gs[0, 1])
    _draw_grouped_bars(ax_b, df_orf59, WETLAB_TCR_COLORS["ORF59-Specific TCR"],
                       ylim=80)
    ax_b.set_title("B", loc="left", fontweight="bold", fontsize=18)

    # Panel C: ORF6 (KS TCR 2)
    ax_c = fig.add_subplot(gs[1, 0])
    _draw_grouped_bars(ax_c, df_orf6_tcr2, _TCR2_COLOR, ylim=80)
    ax_c.set_title("C", loc="left", fontweight="bold", fontsize=18)
    ax_c.text(0.5, 1.08, "ORF6-Specific TCR 1", transform=ax_c.transAxes,
              ha="center", fontsize=13, fontweight="bold", color=_TCR2_COLOR)

    # Panel D: ORF6 (KS TCR 3)
    ax_d = fig.add_subplot(gs[1, 1])
    _draw_grouped_bars(ax_d, df_orf6_tcr3, _TCR3_COLOR, ylim=80)
    ax_d.set_title("D", loc="left", fontweight="bold", fontsize=18)
    ax_d.text(0.5, 1.08, "ORF6-Specific TCR 2", transform=ax_d.transAxes,
              ha="center", fontsize=13, fontweight="bold", color=_TCR3_COLOR)

    out = Path(out_dir) / "Supplementary_Figure_COS7_reporter.pdf"
    fig.savefig(out)
    plt.close(fig)
    return out
