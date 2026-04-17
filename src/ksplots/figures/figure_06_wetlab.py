"""Figure 6 (wetlab) — PEL cell killing.

Reproduces FIGURES_IT.pptx Slide 6:

* **Panel A** — Specific Lysis (%) vs peptide concentration for ORF6 TCRs
  against PEL cell lines (JSC-1, VG-1).
* **Panel B** — mNeonGreen+ (%) bar chart: (-) TPA vs (+) TPA vs OKT3
  for ORF6-Specific TCR 1 (KS TCR 2).
* **Panel C** — IFNγ (pg/mL) bar chart: (-) TPA vs (+) TPA.
* **Panel D** — Same as panel B for ORF6-Specific TCR 2 (KS TCR 3).

Data sources:
  - CRA Prism for dose-response
  - Master Prism for mNeonGreen+ bars
  - ELISA Prism for IFNγ
"""
from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.gridspec import GridSpec

from ..config import (
    PRISM_CRA,
    PRISM_ELISA,
    PRISM_MASTER,
    WETLAB_TCR_COLORS,
)
from ..io.prism import extract_grouped_mean_sem, extract_xy_mean_sem

_TCR2_COLOR = WETLAB_TCR_COLORS["ORF6-Specific TCR 1"]
_TCR3_COLOR = WETLAB_TCR_COLORS["ORF6-Specific TCR 2"]
_NEG_COLOR = "#333333"


def _draw_pel_bars(ax, df, color, ylabel="mNeonGreen+ (%)"):
    """Draw PEL reactivation bar chart."""
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
               color=color, edgecolor="black", linewidth=0.5,
               label=neg_label, error_kw={"linewidth": 0.8},
               alpha=0.4)
        ax.bar(x + bar_w / 2, pos_mean, bar_w, yerr=pos_sem, capsize=3,
               color=color, edgecolor="black", linewidth=0.5,
               label=pos_label, error_kw={"linewidth": 0.8})
    elif len(mean_cols) == 1:
        vals = df[mean_cols[0]].values
        errs = df[sem_cols[0]].values if sem_cols else np.zeros(n)
        ax.bar(x, vals, bar_w * 1.5, yerr=errs, capsize=3,
               color=color, edgecolor="black", linewidth=0.5,
               error_kw={"linewidth": 0.8})

    ax.set_xticks(x)
    ax.set_xticklabels(labels, fontsize=13, rotation=30, ha="right")
    ax.set_ylabel(ylabel, fontsize=14, fontweight="bold")
    ax.tick_params(labelsize=12)
    ax.legend(fontsize=10, frameon=False)
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)


def _draw_pel_titration(ax, df, tcr2_color, tcr3_color):
    """Draw PEL killing dose-response for two TCRs."""
    mean_cols = [c for c in df.columns if c.endswith("_mean")]

    for col, color, label in zip(
        mean_cols[:2],
        [tcr2_color, tcr3_color],
        ["ORF6-Specific TCR 1", "ORF6-Specific TCR 2"],
    ):
        sem_col = col.replace("_mean", "_sem")
        x = df["x"].values
        y = df[col].values
        yerr = df[sem_col].values if sem_col in df.columns else None

        ax.errorbar(x, y, yerr=yerr, fmt="o-", color=color, markersize=6,
                    capsize=3, linewidth=1.5, label=label)

    ax.set_xscale("log")
    ax.set_xlabel("Peptide Concentration (nM)", fontsize=14)
    ax.set_ylabel("Specific Lysis (%)", fontsize=14, fontweight="bold")
    ax.set_ylim(0, 80)
    ax.tick_params(labelsize=12)
    ax.legend(fontsize=11, frameon=False)
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)


def render(out_dir: Path) -> Path:
    # Panel B: mNeonGreen+ bars — KS TCR 2 with JSC-1/VG-1 reactivated
    df_tcr2_react = extract_grouped_mean_sem(
        PRISM_MASTER,
        "JSC-1 and VG-1 reactivated_KS TCR 2_24 hr reactivation_v3",
    )

    # Panel D: same for KS TCR 3
    df_tcr3_react = extract_grouped_mean_sem(
        PRISM_MASTER,
        "JSC-1 and VG-1 reactivated_KS TCR 3_24 hr reactivation_v3",
    )

    # Panel C: IFNγ ELISA — KS TCR 2 PEL no peptide
    df_ifng = extract_grouped_mean_sem(
        PRISM_ELISA,
        "ELISA_16 hr co-culture_KS TCR 2_PEL_Newly transduced CD8s_NO PEPTIDE_12Mar2025",
    )

    fig = plt.figure(figsize=(20, 12))
    gs = GridSpec(2, 2, hspace=0.40, wspace=0.30,
                  left=0.08, right=0.95, top=0.93, bottom=0.10)

    # Panel A: dose-response (try to load)
    ax_a = fig.add_subplot(gs[0, 0])
    try:
        df_titration = extract_xy_mean_sem(
            PRISM_CRA,
            "CRA_KS TCR 2 with JSC-1 AND VG-1_Peptide Titration_12Mar2025",
        )
        _draw_pel_titration(ax_a, df_titration, _TCR2_COLOR, _TCR3_COLOR)
    except KeyError:
        # Try alternate table name
        try:
            df_titration = extract_xy_mean_sem(
                PRISM_CRA,
                "CRA_KS TCR 2 and 3 with JSC-1_Peptide Titration_12Mar2025",
            )
            _draw_pel_titration(ax_a, df_titration, _TCR2_COLOR, _TCR3_COLOR)
        except KeyError:
            ax_a.text(0.5, 0.5, "Panel A: titration data\nnot found",
                      transform=ax_a.transAxes, ha="center", va="center",
                      fontsize=14)
    ax_a.set_title("A", loc="left", fontweight="bold", fontsize=18)

    # Panel B: mNeonGreen+ KS TCR 2
    ax_b = fig.add_subplot(gs[0, 1])
    _draw_pel_bars(ax_b, df_tcr2_react, _TCR2_COLOR)
    ax_b.set_ylim(0, 100)
    ax_b.set_title("B", loc="left", fontweight="bold", fontsize=18)
    ax_b.text(0.5, 1.05, "ORF6-Specific TCR 1", transform=ax_b.transAxes,
              ha="center", fontsize=13, fontweight="bold", color=_TCR2_COLOR)

    # Panel C: IFNγ
    ax_c = fig.add_subplot(gs[1, 0])
    _draw_pel_bars(ax_c, df_ifng, _TCR2_COLOR,
                   ylabel="IFN\u03b3 (pg/mL)")
    ax_c.set_title("C", loc="left", fontweight="bold", fontsize=18)

    # Panel D: mNeonGreen+ KS TCR 3
    ax_d = fig.add_subplot(gs[1, 1])
    _draw_pel_bars(ax_d, df_tcr3_react, _TCR3_COLOR)
    ax_d.set_ylim(0, 100)
    ax_d.set_title("D", loc="left", fontweight="bold", fontsize=18)
    ax_d.text(0.5, 1.05, "ORF6-Specific TCR 2", transform=ax_d.transAxes,
              ha="center", fontsize=13, fontweight="bold", color=_TCR3_COLOR)

    out = Path(out_dir) / "Figure_06_wetlab.pdf"
    fig.savefig(out)
    plt.close(fig)
    return out
