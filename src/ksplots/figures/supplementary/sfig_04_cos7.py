"""Supplementary Figure 4 — COS7 Jurkat reporter assays, ORF6 fine mapping.

Reproduces FIGURES_IT.pptx Slide 4:

* **Panel A** — (placeholder, not yet implemented).
* **Panel B** — ORF6 fine epitope mapping ± HLA-B*45:01:
  ORF6, 1-454, 311-375, 311-339, 329-337, AEQALHIGA peptide, OKT3.
* **Panel C** — ORF6 segment mapping ± HLA-B*45:01:
  ORF6, Segment A, A5, A5c, AEQALHIGA, AEQALHIGA peptide, OKT3.

Data sources:
  - Panel B: ``data/wetlab/data/SpFigure3_PanelC.csv``
  - Panel C: ``data/wetlab/data/SpFigure3_PanelD.csv``
"""
from __future__ import annotations

import csv
import io
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.gridspec import GridSpec
from matplotlib.patches import Patch
from scipy.stats import ttest_ind

from ...config import WETLAB_DIR

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------

_PANEL_CSV_DIR = WETLAB_DIR / "data"

# ---------------------------------------------------------------------------
# Colours
# ---------------------------------------------------------------------------

_NEG_COLOR = "#999999"   # grey for (-) HLA
_OKT3_COLOR = "black"    # black for OKT3 positive control
_TCR1_COLOR = "#D4A017"  # golden yellow for ORF6-Specific TCR 1
_TCR2_COLOR = "#008080"  # teal for ORF6-Specific TCR 2
_HLA_ALLELE = "B*45:01"

# ---------------------------------------------------------------------------
# CSV loading
# ---------------------------------------------------------------------------


def _load_panel_csv(path: Path) -> pd.DataFrame:
    """Load a SpFigure3 panel CSV into a DataFrame with mean/SEM columns."""
    text = path.read_text().strip()
    reader = csv.reader(io.StringIO(text))
    rows = list(reader)

    header = rows[0]
    neg_label = header[1].strip()
    pos_label = header[4].strip()

    records = []
    for row in rows[1:]:
        label = row[0].strip()
        vals = []
        for v in row[1:7]:
            v = v.strip() if v else ""
            vals.append(float(v) if v else np.nan)
        neg_vals = np.array(vals[:3])
        pos_vals = np.array(vals[3:6])

        rec = {
            "label": label,
            f"{neg_label}_mean": np.nanmean(neg_vals),
            f"{neg_label}_sem": (
                np.nanstd(neg_vals, ddof=1) / np.sqrt(np.sum(~np.isnan(neg_vals)))
                if np.sum(~np.isnan(neg_vals)) > 1 else 0.0
            ),
            f"{pos_label}_mean": np.nanmean(pos_vals),
            f"{pos_label}_sem": (
                np.nanstd(pos_vals, ddof=1) / np.sqrt(np.sum(~np.isnan(pos_vals)))
                if np.sum(~np.isnan(pos_vals)) > 1 else 0.0
            ),
            f"{neg_label}_r1": vals[0],
            f"{neg_label}_r2": vals[1],
            f"{neg_label}_r3": vals[2],
            f"{pos_label}_r1": vals[3],
            f"{pos_label}_r2": vals[4],
            f"{pos_label}_r3": vals[5],
        }
        records.append(rec)

    return pd.DataFrame(records)


# ---------------------------------------------------------------------------
# Bar chart drawing
# ---------------------------------------------------------------------------


def _draw_grouped_bars(
    ax: plt.Axes,
    df: pd.DataFrame,
    pos_color: str,
    hla_allele: str = "HLA",
    neg_color: str = _NEG_COLOR,
    ylabel: str = "mNeonGreen+ (%)",
    ylim: float = 80,
) -> None:
    """Draw paired bar chart: (-) HLA vs (+) HLA with significance brackets."""
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

        pos_legend = f"(+) HLA-{hla_allele}"

        is_okt3 = [lbl.upper().strip().startswith("OKT3") for lbl in labels]
        pos_colors = [_OKT3_COLOR if okt else pos_color for okt in is_okt3]

        ax.bar(x - bar_w / 2, neg_mean, bar_w, yerr=neg_sem, capsize=3,
               color=neg_color, edgecolor="black", linewidth=0.5,
               label="(-) HLA", error_kw={"linewidth": 0.8})

        bars_pos = ax.bar(x + bar_w / 2, pos_mean, bar_w, yerr=pos_sem,
                          capsize=3, edgecolor="black", linewidth=0.5,
                          color=pos_colors, error_kw={"linewidth": 0.8})
        for idx, okt in enumerate(is_okt3):
            if not okt:
                bars_pos[idx].set_label(pos_legend)
                break

        # Significance testing
        rep_cols = [c for c in df.columns
                    if c != "label" and not c.endswith("_mean") and not c.endswith("_sem")]
        if rep_cols:
            for i in range(n):
                neg_vals = df.iloc[i][[c for c in rep_cols if c.startswith(neg_label)]].dropna().values
                pos_vals = df.iloc[i][[c for c in rep_cols if c.startswith(pos_label)]].dropna().values
                if len(neg_vals) >= 2 and len(pos_vals) >= 2:
                    _, p = ttest_ind(neg_vals.astype(float), pos_vals.astype(float))
                    if p < 0.001:
                        y_bar = min(
                            max(pos_mean[i] + pos_sem[i], neg_mean[i] + neg_sem[i]) * 1.05,
                            ylim * 0.92,
                        )
                        ax.plot(
                            [x[i] - bar_w / 2, x[i] - bar_w / 2,
                             x[i] + bar_w / 2, x[i] + bar_w / 2],
                            [y_bar, y_bar * 1.03, y_bar * 1.03, y_bar],
                            color="black", lw=1.0, clip_on=False,
                        )
                        ax.text(x[i], y_bar * 1.04, "***",
                                ha="center", va="bottom", fontsize=18, fontweight="bold")

    ax.set_xticks(x)
    ax.set_xticklabels(labels, fontsize=18, rotation=30, ha="right")
    ax.set_ylabel(ylabel, fontsize=20, fontweight="bold")
    ax.set_ylim(0, ylim)
    ax.tick_params(labelsize=18)
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)


# ---------------------------------------------------------------------------
# Render
# ---------------------------------------------------------------------------


def render(out_dir: Path) -> Path:
    fig = plt.figure(figsize=(22, 10))

    gs = GridSpec(
        1, 2, width_ratios=[7, 7],
        left=0.08, right=0.95, top=0.88, bottom=0.22, wspace=0.35,
    )

    # Panel B — ORF6 fine epitope mapping (ORF6-Specific TCR 1)
    ax_b = fig.add_subplot(gs[0, 0])
    df_b = _load_panel_csv(_PANEL_CSV_DIR / "SpFigure3_PanelC.csv")
    _draw_grouped_bars(ax_b, df_b, _TCR1_COLOR,
                       hla_allele="TCR1", ylim=80)
    ax_b.set_title("B", loc="left", fontweight="bold", fontsize=28)

    # Panel C — ORF6 segment mapping (ORF6-Specific TCR 2)
    ax_c = fig.add_subplot(gs[0, 1])
    df_c = _load_panel_csv(_PANEL_CSV_DIR / "SpFigure3_PanelD.csv")
    _draw_grouped_bars(ax_c, df_c, _TCR2_COLOR,
                       hla_allele="TCR2", ylim=80)
    ax_c.set_title("C", loc="left", fontweight="bold", fontsize=28)

    # Shared legend with mathtext subscripts for ORF6 residue range
    _tcr1_label = (
        r"(+) HLA-B*45:01; ORF6$_{329\text{-}337}$-Specific TCR;"
        " CASSIAGHEQFF / CAVAASGGYQKVTF"
    )
    _tcr2_label = (
        r"(+) HLA-B*45:01; ORF6$_{329\text{-}337}$-Specific TCR;"
        " CASSIAGHEQYF / CAVGASGGYQKVTF"
    )

    ordered: list[tuple[str, plt.Artist]] = [
        ("(-) HLA", Patch(facecolor=_NEG_COLOR, edgecolor="black", lw=0.5)),
        (_tcr1_label, Patch(facecolor=_TCR1_COLOR, edgecolor="black", lw=0.5)),
        (_tcr2_label, Patch(facecolor=_TCR2_COLOR, edgecolor="black", lw=0.5)),
        ("OKT3", Patch(facecolor=_OKT3_COLOR, edgecolor="black", lw=0.5)),
    ]

    fig.legend(
        [h for _, h in ordered],
        [l for l, _ in ordered],
        loc="lower center", bbox_to_anchor=(0.5, 0.0),
        ncol=2, fontsize=14, frameon=False,
    )

    out = Path(out_dir) / "Supplementary_Figure_04.pdf"
    fig.savefig(out)
    plt.close(fig)
    return out
