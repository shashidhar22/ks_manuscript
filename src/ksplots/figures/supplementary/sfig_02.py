"""Supplementary Figure 2 — CibersortX immune deconvolution.

Reproduces ``reference/Supplementary_figure_two.pdf``:

* **Panel A** — Heatmap of cell type abundance (%) across samples grouped by
  study and tissue type.
* **Panel B** — Macrophages M2 vs T cells CD8 scatter with per-tissue-type
  regression lines and R/p annotations, faceted by study.

Data sources: ``data/cibersort/lid_*.csv``, ``tso_epidemic.csv``,
``gtex_nes.csv``, ``gtex_ses.csv``.
"""
from __future__ import annotations

import re
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.gridspec import GridSpec
from matplotlib.patches import Rectangle

from ...config import DATA_DIR
from ...io.metadata import load_study_metadata

CIBERSORT_DIR = DATA_DIR / "cibersort"
DROP_COLS = {"P-value", "Correlation", "RMSE", "Absolute score (sig.score)"}

CELL_TYPE_ORDER = [
    "T cells regulatory (Tregs)",
    "T cells gamma delta",
    "T cells follicular helper",
    "T cells CD8",
    "T cells CD4 naive",
    "T cells CD4 memory resting",
    "T cells CD4 memory activated",
    "Plasma cells",
    "NK cells resting",
    "NK cells activated",
    "Neutrophils",
    "Monocytes",
    "Mast cells resting",
    "Mast cells activated",
    "Macrophages M2",
    "Macrophages M1",
    "Macrophages M0",
    "Eosinophils",
    "Dendritic cells resting",
    "Dendritic cells activated",
    "B cells naive",
    "B cells memory",
]

STUDY_GROUPS = [
    ("Uganda", [
        "hippos_endemic.csv", "hippos_epidemic.csv",
    ]),
    ("Tanzania 1", [
        "lid_batchone.csv", "lid_batchtwo.csv",
        "lid_batchthree.csv", "lid_batchfour.csv",
    ]),
    ("Tanzania 2", ["tso_epidemic.csv"]),
    ("Non sun-exposed skin", ["gtex_nes.csv"]),
    ("Sun-exposed skin", ["gtex_ses.csv"]),
]

UGANDA_SUBGROUP_ORDER = ["Endemic KS", "Epidemic KS"]
TANZANIA1_SUBGROUP_ORDER = ["Con-\ntrol", "Ax. Skin", "Ende-\nmic KS", "Epidemic KS"]
TANZANIA2_SUBGROUP_ORDER = ["Ax. Skin", "Epide-\nmic KS"]
CONTROL_SUBGROUP = ["Control"]

# Tissue-type palette for Panel B (matching notebook).
TISSUE_COLORS = {
    "Control": "#67a628",
    "Endemic KS - NAT": "#80b1d3",
    "Endemic KS - Tumor": "#377eb8",
    "Epidemic KS - NAT": "#fb8072",
    "Epidemic KS - Tumor": "#e41a1c",
}


def _load_percentage(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path)
    df = df.rename(columns={df.columns[0]: "sample"})
    df["sample"] = df["sample"].str.replace("_quant", "", regex=False)
    sig = df.get("Absolute score (sig.score)")
    cell_cols = [c for c in df.columns if c not in {"sample"} | DROP_COLS]
    if sig is not None:
        for c in cell_cols:
            df[c] = (df[c] * 100) / sig
    return df[["sample"] + cell_cols]


def _draw_boxed_label(ax, x, y, w, h, text, fontsize=13, fontweight="bold"):
    ax.add_patch(Rectangle((x, y), w, h, transform=ax.transAxes,
                            facecolor="white", edgecolor="black",
                            linewidth=1.0, clip_on=False, zorder=4))
    ax.text(x + w / 2, y + h / 2, text, transform=ax.transAxes,
            ha="center", va="center", fontsize=fontsize,
            fontweight=fontweight, zorder=5)


def _load_data() -> pd.DataFrame:
    """Load all cibersort data with study/tissue_type annotations."""
    meta = load_study_metadata()
    meta_lookup = {}
    for _, row in meta.iterrows():
        rid = row.get("repertoire_id", "")
        if pd.notna(rid):
            meta_lookup[rid] = (
                str(row.get("tissue_type", "")),
                str(row.get("phenotype", "")),
            )

    records: list[dict] = []
    for study_label, files in STUDY_GROUPS:
        for fname in files:
            path = CIBERSORT_DIR / fname
            if not path.exists():
                continue
            df = _load_percentage(path)
            cell_cols = [c for c in df.columns if c != "sample"]
            for _, row in df.iterrows():
                sample = row["sample"]
                srr = re.search(r"SRR\d+", sample)
                lookup_key = srr.group() if srr else sample
                tissue, pheno = meta_lookup.get(lookup_key, ("", ""))

                # Heatmap subgroup (with line-break formatting).
                if study_label == "Uganda":
                    if "endemic" in fname:
                        subgroup = "Endemic KS"
                        tissue_type = "Endemic KS - Tumor"
                    else:
                        subgroup = "Epidemic KS"
                        tissue_type = "Epidemic KS - Tumor"
                elif study_label == "Tanzania 1":
                    if tissue == "Normal":
                        subgroup = "Con-\ntrol"
                    elif tissue == "Tumor" and pheno == "Endemic KS":
                        subgroup = "Ende-\nmic KS"
                    elif tissue == "Tumor" and pheno == "Epidemic KS":
                        subgroup = "Epidemic KS"
                    elif tissue == "NAT":
                        subgroup = "Ax. Skin"
                    else:
                        subgroup = tissue or "Unknown"
                elif study_label == "Tanzania 2":
                    if tissue == "Tumor":
                        subgroup = "Epide-\nmic KS"
                    elif tissue == "NAT":
                        subgroup = "Ax. Skin"
                    else:
                        subgroup = tissue or "Epide-\nmic KS"
                else:
                    subgroup = "Control"

                # Panel B tissue type (no line breaks) — skip if already set (Uganda).
                if study_label != "Uganda":
                    if study_label in ("Non sun-exposed skin", "Sun-exposed skin"):
                        tissue_type = "Control"
                    elif tissue == "Normal":
                        tissue_type = "Control"
                    elif tissue == "Tumor":
                        tissue_type = f"{pheno} - Tumor"
                    elif tissue == "NAT":
                        tissue_type = f"{pheno} - NAT"
                    else:
                        tissue_type = "Control"

                rec = {
                    "study": study_label, "subgroup": subgroup,
                    "tissue_type": tissue_type, "sample": sample,
                }
                for c in cell_cols:
                    rec[c] = row[c]
                records.append(rec)

    return pd.DataFrame(records)


def _draw_panel_a(ax, rdf: pd.DataFrame) -> None:
    cell_cols = [c for c in rdf.columns
                 if c not in {"study", "subgroup", "tissue_type", "sample"}]
    cell_cols_ordered = [c for c in CELL_TYPE_ORDER if c in cell_cols]
    cell_cols_ordered += [c for c in cell_cols if c not in cell_cols_ordered]

    subgroup_orders = {
        "Uganda": UGANDA_SUBGROUP_ORDER,
        "Tanzania 1": TANZANIA1_SUBGROUP_ORDER,
        "Tanzania 2": TANZANIA2_SUBGROUP_ORDER,
        "Non sun-exposed skin": CONTROL_SUBGROUP,
        "Sun-exposed skin": CONTROL_SUBGROUP,
    }

    ordered_indices: list[int] = []
    study_spans: list[tuple[str, int, int]] = []
    subgroup_spans: list[tuple[str, int, int]] = []

    for study_label, _ in STUDY_GROUPS:
        study_start = len(ordered_indices)
        sub_order = subgroup_orders.get(study_label, CONTROL_SUBGROUP)
        study_data = rdf[rdf["study"] == study_label]
        for sg in sub_order:
            sg_start = len(ordered_indices)
            sg_data = study_data[study_data["subgroup"] == sg]
            ordered_indices.extend(sg_data.index.tolist())
            if len(sg_data) > 0:
                subgroup_spans.append((sg, sg_start, len(ordered_indices)))
        study_spans.append((study_label, study_start, len(ordered_indices)))

    mat = rdf.loc[ordered_indices, cell_cols_ordered].values.T
    n_cells, n_samples = mat.shape

    im = ax.imshow(mat, aspect="auto", cmap="plasma", interpolation="nearest",
                   vmin=0, vmax=50)

    # Fine white grid between tiles.
    ax.set_xticks(np.arange(n_samples + 1) - 0.5, minor=True)
    ax.set_yticks(np.arange(n_cells + 1) - 0.5, minor=True)
    ax.grid(which="minor", color="white", linewidth=0.5)
    ax.tick_params(which="minor", length=0)

    ax.set_yticks(np.arange(n_cells))
    ax.set_yticklabels(cell_cols_ordered, fontsize=15, fontweight="bold")
    ax.tick_params(axis="y", length=0, pad=4)
    ax.set_xticks([])
    ax.set_ylabel("Cell type", fontsize=17, fontweight="bold")

    # Thicker white dividers between study/subgroups (on top of grid).
    for _, start, end in study_spans:
        if start > 0:
            ax.axvline(start - 0.5, color="white", lw=3.0, zorder=3)
    for _, start, end in subgroup_spans:
        if start > 0:
            ax.axvline(start - 0.5, color="white", lw=1.8, zorder=3)

    for label, start, end in study_spans:
        x0 = start / n_samples
        x1 = end / n_samples
        _draw_boxed_label(ax, x=x0, y=1.06, w=x1 - x0, h=0.04,
                          text=label, fontsize=16)

    for label, start, end in subgroup_spans:
        x0 = start / n_samples
        x1 = end / n_samples
        _draw_boxed_label(ax, x=x0, y=1.005, w=x1 - x0, h=0.05,
                          text=label, fontsize=14)

    ax.text(-0.08, 1.12, "A", transform=ax.transAxes,
            fontsize=24, fontweight="bold", va="top")
    return im


PANEL_B_CELL_TYPES = [
    "Macrophages M2",
    "T cells CD4 memory activated",
    "T cells CD4 memory resting",
    "T cells CD8",
]

PANEL_B_COLORS = {
    "Endemic KS": "#377eb8",
    "Epidemic KS": "#e41a1c",
}


def _draw_panel_b(axes: list, rdf: pd.DataFrame) -> None:
    """Per-cell-type boxplots: Endemic vs Epidemic KS tumors with stats."""
    import seaborn as sns
    from scipy.stats import mannwhitneyu

    # Filter to tumor samples only.
    tumor = rdf[rdf["tissue_type"].isin(
        ["Endemic KS - Tumor", "Epidemic KS - Tumor"]
    )].copy()
    tumor["phenotype"] = tumor["tissue_type"].str.replace(" - Tumor", "")

    for i, ct in enumerate(PANEL_B_CELL_TYPES):
        ax = axes[i]
        if ct not in tumor.columns:
            ax.set_visible(False)
            continue

        plot_df = tumor[["phenotype", ct]].copy()
        plot_df = plot_df.rename(columns={ct: "value"})

        sns.boxplot(data=plot_df, x="phenotype", y="value",
                    hue="phenotype", order=["Endemic KS", "Epidemic KS"],
                    hue_order=["Endemic KS", "Epidemic KS"],
                    palette=PANEL_B_COLORS, legend=False,
                    fliersize=0, linewidth=0.8, ax=ax, width=0.5)
        sns.stripplot(data=plot_df, x="phenotype", y="value",
                      hue="phenotype", order=["Endemic KS", "Epidemic KS"],
                      hue_order=["Endemic KS", "Epidemic KS"],
                      palette=PANEL_B_COLORS, legend=False,
                      size=5, alpha=0.7, linewidth=0, ax=ax)

        # Mann-Whitney U test.
        endemic = plot_df.loc[plot_df["phenotype"] == "Endemic KS", "value"].dropna()
        epidemic = plot_df.loc[plot_df["phenotype"] == "Epidemic KS", "value"].dropna()
        if len(endemic) >= 2 and len(epidemic) >= 2:
            _, p = mannwhitneyu(endemic, epidemic, alternative="two-sided")
            if p < 0.001:
                p_label = "p < 0.001"
            elif p < 0.01:
                p_label = f"p = {p:.3f}"
            else:
                p_label = f"p = {p:.2f}"
            # Draw significance bracket.
            y_max = plot_df["value"].max()
            y_bar = y_max * 1.08
            ax.plot([0, 0, 1, 1], [y_bar, y_bar * 1.02, y_bar * 1.02, y_bar],
                    color="black", lw=1.0, clip_on=False)
            ax.text(0.5, y_bar * 1.04, p_label, ha="center", va="bottom",
                    fontsize=14, fontstyle="italic")

        ax.set_title(ct, fontsize=16, fontweight="bold", pad=14)
        ax.set_ylabel("Abundance (%)" if i == 0 else "", fontsize=17)
        ax.set_xlabel("")
        ax.tick_params(axis="x", labelsize=15)
        ax.tick_params(axis="y", labelsize=15)
        for spine in ("top", "right"):
            ax.spines[spine].set_visible(False)

    # Shared y-axis limits across all panels.
    y_max_all = max(ax.get_ylim()[1] for ax in axes if ax.get_visible())
    for ax in axes:
        if ax.get_visible():
            ax.set_ylim(0, y_max_all)


def render(out_dir: Path) -> Path:
    rdf = _load_data()

    fig = plt.figure(figsize=(34, 15))
    gs = GridSpec(2, 5, height_ratios=[1, 0.7],
                  width_ratios=[1, 1, 1, 1, 0.08],
                  hspace=0.40, wspace=0.3,
                  left=0.10, right=0.93, top=0.94, bottom=0.06)

    # Panel A: heatmap spanning all 4 data columns.
    ax_a = fig.add_subplot(gs[0, :4])
    cax = fig.add_subplot(gs[0, 4])
    im = _draw_panel_a(ax_a, rdf)
    cb = fig.colorbar(im, cax=cax)
    cb.set_label("Cell type\nabundance", fontsize=16)
    cb.set_ticks([10, 20, 30, 40, 50])
    cb.ax.tick_params(labelsize=14)

    # Panel B: 4 boxplots with shared y-axis.
    axes_b = [fig.add_subplot(gs[1, i]) for i in range(4)]
    _draw_panel_b(axes_b, rdf)
    axes_b[0].text(-0.3, 1.12, "B", transform=axes_b[0].transAxes,
                   fontsize=24, fontweight="bold", va="top")

    out = Path(out_dir) / "Supplementary_Figure_02.pdf"
    fig.savefig(out)
    plt.close(fig)
    return out
