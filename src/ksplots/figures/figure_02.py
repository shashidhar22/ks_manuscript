"""Figure 2 — Cumulative TCR frequency by GLIPH-cluster category and cohort.

Derived from panel C of ``data/figures/Figure_04.pdf``, but the per-pathogen
breakdown is collapsed into three categories:

* **Clustered Knowns**   — productive CDR3 matches an entry in the known
  pathogen TCR database (CMV, EBV, HCV, HIV-1, HSV-2, Influenza A/B, MTB,
  SARS-CoV-2, etc., including multi-pathogen).
* **Clustered Unknown**  — CDR3 has no DB match but belongs to at least one
  high-confidence GLIPH cluster.
* **Unclustered Unknown** — CDR3 has no DB match and is not part of any
  GLIPH cluster.

Per repertoire we sum ``duplicate_frequency`` across all CDR3s in each
category, then plot the distribution per cohort (Endemic/Epidemic × NAT/Tumor).
"""
from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import pandas as pd
import seaborn as sns

from ..config import COHORT_COLORS
from ..io.metadata import cohort_label, load_study_metadata
from ..io.tcr import load_gliph_clusters, load_known_tcr_database, load_raw_productive

CATEGORIES = ["Clustered\nKnowns", "Clustered\nUnknowns", "Unclustered\nUnknowns"]


def _build_table() -> pd.DataFrame:
    cols = ["repertoire_id", "junction_aa", "duplicate_frequency"]
    tcr = load_raw_productive(columns=cols)
    tcr = tcr.dropna(subset=["junction_aa", "duplicate_frequency"])
    tcr = tcr[tcr["duplicate_frequency"] > 0]

    meta = load_study_metadata()
    meta = meta.assign(cohort_label=cohort_label(meta))
    cohort = meta[["trb_repertoire_id", "cohort_label"]].rename(
        columns={"trb_repertoire_id": "repertoire_id"}
    )
    cohort = cohort[cohort["cohort_label"].isin(COHORT_COLORS)]
    tcr = tcr.merge(cohort, on="repertoire_id", how="inner")

    db = load_known_tcr_database()
    known = set(db["trb_cdr3_aa"].dropna().unique())

    gliph = load_gliph_clusters().dropna(subset=["pattern", "TcRb"]).copy()
    # High-confidence GLIPH clusters only (matching notebook logic):
    # vb_score ≤ 0.01, ≥3 unique CDR3s, ≥3 unique patients, pattern length > 2,
    # pattern != "single".
    gliph["patient_id"] = gliph["Sample"].astype(str).str.extract(
        r"(008_\d+)", expand=False
    )
    pattern_stats = gliph.groupby("pattern").agg(
        n_cdr=("TcRb", "nunique"),
        n_pts=("patient_id", "nunique"),
        vb=("vb_score", "min"),
    )
    hc_patterns = pattern_stats[
        (pattern_stats["vb"] <= 0.01)
        & (pattern_stats["n_cdr"] >= 3)
        & (pattern_stats["n_pts"] >= 3)
    ].index
    hc_patterns = [
        p for p in hc_patterns
        if isinstance(p, str) and p != "single" and len(p) > 2
    ]
    clustered = set(
        gliph.loc[gliph["pattern"].isin(hc_patterns), "TcRb"].unique()
    )

    in_db = tcr["junction_aa"].isin(known)
    in_cluster = tcr["junction_aa"].isin(clustered)

    # Notebook logic: only CDR3s in HC clusters get a "Clustered" label.
    # Clustered Knowns  = in HC cluster AND in known DB
    # Clustered Unknown = in HC cluster AND NOT in known DB
    # Unclustered Unknown = not in any HC cluster (regardless of DB match)
    cat = pd.Series("Unclustered\nUnknowns", index=tcr.index)
    cat[in_cluster & ~in_db] = "Clustered\nUnknowns"
    cat[in_cluster & in_db] = "Clustered\nKnowns"
    tcr = tcr.assign(category=cat)

    summed = (
        tcr.groupby(["repertoire_id", "cohort_label", "category"], as_index=False)["duplicate_frequency"]
        .sum()
        .rename(columns={"duplicate_frequency": "frequency"})
    )
    return summed


def render(out_dir: Path) -> Path:
    df = _build_table()
    df = df[df["category"].isin(CATEGORIES)].copy()
    df["category"] = pd.Categorical(df["category"], categories=CATEGORIES, ordered=True)

    # Rename NAT → Ax. Skin for display
    df["cohort_label"] = df["cohort_label"].str.replace("NAT", "Ax. Skin", regex=False)
    renamed_colors = {k.replace("NAT", "Ax. Skin"): v for k, v in COHORT_COLORS.items()}
    palette = {k: v for k, v in renamed_colors.items() if k in df["cohort_label"].unique()}
    hue_order = [k for k in renamed_colors if k in palette]

    fig, ax = plt.subplots(figsize=(10, 6))
    # Convert frequency to percentage for plotting
    df["pct"] = df["frequency"] * 100

    sns.boxplot(data=df, x="category", y="pct", hue="cohort_label",
                hue_order=hue_order, palette=palette, fliersize=0,
                linewidth=1.5, ax=ax)
    sns.stripplot(data=df, x="category", y="pct", hue="cohort_label",
                  hue_order=hue_order, palette=palette, dodge=True,
                  size=2.8, alpha=0.6, linewidth=0, ax=ax)
    ax.set_yscale("log")
    ax.yaxis.set_major_locator(mticker.FixedLocator([0.01, 0.1, 1, 10, 100]))
    ax.yaxis.set_major_formatter(mticker.FuncFormatter(
        lambda v, _: f"{v:g}"
    ))
    ax.yaxis.set_minor_formatter(mticker.NullFormatter())
    ax.set_ylabel("Cumulative frequency (%)", fontsize=16)
    ax.set_xlabel("", fontsize=16)
    ax.tick_params(axis="both", labelsize=16)

    guide_lines = {1: "1%", 25: "25%", 75: "75%", 90: "90%"}
    for y_pct, label in guide_lines.items():
        ax.axhline(y_pct, ls="--", lw=1, color="grey")
        ax.text(ax.get_xlim()[1], y_pct, f"  {label}", va="center", ha="left",
                fontsize=12, color="grey", clip_on=False)

    handles, labels = ax.get_legend_handles_labels()
    n = len(hue_order)
    ax.legend(handles[:n], labels[:n], title="Cohort", loc="lower right",
              frameon=True, fontsize=14, title_fontsize=14)

    fig.tight_layout()
    out = Path(out_dir) / "Figure_02.pdf"
    fig.savefig(out)
    plt.close(fig)
    return out
