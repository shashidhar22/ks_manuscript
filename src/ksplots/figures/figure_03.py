"""Figure 3 — Temporal clone tracking (reproduces Figure_five_2026.pdf).

Two-panel alluvial of TCR clone frequencies across time-course samples for
two patients:

* Panel A — Patient 008_141 (Endemic KS, blue): NAT skin V01 → tumors V01,
  V05 (×3), V08 (×2).
* Panel B — Patient 008_008 (Epidemic KS, red): tumors V01 (×2), V05 (×3),
  V08 (×2).

Each sample is rendered as a vertical bar of stacked clones (sorted by
frequency, top clones largest at the bottom). Coloured ribbons connect the
same CDR3 across consecutive samples; clones present in a high-confidence
GLIPH "Clustered Unknown" cluster are highlighted in the cohort colour, and
all other clones are shown in pale grey.
"""
from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from ..io.tcr import load_gliph_clusters, load_known_tcr_database, load_raw_productive

# Patient → ordered (repertoire_id, display label)
PATIENTS: dict[str, dict] = {
    "008_141": {
        "color": "#377eb8",  # endemic blue
        "title": "A",
        "samples": [
            ("008_141_A", "Ax. Skin\nWeek 0"),
            ("008_141_D", "Tumor D\nWeek 0"),
            ("008_141_E", "Tumor E\nWeek 12"),
            ("008_141_F", "Tumor F\nWeek 12"),
            ("008_141_G", "Tumor G\nWeek 12"),
            ("008_141_H", "Tumor H\nWeek 20"),
            ("008_141_I", "Tumor I\nWeek 20"),
        ],
    },
    "008_008": {
        "color": "#e41a1c",  # epidemic red
        "title": "B",
        "samples": [
            ("008_008_C", "Tumor C\nWeek 0"),
            ("008_008_D", "Tumor D\nWeek 0"),
            ("008_008_E", "Tumor E\nWeek 12"),
            ("008_008_F", "Tumor F\nWeek 12"),
            ("008_008_G", "Tumor G\nWeek 12"),
            ("008_008_H", "Tumor H\nWeek 22"),
            ("008_008_I", "Tumor I\nWeek 22"),
        ],
    },
}

TOP_PER_SAMPLE = 500


def _highlighted_cdr3s() -> set[str]:
    """Reproduce ``kstme_clustered_unknown_table$junction_aa`` from the notebook.

    A CDR3 is "Clustered Unknown" iff:

    1. It belongs to a *high-confidence* GLIPH pattern, where the pattern
       satisfies all of:
         * ``vb_score <= 0.01``
         * ``number_unique_cdr3 >= 3``
         * unique patients (``008_NNN`` extracted from ``Sample``) ``>= 3``
         * ``len(pattern) > 2`` and pattern != "single"
    2. It does *not* appear in the known pathogen TCR database (otherwise it
       would be labelled with that pathogen).
    """
    gliph = load_gliph_clusters().dropna(subset=["pattern", "TcRb"]).copy()

    # Extract patient_id (e.g. 008_141) from "008_141_B:Endemic KS - Tumor".
    sample = gliph["Sample"].astype(str)
    gliph["patient_id"] = sample.str.extract(r"(008_\d+)", expand=False)

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
    hc_patterns = [p for p in hc_patterns if isinstance(p, str)
                   and p != "single" and len(p) > 2]

    known = set(load_known_tcr_database()["trb_cdr3_aa"].dropna().unique())

    hc = gliph[gliph["pattern"].isin(hc_patterns)]
    # Drop any pattern that contains *any* CDR3 in the known pathogen DB —
    # such a pattern is annotated with that pathology, not "Clustered Unknown".
    pattern_has_known = hc.groupby("pattern")["TcRb"].apply(
        lambda s: bool(set(s) & known)
    )
    unknown_patterns = pattern_has_known[~pattern_has_known].index
    hc_cdr3s = set(hc.loc[hc["pattern"].isin(unknown_patterns), "TcRb"].unique())
    return hc_cdr3s - known


def _per_sample(tcr: pd.DataFrame, repertoire_id: str) -> pd.DataFrame:
    sub = tcr[tcr["repertoire_id"] == repertoire_id]
    sub = (
        sub.groupby("junction_aa", as_index=False)["duplicate_frequency"]
        .sum()
        .sort_values("duplicate_frequency", ascending=False)
        .head(TOP_PER_SAMPLE)
        .reset_index(drop=True)
    )
    return sub


def _draw_panel(ax, panel: dict, tcr: pd.DataFrame, highlight_set: set[str]) -> None:
    samples = panel["samples"]
    color = panel["color"]
    n = len(samples)

    bar_w = 0.22
    x_pos = np.arange(n, dtype=float)

    # 1) For each sample build (cdr3 -> (y_low, y_high)) for stacked clones.
    #    To minimise ribbon crossings, carry forward the stacking order from
    #    the previous sample: shared clones keep their relative order, then
    #    new clones are appended sorted by frequency (largest first → top).
    layouts: list[dict[str, tuple[float, float]]] = []
    totals: list[float] = []
    prev_order: list[str] = []
    for rep, _ in samples:
        sub = _per_sample(tcr, rep)
        freq = dict(zip(sub["junction_aa"], sub["duplicate_frequency"]))
        present = set(freq)
        # Clones shared with previous sample keep their relative order.
        ordered: list[str] = [c for c in prev_order if c in present]
        # New clones (not in previous sample) sorted by frequency desc.
        new = sorted(present - set(ordered), key=lambda c: freq[c], reverse=True)
        ordered.extend(new)
        layout: dict[str, tuple[float, float]] = {}
        y = 0.0
        for cdr3 in ordered:
            h = float(freq[cdr3])
            layout[cdr3] = (y, y + h)
            y += h
        layouts.append(layout)
        totals.append(y)
        prev_order = ordered

    ymax = max(totals) if totals else 1.0

    # 2) Draw clone bars (stacked rects) — highlight clones in cohort colour.
    # Each CDR3 gets a visible black border to recreate the LymphoSeq2
    # "ladder of rungs" look from the reference figure.
    for i, layout in enumerate(layouts):
        x = x_pos[i] - bar_w / 2
        for cdr3, (y0, y1) in layout.items():
            face = color if cdr3 in highlight_set else "white"
            ax.add_patch(plt.Rectangle((x, y0), bar_w, y1 - y0,
                                        facecolor=face, edgecolor="black",
                                        linewidth=0.35))

    # 3) Draw alluvial ribbons between consecutive samples.
    for i in range(n - 1):
        a, b = layouts[i], layouts[i + 1]
        shared = set(a) & set(b)
        # Stable order so ribbons cross cleanly: by mid-y in left bar.
        ordered = sorted(shared, key=lambda c: (a[c][0] + a[c][1]) / 2)
        x_l = x_pos[i] + bar_w / 2
        x_r = x_pos[i + 1] - bar_w / 2
        for cdr3 in ordered:
            y_l0, y_l1 = a[cdr3]
            y_r0, y_r1 = b[cdr3]
            if cdr3 in highlight_set:
                rib_color = color
                alpha = 0.7
            else:
                rib_color = "#cccccc"
                alpha = 0.45
            xs = np.linspace(x_l, x_r, 30)
            t = (xs - x_l) / (x_r - x_l)
            smooth = 3 * t**2 - 2 * t**3
            y_lower = y_l0 + (y_r0 - y_l0) * smooth
            y_upper = y_l1 + (y_r1 - y_l1) * smooth
            ax.fill_between(xs, y_lower, y_upper, color=rib_color, alpha=alpha,
                            linewidth=0)

    ax.set_xlim(-0.6, n - 0.4)
    ax.set_ylim(-0.03 * ymax, ymax * 1.02)
    ax.set_xticks(x_pos)
    ax.set_xticklabels([s[1] for s in samples], fontsize=12, fontweight="bold")
    ax.set_xlabel("Sample", fontsize=14)
    ax.set_ylabel("Frequency of CDR3 sequence", fontsize=14)
    ax.tick_params(axis="y", labelsize=11)
    ax.set_title(panel["title"], loc="left", fontweight="bold", fontsize=20)
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)


def render(out_dir: Path) -> Path:
    needed_reps = [s for p in PATIENTS.values() for s, _ in p["samples"]]
    tcr = load_raw_productive(columns=["repertoire_id", "junction_aa", "duplicate_frequency"])
    tcr = tcr[tcr["repertoire_id"].isin(needed_reps)].dropna(
        subset=["junction_aa", "duplicate_frequency"]
    )

    highlight = _highlighted_cdr3s()

    fig, axes = plt.subplots(1, 2, figsize=(12, 8))
    _draw_panel(axes[0], PATIENTS["008_141"], tcr, highlight)
    _draw_panel(axes[1], PATIENTS["008_008"], tcr, highlight)
    fig.tight_layout(w_pad=3)
    out = Path(out_dir) / "Figure_03.pdf"
    fig.savefig(out)
    plt.close(fig)
    return out
