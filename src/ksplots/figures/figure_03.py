"""Figure 3 — Temporal clone tracking (reproduces Figure_five_2026.pdf).

Two-panel alluvial of TCR clone frequencies across time-course samples for
two patients, each panel topped by a schematic body showing where the serial
biopsies were taken (one coloured dot per sample, matched to the coloured dot
under each sample's bar):

* Panel A — Patient 008_141 (Endemic KS, blue): NAT skin V01 -> tumors V01,
  V05 (x3), V08 (x2).
* Panel B — Patient 008_008 (Epidemic KS, red): tumors V01 (x2), V05 (x3),
  V08 (x2).

Each sample is rendered as a vertical bar of stacked clones (sorted by
frequency, top clones largest at the bottom). Coloured ribbons connect the
same CDR3 across consecutive samples; clones present in a high-confidence
GLIPH "Clustered Unknown" cluster are highlighted in the cohort colour, and
all other clones are shown in pale grey.

The body schematic is a generic, publication-licensable human-body outline,
committed as a packaged asset (``src/ksplots/assets/fig_03_body.png``) and
drawn edge-to-edge above each panel; the biopsy-site dots are drawn on a
transparent full-figure overlay so they sit on top of the body. This replaces
the figure's former Illustrator post-processing, so Figure 3 is now fully
reproducible from the repo.
"""
from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from ..io.tcr import load_gliph_clusters, load_known_tcr_database, load_raw_productive

# Packaged body-outline asset (RGBA PNG with a transparent background).
_BODY_ASSET = Path(__file__).resolve().parents[1] / "assets" / "fig_03_body.png"

# Patient -> cohort colour + ordered samples. Each sample is
#   (repertoire_id, display label, dot colour, (body_frac_x, body_frac_y))
# where the body fractions place that sample's biopsy-site dot on the body
# outline (x,y as a fraction of the body image box, y measured from the top).
PATIENTS: dict[str, dict] = {
    "008_141": {
        "color": "#377eb8",  # endemic blue
        "samples": [
            ("008_141_A", "Ax. Skin\nWeek 0",  (0.353, 0.667, 0.275), (0.373, 0.294)),
            ("008_141_D", "Tumor D\nWeek 0",   (0.957, 0.635, 0.380), (0.338, 0.594)),
            ("008_141_E", "Tumor E\nWeek 12",  (0.902, 0.224, 0.275), (0.384, 0.638)),
            ("008_141_F", "Tumor F\nWeek 12",  (0.710, 0.090, 0.620), (0.384, 0.535)),
            ("008_141_G", "Tumor G\nWeek 12",  (0.608, 0.365, 0.898), (0.384, 0.711)),
            ("008_141_H", "Tumor H\nWeek 20",  (0.263, 0.380, 0.933), (0.451, 0.587)),
            ("008_141_I", "Tumor I\nWeek 20",  (0.106, 0.122, 0.541), (0.384, 0.783)),
        ],
    },
    "008_008": {
        "color": "#e41a1c",  # epidemic red
        "samples": [
            ("008_008_C", "Tumor C\nWeek 0",   (0.831, 0.722, 0.173), (0.613, 0.507)),
            ("008_008_D", "Tumor D\nWeek 0",   (0.957, 0.635, 0.380), (0.334, 0.559)),
            ("008_008_E", "Tumor E\nWeek 12",  (0.902, 0.224, 0.275), (0.613, 0.766)),
            ("008_008_F", "Tumor F\nWeek 12",  (0.710, 0.090, 0.620), (0.368, 0.606)),
            ("008_008_G", "Tumor G\nWeek 12",  (0.608, 0.365, 0.898), (0.369, 0.507)),
            ("008_008_H", "Tumor H\nWeek 22",  (0.263, 0.380, 0.933), (0.613, 0.703)),
            ("008_008_I", "Tumor I\nWeek 22",  (0.106, 0.122, 0.541), (0.648, 0.605)),
        ],
    },
}

TOP_PER_SAMPLE = 500

# --- Composite layout (figure-fraction coordinates) --------------------------
_FIG_W, _FIG_H = 12.0, 10.4
_PANEL_RECT = {                       # [x0, y0, w, h] of each alluvial axes
    "008_141": [0.085, 0.110, 0.400, 0.580],
    "008_008": [0.565, 0.110, 0.400, 0.580],
}
_BODY_BY0 = 0.700                     # bottom of the body images
_BODY_BH = 0.270                      # height of the body images
_BODY_DOT_S = 120                     # biopsy-site dot area (points^2)
_SAMPLE_DOT_S = 70                    # per-bar colour-key dot area (points^2)
_XTICK_PAD = 24                       # push x labels below the colour-key dots
_PANEL_LETTER_SIZE = 34
_PANEL_LETTER_Y = 0.985
_PANEL_LETTER_X = {"008_141": 0.020, "008_008": 0.500}


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
    for rep, *_ in samples:
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

    # 4) Colour-key dot under each bar (matched to the body-site dot above).
    for i, s in enumerate(samples):
        ax.scatter([x_pos[i]], [-0.045 * ymax], s=_SAMPLE_DOT_S, c=[s[2]],
                   edgecolors="none", clip_on=False, zorder=5)

    ax.set_xticks(x_pos)
    ax.set_xticklabels([s[1] for s in samples], fontsize=12, fontweight="bold")
    ax.set_xlabel("Sample", fontsize=14)
    ax.set_ylabel("Frequency of CDR3 sequence", fontsize=14)
    ax.tick_params(axis="x", pad=_XTICK_PAD)
    ax.tick_params(axis="y", labelsize=11)
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)


def render(out_dir: Path) -> Path:
    needed_reps = [s[0] for p in PATIENTS.values() for s in p["samples"]]
    tcr = load_raw_productive(columns=["repertoire_id", "junction_aa", "duplicate_frequency"])
    tcr = tcr[tcr["repertoire_id"].isin(needed_reps)].dropna(
        subset=["junction_aa", "duplicate_frequency"]
    )

    highlight = _highlighted_cdr3s()
    body_img = plt.imread(str(_BODY_ASSET))

    fig = plt.figure(figsize=(_FIG_W, _FIG_H))

    # Alluvial panels.
    for key, rect in _PANEL_RECT.items():
        ax = fig.add_axes(rect)
        _draw_panel(ax, PATIENTS[key], tcr, highlight)

    # Body outline centred above each panel (aspect preserved, no distortion).
    asp = body_img.shape[1] / body_img.shape[0]          # width / height
    body_bw = _BODY_BH * asp * _FIG_H / _FIG_W           # figure-fraction width
    body_rects: dict[str, list[float]] = {}
    for key, rect in _PANEL_RECT.items():
        cx = rect[0] + rect[2] / 2
        br = [cx - body_bw / 2, _BODY_BY0, body_bw, _BODY_BH]
        body_rects[key] = br
        bax = fig.add_axes(br)
        bax.imshow(body_img)
        bax.axis("off")

    # Biopsy-site dots on a transparent overlay so they sit on top of the body.
    overlay = fig.add_axes([0, 0, 1, 1])
    overlay.set_xlim(0, 1)
    overlay.set_ylim(0, 1)
    overlay.axis("off")
    overlay.patch.set_alpha(0)
    overlay.set_zorder(20)
    for key, br in body_rects.items():
        for s in PATIENTS[key]["samples"]:
            fx, fy = s[3]
            dx = br[0] + fx * br[2]
            dy = br[1] + (1.0 - fy) * br[3]             # fy is from the top
            overlay.scatter([dx], [dy], s=_BODY_DOT_S, c=[s[2]],
                            edgecolors="none", zorder=21)

    # Panel letters.
    for key, letter in (("008_141", "A"), ("008_008", "B")):
        fig.text(_PANEL_LETTER_X[key], _PANEL_LETTER_Y, letter,
                 fontsize=_PANEL_LETTER_SIZE, fontweight="bold",
                 va="top", ha="left")

    out = Path(out_dir) / "Figure_03.pdf"
    fig.savefig(out)
    plt.close(fig)
    return out
