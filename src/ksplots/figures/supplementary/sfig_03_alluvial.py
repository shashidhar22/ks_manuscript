"""Supplementary Figure 3 — Clone tracking alluvials + COS7 reporter assays.

* **Panel A** — TRB clone tracking alluvial for subject 008_098 (HIV-, CR).
* **Panel B** — TRB clone tracking alluvial for subject 008_002 (HIV+, PR).
* **Panels C–F** — COS7 Jurkat reporter bar charts (±HLA restriction).

Data sources:
  - Panels A/B: immunoSEQ TSV files in ``data/wetlab/R Projects/ImmunoSEQ_KS/``
  - Panels C–F: ``data/wetlab/data/SpFigure3_Panel{C,D,E,F}.csv``
"""
from __future__ import annotations

import csv
import io
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.gridspec import GridSpec
from scipy.stats import ttest_ind

from ...config import WETLAB_DIR, WETLAB_TCR_COLORS

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------

_IMMUNOSEQ_DIR = WETLAB_DIR / "R Projects" / "ImmunoSEQ_KS"
_PANEL_CSV_DIR = WETLAB_DIR / "data"

# ---------------------------------------------------------------------------
# LymphoSeq2-equivalent functions
# ---------------------------------------------------------------------------


def _load_immunoseq(path: Path) -> pd.DataFrame:
    """Read a single immunoSEQ TSV file, standardising column names to MiAIRR.

    Equivalent to LymphoSeq2's ``readImmunoSeq`` for a single file.
    """
    df = pd.read_csv(path, sep="\t")
    df = df.rename(columns={
        "aminoAcid": "junction_aa",
        "count (templates)": "duplicate_count",
        "frequencyCount (%)": "duplicate_frequency_pct",
    })
    df["repertoire_id"] = path.stem
    return df


def _productive_seq(df: pd.DataFrame) -> pd.DataFrame:
    """Filter productive sequences and aggregate by junction_aa.

    Reimplements LymphoSeq2's ``productiveSeq(aggregate="junction_aa")``:

    1. Keep rows where ``junction_aa`` is non-empty and contains no stop
       codon (``*``).  This matches the R logic:
       ``aminoAcid != "" & !grepl("\\\\*", aminoAcid)``
    2. Aggregate by ``junction_aa``, summing ``duplicate_count`` across all
       nucleotide variants of the same amino-acid CDR3.
    3. Recalculate ``duplicate_frequency`` as
       ``duplicate_count / total_productive_count`` (0-1 scale).
    """
    prod = df.copy()
    # Filter: non-null, non-empty, no stop codons
    prod = prod.dropna(subset=["junction_aa"])
    prod = prod[prod["junction_aa"].astype(str).str.strip() != ""]
    prod = prod[~prod["junction_aa"].astype(str).str.contains(r"\*", regex=True)]

    # Aggregate by junction_aa — sum template counts
    agg = (
        prod.groupby("junction_aa", as_index=False)["duplicate_count"]
        .sum()
    )

    # Recalculate frequency from counts (normalised to 0-1)
    total = agg["duplicate_count"].sum()
    agg["duplicate_frequency"] = agg["duplicate_count"] / total if total > 0 else 0.0

    # Sort descending by frequency
    agg = agg.sort_values("duplicate_frequency", ascending=False).reset_index(drop=True)
    return agg


def _top_seqs(df: pd.DataFrame, top: int = 100) -> pd.DataFrame:
    """Return top *N* sequences by frequency.

    Reimplements LymphoSeq2's ``topSeqs(top=N)``.
    """
    return df.head(top).copy()


def _clone_track(
    sample_dfs: dict[str, pd.DataFrame],
    sample_list: list[str],
) -> pd.DataFrame:
    """Build a clone-tracking table across samples.

    Reimplements LymphoSeq2's ``cloneTrack``:

    1. Concatenate all sample DataFrames.
    2. Filter to *sample_list* and the union of all junction_aa present.
    3. Add a ``seen`` column = number of distinct samples each CDR3 appears in.
    """
    frames = []
    for sid in sample_list:
        df = sample_dfs[sid].copy()
        df["repertoire_id"] = sid
        frames.append(df[["repertoire_id", "junction_aa", "duplicate_frequency"]])
    combined = pd.concat(frames, ignore_index=True)

    # Count how many samples each CDR3 appears in
    seen = combined.groupby("junction_aa")["repertoire_id"].nunique().rename("seen")
    combined = combined.merge(seen, on="junction_aa")
    return combined


# ---------------------------------------------------------------------------
# Alluvial drawing (adapted from figure_03.py with multi-colour highlights)
# ---------------------------------------------------------------------------


def _draw_alluvial(
    ax: plt.Axes,
    sample_files: list[tuple[str, str, Path]],
    highlight_map: dict[str, str],
    panel_letter: str,
    annotation: str,
    top_n: int = 100,
    annotation_right: bool = False,
) -> None:
    """Draw an alluvial clone-tracking plot.

    Parameters
    ----------
    sample_files : list of (sample_id, display_label, tsv_path)
    highlight_map : CDR3 amino-acid sequence -> colour string
    """
    # --- Load, filter, aggregate per sample ---
    sample_dfs: dict[str, pd.DataFrame] = {}
    sample_order: list[str] = []
    for sid, _label, fpath in sample_files:
        raw = _load_immunoseq(fpath)
        prod = _productive_seq(raw)
        top = _top_seqs(prod, top_n)
        sample_dfs[sid] = top
        sample_order.append(sid)

    n = len(sample_order)
    bar_w = 0.22
    x_pos = np.arange(n, dtype=float)

    # --- Build stacked layouts per sample ---
    layouts: list[dict[str, tuple[float, float]]] = []
    totals: list[float] = []
    prev_order: list[str] = []

    for sid in sample_order:
        sub = sample_dfs[sid]
        freq = dict(zip(sub["junction_aa"], sub["duplicate_frequency"]))
        present = set(freq)

        # Shared clones keep their relative order from previous sample
        ordered: list[str] = [c for c in prev_order if c in present]
        # New clones sorted by frequency descending
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

    # --- Draw stacked clone bars ---
    for i, layout in enumerate(layouts):
        x = x_pos[i] - bar_w / 2
        for cdr3, (y0, y1) in layout.items():
            face = highlight_map.get(cdr3, "white")
            ax.add_patch(plt.Rectangle(
                (x, y0), bar_w, y1 - y0,
                facecolor=face, edgecolor="black", linewidth=0.35,
            ))

    # --- Draw alluvial ribbons ---
    for i in range(n - 1):
        a, b = layouts[i], layouts[i + 1]
        shared = set(a) & set(b)
        ordered = sorted(shared, key=lambda c: (a[c][0] + a[c][1]) / 2)
        x_l = x_pos[i] + bar_w / 2
        x_r = x_pos[i + 1] - bar_w / 2
        for cdr3 in ordered:
            y_l0, y_l1 = a[cdr3]
            y_r0, y_r1 = b[cdr3]
            if cdr3 in highlight_map:
                rib_color = highlight_map[cdr3]
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

    # --- Axes ---
    ax.set_xlim(-0.6, n - 0.4)
    ax.set_ylim(-0.03 * ymax, ymax * 1.02)
    ax.set_xticks(x_pos)
    ax.set_xticklabels(
        [lbl for _, lbl, _ in sample_files],
        fontsize=20, fontweight="bold",
    )
    ax.set_ylabel("TCRB sequence frequency", fontsize=24, fontweight="bold")
    ax.tick_params(axis="y", labelsize=18)
    ax.set_title(panel_letter, loc="left", fontweight="bold", fontsize=28)
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)

    # --- Annotation text ---
    if annotation_right:
        ax.text(
            0.98, 0.98, annotation, transform=ax.transAxes,
            fontsize=20, fontweight="bold", va="top", ha="right",
        )
    else:
        ax.text(
            0.02, 0.98, annotation, transform=ax.transAxes,
            fontsize=20, fontweight="bold", va="top", ha="left",
        )

    # --- Legend for highlighted clones ---
    if highlight_map:
        from matplotlib.patches import Patch
        handles = [
            Patch(facecolor=color, edgecolor="black", linewidth=0.5, label=cdr3)
            for cdr3, color in highlight_map.items()
        ]
        ncol = len(handles)  # single line
        ax.legend(
            handles=handles, fontsize=18, loc="lower center",
            bbox_to_anchor=(0.5, -0.18), ncol=ncol,
            frameon=False, title="TCRB CDR3 Sequence",
            title_fontproperties={"weight": "bold", "size": 20},
        )


# ---------------------------------------------------------------------------
# CSV bar-chart helpers (Panels C–F)
# ---------------------------------------------------------------------------

def _load_panel_csv(path: Path) -> pd.DataFrame:
    """Load a SpFigure4 panel CSV into a DataFrame with mean/SEM columns.

    CSV format: header row with condition labels (first cell empty, then
    3x negative HLA, 3x positive HLA), data rows with condition name
    followed by 6 values.
    """
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
            # Raw replicates for significance testing
            f"{neg_label}_r1": vals[0],
            f"{neg_label}_r2": vals[1],
            f"{neg_label}_r3": vals[2],
            f"{pos_label}_r1": vals[3],
            f"{pos_label}_r2": vals[4],
            f"{pos_label}_r3": vals[5],
        }
        records.append(rec)

    return pd.DataFrame(records)


_NEG_COLOR = "#999999"   # grey for (-) HLA
_OKT3_COLOR = "black"    # black for OKT3 positive control


def _draw_grouped_bars(
    ax: plt.Axes,
    df: pd.DataFrame,
    pos_color: str,
    hla_allele: str = "HLA",
    neg_color: str = _NEG_COLOR,
    ylabel: str = "mNeonGreen+ (%)",
    ylim: float = 80,
) -> None:
    """Draw paired bar chart: (-) HLA vs (+) HLA with significance brackets.

    ORF condition bars use *pos_color*; OKT3 positive-control bars use black.
    All (-) HLA bars are grey.
    """
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

        # Determine per-bar colours: OKT3 rows get black, others get pos_color
        is_okt3 = [lbl.upper().strip().startswith("OKT3") for lbl in labels]
        pos_colors = [_OKT3_COLOR if okt else pos_color for okt in is_okt3]

        # Draw (-) HLA bars — grey for all
        ax.bar(x - bar_w / 2, neg_mean, bar_w, yerr=neg_sem, capsize=3,
               color=neg_color, edgecolor="black", linewidth=0.5,
               label="(-) HLA", error_kw={"linewidth": 0.8})

        # Draw (+) HLA bars — per-bar colour
        bars_pos = ax.bar(x + bar_w / 2, pos_mean, bar_w, yerr=pos_sem,
                          capsize=3, edgecolor="black", linewidth=0.5,
                          color=pos_colors, error_kw={"linewidth": 0.8})
        # Attach a legend label only to the first non-OKT3 bar
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
    ax.set_xticklabels(labels, fontsize=18, rotation=0, ha="center")
    ax.set_ylabel(ylabel, fontsize=20, fontweight="bold")
    ax.set_ylim(0, ylim)
    ax.tick_params(labelsize=18)
    # Per-panel legend suppressed; shared legend drawn in render().
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)


# ---------------------------------------------------------------------------
# Panel configurations
# ---------------------------------------------------------------------------

_SAMPLES_098 = [
    ("008_098_A", "Ax. Skin", _IMMUNOSEQ_DIR / "008_098_A.tsv"),
    ("008_098_B", "Tumor B", _IMMUNOSEQ_DIR / "008_098_B.tsv"),
    ("008_098_C", "Tumor C", _IMMUNOSEQ_DIR / "008_098_C.tsv"),
]

_SAMPLES_002 = [
    ("008_002_A", "Ax. Skin", _IMMUNOSEQ_DIR / "008_002_A.tsv"),
    ("008_002_B", "Tumor B", _IMMUNOSEQ_DIR / "008_002_B.tsv"),
    ("008_002_C", "Tumor C", _IMMUNOSEQ_DIR / "008_002_C.tsv"),
    ("008_002_D", "Tumor D", _IMMUNOSEQ_DIR / "008_002_D.tsv"),
    ("008_002_H", "Tumor H", _IMMUNOSEQ_DIR / "008_002_H.tsv"),
    ("008_002_I", "Tumor I", _IMMUNOSEQ_DIR / "008_002_I.tsv"),
]

_HIGHLIGHT_098: dict[str, str] = {
    "CASSTGVYGYTF": "royalblue",
    "CASSIAGHEQFF": "red",
}

_HIGHLIGHT_002: dict[str, str] = {
    "CAWNLGDSNQPQHF": "purple",
}

# (csv_name, color, panel_letter, keep_rows or None, hla_allele)
_BAR_PANELS = [
    ("SpFigure3_PanelC.csv", "#D4A017", "C",       # golden yellow — B*45:01
     ["ORF6", "OKT3"], "B*45:01"),
    ("SpFigure3_PanelD.csv", "#D4A017", "D",        # golden yellow — B*45:01 (same allele)
     ["ORF6", "OKT3"], "B*45:01"),
    ("SpFigure3_PanelE.csv", "#4DAF4A", "E",        # green — A*66:01
     None, "A*66:01"),
    ("SpFigure3_PanelF.csv", "#FF7F00", "F",        # orange — B*57:03
     ["ORF59", "OKT3"], "B*57:03"),
]


# ---------------------------------------------------------------------------
# Render
# ---------------------------------------------------------------------------


def render(out_dir: Path) -> Path:
    fig = plt.figure(figsize=(24, 16))

    # Top row: alluvial plots
    gs_top = GridSpec(
        1, 2, width_ratios=[1, 2],
        left=0.05, right=0.95, top=0.95, bottom=0.55, wspace=0.25,
    )
    ax_a = fig.add_subplot(gs_top[0, 0])
    ax_b = fig.add_subplot(gs_top[0, 1])

    _draw_alluvial(
        ax_a, _SAMPLES_098, _HIGHLIGHT_098,
        panel_letter="A",
        annotation="008_098 HIV\u2013\nComplete Response",
    )
    _draw_alluvial(
        ax_b, _SAMPLES_002, _HIGHLIGHT_002,
        panel_letter="B",
        annotation="008_002 HIV+\nPartial Response",
        annotation_right=True,
    )

    # Bottom row: bar charts
    gs_bot = GridSpec(
        1, 4, width_ratios=[2, 2, 2, 2],
        left=0.05, right=0.95, top=0.42, bottom=0.10, wspace=0.45,
    )

    bar_axes = []
    for i, (csv_name, color, letter, keep, hla) in enumerate(_BAR_PANELS):
        ax = fig.add_subplot(gs_bot[0, i])
        df = _load_panel_csv(_PANEL_CSV_DIR / csv_name)
        if keep is not None:
            keep_upper = {k.upper() for k in keep}
            def _match(lbl: str) -> bool:
                s = lbl.upper().strip()
                return s in keep_upper or (
                    "OKT3" in keep_upper and s.startswith("OKT3")
                )
            df = df[df["label"].apply(_match)].reset_index(drop=True)
        _draw_grouped_bars(ax, df, color, hla_allele=hla, ylim=80)
        ax.set_title(letter, loc="left", fontweight="bold", fontsize=28)
        bar_axes.append(ax)

    # Shared legend for all bar panels — ordered: (-) HLA, (+) alleles, OKT3
    from matplotlib.patches import Patch
    raw_handles: dict[str, plt.Artist] = {}
    for ax in bar_axes:
        for h, l in zip(*ax.get_legend_handles_labels()):
            if l not in raw_handles:
                raw_handles[l] = h

    ordered_legend: list[tuple[str, plt.Artist]] = []
    # (-) HLA first
    if "(-) HLA" in raw_handles:
        ordered_legend.append(("(-) HLA", raw_handles.pop("(-) HLA")))
    # (+) HLA entries in panel order
    for _, _, _, _, hla in _BAR_PANELS:
        key = f"(+) HLA-{hla}"
        if key in raw_handles:
            ordered_legend.append((key, raw_handles.pop(key)))
    # OKT3 last
    ordered_legend.append(("OKT3", Patch(facecolor=_OKT3_COLOR,
                                         edgecolor="black", linewidth=0.5)))

    if ordered_legend:
        fig.legend(
            [h for _, h in ordered_legend],
            [l for l, _ in ordered_legend],
            loc="lower center", bbox_to_anchor=(0.5, 0.01),
            ncol=len(ordered_legend), fontsize=18, frameon=False,
        )

    out = Path(out_dir) / "Supplementary_Figure_03.pdf"
    fig.savefig(out)
    plt.close(fig)
    return out
