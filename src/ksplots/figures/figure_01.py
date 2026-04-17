"""Figure 1 — KSHV and HIV viral expression heatmaps.

Reproduces ``reference/Figure_one_2026.pdf`` from the nf-core/rnaseq viral
re-run output:

    data/rnaseq/results_viral/heatmaps/{kshv,hiv1}_salmon_tpm_heatmap.tsv

Layout
------
* Two panels side by side: **A. KSHV** (large) and **B. HIV** (narrow).
* **Rows = samples**, grouped top → bottom by 9 cohorts:
  Uganda Endemic / Epidemic, Tanzania 1 Endemic ± LA-Skin, Tanzania 1
  Epidemic ± LA-Skin, Tanzania 1 Non-KS Control, Tanzania 2 Epidemic ± LA-Skin.
  Each cohort gets a black-bordered name box on the far left and per-sample
  text labels next to the heatmap.
* **Columns = viral genes**, ordered by lytic stage. Stages are wrapped in
  black-bordered banner boxes above each gene group; gene names sit at the
  bottom rotated 90°.
* Two stacked colorbars (KSHV log2(TPM+1), HIV log2(TPM+1)) on the right.
"""
from __future__ import annotations

import re
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.gridspec import GridSpec
from matplotlib.patches import Rectangle

from ..io.viral import load_kshv_gene_map, load_viral_tpm

# (sample-name prefix, display label) — order matches Figure_one_2026.pdf.
COHORTS: list[tuple[str, str]] = [
    ("endemic_008",                   "Uganda\nEndemic KS Lesions"),
    ("epidemic_008",                  "Uganda\nEpidemic KS Lesions"),
    ("lidenge_lesion_endemic",        "Tanzania 1\nEndemic KS Lesions"),
    ("lidenge_control_skin_endemic",  "Tanzania 1 Endemic\nAxillary Skin"),
    ("lidenge_lesion_epidemic",       "Tanzania 1\nEpidemic KS Lesions"),
    ("lidenge_control_skin_epidemic", "Tanzania 1 Epidemic\nAxillary Skin"),
    ("lidenge_normal_skin",           "Tanzania 1 Non-KS\nControl Skin"),
    ("tso_lesion",                    "Tanzania 2\nEpidemic KS Lesions"),
    ("tso_control",                   "Tanzania 2 Epidemic\nAxillary Skin"),
]

# Stage order matching Figure_one_2026.pdf.
KSHV_STAGE_ORDER = ["La", "IE", "E1", "E2", "E3", "L4", "L5", "ncRNA"]
HIV_STAGE_ORDER = ["HIV"]

# Authoritative gene → stage mapping from the manuscript notebook, sourced from:
#   https://journals.plos.org/plospathogens/article?id=10.1371/journal.ppat.1001013
#   https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3894221/
# Display-name keys match _gene_label() output.
KSHV_GENE_STAGE: dict[str, str] = {
    "vIRF-3": "La", "ORF73": "La", "ORF71": "La", "ORF72": "La", "K12": "La",
    "K4.2": "IE", "ORF45": "IE", "ORF48": "IE", "ORF50": "IE", "K8": "IE",
    "ORF70": "E1", "ORF10": "E1", "ORF56": "E1", "ORF11": "E1", "K3": "E1",
    "ORF40": "E1", "ORF59": "E1", "K15": "E1", "ORF54": "E1", "K4": "E1",
    "K4.1": "E1", "K15.1": "E1",
    "ORF2": "E2", "K1": "E2", "K2": "E2", "ORF46": "E2", "ORF9": "E2",
    "K5": "E2", "K7": "E2", "K6": "E2", "vIRF-1": "E2", "K14": "E2",
    "vIRF-2": "E2", "ORF49": "E2", "ORF74": "E2", "ORF7": "E2", "ORF6": "E2",
    "ORF17.5": "E2", "K14.1": "E2",
    "ORF44": "E3", "ORF31": "E3", "ORF19": "E3", "ORF37": "E3", "ORF57": "E3",
    "ORF17": "E3", "ORF36": "E3", "ORF29": "E3", "ORF21": "E3", "ORF61": "E3",
    "ORF60": "E3", "ORF66": "E3", "ORF16": "E3", "vIRF-4": "E3", "ORF69": "E3",
    "ORF23": "L4", "ORF22": "L4", "ORF25": "L4", "ORF20": "L4", "ORF26": "L4",
    "ORF24": "L4", "ORF27": "L4", "ORF30": "L4", "ORF28": "L4", "ORF34": "L4",
    "ORF18": "L4", "ORF35": "L4", "ORF63": "L4", "ORF62": "L4", "ORF65": "L4",
    "ORF32": "L4", "ORF38": "L4", "ORF33": "L4", "ORF43": "L4", "ORF64": "L4",
    "ORF67": "L4", "ORF68": "L4",
    "ORF42": "L5", "ORF4": "L5", "ORF39": "L5", "K8.1": "L5", "ORF8": "L5",
    "ORF47": "L5", "ORF75": "L5", "ORF58": "L5", "ORF52": "L5", "ORF53": "L5",
    "ORF55": "L5", "ORF67A": "L5",
    "T1.1": "ncRNA",
    "miR-K1": "ncRNA", "miR-K2": "ncRNA", "miR-K3": "ncRNA", "miR-K4": "ncRNA",
    "miR-K5": "ncRNA", "miR-K6": "ncRNA", "miR-K7": "ncRNA", "miR-K8": "ncRNA",
    "miR-K9": "ncRNA",
}

# Exact display order from Figure_one_2026.pdf (left → right).
KSHV_GENE_ORDER: list[str] = [
    # La
    "vIRF-3", "ORF71", "ORF73", "ORF72", "K12",
    # IE
    "K4.2", "ORF45", "ORF48", "ORF50", "K8",
    # E1
    "ORF70", "ORF10", "ORF56", "ORF11", "K3", "ORF40", "ORF59", "K15",
    "K15.1", "ORF54", "K4", "K4.1",
    # E2
    "ORF2", "K1", "K2", "ORF46", "ORF9", "K5", "K7", "K6",
    "vIRF-1", "vIRF-2", "ORF49", "ORF74", "ORF7", "ORF6", "ORF17.5", "K14", "K14.1",
    # E3
    "ORF44", "ORF31", "ORF19", "ORF37", "ORF57",
    "ORF17", "ORF36", "ORF29", "ORF21", "ORF61", "ORF60", "ORF66", "ORF16",
    "vIRF-4", "ORF69",
    "ORF23", "ORF22", "ORF25", "ORF20", "ORF26", "ORF24", "ORF27", "ORF30",
    "ORF28", "ORF34", "ORF18", "ORF35", "ORF63", "ORF62", "ORF65", "ORF32",
    "ORF38", "ORF33", "ORF43", "ORF64", "ORF67", "ORF68",
    # L4  (reference figure continues without a visible break from E3→L4 here;
    #      the stage banner covers the right set)
    # — see note below: the reference banner boundaries are extracted separately
    # L5
    "ORF75", "ORF42",
    "ORF4", "ORF39", "K8.1", "ORF8", "ORF47", "ORF58", "ORF52", "ORF53",
    "ORF55", "ORF67A",
    # ncRNA
    "T1.1",
    "miR-K9", "miR-K8", "miR-K7", "miR-K6", "miR-K5", "miR-K4",
    "miR-K3", "miR-K2", "miR-K1",
]

# Override table for the 19 genes the salmon TSV flags "Unclassified".
# Maps gene_id (from the TSV index) → (display name, stage).
ORPHAN_OVERRIDES: dict[str, tuple[str, str]] = {
    "KSHV_HHV8GK18_gp77": ("ORF68",  "L4"),
    "KSHV_HHV8GK18_gp78": ("ORF69",  "E3"),
    "KSHV_HHV8GK18_gp79": ("K12",    "La"),
    "KSHV_HHV8GK18_gp80": ("ORF71",  "La"),
    "KSHV_HHV8GK18_gp81": ("ORF73",  "La"),
    "KSHV_HHV8GK18_gp82": ("ORF72",  "La"),
    "KSHV_HHV8GK18_gp83": ("K14",    "E2"),
    "KSHV_HHV8GK18_gp84": ("K14.1",  "E2"),
    "KSHV_HHV8GK18_gp86": ("K15.1",  "E1"),
    "KSHV_HHV8_gs01":     ("T1.1",   "ncRNA"),
    "KSHV_HHV8_gs02":     ("miR-K1", "ncRNA"),
    "KSHV_HHV8_gs03":     ("miR-K2", "ncRNA"),
    "KSHV_HHV8_gs04":     ("miR-K3", "ncRNA"),
    "KSHV_HHV8_gs05":     ("miR-K4", "ncRNA"),
    "KSHV_HHV8_gs06":     ("miR-K5", "ncRNA"),
    "KSHV_HHV8_gs07":     ("miR-K6", "ncRNA"),
    "KSHV_HHV8_gs08":     ("miR-K7", "ncRNA"),
    "KSHV_HHV8_gs09":     ("miR-K8", "ncRNA"),
    "KSHV_HHV8_gs10":     ("miR-K9", "ncRNA"),
}

# Stage corrections for NAMED genes whose stage in the source TSV doesn't
# match the notebook's kshv_groups mapping.
NAMED_STAGE_OVERRIDES: dict[str, str] = {
    # No overrides needed — KSHV_GENE_STAGE is authoritative
}


def _classify(sample: str) -> tuple[int, str] | None:
    for i, (prefix, label) in enumerate(COHORTS):
        if sample.startswith(prefix):
            return i, label
    return None


def _sample_label(sample: str) -> str:
    """Convert a column name to the per-row text label shown in the figure."""
    if sample.startswith("endemic_008") or sample.startswith("epidemic_008"):
        # endemic_008-061-C  → 008_061 C
        # epidemic_008-001-B-12 → 008_001 B   (drop trailing batch number)
        body = sample.split("_", 1)[1]               # 008-061-C  /  008-001-B-12
        parts = body.split("-")
        if len(parts) >= 3:
            return f"008_{parts[1]} {parts[2]}"
        return body.replace("-", "_")
    if sample.startswith("lidenge_lesion_endemic_") or sample.startswith("lidenge_control_skin_endemic_"):
        return sample.rsplit("_", 1)[-1].upper()
    if sample.startswith("lidenge_lesion_epidemic_") or sample.startswith("lidenge_control_skin_epidemic_"):
        return sample.rsplit("_", 1)[-1].upper()
    if sample.startswith("lidenge_normal_skin_"):
        return sample.rsplit("_", 1)[-1].upper()
    if sample.startswith("tso_lesion_") or sample.startswith("tso_control_"):
        return sample.rsplit("_", 1)[-1].upper()
    return sample


def _ordered_samples(columns: list[str]) -> tuple[list[str], list[tuple[str, int, int]]]:
    rows = []
    for c in columns:
        hit = _classify(c)
        if hit is not None:
            rows.append((hit[0], hit[1], c))
    rows.sort(key=lambda r: (r[0], r[2]))
    ordered = [r[2] for r in rows]
    spans: list[tuple[str, int, int]] = []
    if not rows:
        return ordered, spans
    cur, start = rows[0][1], 0
    for i, r in enumerate(rows[1:], start=1):
        if r[1] != cur:
            spans.append((cur, start, i))
            cur, start = r[1], i
    spans.append((cur, start, len(rows)))
    return ordered, spans


def _normalise_stage(value, gene_name: str, display_label: str) -> str:
    """Map salmon-table lytic_stage to a figure stage.

    Priority: ORPHAN_OVERRIDES → KSHV_GENE_STAGE (by display name) → TSV value.
    """
    override = ORPHAN_OVERRIDES.get(gene_name)
    if override is not None:
        return override[1]
    stage = KSHV_GENE_STAGE.get(display_label)
    if stage is not None:
        return stage
    if pd.isna(value) or value == "" or value == "Unclassified":
        return "ncRNA"
    return str(value)


def _gene_label(gene: str, name_lookup: dict[str, str]) -> str:
    """Display label.

    Priority: ORPHAN_OVERRIDES → GFF ``gene=`` attribute → locus_tag with
    KSHV_/HIV1_ prefix stripped.
    """
    override = ORPHAN_OVERRIDES.get(gene)
    if override is not None:
        return override[0]  # display name
    nice = name_lookup.get(gene)
    if nice and not nice.startswith("HHV8"):
        return nice
    return re.sub(r"^(KSHV|HIV1)_", "", gene)


def _gene_order(
    stage: pd.Series,
    stage_order: list[str],
    label_to_gene: dict[str, str],
    explicit_order: list[str] | None = None,
) -> tuple[list[str], list[tuple[str, int, int]]]:
    """Return (ordered gene_ids, stage spans).

    If *explicit_order* is given (list of display labels), genes are placed in
    that exact order.  Genes in the data but missing from the list are appended
    at the end, sorted by stage then name.
    """
    gene_ids = list(stage.index)

    if explicit_order is not None:
        # Build label→gene_id reverse map
        label2id = {v: k for k, v in label_to_gene.items()}
        ordered_genes: list[str] = []
        seen: set[str] = set()
        for label in explicit_order:
            gid = label2id.get(label)
            if gid is not None and gid in stage.index and gid not in seen:
                ordered_genes.append(gid)
                seen.add(gid)
        # Append any remaining genes not in the explicit order
        remaining = [g for g in gene_ids if g not in seen]
        remaining.sort(key=lambda g: (
            stage_order.index(stage.loc[g]) if stage.loc[g] in stage_order else 999,
            g,
        ))
        ordered_genes.extend(remaining)
    else:
        cats = [s for s in stage_order if s in set(stage)] + sorted(set(stage) - set(stage_order))
        df = pd.DataFrame({
            "gene": stage.index,
            "stage": pd.Categorical(stage.values, categories=cats, ordered=True),
        })
        df = df.sort_values(["stage", "gene"], kind="stable")
        ordered_genes = df["gene"].tolist()

    # Build stage spans
    spans: list[tuple[str, int, int]] = []
    cur, start = None, 0
    for i, g in enumerate(ordered_genes):
        s = stage.loc[g]
        if s != cur:
            if cur is not None:
                spans.append((cur, start, i))
            cur, start = s, i
    if cur is not None:
        spans.append((cur, start, len(ordered_genes)))
    return ordered_genes, spans


def _prepare(
    df: pd.DataFrame,
    stage_order: list[str],
    default_stage: str,
    name_lookup: dict[str, str],
    explicit_order: list[str] | None = None,
) -> tuple[pd.DataFrame, list, list]:
    df = df.copy()
    if "lytic_stage" in df.columns:
        raw_stage = df["lytic_stage"]
        df = df.drop(columns="lytic_stage")
    else:
        raw_stage = pd.Series([default_stage] * len(df), index=df.index)

    # Build gene_id → display label map for every gene in the data
    label_map = {g: _gene_label(g, name_lookup) for g in df.index}

    stage = pd.Series(
        [_normalise_stage(v, g, label_map[g]) for g, v in zip(df.index, raw_stage)],
        index=df.index,
    )
    if default_stage == "HIV":
        stage[:] = "HIV"

    sample_cols, row_spans = _ordered_samples(list(df.columns))
    gene_rows, col_spans = _gene_order(
        stage, stage_order, label_to_gene=label_map, explicit_order=explicit_order,
    )
    mat = df.loc[gene_rows, sample_cols].T  # rows=samples, cols=genes
    return mat, row_spans, col_spans


def _draw_boxed_label(ax, x: float, y: float, w: float, h: float, text: str,
                      fontsize: int = 7, fontweight: str = "bold") -> None:
    """Draw a black-bordered rectangle with centred text inside (axes coords)."""
    ax.add_patch(Rectangle((x, y), w, h, transform=ax.transAxes,
                           facecolor="white", edgecolor="black",
                           linewidth=1.0, clip_on=False, zorder=4))
    ax.text(x + w / 2, y + h / 2, text, transform=ax.transAxes,
            ha="center", va="center", fontsize=fontsize,
            fontweight=fontweight, zorder=5)


def _draw_panel(fig, ax, mat: pd.DataFrame, row_spans, col_spans, vmax: float,
                panel_title: str, show_row_labels: bool, cohort_box_width: float,
                name_lookup: dict[str, str]):
    arr = np.log2(mat.to_numpy(dtype=float) + 1.0)
    im = ax.imshow(arr, aspect="auto", cmap="viridis", interpolation="nearest",
                   vmin=0, vmax=vmax)

    # Gene names along the bottom.
    gene_labels = [_gene_label(g, name_lookup) for g in mat.columns]
    ax.set_xticks(np.arange(len(gene_labels)))
    ax.set_xticklabels(gene_labels, fontsize=15, rotation=90)
    ax.tick_params(axis="x", length=2, pad=2)

    # Sample names on the inside-left of the matrix (only on KSHV panel).
    if show_row_labels:
        sample_labels = [_sample_label(s) for s in mat.index]
        ax.set_yticks(np.arange(len(sample_labels)))
        ax.set_yticklabels(sample_labels, fontsize=13)
        ax.tick_params(axis="y", length=0, pad=4)
        for lbl in ax.get_yticklabels():
            lbl.set_clip_on(False)
    else:
        ax.set_yticks([])

    # White dividers between row groups (cohorts) and column groups (stages).
    for _, start, end in row_spans:
        if start > 0:
            ax.axhline(start - 0.5, color="white", lw=1.4)
    for _, start, end in col_spans:
        if start > 0:
            ax.axvline(start - 0.5, color="white", lw=1.4)

    # Cohort boxes on the far left (only on KSHV panel).
    if show_row_labels:
        n_rows = len(mat.index)
        for label, start, end in row_spans:
            y0 = 1.0 - end / n_rows
            y1 = 1.0 - start / n_rows
            _draw_boxed_label(
                ax,
                x=-cohort_box_width - 0.07,
                y=y0,
                w=cohort_box_width,
                h=y1 - y0,
                text=f"{label}\n(n = {end - start})",
                fontsize=20,
            )

    # Lytic-stage boxes above each column group.
    n_cols = len(mat.columns)
    for label, start, end in col_spans:
        x0 = start / n_cols
        x1 = end / n_cols
        _draw_boxed_label(
            ax,
            x=x0,
            y=1.005,
            w=x1 - x0,
            h=0.018,
            text=label,
            fontsize=20,
        )

    ax.set_title(panel_title, loc="left", fontweight="bold", fontsize=26, pad=58)
    return im


def render(out_dir: Path) -> Path:
    kshv = load_viral_tpm("kshv")
    hiv = load_viral_tpm("hiv1")

    gene_map = load_kshv_gene_map()
    name_lookup = gene_map["name"].to_dict()

    kshv_mat, kshv_rows, kshv_cols = _prepare(
        kshv, KSHV_STAGE_ORDER, default_stage="ncRNA",
        name_lookup=name_lookup, explicit_order=KSHV_GENE_ORDER,
    )
    hiv_mat, hiv_rows, hiv_cols = _prepare(
        hiv, HIV_STAGE_ORDER, default_stage="HIV",
        name_lookup=name_lookup,
    )

    n_kshv_genes = kshv_mat.shape[1]
    n_hiv_genes = hiv_mat.shape[1]

    # Portrait orientation matching reference/Figure_one_2026.pdf (~32 × 40 in).
    fig = plt.figure(figsize=(34, 32))
    gs = GridSpec(
        2, 3,
        width_ratios=[n_kshv_genes, max(n_hiv_genes, 6), 1.0],
        height_ratios=[1, 1],
        wspace=0.04, hspace=0.6, figure=fig,
        left=0.24, right=0.96, top=0.94, bottom=0.06,
    )
    ax_k = fig.add_subplot(gs[:, 0])
    ax_h = fig.add_subplot(gs[:, 1])
    cax_k = fig.add_subplot(gs[0, 2])
    cax_h = fig.add_subplot(gs[1, 2])

    im_k = _draw_panel(fig, ax_k, kshv_mat, kshv_rows, kshv_cols, vmax=5.0,
                       panel_title="A.  KSHV gene\n      expression",
                       show_row_labels=True, cohort_box_width=0.20,
                       name_lookup=name_lookup)
    im_h = _draw_panel(fig, ax_h, hiv_mat, hiv_rows, hiv_cols, vmax=0.4,
                       panel_title="B.  HIV gene\n      expression",
                       show_row_labels=False, cohort_box_width=0.0,
                       name_lookup=name_lookup)

    cb_k = fig.colorbar(im_k, cax=cax_k)
    cb_k.set_label("KSHV log2(TPM + 1)", fontsize=17)
    cb_k.ax.tick_params(labelsize=15)
    cb_h = fig.colorbar(im_h, cax=cax_h)
    cb_h.set_label("HIV log2(TPM + 1)", fontsize=17)
    cb_h.ax.tick_params(labelsize=15)

    out = Path(out_dir) / "Figure_01.pdf"
    fig.savefig(out)
    plt.close(fig)
    return out
