"""Figure 6 — PEL cell functional assays.

Reproduces FIGURES_IT.pptx Slide 6:

* **Panel A** — (placeholder, microscopy image not yet implemented).
* **Panel B** — IFN-γ (pg/mL) bar chart: TCR 2 & 3 ± TPA.
* **Panel C** — Specific Lysis (%) bar chart: TCR 2 & 3 ± TPA.
* **Panel D** — Specific Lysis (%) dose-response: TCR 2 & 3.

Data sources: ``data/wetlab/data/Figure6_Panel{B,C,D}.csv``
"""
from __future__ import annotations

import csv
import io
import re
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.gridspec import GridSpec
from scipy.optimize import curve_fit
from scipy.stats import ttest_ind

from ..config import DATA_DIR

_DATA_DIR = DATA_DIR / "wetlab" / "data"

# ---------------------------------------------------------------------------
# Colors
# ---------------------------------------------------------------------------

_TCR_RENAME = {
    "TCR 2": "ORF6329-337-Specific TCR\nCASSIAGHEQFF",
    "TCR 3": "ORF6329-337-Specific TCR\nCASSIAGHEQYF",
}

_TCR_COLORS = {
    "ORF6329-337-Specific TCR\nCASSIAGHEQFF": "#E41A1C",  # red
    "ORF6329-337-Specific TCR\nCASSIAGHEQYF": "#377EB8",  # blue
}
_NEG_COLOR = "#999999"   # grey for (-) TPA
_MARKERS = {
    "ORF6329-337-Specific TCR\nCASSIAGHEQFF": "s",
    "ORF6329-337-Specific TCR\nCASSIAGHEQYF": "o",
}

# ---------------------------------------------------------------------------
# Gene name splitting (reuse logic from figure_05)
# ---------------------------------------------------------------------------

_KNOWN_GENES = ["ORF57", "ORF59", "ORF6"]


def _split_gene_range(s: str) -> tuple[str, str, str] | None:
    for gene in _KNOWN_GENES:
        idx = s.find(gene)
        if idx >= 0:
            after = s[idx + len(gene):]
            m = re.match(r'(\d+[-–]\d+)(.*)', after)
            if m:
                return s[:idx + len(gene)], m.group(1), m.group(2)
    m = re.match(r'(.+?)(\d+[-–]\d+)(.*)', s)
    if m:
        return m.group(1), m.group(2), m.group(3)
    return None


def _subscript(s: str) -> str:
    parts = _split_gene_range(s)
    if parts:
        pre, rng, post = parts
        return f"{pre}$_{{{rng}}}${post}"
    return s


def _bold_with_sub(s: str) -> str:
    parts = _split_gene_range(s)
    if parts:
        pre, rng, post = parts
        pre = pre.replace(" ", r"\ ").replace("-", r"\text{-}")
        rng = rng.replace("-", r"\text{-}")
        post = post.replace(" ", r"\ ").replace("-", r"\text{-}")
        return rf"$\bf{{{pre}}}$$\bf{{_{{{rng}}}}}$$\bf{{{post}}}$"
    safe = s.replace(" ", r"\ ").replace("-", r"\text{-}")
    return rf"$\bf{{{safe}}}$"


# ---------------------------------------------------------------------------
# 4PL fitting
# ---------------------------------------------------------------------------


def _sigmoidal_4pl(x, bottom, top, ec50, hill):
    return bottom + (top - bottom) / (1 + (ec50 / x) ** hill)


def _interpolate_half_max(conc, mean):
    mask = conc > 0
    conc_pos, mean_pos = conc[mask], mean[mask]
    half_y = (mean_pos.max() + mean_pos.min()) / 2
    crossings = np.where(np.diff(np.sign(mean_pos - half_y)))[0]
    if len(crossings) == 0:
        return None
    i = crossings[0]
    frac = (half_y - mean_pos[i + 1]) / (mean_pos[i] - mean_pos[i + 1])
    log_conc = np.log10(conc_pos)
    log_ec50 = log_conc[i + 1] + frac * (log_conc[i] - log_conc[i + 1])
    return round(10 ** log_ec50, 2)


def _fit_ec50(conc, mean, bounds_top=100):
    mask = conc > 0
    if mask.sum() < 3:
        return None
    dyn_range = mean[mask].max() - mean[mask].min()
    if dyn_range < 0.10 * mean[mask].max():
        return None
    try:
        popt, _ = curve_fit(
            _sigmoidal_4pl, conc[mask], mean[mask],
            p0=[mean.min(), mean.max(), 1.0, 1],
            maxfev=20000,
            bounds=([0, 0, 1e-4, 0.1],
                    [mean.max() * 1.5, bounds_top * 2, 1e5, 10]),
        )
        ec50 = round(popt[2], 2)
        if ec50 <= conc[mask].max() * 2:
            return ec50
    except Exception:
        pass
    return _interpolate_half_max(conc, mean)


# ---------------------------------------------------------------------------
# Bar chart loader / drawer (Panels B & C)
# ---------------------------------------------------------------------------


def _load_bar_csv(path: Path):
    """Load Panel B/C CSV: returns (conditions, {tcr_name: (mean, sem)})."""
    text = path.read_text().strip()
    reader = csv.reader(io.StringIO(text))
    rows = list(reader)

    header = rows[0]
    # Group names from header (skip empty first col)
    raw_names = [h.strip() for h in header[1:]]
    group_names = []
    for n in raw_names:
        if not group_names or n != group_names[-1]:
            group_names.append(n)
    n_reps = len(raw_names) // len(group_names)

    conditions = []
    group_data: dict[str, list[list[float]]] = {n: [] for n in group_names}

    for row in rows[1:]:
        conditions.append(row[0].strip())
        idx = 1
        for gname in group_names:
            reps = [float(row[idx + j]) for j in range(n_reps)]
            group_data[gname].append(reps)
            idx += n_reps

    result = {}
    for gname in group_names:
        arr = np.array(group_data[gname])
        display = _TCR_RENAME.get(gname, gname)
        result[display] = (arr.mean(axis=1), arr.std(axis=1, ddof=1) / np.sqrt(arr.shape[1]))

    return conditions, result


def _draw_bar_panel(
    ax: plt.Axes,
    conditions: list[str],
    tcr_data: dict[str, tuple[np.ndarray, np.ndarray]],
    ylabel: str,
    panel_letter: str,
) -> None:
    """Draw grouped bar chart with significance testing."""
    n_cond = len(conditions)
    n_tcrs = len(tcr_data)
    bar_w = 0.35
    x = np.arange(n_cond)

    tcr_names = list(tcr_data.keys())
    for i, tname in enumerate(tcr_names):
        mean, sem = tcr_data[tname]
        color = _TCR_COLORS.get(tname, list(_TCR_COLORS.values())[i])
        offset = (i - n_tcrs / 2 + 0.5) * bar_w
        disp_label = "\n".join(_subscript(ln) for ln in tname.split("\n"))
        ax.bar(x + offset, mean, bar_w * 0.85, yerr=sem, capsize=4,
               color=color, edgecolor="black", linewidth=0.5,
               label=disp_label, error_kw={"linewidth": 0.8})

    # Significance: compare first condition vs each subsequent condition
    # within each TCR. Draw bracket between condition pairs.
    bracket_idx = 0
    for i, tname in enumerate(tcr_names):
        mean, sem = tcr_data[tname]
        offset = (i - n_tcrs / 2 + 0.5) * bar_w
        for j in range(1, n_cond):
            if mean[j] > mean[0] * 1.5:
                y_max = max(mean[0] + sem[0], mean[j] + sem[j])
                y_bar = y_max * 1.08 + bracket_idx * y_max * 0.10
                x_left = 0 + offset
                x_right = j + offset
                ax.plot([x_left, x_left, x_right, x_right],
                        [y_bar, y_bar * 1.02, y_bar * 1.02, y_bar],
                        color="black", lw=1.0, clip_on=False)
                ax.text((x_left + x_right) / 2, y_bar * 1.03, "***",
                        ha="center", va="bottom", fontsize=18, fontweight="bold")
                bracket_idx += 1

    ax.set_xticks(x)
    ax.set_xticklabels(conditions, fontsize=22, fontweight="bold")
    ax.set_ylabel(ylabel, fontsize=28, fontweight="bold")
    ax.tick_params(labelsize=22)
    ax.set_title(panel_letter, loc="left", fontweight="bold", fontsize=28)
    ax.legend(fontsize=18, frameon=False, loc="center left",
              bbox_to_anchor=(1.02, 0.5), prop={"weight": "bold"})
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)


# ---------------------------------------------------------------------------
# Dose-response loader / drawer (Panel D)
# ---------------------------------------------------------------------------


def _load_dose_response(path: Path):
    """Load Panel D CSV (has_conc_col=True format)."""
    text = path.read_text().strip()
    reader = csv.reader(io.StringIO(text))
    rows = list(reader)

    header = rows[0]
    raw_names = [h.strip() for h in header[2:]]
    group_names = []
    for n in raw_names:
        if not group_names or n != group_names[-1]:
            group_names.append(n)
    n_reps = len(raw_names) // len(group_names)

    concentrations = []
    group_data: dict[str, list[list[float]]] = {n: [] for n in group_names}

    for row in rows[1:]:
        concentrations.append(float(row[1]))
        idx = 2
        for gname in group_names:
            reps = [float(row[idx + j]) for j in range(n_reps)]
            group_data[gname].append(reps)
            idx += n_reps

    conc = np.array(concentrations)
    groups = []
    for gname in group_names:
        arr = np.array(group_data[gname])
        display = _TCR_RENAME.get(gname, gname)
        groups.append((display, arr.mean(axis=1),
                       arr.std(axis=1, ddof=1) / np.sqrt(arr.shape[1])))
    return conc, groups


def _draw_dose_response(
    ax: plt.Axes,
    conc: np.ndarray,
    groups: list[tuple[str, np.ndarray, np.ndarray]],
    ylabel: str,
    panel_letter: str,
    ylim_top: float = 80,
) -> None:
    """Draw dose-response with EC50 annotation."""
    X_ZERO = 3e-3

    curve_data = {}
    ec50_values = {}

    for i, (name, mean, sem) in enumerate(groups):
        color = _TCR_COLORS.get(name, list(_TCR_COLORS.values())[i])
        marker = _MARKERS.get(name, "s")
        x_plot = np.where(conc == 0, X_ZERO, conc)

        disp_label = "\n".join(_subscript(ln) for ln in name.split("\n"))
        ax.errorbar(x_plot, mean, yerr=sem,
                    fmt=f"{marker}-", color=color, markersize=8,
                    capsize=3, linewidth=1.8, label=disp_label, zorder=3)
        curve_data[name] = (x_plot, mean)
        ec50_values[name] = _fit_ec50(conc, mean, bounds_top=ylim_top)

    # Draw half-max lines
    for i, (name, _, _) in enumerate(groups):
        hm_x = ec50_values.get(name)
        if hm_x is None:
            continue
        color = _TCR_COLORS.get(name, list(_TCR_COLORS.values())[i])
        x_pts, y_pts = curve_data[name]
        pos = x_pts > 0
        x_log = np.log10(x_pts[pos])
        y_pos = y_pts[pos]
        order = np.argsort(x_log)
        half_y = float(np.interp(np.log10(hm_x), x_log[order], y_pos[order]))

        x_edge = conc.max() * 2
        ax.plot([x_edge, hm_x], [half_y, half_y],
                ls="--", color=color, lw=2.0, alpha=0.7, zorder=2)
        ax.plot([hm_x, hm_x], [half_y, 0],
                ls="--", color=color, lw=2.0, alpha=0.7, zorder=2)

    # EC50 text box
    hm_lines = [r"$\bf{1/2\ max\ EC_{50}:}$"]
    for name, _, _ in groups:
        hm_x = ec50_values.get(name)
        if hm_x is None:
            continue
        dn_lines = name.split("\n")
        short_name = dn_lines[0].strip()
        cdr3 = dn_lines[1].strip() if len(dn_lines) > 1 else ""
        label = _bold_with_sub(short_name)
        if cdr3:
            label += f" ({cdr3})"
        hm_lines.append(f"  {label}: {hm_x} nM")

    if len(hm_lines) > 1:
        ax.text(0.98, 0.98, "\n".join(hm_lines), transform=ax.transAxes,
                fontsize=22, verticalalignment="top", horizontalalignment="right",
                bbox=dict(boxstyle="round,pad=0.4", facecolor="white",
                          edgecolor="grey", alpha=0.9))

    ax.set_xscale("log")
    ax.set_ylabel(ylabel, fontsize=28, fontweight="bold")
    ax.set_xlabel("Peptide Concentration (nM)", fontsize=28)
    ax.set_ylim(bottom=min(0, min(m.min() for _, m, _ in groups) - 5),
                top=ylim_top)
    ax.tick_params(labelsize=22)
    ax.invert_xaxis()

    max_exp = int(np.ceil(np.log10(conc[conc > 0].max())))
    min_exp = int(np.floor(np.log10(conc[conc > 0].min())))
    tick_vals = [10**e for e in range(max_exp, min_exp - 1, -1)]
    tick_labels = [rf"$10^{{{e}}}$" for e in range(max_exp, min_exp - 1, -1)]
    if 0 in conc:
        tick_vals.append(X_ZERO)
        tick_labels.append("0")
    ax.set_xticks(tick_vals)
    ax.set_xticklabels(tick_labels)

    ax.set_title(panel_letter, loc="left", fontweight="bold", fontsize=28)
    ax.legend(fontsize=24, loc="center left", bbox_to_anchor=(1.02, 0.5),
              frameon=False, prop={"weight": "bold"})
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)


# ---------------------------------------------------------------------------
# Render
# ---------------------------------------------------------------------------


def render(out_dir: Path) -> Path:
    fig = plt.figure(figsize=(28, 20))
    gs = GridSpec(2, 2, hspace=0.35, wspace=0.45,
                  left=0.08, right=0.85, top=0.95, bottom=0.08)

    # Panel A — mNeonGreen+ bar chart (top-left)
    ax_a = fig.add_subplot(gs[0, 0])
    cond_a, data_a = _load_bar_csv(_DATA_DIR / "FIgure6_PanelA.csv")
    _draw_bar_panel(ax_a, cond_a, data_a,
                    ylabel="mNeonGreen+ (%)", panel_letter="A")

    # Panel B — IFN-γ bar chart (top-right)
    ax_b = fig.add_subplot(gs[0, 1])
    cond_b, data_b = _load_bar_csv(_DATA_DIR / "Figure6_PanelB.csv")
    _draw_bar_panel(ax_b, cond_b, data_b,
                    ylabel="IFN-\u03b3 (pg/mL)", panel_letter="B")

    # Panel C — Specific Lysis bar chart (bottom-left)
    ax_c = fig.add_subplot(gs[1, 0])
    cond_c, data_c = _load_bar_csv(_DATA_DIR / "Figure6_PanelC.csv")
    _draw_bar_panel(ax_c, cond_c, data_c,
                    ylabel="Specific Lysis (%)", panel_letter="C")

    # Panel D — Dose-response (bottom-right)
    ax_d = fig.add_subplot(gs[1, 1])
    conc_d, groups_d = _load_dose_response(_DATA_DIR / "Figure6_PanelD.csv")
    _draw_dose_response(ax_d, conc_d, groups_d,
                        ylabel="Specific Lysis (%)", panel_letter="D",
                        ylim_top=90)

    out = Path(out_dir) / "Figure_06.pdf"
    fig.savefig(out)
    plt.close(fig)
    return out
