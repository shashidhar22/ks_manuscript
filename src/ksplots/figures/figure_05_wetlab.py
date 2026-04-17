"""Figure 5 — High-avidity TCR killing against peptide-pulsed LCLs.

Reproduces FIGURES_IT.pptx Slide 2:

* **Panel A** — mNeonGreen+ (%) dose-response for KS TCR 2 & 3 (ORF6).
* **Panel B** — CRA specific lysis dose-response for KS TCR 2, 3, 21, 24.
* **Panel C** — CRA dose-response for ORF57 peptide variants (KS TCR 21).
* **Panel D** — CRA dose-response for ORF59 peptide (KS TCR 24).

Data sources: ``data/wetlab/data/Figure5_Panel{A,B,C,D}.csv``
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

from ..config import DATA_DIR

_DATA_DIR = DATA_DIR / "wetlab" / "data"

# ---------------------------------------------------------------------------
# Shared palette and markers
# ---------------------------------------------------------------------------

_COLORS = {
    "TCR1": "#E41A1C",   # red
    "TCR2": "#377EB8",   # blue
    "TCR3": "#4DAF4A",   # green
    "TCR4": "#FF7F00",   # orange
    "TCR5": "#984EA3",   # purple
}

_MARKERS = ["s", "o", "^", "v", "D", "<", ">", "p"]


# ---------------------------------------------------------------------------
# 4PL fitting
# ---------------------------------------------------------------------------


def _sigmoidal_4pl(x, bottom, top, ec50, hill):
    """4-parameter logistic."""
    return bottom + (top - bottom) / (1 + (ec50 / x) ** hill)


def _interpolate_half_max(conc: np.ndarray, mean: np.ndarray) -> float | None:
    """Estimate EC50 by interpolating where the curve crosses the observed
    half-max value (midpoint between observed min and max)."""
    mask = conc > 0
    conc_pos = conc[mask]
    mean_pos = mean[mask]
    half_y = (mean_pos.max() + mean_pos.min()) / 2

    # Find crossing point
    crossings = np.where(np.diff(np.sign(mean_pos - half_y)))[0]
    if len(crossings) == 0:
        return None
    i = crossings[0]
    # Linear interpolation in log-concentration space
    frac = (half_y - mean_pos[i + 1]) / (mean_pos[i] - mean_pos[i + 1])
    log_conc = np.log10(conc_pos)
    log_ec50 = log_conc[i + 1] + frac * (log_conc[i] - log_conc[i + 1])
    return round(10 ** log_ec50, 2)


def _fit_ec50(conc: np.ndarray, mean: np.ndarray,
              bounds_top: float = 100) -> float | None:
    """Fit 4PL and return EC50, or None on failure.

    Returns None if the curve has insufficient dynamic range (<10% of
    the max value).  Falls back to observed half-max interpolation if
    the 4PL EC50 exceeds the tested concentration range (e.g. when the
    curve has not plateaued).
    """
    mask = conc > 0
    if mask.sum() < 3:
        return None

    # Skip flat / non-responsive curves
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
        conc_pos = conc[mask]
        if ec50 <= conc_pos.max() * 2:
            return ec50
    except Exception:
        pass

    # Fallback: interpolate from observed half-max
    return _interpolate_half_max(conc, mean)


# ---------------------------------------------------------------------------
# Subscript helper
# ---------------------------------------------------------------------------


_KNOWN_GENES = ["ORF57", "ORF59", "ORF6"]  # longest first for matching


def _split_gene_range(s: str) -> tuple[str, str, str] | None:
    """Split 'ORF57382-390' into ('ORF57', '382-390', rest).

    Uses a list of known gene names to correctly split the boundary
    between gene number and residue range (e.g. ORF6 vs ORF63).
    """
    for gene in _KNOWN_GENES:
        idx = s.find(gene)
        if idx >= 0:
            after = s[idx + len(gene):]
            m = re.match(r'(\d+[-–]\d+)(.*)', after)
            if m:
                pre = s[:idx + len(gene)]
                return pre, m.group(1), m.group(2)
    # Fallback: generic pattern for non-ORF names (e.g. KS TCR 2)
    m = re.match(r'(.+?)(\d+[-–]\d+)(.*)', s)
    if m:
        return m.group(1), m.group(2), m.group(3)
    return None


def _subscript(s: str) -> str:
    """Convert residue ranges to subscript, e.g. ORF57382-390 -> ORF57$_{382-390}$."""
    parts = _split_gene_range(s)
    if parts:
        pre, rng, post = parts
        return f"{pre}$_{{{rng}}}${post}"
    return s


def _bold_with_sub(s: str) -> str:
    """Bold text with subscript ranges for annotation box."""
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
# Generic CSV loader
# ---------------------------------------------------------------------------


def _load_dose_response(path: Path, has_conc_col: bool = False):
    """Load a dose-response CSV.

    Returns (concentrations, groups) where groups is a list of
    (name, mean, sem) tuples.

    If *has_conc_col* is True, col1 is the numeric concentration
    (Panel A format: col0=label, col1=conc, col2+=data).
    Otherwise col0 is the concentration label (Panels B-D format).
    """
    text = path.read_text().strip().replace("\r\n", "\n")
    reader = csv.reader(io.StringIO(text))
    rows = list(reader)

    header = rows[0]

    # Parse group names from header — multi-line names joined
    if has_conc_col:
        data_start = 2  # skip label + conc columns
    else:
        data_start = 1  # skip label column

    raw_names = header[data_start:]
    # Group consecutive identical names (3 replicates each)
    group_names: list[str] = []
    for name in raw_names:
        name = name.strip().replace("\r", " ").replace("\n", "\n")
        if not group_names or name != group_names[-1]:
            group_names.append(name)

    # Count replicates per group
    n_reps = len(raw_names) // len(group_names)

    # Parse data rows
    concentrations = []
    group_data: dict[str, list[list[float]]] = {n: [] for n in group_names}

    for row in rows[1:]:
        if not row or all(v.strip() == "" for v in row):
            continue

        if has_conc_col:
            conc_val = float(row[1])
        else:
            # Parse concentration from label like "10³", "10²", "10-1", "0"
            label = row[0].strip()
            if label == "0":
                conc_val = 0.0
            else:
                # Labels are "10^exponent" written as e.g. "103" (=10^3),
                # "10-1" (=10^-1), "100" (=10^0).  Detect this pattern
                # first before falling back to plain float.
                label_clean = label.replace("−", "-").replace("–", "-")
                m = re.match(r'^10([-]?\d+)$', label_clean)
                if m:
                    conc_val = 10 ** int(m.group(1))
                else:
                    conc_val = float(label_clean)

        concentrations.append(conc_val)

        idx = data_start
        for gname in group_names:
            reps = []
            for j in range(n_reps):
                v = row[idx + j].strip() if idx + j < len(row) else ""
                reps.append(float(v) if v else np.nan)
            group_data[gname].append(reps)
            idx += n_reps

    conc = np.array(concentrations)
    groups = []
    for gname in group_names:
        arr = np.array(group_data[gname])
        mean = np.nanmean(arr, axis=1)
        sem = np.nanstd(arr, axis=1, ddof=1) / np.sqrt(
            np.sum(~np.isnan(arr), axis=1)
        )
        groups.append((gname, mean, sem))

    return conc, groups


# ---------------------------------------------------------------------------
# Dose-response drawing
# ---------------------------------------------------------------------------


def _draw_dose_response(
    ax: plt.Axes,
    conc: np.ndarray,
    groups: list[tuple[str, np.ndarray, np.ndarray]],
    ylabel: str,
    panel_letter: str,
    ylim_top: float = 100,
    invert_x: bool = True,
    rename: dict[str, str] | None = None,
    ec50_x: float = 0.98,
    show_ec50: bool = True,
) -> None:
    """Draw dose-response curves with EC50 annotation."""
    X_ZERO = 3e-3  # pseudo-x for concentration = 0

    curve_data = {}
    ec50_values = {}

    for i, (name, mean, sem) in enumerate(groups):
        color = list(_COLORS.values())[i % len(_COLORS)]
        marker = _MARKERS[i % len(_MARKERS)]
        x_plot = np.where(conc == 0, X_ZERO, conc)

        # Apply rename mapping if provided
        display_name = rename.get(name, name) if rename else name

        # Build display name with subscripts
        disp_lines = display_name.split("\n")
        disp_name = "\n".join(_subscript(ln) for ln in disp_lines)

        ax.errorbar(x_plot, mean, yerr=sem,
                    fmt=f"{marker}-", color=color, markersize=8,
                    capsize=3, linewidth=1.8, label=disp_name, zorder=3)
        curve_data[name] = (x_plot, mean)
        if show_ec50:
            ec50_values[name] = _fit_ec50(conc, mean, bounds_top=ylim_top)

    if show_ec50:
        # Draw 1/2 max marker lines
        for i, (name, _, _) in enumerate(groups):
            hm_x = ec50_values.get(name)
            if hm_x is None:
                continue
            color = list(_COLORS.values())[i % len(_COLORS)]
            x_pts, y_pts = curve_data[name]

            pos = x_pts > 0
            x_log = np.log10(x_pts[pos])
            y_pos = y_pts[pos]
            order = np.argsort(x_log)
            half_y = float(np.interp(np.log10(hm_x), x_log[order], y_pos[order]))

            x_edge = conc.max() * 2 if not invert_x else conc.max() * 2
            ax.plot([x_edge, hm_x], [half_y, half_y],
                    ls="--", color=color, lw=2.0, alpha=0.7, zorder=2)
            ax.plot([hm_x, hm_x], [half_y, 0],
                    ls="--", color=color, lw=2.0, alpha=0.7, zorder=2)

        # Text box with EC50 values
        hm_lines = [r"$\bf{1/2\ max\ EC_{50}:}$"]
        for name, _, _ in groups:
            hm_x = ec50_values.get(name)
            if hm_x is None:
                continue
            display_name = rename.get(name, name) if rename else name
            dn_lines = display_name.split("\n")
            short_name = dn_lines[0].strip()
            cdr3 = dn_lines[1].strip() if len(dn_lines) > 1 else ""
            label = _bold_with_sub(short_name)
            if cdr3:
                label += f" ({cdr3})"
            hm_lines.append(f"  {label}: {hm_x} nM")

        if len(hm_lines) > 1:
            ax.text(ec50_x, 0.98, "\n".join(hm_lines), transform=ax.transAxes,
                    fontsize=22, verticalalignment="top", horizontalalignment="right",
                    bbox=dict(boxstyle="round,pad=0.4", facecolor="white",
                              edgecolor="grey", alpha=0.9))

    # Axes
    ax.set_xscale("log")
    ax.set_ylabel(ylabel, fontsize=28, fontweight="bold")
    ax.set_xlabel("Peptide Concentration (nM)", fontsize=28)
    ax.set_ylim(0, ylim_top)
    ax.tick_params(labelsize=22)

    if invert_x:
        ax.invert_xaxis()

    # X-axis ticks matching the data range
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


_TCR_RENAME = {
    "KS TCR 2 \nCASSIAGHEQFF": "ORF6329-337-Specific TCR\nCASSIAGHEQFF",
    "KS TCR 3 \nCASSIAGHEQYF": "ORF6329-337-Specific TCR\nCASSIAGHEQYF",
    "KS TCR 2\nCASSIAGHEQFF": "ORF6329-337-Specific TCR\nCASSIAGHEQFF",
    "KS TCR 3\nCASSIAGHEQYF": "ORF6329-337-Specific TCR\nCASSIAGHEQYF",
}


def render(out_dir: Path) -> Path:
    fig = plt.figure(figsize=(28, 20))
    gs = GridSpec(2, 2, hspace=0.35, wspace=0.55,
                  left=0.07, right=0.82, top=0.95, bottom=0.08)

    # Panel A — Specific Lysis dose-response (ORF6, 2 TCRs)
    ax_a = fig.add_subplot(gs[0, 0])
    conc_a, groups_a = _load_dose_response(
        _DATA_DIR / "Figure5_PanelA.csv", has_conc_col=True
    )
    _draw_dose_response(ax_a, conc_a, groups_a,
                        ylabel="Specific Lysis (%)", panel_letter="A",
                        ylim_top=80, rename=_TCR_RENAME, ec50_x=1.15)

    # Panel B — IFN-gamma dose-response
    ax_b = fig.add_subplot(gs[0, 1])
    conc_b, groups_b = _load_dose_response(
        _DATA_DIR / "Figure5_PanelB.csv", has_conc_col=False
    )
    _draw_dose_response(ax_b, conc_b, groups_b,
                        ylabel="IFN-\u03b3 (pg/mL)", panel_letter="B",
                        ylim_top=max(m.max() for _, m, _ in groups_b) * 1.2,
                        rename=_TCR_RENAME, show_ec50=False)

    # Panel C — mNeonGreen dose-response (ORF57 peptide variants)
    ax_c = fig.add_subplot(gs[1, 0])
    conc_c, groups_c = _load_dose_response(
        _DATA_DIR / "Figure5_PanelC.csv", has_conc_col=False
    )
    _draw_dose_response(ax_c, conc_c, groups_c,
                        ylabel="mNeonGreen+ (%)", panel_letter="C",
                        ylim_top=100)

    # Panel D — mNeonGreen dose-response (ORF59)
    ax_d = fig.add_subplot(gs[1, 1])
    conc_d, groups_d = _load_dose_response(
        _DATA_DIR / "Figure5_PanelD.csv", has_conc_col=False
    )
    _draw_dose_response(ax_d, conc_d, groups_d,
                        ylabel="mNeonGreen+ (%)", panel_letter="D",
                        ylim_top=80)

    out = Path(out_dir) / "Figure_05.pdf"
    fig.savefig(out)
    plt.close(fig)
    return out
