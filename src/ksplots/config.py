"""Paths and shared constants."""
from __future__ import annotations

from pathlib import Path

import matplotlib as mpl
import matplotlib.font_manager as fm

# Use Liberation Sans (metrically identical to Arial) as the default font.
mpl.rcParams["font.family"] = "sans-serif"
mpl.rcParams["font.sans-serif"] = ["Liberation Sans", "Arial", "DejaVu Sans"]
# Rebuild the font cache so the new font is picked up.
fm._load_fontmanager(try_read_cache=False)

REPO_ROOT = Path(__file__).resolve().parents[2]
DATA_DIR = REPO_ROOT / "data"
FIGURES_DIR = DATA_DIR / "figures"

# RNA-seq lives via the data/rnaseq symlink (points to /scratch/...).
RNASEQ_DIR = DATA_DIR / "rnaseq"
VIRAL_RESULTS_DIR = RNASEQ_DIR / "results_viral"
VIRAL_HEATMAPS_DIR = VIRAL_RESULTS_DIR / "heatmaps"

# Cohort palette used across figures.
COHORT_COLORS = {
    "Endemic KS - NAT": "#80b1d3",
    "Endemic KS - Tumor": "#377eb8",
    "Epidemic KS - NAT": "#fb9a99",
    "Epidemic KS - Tumor": "#e41a1c",
    "Control Skin": "#999999",
}

# Wetlab TCR functional assay palette — blues for KSHV TCRs, reds for HIV TCRs.
WETLAB_TCR_COLORS = {
    # KSHV-specific TCRs (blues — matching Endemic KS palette)
    "ORF6-Specific TCR 1": "#377eb8",
    "ORF6-Specific TCR 2": "#80b1d3",
    "ORF57-Specific TCR": "#1b4f72",
    "ORF59-Specific TCR": "#5dade2",
    # HIV-specific TCRs (reds/warm — matching Epidemic KS palette)
    "Nef-Specific TCR 1": "#e41a1c",
    "Nef-Specific TCR 2": "#fb9a99",
    "Vpr-Specific TCR 1": "#d35400",
    "Vpr-Specific TCR 2": "#f39c12",
    "Pol-Specific TCR": "#8e44ad",
}

# Prism file paths for wetlab data.
WETLAB_DIR = DATA_DIR / "wetlab"
PRISM_CRA = WETLAB_DIR / "Data_CRAs" / "Data Analysis_CRA.prism"
PRISM_ELISA = WETLAB_DIR / "Data_ELISAs" / "Data Analysis_ELISAs.prism"
PRISM_CYTOKINE = (
    WETLAB_DIR / "Data_Cytokine analysis" / "Cytokine Analysis_IML"
    / "Cytokine Analysis.prism"
)
PRISM_MASTER = (
    WETLAB_DIR / "Experiments" / "Experiments"
    / "A. DATA ANALYSIS_GRAPHPAD FILE" / "Data Analysis.prism"
)

PATHOGEN_ORDER = [
    "CMV",
    "EBV",
    "HCV",
    "HIV-1",
    "HSV-2",
    "Influenza A",
    "Influenza B",
    "M.tuberculosis",
    "SARS-CoV-2",
    "Other",
    "Multi-pathogen",
    "Unknown",
]
