"""CibersortX cell-type fraction loaders."""
from __future__ import annotations

from pathlib import Path

import pandas as pd

from ..config import DATA_DIR

CIBERSORT_DIR = DATA_DIR / "cibersort"

COHORT_FILES = {
    "Hippos Endemic": "hippos_endemic.csv",
    "Hippos Epidemic": "hippos_epidemic.csv",
    "Lidenge batch 1": "lid_batchone.csv",
    "Lidenge batch 2": "lid_batchtwo.csv",
    "Lidenge batch 3": "lid_batchthree.csv",
    "Lidenge batch 4": "lid_batchfour.csv",
    "Tso Epidemic": "tso_epidemic.csv",
}

DROP_COLS = {"P-value", "Correlation", "RMSE", "Absolute score (sig.score)"}


def load_cibersort_long() -> pd.DataFrame:
    rows = []
    for label, fname in COHORT_FILES.items():
        path = CIBERSORT_DIR / fname
        if not path.exists():
            continue
        df = pd.read_csv(path)
        df = df.rename(columns={df.columns[0]: "sample"})
        cell_cols = [c for c in df.columns if c not in {"sample"} | DROP_COLS]
        long = df.melt(id_vars=["sample"], value_vars=cell_cols,
                       var_name="cell_type", value_name="fraction")
        long["cohort"] = label
        rows.append(long)
    return pd.concat(rows, ignore_index=True) if rows else pd.DataFrame()
