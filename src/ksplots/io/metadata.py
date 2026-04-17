"""Study metadata loader."""
from __future__ import annotations

import pandas as pd

from ..config import DATA_DIR


def load_study_metadata() -> pd.DataFrame:
    """Load KSTME study metadata as a DataFrame."""
    df = pd.read_csv(DATA_DIR / "KSTME_study_metadata.csv", low_memory=False)
    return df


def cohort_label(df: pd.DataFrame) -> pd.Series:
    """Return an 'Endemic/Epidemic KS - NAT/Tumor' style label per row."""
    phenotype = df["phenotype"].fillna("")
    tissue = df["tissue_type"].fillna("")
    return phenotype.str.strip() + " - " + tissue.str.strip()
