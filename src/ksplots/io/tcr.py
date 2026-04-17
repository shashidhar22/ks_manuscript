"""TCR repertoire and pathogen-database loaders."""
from __future__ import annotations

from functools import lru_cache
from pathlib import Path

import pandas as pd
import pyreadr

from ..config import DATA_DIR

AIRRSEQ_DIR = DATA_DIR / "airrseq"
TCR_DB_PATH = DATA_DIR / "tcr_database" / "known_tcr_database.tsv"


@lru_cache(maxsize=4)
def _load_tcr_rda() -> dict:
    """Load the small kstme_paper_tcr_tables.rda (~98K) once and cache."""
    return pyreadr.read_r(str(AIRRSEQ_DIR / "kstme_paper_tcr_tables.rda"))


def load_tcr_summary() -> pd.DataFrame:
    return _load_tcr_rda()["annotated_summary"]


def load_renyi_table() -> pd.DataFrame:
    return _load_tcr_rda()["study_renyi_table"]


# The big raw tables live as parquet (converted once from
# kstme_paper_raw_tables.rda; see scripts/convert_rda_to_parquet.py).
RAW_PARQUETS = {
    "raw": "study_raw_table.parquet",
    "nprod": "study_nprod_table.parquet",
    "annotated_nprod": "study_annotated_nprod_table.parquet",
    "amino_acid": "study_amino_acid_table.parquet",
}


def load_raw_table(name: str = "annotated_nprod", columns: list[str] | None = None) -> pd.DataFrame:
    """Load one of the slim parquet exports of kstme_paper_raw_tables.rda."""
    path = AIRRSEQ_DIR / RAW_PARQUETS[name]
    return pd.read_parquet(path, columns=columns)


def load_raw_productive(columns: list[str] | None = None) -> pd.DataFrame:
    """Default productive-aa table used by figures (annotated_nprod)."""
    return load_raw_table("annotated_nprod", columns=columns)


def load_known_tcr_database() -> pd.DataFrame:
    """Pathogen-annotated TCR reference database."""
    return pd.read_csv(TCR_DB_PATH, sep="\t")


def load_gliph_clusters() -> pd.DataFrame:
    df = pd.read_csv(AIRRSEQ_DIR / "kstme_gliph_clusters.csv")
    df.columns = [c.strip() for c in df.columns]
    return df
