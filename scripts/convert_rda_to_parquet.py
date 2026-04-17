"""One-time conversion of kstme_paper_raw_tables.rda to slim parquet files.

The original .rda is ~575 MB and balloons to ~23 GB in RAM when loaded via
pyreadr (which holds the full study_raw_table at 6.25M × 145 cols). We only
need a handful of columns for figure generation, so this script extracts them
once and writes one parquet per inner data frame.

Run once after cloning the repo:

    python scripts/convert_rda_to_parquet.py
"""
from __future__ import annotations

import gc
from pathlib import Path

import pyreadr

AIRRSEQ_DIR = Path(__file__).resolve().parents[1] / "data" / "airrseq"
RDA = AIRRSEQ_DIR / "kstme_paper_raw_tables.rda"
KEEP = [
    "repertoire_id",
    "junction",
    "junction_aa",
    "duplicate_count",
    "duplicate_frequency",
    "v_call",
    "j_call",
    "productive",
]


def main() -> None:
    print(f"Loading {RDA} (this can take ~10 min and ~25 GB RAM)...", flush=True)
    objs = pyreadr.read_r(str(RDA))
    for key in list(objs):
        df = objs[key]
        cols = [c for c in KEEP if c in df.columns]
        slim = df[cols] if cols else df
        out = AIRRSEQ_DIR / f"{key}.parquet"
        slim.to_parquet(out, index=False)
        print(f"  wrote {out} {slim.shape}", flush=True)
        del objs[key], df, slim
        gc.collect()


if __name__ == "__main__":
    main()
