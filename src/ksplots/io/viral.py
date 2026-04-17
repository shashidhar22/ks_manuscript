"""HIV/KSHV viral expression heatmap loader (Salmon TPM)."""
from __future__ import annotations

from functools import lru_cache

import pandas as pd

from ..config import RNASEQ_DIR, VIRAL_HEATMAPS_DIR

KSHV_GENE_MAP_PATH = RNASEQ_DIR / "references" / "kshv_gene_map.tsv"


def load_viral_tpm(virus: str) -> pd.DataFrame:
    """Load Salmon TPM heatmap matrix (genes × samples) for 'hiv1' or 'kshv'."""
    path = VIRAL_HEATMAPS_DIR / f"{virus}_salmon_tpm_heatmap.tsv"
    df = pd.read_csv(path, sep="\t", index_col=0)
    return df


@lru_cache(maxsize=1)
def load_kshv_gene_map() -> pd.DataFrame:
    """Locus_tag → (common name, biotype) lookup parsed from the KSHV RefSeq GFF.

    Built once with::

        awk -F'\\t' '$3=="gene" && $1=="NC_009333.1" { ... }' \\
            data/rnaseq/references/combined_human_viral/viral_raw/kshv_dataset/\\
            ncbi_dataset/data/GCF_000838265.1/genomic.gff > kshv_gene_map.tsv

    Columns: ``gene_id`` (e.g. ``KSHV_HHV8GK18_gp01``), ``name`` (e.g. ``K1``;
    falls back to the locus tag when the GFF has no ``gene=`` attribute) and
    ``biotype`` (``protein_coding`` / ``misc_RNA`` / …).
    """
    df = pd.read_csv(KSHV_GENE_MAP_PATH, sep="\t",
                     names=["gene_id", "name", "biotype"])
    return df.set_index("gene_id")
