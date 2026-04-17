# T-cells specific for KSHV and HIV migrate to Kaposi sarcoma tumors and persist over time

This repository contains the code and data to reproduce all figures from our
manuscript titled "T-cells specific for KSHV and HIV migrate to Kaposi sarcoma
tumors and persist over time".

## Repository Structure

```
ks_manuscript/
  src/ksplots/              Python figure-generation package
    config.py               Paths, palettes, shared constants
    cli.py                  CLI entry point (ksplots command)
    io/                     Data loaders (TCR, CIBERSORTx, Prism, metadata)
    plotting/               Shared matplotlib theme
    figures/                Figure modules (one per figure)
      figure_01.py          Figure 1: KSHV/HIV viral expression heatmaps
      figure_02.py          Figure 2: TCR repertoire diversity & pathogen annotation
      figure_03.py          Figure 3: Temporal clone tracking alluvials
      figure_04_wetlab.py   Figure 4: HIV-specific TCR functional validation
      figure_05_wetlab.py   Figure 5: High-avidity killing dose-response
      figure_06_pel.py      Figure 6: PEL cell functional assays
      supplementary/
        sfig_01.py          Supp. Figure 1: NanoString targeted expression
        sfig_02.py          Supp. Figure 2: CIBERSORTx immune deconvolution
        sfig_03_alluvial.py Supp. Figure 3: Clone tracking + COS7 reporter assays
        sfig_04_cos7.py     Supp. Figure 4: ORF6 fine epitope mapping
  zenodo_data/              All source data required for figure generation
  data/figures/             Rendered figure PDFs (output directory)
  legacy/                   Original R notebooks and Rmd (archived)
  scripts/                  Utility scripts
  tests/                    Smoke tests
  pyproject.toml            Package metadata and dependencies
```

## Quick Start

### 1. Install

```bash
pip install -e .
```

### 2. One-time data setup

If starting from the raw R data archive, convert to parquet:

```bash
python scripts/convert_rda_to_parquet.py
```

If using the Zenodo data deposit, ensure `zenodo_data/` contents are placed
under `data/` following the structure described in `zenodo_data/README.md`.

### 3. Render figures

```bash
ksplots list               # list all available figures
ksplots figure 1           # render a single figure
ksplots figure s3          # supplementary figures use 's' prefix
ksplots all                # render all figures
```

Output PDFs are written to `data/figures/` by default. Use `--out` to
specify an alternative output directory.

## Available Figures

| Key  | Output File | Description |
|------|-------------|-------------|
| `1`  | Figure_01.pdf | KSHV and HIV viral gene expression heatmaps |
| `2`  | Figure_02.pdf | TCR repertoire diversity and pathogen annotation |
| `3`  | Figure_03.pdf | Temporal clone tracking alluvials |
| `4`  | Figure_04.pdf | HIV-specific TCR functional validation |
| `5`  | Figure_05.pdf | High-avidity TCR killing dose-response |
| `6`  | Figure_06.pdf | PEL cell functional assays |
| `s1` | Supplementary_Figure_01.pdf | NanoString targeted viral expression |
| `s2` | Supplementary_Figure_02.pdf | CIBERSORTx immune deconvolution |
| `s3` | Supplementary_Figure_03.pdf | Clone tracking alluvials + COS7 reporter assays |
| `s4` | Supplementary_Figure_04.pdf | ORF6 fine epitope mapping (two TCRs) |

## Dependencies

- Python >= 3.10
- matplotlib, numpy, pandas, seaborn, scipy
- scanpy, anndata (single-cell)
- pyreadr (R data files)
- openpyxl (Excel files)
- See `pyproject.toml` for the full list.

## Data Availability

Source data for all figures are deposited on Zenodo. The `zenodo_data/`
directory contains a self-describing README with file descriptions and
a figure-to-data mapping.

## Abstract

Kaposi sarcoma-associated herpesvirus (KSHV) is the etiologic agent of Kaposi
sarcoma (KS), which causes significant morbidity and mortality worldwide,
particularly in people living with HIV (PLWH) and in sub-Saharan Africa where
KSHV seroprevalence is high. Postulating that T-cells specific for KSHV and HIV
would be attracted to KS tumors, we performed transcriptional profiling and
T-cell receptor (TCR) repertoire analysis of tumor biopsies from 144 Ugandan
adults with KS, 106 of whom were also living with HIV. We show that CD8+ T-cells
and M2-polarized macrophages are the most common immune cells in KS tumors. The
TCR repertoire of T-cells associated with KS tumors is shared across spatially
and temporally distinct tumors from the same individual. Clusters of T-cells with
predicted shared specificity for uncharacterized antigens, potentially encoded by
KSHV or HIV, comprise ~25% of the T-cells in KS tumors. Single-cell
RNA-sequencing of blood from a subset of 9 adults captured 4,283 unique TCRs
carried in 14,698 putative KSHV- or HIV-specific T-cells, which carried an
antigen-experienced effector phenotype. T-cells engineered to express a
representative sample of these TCRs showed high-avidity recognition of KSHV- or
HIV-encoded antigens. These results suggest that a polyspecific, high-avidity
KSHV- and HIV-specific T-cell response, potentially inhibited by M2 macrophages,
migrates to and localizes with KS tumors. Further analysis of KSHV- and
HIV-specific T-cells in KS tumors will provide insight into the pathogenesis of
KS and could guide the development of specific immune therapy based on adoptive
transfer or vaccination.
