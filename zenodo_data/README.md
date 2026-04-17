# KS Tumor Microenvironment Manuscript - Source Data

Source data for all figures in the KS tumor microenvironment manuscript.
These files are the minimal set required to reproduce every main and
supplementary figure using the `ksplots` Python package included in
the companion code repository.

## Directory Structure

```
zenodo_data/
  rnaseq/                  RNA-seq viral expression data (Figures 1)
  airrseq/                 AIRR-seq TCR repertoire data (Figures 2, 3)
  cibersort/               CIBERSORTx immune deconvolution (Supp. Figure 2)
  tcr_database/            Known pathogen-associated TCR reference (Figures 2, 3)
  nanostring/              NanoString targeted expression (Supp. Figure 1)
  metadata/                Study-wide sample metadata
  wetlab/
    panels/                Functional assay summary data (Figures 4-6, Supp. 3-4)
    immunoseq/             ImmunoSEQ TCR repertoire files (Supp. Figure 3)
```

## File Descriptions

### rnaseq/
| File | Description |
|------|-------------|
| kshv_salmon_tpm_heatmap.tsv | KSHV gene expression (TPM) across all samples |
| hiv1_salmon_tpm_heatmap.tsv | HIV-1 gene expression (TPM) across all samples |
| kshv_gene_map.tsv | KSHV gene ID to name and lytic stage mapping |

### airrseq/
| File | Description |
|------|-------------|
| study_annotated_nprod_table.parquet | Productive TCR sequences with cohort annotations |
| kstme_gliph_clusters.csv | GLIPH2 clustering results for TCR specificity grouping |
| kstme_paper_tcr_tables.rda | R data archive with TCR summary and Renyi diversity tables |

### cibersort/
| File | Description |
|------|-------------|
| hippos_endemic.csv | CIBERSORTx results, Uganda endemic KS tumors |
| hippos_epidemic.csv | CIBERSORTx results, Uganda epidemic KS tumors |
| lid_batchone.csv - lid_batchfour.csv | CIBERSORTx results, Tanzania 1 cohort (4 batches) |
| tso_epidemic.csv | CIBERSORTx results, Tanzania 2 epidemic KS |
| gtex_nes.csv | CIBERSORTx results, GTEx non-sun-exposed skin controls |
| gtex_ses.csv | CIBERSORTx results, GTEx sun-exposed skin controls |

### tcr_database/
| File | Description |
|------|-------------|
| known_tcr_database.tsv | Curated database of pathogen-associated TCR CDR3 sequences |

### nanostring/
| File | Description |
|------|-------------|
| nanostring_counts.xlsx | NanoString nCounter gene expression counts (Sheet 4 used) |

### metadata/
| File | Description |
|------|-------------|
| KSTME_study_metadata.csv | Per-sample metadata: cohort, tissue type, HLA alleles, clinical response |

### wetlab/panels/
| File | Figure | Description |
|------|--------|-------------|
| Figure4_PanelA.csv | Fig 4A | HLA restriction bar chart data |
| Figure4_PanelB.csv | Fig 4B | CRA dose-response data |
| Figure5_PanelA.csv | Fig 5A | Specific lysis dose-response (ORF6 TCRs) |
| Figure5_PanelB.csv | Fig 5B | IFN-gamma dose-response (ORF6 TCRs) |
| Figure5_PanelC.csv | Fig 5C | mNeonGreen dose-response (ORF57 peptide variants) |
| Figure5_PanelD.csv | Fig 5D | mNeonGreen dose-response (ORF59 peptide) |
| FIgure6_PanelA.csv | Fig 6A | mNeonGreen+ bar chart (PEL, +/- TPA) |
| Figure6_PanelB.csv | Fig 6B | IFN-gamma bar chart (PEL, +/- TPA) |
| Figure6_PanelC.csv | Fig 6C | Specific lysis bar chart (PEL, +/- TPA) |
| Figure6_PanelD.csv | Fig 6D | Specific lysis dose-response (PEL) |
| SpFigure3_PanelC.csv | Supp 3C, 4B | ORF6 fine epitope mapping (+/- HLA-B*45:01) |
| SpFigure3_PanelD.csv | Supp 3D, 4C | ORF6 segment mapping (+/- HLA-B*45:01) |
| SpFigure3_PanelE.csv | Supp 3E | ORF57 (+/- HLA-A*66:01) |
| SpFigure3_PanelF.csv | Supp 3F | ORF59 (+/- HLA-B*57:03) |

### wetlab/immunoseq/
| File | Description |
|------|-------------|
| 008_098_{A,B,C}.tsv | ImmunoSEQ TRB repertoires, subject 008_098 (HIV-, CR) |
| 008_002_{A,B,C,D,H,I}.tsv | ImmunoSEQ TRB repertoires, subject 008_002 (HIV+, PR) |

## Figure-to-Data Mapping

| Figure | Data Sources |
|--------|-------------|
| Figure 1 | rnaseq/ |
| Figure 2 | airrseq/, tcr_database/, metadata/ |
| Figure 3 | airrseq/, tcr_database/ |
| Figure 4 | wetlab/panels/Figure4_* |
| Figure 5 | wetlab/panels/Figure5_* |
| Figure 6 | wetlab/panels/FIgure6_*, Figure6_* |
| Supp. Figure 1 | nanostring/ |
| Supp. Figure 2 | cibersort/, metadata/ |
| Supp. Figure 3 | wetlab/immunoseq/, wetlab/panels/SpFigure3_* |
| Supp. Figure 4 | wetlab/panels/SpFigure3_PanelC.csv, SpFigure3_PanelD.csv |

## Reproducing Figures

```bash
pip install -e .          # install ksplots from the code repository
ksplots list              # show available figures
ksplots figure 1          # render Figure 1
ksplots all               # render all figures
```

Output PDFs are written to `data/figures/` by default.
