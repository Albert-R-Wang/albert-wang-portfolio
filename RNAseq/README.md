# RNA-seq Analysis

This folder contains RNA-seq analysis workflows and scripts developed across multiple research projects in cancer biology and computational genomics.

The workflows included here cover RNA-seq data processing, differential expression analysis, visualization, and downstream biological interpretation across bulk and single-cell datasets.

## Contents

- `DESeq2_GO_analysis.R`: A complete RNA-seq differential expression analysis using the DESeq2 package. This script also generates volcano plots and VST-normalized expression data, and conducts gene ontology (GO) enrichment analysis.
- `TCGAanalysis_external-gene-list.R`: This R script implements an RNA-seq analysis pipeline using The Cancer Genome Atlas Lung Adenocarcinoma (TCGA-LUAD) dataset, integrating an external lung–brain metastasis (LBM) gene signature to generate annotated heatmaps and perform clustering and survival analysis.
- `TCGA_DESeq2_analysis.Rmd`: This R Markdown document presents an RNA-seq differential expression analysis pipeline for TCGA-LUAD data using DESeq2. It covers data acquisition via TCGAbiolinks, preprocessing, identification of differentially expressed genes (DEGs), visualization, and biological interpretation.
  - `TCGA_DESeq2_analysis.pdf`: Example output of the TCGA DESeq2 analysis R Markdown document.
- `NicheNet_LBM.R`: A tailored NicheNet pipeline for the LBM project, this script identifies and visualizes ligand–target regulatory networks driving intercellular communication between defined sender and receiver cell populations.
- `scRNAseq_SeuratV5_analysis.R`: This R script provides a comprehensive pipeline for single-cell RNA-seq analysis using Seurat (v5).

## Notes

- The code is provided as-is for reference and educational purposes.
- Project-specific data and raw input files are not included due to data sharing restrictions.
- Example scripts and templates are included to demonstrate key analysis steps.
