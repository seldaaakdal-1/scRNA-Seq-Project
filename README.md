# scRNA-seq Analysis Pipeline (Seurat + DoubletFinder + SingleR)

This repository contains a single-cell RNA-seq (scRNA-seq) analysis pipeline implemented in R using Seurat, DoubletFinder, and SingleR. The workflow includes preprocessing, quality control, clustering, doublet detection, and cell type annotation.

## Workflow

- Load raw scRNA-seq data (MTX format)
- Quality control and filtering (gene count, mitochondrial content)
- Normalization and feature selection
- Dimensionality reduction (PCA, UMAP)
- Clustering of cells
- Doublet detection with DoubletFinder
- Removal of doublets and re-analysis
- Cell type annotation using SingleR
- Marker gene analysis

## Tools

- Seurat
- DoubletFinder
- SingleR
- celldex
- tidyverse

## Output

- UMAP plots
- Cell clusters
- Doublet classification
- Cell type annotations
- Marker gene visualizations

## Note

This pipeline is based on a breast cancer scRNA-seq dataset but can be adapted to other datasets.

