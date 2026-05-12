################## Single-Cell RNA-Seq Analysis with Seurat  #############

########### Seurat, DoubletFinder and SingleR scRNA-seq Pipeline ###################

## Standard Preprocessing Workflow

# Install and load required packages for single-cell RNA-seq analysis.

install.packages("devtools")
devtools::install_github('chris-mcginnis-ucsf/DoubletFinder')

library(DoubletFinder)
library(Seurat)
library(tidyverse)

# Load raw count matrix data in MTX format using ReadMtx function

cts <- ReadMtx(mtx = "GSE202501_RAW/GSM6123277_matrix.mtx.gz",
               features = "GSE202501_RAW/GSM6123277_features.tsv.gz",
               cells = "GSE202501_RAW/GSM6123277_barcodes.tsv.gz")

# Create a Seurat object with basic filtering of low-quality cells.

breastcancer <- CreateSeuratObject(counts = cts, project = "breastcancer", min.cells = 3, min.features = 200)

View(breastcancer@meta.data)

# Calculate mitochondrial gene percentage as a quality control metric. 

breastcancer$percent.mt <- PercentageFeatureSet(breastcancer, pattern = "^MT-")

# Visualize QC metrics before filtering.

VlnPlot(breastcancer, features = c("nCount_RNA", "nFeature_RNA", "percent.mt"))

# Filter out low-quality cells based on gene count and mitochondrial content.

breastcancer <- subset(breastcancer, subset = nFeature_RNA < 12500 & nFeature_RNA > 3000 & percent.mt < 20)

# Visualize QC metrics after filtering.

VlnPlot(breastcancer, features = c("nCount_RNA", "nFeature_RNA", "percent.mt"))

# Normalize gene expression data across all cells.

breastcancer <- NormalizeData(breastcancer)

# Identify highly variable features for downstream analysis.

breastcancer <- FindVariableFeatures(breastcancer)

# Scale and center the data prior to dimensional reduction.

breastcancer <- ScaleData(breastcancer)

# Perform Principal Component Analysis (PCA) for dimensionality reduction.

breastcancer <- RunPCA(breastcancer)

# Visualize principal components to select significant dimensions.

ElbowPlot(breastcancer)

# Construct nearest neighbor graph using selected principal components.

breastcancer <- FindNeighbors(breastcancer, dims = 1:15)

# Cluster cells based on gene expression similarity.

breastcancer <- FindClusters(breastcancer, resolution = 0.7)

# Perform UMAP for non-linear visualization of cell clusters.

breastcancer <- RunUMAP(breastcancer, dims = 1:15)

# Visualize clustered cells in reduced UMAP space.

DimPlot(breastcancer, reduction = "umap")

####  Doublet Detection ######

# Perform parameter sweep to evaluate optimal pK value for doublet detection.

sweep.res.list_cancer <- paramSweep(breastcancer, PCs = 1:15, sct = FALSE)

# Summarize sweep results to evaluate doublet detection performance.

sweep.stats_cancer <- summarizeSweep(sweep.res.list_cancer, GT = FALSE)

# Identify optimal pK value based on BCmetric performance.

bcmvn_cancer <- find.pK(sweep.stats_cancer)

# # Visualize BCmetric across different pK values.

ggplot(bcmvn_cancer, aes(pK, BCmetric, group = 1)) +
  geom_point() +
  geom_line()

# Select the pK value with the highest BCmetric score.

pK <- bcmvn_cancer %>% 
  filter(BCmetric == max(BCmetric)) %>% 
  select(pK) 

# Convert selected pK value into numeric format for analysis.

pK <- as.numeric(as.character(pK[[1]]))

# Extract cluster annotations for homotypic doublet estimation.

annotations <- breastcancer@meta.data$seurat_clusters 

# Estimate proportion of homotypic doublets based on cluster identity.

homotypic.prop <- modelHomotypic(annotations)         

# Estimate expected number of doublets based on dataset size.

nExp_poi <- round(0.076*nrow(breastcancer@meta.data)) 

# Adjust expected doublet number using homotypic proportion.

nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop)) 

# Run DoubletFinder to identify potential doublet cells.

breastcancer <- doubletFinder(breastcancer,           
                              PCs = 1:15,
                              pK=pK,
                              pN=0.25,
                              nExp = nExp_poi.adj)
View(breastcancer@meta.data)
head(breastcancer@meta.data)

# Visualize predicted singlets and doublets in UMAP space.

DimPlot(breastcancer, reduction = "umap", group.by = "DF.classifications_0.25_0.09_283")

# Filter out predicted doublet cells and retain only singlets.

breastcancerfiltered <- subset(breastcancer, subset = DF.classifications_0.25_0.09_283 == "Singlet")

# Re-run full preprocessing pipeline after doublet removal.

breastcancerfiltered <- breastcancerfiltered %>% 
  NormalizeData() %>% 
  FindVariableFeatures() %>% 
  ScaleData() %>% 
  RunPCA() %>% 
  FindNeighbors(dims = 1:15) %>% 
  FindClusters() %>% 
  RunUMAP(dims = 1:15)

# Visualize final clustered single-cell data.

DimPlot(breastcancerfiltered, reduction = "umap", label=T)

### Cell Type Annotation with SingleR ###

# # Install and load SingleR and reference datasets for cell type annotation.

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("SingleR")

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("celldex")

library(SingleR)
library(celldex)

# Load reference dataset containing known human primary cell types.

ref<- celldex::HumanPrimaryCellAtlasData()

View(as.data.frame(colData(ref)))

# Extract normalized expression data from Seurat object for annotation.

breast_counts <- GetAssayData(object = breastcancerfiltered, layer = "data") 

# Perform cell type prediction using SingleR against reference dataset.

pred <- SingleR(test = breast_counts,
                ref = ref,
                labels = ref$label.main)
pred

# Assign predicted cell types back into Seurat metadata.

breastcancerfiltered$singleR <- pred$labels

head(breastcancerfiltered@meta.data)

View(breastcancerfiltered@meta.data)

# Visualize annotated cell types on UMAP plot.

a1<- DimPlot(breastcancerfiltered, reduction = "umap", group.by = "SingleR", label=T, repel=T)

# Define SingleR-based identities for downstream marker analysis.

Idents(breastcancerfiltered) <- breastcancerfiltered$singleR

# Identify cluster-specific marker genes.

breast.markers <- FindAllMarkers(breastcancerfiltered, only.pos = TRUE)

# Visualize selected marker gene expression across cell populations.

VlnPlot(breastcancerfiltered, features = c("MS4A1", "CD79A"))

# Generate feature plots for selected immune-related genes.

b1<- FeaturePlot(breastcancerfiltered, features = c("GNLY", "CD3E", "CD14", "FCGR3A", "LYZ",
                                               "CD8A"), label=T, repel = T)

# Save visualization outputs as image files for reporting or publication.

ggsave("a1.png", width = 10, height = 8)
ggsave("b1.png", width = 10, height = 8)

############################################################################################################
############################################################################################################









