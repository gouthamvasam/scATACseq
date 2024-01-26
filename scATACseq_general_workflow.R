#!/usr/bin/env Rscript

# Load necessary libraries and install them if they are not already installed
required_packages <- c("ggplot2", "cowplot")
bioconductor_packages <- c("Signac", "Seurat", "GenomicRanges", "EnsDb.Hsapiens.v86")

# Install CRAN packages if they are not installed
for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
  }
}

# Install Bioconductor Manager if it's not installed
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

# Install Bioconductor packages if they are not installed
BiocManager::install(bioconductor_packages)

# Load necessary libraries
library(Signac)
library(Seurat)
library(GenomicRanges)
library(EnsDb.Hsapiens.v86)
library(ggplot2)
library(cowplot)

# Set paths to input and output directories
input_dir <- '/path/to/cellranger_atac/output'
output_dir <- '/path/to/analysis/output'

# Ensure output directory exists
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Load the Cell Ranger ATAC output as a Seurat object
counts <- Read10X_h5(filename = file.path(input_dir, 'filtered_peak_bc_matrix.h5'))
metadata <- read.csv(file.path(input_dir, 'metadata.csv'), header = TRUE, row.names = 1)

# Create a Seurat object
seurat_obj <- CreateSeuratObject(counts = counts, project = 'ATAC', meta.data = metadata)

# Perform quality control
seurat_obj <- subset(seurat_obj, subset = nCount_ATAC > 3000 & nCount_ATAC < 50000 & percent.mt < 5)

# Normalize the data
seurat_obj <- NormalizeData(seurat_obj, normalization.method = 'LogNormalize', scale.factor = 10000)

# Identify the most significant peaks
seurat_obj <- FindTopFeatures(seurat_obj, min.cutoff = 'q0.5')

# Perform linear dimensional reduction
seurat_obj <- RunTFIDF(seurat_obj)
seurat_obj <- FindTopFeatures(seurat_obj, min.cutoff = 20)
seurat_obj <- RunSVD(seurat_obj)
seurat_obj <- RunUMAP(seurat_obj, dims = 2:30)

# Perform graph-based clustering
seurat_obj <- FindNeighbors(seurat_obj, reduction = 'lsi', dims = 2:30)
seurat_obj <- FindClusters(seurat_obj, resolution = 0.5)

# Find differentially accessible peaks between clusters
da_peaks <- FindMarkers(seurat_obj, ident.1 = 'cluster1', ident.2 = 'cluster2', min.pct = 0.2, test.use = 'LR', latent.vars = 'nCount_ATAC')

# Annotate the clusters
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
genome(annotations) <- 'hg38'
Annotation(seurat_obj) <- annotations

# Visualize the results
p1 <- DimPlot(seurat_obj, reduction = 'umap', group.by = 'clusters')
p2 <- FeaturePlot(seurat_obj, features = head(rownames(da_peaks), 5))

# Save plots
ggsave(file.path(output_dir, 'umap_clusters.pdf'), p1)
ggsave(file.path(output_dir, 'feature_plots.pdf'), p2)

# Save the Seurat object
saveRDS(seurat_obj, file.path(output_dir, 'seurat_obj.rds'))