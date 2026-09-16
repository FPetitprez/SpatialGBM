#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(tidyverse)
})

set.seed(22)

# ==============================================================================
# Functions
# ==============================================================================

is_outlier <- function(seurat_obj, metric, nmads) {
  values <- seurat_obj[[metric, drop = TRUE]]
  
  values < median(values) - nmads * mad(values) |
    values > median(values) + nmads * mad(values)
}


remove_outliers <- function(
    seurat_obj,
    mad_nfeature = 5,
    mad_ncount = 5,
    mad_mito = 3
) {
  
  keep <- !is_outlier(seurat_obj, "nFeature_RNA", mad_nfeature) &
    !is_outlier(seurat_obj, "nCount_RNA", mad_ncount)
  
  if ("percent.mt" %in% colnames(seurat_obj@meta.data)) {
    keep <- keep &
      !is_outlier(seurat_obj, "percent.mt", mad_mito)
  }
  
  cells_to_keep <- colnames(seurat_obj)[keep]
  
  message(
    "Cells retained after QC: ",
    length(cells_to_keep),
    " / ",
    ncol(seurat_obj),
    " (",
    round(100 * length(cells_to_keep) / ncol(seurat_obj), 1),
    "%)"
  )
  
  subset(seurat_obj, cells = cells_to_keep)
}


preprocess_reference <- function(
    seurat_obj,
    n_variable_features = 1500,
    n_pcs = 100
) {
  
  seurat_obj <- NormalizeData(
    seurat_obj,
    normalization.method = "LogNormalize",
    scale.factor = 10000
  )
  
  seurat_obj <- FindVariableFeatures(
    seurat_obj,
    selection.method = "vst",
    nfeatures = n_variable_features
  )
  
  seurat_obj <- ScaleData(
    seurat_obj,
    features = rownames(seurat_obj)
  )
  
  seurat_obj <- RunPCA(
    seurat_obj,
    features = VariableFeatures(seurat_obj),
    npcs = n_pcs
  )
  
  seurat_obj
}


# ==============================================================================
# Create reference
# ==============================================================================

GBM_data <- CreateSeuratObject(
  counts = spatial_counts,
  meta.data = spatial_meta
)

GBM_data <- remove_outliers(
  GBM_data,
  mad_nfeature = 5,
  mad_ncount = 5,
  mad_mito = 3
)

GBM_data <- preprocess_reference(
  GBM_data,
  n_variable_features = 1500,
  n_pcs = 100
)

GBM_data <- FindNeighbors(
  GBM_data,
  dims = 1:40
)

GBM_data <- FindClusters(
  GBM_data,
  resolution = 0.1
)

GBM_data <- RunUMAP(
  GBM_data,
  dims = 1:40
)


# ==============================================================================
# Niche-based cluster filtering
# ==============================================================================

cluster_niche_map <- c(
  "0"  = "N4",
  "1"  = "N4",
  "2"  = "N4",
  "3"  = "N6",
  "4"  = "N3",
  "5"  = "N2",
  "6"  = "N4",
  "7"  = "N3",
  "8"  = "N5",
  "9"  = "N5",
  "10" = "N2",
  "11" = "N1",
  "12" = "N6",
  "13" = "N6",
  "14" = "N6",
  "15" = "N1",
  "16" = "N1",
  "17" = "N6"
)

meta <- GBM_data@meta.data

expected_niche <- cluster_niche_map[
  as.character(meta$seurat_clusters)
]

cells_to_keep <- rownames(meta)[
  !is.na(expected_niche) &
    as.character(meta$niche) == expected_niche
]

reference <- subset(
  GBM_data,
  cells = cells_to_keep
)

message(
  "Cells retained after niche filtering: ",
  ncol(reference)
)


# ==============================================================================
# Export DWLS reference
# ==============================================================================

reference_data <- list(
  counts = reference[["RNA"]]$counts,
  annotations = reference$niche,
  batch_ids = reference$sample
)

saveRDS(
  reference_data,
  "data/processed/dwls_reference.rds"
)

saveRDS(
  reference,
  "data/processed/dwls_reference_seurat.rds"
)