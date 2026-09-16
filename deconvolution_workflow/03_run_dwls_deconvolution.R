#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(tidyverse)
  library(Seurat)
  library(omnideconv)
})


# ==============================================================================
# Functions
# ==============================================================================

subset_cells <- function(cell_matrix, annotations, batch_ids,
                         num_cells, seed) {
  
  seurat_obj <- CreateSeuratObject(
    counts = cell_matrix,
    assay = "RNA"
  )
  
  seurat_obj <- AddMetaData(
    seurat_obj,
    metadata = annotations,
    col.name = "cell_type"
  )
  
  seurat_obj <- AddMetaData(
    seurat_obj,
    metadata = batch_ids,
    col.name = "batch_id"
  )
  
  set.seed(seed)
  
  sampled_metadata <- seurat_obj@meta.data %>%
    rownames_to_column("barcode") %>%
    group_by(cell_type) %>%
    nest() %>%
    mutate(n = map_dbl(data, nrow)) %>%
    mutate(n = pmin(n, num_cells)) %>%
    ungroup() %>%
    mutate(samp = map2(data, n, sample_n)) %>%
    select(-data) %>%
    unnest(samp)
  
  sampled_obj <- subset(
    seurat_obj,
    cells = sampled_metadata$barcode
  )
  
  list(
    data = as.matrix(sampled_obj[["RNA"]]$counts),
    annotations = sampled_obj@meta.data$cell_type,
    batch_ids = sampled_obj@meta.data$batch_id
  )
}


run_dwls <- function(bulk_data, signature) {
  
  omnideconv::deconvolute_dwls(
    as.matrix(bulk_data),
    signature,
    "DampenedWLS"
  )
}


# ==============================================================================
# Load data
# ==============================================================================

reference_data <- readRDS(
  "data/processed/dwls_reference.rds"
)

pseudobulk_cpm <- readRDS(
  "data/processed/spatial_pseudobulk_cpm.rds"
)

tcga_tpm <- readRDS(
  "data/raw/20240430_TCGA_GBM_tpm.rds"
)


# ==============================================================================
# Subsample single-cell reference
# ==============================================================================

reference_subset <- subset_cells(
  cell_matrix = reference_data$counts,
  annotations = reference_data$annotations,
  batch_ids = reference_data$batch_ids,
  num_cells = 5000,
  seed = 22
)


# ==============================================================================
# Build DWLS signature
# ==============================================================================

dwls_signature <- omnideconv::build_model_dwls(
  single_cell_object = as.matrix(reference_subset$data),
  cell_type_annotations = as.vector(reference_subset$annotations),
  "mast_optimized",
  ncores = 5
)

saveRDS(
  dwls_signature,
  "results/models/dwls_signature.rds"
)


# ==============================================================================
# Spatial pseudobulk deconvolution
# ==============================================================================

dwls_pseudobulk <- run_dwls(
  bulk_data = pseudobulk_cpm,
  signature = dwls_signature
)

saveRDS(
  dwls_pseudobulk,
  "results/deconvolution/dwls_results_pseudobulk.rds"
)


# ==============================================================================
# TCGA GBM deconvolution
# ==============================================================================

dwls_tcga <- run_dwls(
  bulk_data = tcga_tpm,
  signature = dwls_signature
)

saveRDS(
  dwls_tcga,
  "results/deconvolution/dwls_results_tcga.rds"
)

message("DWLS deconvolution completed.")