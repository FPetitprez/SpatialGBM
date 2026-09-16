#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Matrix)
})

# ==============================================================================
# Functions
# ==============================================================================

aggregate_spatial_counts <- function(counts, metadata, sample_column = "sample") {
  
  stopifnot(!is.null(colnames(counts)))
  
  metadata <- as.data.frame(metadata)
  
  idx <- match(
    colnames(counts),
    rownames(metadata)
  )
  
  if (anyNA(idx)) {
    missing <- colnames(counts)[is.na(idx)]
    
    stop(
      "Spots in count matrix missing from metadata: ",
      paste(head(missing, 5), collapse = ", ")
    )
  }
  
  metadata <- metadata[idx, , drop = FALSE]
  
  sample_ids <- factor(metadata[[sample_column]])
  
  design <- sparse.model.matrix(
    ~ 0 + sample_ids
  )
  
  counts_by_sample <- counts %*% design
  
  colnames(counts_by_sample) <- levels(sample_ids)
  
  counts_by_sample
}


counts_to_cpm <- function(counts) {
  
  libsize <- Matrix::colSums(counts)
  
  if (any(libsize == 0)) {
    warning("Samples with zero library size detected.")
  }
  
  scaling_factor <- ifelse(
    libsize == 0,
    0,
    1e6 / libsize
  )
  
  t(t(counts) * scaling_factor)
}


# ==============================================================================
# Load data
# ==============================================================================

spatial_counts <- readRDS(
  "data/raw/20240430_SpatialGBMcounts.rds"
)

spatial_metadata <- readRDS(
  "data/raw/20240430_SpatialGBMmetadata.rds"
)


# ==============================================================================
# Aggregate spots into sample-level pseudobulk
# ==============================================================================

pseudobulk_counts <- aggregate_spatial_counts(
  counts = spatial_counts,
  metadata = spatial_metadata,
  sample_column = "sample"
)

pseudobulk_cpm <- counts_to_cpm(
  pseudobulk_counts
)


# ==============================================================================
# Save
# ==============================================================================

saveRDS(
  pseudobulk_counts,
  "data/processed/spatial_pseudobulk_counts.rds"
)

saveRDS(
  pseudobulk_cpm,
  "data/processed/spatial_pseudobulk_cpm.rds"
)