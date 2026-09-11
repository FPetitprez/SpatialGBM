##############################################################
#                                                            #
# Run RCTD using the Ruiz-Moreno atlas on SpatialGBM samples #
#                                                            #
##############################################################

## This script runs RCTD using different resolutions of annotations in the atlas by Ruiz-Moreno et al. (bioRXiv 2022) on SpatialGBM samples and adds the results in the SPATA2 object.


############################
# Library and data loading #
############################


library(dplyr)
library(spacexr)
library(Matrix)
library(doParallel)
library(Seurat)
library(ggplot2)
library(SPATA2)
library(tictoc)
library(colorspace)


Ruiz_Moreno_Seurat <- readRDS("cellxGene_Seurat.rds") ## This is the file from the Ruiz-Moreno CellxGene instance.

############################################
# Create deconvolution references for RCTD #
############################################

RM_ref_counts <- Ruiz_Moreno_Seurat@assays$RNA$counts

rownames(RM_ref_counts) <- as.character(Ruiz_Moreno_Seurat@assays$RNA@meta.features$feature_name)


RM_cellTypes_large <- Ruiz_Moreno_Seurat$annotation_level_3
RM_cellTypes_large <- recode_factor(RM_cellTypes_large,
                                   'CD4/CD8' = "T cell",
                                   'Mural cell' = "Vascular",
                                   'Endothelial' = "Vascular")


RM_cellTypes_fine <- Ruiz_Moreno_Seurat$annotation_level_4

RM_cellTypes_fine <- recode_factor(RM_cellTypes_fine,
                                  'TAM-BDM hypoxia/MES' = "TAM-BDM hypoxia-MES",
                                  # Simplify vascular cell types
                                  'Pericyte' = "Vascular",
                                  'Perivascular fibroblast' = "Vascular",
                                  'Endo arterial' = "Vascular",
                                  'Endo capilar' = "Vascular",
                                  'SMC' = "Vascular",
                                  'SMC COL' = "Vascular",
                                  'SMC prolif' = "Vascular",
                                  'Scavenging endothelial' = "Vascular",
                                  'Scavenging pericyte' = "Vascular",
                                  'Tip-like' = "Vascular",
                                  'VLMC' = "Vascular",
                                  # Simplify tumor cell types
                                  'AC-like Prolif' = "AC-like",
                                  'MES-like hypoxia independent' = "MES-like",
                                  'MES-like hypoxia/MHC' = "MES-like",
                                  'NPC-like OPC' = "NPC-like",
                                  'NPC-like Prolif' = "NPC-like",
                                  'NPC-like neural' = "NPC-like",
                                  'OPC-like Prolif' = "OPC-like")


#populations to remove (less than 100 cells in atlas) and poorly annotated DC and T subsets
popToRemove <- unique(c(names(table(RM_cellTypes_large))[table(RM_cellTypes_large)<100],names(table(RM_cellTypes_fine))[table(RM_cellTypes_fine)<100], "DC1", "DC2", "DC3", "Stress sig", "Prolif T"))

Ruiz_Moreno_reference_large <-Reference(RM_ref_counts[,!RM_cellTypes_large %in% popToRemove],RM_cellTypes_large[!RM_cellTypes_large %in% popToRemove])
Ruiz_Moreno_reference_fine <- Reference(RM_ref_counts[,!RM_cellTypes_fine %in% popToRemove],RM_cellTypes_fine[!RM_cellTypes_fine %in% popToRemove])


#################
# Deconvolution #
#################



runRCTDonSPATA2 <- function(processedObject, sampleName, samplePath){
  
  spata2object <- processedObject
  
  date <- fp.strsplit(as.character(Sys.time())," ")
  
  # Open a file to send messages to
  zz <- file(paste(date,sampleName,"RCTD_log.Rout",sep="_"), open = "wt")
  # Divert messages to that file
  sink(zz, type = "message")
  
  tic(paste("Done processing sample",sampleName))
  
  
  ## Step 4: Global populations
  
  tic("RCTD for large populations")
  
  spatialRNA <- spacexr::read.VisiumSpatialRNA(paste(samplePath,"outs",sep=""))
  
  globalRCTD <- spacexr::create.RCTD(spatialRNA,Ruiz_Moreno_reference_large, max_cores = 12)
  globalRCTD <- spacexr::run.RCTD(globalRCTD,doublet_mode = "full")
  
  results <- globalRCTD@results$weights
  fractions <- data.frame(spacexr::normalize_weights(results))
  colnames(fractions) <- paste("RM_large_RCTD_",colnames(fractions),sep="")
  feats <- colnames(fractions)
  fractions$barcodes <- row.names(fractions)
  
  spata2object <- addFeatures(spata2object,fractions,colnames(fractions))
  
  toc()
  
  
  ## Step 5: Large immune populations
  
  tic("RCTD for fine populations")
  
  fineImmuneRCTD <- spacexr::create.RCTD(spatialRNA,Ruiz_Moreno_reference_fine, max_cores = 12)
  fineImmuneRCTD <- spacexr::run.RCTD(fineImmuneRCTD,doublet_mode = "full")
  
  results <- fineImmuneRCTD@results$weights
  fractions <- data.frame(spacexr::normalize_weights(results))
  colnames(fractions) <- paste("RM_fine_RCTD_",colnames(fractions),sep="")
  feats <- colnames(fractions)
  fractions$barcodes <- row.names(fractions)
  
  spata2object <- addFeatures(spata2object,fractions,colnames(fractions))
  
  toc()

  sink(type="message")
  
  toc()
  
  gc()
  
  return(spata2object)
  
}




#####################
# Load RCTD objects #
#####################


#samplesNames <-  ## char of all sample names


for(i in 1:length(samplesNames)){
  print(samplesNames[i])
  load(paste0("path/to/SPATA2/RdataFiles/",samplesNames[i],".RData"))
  eval(parse(text = paste0("res <- ",samplesNames[i])))
  res <- runRCTDonSPATA2(res,samplesNames[i],folders[i])
  #res <- updateSpataObject(res) uncomment if run on an older version of SPATA2
  gdata::mv("res",gsub("-","_",samplesNames[i]))
  eval(parse(text = paste0("save(",gsub("-","_",samplesNames[i]),", file ='",samplesNames[i],".RData')")))
}











