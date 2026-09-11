######################
# Manuscript figures #
######################


## This script assumes the All_RCTD_large and All_RCTD_largeobject are loaded (see script 01 to 04)



source("./SpatialGBM_utils.R")


library("spacexr")
library("Matrix")
library("doParallel")
library("Seurat")
library("ggplot2")
library("ggpubr")
library("SPATA2")
library("colorspace")
library("patchwork")
library("ComplexHeatmap")
library("circlize")
library("dplyr")
library("CellChat")
library("forestplot")


##################################################################################################################################################################################################################################################

############
# Figure 1 #
############

## Figure 1e: spatial heatmap on 2 replicates

samplesNames = c("s08_E53","s16_E53")


for(s in samplesNames[]){

  print(s)

  # load data
  load(paste0("path/to/RData/files/",s,".RData"))
  
  # add nuclei count
  eval(parse(text = paste0(s,"<- addNucleiCount(",s,",'path/to/cell/counts",s,"/",s,"_detectedCells.txt')")))

}


cancerPops <- c("RM_large_RCTD_MES.like","RM_large_RCTD_AC.like","RM_large_RCTD_NPC.like","RM_large_RCTD_OPC.like")
TAMsPops <- c("RM_large_RCTD_TAM.BDM","RM_large_RCTD_TAM.MG")
OtherPops <- c("RM_large_RCTD_Vascular","RM_large_RCTD_Oligodendrocyte")

SpatialGBM.plotPops(s16_E53,cancerPops, ncol = 4, pt_size = 2.5)
SpatialGBM.plotPops(s16_E53,TAMsPops, pt_size = 2.5) + SpatialGBM.plotPops(s16_E53,OtherPops, pt_size = 2.5)

SpatialGBM.plotPops(s08_E53,cancerPops, ncol = 4, pt_size = 3)
SpatialGBM.plotPops(s08_E53,TAMsPops, pt_size = 3) + SpatialGBM.plotPops(s08_E53,OtherPops, pt_size = 3)

## Figure 1d
SpatialGBM.plotPops(s16_E53,"nuclei.count", as.cell.counts = FALSE, pt_size = 2.5)


##################################################################################################################################################################################################################################################


############
# Figure 2 #
############

## Figure 2a-b: heatmaps of cancer and TAMs in other samples

samplesNames <- c("UKF260_T","s04_E57")

for(s in samplesNames[]){
  
  print(s)
  
  # load data
  load(paste0("path/to/RData/files/",s,".RData"))
  
  # add nuclei count
  eval(parse(text = paste0(s,"<- addNucleiCount(",s,",'path/to/cell/counts",s,"/",s,"_detectedCells.txt')")))
  
}



SpatialGBM.plotPops(UKF260_T,cancerPops, ncol = 4, pt_size = 1.5)
SpatialGBM.plotPops(UKF260_T,TAMsPops, pt_size = 1.5)

SpatialGBM.plotPops(s04_E57,cancerPops, ncol = 4, pt_size = 1.5)
SpatialGBM.plotPops(s04_E57,TAMsPops, pt_size = 1.5)


## 2c: correlation heatmap

large_correlation = cor(All_RCTD_large[,3:20], method = "pearson")

corrHM_large = Heatmap(matrix = large_correlation,
                       col = colorRamp2(c(-1,-.25,0,.25,1),c("#2166ac","#92c5de","#f7f7f7","#f4a582","#b2182b")),
                       show_column_names = TRUE, show_row_names = TRUE, name = "Correlation",column_split = 3, row_split = 3,
                       cluster_rows = TRUE, cluster_columns = TRUE)
draw(corrHM_large)



##################################################################################################################################################################################################################################################


############
# Figure 3 #
############


## 3a: UMAP

NicheColorCode <- c("#cc79a7","#ee5e00","#e69f00","#f0e442","#009e73","#56b4e9")
names(NicheColorCode) <- paste0("N",1:6)

umap_Niche = ggplot(All_RCTD_large[sample(1:nrow(All_RCTD_large),replace=FALSE),],aes(x = umap.x, y = umap.y, colour = Niche)) +
  geom_point(size = .7) +
  scale_color_manual(values = NicheColorCode) +
  guides(colour = guide_legend(override.aes = list(size=5))) +
  labs(colour = "Niche", x = "UMAP1", y = "UMAP2") +
  theme_classic(base_size = 30)
umap_Niche


## 3c-d surface plots of niches

samplesNames = c("s08_E53","s16_E53","s03_E57","s04_E57","s09_E34","s10_E34","UKF259_T","UKF260_T","UKF248_T","UKF275_T","UKF255_T","UKF243_T")

for(s in samplesNames[]){
  
  print(s)
  
  # load data
  load(paste0("path/to/RData/files/",s,".RData"))
  
  # add nuclei count
  eval(parse(text = paste0(s,"<- addNucleiCount(",s,",'path/to/cell/counts",s,"/",s,"_detectedCells.txt')")))
  
}



for(s in samplesNames){
  
  clusterDf = All_RCTD_large[All_RCTD_large$sample == s, c("barcodes","Niche","tumour")]
  clusterDf$Niche = factor(clusterDf$Niche, levels = paste0("N",1:6))
  
  #eval(parse(text = paste0(s," <- SPATA2::updateSpataObject(",s,",method = 'VisiumSmall')")))
  eval(parse(text = paste0(s," <- addFeatures(",s,",feature_df = clusterDf, overwrite = TRUE)")))
  
}


plotSurface(s09_E34,"Niche", pt_size = 1.2, img_alpha = .5) + scale_color_manual(values = NicheColorCode)

p1 = plotSurface(s16_E53,"Niche", pt_size = 1.5, img_alpha = .5) + scale_color_manual(values = NicheColorCode, na.value = alpha("lightgrey",.5)) + plotSurface(s08_E53,"Niche", pt_size = 2, img_alpha = .5) + scale_color_manual(values = NicheColorCode, na.value = alpha("lightgrey",.5))
p2 = plotSurface(s03_E57,"Niche", pt_size = 1, img_alpha = .5) + scale_color_manual(values = NicheColorCode, na.value = alpha("lightgrey",.5)) + plotSurface(s04_E57,"Niche", pt_size = 1, img_alpha = .5) + scale_color_manual(values = NicheColorCode, na.value = alpha("lightgrey",.5))
p3 = plotSurface(s09_E34,"Niche", pt_size = 1, img_alpha = .5) + scale_color_manual(values = NicheColorCode, na.value = alpha("lightgrey",.5)) + plotSurface(s10_E34,"Niche", pt_size = 1, img_alpha = .5) + scale_color_manual(values = NicheColorCode, na.value = alpha("lightgrey",.5))

p4 = plotSurface(UKF260_T,"Niche", pt_size = 1.25, img_alpha = .5) + scale_color_manual(values = NicheColorCode, na.value = alpha("lightgrey",.5)) + plotSurface(UKF259_T,"Niche", pt_size = 1.25, img_alpha = .5) + scale_color_manual(values = NicheColorCode, na.value = alpha("lightgrey",.5))
p5 = plotSurface(UKF275_T,"Niche", pt_size = 1.25, img_alpha = .5) + scale_color_manual(values = NicheColorCode, na.value = alpha("lightgrey",.5)) + plotSurface(UKF243_T,"Niche", pt_size = 1, img_alpha = .5) + scale_color_manual(values = NicheColorCode, na.value = alpha("lightgrey",.5))
p6 = plotSurface(UKF248_T,"Niche", pt_size = 1, img_alpha = .5) + scale_color_manual(values = NicheColorCode, na.value = alpha("lightgrey",.5)) + plotSurface(UKF255_T,"Niche", pt_size = 1.5, img_alpha = .5) + scale_color_manual(values = NicheColorCode, na.value = alpha("lightgrey",.5))


(p1 / p2 / p3)
(p4 / p5 / p6)


##################################################################################################################################################################################################################################################


############
# Figure 4 #
############


## Fig. 4a: spatial transcriptomics: cell types

samplesNames = c("E24","E26")

for(s in samplesNames[]){
  
  print(s)
  
  # load data
  load(paste0("path/to/RData/files/",s,".RData"))
  
  # add nuclei count
  eval(parse(text = paste0(s,"<- addNucleiCount(",s,",'path/to/cell/counts",s,"/",s,"_detectedCells.txt')")))
  
}


SpatialGBM.plotPops(E24,cancerPops, ncol = 4, pt_size = 1.5) / SpatialGBM.plotPops(E24,TAMsPops, pt_size = 1.5)
SpatialGBM.plotPops(E26,cancerPops, ncol = 4, pt_size = 4) / SpatialGBM.plotPops(E26,TAMsPops, pt_size = 4)


## Fig. 4b: niches

centroids <- aggregate(All_RCTD_large[,3:20],by = All_RCTD_large$Niche, FUN = median)
centroids <- t(log10(1+centroids))


### predict by centroid (renormalised) on large cell types
largeE24 <- E24@meta_obs[,grep("RM_large",colnames(E24@meta_obs))] |> as.data.frame()
rownames(largeE24) <- E24@meta_obs$barcodes
colnames(largeE24) <- gsub("."," ",fp.strsplit(colnames(largeE24),"RCTD_",2),fixed=TRUE)
largeE24 <- log10(1+largeE24 * E24@meta_obs$nuclei.count) ## multiply by nuclei count and log-tranform

niche_E24 <- apply(largeE24,1,function(x){
  distances <- dist(rbind(x,t(centroids)),method = "manhattan") |> as.matrix()
  return(names(which.min(distances[-1,1])))
}) |> unlist()

niches_E24 <- data.frame(barcodes = names(niche_E24), Niche = niche_E24)
E24 <- addMetaDataObs(E24,niches_E24,overwrite = TRUE)

pE24 <- plotSurface(E24,"Niche", pt_size = 1.2, img_alpha = .5) + scale_color_manual(values = NicheColorCode, na.value = alpha("lightgrey",.5))


### predict by centroid (renormalised) on large cell types
largeE26 <- E26@meta_obs[,grep("RM_large",colnames(E26@meta_obs))] |> as.data.frame()
rownames(largeE26) <- E26@meta_obs$barcodes
colnames(largeE26) <- gsub("."," ",fp.strsplit(colnames(largeE26),"RCTD_",2),fixed=TRUE)
largeE26 <- log10(1+largeE26 * E26@meta_obs$nuclei.count) ## multiply by nuclei count and log-tranform

niche_E26 <- apply(largeE26,1,function(x){
  distances <- dist(rbind(x,t(centroids)),method = "manhattan") |> as.matrix()
  return(names(which.min(distances[-1,1])))
}) |> unlist()

niches_E26 <- data.frame(barcodes = names(niche_E26), Niche = niche_E26)
E26 <- addMetaDataObs(E26,niches_E26,overwrite = TRUE)

pE26 <- plotSurface(E26,"Niche", pt_size = 1.2, img_alpha = .5) + scale_color_manual(values = NicheColorCode, na.value = alpha("lightgrey",.5))



## Fig. 4d

cellData <- data.frame() # loaded the dataset containing, for each cell of each sample, with a column tumour with sample name, and a column Markers with all markers each cell is positive for (space-separated). Columns Centroid.X.µm and Centroid.Y.µm contain the cell coordinates. Column Cell.Area.µm.2 contains each cell's area.

cellData$avg.diameter <- 2*sqrt(cellData$Cell.Area.µm.2/pi)

markersPos <- strsplit(cellData$Markers, " ", fixed=TRUE)

for(marker in c("CD14","CD44","IBA1","TMEM119","SOX2")){
  for(i in 1:nrow(cellData)){
    cellData[i,marker] = marker %in% markersPos[[i]]
  }
}


# Cell typing
cellData$cell.type <- "Uncharacterised"

cellData$cell.type <- ifelse(cellData$SOX2 & cellData$CD44, "Cancer - MES-like", cellData$cell.type)
cellData$cell.type <- ifelse(cellData$SOX2 & !cellData$CD44, "Cancer - non-MES-like", cellData$cell.type)
cellData$cell.type <- ifelse(!cellData$SOX2 & (cellData$IBA1 | cellData$CD14 | cellData$TMEM119) & cellData$TMEM119, "TAM-MG", cellData$cell.type)
cellData$cell.type <- ifelse(!cellData$SOX2 & (cellData$IBA1 | cellData$CD14) & !cellData$TMEM119, "TAM-BDM", cellData$cell.type)

# Plot cell types
colorCode <- c("Cancer - MES-like" = "#e31a1c", "Cancer - non-MES-like" = "#fb9a99", "TAM-BDM" = "#33a02c", "TAM-MG" = "#b2df8a", "Uncharacterised" = "lightgrey")

cellData <- rbind(subset(cellData,cell.type == "Uncharacterised"),subset(cellData,cell.type != "Uncharacterised")) ## Placing uncharactarised cels at the top of the table will bring them to the bottom of the plot to make it clearer.

ggplot(cellData,aes(x = Centroid.X.µm, y = Centroid.Y.µm, colour = cell.type, size = avg.diameter)) +
  rasterise(
    geom_point(), 
    dpi = 300
  ) +
  facet_wrap(~tumour) +
  scale_size_continuous(range = c(0, 1)) +
  scale_colour_manual(values = colorCode) +
  theme_classic(base_size = 16)



## Fig. 4e: distances

# Initialize new columns for the results with NA values
cellData$dist_to_nearest_MES <- NA
cellData$dist_to_nearest_nonMES <- NA

for (samp in unique(cellData$tumour)) {
  # Subset data for the current sample to avoid repeated filtering
  sample_data <- cellData[cellData$tumour == samp, ]
  
  # Identify indices for the target cell type (Cancer - MES-like)
  idx_mac <- which(substr(sample_data$cell.type,1,3) == "TAM")
  
  # Proceed only if MES-like cells exist in this sample
  if (length(idx_mac) > 0) {
    
    # Extract coordinates for MES-like cells
    mac_x <- sample_data$Centroid.X.µm[idx_mac]
    mac_y <- sample_data$Centroid.Y.µm[idx_mac]
    
    # --- Calculation for TAM-BDM ---
    idx_cancer_mes <- which(sample_data$cell.type == "Cancer - MES-like")
    
    if (length(idx_cancer_mes) > 0) {
      cancer_mes_x <- sample_data$Centroid.X.µm[idx_cancer_mes]
      cancer_mes_y <- sample_data$Centroid.Y.µm[idx_cancer_mes]
      
      # Calculate distance matrix (rows: MES, cols: TAM-BDM)
      diff_x_mes <- outer(mac_x, cancer_mes_x, "-")
      diff_y_mes <- outer(mac_y, cancer_mes_y, "-")
      dist_matrix_mes <- sqrt(diff_x_mes^2 + diff_y_mes^2)
      
      # Find minimum distance for each MES cell
      min_dists_mes <- apply(dist_matrix_mes, 1, min)
      
      # Map back to original row indices
      original_row_indices_mac <- which(cellData$tumour == samp)[idx_mac]
      cellData$dist_to_nearest_MES[original_row_indices_mac] <- min_dists_mes
    }
    
    # --- Calculation for TAM-MG ---
    idx_cancer_nmes <- which(sample_data$cell.type == "Cancer - non-MES-like")
    
    if (length(idx_cancer_nmes) > 0) {
      cancer_nmes_x <- sample_data$Centroid.X.µm[idx_cancer_nmes]
      cancer_nmes_y <- sample_data$Centroid.Y.µm[idx_cancer_nmes]
      
      # Calculate distance matrix (rows: MES, cols: TAM-MG)
      diff_x_nmes <- outer(mac_x, cancer_nmes_x, "-")
      diff_y_nmes <- outer(mac_y, cancer_nmes_y, "-")
      dist_matrix_nmes <- sqrt(diff_x_nmes^2 + diff_y_nmes^2)
      
      # Find minimum distance for each MES cell
      min_dists_nmes <- apply(dist_matrix_nmes, 1, min)
      
      # Map back to original row indices (same MES indices as above)
      original_row_indices_mac <- which(cellData$tumour == samp)[idx_mac]
      cellData$dist_to_nearest_nonMES[original_row_indices_mac] <- min_dists_nmes
    }
  }
}

p_MES <- ggplot(subset(cellData, substr(cell.type,1,3) == "TAM" & tumour !="E20"), aes(x = cell.type, y = dist_to_nearest_MES, fill = cell.type)) +
  geom_violin() +
  geom_hline(yintercept = 40, lty=2) +
  geom_boxplot(fill="white",alpha = 0.5, width = .3) +
  #facet_wrap(~tumour) +
  #geom_jitter(width = .1, alpha = 0.5) +
  theme(legend.position = "none") +
  labs(title = "Distance to neareset MES cancer", y = "Distance to nearest MES (µm)", x = NULL) +
  ggpubr::stat_pwc(method = "dunn.test",p.adjust.method = "BH",label = "p.format", hide.ns = FALSE,p.digits = 3) +
  scale_y_log10() +
  theme_bw(base_size = 14) +
  theme(legend.position = "none")

p_nonMES <- ggplot(subset(cellData, substr(cell.type,1,3) == "TAM" & tumour !="E20"), aes(x = cell.type, y = dist_to_nearest_nonMES, fill = cell.type)) +
  geom_violin() +
  geom_hline(yintercept = 40, lty=2) +
  geom_boxplot(fill="white",alpha = 0.5, width = .3) +
  #facet_wrap(~tumour) +
  #geom_jitter(width = .1, alpha = 0.5) +
  theme(legend.position = "none") +
  labs(title = "Distance to neareset non-MES cancer", y = "Distance to nearest non-MES (µm)", x = NULL) +
  ggpubr::stat_pwc(method = "dunn.test",p.adjust.method = "BH",label = "p.format", hide.ns = FALSE,p.digits = 3) +
  scale_y_log10() +
  theme_bw(base_size = 14) +
  theme(legend.position = "none")


##################################################################################################################################################################################################################################################


############
# Figure 5 #
############


# seuratGBM contains a seurat object initialised with the concatenated counts and metadata of all samples

Idents(seuratGBM) <- meta$seuratGBM$Niche
markerGenes <- FindAllMarkers(seuratGBM)
seuratGBM <- NormalizeData(seuratGBM,normalization.method = "LogNormalize",verbose = TRUE)
seuratGBM <- subset(seuratGBM, downsample = 2000)

# volcanoplots
volcanoPlotsList <- lapply(paste0("N",1:6), function(n){
  niche_markers <- FindMarkers(seuratGBM,ident.1 = n)
  p <- EnhancedVolcano::EnhancedVolcano(niche_markers,x = "avg_log2FC",y = "p_val_adj",lab = rownames(niche_markers),title = n,subtitle = NULL,
                                 drawConnectors = TRUE,typeConnectors = "open", arrowheads = FALSE)
  return(p)})
names(volcanoPlotsList) <- paste0("N",1:6)

# enrichment plots
enrichmentPosPlots <- lapply(paste0("N",1:6),function(n){
  enrichmentPlot(subset(markerGenes,Niche == "n" & avg_log2FC > 0)[,"gene"],"path/to/MSigDB/Hallmarks/gmt/file")
})
enrichmentNegPlots <- lapply(paste0("N",1:6),function(n){
  enrichmentPlot(subset(markerGenes,Niche == "n" & avg_log2FC < 0)[,"gene"],"path/to/MSigDB/Hallmarks/gmt/file")
})



##################################################################################################################################################################################################################################################



############
# Figure 6 #
############

## This part requires the estimates of deconvolution results, and the TCGA annotation data (see script 05) 


## Fig 6a: DWLS in TCGA
DWLS_results <- readRDS("/path/to/DWLS/results") ## output from decenvolution script

DWLS_results <- DWLS_results[order(DWLS_results[,"N1"]),order(colnames(DWLS_results))] # reorder by estimate for N1

barplot(t(DWLS_results),horiz = TRUE, col = NicheColorCode,border = FALSE)



## Fig 6b: forestplot

DWLS_results <- DWLS_results[rownames(annot_TCGA_GBM),]
annot_TCGA_GBM <- cbind(annot_TCGA_GBM,DWLS_results)

censorTime <- 12
annot_TCGA_GBM$OS.event.censored <- ifelse(annot_TCGA_GBM$OS.time > censorTime, 0, annot_TCGA_GBM$OS.event)
annot_TCGA_GBM$OS.time.censored <- ifelse(annot_TCGA_GBM$OS.time > censorTime, censorTime, annot_TCGA_GBM$OS.time)



coxphTable <- function(df, eventName, timeName, variables){
  
  
  event <- df[,eventName]
  time <- df[,timeName]
  surv <- Surv(time = time, event = event)
  
  res <- data.frame(variable = "variables", level = NA, n = NA, n.events = NA, HR = NA, CI.95.low = NA, CI.95.high = NA, p.value = NA)[-1,]
  
  for(var in variables){
    
    res.loc <- data.frame(variable = var, level = NA, n = NA, n.events = NA, HR = NA, CI.95.low = NA, CI.95.high = NA, p.value = NA)
    
    cox <- coxph(surv ~ df[,var])
    
    res.loc[,"n"] <- cox$n
    res.loc[,"n.events"] <- cox$nevent
    res.loc[,"HR"] <- exp(cox$coefficients[1])
    res.loc[,"CI.95.low"] <- exp(confint(cox,level = .95))[1]
    res.loc[,"CI.95.high"] <- exp(confint(cox,level = .95))[2]
    res.loc$p.value <- summary(cox)$coefficients[5]
    
    res <- rbind(res,res.loc)
    
  }
  
  return(res)
  
}


coxTable = coxphTable(annot = annot_TCGA_GBM,event = "OS.event.censored",delay = "OS.time.censored",variables = paste0("N",1:6))

forestTable = data.frame(annotation = coxTable$variable,
                         mean = coxTable$HR,
                         lower = coxTable$CI.95.low,
                         upper = coxTable$CI.95.high,
                         HR_CI = paste0(coxTable$H.R.," (",coxTable$CI.95.low," - ",coxTable$CI.95.high,")"),
                         p_value = coxTable$p.value)

forestTable |>
  forestplot(labeltext = c(annotation, HR_CI, p_value),
             xlog = TRUE, boxsize = .2,lwd.ci = 5) |>
  fp_add_header(annotation = c("", "Study"),
                HR_CI = c("", "Hazard Ratio (95% CI)"),
                p_value = c("", "p-value"))



## Fig. 6c-d

discretize <- function(x,n){
  br <- quantile(x, probs = seq(0, 1, length.out = n + 1), na.rm = TRUE)
  return(as.integer(cut(x,br, labels=FALSE, include.lowest=TRUE)))
}

for(niche in paste0("N",1:6)){
  annot_TCGA_GBM[,paste(niche,"Disc",sep=".")] = paste0("Q",discretize(annot_TCGA_GBM[,niche],4))
}



kmPlot <- function(time,event, stratification, colors){
  
  #survival model
  surv = Surv(time,event)
  fit = survfit(surv~stratification)
  diff = survdiff(surv~stratification)
  
  #p-value for difference between survival curves
  p = signif(pchisq(diff$chisq,length(diff$obs)-1,lower.tail=F), digits=3)
  
  print(paste0("p = ",p))
  
  #plot
  plot(fit, col=colors)
  
  #plot legend
  n = fit$n
  strata = names(fit$strata)
  labels = c()
  for(i in 1:length(n)){
    namesplit = strsplit(strata[i],"=")
    name = namesplit[[1]][length(namesplit[[1]])]
    labels[i] = paste(name, " (n=",n[i],")",sep="")
  }
  legend(10,0.4,legend=labels, lty=1, col=colors, bty="n")
  
}


for(niche in paste0("N",c(1,5))){
  kmPlot(time = annot_TCGA_GBM$OS.time.censored, event = annot_TCGA_GBM$OS.event.censored,stratification = annot_TCGA_GBM[,paste(niche,"Disc",sep=".")], colors = rev(brewer.pal(4,"RdBu")))
}





##################################################################################################################################################################################################################################################



############
# Figure 7 #
############


## This section assumes access to the deconvolution results, and that annot_PRJNA482620 is loaded (see script 06)

deconvolutionResults <- readRDS("dwls_results_PRJA482620.rds")


deconvolutionResults <- deconvolutionResults[rownames(annot_PRJNA482620),]
for(niche in paste0("N",1:6)){
  annot_PRJNA482620[,niche] <- deconvolutionResults[,niche]
}

annot_PRJNA482620 <- annot_PRJNA482620[!is.na(annot_PRJNA482620$PATIENT_ID),]

## Fig. 7a

nichePlots <- lapply(paste0("N",1:6),function(niche){
  p <- ggplot(annot_PRJNA482620[annot_PRJNA482620$Immediately_Pre,], aes(x = `Response to PD1?`, y = annot_PRJNA482620[annot_PRJNA482620$Immediately_Pre,niche]), fill = annot_PRJNA482620[annot_PRJNA482620$Immediately_Pre,"Response to PD1?"]) +
    geom_violin() + scale_fill_manual(values = c("#fb9a99","#b2df8a")) +
    geom_boxplot(width = 0.25, fill = "white", alpha = .7) +
    geom_jitter(colour = "black", width = 0.05, alpha = .3) +
    ggpubr::stat_pwc(method = "wilcox_test") +
    labs(y = paste0(niche," (pretreatment)")) +
    theme(legend.position = "none") +
    theme_bw()
  return(p)
})

do.call("grid.arrange", c(nichePlots, ncol=3))


## Fig. 7b

nichePlots_evolution <- lapply(paste0("N",1:6),function(niche){
  p <- ggpaired(annot_PRJNA482620[annot_PRJNA482620$Pre_Post_aPD1=="Post" | annot_PRJNA482620$Immediately_Pre,], x = "Pre_Post_aPD1", y = niche, id = "PATIENT_ID",
                facet.by = "Response to PD1?",
                color = "Pre_Post_aPD1", line.color = "gray", line.size = 0.4,ggtheme = theme_bw(), palette = brewer.pal(4,"Set1")) +
    labs(x = "Time point", y = niche) +
    theme(legend.position="none") +
    stat_pwc(method="wilcox_test")
  return(p)
})

do.call("grid.arrange", c(nichePlots_evolution, ncol=3))






