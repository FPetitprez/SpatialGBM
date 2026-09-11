#########################
# Supplementary figures #
#########################



## This script assumes the All_RCTD_large and All_RCTD_largeobject are loaded (see script 01 to 04)
## Only items differing or not covered in the main figures script are detailed here. 


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
library("future")
library("CellChat")
library("forestplot")






##########
# Fig S1 #
##########

tumourList <- c("E31","E34","E37","E51","E53","E55","E57")
samplesNames <- unique(All_RCTD_large$sample)


# Samples Color Code
samplesColorCode <- c("s06-E31" = "#ffff33", "s14-E31" = "#ffff99",
                     "s09-E34" = "#555555","s10-E34" = "#a0a0a0", "s11-E37" = "#ff7f00", "s12-E37" = "#fdbf6f",
                     "s07-E51" = "#e31a1c", "s15-E51" = "#fb9a99", "s08-E53" = "#33a02c", "s16-E53" = "#b2df8a",
                     "s01-E55" = "#1f78b4", "s02-E55" = "#a6cee3", "s03-E57" = "#6a3d9a", "s04-E57" = "#cab2d6",
                     "UKF242_T" = "#1b9e77", "UKF243_T" = "#8dd3c7", "UKF248_T" = "#e7298a",
                     "UKF255_T" = "#002d9c", "UKF259_T" = "#00bdff", "UKF260_T" = "#0b84a5", "UKF275_T" = "#9f1853",
                     "UKF296_T" = "black", "UKF304_T" = "#423300", "UKF313_T" = "#005700", "UKF334_T" = "#858500")




KSpvals_corr <- lapply(colnames(All_RCTD_large)[3:20],function(pop){
  KSpvals <- lapply(tumourList,function(tum){
    samples <- samplesNames[grep(tum,samplesNames)]
    ks <- ks.test(All_RCTD_large[All_RCTD_large$sample==samples[1],pop],All_RCTD_large[All_RCTD_large$sample==samples[2],pop])
    return(ks$p.value)
  }) |> unlist()
  
  names(KSpvals) <- tumourList
  return(KSpvals)})

allKS <- unlist(KSpvals_corr) |> p.adjust(method = "bonferroni")
for(i in 1:18){
  KSpvals_corr[[i]] <- allKS[((i-1)*length(tumourList)+1):(i*length(tumourList))]
}
names(KSpvals_corr) <- colnames(All_RCTD_large)[3:20]


for(pop in c("MES like","TAM MG")){
  distribList = lapply(samplesNames[substr(samplesNames,1,1)=="s"], function(s){
    textAnnot <- grobTree(textGrob(paste0("KS p=",signif(KSpvals_corr[[pop]][fp.strsplit(s,"-",2)],2)), x=0.95,  y=0.9, hjust=1))
    
    popName = sym(pop)
    
    p = ggplot(All_RCTD_large[All_RCTD_large$sample==s,], aes(x = !!popName)) +
      geom_histogram(fill = samplesColorCode[gsub("-","_",s,fixed=TRUE)]) +
      labs(title = s) +
      theme_minimal() +
      theme(legend.position = "none")
    
    return(p+annotation_custom(textAnnot))
  })
  
  do.call("grid.arrange", c(distribList, ncol=2))
  
}


##################################################################################################################################################################################################################################################

###########
# Fig. S3 #
###########


umap_Sample = ggplot(All_RCTD_large[sample(1:nrow(All_RCTD_large),replace=FALSE),],aes(x = umap.x, y = umap.y, colour = sample)) +
  geom_point(size = .7) +
  scale_color_manual(values = samplesColorCode) +
  guides(colour = guide_legend(override.aes = list(size=5))) +
  labs(colour = "Sample", x = "UMAP1", y = "UMAP2") +
  theme_classic(base_size = 30)
umap_Sample



umap_CellCount = ggplot(All_RCTD_large[sample(1:nrow(All_RCTD_large),replace=FALSE),],aes(x = umap.x, y = umap.y, colour = CellCount)) +
  ggrastr::rasterise(geom_point(size = .7)) +
  scale_color_viridis_c(option = "turbo",limits = c(0,25)) +
  labs(colour = "Cell count", x = "UMAP1", y = "UMAP2") +
  theme_classic(base_size = 20)
umap_CellCount



umap_MESlike = ggplot(All_RCTD_large[sample(1:nrow(All_RCTD_large),replace=FALSE),],aes(x = umap.x, y = umap.y, colour = `MES like`)) +
  ggrastr::rasterise(geom_point(size = .5)) +
  scale_color_viridis_c(option = "turbo",limits = c(0,15)) +
  #guides(colour = guide_legend(override.aes = list(size=5))) +
  labs(colour = "MES-like cells", x = "UMAP1", y = "UMAP2") +
  theme_classic(base_size = 15)

umap_AClike = ggplot(All_RCTD_large[sample(1:nrow(All_RCTD_large),replace=FALSE),],aes(x = umap.x, y = umap.y, colour = `AC like`)) +
  ggrastr::rasterise(geom_point(size = .5)) +
  scale_color_viridis_c(option = "turbo",limits = c(0,15)) +
  #guides(colour = guide_legend(override.aes = list(size=5))) +
  labs(colour = "AC-like cells", x = "UMAP1", y = "UMAP2") +
  theme_classic(base_size = 15)

umap_OPClike = ggplot(All_RCTD_large[sample(1:nrow(All_RCTD_large),replace=FALSE),],aes(x = umap.x, y = umap.y, colour = `OPC like`)) +
  ggrastr::rasterise(geom_point(size = .5)) +
  scale_color_viridis_c(option = "turbo",limits = c(0,5)) +
  #guides(colour = guide_legend(override.aes = list(size=5))) +
  labs(colour = "OPC-like cells", x = "UMAP1", y = "UMAP2") +
  theme_classic(base_size = 15)

umap_NPClike = ggplot(All_RCTD_large[sample(1:nrow(All_RCTD_large),replace=FALSE),],aes(x = umap.x, y = umap.y, colour = `NPC like`)) +
  ggrastr::rasterise(geom_point(size = .5)) +
  scale_color_viridis_c(option = "turbo",limits = c(0,5)) +
  #guides(colour = guide_legend(override.aes = list(size=5))) +
  labs(colour = "NPC-like cells", x = "UMAP1", y = "UMAP2") +
  theme_classic(base_size = 15)

(umap_MESlike + umap_AClike) / (umap_OPClike + umap_NPClike)


umap_TAM_BDM = ggplot(All_RCTD_large[sample(1:nrow(All_RCTD_large),replace=FALSE),],aes(x = umap.x, y = umap.y, colour = `TAM BDM`)) +
  ggrastr::rasterise(geom_point(size = .5)) +
  scale_color_viridis_c(option = "turbo",limits = c(0,3)) +
  #guides(colour = guide_legend(override.aes = list(size=5))) +
  labs(colour = "TAM BDM", x = "UMAP1", y = "UMAP2") +
  theme_classic(base_size = 15)

umap_TAM_MG = ggplot(All_RCTD_large[sample(1:nrow(All_RCTD_large),replace=FALSE),],aes(x = umap.x, y = umap.y, colour = `TAM MG`)) +
  ggrastr::rasterise(geom_point(size = .5)) +
  scale_color_viridis_c(option = "turbo",limits = c(0,3)) +
  #guides(colour = guide_legend(override.aes = list(size=5))) +
  labs(colour = "TAM MG", x = "UMAP1", y = "UMAP2") +
  theme_classic(base_size = 15)

umap_TAM_MG + umap_TAM_BDM





##################################################################################################################################################################################################################################################

###########
# Fig. S8 #
###########


## Fig. S8a-c

samplesNames <- c() #put in all available spatial transcriptomics samples

s = samplesNames[1]
# load data
load(paste0("/path/to/RData/files/",s,".RData"))
eval(parse(text = paste0("spata2object <- ",s)))
# spata2object <- updateSpataObject(spata2object, method = "VisiumSmall") # uncomment if RData generated with an older version of SPATA2
eval(parse(text = paste0("rm(",s,")")))
locCounts <- getCountMatrix(spata2object)

locCoords <- getCoordsDf(spata2object)
locMeta <- data.frame(barcode = colnames(locCounts), sample = s, x = locCoords$x, y = locCoords$y)

spatialGBMcounts = locCounts
meta = locMeta

for(s in samplesNames[2:length(samplesNames)]){
  print(s)
  load(paste0("/path/to/RData/files/",s,".RData"))
  eval(parse(text = paste0("spata2object <- ",s)))
  spata2object <- updateSpataObject(spata2object, method = "VisiumSmall")
  eval(parse(text = paste0("rm(",s,")")))
  locCounts <- getCountMatrix(spata2object)
  locCoords <- getCoordsDf(spata2object)
  locMeta <- data.frame(barcode = colnames(locCounts), sample = s, x = locCoords$x, y = locCoords$y)
  
  commonGenes <- intersect(rownames(spatialGBMcounts), rownames(locCounts))
  spatialGBMcounts <- cbind(spatialGBMcounts[commonGenes,],locCounts[commonGenes,])
  meta <- rbind(meta,locMeta)
  gc()
}

identical(colnames(spatialGBMcounts),meta$barcode)

meta <- as.data.frame(meta)
rownames(meta) <- paste(meta$sample,meta$barcode,sep = "_")
colnames(spatialGBMcounts) <- rownames(meta)

meta$Niche <- All_RCTD_fine[rownames(meta),"Niche"]
meta <- meta[!is.na(meta$Niche),]
spatialGBMcounts <- spatialGBMcounts[,rownames(meta)]

spatialGBMmetadata <- meta
colnames(spatialGBMmetadata) <- c("barcode","samples","x","y","Niche")
meta$labels <- meta$Niche

normCounts_SpatialGBM <- CellChat::normalizeData(spatialGBMcounts)
beepr::beep()


## prepare scale factors
spatial.factors <- data.frame(samples = c(), ratio = c(), tol = c())

for(s in samplesNames){
  folder <- paste0("/path/to/preprocessed/data/",s,"/outs/spatial/")
  sf <- jsonlite::fromJSON(txt = file.path(folder, 'scalefactors_json.json'))
  spot.size <- 55
  conversion.factor <- spot.size/sf$spot_diameter_fullres
  spatial.factors.line <- data.frame(samples = s, ratio = conversion.factor, tol = spot.size/2)
  spatial.factors <- rbind(spatial.factors,spatial.factors.line)
}

rownames(spatial.factors) <- gsub("-","_",spatial.factors$samples,fixed=TRUE)
spatial.factors <- spatial.factors[samplesNames,]

colnames(meta)[grep("sample",colnames(meta))] <- "samples"



## create cellchat object
cellchat <- createCellChat(object = normCounts_SpatialGBM, meta = meta, group.by = "labels", datatype = "spatial", coordinates = meta[,c("x","y")],spatial.factors = spatial.factors)
CellChatDB <- CellChatDB.human
cellchat@DB = CellChatDB
cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
future::plan("multisession", workers = 4) 
cellchat <- identifyOverExpressedGenes(cellchat)

options(future.globals.maxSize = 8000 * 1024^2)
cellchat <- identifyOverExpressedInteractions(cellchat)


## Inference of cell-cell communication network
ptm = Sys.time()
cellchat <- computeCommunProb(cellchat, type = "truncatedMean", trim = 0.1, 
                              distance.use = FALSE, interaction.range = 250, scale.distance = NULL,
                              contact.dependent = TRUE, contact.range = 100)

cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)

execution.time = Sys.time() - ptm
print(as.numeric(execution.time, units = "secs"))

for(niche in paste0("N",c(1,3,5))){
  netVisual_bubble(cellchat,sources.use = niche, targets.use = niche)
}



