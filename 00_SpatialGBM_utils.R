####################
# SpatialGBM utils #
####################

require("patchwork")
require("enrichplot")


SpatialGBM.plotRCTDestimates = function(spata2object, resolution = c("large","fine")[1], as.cell.counts = TRUE, pt_size = 2,display_image = TRUE, normalize = FALSE, smooth = FALSE, palette = "turbo", ...){
  
  if(resolution == "large"){estimates = getFeatureNames(spata2object)[grep("RM_large_RCTD_",getFeatureNames(spata2object))]}
  if(resolution == "fine"){estimates = getFeatureNames(spata2object)[grep("RM_fine_RCTD_",getFeatureNames(spata2object))]}
  names(estimates) = gsub("."," ",fp.strsplit(estimates,"_",4),fixed = TRUE)
  
  if(as.cell.counts){
    if(!"nuclei.count" %in% getFeatureNames(spata2object)){stop("Feature nuclei.count not found in SPATA2 object.")}
    nucleiCount = SPATA2::getFeatureVariables(spata2object,features = "nuclei.count")
    for(feat in estimates){
      res = getFeatureVariables(spata2object,feat)
      res[,feat] = res[,feat]*nucleiCount$nuclei.count
      res[which(res[,feat]==0),feat] <- NA
      spata2object = addFeatures(spata2object,res,feat, overwrite = T)
    }
  }
  
  for(i in 1:length(estimates)){
    eval(parse(text = paste0("spata2object <- SPATA2::renameFeatures(spata2object, '",names(estimates)[i], "' = '",estimates[i],"')")))
  }
  
  p = NULL
  
  p = plotSurfaceComparison(object = spata2object, color_by = names(estimates), pt_size = pt_size,display_image = display_image, normalize = normalize, smooth = smooth, ...)
  
  if((!is.null(p)) & palette == "turbo"){p = p+scale_color_viridis_c(option="turbo")}
  
  return(p)
  
}



SpatialGBM.plotAllRCTD = function(spta2object, ...){
  p1 = SpatialGBM.plotRCTDestimates(spta2object,"global", ncol = 10, ...)
  p2 = SpatialGBM.plotRCTDestimates(spta2object,"fine", ncol = 10, ...)
  p1 / p2 + patchwork::plot_layout(heights = c(2,3))
}

fp.strsplit = function(charVector,sep,index = 1){
  return(unlist(lapply(1:length(charVector),function(i){strsplit(charVector[i],sep,fixed=T)[[1]][index]})))
}



addNucleiCount = function(spata2object,detectedcellsPath){
  cellCounts = read.table(detectedcellsPath, sep = "\t", comment.char = "", header = TRUE)
  return(addFeatures(object = spata2object, feature_names = "nuclei.count", feature_df = data.frame(barcodes = cellCounts$Name, nuclei.count = cellCounts$Num.Detections)))
}



SpatialGBM.plotPops <- function(spata2object, populations, as.cell.counts = TRUE, pt_size = 2,display_image = TRUE, normalize = FALSE, smooth = FALSE, palette = "turbo", ...){
  
  names(populations) <- gsub("."," ",fp.strsplit(populations,"_",4),fixed = TRUE)
  
  if(as.cell.counts){
    if(!"nuclei.count" %in% getFeatureNames(spata2object)){stop("Feature nuclei.count not found in SPATA2 object.")}
    nucleiCount <- getMetaDf(spata2object)$nuclei.count
    for(feat in populations){
      res <- getMetaDf(spata2object)
      res[,feat] <- res[,feat]*nucleiCount
      res[which(res[,feat]==0),feat] <- NA
      spata2object <- addFeatures(spata2object,res,feat, overwrite = T)
    }
  }
  
  p <- NULL
  
  p <- plotSurfaceComparison(object = spata2object, color_by = populations, pt_size = pt_size,display_image = display_image, normalize = normalize, smooth = smooth, ...)
  
  if((!is.null(p)) & palette == "turbo"){p = p+scale_color_viridis_c(option="turbo")}
  
  return(p)
  
}


enrichmentPlot <- function(geneList, gmtFilePath = "path/to/hallmarks/gmt/file"){
  gmtFile <- clusterProfiler::read.gmt(gmtF)
  enrichment <- clusterProfiler::enricher(geneList,TERM2GENE = gmtFile)
  plot(dotplot(enrichment,showCategory = ncat,title = plotNames[coll])+labs(x = "Gene Ratio"))
  return(enrichmentList)
} 
