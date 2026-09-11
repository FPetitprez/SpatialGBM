########################################
# SPATA processing for Visium GBM data #
########################################


library(spacexr)
library(Matrix)
library(doParallel)
library(Seurat)
library(ggplot2)
library(SPATA2)
library(tictoc)
library(colorspace)




########################
# Processing functions #
########################




processSPATA2object = function(sampleName = "s15_E51", samplePath = "path/to/sample/files"){
  
  date = fp.strsplit(as.character(Sys.time())," ")
  
  # Open a file to send messages to
  zz <- file(paste(date,sampleName,"SPATA2_log.Rout",sep="_"), open = "wt")
  # Divert messages to that file
  sink(zz, type = "message")
  
  tic(paste("Done processing sample",sampleName))
  
  ## Step 1: load data
  
  spata2object = initiateSpataObject_10X(samplePath, sampleName)
  
  
  ## Step 2: autoencoder
  
  tic("Autoencoder")
  
  spata2object = runAutoencoderAssessment(spata2object,
                                          activations = c("relu", "selu", "sigmoid"), 
                                          bottlenecks = c(32, 40, 48, 56, 64),
                                          epochs = 20, 
                                          layers = c(128, 64, 32), 
                                          dropout = 0.1)
  
  assessment = getAutoencoderAssessment(spata2object)$df
  assessment = assessment[order(assessment$total_var,decreasing = TRUE),]
  opt.activation = assessment[1,"activation"]
  opt.bottleneck = as.numeric(as.character(assessment[1,"bottleneck"]))
  
  spata2object = runAutoencoderDenoising(
    object = spata2object, 
    activation = opt.activation, 
    bottleneck = opt.bottleneck, 
    epochs = 20, 
    layers = c(128, 64, 32), 
    dropout = 0.1
  )
  
  toc()
  
  ## Step 3: CNV
  
  tic("CNV")
  
  
  CNVdir = paste(sampleName,"_CNV",sep="")
  
  if (!paste(".",CNVdir,sep="/") %in% list.dirs()){system(paste("mkdir",CNVdir,sep=" "))}
  
  spata2object <- runCnvAnalysis(
    object = spata2object,
    directory_cnv_folder = CNVdir, # example directory
    cnv_prefix = "Chr"
  )
  
  toc()
  
  
  toc()
  
  return(spata2object)
  
}




######################
# Run on all samples #
######################

folders = c("path/to/data")

samplesNames = c() # char of all samples names

for(i in 1:length(samplesNames)){
  res <- processSPATA2object(sampleName = samplesNames[i], samplePath = folders[i])
  res <- runRCTDonSPATA2(res,samplesNames[i],folders[i])
  gdata::mv("res",gsub("-","_",samplesNames[i]))
  eval(parse(text = paste0("save(",samplesNames[i],",file = paste0('./RdataFiles/',",samplesNames[i],",'.RData'))")))
}
