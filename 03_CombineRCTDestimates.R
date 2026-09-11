##############################
# Combine all RCTD estimates #
##############################



samplesNames <- c() ## char of all samples names


s <- samplesNames[1]



# load data
load(paste0("path/to/RData/files/",s,".RData"))


# add nuclei count
eval(parse(text = paste0(s,"<- addNucleiCount(",s,",'../../2023_04_cellDetection/",s,"/",s,"_detectedCells.txt')")))

spata2object <- eval(parse(text = paste0("spata2object <- ",s)))


large_estimates <- getFeatureNames(spata2object)[grep("RM_large_RCTD_",getFeatureNames(spata2object))]
names(large_estimates) <- gsub("."," ",fp.strsplit(large_estimates,"_",4),fixed = TRUE)

nucleiCount <- SPATA2::getFeatureVariables(spata2object,features = "nuclei.count")

All_RCTD_large <- getFeatureVariables(spata2object,large_estimates)
All_RCTD_large[3:ncol(All_RCTD_large)] = apply(All_RCTD_large[3:ncol(All_RCTD_large)],2,function(x){x*nucleiCount$nuclei.count})
colnames(All_RCTD_large) = c(colnames(All_RCTD_large)[1:2],names(large_estimates))

All_RCTD_fine <- getFeatureVariables(spata2object,fine_estimates)
All_RCTD_fine[3:ncol(All_RCTD_fine)] = apply(All_RCTD_fine[3:ncol(All_RCTD_fine)],2,function(x){x*nucleiCount$nuclei.count})
colnames(All_RCTD_fine) = c(colnames(All_RCTD_fine)[1:2],names(fine_estimates))


## loop on all other samples

for(s in samplesNames[2:length(samplesNames)]){
  
  print(s)
  
  # load data
  load(paste0("~/Documents/UoE/2020_10_spatialGBM/2022_scRNAseqAtlas/Ruiz_Moreno_2022/RData files/",s,".RData"))
  
  # add nuclei count
  eval(parse(text = paste0(s,"<- addNucleiCount(",s,",'../../2023_04_cellDetection/",s,"/",s,"_detectedCells.txt')")))
  
  eval(parse(text = paste0("spata2object <- ",s)))
  
  nucleiCount <- SPATA2::getFeatureVariables(spata2object,features = "nuclei.count")
  
  Loc_RCTD_large <- getFeatureVariables(spata2object,large_estimates)
  Loc_RCTD_large[3:ncol(Loc_RCTD_large)] = apply(Loc_RCTD_large[3:ncol(Loc_RCTD_large)],2,function(x){x*nucleiCount$nuclei.count})
  colnames(Loc_RCTD_large) <- c(colnames(Loc_RCTD_large)[1:2],names(large_estimates))
  
  Loc_RCTD_fine <- getFeatureVariables(spata2object,fine_estimates)
  Loc_RCTD_fine[3:ncol(Loc_RCTD_fine)] = apply(Loc_RCTD_fine[3:ncol(Loc_RCTD_fine)],2,function(x){x*nucleiCount$nuclei.count})
  colnames(Loc_RCTD_fine) <- c(colnames(Loc_RCTD_fine)[1:2],names(fine_estimates))
  
  All_RCTD_large <- rbind(All_RCTD_large, Loc_RCTD_large)
  All_RCTD_fine <- rbind(All_RCTD_fine, Loc_RCTD_fine)
  
  eval(parse(text = paste0("rm(",s,")")))
  
}

All_RCTD_fine = All_RCTD_fine[!is.na(All_RCTD_fine$Vascular),]
All_RCTD_large = All_RCTD_large[!is.na(All_RCTD_large$Vascular),]

#remove empty spots
All_RCTD_fine = All_RCTD_fine[apply(All_RCTD_fine[,3:35],1,sum)>0,]
All_RCTD_large = All_RCTD_large[apply(All_RCTD_large[,3:20],1,sum)>0,]

All_RCTD_large$CellCount = apply(All_RCTD_large[,3:20],1,sum)
All_RCTD_fine$CellCount = apply(All_RCTD_large[,3:20],1,sum)

