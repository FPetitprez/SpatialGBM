#######################
# UMAP and clustering #
#######################

## This script asusmes that All_RCTD_fine is loaded (script 03)


set.seed(42)
library("factoextra")
library("cluster")


####################
# UMAP coordinates #
####################


logCellCounts <- log10(1+All_RCTD_fine[,3:35])

umapCoordinates <- umap::umap(logCellCounts,)


coords <- umapCoordinates$layout
All_RCTD_fine$umap.x <- coords[,1]
All_RCTD_fine$umap.y <- coords[,2]

All_RCTD_large$umap.x <- All_RCTD_fine$umap.x
All_RCTD_large$umap.y <- All_RCTD_fine$umap.y



######################
# k means clustering #
######################

method1 <- fviz_nbclust(logCellCounts, method = "silhouette", k.max = 15, verbose = TRUE, FUNcluster = kmeans, print.summary = TRUE)
method2 <- fviz_nbclust(logCellCounts, pam, method = "gap_stat", k.max = 15,nboot = 20,print.summary = TRUE)

k6 <- cluster::pam(logCellCounts, k = 6)

All_RCTD_fine$cluster <- paste0("C",k6$clustering)
All_RCTD_fine$Niche <- c("C1"="N4", "C2"="N6","C3"="N3","C4"="N2","C5"="N1","C6"="N5")[All_RCTD_fine$cluster] # renumbering the clusters into niches

All_RCTD_large$cluster <- paste0("C",k6$clustering)
All_RCTD_large$Niche <- c("C1"="N4", "C2"="N6","C3"="N3","C4"="N2","C5"="N1","C6"="N5")[All_RCTD_large$cluster] # renumbering the clusters into niches

