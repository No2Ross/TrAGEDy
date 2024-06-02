#Trade seq and STACAS do not get along 
library(reticulate)
use_python("/Users/rosslaidlaw/mambaforge/bin/python")
reticulate::import("phate")
library(gridExtra)
library(grid)
library(ggplot2)
library(lattice)
library(phateR)
library(dplyr)
library(Seurat)
library(rgl)
library(ggplot2)
library(ggpubr)
library(cowplot)
library(dtw)
library(scater)
library(scran)
library(clustree)
library(RColorBrewer)
library(SingleCellExperiment)
library(tradeSeq)
library(slingshot)
library(stats)
library(stringi)
library(reshape2)
library(stringr)
library(mclust)
setwd("~/R")
source("Scripts/methods/own_method_functions.R")

emma_obj <- readRDS("/Users/rosslaidlaw/R/TrAGEDy_V2/Tbrucei_Zc3h20_KO/objects/WT_ZC3H20_with_phate.rds")

emma_obj <- AddMetaData(emma_obj, emma_obj@active.ident, col.name = "cell_type")

split_output <- SplitObject(emma_obj, "orig.ident")

#Merge the A and B subsets into one
split_output$WT_01 <- RenameIdents(split_output$WT_01, "LS A.1" = "LS A", "LS A.2"= "LS A", "LS B.1"="LS B", "LS B.2" ="LS B")
split_output$WT_02 <- RenameIdents(split_output$WT_02, "LS A.1" = "LS A", "LS A.2"= "LS A", "LS B.1"="LS B", "LS B.2" ="LS B")
split_output$ZC3H20_KO <- RenameIdents(split_output$ZC3H20_KO, "LS A.1" = "LS A", "LS A.2"= "LS A", "LS B.1"="LS B", "LS B.2" ="LS B")


markers_list <- list()

DimPlot(split_output$WT_01, group.by = "seurat_clusters")

#Get feature space for alignment
features <- c()
for (i in 1:length(split_output)){
  
  current <- split_output[[i]]
  
  DefaultAssay(current) <- "RNA"
  
  markers_ <- FindAllMarkers(current, only.pos = T, logfc.threshold = 0.75, min.pct = 0.25)
  
  markers <- subset(markers_, p_val_adj < 0.05)
  
  markers %>%
    group_by(cluster) %>%
    top_n(n = 200, wt = avg_log2FC) -> markers
  
  markers <- markers$gene
  markers <- unique(markers)
  
  markers_list[[i]] <- markers_
  
  features <- append(features, markers)
  
  rm(current)
  gc()
  
}


features <- unique(features)
features <- features[which(features %in% row.names(split_output[[1]]))]
features <- features[which(features %in% row.names(split_output[[2]]))]
features <- features[which(features %in% row.names(split_output[[3]]))]

write.csv(features, "/Users/rosslaidlaw/R/TrAGEDy_V2/Tbrucei_Zc3h20_KO/WT_ZC3H20_KO_feature_space.csv")

split_output$WT_01@misc$feature_space <- features
split_output$WT_02@misc$feature_space <- features
split_output$ZC3H20_KO@misc$feature_space <- features

split_output$WT_01 <- AddMetaData(split_output$WT_01, Idents(split_output$WT_01), "cell_type")
split_output$WT_02 <- AddMetaData(split_output$WT_02, Idents(split_output$WT_02), "cell_type")
split_output$ZC3H20_KO <- AddMetaData(split_output$ZC3H20_KO, Idents(split_output$ZC3H20_KO), "cell_type")

saveRDS(split_output$WT_01,"/Users/rosslaidlaw/R/TrAGEDy_V2/Tbrucei_Zc3h20_KO/objects/WT_01_raw_brucei.rds")
saveRDS(split_output$WT_02,"/Users/rosslaidlaw/R/TrAGEDy_V2/Tbrucei_Zc3h20_KO/objects/WT_02_raw_brucei.rds")
saveRDS(split_output$ZC3H20_KO, "/Users/rosslaidlaw/R/TrAGEDy_V2/Tbrucei_Zc3h20_KO/objects/KO_raw_brucei.rds")

