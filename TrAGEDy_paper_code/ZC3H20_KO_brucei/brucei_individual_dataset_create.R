#Trade seq and STACAS do not get along 
library(reticulate)
use_python("/usr/local/bin/python3.8")
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
library(STACAS)
library(tradeSeq)
library(monocle)
library(slingshot)
library(stats)
library(stringi)
library(reshape2)
library(stringr)
library(mclust)
setwd("~/R")
source("Scripts/methods/own_method_functions.R")

emma_obj <- readRDS("emma_KO_WT_brucei_data/WT_KO_int_phate/WT_ZC3H20_with_phate.rds")

emma_obj <- AddMetaData(emma_obj, emma_obj@active.ident, col.name = "cell_type")

split_output <- SplitObject(emma_obj, "orig.ident")

#Merge the A and B subsets into one
split_output$WT_01 <- RenameIdents(split_output$WT_01, "LS A.1" = "LS A", "LS A.2"= "LS A", "LS B.1"="LS B", "LS B.2" ="LS B")
split_output$WT_02 <- RenameIdents(split_output$WT_02, "LS A.1" = "LS A", "LS A.2"= "LS A", "LS B.1"="LS B", "LS B.2" ="LS B")
split_output$ZC3H20_KO <- RenameIdents(split_output$ZC3H20_KO, "LS A.1" = "LS A", "LS A.2"= "LS A", "LS B.1"="LS B", "LS B.2" ="LS B")


for (i in 1:length(split_output)){
  current <- split_output[[i]]
  DefaultAssay(current) <- "RNA"
  
  current <- FindVariableFeatures(current, selection.method = "vst", nfeatures = 2000)
  all.genes <- rownames(current)
  current <- ScaleData(current, features = all.genes)
  
  current <- RunPCA(current, features = VariableFeatures(object = current))
  
  current <- RunUMAP(current, dims = 1:10, n.components = 3)
  
  split_output[[i]] <- current
  
}

#Get feature space for alignment
features <- c()
for (i in 1:length(split_output)){
  
  current <- split_output[[i]]
  DefaultAssay(current) <- "RNA"
  
  markers <- FindAllMarkers(current)
  
  markers %>%
    group_by(cluster) %>%
    top_n(n = 300, wt = avg_log2FC) -> markers
  
  markers <- markers$gene
  markers <- unique(markers)
  
  split_output[[i]] <- current
  features <- append(features, markers)
  
  rm(current)
  gc()
  
}

features <- unique(features)
features <- features[which(features %in% row.names(split_output[[1]]))]
features <- features[which(features %in% row.names(split_output[[2]]))]
features <- features[which(features %in% row.names(split_output[[3]]))]

split_output$WT_01@misc$feature_space <- features
split_output$WT_02@misc$feature_space <- features
split_output$ZC3H20_KO@misc$feature_space <- features

split_output$WT_01 <- AddMetaData(split_output$WT_01, Idents(split_output$WT_01), "cell_type")
split_output$WT_02 <- AddMetaData(split_output$WT_02, Idents(split_output$WT_02), "cell_type")
split_output$ZC3H20_KO <- AddMetaData(split_output$ZC3H20_KO, Idents(split_output$ZC3H20_KO), "cell_type")

