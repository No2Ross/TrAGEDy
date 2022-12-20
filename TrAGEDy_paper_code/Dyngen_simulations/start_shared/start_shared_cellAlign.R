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
library(clustree)
library(STACAS)
library(monocle)
library(slingshot)
library(stats)
library(stringi)
library(reshape2)
library(stringr)
library(mclust)
library(cellAlign)
library(coin)
library(biomaRt)
library(metap)
setwd("~/R")
source("Scripts/methods/own_method_functions.R")

WT_obj <- readRDS("TrAGEDy_results/simulated_data/start_shared/start_shared_wt.rds")
KO_obj <- readRDS("TrAGEDy_results/simulated_data/start_shared/start_shared_ko.rds")

output <- list(WT_obj, KO_obj)

for (i in 1:length(output)){
  current <- output[[i]]
  
  current <- NormalizeData(current)
  
  current@assays$RNA@var.features <- row.names(current)
  
  all.genes <- rownames(current)
  current <- ScaleData(current, features = all.genes)
  
  current <- RunPCA(current, features = VariableFeatures(object = current))
  
  print(ElbowPlot(current, ndims = 50))
  
  current <- FindNeighbors(current, dims = 1:15)
  current <- FindClusters(current, resolution = 0.5)
  
  current <- RunUMAP(current, dims = 1:15, n.components = 2)
  
  print(DimPlot(current, reduction = "pca"))
  
  
  output[[i]] <- current
  
}

WT_obj <- output[[1]]
KO_obj <- output[[2]]

WT_sce <- as.SingleCellExperiment(WT_obj, assay = "RNA")
start_cluster_1 <- levels(WT_obj$seurat_clusters)[which(table(WT_obj$seurat_clusters[which(WT_obj$sim_time == 0)]) == max(table(WT_obj$seurat_clusters[which(WT_obj$sim_time == 0)])))]
WT_sce <- slingshot(WT_sce, reducedDim = 'UMAP', clusterLabels = WT_sce@colData@listData[["RNA_snn_res.0.5"]], start.clus = start_cluster_1)
WT_sling <- SlingshotDataSet(WT_sce)

KO_sce <- as.SingleCellExperiment(KO_obj, assay = "RNA")

start_cluster_2 <- levels(KO_obj$seurat_clusters)[which(table(KO_obj$seurat_clusters[which(KO_obj$sim_time == 0)]) == max(table(KO_obj$seurat_clusters[which(KO_obj$sim_time == 0)])))]

KO_sce <- slingshot(KO_sce, reducedDim = 'UMAP', clusterLabels = KO_sce@colData@listData[["RNA_snn_res.0.5"]], start.clus = start_cluster_2)
KO_sling <- SlingshotDataSet(KO_sce)

WT.markers <- FindAllMarkers(WT_obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
WT.markers_top15 <- WT.markers %>%
  group_by(cluster) %>%
  slice_max(n = 100, order_by = avg_log2FC)


KO.markers <- FindAllMarkers(KO_obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
KO.markers_top15 <- KO.markers %>%
  group_by(cluster) %>%
  slice_max(n = 100, order_by = avg_log2FC)

features <- unique(c(KO.markers_top15$gene, WT.markers_top15$gene))

pseudo_end <- min(c(max(KO_sce$slingPseudotime_1, WT_sce$slingPseudotime_1, na.rm = T)))
window <- pseudo_end / 45

WT_data <- as.matrix(WT_sce@assays@data$counts)
KO_data <- as.matrix(KO_sce@assays@data$counts)

WT_pseudo <- as.numeric(WT_sce$slingPseudotime_1)
KO_pseudo <- as.numeric(KO_sce$slingPseudotime_1)

numPts = 50
interGlobalWT = cellAlign::interWeights(expDataBatch = WT_data, trajCond = WT_pseudo,
                                        winSz = window, numPts = numPts)
interGlobalKO = cellAlign::interWeights(expDataBatch = KO_data, trajCond = KO_pseudo,
                                        winSz = window, numPts = numPts)

interScaledGlobalWT = cellAlign::scaleInterpolate(interGlobalWT)
interScaledGlobalKO = cellAlign::scaleInterpolate(interGlobalKO)


interScaledGlobalWT$scaledData <- interScaledGlobalWT$scaledData[ which( row.names(interScaledGlobalWT$scaledData) %in% features), ]
interScaledGlobalKO$scaledData <- interScaledGlobalKO$scaledData[ which( row.names(interScaledGlobalKO$scaledData) %in% features), ]


interScaledGlobalWT$scaledError <- interScaledGlobalWT$scaledError[ which( row.names(interScaledGlobalWT$scaledError) %in% features), ]
interScaledGlobalKO$scaledError <- interScaledGlobalKO$scaledError[ which( row.names(interScaledGlobalKO$scaledError) %in% features), ]


alignment = globalAlign(interScaledGlobalWT$scaledData, interScaledGlobalKO$scaledData,
                        scores = list(query = interScaledGlobalWT$traj, 
                                      ref = interScaledGlobalKO$traj),
                        sigCalc = F, numPerm = 20, dist.method = "Euclidean", normDist = F)
plotAlign(alignment)

alignment = globalAlign(interGlobalWT$interpolatedVals[ which( row.names(interGlobalWT$interpolatedVals) %in% features), ], interGlobalKO$interpolatedVals[ which( row.names(interGlobalKO$interpolatedVals) %in% features), ],
                        scores = list(query = interGlobalWT$traj, 
                                      ref = interGlobalKO$traj),
                        sigCalc = F, numPerm = 20, dist.method = "correlation")
plotAlign(alignment)

