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
library(coin)
library(biomaRt)
library(metap)
setwd("~/R")
source("TrAGEDy_results/TrAGEDy_functions.R")

c_obj <- readRDS("TrAGEDy_results/simulated_data/not_shared/negative_control_1.rds")
d_obj <- readRDS("TrAGEDy_results/simulated_data/not_shared/negative_control_2.rds")

output <- list(c_obj, d_obj)

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
  
  
  print(DimPlot(current, reduction = "umap"))
  
  
  output[[i]] <- current
  
}

c_obj <- output[[1]]
d_obj <- output[[2]]

#remove outlier in the d object
d_obj <- subset(d_obj, seurat_clusters != "7")

c_sce <- as.SingleCellExperiment(c_obj, assay = "RNA")
start_cluster_1 <- levels(c_obj$seurat_clusters)[which(table(c_obj$seurat_clusters[which(c_obj$sim_time == 0)]) == max(table(c_obj$seurat_clusters[which(c_obj$sim_time == 0)])))]
c_sce <- slingshot(c_sce, reducedDim = 'UMAP', clusterLabels = c_sce@colData@listData[["RNA_snn_res.0.5"]], start.clus = start_cluster_1)


d_sce <- as.SingleCellExperiment(d_obj, assay = "RNA")
start_cluster_2 <- levels(d_obj$seurat_clusters)[which(table(d_obj$seurat_clusters[which(d_obj$sim_time == 0)]) == max(table(d_obj$seurat_clusters[which(d_obj$sim_time == 0)])))]
d_sce <- slingshot(d_sce, reducedDim = 'UMAP', clusterLabels = d_sce@colData@listData[["RNA_snn_res.0.5"]], start.clus = start_cluster_2)
d_sling <- SlingshotDataSet(d_sce)

WT.markers <- FindAllMarkers(c_obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
WT.markers_top15 <- WT.markers %>%
  group_by(cluster) %>%
  slice_max(n = 100, order_by = avg_log2FC)


KO.markers <- FindAllMarkers(d_obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
KO.markers_top15 <- KO.markers %>%
  group_by(cluster) %>%
  slice_max(n = 100, order_by = avg_log2FC)

features <- unique(c(KO.markers_top15$gene, WT.markers_top15$gene))
# features <- row.names(d_obj)

#c_sce$slingPseudotime_1 <- c_sce$sim_time
#d_sce$slingPseudotime_1 <- d_sce$sim_time

pseudo_end <- min(c(max(d_sce$slingPseudotime_1, c_sce$slingPseudotime_1, na.rm = T)))
window <- pseudo_end / 45

#3: create metacells across pseudotime 
c_cell_pseudotime <- matrix(c_sce$slingPseudotime_1, dimnames =list(c_sce@colData@rownames))
d_cell_pseudotime <- matrix(d_sce$slingPseudotime_1, dimnames =list(d_sce@colData@rownames))
c_ID <- data.frame(c_sce$seurat_clusters, row.names =c_sce@colData@rownames)
d_ID <- data.frame(d_sce$seurat_clusters, row.names =d_sce@colData@rownames)

c_tree <- nodePseudotime(c_cell_pseudotime,c_ID, 50, "WT")
d_tree <- nodePseudotime(d_cell_pseudotime,d_ID, 50, "KO")

#cellalign node exp mtx - not scaled with cellalign way
d_cell_pseudo <- data.frame("ID" = d_sce@colData@rownames, "pseudo" = d_sce$slingPseudotime_1)
d_node_pseudo <- data.frame("ID" = row.names(d_tree$pseudotime), "pseudo" = d_tree$pseudotime$pseudotime)

c_cell_pseudo <- data.frame("ID" = c_sce@colData@rownames, "pseudo" = c_sce$slingPseudotime_1)
c_node_pseudo <- data.frame("ID" = row.names(c_tree$pseudotime), "pseudo" = c_tree$pseudotime$pseudotime)

d_node_pseudotime <- matrix(d_tree$pseudotime$pseudotime , dimnames = list(row.names(d_tree$pseudotime)), )
c_node_pseudotime <- matrix(c_tree$pseudotime$pseudotime , dimnames = list(row.names(c_tree$pseudotime)), )

d_node_exp_mtx <- nodeExpressionEstimate(d_sce@assays@data@listData$logcounts, d_node_pseudotime, d_cell_pseudotime, window, adjust.window = T)

c_node_exp_mtx <- nodeExpressionEstimate(c_sce@assays@data@listData$logcounts, c_node_pseudotime, c_cell_pseudotime, window, adjust.window = T)

d_node_exp_mtx  <- d_node_exp_mtx[ features, ]
c_node_exp_mtx <- c_node_exp_mtx[ features, ]

#3.5: remove genes from our feature space which are not temporarily regulated in our process

#4.1: Create correlation matrix based on the new feature space we used across the metacells
row.names(as.matrix(d_node_exp_mtx)) == row.names(as.matrix(c_node_exp_mtx))

penalty_mtx_cut <- dis_mtx_calculator(as.matrix(c_node_exp_mtx), as.matrix(d_node_exp_mtx), "spearman")

output_solution <- pathfind(penalty_mtx_cut, cut_type = "minimum", method = "mean")
output_solution_cut <- cut_deviate(output_solution, penalty_mtx_cut, method = "mean")

#output_solution_cut$Status <- rep("cut", length(output_solution_cut$Status))
PlotAlignment(output_solution_cut, penalty_mtx_cut)


