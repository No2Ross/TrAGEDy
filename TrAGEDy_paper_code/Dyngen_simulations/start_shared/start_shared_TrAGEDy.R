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
source("Scripts/methods/own_method_functions.R")

#childhood/quite young (2/3 years old) - Getting out of the cradle (dressed in wellies), went to the front door and opened it and walked to his nanas. 
#Remembers seeing this head bobbed up and down. The street was slopped so the bob was the walk of the person. Nana says it was true

#Midteens - Had a girlfriend at school. He never said anything to his mum and dad. She came round with a christmas present and he brushed it off and said thanks and shut the door

#Young adult - He had flown a couple of times before the cayman islands, never flown transatlantic. Never gone to a really hot country. Plane lansded in the evening in the cayman islands
# He walked out the door and was hit by a wall of heat he hadn't experienced before

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


#WT_sce$slingPseudotime_1 <- WT_sce$sim_time
#KO_sce$slingPseudotime_1 <- KO_sce$sim_time

pseudo_end <- min(c(max(KO_sce$slingPseudotime_1, WT_sce$slingPseudotime_1, na.rm = T)))
window <- pseudo_end / 45

#3: create metacells across pseudotime 
WT_cell_pseudotime <- matrix(WT_sce$slingPseudotime_1, dimnames =list(WT_sce@colData@rownames))
KO_cell_pseudotime <- matrix(KO_sce$slingPseudotime_1, dimnames =list(KO_sce@colData@rownames))
WT_ID <- data.frame(WT_sce$cell_type, row.names =WT_sce@colData@rownames)
KO_ID <- data.frame(KO_sce$cell_type, row.names =KO_sce@colData@rownames)

source("Scripts/methods/own_method_functions.R")
WT_tree <- nodePseudotime(WT_cell_pseudotime,WT_ID, 50, "WT")
KO_tree <- nodePseudotime(KO_cell_pseudotime,KO_ID, 50, "KO")

#cellalign node exp mtx - not scaled with cellalign way
KO_cell_pseudo <- data.frame("ID" = KO_sce@colData@rownames, "pseudo" = KO_sce$slingPseudotime_1)
KO_node_pseudo <- data.frame("ID" = row.names(KO_tree$pseudotime), "pseudo" = KO_tree$pseudotime$pseudotime)

WT_cell_pseudo <- data.frame("ID" = WT_sce@colData@rownames, "pseudo" = WT_sce$slingPseudotime_1)
WT_node_pseudo <- data.frame("ID" = row.names(WT_tree$pseudotime), "pseudo" = WT_tree$pseudotime$pseudotime)

KO_node_pseudotime <- matrix(KO_tree$pseudotime$pseudotime , dimnames = list(row.names(KO_tree$pseudotime)), )
WT_node_pseudotime <- matrix(WT_tree$pseudotime$pseudotime , dimnames = list(row.names(WT_tree$pseudotime)), )

source("Scripts/methods/own_method_functions.R")
KO_node_exp_mtx <- nodeExpressionEstimate(KO_sce@assays@data@listData$logcounts, KO_node_pseudotime, KO_cell_pseudotime, window, adjust.window = T)

WT_node_exp_mtx <- nodeExpressionEstimate(WT_sce@assays@data@listData$logcounts, WT_node_pseudotime, WT_cell_pseudotime, window, adjust.window = T)

KO_node_exp_mtx  <- KO_node_exp_mtx[ features, ]
WT_node_exp_mtx <- WT_node_exp_mtx[ features, ]

#3.5: remove genes from our feature space which are not temporarily regulated in our process

#4.1: Create correlation matrix based on the new feature space we used across the metacells
row.names(as.matrix(KO_node_exp_mtx)) == row.names(as.matrix(WT_node_exp_mtx))
source("Scripts/methods/own_method_functions.R")

penalty_mtx_cut <- dis_mtx_calculator(as.matrix(WT_node_exp_mtx), as.matrix(KO_node_exp_mtx), "spearman")

source("Scripts/methods/own_method_functions.R")
output_solution <- pathfind(penalty_mtx_cut, cut_type = "minimum", method = "mean")
output_solution_cut <- cut_deviate(output_solution, penalty_mtx_cut, method = "mean")

#output_solution_cut$Status <- rep("cut", length(output_solution_cut$Status))
PlotAlignment(output_solution_cut, penalty_mtx_cut)



