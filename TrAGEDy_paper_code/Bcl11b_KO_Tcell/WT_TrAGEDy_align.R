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

genesV2 <- read.csv("path/to/seurat/cellcycle/genes/human_to_mouse_seurat_cellcycle.csv")

cell_cycle <- readxl::read_xlsx("path/to/cellcycle_genes/mgi_cellCycle_GO.xlsx")

WT_1_obj <- readRDS("path/to/annotated/WT/1/object/WT_1_obj.rds")
WT_2_obj <- readRDS("path/to/annotated/WT/2/object/WT_2_obj.rds")
features <- read.delim("path/to/feature/space/feature_space.csv", sep = ",")[,2]

features <- setdiff(features, unique(c(genesV2$MGI.symbol, cell_cycle$Symbol)))

phate_output <- as.matrix(phate(t(WT_1_obj@assays$RNA@data[features,]), ndim=10, seed=1))

phate.reduction <- CreateDimReducObject(embeddings = phate_output,
                                        key="PHATE_",
                                        assay="RNA")

WT_1_obj@reductions$phate <- phate.reduction


start_clus1 <- "Cd34_TSP"

WT_1_sce <- as.SingleCellExperiment(WT_1_obj, assay = "RNA")
WT_1_sce <- slingshot(WT_1_sce, reducedDim = 'PHATE', clusterLabels = WT_1_sce@colData@listData[["cell_type"]], start.clus = start_clus1, use.median = T, dist.method = "simple")
WT_1_sling <- SlingshotDataSet(WT_1_sce)

phate_output <- as.matrix(phate(t(WT_2_obj@assays$RNA@data[features,]), ndim=10, seed=1))

phate.reduction <- CreateDimReducObject(embeddings = phate_output,
                                        key="PHATE_",
                                        assay="RNA")

WT_2_obj@reductions$phate <- phate.reduction

start_clus2 <- "early_T"

WT_2_sce <- as.SingleCellExperiment(WT_2_obj, assay = "RNA")
WT_2_sce <- slingshot(WT_2_sce, reducedDim = 'PHATE', clusterLabels = WT_2_sce@colData@listData[["cell_type"]], start.clus = start_clus2, use.median = T, dist.method = "simple")
WT_2_sling <- SlingshotDataSet(WT_2_sce)


WT_1_sce <- subset(WT_1_sce, , slingPseudotime_1 < 0.1)

pseudo_end <- min(c(max(WT_1_sce$slingPseudotime_1, WT_2_sce$slingPseudotime_1)))
window <- pseudo_end / 90

WT_1_cell_pseudotime <- matrix(WT_1_sce$slingPseudotime_1, dimnames =list(WT_1_sce@colData@rownames))
WT_2_cell_pseudotime <- matrix(WT_2_sce$slingPseudotime_1, dimnames =list(WT_2_sce@colData@rownames))
WT_1_ID <- data.frame(WT_1_sce$cell_type, row.names =WT_1_sce@colData@rownames)
WT_2_ID <- data.frame(WT_2_sce$cell_type, row.names =WT_2_sce@colData@rownames)

WT_1_tree <- nodePseudotime(WT_1_cell_pseudotime,WT_1_ID, 100, "WT1")
WT_2_tree <- nodePseudotime(WT_2_cell_pseudotime,WT_2_ID, 100, "WT2")

WT_2_cell_pseudo <- data.frame("ID" = WT_2_sce@colData@rownames, "pseudo" = WT_2_sce$slingPseudotime_1)
WT_2_node_pseudo <- data.frame("ID" = row.names(WT_2_tree$pseudotime), "pseudo" = WT_2_tree$pseudotime$pseudotime)

WT_1_cell_pseudo <- data.frame("ID" = WT_1_sce@colData@rownames, "pseudo" = WT_1_sce$slingPseudotime_1)
WT_1_node_pseudo <- data.frame("ID" = row.names(WT_1_tree$pseudotime), "pseudo" = WT_1_tree$pseudotime$pseudotime)

WT_2_node_pseudotime <- matrix(WT_2_tree$pseudotime$pseudotime , dimnames = list(row.names(WT_2_tree$pseudotime)), )
WT_1_node_pseudotime <- matrix(WT_1_tree$pseudotime$pseudotime , dimnames = list(row.names(WT_1_tree$pseudotime)), )

WT_2_node_exp_mtx <- nodeExpressionEstimate(WT_2_sce@assays@data@listData$logcounts, WT_2_node_pseudotime, WT_2_cell_pseudotime, window, adjust.window = T)

WT_1_node_exp_mtx <- nodeExpressionEstimate(WT_1_sce@assays@data@listData$logcounts, WT_1_node_pseudotime, WT_1_cell_pseudotime, window, adjust.window = T)

WT_2_node_exp_mtx <- WT_2_node_exp_mtx[[1]]
WT_1_node_exp_mtx <- WT_1_node_exp_mtx[[1]]

WT_2_node_exp_mtx  <- WT_2_node_exp_mtx[ features, ]
WT_1_node_exp_mtx <- WT_1_node_exp_mtx[ features, ]

#Calculate dissimilarity matrix
penalty_mtx_cut <- dis_mtx_calculator(as.matrix(WT_1_node_exp_mtx), as.matrix(WT_2_node_exp_mtx), "spearman")

#Find optimal path
output_solution <- pathfind(penalty_mtx_cut, cut_type = "minimum", method = "mean")
output_solution_cut <- cut_deviate(output_solution, penalty_mtx_cut, method = "mean")

PlotAlignment(output_solution_cut, penalty_mtx_cut)

PlotOutput(WT_1_tree, WT_2_tree, output_solution_cut)

#Align pseudotime of interpolated points
test <- chunk_node(WT_1_tree$pseudotime, WT_2_tree$pseudotime, output_solution_cut)

WT_1_tree_new <- WT_1_tree
WT_1_tree_new$pseudotime <- data.frame(test$condition_1$pseudotime, row.names = row.names(WT_1_tree$pseudotime) )

WT_2_tree_new <- WT_2_tree
WT_2_tree_new$pseudotime <- data.frame(test$condition_2$pseudotime, row.names = row.names(WT_2_tree$pseudotime) )

PlotOutput(WT_1_tree_new, WT_2_tree_new, output_solution_cut)

WT_2_cell_pseudo <- data.frame(WT_2_sce$slingPseudotime_1, row.names = WT_2_sce@colData@rownames)
WT_2_node_pseudo_new <- data.frame(row.names(test$condition_2),test$condition_2$pseudotime, test$condition_2$alignment, row.names =row.names(test$condition_2) )
WT_2_node_pseudo <- data.frame( pseudotime = WT_2_tree$pseudotime, row.names = row.names(WT_2_tree$pseudotime))

WT_2_cell_pseudo_new <- pseudo_cell_align(WT_2_cell_pseudo , WT_2_node_pseudo_new, WT_2_node_pseudo, window)
WT_2_cell_pseudo_new <- WT_2_cell_pseudo_new[order(match(row.names(WT_2_cell_pseudo_new), row.names(WT_2_cell_pseudo))),]

WT_1_cell_pseudo <- data.frame(WT_1_sce$slingPseudotime_1, row.names = WT_1_sce@colData@rownames)
WT_1_node_pseudo_new <- data.frame(row.names(test$condition_1),test$condition_1$pseudotime, test$condition_1$align, row.names =row.names(test$condition_1) )
WT_1_node_pseudo <- data.frame(pseudotime = WT_1_tree$pseudotime, row.names = row.names(WT_1_tree$pseudotime))

WT_1_cell_pseudo_new <- pseudo_cell_align(WT_1_cell_pseudo , WT_1_node_pseudo_new, WT_1_node_pseudo, window)
WT_1_cell_pseudo_new <- WT_1_cell_pseudo_new[order(match(row.names(WT_1_cell_pseudo_new), row.names(WT_1_cell_pseudo))),]

#Add metadata
WT_2_sce$oldPseudotime <- WT_2_sce$slingPseudotime_1
WT_2_sce$newPseudotime <- WT_2_cell_pseudo_new$pseudotime
WT_1_sce$oldPseudotime <- WT_1_sce$slingPseudotime_1
WT_1_sce$newPseudotime <- WT_1_cell_pseudo_new$pseudotime
WT_2_sce$Status <- as.factor(WT_2_cell_pseudo_new$status)
WT_1_sce$Status <- as.factor(WT_1_cell_pseudo_new$status)
WT_1_sce$cell_names <- colnames(WT_1_sce)

WT_1_obj$cell_names <- colnames(WT_1_obj)
WT_1_obj <- subset(WT_1_obj, cell_names %in% WT_1_sce$cell_names)

WT_2_obj$oldPseudotime <- WT_2_sce$slingPseudotime_1
WT_2_obj$newPseudotime <- WT_2_cell_pseudo_new$pseudotime
WT_1_obj$oldPseudotime <- WT_1_sce$slingPseudotime_1
WT_1_obj$newPseudotime <- WT_1_cell_pseudo_new$pseudotime
WT_2_obj$Status <- as.factor(WT_2_cell_pseudo_new$status)
WT_1_obj$Status <- as.factor(WT_1_cell_pseudo_new$status)
