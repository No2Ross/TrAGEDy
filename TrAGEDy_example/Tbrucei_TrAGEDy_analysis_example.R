library(slingshot)
library(RColorBrewer)
library(Seurat)
library(SingleCellExperiment)
library(ggplot2)
library(dplyr)
library(stats)
library(stringr)
library(rgl)

setwd("path/to/tutorial")
source("TrAGEDy_functions.R")

set.seed(42)

WT_1_sce <- readRDS("WT_1_sce.rds")
WT_2_sce <- readRDS("WT_2_sce.rds")

features <- read.csv("WTAlignment_brucei_feature_space.csv", row.names = 1)
features <- features$x

pseudo_end <- max(WT_2_sce$slingPseudotime_1, WT_1_sce$slingPseudotime_1)
window <- pseudo_end / 45

WT_1_tree <- nodePseudotime(WT_1_sce,"slingPseudotime_1","cell_type", 50, "WT1")
WT_2_tree <- nodePseudotime(WT_2_sce,"slingPseudotime_1","cell_type", 50, "WT2")

WT_1_node_exp_mtx <- nodeExpressionEstimate(WT_1_sce@assays@data@listData$logcounts, WT_1_tree, window, adjust.window = T)
WT_2_node_exp_mtx <- nodeExpressionEstimate(WT_2_sce@assays@data@listData$logcounts, WT_2_tree, window, adjust.window = T)

WT_2_node_exp_mtx  <- WT_2_node_exp_mtx[ features, ]
WT_1_node_exp_mtx <- WT_1_node_exp_mtx[ features, ]

penalty_mtx_cut <- dis_mtx_calculator(as.matrix(WT_1_node_exp_mtx), as.matrix(WT_2_node_exp_mtx), "spearman")

path_tragedy <- bootstrap_pathfind(sequence_1 = as.matrix(WT_1_node_exp_mtx), sequence_2 = as.matrix(WT_2_node_exp_mtx)
                        , similarity_method= "spearman", threshold_method = "mean")

path_tragedy_cut <- cut_deviate(path_tragedy[[1]], penalty_mtx_cut, method = "mean")

PlotAlignment(path_tragedy_cut, penalty_mtx_cut)

PlotOutput(WT_1_tree, WT_2_tree, path_tragedy_cut)


alignedPoints <- chunk_node(WT_1_tree$node_pseudotime, WT_2_tree$node_pseudotime, path_tragedy_cut)


WT_1_tree_aligned <- WT_1_tree
WT_1_tree_aligned$node_pseudotime <- alignedPoints$condition_1$pseudotime
names(WT_1_tree_aligned$node_pseudotime) <- row.names(alignedPoints$condition_1)

WT_2_tree_aligned <- WT_2_tree
WT_2_tree_aligned$node_pseudotime <- alignedPoints$condition_2$pseudotime
names(WT_2_tree_aligned$node_pseudotime) <- row.names(alignedPoints$condition_2)

PlotOutput(WT_1_tree_aligned, WT_2_tree_aligned, path_tragedy_cut)

WT_1_cell_pseudo_new <- pseudo_cell_align(WT_1_tree$cell_pseudotime, alignedPoints$condition_1 , WT_1_tree$node_pseudotime, window)
WT_2_cell_pseudo_new <- pseudo_cell_align(WT_2_tree$cell_pseudotime, alignedPoints$condition_2 , WT_2_tree$node_pseudotime, window)

WT_2_sce$oldPseudotime <- WT_2_sce$slingPseudotime_1
WT_2_sce$newPseudotime <- WT_2_cell_pseudo_new$pseudotime
WT_1_sce$oldPseudotime <- WT_1_sce$slingPseudotime_1
WT_1_sce$newPseudotime <- WT_1_cell_pseudo_new$pseudotime
WT_2_sce$Status <- as.factor(WT_2_cell_pseudo_new$status)
WT_1_sce$Status <- as.factor(WT_1_cell_pseudo_new$status)

output <- TrajDE(list(WT_1_sce, WT_2_sce), list(WT_1_tree_aligned, WT_2_tree_aligned), path_tragedy_cut, n_windows = 4, 
                 overlap = 1, p_val = 0.05, min.pct = 0.1, logfc = 0.5, all.genes = F, test_use = "wilcox", correct = T)



