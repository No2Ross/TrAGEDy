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
library(monocle)
library(slingshot)
library(stats)
library(stringi)
library(reshape2)
library(stringr)
library(mclust)
setwd("~/R")
source("Scripts/methods/own_method_functions.R")

WT_01_obj <- readRDS("emma_KO_WT_brucei_data/Emma Annotation/WT_01.rds")
WT_02_obj <- readRDS("emma_KO_WT_brucei_data/Emma Annotation/WT_02.rds")

features <- WT_01_obj@misc$feature_space

DimPlot(WT_01_obj, reduction = "umap")

phate_output <- as.matrix(phate(t(WT_01_obj@assays$RNA@data[features,]), ndim=10, seed=1))

colnames(phate_output)<- paste0("PHATE_", 1:ncol(phate_output))
phate.reduction <- CreateDimReducObject(embeddings = phate_output,
                                        key="PHATE_",
                                        assay="integrated")

WT_01_obj@reductions$phate <- phate.reduction

WT_01_sce <- as.SingleCellExperiment(WT_01_obj, assay = "RNA")
WT_01_sce <- slingshot(WT_01_sce, reducedDim = 'PHATE', clusterLabels = WT_01_sce@colData@listData[["cell_type"]], start.clus = "LS A")
WT_01_sling <- SlingshotDataSet(WT_01_sce)



phate_output <- as.matrix(phate(t(WT_02_obj@assays$RNA@data[features,]), ndim=10, seed=1))
phate.reduction <- CreateDimReducObject(embeddings = phate_output,
                                        key="PHATE_",
                                        assay="integrated")

WT_02_obj@reductions$phate <- phate.reduction
WT_02_sce <- as.SingleCellExperiment(WT_02_obj, assay = "RNA")
WT_02_sce <- slingshot(WT_02_sce, reducedDim = 'PHATE', clusterLabels = WT_02_sce@colData@listData[["cell_type"]], start.clus = "LS A")
WT_02_sling <- SlingshotDataSet(WT_02_sce)


pseudo_end <- max(WT_02_sce$slingPseudotime_1, WT_01_sce$slingPseudotime_1)
window <- pseudo_end / 45

WT_01_cell_pseudotime <- matrix(WT_01_sce$slingPseudotime_1, dimnames =list(WT_01_sce@colData@rownames))
WT_02_cell_pseudotime <- matrix(WT_02_sce$slingPseudotime_1, dimnames =list(WT_02_sce@colData@rownames))
WT_01_ID <- data.frame(WT_01_sce$cell_type, row.names =WT_01_sce@colData@rownames)
WT_02_ID <- data.frame(WT_02_sce$cell_type, row.names =WT_02_sce@colData@rownames)

source("Scripts/methods/own_method_functions.R")
WT_01_tree <- nodePseudotime(WT_01_cell_pseudotime,WT_01_ID, 50, "WT_01")
WT_02_tree <- nodePseudotime(WT_02_cell_pseudotime,WT_02_ID, 50, "WT_02")

WT_02_node_pseudotime <- matrix(WT_02_tree$pseudotime$pseudotime , dimnames = list(row.names(WT_02_tree$pseudotime)), )
WT_01_node_pseudotime <- matrix(WT_01_tree$pseudotime$pseudotime , dimnames = list(row.names(WT_01_tree$pseudotime)), )

source("Scripts/methods/own_method_functions.R")
WT_02_node_exp_mtx <- nodeExpressionEstimate(WT_02_sce@assays@data@listData$logcounts, WT_02_node_pseudotime, WT_02_cell_pseudotime, window, adjust.window = T)

WT_01_node_exp_mtx <- nodeExpressionEstimate(WT_01_sce@assays@data@listData$logcounts, WT_01_node_pseudotime, WT_01_cell_pseudotime, window, adjust.window = T)


WT_02_node_exp_mtx <- WT_02_node_exp_mtx[[1]]
WT_01_node_exp_mtx <- WT_01_node_exp_mtx[[1]]

WT_02_node_exp_mtx  <- WT_02_node_exp_mtx[ features, ]
WT_01_node_exp_mtx <- WT_01_node_exp_mtx[ features, ]

penalty_mtx_cut <- dis_mtx_calculator(as.matrix(WT_01_node_exp_mtx), as.matrix(WT_02_node_exp_mtx), "spearman")

output_solution <- pathfind(penalty_mtx_cut, cut_type = "minimum", method = "mean")
output_solution_cut <- cut_deviate(output_solution, penalty_mtx_cut, method = "mean")

PlotAlignment(output_solution_cut, penalty_mtx_cut)

pseudo_1 <- WT_02_tree$pseudotime
pseudo_2 <- WT_01_tree$pseudotime
alignment <- output_solution_cut
id1 <- WT_02_tree$ID
id2 <- WT_01_tree$ID

PlotOutput(WT_01_tree, WT_02_tree, alignment)

test <- chunk_node(WT_01_tree$pseudotime, WT_02_tree$pseudotime, output_solution_cut)

WT_01_tree_new <- WT_01_tree
WT_01_tree_new$pseudotime <- data.frame(test$condition_1$pseudotime, row.names = row.names(WT_01_tree$pseudotime) )

WT_02_tree_new <- WT_02_tree
WT_02_tree_new$pseudotime <- data.frame(test$condition_2$pseudotime, row.names = row.names(WT_02_tree$pseudotime) )


WT_02_cell_pseudo <- data.frame(WT_02_sce$slingPseudotime_1, row.names = WT_02_sce@colData@rownames)
WT_02_node_pseudo_new <- data.frame(row.names(test$condition_2),test$condition_2$pseudotime, test$condition_2$alignment, row.names =row.names(test$condition_2) )
WT_02_node_pseudo <- data.frame( pseudotime = WT_02_tree$pseudotime, row.names = row.names(WT_02_tree$pseudotime))

WT_02_cell_pseudo_new <- pseudo_cell_align(WT_02_cell_pseudo , WT_02_node_pseudo_new, WT_02_node_pseudo, window)
WT_02_cell_pseudo_new <- WT_02_cell_pseudo_new[order(match(row.names(WT_02_cell_pseudo_new), row.names(WT_02_cell_pseudo))),]

WT_01_cell_pseudo <- data.frame(WT_01_sce$slingPseudotime_1, row.names = WT_01_sce@colData@rownames)
WT_01_node_pseudo_new <- data.frame(row.names(test$condition_1),test$condition_1$pseudotime, test$condition_1$align, row.names =row.names(test$condition_1) )
WT_01_node_pseudo <- data.frame(pseudotime = WT_01_tree$pseudotime, row.names = row.names(WT_01_tree$pseudotime))

WT_01_cell_pseudo_new <- pseudo_cell_align(WT_01_cell_pseudo , WT_01_node_pseudo_new, WT_01_node_pseudo, window)
WT_01_cell_pseudo_new <- WT_01_cell_pseudo_new[order(match(row.names(WT_01_cell_pseudo_new), row.names(WT_01_cell_pseudo))),]

#Add metadata
WT_02_sce$oldPseudotime <- WT_02_sce$slingPseudotime_1
WT_02_sce$newPseudotime <- WT_02_cell_pseudo_new$pseudotime
WT_01_sce$oldPseudotime <- WT_01_sce$slingPseudotime_1
WT_01_sce$newPseudotime <- WT_01_cell_pseudo_new$pseudotime
WT_02_sce$Status <- as.factor(WT_02_cell_pseudo_new$status)
WT_01_sce$Status <- as.factor(WT_01_cell_pseudo_new$status)

WT_02_obj$oldPseudotime <- WT_02_sce$slingPseudotime_1
WT_02_obj$newPseudotime <- WT_02_cell_pseudo_new$pseudotime
WT_01_obj$oldPseudotime <- WT_01_sce$slingPseudotime_1
WT_01_obj$newPseudotime <- WT_01_cell_pseudo_new$pseudotime
WT_02_obj$Status <- as.factor(WT_02_cell_pseudo_new$status)
WT_01_obj$Status <- as.factor(WT_01_cell_pseudo_new$status)




