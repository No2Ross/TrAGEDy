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

KO_1_obj <- readRDS("path/to/annotated/KO/1/object/KO_1_obj.rds")
KO_2_obj <- readRDS("path/to/annotated/KO/2/object/KO_2_obj.rds")
features <- read.delim("path/to/feature/space/feature_space.csv", sep = ",")[,2]

features <- setdiff(features, unique(c(genesV2$MGI.symbol, cell_cycle$Symbol)))

phate_output <- as.matrix(phate(t(KO_1_obj@assays$RNA@data[features,]), ndim=10, seed=1))

phate.reduction <- CreateDimReducObject(embeddings = phate_output,
                                        key="PHATE_",
                                        assay="RNA")

KO_1_obj@reductions$phate <- phate.reduction

DimPlot(KO_1_obj, reduction="phate", label = T, group.by = "cell_type")

start_clus1 <- "Cd34_TSP"

KO_1_sce <- as.SingleCellExperiment(KO_1_obj, assay = "RNA")
KO_1_sce <- slingshot(KO_1_sce, reducedDim = 'PHATE', clusterLabels = KO_1_sce@colData@listData[["cell_type"]], start.clus = start_clus1, use.median = T, dist.method = "simple")
KO_1_sling <- SlingshotDataSet(KO_1_sce)


phate_output <- as.matrix(phate(t(KO_2_obj@assays$RNA@data[features,]), ndim=10, seed=1))

phate.reduction <- CreateDimReducObject(embeddings = phate_output,
                                        key="PHATE_",
                                        assay="RNA")

KO_2_obj@reductions$phate <- phate.reduction

DimPlot(KO_2_obj, reduction="phate", label = T, group.by = "cell_type")


start_clus2 <- "Cd34_TSP"
KO_2_sce <- as.SingleCellExperiment(KO_2_obj, assay = "RNA")
KO_2_sce <- slingshot(KO_2_sce, reducedDim = 'PHATE', clusterLabels = KO_2_sce@colData@listData[["cell_type"]], start.clus = start_clus2, use.median = T, dist.method = "simple")
KO_2_sling <- SlingshotDataSet(KO_2_sce)

pseudo_end <- min(c(max(KO_1_sce$slingPseudotime_1, KO_2_sce$slingPseudotime_1)))
window <- pseudo_end / 90

KO_1_cell_pseudotime <- matrix(KO_1_sce$slingPseudotime_1, dimnames =list(KO_1_sce@colData@rownames))
KO_2_cell_pseudotime <- matrix(KO_2_sce$slingPseudotime_1, dimnames =list(KO_2_sce@colData@rownames))
KO_1_ID <- data.frame(KO_1_sce$cell_type, row.names =KO_1_sce@colData@rownames)
KO_2_ID <- data.frame(KO_2_sce$cell_type, row.names =KO_2_sce@colData@rownames)

KO_1_tree <- nodePseudotime(KO_1_cell_pseudotime,KO_1_ID, 100, "KO1")
KO_2_tree <- nodePseudotime(KO_2_cell_pseudotime,KO_2_ID, 100, "KO2")

KO_2_cell_pseudo <- data.frame("ID" = KO_2_sce@colData@rownames, "pseudo" = KO_2_sce$slingPseudotime_1)
KO_2_node_pseudo <- data.frame("ID" = row.names(KO_2_tree$pseudotime), "pseudo" = KO_2_tree$pseudotime$pseudotime)

KO_1_cell_pseudo <- data.frame("ID" = KO_1_sce@colData@rownames, "pseudo" = KO_1_sce$slingPseudotime_1)
KO_1_node_pseudo <- data.frame("ID" = row.names(KO_1_tree$pseudotime), "pseudo" = KO_1_tree$pseudotime$pseudotime)

KO_2_node_pseudotime <- matrix(KO_2_tree$pseudotime$pseudotime , dimnames = list(row.names(KO_2_tree$pseudotime)), )
KO_1_node_pseudotime <- matrix(KO_1_tree$pseudotime$pseudotime , dimnames = list(row.names(KO_1_tree$pseudotime)), )

KO_2_node_exp_mtx <- nodeExpressionEstimate(KO_2_sce@assays@data@listData$logcounts, KO_2_node_pseudotime, KO_2_cell_pseudotime, window, adjust.window = T)

KO_1_node_exp_mtx <- nodeExpressionEstimate(KO_1_sce@assays@data@listData$logcounts, KO_1_node_pseudotime, KO_1_cell_pseudotime, window, adjust.window = T)

KO_2_node_exp_mtx <- KO_2_node_exp_mtx[[1]]
KO_1_node_exp_mtx <- KO_1_node_exp_mtx[[1]]

#Subset down to feature space
KO_2_node_exp_mtx  <- KO_2_node_exp_mtx[ features, ]
KO_1_node_exp_mtx <- KO_1_node_exp_mtx[ features, ]

#Calculate dissimilarity matrix
penalty_mtx_cut <- dis_mtx_calculator(as.matrix(KO_1_node_exp_mtx), as.matrix(KO_2_node_exp_mtx), "spearman")

#Find optimal path
output_solution <- pathfind(penalty_mtx_cut, cut_type = "minimum", method = "mean")
output_solution_cut <- cut_deviate(output_solution, penalty_mtx_cut, method = "mean")

PlotAlignment(output_solution_cut, penalty_mtx_cut)

PlotOutput(KO_1_tree, KO_2_tree, output_solution_cut)

#align the pseudotime of interpolated points
test <- chunk_node(KO_1_tree$pseudotime, KO_2_tree$pseudotime, output_solution_cut)

KO_1_tree_new <- KO_1_tree
KO_1_tree_new$pseudotime <- data.frame(test$condition_1$pseudotime, row.names = row.names(KO_1_tree$pseudotime) )

KO_2_tree_new <- KO_2_tree
KO_2_tree_new$pseudotime <- data.frame(test$condition_2$pseudotime, row.names = row.names(KO_2_tree$pseudotime) )

PlotOutput(KO_1_tree_new, KO_2_tree_new, output_solution_cut)

KO_2_cell_pseudo <- data.frame(KO_2_sce$slingPseudotime_1, row.names = KO_2_sce@colData@rownames)
KO_2_node_pseudo_new <- data.frame(row.names(test$condition_2),test$condition_2$pseudotime, test$condition_2$alignment, row.names =row.names(test$condition_2) )
KO_2_node_pseudo <- data.frame( pseudotime = KO_2_tree$pseudotime, row.names = row.names(KO_2_tree$pseudotime))

KO_2_cell_pseudo_new <- pseudo_cell_align(KO_2_cell_pseudo , KO_2_node_pseudo_new, KO_2_node_pseudo, window)
KO_2_cell_pseudo_new <- KO_2_cell_pseudo_new[order(match(row.names(KO_2_cell_pseudo_new), row.names(KO_2_cell_pseudo))),]

KO_1_cell_pseudo <- data.frame(KO_1_sce$slingPseudotime_1, row.names = KO_1_sce@colData@rownames)
KO_1_node_pseudo_new <- data.frame(row.names(test$condition_1),test$condition_1$pseudotime, test$condition_1$align, row.names =row.names(test$condition_1) )
KO_1_node_pseudo <- data.frame(pseudotime = KO_1_tree$pseudotime, row.names = row.names(KO_1_tree$pseudotime))

KO_1_cell_pseudo_new <- pseudo_cell_align(KO_1_cell_pseudo , KO_1_node_pseudo_new, KO_1_node_pseudo, window)
KO_1_cell_pseudo_new <- KO_1_cell_pseudo_new[order(match(row.names(KO_1_cell_pseudo_new), row.names(KO_1_cell_pseudo))),]

#Add metadata
KO_2_sce$oldPseudotime <- KO_2_sce$slingPseudotime_1
KO_2_sce$newPseudotime <- KO_2_cell_pseudo_new$pseudotime
KO_1_sce$oldPseudotime <- KO_1_sce$slingPseudotime_1
KO_1_sce$newPseudotime <- KO_1_cell_pseudo_new$pseudotime
KO_2_sce$Status <- as.factor(KO_2_cell_pseudo_new$status)
KO_1_sce$Status <- as.factor(KO_1_cell_pseudo_new$status)

KO_2_obj$oldPseudotime <- KO_2_sce$slingPseudotime_1
KO_2_obj$newPseudotime <- KO_2_cell_pseudo_new$pseudotime
KO_1_obj$oldPseudotime <- KO_1_sce$slingPseudotime_1
KO_1_obj$newPseudotime <- KO_1_cell_pseudo_new$pseudotime
KO_2_obj$Status <- as.factor(KO_2_cell_pseudo_new$status)
KO_1_obj$Status <- as.factor(KO_1_cell_pseudo_new$status)

