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
# library(STACAS)
# library(monocle)
library(slingshot)
library(stats)
library(stringi)
library(reshape2)
library(stringr)
library(mclust)
setwd("~/R")
source("Scripts/methods/own_method_functions.R")

WT_01_obj <- readRDS("/Users/rosslaidlaw/R/TrAGEDy_V2/Tbrucei_Zc3h20_KO/objects/WT_01_raw_brucei.rds")
WT_02_obj <- readRDS("/Users/rosslaidlaw/R/TrAGEDy_V2/Tbrucei_Zc3h20_KO/objects/WT_02_raw_brucei.rds")

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

DimPlot(WT_02_obj, reduction = "phate", group.by = "cell_type")
DimPlot(WT_01_obj, reduction = "phate", group.by = "cell_type")


pseudo_end <- max(WT_02_sce$slingPseudotime_1, WT_01_sce$slingPseudotime_1)
window <- pseudo_end / 45

source("~/R/Scripts/Methods/own_method_functions.R")

#The node name cannot have underscores in it
WT_1_tree <- nodePseudotime(WT_01_sce,"slingPseudotime_1","cell_type", 50, "WT1")
WT_2_tree <- nodePseudotime(WT_02_sce,"slingPseudotime_1","cell_type", 50, "WT2")

WT_2_node_exp_mtx <- nodeExpressionEstimate(WT_02_sce@assays@data@listData$logcounts, WT_2_tree, window, adjust.window = T)
WT_1_node_exp_mtx <- nodeExpressionEstimate(WT_01_sce@assays@data@listData$logcounts, WT_1_tree, window, adjust.window = T)



WT_2_node_exp_mtx  <- WT_2_node_exp_mtx[ features, ]
WT_1_node_exp_mtx <- WT_1_node_exp_mtx[ features, ]

penalty_mtx_cut <- dis_mtx_calculator(as.matrix(WT_1_node_exp_mtx), as.matrix(WT_2_node_exp_mtx), "spearman")

source("~/R/Scripts/Methods/own_method_functions.R")
x <- bootstrap_pathfind(sequence_1 = as.matrix(WT_1_node_exp_mtx), sequence_2 = as.matrix(WT_2_node_exp_mtx)
                        , similarity_method= "spearman", threshold_method = "mean")

output_solution_cut <- cut_deviate(x[[1]], penalty_mtx_cut, method = "mean")

PlotAlignment(output_solution_cut, penalty_mtx_cut)

pdf("/Users/rosslaidlaw/R/TrAGEDy_V2/Tbrucei_Zc3h20_KO/plots/PlotAlignment_WT.pdf", width = 7.75)
PlotAlignment(output_solution_cut, penalty_mtx_cut)
dev.off()


pseudo_1 <- WT_2_tree$pseudotime
pseudo_2 <- WT_1_tree$pseudotime
alignment <- output_solution_cut
id1 <- WT_2_tree$ID
id2 <- WT_1_tree$ID

pdf("/Users/rosslaidlaw/R/TrAGEDy_V2/Tbrucei_Zc3h20_KO/plots/PlotOutput_unaligned_WT.pdf")
PlotOutput(WT_1_tree, WT_2_tree, alignment)
dev.off()


test <- chunk_node(WT_1_tree$node_pseudotime, WT_2_tree$node_pseudotime, output_solution_cut)


WT_1_tree_aligned <- WT_1_tree
WT_1_tree_aligned$node_pseudotime <- test$condition_1$pseudotime
names(WT_1_tree_aligned$node_pseudotime) <- row.names(test$condition_1)

WT_2_tree_aligned <- WT_2_tree
WT_2_tree_aligned$node_pseudotime <- test$condition_2$pseudotime
names(WT_2_tree_aligned$node_pseudotime) <- row.names(test$condition_2)

pdf("/Users/rosslaidlaw/R/TrAGEDy_V2/Tbrucei_Zc3h20_KO/plots/PlotOutput_aligned_WT.pdf")
PlotOutput(WT_1_tree_aligned, WT_2_tree_aligned, alignment)
dev.off()

source("Scripts/methods/own_method_functions.R")
WT_1_cell_pseudo_new <- pseudo_cell_align_(WT_1_tree$cell_pseudotime, test$condition_1 , WT_1_tree$node_pseudotime, window)
WT_2_cell_pseudo_new <- pseudo_cell_align_(WT_2_tree$cell_pseudotime, test$condition_2 , WT_2_tree$node_pseudotime, window)


#Add metadata
WT_02_sce$oldPseudotime <- WT_02_sce$slingPseudotime_1
WT_02_sce$newPseudotime <- WT_2_cell_pseudo_new$pseudotime
WT_01_sce$oldPseudotime <- WT_01_sce$slingPseudotime_1
WT_01_sce$newPseudotime <- WT_1_cell_pseudo_new$pseudotime
WT_02_sce$Status <- as.factor(WT_2_cell_pseudo_new$status)
WT_01_sce$Status <- as.factor(WT_1_cell_pseudo_new$status)

WT_02_obj$oldPseudotime <- WT_02_sce$slingPseudotime_1
WT_02_obj$newPseudotime <- WT_2_cell_pseudo_new$pseudotime

WT_01_obj$newPseudotime <- WT_1_cell_pseudo_new$pseudotime
WT_02_obj$Status <- as.factor(WT_2_cell_pseudo_new$status)
WT_01_obj$Status <- as.factor(WT_1_cell_pseudo_new$status)

saveRDS(WT_01_obj,"/Users/rosslaidlaw/R/TrAGEDy_V2/Tbrucei_Zc3h20_KO/objects/WT_1_align.rds")
saveRDS(WT_02_obj,"/Users/rosslaidlaw/R/TrAGEDy_V2/Tbrucei_Zc3h20_KO/objects/WT_2_align.rds")



