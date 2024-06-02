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
library(clustree)
# library(STACAS)
# library(monocle)
library(slingshot)
library(stats)
library(stringi)
library(reshape2)
library(stringr)
library(mclust)
library(coin)
library(biomaRt)
library(metap)

genesV2 <- read.csv("/Users/rosslaidlaw/R/TrAGEDy_V2/general/human_to_mouse_seurat_cellcycle.csv")

cell_cycle <- readxl::read_xlsx("/Users/rosslaidlaw/R/TrAGEDy_V2/general/mgi_cellCycle_GO.xlsx")

WT_1_obj <- readRDS("/Users/rosslaidlaw/R/TrAGEDy_V2/Tcell_Bcl11b_KO/objects/WT_1_obj.rds")
WT_2_obj <- readRDS("/Users/rosslaidlaw/R/TrAGEDy_V2/Tcell_Bcl11b_KO/objects/WT_2_obj.rds")

features <- read.delim("/Users/rosslaidlaw/R/TrAGEDy_V2/Tcell_Bcl11b_KO/TrAGEDy_feature_space_with_rearrangeTCR.csv", sep = ",")[,2]


features <- setdiff(features, unique(c(genesV2$MGI.symbol, cell_cycle$Symbol)))

phate_output <- as.matrix(phate(t(WT_1_obj@assays$RNA@data[features,]), ndim=10, seed=1))

phate.reduction <- CreateDimReducObject(embeddings = phate_output,
                                        key="PHATE_",
                                        assay="RNA")

WT_1_obj@reductions$phate <- phate.reduction

DimPlot(WT_1_obj, reduction = "phate")

start_clus1 <- "Cd34_TSP"

WT_1_sce <- as.SingleCellExperiment(WT_1_obj, assay = "RNA")
WT_1_sce <- slingshot(WT_1_sce, reducedDim = 'PHATE', clusterLabels = WT_1_sce@colData@listData[["cell_type"]], start.clus = start_clus1, use.median = T, dist.method = "simple")
WT_1_sling <- SlingshotDataSet(WT_1_sce)

plot(WT_1_sce@int_colData$reducedDims@listData$PHATE[,1], WT_1_sce@int_colData$reducedDims@listData$PHATE[,2])

phate_output <- as.matrix(phate(t(WT_2_obj@assays$RNA@data[features,]), ndim=10, seed=1))

phate.reduction <- CreateDimReducObject(embeddings = phate_output,
                                        key="PHATE_",
                                        assay="RNA")

WT_2_obj@reductions$phate <- phate.reduction

DimPlot(WT_2_obj, reduction = "phate")
DimPlot(WT_1_obj, reduction = "phate")

pdf("/Users/rosslaidlaw/R/TrAGEDy_V2/Tcell_Bcl11b_KO/plots/WT_1_outlier_PHATE_colour_celltype")
DimPlot(WT_1_obj, reduction = "phate")
dev.off()

WT_1_obj$oldPseudotime <- WT_1_sce$slingPseudotime_1

pdf("/Users/rosslaidlaw/R/TrAGEDy_V2/Tcell_Bcl11b_KO/plots/WT_1_outlier_PHATE_colour_pseudotime")
FeaturePlot(WT_1_obj, reduction = "phate", feature = "oldPseudotime") +ylim(-0.01, 0.04) + xlim(-0.045, 0.055)
dev.off()

start_clus2 <- "Cd3_lo_T"

WT_2_sce <- as.SingleCellExperiment(WT_2_obj, assay = "RNA")
WT_2_sce <- slingshot(WT_2_sce, reducedDim = 'PHATE', clusterLabels = WT_2_sce@colData@listData[["cell_type"]], start.clus = start_clus2, use.median = T, dist.method = "simple")
WT_2_sling <- SlingshotDataSet(WT_2_sce)

hist(WT_1_sce$slingPseudotime_1)

WT_1_sce <- subset(WT_1_sce, , slingPseudotime_1 < 0.1)

DimPlot(WT_1_obj, reduction = "phate", label = T)
DimPlot(WT_2_obj, reduction = "phate", label = T)

pseudo_end <- min(c(max(WT_1_sce$slingPseudotime_1, WT_2_sce$slingPseudotime_1)))
window <- pseudo_end / 90

source("~/R/Scripts/Methods/own_method_functions.R")

#The node name cannot have underscores in it
WT_1_tree <- nodePseudotime(WT_1_sce,"slingPseudotime_1","cell_type", 100, "WT1")
WT_2_tree <- nodePseudotime(WT_2_sce,"slingPseudotime_1","cell_type", 100, "WT2")


WT_2_node_exp_mtx <- nodeExpressionEstimate(WT_2_sce@assays@data@listData$logcounts, WT_2_tree, window, adjust.window = T)
WT_1_node_exp_mtx <- nodeExpressionEstimate(WT_1_sce@assays@data@listData$logcounts, WT_1_tree, window, adjust.window = T)

WT_2_node_exp_mtx  <- WT_2_node_exp_mtx[ features, ]
WT_1_node_exp_mtx <- WT_1_node_exp_mtx[ features, ]

penalty_mtx_cut <- dis_mtx_calculator(as.matrix(WT_1_node_exp_mtx), as.matrix(WT_2_node_exp_mtx), "spearman")

source("~/R/Scripts/Methods/own_method_functions.R")
x <- bootstrap_pathfind(sequence_1 = as.matrix(WT_1_node_exp_mtx), sequence_2 = as.matrix(WT_2_node_exp_mtx)
                    , similarity_method= "spearman", threshold_method = "mean")


output_solution_cut <- cut_deviate(x[[1]], penalty_mtx_cut, method = "mean")

PlotAlignment(output_solution_cut, penalty_mtx_cut)

pdf("/Users/rosslaidlaw/R/TrAGEDy_V2/Tcell_Bcl11b_KO/plots/PlotAlignment_WT.pdf", width = 7.75)
PlotAlignment(output_solution_cut, penalty_mtx_cut)
dev.off()

pseudo_1 <- WT_2_tree$pseudotime
pseudo_2 <- WT_1_tree$pseudotime
alignment <- output_solution_cut
id1 <- WT_2_tree$ID
id2 <- WT_1_tree$ID

pdf("/Users/rosslaidlaw/R/TrAGEDy_V2/Tcell_Bcl11b_KO/plots/PlotOutput_unaligned_WT.pdf")
PlotOutput(WT_1_tree, WT_2_tree, output_solution_cut)
dev.off()


test <- chunk_node(WT_1_tree$node_pseudotime, WT_2_tree$node_pseudotime, output_solution_cut)


WT_1_tree_aligned <- WT_1_tree
WT_1_tree_aligned$node_pseudotime <- test$condition_1$pseudotime
names(WT_1_tree_aligned$node_pseudotime) <- row.names(test$condition_1)

WT_2_tree_aligned <- WT_2_tree
WT_2_tree_aligned$node_pseudotime <- test$condition_2$pseudotime
names(WT_2_tree_aligned$node_pseudotime) <- row.names(test$condition_2)

pdf("/Users/rosslaidlaw/R/TrAGEDy_V2/Tcell_Bcl11b_KO/plots/PlotOutput_aligned_WT.pdf")
PlotOutput(WT_1_tree_aligned, WT_2_tree_aligned, output_solution_cut)
dev.off()

WT_1_cell_pseudo_new <- pseudo_cell_align_(WT_1_tree$cell_pseudotime, test$condition_1 , WT_1_tree$node_pseudotime, window)
WT_2_cell_pseudo_new <- pseudo_cell_align_(WT_2_tree$cell_pseudotime, test$condition_2 , WT_2_tree$node_pseudotime, window)


#Add metadata
WT_2_sce$oldPseudotime <- WT_2_sce$slingPseudotime_1
WT_2_sce$newPseudotime <- WT_2_cell_pseudo_new$pseudotime
WT_1_sce$oldPseudotime <- WT_1_sce$slingPseudotime_1
WT_1_sce$newPseudotime <- WT_1_cell_pseudo_new$pseudotime
WT_2_sce$Status <- as.factor(WT_2_cell_pseudo_new$status)
WT_1_sce$Status <- as.factor(WT_1_cell_pseudo_new$status)

WT_2_obj$oldPseudotime <- WT_2_sce$slingPseudotime_1
WT_2_obj$newPseudotime <- WT_2_cell_pseudo_new$pseudotime
WT_1_obj <- subset(WT_1_obj, oldPseudotime < 0.1)

WT_1_obj$newPseudotime <- WT_1_cell_pseudo_new$pseudotime
WT_2_obj$Status <- as.factor(WT_2_cell_pseudo_new$status)
WT_1_obj$Status <- as.factor(WT_1_cell_pseudo_new$status)

saveRDS(WT_1_obj,"/Users/rosslaidlaw/R/TrAGEDy_V2/Tcell_Bcl11b_KO/objects/WT_1_align.rds")
saveRDS(WT_2_obj,"/Users/rosslaidlaw/R/TrAGEDy_V2/Tcell_Bcl11b_KO/objects/WT_2_align.rds")

plot(WT_1_obj$oldPseudotime, WT_1_obj$newPseudotime)

FeaturePlot(WT_1_obj, reduction = "phate", feature = "oldPseudotime") +ylim(-0.01, 0.04) + xlim(-0.045, 0.055)



