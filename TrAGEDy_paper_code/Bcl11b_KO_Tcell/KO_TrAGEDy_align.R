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
setwd("~/R")

genesV2 <- read.csv("/Users/rosslaidlaw/R/TrAGEDy_V2/general/human_to_mouse_seurat_cellcycle.csv")

cell_cycle <- readxl::read_xlsx("/Users/rosslaidlaw/R/TrAGEDy_V2/general/mgi_cellCycle_GO.xlsx")

KO_1_obj <- readRDS("/Users/rosslaidlaw/R/TrAGEDy_V2/Tcell_Bcl11b_KO/objects/KO_1_obj.rds")
KO_2_obj <- readRDS("/Users/rosslaidlaw/R/TrAGEDy_V2/Tcell_Bcl11b_KO/objects/KO_2_obj.rds")
features <- read.delim("/Users/rosslaidlaw/R/TrAGEDy_V2/Tcell_Bcl11b_KO/TrAGEDy_feature_space_with_rearrangeTCR.csv", sep = ",")[,2]

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
DimPlot(KO_1_obj, reduction="phate", label = T, group.by = "cell_type")

plot(reducedDims(KO_1_sce)$PHATE, col = brewer.pal(9,'Set1')[KO_1_sce$cell_type], pch=16, asp = 1)
lines(KO_1_sling, lwd=2, type = 'lineages', col = 'black')

FeaturePlot(KO_1_obj, reduction="phate", feature = "Mki67") + xlim(-0.05, 0.07) + ylim(-0.03, 0.04)
FeaturePlot(KO_1_obj, reduction="phate", feature = "Ptcra") + xlim(-0.03, 0.055) + ylim(-0.03, 0.04)



start_clus2 <- "Cd34_TSP"
KO_2_sce <- as.SingleCellExperiment(KO_2_obj, assay = "RNA")
KO_2_sce <- slingshot(KO_2_sce, reducedDim = 'PHATE', clusterLabels = KO_2_sce@colData@listData[["cell_type"]], start.clus = start_clus2, use.median = T, dist.method = "simple")
KO_2_sling <- SlingshotDataSet(KO_2_sce)

KO_1_obj$unalign_pseudotime <- KO_1_sce$slingPseudotime_1
KO_2_obj$unalign_pseudotime <- KO_2_sce$slingPseudotime_1


FeaturePlot(KO_1_obj, reduction="phate", feature = "unalign_pseudotime") + xlim(-0.03, 0.055) + ylim(-0.03, 0.04)
FeaturePlot(KO_2_obj, reduction="phate", feature = "unalign_pseudotime") + xlim(-0.04, 0.05) + ylim(-0.03, 0.02)

FeaturePlot(KO_1_obj, reduction="phate", feature = features[3]) + xlim(-0.03, 0.055) + ylim(-0.03, 0.04)

DimPlot(KO_1_obj, reduction = "phate", label = T)
DimPlot(KO_2_obj, reduction = "phate", label = T)

FeaturePlot(KO_1_obj, reduction="phate", feature = "Atad2") + xlim(-0.03, 0.055) + ylim(-0.03, 0.04)


plot(reducedDims(KO_1_sce)$PHATE, col = brewer.pal(9,'Set1')[KO_1_sce$cell_type], pch=16, asp = 1)
lines(SlingshotDataSet(KO_1_sce), lwd=2, type = 'lineages', col = 'black')

plot(reducedDims(KO_2_sce)$PHATE, col = brewer.pal(9,'Set1')[KO_2_sce$cell_type], pch=16, asp = 1)
lines(SlingshotDataSet(KO_2_sce), lwd=2, type = 'lineages', col = 'black')

pseudo_end <- min(c(max(KO_1_sce$slingPseudotime_1, KO_2_sce$slingPseudotime_1)))
window <- pseudo_end / 90

source("~/R/Scripts/Methods/own_method_functions.R")

#The node name cannot have underscores in it
KO_1_tree <- nodePseudotime(KO_1_sce,"slingPseudotime_1","cell_type", 100, "KO1")
KO_2_tree <- nodePseudotime(KO_2_sce,"slingPseudotime_1","cell_type", 100, "KO2")

KO_2_node_exp_mtx <- nodeExpressionEstimate(KO_2_sce@assays@data@listData$logcounts, KO_2_tree, window, adjust.window = T)
KO_1_node_exp_mtx <- nodeExpressionEstimate(KO_1_sce@assays@data@listData$logcounts, KO_1_tree, window, adjust.window = T)

KO_2_node_exp_mtx  <- KO_2_node_exp_mtx[ features, ]
KO_1_node_exp_mtx <- KO_1_node_exp_mtx[ features, ]

penalty_mtx_cut <- dis_mtx_calculator(as.matrix(KO_1_node_exp_mtx), as.matrix(KO_2_node_exp_mtx), "spearman")

source("~/R/Scripts/Methods/own_method_functions.R")
x <- bootstrap_pathfind(sequence_1 = as.matrix(KO_1_node_exp_mtx), sequence_2 = as.matrix(KO_2_node_exp_mtx)
                        , similarity_method= "spearman", threshold_method = "mean")

output_solution_cut <- cut_deviate(x[[1]], penalty_mtx_cut, method = "mean")

pdf("/Users/rosslaidlaw/R/TrAGEDy_V2/Tcell_Bcl11b_KO/plots/PlotAlignment_KO.pdf", width = 7.75)
PlotAlignment(output_solution_cut, penalty_mtx_cut)
dev.off()

dev.off()
PlotAlignment(output_solution_cut, penalty_mtx_cut)


alignment <- output_solution_cut

pdf("/Users/rosslaidlaw/R/TrAGEDy_V2/Tcell_Bcl11b_KO/plots/Plotoutput_unaligned_KO.pdf")
PlotOutput(KO_1_tree, KO_2_tree, alignment)
dev.off()

test <- chunk_node(KO_1_tree$node_pseudotime, KO_2_tree$node_pseudotime, output_solution_cut)

KO_1_tree_aligned <- KO_1_tree
KO_1_tree_aligned$node_pseudotime <- test$condition_1$pseudotime
names(KO_1_tree_aligned$node_pseudotime) <- row.names(test$condition_1)

KO_2_tree_aligned <- KO_2_tree
KO_2_tree_aligned$node_pseudotime <- test$condition_2$pseudotime
names(KO_2_tree_aligned$node_pseudotime) <- row.names(test$condition_2)

pdf("/Users/rosslaidlaw/R/TrAGEDy_V2/Tcell_Bcl11b_KO/plots/Plotoutput_aligned_KO.pdf")
PlotOutput(KO_1_tree_aligned, KO_2_tree_aligned, alignment)
dev.off()

source("Scripts/methods/own_method_functions.R")
KO_1_cell_pseudo_new <- pseudo_cell_align(KO_1_tree$cell_pseudotime, test$condition_1 , KO_1_tree$node_pseudotime, window)
KO_2_cell_pseudo_new <- pseudo_cell_align(KO_2_tree$cell_pseudotime, test$condition_2 , KO_2_tree$node_pseudotime, window)

PlotOutput(KO_1_tree_aligned, KO_2_tree_aligned, output_solution_cut) + theme(legend.text=element_text(size=15), axis.text.y = element_text(size = 20), axis.text.x = element_text(size = 20), axis.title.y = element_text(size = 20),
                                                                  axis.title.x = element_text(size = 20))

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

saveRDS(KO_1_obj,"/Users/rosslaidlaw/R/TrAGEDy_V2/Tcell_Bcl11b_KO/objects/KO_1_align.rds")
saveRDS(KO_2_obj,"/Users/rosslaidlaw/R/TrAGEDy_V2/Tcell_Bcl11b_KO/objects/KO_2_align.rds")




