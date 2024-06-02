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

genesV2 <- read.csv("/Users/rosslaidlaw/R/TrAGEDy_results/TrAGEDy_paper/WT_Bcl11bKO_tcell/human_to_mouse_seurat_cellcycle.csv")

cell_cycle <- readxl::read_xlsx("/Users/rosslaidlaw/R/TrAGEDy_results/TrAGEDy_paper/WT_Bcl11bKO_tcell/mgi_cellCycle_GO.xlsx")

start_time <- Sys.time()

#merge WT
WT_1_obj <- readRDS("/Users/rosslaidlaw/R/TrAGEDy_V2/Tcell_Bcl11b_KO/objects/WT_1_align.rds")
WT_2_obj <- readRDS("/Users/rosslaidlaw/R/TrAGEDy_V2/Tcell_Bcl11b_KO/objects/WT_2_align.rds")

WT_1_obj$oldPseudotime <- WT_1_obj$newPseudotime
WT_1_obj$newPseudotime <- NULL

WT_2_obj$oldPseudotime <- WT_2_obj$newPseudotime
WT_2_obj$newPseudotime <- NULL

WT_obj <- merge(WT_1_obj, WT_2_obj)

#merge KO
KO_1_obj <- readRDS("/Users/rosslaidlaw/R/TrAGEDy_V2/Tcell_Bcl11b_KO/objects/KO_1_align.rds")
KO_2_obj <- readRDS("/Users/rosslaidlaw/R/TrAGEDy_V2/Tcell_Bcl11b_KO/objects/KO_2_align.rds")

KO_1_obj$oldPseudotime <- KO_1_obj$newPseudotime
KO_1_obj$newPseudotime <- NULL

KO_2_obj$oldPseudotime <- KO_2_obj$newPseudotime
KO_2_obj$newPseudotime <- NULL

KO_obj <- merge(KO_1_obj, KO_2_obj)

DefaultAssay(WT_obj) <- "RNA"
DefaultAssay(KO_obj) <- "RNA"

WT_sce <- as.SingleCellExperiment(WT_obj, assay = "RNA")
KO_sce <- as.SingleCellExperiment(KO_obj, assay = "RNA")

features <- read.delim("/Users/rosslaidlaw/R/TrAGEDy_V2/Tcell_Bcl11b_KO/TrAGEDy_feature_space_with_rearrangeTCR.csv", sep = ",")[,2]

features <- setdiff(features, unique(c(genesV2$MGI.symbol, cell_cycle$Symbol)))

pseudo_end <- max(WT_sce$oldPseudotime, KO_sce$oldPseudotime)
window <- pseudo_end / 90


source("Scripts/methods/own_method_functions.R")
WT_tree <- nodePseudotime(WT_sce,"oldPseudotime","cell_type", 100, "WT")

source("Scripts/methods/own_method_functions.R")
KO_tree <- nodePseudotime(KO_sce,"oldPseudotime","cell_type", 100, "KO")

KO_node_exp_mtx <- nodeExpressionEstimate(KO_sce@assays@data@listData$logcounts, KO_tree, window, adjust.window = T)
WT_node_exp_mtx <- nodeExpressionEstimate(WT_sce@assays@data@listData$logcounts, WT_tree, window, adjust.window = T)

KO_node_exp_mtx  <- KO_node_exp_mtx[ features, ]
WT_node_exp_mtx <- WT_node_exp_mtx[ features, ]

row.names(as.matrix(KO_node_exp_mtx)) == row.names(as.matrix(WT_node_exp_mtx))
source("Scripts/methods/own_method_functions.R")
penalty_mtx_cut <- dis_mtx_calculator(as.matrix(WT_node_exp_mtx), as.matrix(KO_node_exp_mtx), "spearman")

source("Scripts/methods/own_method_functions.R")

x <- bootstrap_pathfind(sequence_1 = as.matrix(WT_node_exp_mtx), sequence_2 = as.matrix(KO_node_exp_mtx)
                        , similarity_method= "spearman", threshold_method = "mean")

output_solution_cut <- cut_deviate(x[[1]], penalty_mtx_cut, method = "mean")

#output_solution_cut$Status <- rep("cut", length(output_solution_cut$Status))
pdf(file = "TrAGEDy_V2/Tcell_Bcl11b_KO/plots/PlotAlignment_heatmap.pdf",
    width = 7.75)
alignment_plot <- PlotAlignment(output_solution_cut, penalty_mtx_cut)
dev.off()

pdf(file = "TrAGEDy_V2/Tcell_Bcl11b_KO/plots/PlotAlignment_dot_unaligned.pdf", width = 9)
PlotOutput(WT_tree, KO_tree, output_solution_cut) + theme(legend.text=element_text(size=15), axis.text.y = element_text(size = 20), axis.text.x = element_text(size = 20), axis.title.y = element_text(size = 20),
                                                          axis.title.x = element_text(size = 20))
dev.off()

test <- chunk_node(WT_tree$node_pseudotime, KO_tree$node_pseudotime, output_solution_cut)


WT_tree_new <- WT_tree
WT_tree_new$node_pseudotime <- test$condition_1$pseudotime
names(WT_tree_new$node_pseudotime) <- names(WT_tree$node_pseudotime)

KO_tree_new <- KO_tree
KO_tree_new$node_pseudotime <- test$condition_2$pseudotime
names(KO_tree_new$node_pseudotime) <- names(KO_tree$node_pseudotime)

pdf(file = "TrAGEDy_V2/Tcell_Bcl11b_KO/plots/PlotAlignment_dot_aligned.pdf", width = 9)
PlotOutput(WT_tree_new, KO_tree_new, output_solution_cut) + theme(legend.text=element_text(size=15), axis.text.y = element_text(size = 20), axis.text.x = element_text(size = 20), axis.title.y = element_text(size = 20),
                                                                  axis.title.x = element_text(size = 20))
dev.off()

source("Scripts/methods/own_method_functions.R")
KO_cell_pseudo_new_new <- pseudo_cell_align(KO_tree$cell_pseudotime, test$condition_2 , KO_tree$node_pseudotime, pseudo_end / 90)
WT_cell_pseudo_new_new <- pseudo_cell_align(WT_tree$cell_pseudotime, test$condition_1 , WT_tree$node_pseudotime, pseudo_end / 90)

WT_tree_new$cell_pseudotime <- WT_cell_pseudo_new_new$pseudotime
names(WT_tree_new$cell_pseudotime) <- row.names(WT_cell_pseudo_new_new)

KO_tree_new$cell_pseudotime <- KO_cell_pseudo_new_new$pseudotime
names(KO_tree_new$cell_pseudotime) <- row.names(KO_cell_pseudo_new_new)

plot(KO_cell_pseudo_new_new$pseudotime, KO_tree$cell_pseudotime)
plot(WT_cell_pseudo_new_new$pseudotime, WT_tree$cell_pseudotime)


KO_sce$newPseudotime <- KO_cell_pseudo_new_new$pseudotime
WT_sce$newPseudotime <- WT_cell_pseudo_new_new$pseudotime
KO_sce$Status <- as.factor(KO_cell_pseudo_new_new$status)
WT_sce$Status <- as.factor(WT_cell_pseudo_new_new$status)

KO_obj$newPseudotime <- KO_cell_pseudo_new_new$pseudotime
WT_obj$newPseudotime <- WT_cell_pseudo_new_new$pseudotime
KO_obj$Status <- as.factor(KO_cell_pseudo_new_new$status)
WT_obj$Status <- as.factor(WT_cell_pseudo_new_new$status)


tragedy_merge <- merge(WT_obj, y = KO_obj)

saveRDS(tragedy_merge, "/Users/rosslaidlaw/R/TrAGEDy_V2/Tcell_Bcl11b_KO/objects/tragedy_merge_obj.rds")

#Find DE genes - TrAGEDy

start_time_tragedy <- Sys.time()
output <- TrajDE(list(WT_sce, KO_sce), list(WT_tree_new, KO_tree_new), output_solution_cut, n_windows = 6, 
                  overlap = 1, p_val = 0.05, min.pct = 0.1, logfc = 0.75, all.genes = F, test_use = "wilcox", correct = T)


end_time_tragedy <- Sys.time()
time_taken_TrAGEDy <- end_time_tragedy - start_time_tragedy

gene_output <- output[[1]]

unique_new <- c()
for (i in 1:length(gene_output)) {
  unique_new <- append(unique_new, gene_output[[i]][,"gene_name"])
  
}

unique_new <- unique(unique_new)

#KO side DE genes
unique_new_KO <- c()
for(i in 1:length(gene_output)){
  
  current <- gene_output[[i]]
  current <- subset(current, logfc < 0)
  unique_new_KO <- append(unique_new_KO, current$gene_name)
  
}
unique_new_KO <- unique(unique_new_KO)


#WT side DE genes
unique_new_WT <- c()
for(i in 1:length(gene_output)){
  
  current <- gene_output[[i]]
  current <- subset(current, logfc > 0)
  unique_new_WT <- append(unique_new_WT, current$gene_name)
  
}

unique_new_WT <- unique(unique_new_WT)

window_genes <- list()

for(i in 1:length(gene_output)){
  
  current <- gene_output[[i]]
  window_genes[[i]] <- current$gene_name
  
}


combined_seurat <- readRDS("/Users/rosslaidlaw/R/TrAGEDy_V2/Tcell_Bcl11b_KO/objects/phate_cell_cycle_UMI_regress_seurat_integrated.rds")

DimPlot(combined_seurat, label = T)

DimPlot(combined_seurat, label = T)

DefaultAssay(combined_seurat) <- "RNA"

combined_seurat$celltype.genotype <- paste(Idents(combined_seurat), combined_seurat$condition, sep = "_")
Idents(combined_seurat) <- "celltype.genotype"

seurat_gene_output <- list()

start_time_seurat <- Sys.time()
for (i in 1:length(unique(combined_seurat$cell_type_seurat))){
  
  current_ident <- unique(as.character(combined_seurat$cell_type_seurat))[i]
  
  current <- FindMarkers(combined_seurat, ident.1 = paste0(current_ident, "_WT"), ident.2=paste0(current_ident, "_Bcl11b_KO"), logfc.threshold = 0.75,
                         min.pct = 0.1)
  
  current <- subset(current, p_val_adj < 0.05)
  
  seurat_gene_output[[current_ident]] <- current
}

end_time_seurat <- Sys.time()
time_taken_seurat <- end_time_seurat - start_time_seurat


seurat_WT_de <- c()
seurat_KO_de <- c()

for (i in 1:length(seurat_gene_output)){
  current <- seurat_gene_output[[i]]
  
  if(length(current[,1]) != 0){
    WT_subset <- subset(current, avg_log2FC > 0)
    KO_subset <- subset(current, avg_log2FC < 0)
    
    current_WT <- row.names(WT_subset)
    
    current_KO <- row.names(KO_subset)
    
    seurat_WT_de <- append(seurat_WT_de, current_WT)
    seurat_KO_de <- append(seurat_KO_de, current_KO)
  }
  
}

seurat_WT_de <- unique(seurat_WT_de)
seurat_KO_de <- unique(seurat_KO_de)
seurat_unique <- unique(c(seurat_KO_de, seurat_WT_de))


TrAGEDy_noSeurat <- setdiff(unique_new, seurat_unique)
seurat_noTrAGEDy <- setdiff(seurat_unique, unique_new)

TrAGEDy_noSeurat_WT <- setdiff(unique_new_WT, seurat_unique)
TrAGEDy_noSeurat_KO <- setdiff(unique_new_KO, seurat_unique)

seurat_noTrAGEDy_WT <- setdiff(seurat_WT_de, unique_new)
seurat_noTrAGEDy_KO <- setdiff(seurat_KO_de, unique_new)

#Get Tradeseq results for original pseudotime
tradeseq_gene_output <- list()
tradeseq_genes_all <- c()

for(i in 1:6){
  current <- read.delim(paste0("/Users/rosslaidlaw/R/TrAGEDy_V2/Tcell_Bcl11b_KO/tables/trade_seq_win",i, ".csv" ), sep = ",")
  tradeseq_gene_output[[i]] <- current
  tradeseq_genes_all <- append(tradeseq_genes_all, current$X)
  
}
tradeseq_genes_all <- unique(tradeseq_genes_all)


TrAGEDy_noSeurat_noTradeseq <- setdiff(unique_new, c(seurat_unique, tradeseq_genes_all))

TrAGEDy_noSeuratTradeSeq_WT <- setdiff(unique_new, c(seurat_unique, tradeseq_genes_all))
TrAGEDy_noSeuratTradeSeq_KO <- setdiff(unique_new, c(seurat_unique, tradeseq_genes_all))

seurat_noTrAGEDy_noTradeseq <- setdiff(seurat_unique, c(unique_new,tradeseq_genes_all))

Tradeseq_noTrAGEDy_noSeurat <- setdiff(tradeseq_genes_all, c(unique_new, seurat_unique))

tradeseq_TrAGEDy_intersect <- intersect(unique_new, tradeseq_genes_all)

length(TrAGEDy_noSeurat_noTradeseq)
length(seurat_noTrAGEDy_noTradeseq)
length(Tradeseq_noTrAGEDy_noSeurat)


#save window information for TrAGEDy
for(i in 1:length(gene_output)){
  current <- gene_output[[i]]
  
  write.csv(current, paste0("/Users/rosslaidlaw/R/TrAGEDy_V2/Tcell_Bcl11b_KO/tables/window_",i, ".csv" ))
  
}

#save cluster information for seurat
for(i in 1:length(seurat_gene_output)){
  current <- seurat_gene_output[[i]]
  
  write.csv(current, paste0("/Users/rosslaidlaw/R/TrAGEDy_V2/Tcell_Bcl11b_KO/tables/cluster_",names(seurat_gene_output)[i], ".csv" ))
  
}

write.csv(Tradeseq_noTrAGEDy_noSeurat, "/Users/rosslaidlaw/R/TrAGEDy_V2/Tcell_Bcl11b_KO/tables/tradeseq_noTrAGEDy_noSeurat.csv")
write.csv(TrAGEDy_noSeurat_noTradeseq, "/Users/rosslaidlaw/R/TrAGEDy_V2/Tcell_Bcl11b_KO/tables/TrAGEDy_noTradeSeq_noSeurat.csv")
write.csv(seurat_noTrAGEDy_noTradeseq, "/Users/rosslaidlaw/R/TrAGEDy_V2/Tcell_Bcl11b_KO/tables/seurat_noTrAGEDy_noTradeseq.csv")

write.csv(unique_new  , "/Users/rosslaidlaw/R/TrAGEDy_V2/Tcell_Bcl11b_KO/tables/TrAGEDy_all.csv")
write.csv(seurat_unique , "/Users/rosslaidlaw/R/TrAGEDy_V2/Tcell_Bcl11b_KO/tables/seurat_all.csv")
write.csv(tradeseq_genes_all , "/Users/rosslaidlaw/R/TrAGEDy_V2/Tcell_Bcl11b_KO/tables/tradeseq_all.csv")

#Get Overlap of DE genes

library(VennDiagram)
library(RColorBrewer)
myCol <- c("#619cff", "#f8766d", "#00ba38")

# Chart
venn.diagram(
  x = list(unique_new, seurat_unique, tradeseq_genes_all),
  category.names = c("TrAGEDy" , "Seurat" , "Tradeseq"),
  filename = "TrAGEDy_V2/Tcell_Bcl11b_KO/plots/venn_diagram_0.5logfc.png",
  lwd = 2,
  lty = 'blank',
  fill = myCol,
  cat.cex = 1
)


#Upset plot creation
for(i in 1:length(gene_output)){
  current <- gene_output[[i]]
}

all_venn <- list()


for(i in 1:length(gene_output)){
  all_venn[[ paste0("Window_", i) ]] <- row.names(gene_output[[i]])
}

for(i in 1:length(seurat_gene_output)){
  all_venn[[names(seurat_gene_output)[i]]] <- row.names(seurat_gene_output[[i]])
}

for(i in 1:length(tradeseq_gene_output)){
  all_venn[[ paste0("tradeseq_knot_", i)  ]] <- tradeseq_gene_output[[i]]$X
}


library(UpSetR)

all_genes <- unique(c(unique_new, seurat_unique, tradeseq_genes_all))

upset_df <- matrix(0, nrow = length(all_genes), ncol = length(all_venn), dimnames = list( all_genes , names(all_venn) ) )

for(i in 1:length(all_genes)){
  gene <- all_genes[i]
  gene_add <- rep(0, length(all_venn))
  
  for(j in 1:length(all_venn)){
    
    if(gene %in% all_venn[[j]]){
      gene_add[j] <- gene_add[j] + 1
    }
    
  }
  
  upset_df[gene,] <- gene_add
  
}

upset_out <- upset(as.data.frame(upset_df), sets = rev(names(all_venn)), keep.order = T, order.by = "freq", text.scale = 2, point.size = 3, nintersects = 10, mb.ratio = c(0.4,0.6))
upset_out

pdf(file = "TrAGEDy_V2/Tcell_Bcl11b_KO/plots/Upset.pdf",
    onefile = TRUE, width = 15, height = 9)
upset_out
dev.off()

#Plot Heatmap
output <- TrajDE_(list(WT_sce, KO_sce), list(WT_tree_new, KO_tree_new), output_solution_cut, n_windows = 6, overlap = 1, p_val = 0.05, min.pct = 0.1, logfc = 0.5, all.genes = F, test_use = "wilcox", correct = T, own.genes = genes_test)

gene_output_all <- output[[1]]


#get top 25 genes from each logfc level direction
logfc_matrix <- matrix(0, nrow = length(gene_output_all[[1]]$gene_name), ncol  = length(gene_output_all))
row.names(logfc_matrix) <- row.names(gene_output_all[[1]])
for (i in 1:length(gene_output_all)){
  current <- gene_output_all[[i]]
  logfc_matrix[,i] <- current$logfc
}

logfc_mtx_mean <- rowMeans(logfc_matrix)
names(logfc_mtx_mean) <- row.names(logfc_matrix)
logfc_mtx_mean <- logfc_mtx_mean[order(logfc_mtx_mean, decreasing = T)]

top100 <- names(logfc_mtx_mean[1:25])

logfc_matrix <- logfc_matrix[names(logfc_mtx_mean[1:25]),]

#Get LogFC heatmap ordered
logfc_matrix_save <- logfc_matrix[ order( match ( names(logfc_mtx_mean[1:25]), row.names(logfc_matrix) ) ), ]



logfc_matrix <- matrix(0, nrow = length(gene_output_all[[1]]$gene_name), ncol  = length(gene_output_all))
row.names(logfc_matrix) <- row.names(gene_output_all[[1]])
for (i in 1:length(gene_output_all)){
  current <- gene_output_all[[i]]
  logfc_matrix[,i] <- current$logfc
}

logfc_mtx_mean <- rowMeans(logfc_matrix)
names(logfc_mtx_mean) <- row.names(logfc_matrix)
logfc_mtx_mean <- logfc_mtx_mean[order(logfc_mtx_mean, decreasing = F)]

top100 <- append(top100, names(logfc_mtx_mean[1:25]))

logfc_matrix <- logfc_matrix[names(logfc_mtx_mean[1:25]),]

logfc_matrix <- logfc_matrix[ order( match ( names(logfc_mtx_mean[1:25]), row.names(logfc_matrix) ) ), ]

logfc_matrix_save <- rbind(logfc_matrix_save, logfc_matrix)

write.csv(logfc_matrix_save, "TrAGEDy_V2/figures/supplementary/TcellHeatmapGenes.csv")

top100 <- unique(top100)

write.csv()

source("Scripts/methods/own_method_functions.R")
PlotWindowHeatmap(gene_output_all, gene_output, top100)

pdf(file = "TrAGEDy_V2/Tcell_Bcl11b_KO/plots/top50_condition_heatmap.pdf",
    onefile = TRUE, width = 9, height = 15)
PlotWindowHeatmap(gene_output_all, gene_output, top100) 
dev.off()

