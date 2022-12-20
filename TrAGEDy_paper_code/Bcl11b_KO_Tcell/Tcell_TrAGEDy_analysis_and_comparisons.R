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

#merge WT
WT_1_obj <- readRDS("path/to/WT/1/WT_1_align.rds")
WT_2_obj <- readRDS("path/to/WT/2/WT_2_align.rds")

WT_1_obj$oldPseudotime <- WT_1_obj$newPseudotime
WT_1_obj$newPseudotime <- NULL

WT_2_obj$oldPseudotime <- WT_2_obj$newPseudotime
WT_2_obj$newPseudotime <- NULL

WT_obj <- merge(WT_1_obj, WT_2_obj)

#merge KO
KO_1_obj <- readRDS("path/to/KO/1/KO_1_align.rds")
KO_2_obj <- readRDS("path/to/KO/2/KO_2_align.rds")

KO_1_obj$oldPseudotime <- KO_1_obj$newPseudotime
KO_1_obj$newPseudotime <- NULL

KO_2_obj$oldPseudotime <- KO_2_obj$newPseudotime
KO_2_obj$newPseudotime <- NULL

KO_obj <- merge(KO_1_obj, KO_2_obj)

Idents(WT_obj)

DefaultAssay(WT_obj) <- "RNA"
DefaultAssay(KO_obj) <- "RNA"

WT_sce <- as.SingleCellExperiment(WT_obj, assay = "RNA")
KO_sce <- as.SingleCellExperiment(KO_obj, assay = "RNA")

features <- read.delim("TrAGEDy_results/Bcl11b_KO/feature_space.csv", sep = ",")[,2]

features <- setdiff(features, unique(c(genesV2$MGI.symbol, cell_cycle$Symbol)))

pseudo_end <- min(c(max(WT_sce$oldPseudotime, KO_sce$oldPseudotime)))
window <- pseudo_end / 90


WT_cell_pseudotime <- matrix(WT_sce$oldPseudotime, dimnames =list(WT_sce@colData@rownames))
KO_cell_pseudotime <- matrix(KO_sce$oldPseudotime, dimnames =list(KO_sce@colData@rownames))
WT_ID <- data.frame(WT_sce$cell_type, row.names =WT_sce@colData@rownames)
KO_ID <- data.frame(KO_sce$cell_type, row.names =KO_sce@colData@rownames)

WT_tree <- nodePseudotime(WT_cell_pseudotime,WT_ID, 100, "WT")
KO_tree <- nodePseudotime(KO_cell_pseudotime,KO_ID, 100, "KO")

KO_cell_pseudo <- data.frame("ID" = KO_sce@colData@rownames, "pseudo" = KO_sce$oldPseudotime)
KO_node_pseudo <- data.frame("ID" = row.names(KO_tree$pseudotime), "pseudo" = KO_tree$pseudotime$pseudotime)

WT_cell_pseudo <- data.frame("ID" = WT_sce@colData@rownames, "pseudo" = WT_sce$oldPseudotime)
WT_node_pseudo <- data.frame("ID" = row.names(WT_tree$pseudotime), "pseudo" = WT_tree$pseudotime$pseudotime)

KO_node_pseudotime <- matrix(KO_tree$pseudotime$pseudotime , dimnames = list(row.names(KO_tree$pseudotime)), )
WT_node_pseudotime <- matrix(WT_tree$pseudotime$pseudotime , dimnames = list(row.names(WT_tree$pseudotime)), )

KO_node_exp_mtx <- nodeExpressionEstimate(KO_sce@assays@data@listData$logcounts, KO_node_pseudotime, KO_cell_pseudotime, window, adjust.window = T)

WT_node_exp_mtx <- nodeExpressionEstimate(WT_sce@assays@data@listData$logcounts, WT_node_pseudotime, WT_cell_pseudotime, window, adjust.window = T)

KO_node_exp_mtx  <- KO_node_exp_mtx[ features, ]
WT_node_exp_mtx <- WT_node_exp_mtx[ features, ]

#Calculate dissimilarity matrix
penalty_mtx_cut <- dis_mtx_calculator(as.matrix(WT_node_exp_mtx), as.matrix(KO_node_exp_mtx), "spearman")

#Find optimal path
output_solution <- pathfind(penalty_mtx_cut, cut_type = "minimum", method = "mean")
output_solution_cut <- cut_deviate(output_solution, penalty_mtx_cut, method = "mean")

PlotAlignment(output_solution_cut, penalty_mtx_cut)

 
PlotOutput(WT_tree, KO_tree, output_solution_cut) + theme(legend.text=element_text(size=15), axis.text.y = element_text(size = 20), axis.text.x = element_text(size = 20), axis.title.y = element_text(size = 20),
                                                          axis.title.x = element_text(size = 20)) 

#Align pseudotime of interpolated points
test <- chunk_node(WT_tree$pseudotime, KO_tree$pseudotime, output_solution_cut)

WT_tree_new <- WT_tree
WT_tree_new$pseudotime <- data.frame(test$condition_1$pseudotime, row.names = row.names(WT_tree$pseudotime) )

KO_tree_new <- KO_tree
KO_tree_new$pseudotime <- data.frame(test$condition_2$pseudotime, row.names = row.names(KO_tree$pseudotime) )


PlotOutput(WT_tree_new, KO_tree_new, output_solution_cut) + theme(legend.text=element_text(size=15), axis.text.y = element_text(size = 20), axis.text.x = element_text(size = 20), axis.title.y = element_text(size = 20),
                                                                  axis.title.x = element_text(size = 20)) 


KO_cell_pseudo <- data.frame(KO_sce$oldPseudotime, row.names = KO_sce@colData@rownames)
KO_node_pseudo_new <- data.frame(row.names(test$condition_2),test$condition_2$pseudotime, test$condition_2$alignment, row.names =row.names(test$condition_2) )
KO_node_pseudo <- data.frame( pseudotime = KO_tree$pseudotime, row.names = row.names(KO_tree$pseudotime))

KO_cell_pseudo_new <- pseudo_cell_align(KO_cell_pseudo , KO_node_pseudo_new, KO_node_pseudo, window)
KO_cell_pseudo_new <- KO_cell_pseudo_new[order(match(row.names(KO_cell_pseudo_new), row.names(KO_cell_pseudo))),]

WT_cell_pseudo <- data.frame(WT_sce$oldPseudotime, row.names = WT_sce@colData@rownames)
WT_node_pseudo_new <- data.frame(row.names(test$condition_1),test$condition_1$pseudotime, test$condition_1$align, row.names =row.names(test$condition_1) )
WT_node_pseudo <- data.frame(pseudotime = WT_tree$pseudotime, row.names = row.names(WT_tree$pseudotime))

WT_cell_pseudo_new <- pseudo_cell_align(WT_cell_pseudo , WT_node_pseudo_new, WT_node_pseudo, window)
WT_cell_pseudo_new <- WT_cell_pseudo_new[order(match(row.names(WT_cell_pseudo_new), row.names(WT_cell_pseudo))),]

#Add metadata
KO_sce$oldPseudotime <- KO_sce$oldPseudotime
KO_sce$newPseudotime <- KO_cell_pseudo_new$pseudotime
WT_sce$oldPseudotime <- WT_sce$oldPseudotime
WT_sce$newPseudotime <- WT_cell_pseudo_new$pseudotime
KO_sce$Status <- as.factor(KO_cell_pseudo_new$status)
WT_sce$Status <- as.factor(WT_cell_pseudo_new$status)

KO_obj$oldPseudotime <- KO_sce$oldPseudotime
KO_obj$newPseudotime <- KO_cell_pseudo_new$pseudotime
WT_obj$oldPseudotime <- WT_sce$oldPseudotime
WT_obj$newPseudotime <- WT_cell_pseudo_new$pseudotime
KO_obj$Status <- as.factor(KO_cell_pseudo_new$status)
WT_obj$Status <- as.factor(WT_cell_pseudo_new$status)

#Find DE genes
output <- TrajDE(list(WT_sce, KO_sce), list(WT_node_pseudo, KO_node_pseudo), output_solution_cut, n_windows = 7, overlap = 1, p_val = 0.05, min.pct = 0.1, logfc = 0.5, all.genes = F, test_use = "wilcox", correct = T)

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

combined_seurat <- readRDS("path/to/seurat/integration/cell_cycle_UMI_regress_seurat_integrated.rds")

DimPlot(combined_seurat, label = T)

combined_seurat <- subset(combined_seurat, cell_type_seurat != "Interferon_response" & cell_type_seurat != "Rora_T" & 
                            cell_type_seurat != "Outlier" & cell_type_seurat != "early_Rora_T")

DimPlot(combined_seurat, label = T)

DefaultAssay(combined_seurat) <- "RNA"

combined_seurat$celltype.genotype <- paste(Idents(combined_seurat), combined_seurat$condition, sep = "_")
Idents(combined_seurat) <- "celltype.genotype"

seurat_gene_output <- list()

for (i in 1:length(unique(combined_seurat$cell_type_seurat))){
  
  current_ident <- unique(as.character(combined_seurat$cell_type_seurat))[i]
  
  current <- FindMarkers(combined_seurat, ident.1 = paste0(current_ident, "_WT"), ident.2=paste0(current_ident, "_Bcl11b_KO"), logfc.threshold  =0.5)
  
  seurat_gene_output[[current_ident]] <- current
}


seurat_WT_de <- c()
seurat_KO_de <- c()

for (i in 1:length(seurat_gene_output)){
  current <- seurat_gene_output[[i]]
  
  if(length(current[,1]) != 0){
    WT_subset <- subset(current, avg_log2FC > 0 & p_val_adj < 0.05)
    KO_subset <- subset(current, avg_log2FC < 0 & p_val_adj < 0.05)
    
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

TrAGEDy_noSeurat_WT <- setdiff(unique_new_WT, seurat_unique)
TrAGEDy_noSeurat_KO <- setdiff(unique_new_KO, seurat_unique)

seurat_noTrAGEDy_WT <- setdiff(seurat_WT_de, unique_new)
seurat_noTrAGEDy_KO <- setdiff(seurat_KO_de, unique_new)

tradeseq_gene_output <- list()
tradeseq_genes_all <- c()

for(i in 1:7){
  current <- read.delim(paste0("path/to/tradeseq/files/trade_seq_win",i, ".csv" ), sep = ",")
  tradeseq_gene_output[[i]] <- current
  tradeseq_genes_all <- append(tradeseq_genes_all, current$X)
  
}
tradeseq_genes_all <- unique(tradeseq_genes_all)

TrAGEDy_noSeurat_noTradeseq <- setdiff(unique_new, c(seurat_unique, tradeseq_genes_all))

TrAGEDy_noSeuratTradeSeq_WT <- setdiff(unique_new_WT, c(seurat_unique, tradeseq_genes_all))
TrAGEDy_noSeuratTradeSeq_KO <- setdiff(unique_new_KO, c(seurat_unique, tradeseq_genes_all))

seurat_noTrAGEDy_noTradeseq <- setdiff(seurat_unique, c(unique_new,tradeseq_genes_all))

Tradeseq_noTrAGEDy_noSeurat <- setdiff(tradeseq_genes_all, c(unique_new, seurat_unique))

tradeseq_TrAGEDy_intersect <- intersect(unique_new, tradeseq_genes_all)

#save window information for TrAGEDy
for(i in 1:length(gene_output)){
  current <- gene_output[[i]]
  
}

#save cluster information for seurat
for(i in 1:length(seurat_gene_output)){
  current <- seurat_gene_output[[i]]
  
}

#Upset plot creation
for(i in 1:length(gene_output)){
  current <- gene_output[[i]]
  print(g %in% current$gene_name)
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

upset(as.data.frame(upset_df), sets = rev(names(all_venn)), keep.order = T, order.by = "freq", text.scale = 2, point.size = 3, nintersects = 10, mb.ratio = c(0.4,0.6))

#Plot Heatmap
output_all <- TrajDE(list(WT_sce, KO_sce), list(WT_node_pseudo, KO_node_pseudo), output_solution_cut, n_windows = 7, overlap = 1, p_val = 0.05, min.pct = 0.1, logfc = 0.5, all.genes = F, test_use = "wilcox", correct = T, own.genes = unique_new)

gene_output_all <- output_all[[1]]


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

top100 <- unique(top100)

source("Scripts/methods/own_method_functions.R")
PlotWindowHeatmap(gene_output_all, gene_output, top100)



