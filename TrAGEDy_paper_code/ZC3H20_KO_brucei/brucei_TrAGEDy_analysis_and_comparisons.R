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
library(coin)
setwd("~/R")
source("Scripts/methods/own_method_functions.R")

emma_tradeseq_de <- read.delim("Emma_data_DE/emma_tradeseq_de.csv", sep = ",")

WT_01_obj <- readRDS("emma_KO_WT_brucei_data/Emma Annotation/WT_01_aligned.rds")
WT_02_obj <- readRDS("emma_KO_WT_brucei_data/Emma Annotation/WT_02_aligned.rds")
KO_obj <- readRDS("emma_KO_WT_brucei_data/Emma annotation/ZC3H20_KO.rds")

WT_obj <- merge(WT_01_obj, WT_02_obj)

features <- read.delim("/Users/rosslaidlaw/R/TrAGEDy_results/TrAGEDy_example/WT_ZC3H20_KO_feature_space.csv", sep = ",")[,2]

KO_obj <- subset(KO_obj, cell_type %in% c("LS A", "LS B"))
Idents(KO_obj) <- "cell_type"

WT_sce <- as.SingleCellExperiment(WT_obj, assay = "RNA")

#Use KO clusters
DimPlot(KO_obj, reduction = "phate")
DimPlot(KO_obj, reduction = "phate", group.by = "cell_type")

phate_output <- as.matrix(phate(t(KO_obj@assays$RNA@data[features,]), ndim=10, seed=1))
phate.reduction <- CreateDimReducObject(embeddings = phate_output,
                                        key="PHATE_",
                                        assay="integrated")

KO_obj@reductions$phate <- phate.reduction
DimPlot(KO_obj, reduction="phate", group.by = "cell_type")

KO_sce <- as.SingleCellExperiment(KO_obj, assay = "RNA")
KO_sce <- slingshot(KO_sce, reducedDim = 'PHATE', clusterLabels = KO_sce@colData@listData[["cell_type"]], start.clus = "LS A")
KO_sling <- SlingshotDataSet(KO_sce)

WT_sce$slingPseudotime_1 <- WT_sce$newPseudotime
WT_sce$newPseudotime <- NULL

pseudo_end <- min(c(max(KO_sce$slingPseudotime_1, WT_sce$slingPseudotime_1)))
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

penalty_mtx_cut <- dis_mtx_calculator(as.matrix(WT_node_exp_mtx), as.matrix(KO_node_exp_mtx), "spearman")

output_solution <- pathfind(penalty_mtx_cut, cut_type = "minimum", method = "mean")
output_solution_cut <- cut_deviate(output_solution, penalty_mtx_cut, method = "mean")

#output_solution_cut$Status <- rep("cut", length(output_solution_cut$Status))
alignment_plot <- PlotAlignment(output_solution_cut, penalty_mtx_cut)


PlotOutput(WT_tree, KO_tree, output_solution_cut) + theme(legend.text=element_text(size=15), axis.text.y = element_text(size = 20), axis.text.x = element_text(size = 20), axis.title.y = element_text(size = 20),
                                                          axis.title.x = element_text(size = 20)) 


test <- chunk_node(WT_tree$pseudotime, KO_tree$pseudotime, output_solution_cut)

WT_tree_new <- WT_tree
WT_tree_new$pseudotime <- data.frame(test$condition_1$pseudotime, row.names = row.names(WT_tree$pseudotime) )

KO_tree_new <- KO_tree
KO_tree_new$pseudotime <- data.frame(test$condition_2$pseudotime, row.names = row.names(KO_tree$pseudotime) )

PlotOutput(WT_tree_new, KO_tree_new, output_solution_cut) + theme(legend.text=element_text(size=15), axis.text.y = element_text(size = 20), axis.text.x = element_text(size = 20), axis.title.y = element_text(size = 20),
                                                                  axis.title.x = element_text(size = 20)) 

#Will run into problem when both align from the beginning and we have to move cells backwards that are already at zero 
KO_cell_pseudo <- data.frame(KO_sce$slingPseudotime_1, row.names = KO_sce@colData@rownames)
KO_node_pseudo_new <- data.frame(row.names(test$condition_2),test$condition_2$pseudotime, test$condition_2$alignment, row.names =row.names(test$condition_2) )
KO_node_pseudo <- data.frame( pseudotime = KO_tree$pseudotime, row.names = row.names(KO_tree$pseudotime))

KO_cell_pseudo_new <- pseudo_cell_align(KO_cell_pseudo , KO_node_pseudo_new, KO_node_pseudo, window)
KO_cell_pseudo_new <- KO_cell_pseudo_new[order(match(row.names(KO_cell_pseudo_new), row.names(KO_cell_pseudo))),]

WT_cell_pseudo <- data.frame(WT_sce$slingPseudotime_1, row.names = WT_sce@colData@rownames)
WT_node_pseudo_new <- data.frame(row.names(test$condition_1),test$condition_1$pseudotime, test$condition_1$align, row.names =row.names(test$condition_1) )
WT_node_pseudo <- data.frame(pseudotime = WT_tree$pseudotime, row.names = row.names(WT_tree$pseudotime))

WT_cell_pseudo_new <- pseudo_cell_align(WT_cell_pseudo , WT_node_pseudo_new, WT_node_pseudo, window)
WT_cell_pseudo_new <- WT_cell_pseudo_new[order(match(row.names(WT_cell_pseudo_new), row.names(WT_cell_pseudo))),]


KO_sce$oldPseudotime <- KO_sce$slingPseudotime_1
KO_sce$newPseudotime <- KO_cell_pseudo_new$pseudotime
WT_sce$oldPseudotime <- WT_sce$slingPseudotime_1
WT_sce$newPseudotime <- WT_cell_pseudo_new$pseudotime
KO_sce$Status <- as.factor(KO_cell_pseudo_new$status)
WT_sce$Status <- as.factor(WT_cell_pseudo_new$status)

KO_obj$oldPseudotime <- KO_sce$slingPseudotime_1
KO_obj$newPseudotime <- KO_cell_pseudo_new$pseudotime
WT_obj$oldPseudotime <- WT_sce$slingPseudotime_1
WT_obj$newPseudotime <- WT_cell_pseudo_new$pseudotime
KO_obj$Status <- as.factor(KO_cell_pseudo_new$status)
WT_obj$Status <- as.factor(WT_cell_pseudo_new$status)

source("Scripts/methods/own_method_functions.R")
output <- TrajDE(list(WT_sce, KO_sce), list(WT_node_pseudo, KO_node_pseudo), output_solution_cut, n_windows = 4, overlap = 1, p_val = 0.05, min.pct = 0.1, logfc = 0.5, all.genes = F, test_use = "wilcox", correct = T)

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

#create data frame


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

upset(as.data.frame(upset_df), sets = rev(names(all_venn)), keep.order = T, order.by = "freq", text.scale = 2, point.size = 3, nintersects = 10, mb.ratio = c(0.6,0.4))


output_all <- TrajDE(list(WT_sce, KO_sce), list(WT_node_pseudo, KO_node_pseudo), path_cut, n_windows = 4, overlap = 1, p_val = 0.05, min.pct = 0.1, logfc = 0.5, all.genes = F, test_use = "wilcox", correct = T, own.genes = all_genes)

gene_output_all <- output_all[[1]]

#The genes in our dataset all start with the prefix "Tbrucei---"
#To make plotting easier, we will strip them from the genes

for (i in 1:length(gene_output)){
  gene_output[[i]][, "gene_name"] <- str_replace(gene_output[[i]][, "gene_name"], "Tbrucei---", "")
  row.names(gene_output[[i]]) <- str_replace(row.names(gene_output[[i]]), "Tbrucei---", "")
}

for (i in 1:length(gene_output_all)){
  gene_output_all[[i]][, "gene_name"] <- str_replace(gene_output_all[[i]][, "gene_name"], "Tbrucei---", "")
  row.names(gene_output_all[[i]]) <- str_replace(row.names(gene_output_all[[i]]), "Tbrucei---", "")
}

#We want to make a heatmap of the top 15 genes (in terms of absolute log2 fold change) that are only significant in one window of comparison

logfc_mtx <- matrix(0, ncol = length(gene_output_all), nrow = length(gene_output_all[[1]]$gene_name), dimnames = list( gene_output_all[[1]]$gene_name , paste0("window_", seq(1,4,1))))
for(i in 1:length(gene_output_all)){
  logfc_mtx[,i] <- gene_output_all[[i]]$logfc
  
}

pval_mtx <- matrix(0, ncol = length(gene_output_all), nrow = length(gene_output_all[[1]]$gene_name), dimnames = list( gene_output_all[[1]]$gene_name , paste0("window_", seq(1,4,1))))
for(i in 1:length(gene_output_all)){
  pval_mtx[which(gene_output_all[[i]]$adj_pval < 0.05 ),i] <- 1
}

#subset genes which are only significant in one window
pval_mtx <- pval_mtx[ which(rowSums(pval_mtx) == 1), ]

#sort the genes only significant in one window by log2fc
top15 <- c()
for(i in 1:length(pval_mtx[1,])){
  sig_genes <- row.names(pval_mtx[which(pval_mtx[,i] == 1),])
  
  current_logfc <- logfc_mtx[sig_genes,i]
  current_logfc <- current_logfc[order(abs(current_logfc), decreasing = T)]
  
  if(length(current_logfc < 15)){
    top15 <- append(top15, names(current_logfc))
    
  }
  
  else{
    top15 <- append(top15, names(current_logfc)[1:15])
    
  }
  
}

top15 <- unique(top15)

top15 <- str_replace(top15, "Tbrucei---", "")

#Finally, we can plot the resulting gene list
PlotWindowHeatmap(gene_output_all, gene_output, top15)

