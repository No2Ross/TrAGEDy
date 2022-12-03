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

WT_sce <- readRDS("TrAGEDy_results/WT_vs_KO_brucei/github_files/WT_sce_merged.rds")
KO_sce <- readRDS("TrAGEDy_results/WT_vs_KO_brucei/github_files/KO_sce.rds")

features <- read.delim("TrAGEDy_results/WT_vs_KO_brucei/github_files/WT_ZC3H20_KO_feature_space.csv", sep = ",")[,2]

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

penalty_mtx <- dis_mtx_calculator(as.matrix(WT_node_exp_mtx), as.matrix(KO_node_exp_mtx), "spearman")

source("Scripts/methods/own_method_functions.R")
path_uncut <- pathfind(penalty_mtx, cut_type = "minimum", method = "mean")
path_cut <- cut_deviate(path_uncut, penalty_mtx, method = "mean")

#path_cut$Status <- rep("cut", length(path_cut$Status))
PlotAlignment(path_cut, penalty_mtx)

source("Scripts/methods/own_method_functions.R")
PlotOutput(WT_tree, KO_tree, path_cut) 

source("Scripts/methods/own_method_functions.R")
test <- chunk_node(WT_tree$pseudotime, KO_tree$pseudotime, path_cut)

WT_tree_new <- WT_tree
WT_tree_new$pseudotime <- data.frame(test$condition_1$pseudotime, row.names = row.names(WT_tree$pseudotime) )

KO_tree_new <- KO_tree
KO_tree_new$pseudotime <- data.frame(test$condition_2$pseudotime, row.names = row.names(KO_tree$pseudotime) )


PlotOutput(WT_tree_new, KO_tree_new , path_cut)

#Will run into problem when both align from the beginning and we have to move cells backwards that are already at zero 
source("Scripts/methods/own_method_functions.R")
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

#Find fold change of condition 1 in relation to condition 2


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
output <- TrajDE(list(WT_sce, KO_sce), list(WT_node_pseudo, KO_node_pseudo), path_cut, n_windows = 4, overlap = 1, p_val = 0.05, min.pct = 0.1, logfc = 0.5, all.genes = F, test_use = "wilcox", correct = T)

gene_output <- output[[1]]

all_genes <- c()
for(i in gene_output){
  all_genes <- append(all_genes,i$gene_name)
}

all_genes <- unique(all_genes)

#Plot Heatmap
output_all <- TrajDE(list(WT_sce, KO_sce), list(WT_node_pseudo, KO_node_pseudo), path_cut, n_windows = 4, overlap = 1, p_val = 0.05, min.pct = 0.1, logfc = 0.5, all.genes = F, test_use = "wilcox", correct = T, own.genes = all_genes)

gene_output_all <- output_all[[1]]

for (i in 1:length(gene_output)){
  gene_output[[i]][, "gene_name"] <- str_replace(gene_output[[i]][, "gene_name"], "Tbrucei---", "")
  row.names(gene_output[[i]]) <- str_replace(row.names(gene_output[[i]]), "Tbrucei---", "")
}

for (i in 1:length(gene_output_all)){
  gene_output_all[[i]][, "gene_name"] <- str_replace(gene_output_all[[i]][, "gene_name"], "Tbrucei---", "")
  row.names(gene_output_all[[i]]) <- str_replace(row.names(gene_output_all[[i]]), "Tbrucei---", "")
}

#join output list into logfc matrix

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

#sort the genes only significant in one window by logfc

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


source("Scripts/methods/own_method_functions.R")
PlotWindowHeatmap(gene_output_all, gene_output, top15)

