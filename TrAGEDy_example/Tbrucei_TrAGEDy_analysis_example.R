library(slingshot)
library(RColorBrewer)
library(Seurat)
library(SingleCellExperiment)
library(phateR)
library(reticulate)
use_python("/usr/local/bin/python3.8")
reticulate::import("phate")
library(ggplot2)
library(dplyr)
library(stats)
library(stringr)
library(rgl)
source("/path/to/TrAGEDy_functions.R")

WT_sce <- readRDS("path/to/WT_sce_merged.rds")
KO_sce <- readRDS("path/to/KO_sce.rds")

features <- read.delim("path/to/WT_ZC3H20_KO_feature_space.csv", sep = ",")[,2]

#Get the window size in which to estimate gene expression values for the Interpolated points
pseudo_end <- min(c(max(KO_sce$slingPseudotime_1, WT_sce$slingPseudotime_1)))
window <- pseudo_end / 45

WT_cell_pseudotime <- matrix(WT_sce$slingPseudotime_1, dimnames =list(WT_sce@colData@rownames))
KO_cell_pseudotime <- matrix(KO_sce$slingPseudotime_1, dimnames =list(KO_sce@colData@rownames))
WT_ID <- data.frame(WT_sce$cell_type, row.names =WT_sce@colData@rownames)
KO_ID <- data.frame(KO_sce$cell_type, row.names =KO_sce@colData@rownames)

#Create Interpolated points across pseudotime 
WT_tree <- nodePseudotime(WT_cell_pseudotime,WT_ID, 50, "WT")
KO_tree <- nodePseudotime(KO_cell_pseudotime,KO_ID, 50, "KO")

KO_cell_pseudo <- data.frame("ID" = KO_sce@colData@rownames, "pseudo" = KO_sce$slingPseudotime_1)
KO_node_pseudo <- data.frame("ID" = row.names(KO_tree$pseudotime), "pseudo" = KO_tree$pseudotime$pseudotime)

WT_cell_pseudo <- data.frame("ID" = WT_sce@colData@rownames, "pseudo" = WT_sce$slingPseudotime_1)
WT_node_pseudo <- data.frame("ID" = row.names(WT_tree$pseudotime), "pseudo" = WT_tree$pseudotime$pseudotime)

KO_node_pseudotime <- matrix(KO_tree$pseudotime$pseudotime , dimnames = list(row.names(KO_tree$pseudotime)), )
WT_node_pseudotime <- matrix(WT_tree$pseudotime$pseudotime , dimnames = list(row.names(WT_tree$pseudotime)), )

#Get gene expression values for the interpolated points
KO_node_exp_mtx <- nodeExpressionEstimate(KO_sce@assays@data@listData$logcounts, KO_node_pseudotime, KO_cell_pseudotime, window, adjust.window = T)
WT_node_exp_mtx <- nodeExpressionEstimate(WT_sce@assays@data@listData$logcounts, WT_node_pseudotime, WT_cell_pseudotime, window, adjust.window = T)

#Subset interpolated poinbts to only contain genes present in your feature space
KO_node_exp_mtx  <- KO_node_exp_mtx[ features, ]
WT_node_exp_mtx <- WT_node_exp_mtx[ features, ]

#Calculate dissimilarity between the conditions
penalty_mtx <- dis_mtx_calculator(as.matrix(WT_node_exp_mtx), as.matrix(KO_node_exp_mtx), "spearman")

#Find optimal path through the dissimilarity matrix then cut and matches that have high dissimilarity
path_uncut <- pathfind(penalty_mtx, cut_type = "minimum", method = "mean")
path_cut <- cut_deviate(path_uncut, penalty_mtx, method = "mean")

#Visualise matching between the different interpolated points
PlotAlignment(path_cut, penalty_mtx)
PlotOutput(WT_tree, KO_tree, path_cut) 

#Adjust the pseudotime of the interpolated points based on the optimal path
adjustedPseudotime <- chunk_node(WT_tree$pseudotime, KO_tree$pseudotime, path_cut)

WT_tree_new <- WT_tree
WT_tree_new$pseudotime <- data.frame(adjustedPseudotime$condition_1$pseudotime, row.names = row.names(WT_tree$pseudotime) )

KO_tree_new <- KO_tree
KO_tree_new$pseudotime <- data.frame(adjustedPseudotime$condition_2$pseudotime, row.names = row.names(KO_tree$pseudotime) )

#Visualise the shift in pseudotime of interpolated points
PlotOutput(WT_tree_new, KO_tree_new , path_cut)

#Use the pseudotime of the interpolated points to change the pseudotime of the individual cells
KO_cell_pseudo <- data.frame(KO_sce$slingPseudotime_1, row.names = KO_sce@colData@rownames)
KO_node_pseudo_new <- data.frame(row.names(adjustedPseudotime$condition_2),adjustedPseudotime$condition_2$pseudotime, adjustedPseudotime$condition_2$alignment, row.names =row.names(adjustedPseudotime$condition_2) )
KO_node_pseudo <- data.frame( pseudotime = KO_tree$pseudotime, row.names = row.names(KO_tree$pseudotime))

KO_cell_pseudo_new <- pseudo_cell_align(KO_cell_pseudo , KO_node_pseudo_new, KO_node_pseudo, window)
KO_cell_pseudo_new <- KO_cell_pseudo_new[order(match(row.names(KO_cell_pseudo_new), row.names(KO_cell_pseudo))),]

WT_cell_pseudo <- data.frame(WT_sce$slingPseudotime_1, row.names = WT_sce@colData@rownames)
WT_node_pseudo_new <- data.frame(row.names(adjustedPseudotime$condition_1),adjustedPseudotime$condition_1$pseudotime, adjustedPseudotime$condition_1$align, row.names =row.names(adjustedPseudotime$condition_1) )
WT_node_pseudo <- data.frame(pseudotime = WT_tree$pseudotime, row.names = row.names(WT_tree$pseudotime))

WT_cell_pseudo_new <- pseudo_cell_align(WT_cell_pseudo , WT_node_pseudo_new, WT_node_pseudo, window)
WT_cell_pseudo_new <- WT_cell_pseudo_new[order(match(row.names(WT_cell_pseudo_new), row.names(WT_cell_pseudo))),]

#Add new metadata back into the objects
KO_sce$oldPseudotime <- KO_sce$slingPseudotime_1
KO_sce$newPseudotime <- KO_cell_pseudo_new$pseudotime
WT_sce$oldPseudotime <- WT_sce$slingPseudotime_1
WT_sce$newPseudotime <- WT_cell_pseudo_new$pseudotime
KO_sce$Status <- as.factor(KO_cell_pseudo_new$status)
WT_sce$Status <- as.factor(WT_cell_pseudo_new$status)

#Perform TrajDE to identify differentially expressed (DE) genes between the conditions
#We chose 4 windows in order to keep resuls comparable with Seurat, which compared 4 clusters
output <- TrajDE(list(WT_sce, KO_sce), list(WT_node_pseudo, KO_node_pseudo), path_cut, n_windows = 4, overlap = 1, p_val = 0.05, min.pct = 0.1, logfc = 0.5, all.genes = F, test_use = "wilcox", correct = T)

gene_output <- output[[1]]

all_genes <- c()
for(i in gene_output){
  all_genes <- append(all_genes,i$gene_name)
}

all_genes <- unique(all_genes)

#We can then visualise Log2 fold change and significance for our genes in a heatmap.

#First, we need to get statistics for all the genes we want to plot by supplying them in the 'own.genes' parameter

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

