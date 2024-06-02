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
library(coin)
setwd("~/R")
source("Scripts/methods/own_method_functions.R")

emma_tradeseq_de <- read.delim("Emma_data_DE/emma_tradeseq_de.csv", sep = ",")

WT_01_obj <- readRDS("/Users/rosslaidlaw/R/TrAGEDy_V2/Tbrucei_Zc3h20_KO/objects/WT_1_align.rds")
WT_02_obj <- readRDS("/Users/rosslaidlaw/R/TrAGEDy_V2/Tbrucei_Zc3h20_KO/objects/WT_2_align.rds")
KO_obj <- readRDS("/Users/rosslaidlaw/R/TrAGEDy_V2/Tbrucei_Zc3h20_KO/objects/KO_raw_brucei.rds")

DimPlot(WT_01_obj, group.by = "seurat_clusters")

WT_obj <- merge(WT_01_obj, WT_02_obj)

features <- read.table("/Users/rosslaidlaw/R/TrAGEDy_V2/Tbrucei_Zc3h20_KO/WT_ZC3H20_KO_feature_space.csv", sep = ",", col.names = 1)
features <- features$X1

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

#The node name cannot have underscores in it
WT_tree <- nodePseudotime(WT_sce,"slingPseudotime_1","cell_type", 50, "WT")
KO_tree <- nodePseudotime(KO_sce,"slingPseudotime_1","cell_type", 50, "KO")

WT_node_exp_mtx <- nodeExpressionEstimate(WT_sce@assays@data@listData$logcounts, WT_tree, window, adjust.window = T)
KO_node_exp_mtx <- nodeExpressionEstimate(KO_sce@assays@data@listData$logcounts, KO_tree, window, adjust.window = T)


WT_node_exp_mtx  <- WT_node_exp_mtx[ features, ]
KO_node_exp_mtx <- KO_node_exp_mtx[ features, ]

penalty_mtx_cut <- dis_mtx_calculator(as.matrix(WT_node_exp_mtx), as.matrix(KO_node_exp_mtx), "spearman")

source("~/R/Scripts/Methods/own_method_functions.R")
x <- bootstrap_pathfind(sequence_1 = as.matrix(WT_node_exp_mtx), sequence_2 = as.matrix(KO_node_exp_mtx)
                        , similarity_method= "spearman", threshold_method = "mean")

output_solution_cut <- cut_deviate(x[[1]], penalty_mtx_cut, method = "mean")

PlotAlignment(output_solution_cut, penalty_mtx_cut)

pdf("/Users/rosslaidlaw/R/TrAGEDy_V2/Tbrucei_Zc3h20_KO/plots/TrAGEDy_heatmap_align.pdf", width = 7.75)
PlotAlignment(output_solution_cut, penalty_mtx_cut)
dev.off()

pdf("/Users/rosslaidlaw/R/TrAGEDy_V2/Tbrucei_Zc3h20_KO/plots/TrAGEDy_dotplot_noalign.pdf")
PlotOutput(WT_tree, KO_tree, output_solution_cut) + theme(legend.text=element_text(size=15), axis.text.y = element_text(size = 20), axis.text.x = element_text(size = 20), axis.title.y = element_text(size = 20),
                                                          axis.title.x = element_text(size = 20))
dev.off()


source("Scripts/methods/own_method_functions.R")
test <- chunk_node(WT_tree$node_pseudotime, KO_tree$node_pseudotime, output_solution_cut)


WT_tree_new <- WT_tree
WT_tree_new$node_pseudotime <- test$condition_1$pseudotime
names(WT_tree_new$node_pseudotime) <- names(WT_tree$node_pseudotime)

KO_tree_new <- KO_tree
KO_tree_new$node_pseudotime <- test$condition_2$pseudotime
names(KO_tree_new$node_pseudotime) <- names(KO_tree$node_pseudotime)

pdf("/Users/rosslaidlaw/R/TrAGEDy_V2/Tbrucei_Zc3h20_KO/plots/TrAGEDy_dotplot_align.pdf")
PlotOutput(WT_tree_new, KO_tree_new, output_solution_cut) + theme(legend.text=element_text(size=15), axis.text.y = element_text(size = 20), axis.text.x = element_text(size = 20), axis.title.y = element_text(size = 20),
                                                                  axis.title.x = element_text(size = 20))
dev.off()


source("Scripts/methods/own_method_functions.R")
KO_cell_pseudo_new_new <- pseudo_cell_align_(KO_tree$cell_pseudotime, test$condition_2 , KO_tree$node_pseudotime, window)
WT_cell_pseudo_new_new <- pseudo_cell_align_(WT_tree$cell_pseudotime, test$condition_1 , WT_tree$node_pseudotime, window)

KO_sce$oldPseudotime <- KO_sce$slingPseudotime_1
KO_sce$newPseudotime <- KO_cell_pseudo_new_new$pseudotime
WT_sce$oldPseudotime <- WT_sce$slingPseudotime_1
WT_sce$newPseudotime <- WT_cell_pseudo_new_new$pseudotime
KO_sce$Status <- as.factor(KO_cell_pseudo_new_new$status)
WT_sce$Status <- as.factor(WT_cell_pseudo_new_new$status)

KO_obj$oldPseudotime <- KO_sce$slingPseudotime_1
KO_obj$newPseudotime <- KO_cell_pseudo_new_new$pseudotime
WT_obj$oldPseudotime <- WT_sce$slingPseudotime_1
WT_obj$newPseudotime <- WT_cell_pseudo_new_new$pseudotime
KO_obj$Status <- as.factor(KO_cell_pseudo_new_new$status)
WT_obj$Status <- as.factor(WT_cell_pseudo_new_new$status)

WT_tree_new$cell_pseudotime <- WT_cell_pseudo_new_new$pseudotime
names(WT_tree_new$cell_pseudotime) <- row.names(WT_cell_pseudo_new_new)

KO_tree_new$cell_pseudotime <- KO_cell_pseudo_new_new$pseudotime
names(KO_tree_new$cell_pseudotime) <- row.names(KO_cell_pseudo_new_new)

tragedy_merge <- merge(WT_obj, y = KO_obj)

saveRDS(tragedy_merge, "/Users/rosslaidlaw/R/TrAGEDy_V2/Tbrucei_Zc3h20_KO/objects/tragedy_merge_obj.rds")


tragedy_start_time <- Sys.time()
source("Scripts/methods/own_method_functions.R")
output <- TrajDE_(list(WT_sce, KO_sce), list(WT_tree_new, KO_tree_new), output_solution_cut, n_windows = 4, overlap = 1, p_val = 0.05, min.pct = 0.1, logfc = 0.5, all.genes = F, test_use = "wilcox", correct = T)
tragedy_end_time <- Sys.time()
time_taken_TrAGEDy <- tragedy_end_time - tragedy_start_time

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

combined_seurat <- readRDS("/Users/rosslaidlaw/R/TrAGEDy_V2/Tbrucei_Zc3h20_KO/objects/WT_ZC3H20_with_phate.rds")

DimPlot(combined_seurat, label = T)

combined_seurat$celltype <- Idents(combined_seurat)

combined_seurat <- subset(combined_seurat, celltype != "SS A" & celltype != "SS B")

DimPlot(combined_seurat, label = T)

combined_seurat$celltype.genotype <- paste(Idents(combined_seurat), combined_seurat$line, sep = "_")
Idents(combined_seurat) <- "celltype.genotype"


seurat_gene_output <- list()

seurat_start_time <- Sys.time()
for (i in 1:length(unique(combined_seurat$celltype))){
  
  current_ident <- unique(as.character(combined_seurat$celltype))[i]
  
  current <- FindMarkers(combined_seurat, ident.1 = paste0(current_ident, "_WT"), ident.2=paste0(current_ident, "_ZC3H20_KO"), logfc.threshold  =0.5, min.pct = 0.1)
  
  current <- subset(current, p_val_adj < 0.05)
  
  seurat_gene_output[[current_ident]] <- current
}
seurat_end_time <- Sys.time()
time_taken_seurat <- seurat_end_time - seurat_start_time


seurat_WT_de <- c()
seurat_KO_de <- c()

for (i in 1:length(seurat_gene_output)){
  current <- seurat_gene_output[[i]]
  
  if(length(current[,1]) != 0){
    WT_subset <- subset(current, avg_log2FC > 0 )
    KO_subset <- subset(current, avg_log2FC < 0 )
    
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

#open tradeseq DE genes
tradeseq_gene_output <- list()
tradeseq_genes_all <- c()

for(i in 1:4){
  current <- read.delim(paste0("/Users/rosslaidlaw/R/TrAGEDy_V2/Tbrucei_Zc3h20_KO/tables/trade_seq_win",i, ".csv" ), sep = ",")
  tradeseq_gene_output[[i]] <- current
  tradeseq_genes_all <- append(tradeseq_genes_all, current$X)
  
}

tradeseq_genes_all <- unique(tradeseq_genes_all)

TrAGEDy_noSeurat_noTradeseq <- str_replace(setdiff(unique_new, c(seurat_unique, tradeseq_genes_all)), "Tbrucei---", "")

TrAGEDy_noSeuratTradeSeq_WT <- str_replace(setdiff(unique_new_WT, c(seurat_unique, tradeseq_genes_all)), "Tbrucei---", "" )
TrAGEDy_noSeuratTradeSeq_KO <- str_replace(setdiff(unique_new_KO, c(seurat_unique, tradeseq_genes_all)), "Tbrucei---", "")

seurat_noTrAGEDy_noTradeseq <- str_replace(setdiff(seurat_unique, c(unique_new, tradeseq_genes_all)), "Tbrucei---", "")

Tradeseq_noTrAGEDy_noSeurat <- str_replace(setdiff(tradeseq_genes_all, c(unique_new, seurat_unique)), "Tbrucei---", "")

#save window information for TrAGEDy
for(i in 1:length(gene_output)){
  current <- gene_output[[i]]
  
  write.csv(current, paste0("/Users/rosslaidlaw/R/TrAGEDy_V2/Tbrucei_Zc3h20_KO/tables/window_",i, ".csv" ))
  
}

#save cluster information for seurat
for(i in 1:length(seurat_gene_output)){
  current <- seurat_gene_output[[i]]
  
  write.csv(current, paste0("/Users/rosslaidlaw/R/TrAGEDy_V2/Tbrucei_Zc3h20_KO/tables/cluster_",names(seurat_gene_output)[i], ".csv" ))
  
}

write.csv(Tradeseq_noTrAGEDy_noSeurat, "/Users/rosslaidlaw/R/TrAGEDy_V2/Tbrucei_Zc3h20_KO/tables/tradeseq_noTrAGEDy_noSeurat.csv")
write.csv(TrAGEDy_noSeurat_noTradeseq, "/Users/rosslaidlaw/R/TrAGEDy_V2/Tbrucei_Zc3h20_KO/tables/TrAGEDy_noTradeSeq_noSeurat.csv")
write.csv(seurat_noTrAGEDy_noTradeseq, "/Users/rosslaidlaw/R/TrAGEDy_V2/Tbrucei_Zc3h20_KO/tables/seurat_noTrAGEDy_noTradeseq.csv")


write.csv(str_replace_all(unique_new , "Tbrucei---", "") , "/Users/rosslaidlaw/R/TrAGEDy_V2/Tbrucei_Zc3h20_KO/tables/TrAGEDy_all.csv")
write.csv(str_replace_all(seurat_unique , "Tbrucei---", ""), "/Users/rosslaidlaw/R/TrAGEDy_V2/Tbrucei_Zc3h20_KO/tables/seurat_all.csv")
write.csv(str_replace_all(tradeseq_genes_all , "Tbrucei---", ""), "/Users/rosslaidlaw/R/TrAGEDy_V2/Tbrucei_Zc3h20_KO/tables/tradeseq_all.csv")


#Get venn diagram of intersection of DE gene lists

library(VennDiagram)
library(RColorBrewer)
myCol <- c("#619cff", "#f8766d", "#00ba38")

# Chart
venn.diagram(
  x = list(unique_new, seurat_unique, tradeseq_genes_all),
  category.names = c("TrAGEDy" , "Seurat" , "Tradeseq"),
  filename = "TrAGEDy_V2/Tbrucei_Zc3h20_KO/plots/venn_diagram_0.5logfc.png",
  lwd = 2,
  lty = 'blank',
  fill = myCol,
  cat.cex = 1
)

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

#Create Upset plot
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

upset_out <- upset(as.data.frame(upset_df), sets = rev(names(all_venn)), keep.order = T, order.by = "freq", text.scale = 2, point.size = 2.2, nintersects = 6)
upset_out
upset_out <- upset(as.data.frame(upset_df), sets = rev(names(all_venn)), keep.order = T, order.by = "freq", text.scale = 2, point.size = 3, nintersects = 10, mb.ratio = c(0.4,0.6))

pdf(file = "TrAGEDy_V2/Tbrucei_Zc3h20_KO/plots/Upset.pdf",
    onefile = TRUE, width = 12, height = 10)
upset_out
dev.off()


#Plot Heatmap
genes_test <- unique_new

output <- TrajDE_(list(WT_sce, KO_sce), list(WT_tree_new, KO_tree_new), output_solution_cut, n_windows = 4, overlap = 1, p_val = 0.05, min.pct = 0.1, logfc = 0.5, all.genes = F, test_use = "wilcox", correct = T, own.genes = genes_test)

gene_output_all <- output[[1]]

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
  
  if(length(current_logfc) < 15){
    top15 <- append(top15, names(current_logfc))
    
  }
  
  else{
    top15 <- append(top15, names(current_logfc)[1:15])
    
  }
  
}

top15 <- unique(top15)

top15 <- str_replace(top15, "Tbrucei---", "")

logfc_mtx <- logfc_mtx[top15,]

#Finally, we can plot the resulting gene list
PlotWindowHeatmap(gene_output_all, gene_output, top15)

write.csv(logfc_mtx, "TrAGEDy_V2/figures/supplementary/BruceiHeatmapGenes.csv")

pdf(file = "TrAGEDy_V2/Tbrucei_Zc3h20_KO/plots/top15_heatmap.pdf",
    onefile = TRUE, width = 8, height = 11)
PlotWindowHeatmap(gene_output_all, gene_output, top15)
dev.off()

