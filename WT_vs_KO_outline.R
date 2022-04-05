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
library(clustree)
library(RColorBrewer)
library(SingleCellExperiment)
library(slingshot)
library(stats)
library(stringi)
library(stringr)

#load in seurat objects
WT_obj <- readRDS("file/path/object.rds")
KO_obj <- readRDS("file/path/object.rds")

#this is where you select the features that defin both processes. 
#Usually what i do is find the DE genes for each of the clusters for the two conditions and find the intersection of these two lists
features <- intersect(obj_1_de_genes, obj_2_de_genes)

#Trajectory inference stages
phate_output <- as.matrix(phate(t(WT_obj@assays$RNA@data[features,]), ndim=2, seed=1))
colnames(phate_output)<- paste0("PHATE_", 1:ncol(phate_output))
phate.reduction <- CreateDimReducObject(embeddings = phate_output,
                                        key="PHATE_",
                                        assay="integrated")

WT_obj@reductions$phate <- phate.reduction
DimPlot(WT_obj, reduction="phate")
WT_sce <- as.SingleCellExperiment(WT_obj, assay = "RNA")
WT_sce <- slingshot(WT_sce, reducedDim = 'PHATE', clusterLabels = WT_sce@colData@listData[["insert metadata name for cluster"]], start.clus = "insert starting cluster name")
WT_sling <- SlingshotDataSet(WT_sce)

phate_output <- as.matrix(phate(t(KO_obj@assays$RNA@data[features,]), ndim=2, seed=1))
phate.reduction <- CreateDimReducObject(embeddings = phate_output,
                                        key="PHATE_",
                                        assay="integrated")

KO_obj@reductions$phate <- phate.reduction
DimPlot(KO_obj, reduction="phate")
KO_sce <- as.SingleCellExperiment(KO_obj, assay = "RNA")
KO_sce <- slingshot(KO_sce, reducedDim = 'PHATE', clusterLabels = KO_sce@colData@listData[["insert metadata name for cluster"]], start.clus = "insert starting cluster name")
KO_sling <- SlingshotDataSet(KO_sce)


#Create meta nodes
source("Scripts/methods/own_method_functions.R")
wt_tree <- slingshot_node_maker(WT_sce, WT_sling, c(n_nodes) ,"X", dims=2, "cluster")
ko_tree <- slingshot_node_maker(KO_sce, KO_sling, c(n_nodes) ,"Y", dims=2, "cluster")

#not a set in stone decision, can change what you want the window to be
pseudo_end <- min(c(max(KO_sce$slingPseudotime_1, WT_sce$slingPseudotime_1)))
window <- pseudo_end / (n_nodes - (n_nodes/10) )


source("Scripts/methods/own_method_functions.R")
WT_pseudobulk <- slingshot_pseudo_bulk(WT_sce@assays@data@listData$counts, WT_sce, WT_sling, wt_tree)
KO_pseudobulk <- slingshot_pseudo_bulk(KO_sce@assays@data@listData$counts, KO_sce, KO_sling, ko_tree)

# remove features which aren't in the present feature list and order 
WT_pseudobulk$lineage_1_pseudobulk <- WT_pseudobulk$lineage_1_pseudobulk[ which( row.names(WT_pseudobulk$lineage_1_pseudobulk) %in% features ), ]
KO_pseudobulk$lineage_1_pseudobulk <- KO_pseudobulk$lineage_1_pseudobulk[ which( row.names(KO_pseudobulk$lineage_1_pseudobulk) %in% features ), ]

#cellalign node exp mtx - not scaled with cellalign way
WT_cell_pseudo <- data.frame(WT_sce@colData@rownames, WT_sce$slingPseudotime_1)
WT_node_pseudo <- data.frame(row.names(wt_tree$lineage_1_pseudotime),wt_tree$lineage_1_pseudotime$pseudotime)
colnames(WT_cell_pseudo) <- c("ID", "pseudo")
colnames(WT_node_pseudo) <- c("ID", "pseudo")

KO_cell_pseudo <- data.frame(KO_sce@colData@rownames, KO_sce$slingPseudotime_1)
KO_node_pseudo <- data.frame(row.names(ko_tree$lineage_1_pseudotime),ko_tree$lineage_1_pseudotime$pseudotime)
colnames(KO_cell_pseudo) <- c("ID", "pseudo")
colnames(KO_node_pseudo) <- c("ID", "pseudo")

WT_node_exp_mtx <- as.matrix(cell_align_node_GEV(WT_sce@assays@data@listData$counts, WT_node_pseudo, WT_cell_pseudo , window))
KO_node_exp_mtx <- as.matrix(cell_align_node_GEV(KO_sce@assays@data@listData$counts, KO_node_pseudo, KO_cell_pseudo , window))

WT_node_exp_mtx <- WT_node_exp_mtx[[1]]
KO_node_exp_mtx <- KO_node_exp_mtx[[1]]

WT_node_exp_mtx  <- WT_node_exp_mtx[ features, ]
KO_node_exp_mtx <- KO_node_exp_mtx[ features, ]

WT_pseudobulk$lineage_1_pseudobulk <- WT_node_exp_mtx
KO_pseudobulk$lineage_1_pseudobulk <- KO_node_exp_mtx

for (lineage in 1:1){
  penalty_mtx <- dis_mtx_calculator(as.matrix(WT_pseudobulk$lineage_1_pseudobulk), as.matrix(KO_pseudobulk$lineage_1_pseudobulk), "spearman")
  
}

penalty_mtx_cut <- penalty_mtx

#don't need to run this, just visulisation
mtx <- matrix(nrow=1, ncol=3)
for (i in 1:length(penalty_mtx_cut[,1])){
  
  for (j in 1:length(penalty_mtx_cut[1,])){
    mtx <- rbind(mtx, c(i , j , penalty_mtx_cut[i,j]) )
    
  }
}
mtx <- as.data.frame(mtx[-1,])
colnames(mtx) <- c("KO", "WT", "Score")
options(rgl.useNULL=TRUE)
options(rgl.printRglwidget = TRUE)
plot3d(mtx$KO, mtx$WT, mtx$Score)
zlim <- range(penalty_mtx_cut)
zlen <- zlim[2] - zlim[1] + 1

colorlut <- terrain.colors(zlen) # height color lookup table
colorlut <- "#0066ff"

col <- colorlut[ penalty_mtx_cut - zlim[1] + 1 ]
col_problem <- penalty_mtx_cut - zlim[1] + 1
col_notproblem <- penalty_mtx_cut-zlim[1] + 1

Y_length <- length(penalty_mtx_cut[,1])
X_length <- length(penalty_mtx_cut[1,])

surface3d(seq(1,(Y_length),1), seq(1,(X_length),1), penalty_mtx_cut, col = col)


#find path through matrix
output_solution <- pathfind(penalty_mtx_cut, cut_type = "minimum")
output_solution_cut <- cut_deviate(output_solution, penalty_mtx_cut)

#align pseudotime of nodes
test <- chunk_node(wt_tree$lineage_1_pseudotime, ko_tree$lineage_1_pseudotime, output_solution_cut)

#get new pseudotimes for cells and plot
WT_cell_pseudo <- data.frame(WT_sce$slingPseudotime_1, row.names = WT_sce@colData@rownames)
WT_node_pseudo_new <- data.frame(row.names(test$condition_1),test$condition_1$pseudotime, test$condition_1$alignment, row.names =row.names(test$condition_1) )
WT_node_pseudo <- data.frame( wt_tree$lineage_1_pseudotime$pseudotime, row.names = row.names(wt_tree$lineage_1_pseudotime))

WT_cell_pseudo_new <- pseudo_cell_align(WT_cell_pseudo , WT_node_pseudo_new, WT_node_pseudo, window)
WT_cell_pseudo_oldV <- pseudo_align_cells(WT_cell_pseudo , WT_node_pseudo_new, WT_node_pseudo)
WT_cell_pseudo_oldV <- WT_cell_pseudo_oldV[order(match(WT_cell_pseudo_oldV[,1], WT_cell_pseudo[,1])),]
WT_cell_pseudo_new <- WT_cell_pseudo_new[order(match(row.names(WT_cell_pseudo_new), row.names(WT_cell_pseudo))),]

KO_cell_pseudo <- data.frame(KO_sce$slingPseudotime_1, row.names = KO_sce@colData@rownames)
KO_node_pseudo_new <- data.frame(row.names(test$condition_2),test$condition_2$pseudotime, test$condition_2$align, row.names =row.names(test$condition_2) )
KO_node_pseudo <- data.frame(ko_tree$lineage_1_pseudotime$pseudotime, row.names = row.names(ko_tree$lineage_1_pseudotime))

KO_cell_pseudo_new <- pseudo_cell_align(KO_cell_pseudo , KO_node_pseudo_new, KO_node_pseudo, window)
KO_cell_pseudo_new <- KO_cell_pseudo_new[order(match(row.names(KO_cell_pseudo_new), row.names(KO_cell_pseudo))),]
KO_cell_pseudo_oldV <- pseudo_align_cells(KO_cell_pseudo , KO_node_pseudo_new, KO_node_pseudo)
KO_cell_pseudo_oldV <- KO_cell_pseudo_oldV[order(match(KO_cell_pseudo_oldV[,1], KO_cell_pseudo[,1])),]

#Correlation of original and new pseudotime values
cor(cbind(WT_cell_pseudo_new[,1], WT_cell_pseudo[,1]), method = "spearman")
cor(cbind(KO_cell_pseudo_new[,1], KO_cell_pseudo[,1]), method = "spearman")

WT_frame <- data.frame("cell_type" = WT_sce$cell_type, "old_pseudo" = WT_sce$slingPseudotime_1, "new_pseudo" = WT_cell_pseudo_new$pseudotime)

KO_frame <- data.frame("cell_type" = KO_sce$cell_type, "old_pseudo" = KO_sce$slingPseudotime_1, "new_pseudo" = KO_cell_pseudo_new$pseudotime)

plot_frame <- data.frame("origin" = c( rep("WT_new", length(WT_cell_pseudo_new$pseudotime) ) , rep("KO_new", length(KO_cell_pseudo_new$pseudotime)),rep("WT_old", length(WT_sce$slingPseudotime_1)) , 
                                       rep("KO_old", length(KO_sce$slingPseudotime_1))  ),
                         "pseudotime" = c(WT_cell_pseudo_new$pseudotime, KO_cell_pseudo_new$pseudotime,WT_sce$slingPseudotime_1, KO_sce$slingPseudotime_1),
                         
                         "cell_type" = c(WT_sce$cell_type, KO_sce$cell_type , WT_sce$cell_type, KO_sce$cell_type)
)
plot_frame$origin <- factor(plot_frame$origin, c("WT_old", "KO_old", "WT_new", "KO_new"))
ggplot() + geom_point(data = plot_frame, aes(pseudotime, origin, col = cell_type))




