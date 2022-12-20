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

Tcell.combined <- readRDS( "path/to/combined/object/cell_cycle_UMI_regress_seurat_integrated.rds")

Tcell.combined$sample_id <- as.factor(Tcell.combined$sample_id)

Tcell.combined[["percent.mt"]] <- PercentageFeatureSet(Tcell.combined, pattern = "^mt-")

DefaultAssay(Tcell.combined) <- "RNA"

KO_obj <- subset(Tcell.combined, condition == "Bcl11b_KO")

WT_obj <- subset(Tcell.combined, condition == "WT")

DimPlot(KO_obj, reduction = "umap", label = T, group.by = "cell_type")

KO_split <- SplitObject(KO_obj, "orig.ident")

WT_split <- SplitObject(WT_obj, "orig.ident")

WT_1_obj <- WT_split[[1]]
WT_2_obj <- WT_split[[2]]

KO_1_obj <- KO_split[[1]]

DefaultAssay(KO_1_obj) <- "RNA"

KO_2_obj <- KO_split[[2]]

#Get combined variable features list
KO_1_obj <- NormalizeData(KO_1_obj)
KO_2_obj <- NormalizeData(KO_2_obj)
WT_2_obj <- NormalizeData(WT_2_obj)
WT_1_obj <- NormalizeData(WT_1_obj)

KO_1_obj <- FindVariableFeatures(KO_1_obj, selection.method = "vst", nfeatures = 3000)
KO_2_obj <- FindVariableFeatures(KO_2_obj, selection.method = "vst", nfeatures = 3000)

WT_1_obj <- FindVariableFeatures(WT_1_obj, selection.method = "vst", nfeatures = 3000)
WT_2_obj <- FindVariableFeatures(WT_2_obj, selection.method = "vst", nfeatures = 3000)

s.genes <- c(str_to_title(tolower(cc.genes$s.genes)), "Cenpu")
g2m.genes <- c(str_to_title(tolower(cc.genes$g2m.genes)), "Pimreg", "Jpt1")

DefaultAssay(KO_1_obj) <- "RNA"

KO_1_obj <- CellCycleScoring(KO_1_obj, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

DefaultAssay(KO_2_obj) <- "RNA"


KO_2_obj <- NormalizeData(KO_2_obj)

KO_2_obj <- CellCycleScoring(KO_2_obj, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

KO_1_obj <- ScaleData(KO_1_obj, verbose = FALSE , vars.to.regress = c(genesV2$MGI.symbol),do.scale = T, do.center = T, features = row.names(KO_1_obj))
KO_1_obj <- RunPCA(KO_1_obj, npcs = 50, verbose = FALSE)
ElbowPlot(KO_1_obj, ndims = 50)
KO_1_obj <- RunUMAP(KO_1_obj, reduction = "pca", dims = 1:20)
KO_1_obj <- FindNeighbors(KO_1_obj, dims = 1:20)

res_options <- c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0)

for (i in res_options){
  KO_1_obj <- FindClusters(KO_1_obj, res = i)
}

clustree(KO_1_obj, prefix = "RNA_snn_res.", node_size = 10,
         node_alpha = 0.8)

KO_1_obj <- FindClusters(KO_1_obj, res = 0.7)

KO_1_obj <- RenameIdents(KO_1_obj, "0" = "rearrange_TCR", "1" = "Cd3_lo_T", "2" = "Cd3_T", "3" = "Cd34_TSP", "4" = "Rora_T")

KO_1_obj <- AddMetaData(KO_1_obj, KO_1_obj@active.ident, "cell_type")

KO_1_obj <- subset(KO_1_obj, cell_type != "Outlier" & cell_type != "Rora_T")

KO_2_obj <- ScaleData(KO_2_obj, verbose = FALSE ,vars.to.regress = c(genesV2$MGI.symbol), do.scale = T, do.center = T, features = row.names(KO_2_obj))
KO_2_obj <- RunPCA(KO_2_obj, npcs = 50, verbose = FALSE)
ElbowPlot(KO_2_obj, ndims = 50)
KO_2_obj <- RunUMAP(KO_2_obj, reduction = "pca", dims = 1:30)
KO_2_obj <- FindNeighbors(KO_2_obj, dims = 1:30)

res_options <- c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0)

for (i in res_options){
  KO_2_obj <- FindClusters(KO_2_obj, res = i)
}

clustree(KO_2_obj, prefix = "RNA_snn_res.", node_size = 10,
         node_alpha = 0.8)

KO_2_obj <- FindClusters(KO_2_obj, res = 0.6)

KO_2_obj <- RenameIdents(KO_2_obj, "0" = "Cd3_T0", "1" = "Cd3_lo_T", "2" = "Cd3_T2", "3" = "early_T", "4" = "rearrange_TCR", 
                         "5" = "Rora_T","6" = "Cd34_TSP", "7" = "Interferon_response")

KO_2_obj <- AddMetaData(KO_2_obj, KO_2_obj@active.ident, "cell_type")

KO_2_obj <- subset(KO_2_obj,  cell_type != "Interferon_response" & cell_type != "Rora_T")

#WT object
Tcell.combined <- readRDS( "zhou_2022/cell_cycle_UMI_regress_seurat_integrated.rds")

Tcell.combined$sample_id <- as.factor(Tcell.combined$sample_id)

DefaultAssay(Tcell.combined) <- "RNA"

WT_obj <- subset(Tcell.combined, condition == "WT")

DimPlot(WT_obj, reduction = "umap", label = T, group.by = "cell_type")

WT_split <- SplitObject(WT_obj, "orig.ident")

table(WT_obj$cell_type)

WT_1_obj <- WT_split[[1]]

DefaultAssay(WT_1_obj) <- "RNA"


WT_1_obj <- NormalizeData(WT_1_obj)

s.genes <- c(str_to_title(tolower(cc.genes$s.genes)), "Cenpu")
g2m.genes <- c(str_to_title(tolower(cc.genes$g2m.genes)), "Pimreg", "Jpt1")

DefaultAssay(WT_1_obj) <- "RNA"

WT_1_obj <- CellCycleScoring(WT_1_obj, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)


WT_2_obj <- WT_split[[2]]

DefaultAssay(WT_2_obj) <- "RNA"


WT_2_obj <- NormalizeData(WT_2_obj)


WT_2_obj <- CellCycleScoring(WT_2_obj, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)


WT_1_obj <- FindVariableFeatures(WT_1_obj, selection.method = "vst", nfeatures = 3000)
WT_2_obj <- FindVariableFeatures(WT_2_obj, selection.method = "vst", nfeatures = 3000)


WT_1_obj <- ScaleData(WT_1_obj, verbose = FALSE ,vars.to.regress = c(genesV2$MGI.symbol), do.scale = T, do.center = T, features = row.names(WT_1_obj))
WT_1_obj <- RunPCA(WT_1_obj, npcs = 50, verbose = FALSE)
ElbowPlot(WT_1_obj, ndims = 50)
WT_1_obj <- RunUMAP(WT_1_obj, reduction = "pca", dims = 1:20)
WT_1_obj <- FindNeighbors(WT_1_obj, dims = 1:20)

res_options <- c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0)

for (i in res_options){
  WT_1_obj <- FindClusters(WT_1_obj, res = i)
}

clustree(WT_1_obj, prefix = "RNA_snn_res.", node_size = 10,
         node_alpha = 0.8)

WT_1_obj <- FindClusters(WT_1_obj, res = 0.7)

WT_1_obj <- RenameIdents(WT_1_obj, "0" = "rearrange_TCR", "1" = "Cd3_lo_T", "2" = "Cd3_T", "3" = "early_T", "4" = "Outlier", "5" = "Cd34_TSP")


WT_1_obj <- AddMetaData(WT_1_obj, WT_1_obj@active.ident, "cell_type")

WT_1_obj <- subset(WT_1_obj, cell_type != "Outlier" & cell_type != "Interferon_response" )


WT_2_obj <- ScaleData(WT_2_obj, verbose = FALSE,  vars.to.regress = c(genesV2$MGI.symbol), do.scale = T, do.center = T, features = row.names(WT_2_obj))
WT_2_obj <- RunPCA(WT_2_obj, npcs = 50, verbose = FALSE)
ElbowPlot(WT_2_obj, ndims = 50)
WT_2_obj <- RunUMAP(WT_2_obj, reduction = "pca", dims = 1:20)
WT_2_obj <- FindNeighbors(WT_2_obj, dims = 1:20)

res_options <- c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0)

for (i in res_options){
  WT_2_obj <- FindClusters(WT_2_obj, res = i)
}

clustree(WT_2_obj, prefix = "RNA_snn_res.", node_size = 10,
         node_alpha = 0.8)

WT_2_obj <- FindClusters(WT_2_obj, res = 0.6)
DimPlot(WT_2_obj, reduction = "umap", label = T)

WT_2_obj <- RenameIdents(WT_2_obj, "0" = "Cd3_lo_T", "1" = "rearrange_TCR", "2" = "Cd3_T", "3" = "early_T", 
                         "4" = "signaling_T", "5" = "Rora_T")

WT_2_obj <- AddMetaData(WT_2_obj, WT_2_obj@active.ident, "cell_type")

WT_2_obj <- subset(WT_2_obj, cell_type != "Rora_T")

#do feature space again
markers_WT2 <- FindAllMarkers(WT_2_obj, only.pos = T)

markers_WT2 %>%
  group_by(cluster) %>%
  top_n(n = 250, wt = avg_log2FC) -> markers_WT2

markers_WT2 <- markers_WT2$gene
markers_WT2 <- unique(markers_WT2)

markers_WT1 <- FindAllMarkers(WT_2_obj, only.pos = T)

markers_WT1 %>%
  group_by(cluster) %>%
  top_n(n = 250, wt = avg_log2FC) -> markers_WT1

markers_WT1 <- markers_WT1$gene
markers_WT1 <- unique(markers_WT1)

markers_KO2 <- FindAllMarkers(KO_2_obj, only.pos = T)

markers_KO2 %>%
  group_by(cluster) %>%
  top_n(n = 250, KO = avg_log2FC) -> markers_KO2

markers_KO2 <- markers_KO2$gene
markers_KO2 <- unique(markers_KO2)

markers_KO1 <- FindAllMarkers(KO_2_obj, only.pos = T)

markers_KO1 %>%
  group_by(cluster) %>%
  top_n(n = 250, KO = avg_log2FC) -> markers_KO1

markers_KO1 <- markers_KO1$gene
markers_KO1 <- unique(markers_KO1)

feature_space <- unique(c(markers_WT1, markers_WT2, markers_KO1, markers_KO2))

write.csv(feature_space, "TrAGEDy_results/Bcl11b_KO/feature_space.csv")


