# library(reticulate)
# use_python("/usr/local/bin/python3.8")
# reticulate::import("phate")
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
#library(STACAS)
#library(monocle)
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

Tcell.combined <- readRDS("/Users/rosslaidlaw/R/TrAGEDy_V2/Tcell_Bcl11b_KO/objects/WZ_Bcl11bKO_tworuns_integrated.Rds")

DefaultAssay(Tcell.combined) <- "RNA"

WT_ff_meta <- data.frame(row.names = row.names(Tcell.combined@meta.data), old = Tcell.combined$sample_id, new = rep("WT", length(colnames(Tcell.combined))))

WT_ff_meta$new[which(str_detect(WT_ff_meta$old, "FF") == T)] <- "Bcl11b_KO"

Tcell.combined <- AddMetaData(Tcell.combined, WT_ff_meta$new, "condition")

day_meta <- data.frame(row.names = row.names(Tcell.combined@meta.data), old = Tcell.combined$sample_id, new = rep("D10", length(colnames(Tcell.combined))))

day_meta$new[which(str_detect(day_meta$old, "13") == T)] <- "D13"

Tcell.combined <- AddMetaData(Tcell.combined, day_meta$new, "day")

DefaultAssay(Tcell.combined) <- "RNA"

#remove no cre control
Tcell.combined <- subset(Tcell.combined, sample_id != "D10_FF_NoCre_rep1")

Tcell.combined$sample_id <- as.factor(Tcell.combined$sample_id)

Tcell.combined[["percent.mt"]] <- PercentageFeatureSet(Tcell.combined, pattern = "^mt-")

DefaultAssay(Tcell.combined) <- "RNA"

KO_obj <- subset(Tcell.combined, condition == "Bcl11b_KO")

WT_obj <- subset(Tcell.combined, condition == "WT")

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

KO_1_obj <- FindVariableFeatures(KO_1_obj, selection.method = "vst", nfeatures = 2000)
KO_2_obj <- FindVariableFeatures(KO_2_obj, selection.method = "vst", nfeatures = 2000)

WT_1_obj <- FindVariableFeatures(WT_1_obj, selection.method = "vst", nfeatures = 2000)
WT_2_obj <- FindVariableFeatures(WT_2_obj, selection.method = "vst", nfeatures = 2000)

var_features <- SelectIntegrationFeatures(list(WT_1_obj, WT_2_obj, KO_1_obj, KO_2_obj), nfeatures = 2000)

s.genes <- c(str_to_title(tolower(cc.genes$s.genes)), "Cenpu")
g2m.genes <- c(str_to_title(tolower(cc.genes$g2m.genes)), "Pimreg", "Jpt1")

DefaultAssay(KO_1_obj) <- "RNA"

KO_1_obj <- CellCycleScoring(KO_1_obj, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

DefaultAssay(KO_2_obj) <- "RNA"


KO_2_obj <- NormalizeData(KO_2_obj)

KO_2_obj <- CellCycleScoring(KO_2_obj, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

KO_1_obj <- ScaleData(KO_1_obj, verbose = FALSE , vars.to.regress = c(genesV2$MGI.symbol),do.scale = T, do.center = T, features = var_features)
KO_1_obj <- RunPCA(KO_1_obj, npcs = 50, verbose = FALSE, features = var_features)
ElbowPlot(KO_1_obj, ndims = 50)
KO_1_obj <- RunUMAP(KO_1_obj, reduction = "pca", dims = 1:25)
KO_1_obj <- FindNeighbors(KO_1_obj, dims = 1:25)

KO_1_obj <- FindClusters(KO_1_obj, res = 0.6)
DimPlot(KO_1_obj, label = T)

FeaturePlot(KO_1_obj, features = c("Cd34", "Zap70", "Cd3e", "Cd3d", "Ptcra", "Rora"))
FeaturePlot(KO_1_obj, features = c("Rora"), split.by = "seurat_clusters")
VlnPlot(KO_1_obj, features = c("Cd34", "Rag1", "Cd3e", "Cd3d", "Ptcra", "Rora"))
FeaturePlot(KO_1_obj, features = c("Id2"))

#KO genes
VlnPlot(KO_1_obj, features = c("Il2rb", "Sox5", "Id2","Rora"))

DimPlot(KO_1_obj, label = T)
DimPlot(KO_1_obj, label = T, reduction = "pca")

KO_1_markers <- FindAllMarkers(KO_1_obj, logfc.threshold = 0.5, only.pos = T)

KO_1_obj <- RenameIdents(KO_1_obj, "0" = "rearrange_TCR_T", "1" = "Cd3_lo_T", "2" = "Rora_T", "3" = "Cd3_T", "4" = "Cd34_TSP")

KO_1_obj <- AddMetaData(KO_1_obj, KO_1_obj@active.ident, "cell_type")

KO_1_obj <- subset(KO_1_obj, cell_type != "Rora_T")

saveRDS(KO_1_obj, "/Users/rosslaidlaw/R/TrAGEDy_V2/Tcell_Bcl11b_KO/objects/KO_1_obj.rds")

KO_2_obj <- ScaleData(KO_2_obj, verbose = FALSE , vars.to.regress = c(genesV2$MGI.symbol),do.scale = T, do.center = T, features = var_features)
KO_2_obj <- RunPCA(KO_2_obj, npcs = 50, verbose = FALSE, features = var_features)

ElbowPlot(KO_2_obj, ndims = 50)
KO_2_obj <- RunUMAP(KO_2_obj, reduction = "pca", dims = 1:30)
KO_2_obj <- FindNeighbors(KO_2_obj, dims = 1:30)


KO_2_obj <- FindClusters(KO_2_obj, res = 0.5)

DimPlot(KO_2_obj, label = T)
DimPlot(KO_2_obj, reduction = "pca", label = T)
DimPlot(KO_2_obj, reduction = "pca", label = T, dims = c(2,3))

FeaturePlot(KO_2_obj, features = c("Cd34", "Zap70", "Cd3e", "Cd3d", "Ptcra", "Rora"))
VlnPlot(KO_2_obj, features = c("Cd34", "Zap70", "Cd3e", "Rag1", "Ptcra", "Rora"))

#KO genes
VlnPlot(KO_2_obj, features = c("Il2rb", "Sox5", "Id2","Rora"))

KO_2_markers <- FindAllMarkers(KO_2_obj, logfc.threshold = 0.5, only.pos = T)


#Early T: Doesn't express CD3, does express Cd34, doesn't express Cd74
#Cd34_TSP: Doesn't express CD3, does express Cd34, does express Cd74
KO_2_obj <- RenameIdents(KO_2_obj, "0" = "Cd3_hi_T", "1" = "Cd3_int_T", "2" = "Cd3_lo_T", "3" = "rearrange_TCR_T", "4" = "Rora_T", 
                         "5" = "Cd34_TSP", "6" = "Interferon_response")


KO_2_obj <- AddMetaData(KO_2_obj, KO_2_obj@active.ident, "cell_type")

KO_2_obj <- subset(KO_2_obj, cell_type != "Interferon_response" & cell_type != "Rora_T")

saveRDS(KO_2_obj, "/Users/rosslaidlaw/R/TrAGEDy_V2/Tcell_Bcl11b_KO/objects/KO_2_obj.rds")


#WT analysis

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


WT_1_obj <- ScaleData(WT_1_obj, verbose = FALSE , vars.to.regress = c(genesV2$MGI.symbol),do.scale = T, do.center = T, features = var_features)
WT_1_obj <- RunPCA(WT_1_obj, npcs = 50, verbose = FALSE, features = var_features)
ElbowPlot(WT_1_obj, ndims = 50)
WT_1_obj <- RunUMAP(WT_1_obj, reduction = "pca", dims = 1:22)
WT_1_obj <- FindNeighbors(WT_1_obj, dims = 1:22)

WT_1_obj <- FindClusters(WT_1_obj, res = 0.6)

DimPlot(WT_1_obj, label = T)

WT_1_markers <- FindAllMarkers(WT_1_obj, logfc.threshold = 0.5, only.pos = T)

FeaturePlot(WT_1_obj, features = c("Cd34", "Zap70", "Cd3e", "Cd3d", "Ptcra", "Rora"))
VlnPlot(WT_1_obj, features = c("Cd34", "Rag1", "Cd3e", "Cd3d", "Ptcra", "Rora"))

WT_1_obj <- RenameIdents(WT_1_obj, "0" = "rearrange_TCR_T", "1" = "Cd3_int_T", "2" = "Cd3_hi_T", "3" = "Cd3_lo_T", 
                         "4" = "Rora_T", "5" = "Cd3_T_MT_high", "6" = "Cd34_TSP")

WT_1_obj <- AddMetaData(WT_1_obj, WT_1_obj@active.ident, "cell_type")

WT_1_obj <- subset(WT_1_obj,  cell_type != "Cd3_T_MT_high" & cell_type != "Rora_T" )

saveRDS(WT_1_obj, "/Users/rosslaidlaw/R/TrAGEDy_V2/Tcell_Bcl11b_KO/objects/WT_1_obj.rds")



WT_2_obj <- ScaleData(WT_2_obj, verbose = FALSE , vars.to.regress = c(genesV2$MGI.symbol),do.scale = T, do.center = T, features = var_features)
WT_2_obj <- RunPCA(WT_2_obj, npcs = 50, verbose = FALSE, features = var_features)

ElbowPlot(WT_2_obj, ndims = 50)
WT_2_obj <- RunUMAP(WT_2_obj, reduction = "pca", dims = 1:20)
WT_2_obj <- FindNeighbors(WT_2_obj, dims = 1:20)


WT_2_obj <- FindClusters(WT_2_obj, res = 0.6)
DimPlot(WT_2_obj, reduction = "umap", label = T)

FeaturePlot(WT_2_obj, features = c("Cd34", "Zap70", "Cd3e", "Cd3d", "Ptcra", "Rora", "Trac"))
VlnPlot(WT_2_obj, features = c("Cd34", "Zap70", "Cd3e", "Cd3d", "Ptcra", "Rag1"))
VlnPlot(WT_2_obj, features = c("Cd4", "Cd8a", "Rora", "Trac"))

WT_2_obj <- RenameIdents(WT_2_obj, "0" = "Cd3_hi_T", "1" = "Cd3_int_T", "2" = "rearrange_TCR_T", "3" = "Cd3_lo_T", 
                         "4" = "TCRAB_T", "5" = "Rora_T")

WT_2_obj <- AddMetaData(WT_2_obj, WT_2_obj@active.ident, "cell_type")

WT_2_obj <- subset(WT_2_obj, cell_type != "Rora_T")

saveRDS(WT_2_obj, "/Users/rosslaidlaw/R/TrAGEDy_V2/Tcell_Bcl11b_KO/objects/WT_2_obj.rds")

DimPlot(WT_1_obj, label = T)
DimPlot(WT_2_obj, label = T)
DimPlot(KO_1_obj, label = T)
DimPlot(KO_2_obj, label = T)

write.csv(var_features,"/Users/rosslaidlaw/R/TrAGEDy_V2/Tcell_Bcl11b_KO/selectIntegrationFeatures_TrAGEDy.csv")

#do feature space again
markers_WT2 <- FindAllMarkers(WT_2_obj, only.pos = T, logfc.threshold = 0.5)

markers_WT2 <- subset(markers_WT2, p_val_adj < 0.05)
markers_WT2 <- markers_WT2$gene
markers_WT2 <- unique(markers_WT2)

markers_WT1 <- FindAllMarkers(WT_1_obj, only.pos = T, logfc.threshold = 0.5)

markers_WT1 <- subset(markers_WT1, p_val_adj < 0.05)
markers_WT1 <- markers_WT1$gene
markers_WT1 <- unique(markers_WT1)

markers_KO2 <- FindAllMarkers(KO_2_obj, only.pos = T, logfc.threshold = 0.5)

markers_KO2 <- subset(markers_KO2, p_val_adj < 0.05)
markers_KO2 <- markers_KO2$gene
markers_KO2 <- unique(markers_KO2)

markers_KO1 <- FindAllMarkers(KO_1_obj, only.pos = T, logfc.threshold = 0.5)

markers_KO1 <- subset(markers_KO1, p_val_adj < 0.05)
markers_KO1 <- markers_KO1$gene
markers_KO1 <- unique(markers_KO1)

feature_space <- unique(c(markers_WT1, markers_WT2, markers_KO1, markers_KO2))

write.csv(feature_space, "/Users/rosslaidlaw/R/TrAGEDy_V2/Tcell_Bcl11b_KO/TrAGEDy_feature_space_with_rearrangeTCR.csv")

