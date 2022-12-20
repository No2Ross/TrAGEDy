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

int_obj <- readRDS("path/to/Tcell/dataset/WZ_Bcl11bKO_tworuns_integrated.rds")

genesV2 <- read.csv("path/to/seurat/cellcycle/genes/human_to_mouse_seurat_cellcycle.csv")

DefaultAssay(int_obj) <- "RNA"

WT_ff_meta <- data.frame(row.names = row.names(int_obj@meta.data), old = int_obj$sample_id, new = rep("WT", length(colnames(int_obj))))

WT_ff_meta$new[which(str_detect(WT_ff_meta$old, "FF") == T)] <- "Bcl11b_KO"

int_obj <- AddMetaData(int_obj, WT_ff_meta$new, "condition")

day_meta <- data.frame(row.names = row.names(int_obj@meta.data), old = int_obj$sample_id, new = rep("D10", length(colnames(int_obj))))

day_meta$new[which(str_detect(day_meta$old, "13") == T)] <- "D13"

int_obj <- AddMetaData(int_obj, day_meta$new, "day")

DefaultAssay(int_obj) <- "RNA"

#remove no cre control
int_obj <- subset(int_obj, sample_id != "D10_FF_NoCre_rep1")

FeaturePlot(int_obj, reduction = "umap", features = c("Ptcra", "Rag1"), split.by = "condition")
FeaturePlot(int_obj, reduction = "umap", features = "Trac", split.by = "condition")
FeaturePlot(int_obj, reduction = "umap", features = "Rora", split.by = "condition")

VlnPlot(int_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

#int_obj <- subset(int_obj, nCount_RNA < 20000)

WT_obj <- subset(int_obj, condition == "WT")
KO_obj <- subset(int_obj, condition == "Bcl11b_KO")

FeaturePlot(WT_obj, feature = "Zap70", split.by = "day")

WT.list <- SplitObject(WT_obj, "orig.ident")

out.list <- list()

for (i in 1:length(WT.list)){
  current <- WT.list[[i]]
  
  s.genes <- c(str_to_title(tolower(cc.genes$s.genes)), "Cenpu")
  g2m.genes <- c(str_to_title(tolower(cc.genes$g2m.genes)), "Pimreg", "Jpt1")
  
  DefaultAssay(current) <- "RNA"
  
  current <- CellCycleScoring(current, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
  
  
  current.list <- SplitObject(current, "sample_id")
  
  current.list <- lapply(X = current.list, FUN = function(x) {
    x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 3000)
    x <- NormalizeData(x)
  })
  
  features <- SelectIntegrationFeatures(object.list = current.list, nfeatures = 3000)
  
  WT.anchors <- FindIntegrationAnchors(object.list = current.list, anchor.features = features)
  
  # this command creates an 'integrated' data assay
  WT.combined <- IntegrateData(anchorset = WT.anchors)
  
  out.list[[i]] <- WT.combined
  
  
}

features <- SelectIntegrationFeatures(object.list = out.list, nfeatures = 3000)


WT.anchors <- FindIntegrationAnchors(object.list = out.list, anchor.features = features)

# this command creates an 'integrated' data assay
WT.combined <- IntegrateData(anchorset = WT.anchors)


DefaultAssay(WT.combined) <- "integrated"

KO.list <- SplitObject(KO_obj, "orig.ident")

out.list <- list()

for (i in 1:length(KO.list)){
  current <- KO.list[[i]]
  
  s.genes <- c(str_to_title(tolower(cc.genes$s.genes)), "Cenpu")
  g2m.genes <- c(str_to_title(tolower(cc.genes$g2m.genes)), "Pimreg", "Jpt1")
  
  DefaultAssay(current) <- "RNA"
  
  current <- CellCycleScoring(current, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
  
  
  current.list <- SplitObject(current, "sample_id")
  
  current.list <- lapply(X = current.list, FUN = function(x) {
    x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 3000)
    x <- NormalizeData(x)
  })
  
  features <- SelectIntegrationFeatures(object.list = current.list, nfeatures = 3000)
  
  
  KO.anchors <- FindIntegrationAnchors(object.list = current.list, anchor.features = features)
  
  rm(current)
  rm(current.list)
  gc()
  
  # this command creates an 'integrated' data assay
  KO.combined <- IntegrateData(anchorset = KO.anchors)
  
  out.list[[i]] <- KO.combined
  
  
}

features <- SelectIntegrationFeatures(object.list = out.list, nfeatures = 3000)


KO.anchors <- FindIntegrationAnchors(object.list = out.list, anchor.features = features)

KO.combined <- IntegrateData(anchorset = KO.anchors)

DefaultAssay(KO.combined) <- "integrated"

Tcell.list <- list(KO.combined, WT.combined)
names(Tcell.list) <- c("Bcl11b_KO", "WT")

features <- SelectIntegrationFeatures(object.list = Tcell.list, nfeatures = 3000)

Tcell.anchors <- FindIntegrationAnchors(object.list = Tcell.list, anchor.features = features)

Tcell.combined <- IntegrateData(anchorset = Tcell.anchors)


DefaultAssay(Tcell.combined) <- "integrated"

Tcell.combined <- ScaleData(Tcell.combined, verbose = FALSE, vars.to.regress = c(genesV2$MGI.symbol))
Tcell.combined <- RunPCA(Tcell.combined, npcs = 50, verbose = FALSE)
ElbowPlot(Tcell.combined, ndims = 50)
Tcell.combined <- RunUMAP(Tcell.combined, reduction = "pca", dims = 1:25, n.components = 2)
Tcell.combined <- FindNeighbors(Tcell.combined, dims = 1:25)

res_options <- c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0)

for (i in res_options){
  Tcell.combined <- FindClusters(Tcell.combined, res = i)
}

clustree(Tcell.combined, prefix = "integrated_snn_res.", node_size = 10,
         node_alpha = 0.8)

Tcell.combined <- FindClusters(Tcell.combined, res = 0.6)

DefaultAssay(Tcell.combined) <- "RNA"

Tcell.combined <- RenameIdents(Tcell.combined, "0" = "Cd3_T", "1" = "Cd3_lo_T", "2" = "early_T", "3" = "rearrange_TCR_KO", "4" = "rearrange_TCR_WT", "5" = "early_Rora_T","6" = "Cd34_TSP", "7" = "signalling_T",
                               "8" = "Rora_T", "9" = "Outlier", "10" = "Interferon_response")

Tcell.combined <- AddMetaData(Tcell.combined, Tcell.combined@active.ident, "cell_type_seurat")

DimPlot(Tcell.combined, reduction = "umap", label = T, split.by = "condition", group.by = "cell_type_seurat")

Tcell.combined <- subset(Tcell.combined, cell_type_seurat != "Interferon_response" & cell_type_seurat != "Rora_T" &
                           cell_type_seurat != "early_Rora_T" & cell_type_seurat != "Outlier" & cell_type_seurat != "signalling_T")


DefaultAssay(Tcell.combined) <- "RNA"


#Get pseudotime information for TradeSeq analysis
phate_output <- as.matrix(phate(t(Tcell.combined@assays$integrated@scale.data), ndim=10, seed=1))
phate.reduction <- CreateDimReducObject(embeddings = phate_output,
                                        key="PHATE_",
                                        assay="RNA")

Tcell.combined@reductions$phate <- phate.reduction

DimPlot(Tcell.combined, reduction = "phate", label = T, split.by = "condition", group.by = "cell_type_seurat")

combined_sce <- as.SingleCellExperiment(Tcell.combined, assay = "RNA")
combined_sce <- slingshot(combined_sce, reducedDim = 'PHATE', clusterLabels = combined_sce@colData@listData[["cell_type_seurat"]], start.clus = "Cd34_TSP")
combined_sling <- SlingshotDataSet(combined_sce)

combined_sce$pseudotime <- combined_sce$slingPseudotime_1

combined_sce <- subset(combined_sce, , is.na(slingPseudotime_1) == F)


