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

genesV2 <- read.csv("/Users/rosslaidlaw/R/TrAGEDy_V2/general/human_to_mouse_seurat_cellcycle.csv")
cell_cycle <- readxl::read_xlsx("/Users/rosslaidlaw/R/TrAGEDy_V2/general/mgi_cellCycle_GO.xlsx")

int_obj <- readRDS("/Users/rosslaidlaw/R/TrAGEDy_V2/Tcell_Bcl11b_KO/objects/WZ_Bcl11bKO_tworuns_integrated.Rds")

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

DimPlot(int_obj, group.by = "seurat_clusters")

FeaturePlot(int_obj, reduction = "umap", features = c("Ptcra", "Rag1"), split.by = "condition")
FeaturePlot(int_obj, reduction = "umap", features = "Trac", split.by = "condition")
FeaturePlot(int_obj, reduction = "umap", features = "Zap70", split.by = "condition")
FeaturePlot(int_obj, reduction = "umap", features = c("Cd8a", "Zap70"), split.by = "condition")

VlnPlot(int_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)


int_obj[["RNA"]] <- split(int_obj[["RNA"]], f = int_obj$orig.ident)

int_obj <- NormalizeData(int_obj)
int_obj <- FindVariableFeatures(int_obj)
int_obj <- ScaleData(int_obj, vars.to.regress = c("nCount_RNA", genesV2$MGI.symbol))
int_obj <- RunPCA(int_obj)

int_obj <- IntegrateLayers(object = int_obj, method = CCAIntegration, orig.reduction = "pca", new.reduction = "integrated.cca",
                        verbose = FALSE)

DefaultAssay(int_obj) <- "integrated"

int_obj <- ScaleData(int_obj, vars.to.regress = c("nCount_RNA", genesV2$MGI.symbol))

DefaultAssay(int_obj) <- "RNA"

int_obj[["RNA"]] <- JoinLayers(int_obj[["RNA"]])
ElbowPlot(int_obj, ndims = 50)

int_obj <- FindNeighbors(int_obj, reduction = "integrated.cca", dims = 1:25)
int_obj <- FindClusters(int_obj, resolution = 0.6)
int_obj <- RunUMAP(int_obj, reduction = "integrated.cca",dims = 1:25)


DimPlot(int_obj, label = T)
DimPlot(int_obj, reduction = "integrated.cca")


DefaultAssay(int_obj) <- "RNA"

combined_markers <- FindAllMarkers(int_obj, logfc.threshold = 0.5, only.pos = T)
combined_markers <- subset(combined_markers, p_val_adj < 0.05)

FeaturePlot(int_obj, features = c("Cd34", "Zap70", "Cd3e", "Trac", "Ptcra", "Rora"))
VlnPlot(int_obj, features = c("Cd34", "Rag1", "Cd3e", "Trac", "Ptcra", "Rora"))

VlnPlot(int_obj, features = c("Tcrg-C2", "Tcrg-C4", "Dntt"))


VlnPlot(int_obj, features = c("mt-Nd1", "mt-Nd3"))


int_obj <- RenameIdents(int_obj, "0" = "Cd3_hi_T", "1" = "Cd3_int_T_1", "2" = "rearrange_TCR_T", "3" = "Cd34_lo_Cd3_lo_T", "4" = "Cd3_int_T_2",
                               "5" = "Cd3_lo_T", "6" = "TCR_AB_T",  "7" = "Rora_T", "8" = "Cd34_TSP",
                        "9" = "Interferon_response")

int_obj$cell_type_seurat <- int_obj@active.ident

DimPlot(int_obj, label = T)

int_obj <- subset(int_obj, cell_type_seurat != "Interferon_response")

#Get pseudotime information for TradeSeq analysis
phate_output <- as.matrix(phate(t(int_obj@assays$integrated@scale.data), ndim=15, seed=1))
phate.reduction <- CreateDimReducObject(embeddings = phate_output,
                                        key="PHATE_",
                                        assay="RNA")

int_obj@reductions$phate <- phate.reduction

DimPlot(int_obj, reduction = "phate", label = T, split.by = "condition", group.by = "cell_type_seurat")
FeaturePlot(int_obj, reduction="phate", feature = "Mki67") + xlim(-0.04, 0.04) + ylim(-0.03, 0.04)


combined_sce <- as.SingleCellExperiment(int_obj, assay = "RNA")
combined_sce <- slingshot(combined_sce, reducedDim = 'PHATE', clusterLabels = combined_sce@colData@listData[["cell_type_seurat"]], start.clus = "Cd34_TSP")
combined_sling <- SlingshotDataSet(combined_sce)

combined_sce$pseudotime <- combined_sce$slingPseudotime_1
int_obj$pseudotime <- combined_sce$pseudotime 

combined_sce <- subset(combined_sce, ,is.na(slingPseudotime_1) == F)
combined_sce <- subset(combined_sce, ,cell_type_seurat != "Cd3_lo_T" & cell_type_seurat != "Cd3_int_T_2")

kept_cells <- colnames(combined_sce)

keep_cells <- data.frame("keep_cells" = rep("no", dim(int_obj)[2]),
                         row.names = colnames(int_obj))
keep_cells[kept_cells,] <- "yes"

int_obj$keep_cells <- keep_cells$keep_cells

int_obj <- subset(int_obj, keep_cells == "yes")

saveRDS(combined_sce, "/Users/rosslaidlaw/R/TrAGEDy_V2/Tcell_Bcl11b_KO/objects/phate_cell_cycle_UMI_regress_sce_integrated.rds")
saveRDS(int_obj, "/Users/rosslaidlaw/R/TrAGEDy_V2/Tcell_Bcl11b_KO/objects/phate_cell_cycle_UMI_regress_seurat_integrated.rds")





