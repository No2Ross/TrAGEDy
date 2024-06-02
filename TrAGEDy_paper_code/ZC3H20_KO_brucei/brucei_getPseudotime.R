library(Seurat)
library(ggplot2)
library(RColorBrewer)
library(SingleCellExperiment)
library(slingshot)
library(tradeSeq)


tryp_integrated <- readRDS("/Users/rosslaidlaw/R/TrAGEDy_V2/Tbrucei_Zc3h20_KO/objects/WT_ZC3H20_with_phate.rds")

DimPlot(tryp_integrated, reduction = "phate")

tryp_integrated <- AddMetaData(tryp_integrated, tryp_integrated@active.ident, "cell_type")

#subset out stumpy because the Trajectory runs from LS.A1 to LS B.2 and doesn't go to stumpy
tryp_integrated <- subset(tryp_integrated, cell_type != "SS A" & cell_type != "SS B")

DimPlot(tryp_integrated, reduction = "phate", group.by = "cell_type")

tryp_integrated <- AddMetaData(tryp_integrated, tryp_integrated@active.ident, "cell_type")

tryp_sce <- as.SingleCellExperiment(tryp_integrated, assay = "RNA")
tryp_sce <- slingshot(tryp_sce, reducedDim = 'PHATE', clusterLabels = tryp_sce@colData@listData[["cell_type"]], start.clus = "LS A.1")
tryp_sling <- SlingshotDataSet(tryp_sce)

saveRDS(tryp_sce, "/Users/rosslaidlaw/R/TrAGEDy_V2/Tbrucei_Zc3h20_KO/objects/WT_ZC3H20_with_phate_sce.rds")

