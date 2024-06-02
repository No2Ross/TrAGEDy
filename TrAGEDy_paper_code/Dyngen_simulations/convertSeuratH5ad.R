library(Seurat)
library(zellkonverter)

path <- "/Users/rosslaidlaw/R/TrAGEDy_V2/simulated_datasets/objects/"

datasets <- c("diverge_converge_c",
              "diverge_converge_d",
              "negative_control_1",
              "negative_control_2",
              "start_shared_wt",
              "start_shared_ko")

for (i in datasets){
  
  current <- readRDS(paste0(path, i, ".rds"))
  
  current@assays$spliced <- NULL
  current@assays$unspliced <- NULL
  current@assays$protein <- NULL
  
  current_sce <- as.SingleCellExperiment(current)
  
  writeH5AD(current_sce, file = paste0(path, i, ".h5ad"))
  
}





