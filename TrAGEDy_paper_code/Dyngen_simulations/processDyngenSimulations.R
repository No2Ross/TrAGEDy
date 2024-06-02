library(Seurat)
library(SingleCellExperiment)
library(slingshot)
library(zellkonverter)
library(stringr)
library(dplyr)



processDatasets <- function(dataset, dim){
  dataset <- NormalizeData(dataset)
  
  dataset@assays$RNA@var.features <- row.names(dataset)
  
  all.genes <- rownames(dataset)
  dataset <- ScaleData(dataset, features = all.genes)
  
  dataset <- RunPCA(dataset, features = VariableFeatures(object = dataset))
  
  dataset <- FindNeighbors(dataset, dims = 1:dim)
  dataset <- FindClusters(dataset, resolution = 0.2)
  
  dataset <- RunUMAP(dataset, dims = 1:dim, n.components = 2)
  
  return(dataset)
}

path <- "/Users/rosslaidlaw/R/TrAGEDy_V2/simulated_datasets/objects/"

datasets <- c("diverge_converge_c",
              "diverge_converge_d",
              "negative_control_1",
              "negative_control_2",
              "start_shared_wt",
              "start_shared_ko")

counter <- 0
feature_space <- c()
dim <- 5

for (i in datasets){
  
  current <- readRDS(paste0(path, i, ".rds"))
  
  current@assays$spliced <- NULL
  current@assays$unspliced <- NULL
  current@assays$protein <- NULL
  
  state_dataframe <- as.data.frame(current@misc$traj_progressions)
  state_dataframe$id <- state_dataframe$cell_id
  
  state_dataframe[which(state_dataframe$percentage >= 0 & state_dataframe$percentage <0.5),"id"] <- state_dataframe[which(state_dataframe$percentage >= 0 & state_dataframe$percentage <0.5),"from"]
  state_dataframe[which(state_dataframe$percentage > 0.51 & state_dataframe$percentage <=1),"id"] <- state_dataframe[which(state_dataframe$percentage > 0.51 & state_dataframe$percentage <=1),"to"]
  
  current$state <- state_dataframe$id
  
  current <- processDatasets(current, dim)
  
  print(DimPlot(current, label = T))
  
  markers <- FindAllMarkers(current, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.1)
  markers_top100 <- markers %>%
    group_by(cluster) %>%
    slice_max(n = 100, order_by = avg_log2FC)
  
  current_sce <- as.SingleCellExperiment(current)
  
  start_cluster <- levels(current$seurat_clusters)[which(table(current$seurat_clusters[which(current$sim_time == 0)]) == max(table(current$seurat_clusters[which(current$sim_time == 0)])))]
  
  current_sce <- slingshot(current_sce, reducedDim = 'UMAP', clusterLabels = current_sce@colData@listData[["RNA_snn_res.0.2"]], start.clus = start_cluster)
  
  current$slingPseudotime <- current_sce$slingPseudotime_1
  
  current_sce_save <- as.SingleCellExperiment(current)
  
  feature_space <- append(feature_space, markers_top100$gene)
  
  counter <- counter + 1
  
  #When we've processed a pair of datasets, save the feature space
  if(counter == 2){
    
    feature_space <- unique(feature_space)
    
    file_name <- paste0( str_split(i, "_")[[1]][1] , "_", str_split(i, "_")[[1]][2], "_feature_space" )
    
    write.csv(feature_space, paste0(path, file_name,".csv" ))
    
    feature_space <- c()
    counter <- 0
    
  }
  
  saveRDS(current, file = paste0(path, i, "_processed.rds"))
  writeH5AD(current_sce_save, file = paste0(path, i, "_processed.h5ad"))
  
}

