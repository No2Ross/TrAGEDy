library(Seurat)
library(slingshot)
library(SingleCellExperiment)
library(cellAlign)

c_obj <- readRDS("TrAGEDy_results/simulated_data/not_shared/negative_control_1.rds")
d_obj <- readRDS("TrAGEDy_results/simulated_data/not_shared/negative_control_2.rds")

output <- list(c_obj, d_obj)

for (i in 1:length(output)){
  current <- output[[i]]
  
  current <- NormalizeData(current)
  
  current <- FindVariableFeatures(current, selection.method = "vst", nfeatures = 2000)
  
  all.genes <- rownames(current)
  current <- ScaleData(current, features = all.genes)
  
  current <- RunPCA(current, features = VariableFeatures(object = current))
  
  #print(ElbowPlot(current, ndims = 50))
  
  current <- FindNeighbors(current, dims = 1:15)
  current <- FindClusters(current, resolution = 0.5)
  
  #print(DimPlot(current, reduction = "pca"))
  current <- RunUMAP(current, dims = 1:15, n.components = 2)
  
  output[[i]] <- current
  
}

c_obj <- output[[1]]
d_obj <- output[[2]]


FeaturePlot(c_obj, reduction = "pca", features  = "sim_time")
print(DimPlot(c_obj, reduction = "pca"))
#d_obj <- subset(d_obj, sim_time < 100)
FeaturePlot(d_obj, reduction = "pca", features = "sim_time")
print(DimPlot(d_obj, reduction = "pca"))


print(DimPlot(c_obj, reduction = "umap"))
FeaturePlot(c_obj, reduction = "umap", features = "sim_time")

print(DimPlot(d_obj, reduction = "umap"))
FeaturePlot(d_obj, reduction = "umap", features = "sim_time")


c.markers <- FindAllMarkers(c_obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
c.markers_top15 <- c.markers %>%
  group_by(cluster) %>%
  slice_max(n = 100, order_by = avg_log2FC)

d.markers <- FindAllMarkers(d_obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
d.markers_top15 <- d.markers %>%
  group_by(cluster) %>%
  slice_max(n = 100, order_by = avg_log2FC)

features <- unique(d.markers_top15$gene, c.markers_top15$gene)

c_obj@assays$RNA@var.features <- features
d_obj@assays$RNA@var.features <- features

plot(c_sce$sim_time, c_sce@assays@data$counts["Target234",])
plot(d_sce$sim_time, d_sce@assays@data$counts["Target234",])


features <- c_obj@assays$RNA@var.features

#2: create trajectories seperately
c_sce <- as.SingleCellExperiment(c_obj, assay = "RNA")
start_cluster_1 <- levels(c_obj$seurat_clusters)[which(table(c_obj$seurat_clusters[which(c_obj$sim_time == 0)]) == max(table(c_obj$seurat_clusters[which(c_obj$sim_time == 0)])))]
c_sce <- slingshot(c_sce, reducedDim = 'UMAP', clusterLabels = c_sce@colData@listData[["RNA_snn_res.0.5"]], start.clus = start_cluster_1)
c_sling <- SlingshotDataSet(c_sce)

DimPlot(c_obj, reduction = "phate")
DimPlot(c_obj, reduction = "pca")
FeaturePlot(c_obj, reduction = "phate", features = "sim_time") +xlim(-0.03, 0.03) + ylim(-0.025, 0.025)
FeaturePlot(c_obj, reduction = "pca", features = "sim_time")
ggplot() + geom_point(aes(c_sce@int_colData@listData[["reducedDims"]]@listData[["PCA"]][,1], c_sce@int_colData@listData[["reducedDims"]]@listData[["PCA"]][,2], col = c_sce$slingPseudotime_1) )

plot3d(c_obj@reductions$phate@cell.embeddings[,1],c_obj@reductions$phate@cell.embeddings[,2], c_obj@reductions$phate@cell.embeddings[,3] )

require(scales)

# Create vector with levels of object@ident
identities <- levels(c_obj$RNA_snn_res.0.5)

# Create vector of default ggplot2 colors
my_color_palette <- hue_pal()(length(identities))

plot(reducedDims(c_sce)$PHATE, col = my_color_palette[c_sce$RNA_snn_res.0.5], pch = 16, cex = 0.5, bty='l')
lines(SlingshotDataSet(c_sce), col = "black", lwd = 2)

d_sce <- as.SingleCellExperiment(d_obj, assay = "RNA")
start_cluster_2 <- levels(d_obj$seurat_clusters)[which(table(d_obj$seurat_clusters[which(d_obj$sim_time == 0)]) == max(table(d_obj$seurat_clusters[which(d_obj$sim_time == 0)])))]

d_sce <- slingshot(d_sce, reducedDim = 'UMAP', clusterLabels = d_sce@colData@listData[["RNA_snn_res.0.5"]], start.clus = start_cluster_2)
d_sling <- SlingshotDataSet(d_sce)
DimPlot(d_obj, reduction = "phate")
DimPlot(d_obj, reduction = "umap")
DimPlot(d_obj, reduction = "pca")

FeaturePlot(d_obj, reduction = "phate", features = "sim_time") +xlim(-0.045, 0.045) + ylim(-0.025, 0.025)
ggplot() + geom_point(aes(d_sce@int_colData@listData[["reducedDims"]]@listData[["PCA"]][,1], d_sce@int_colData@listData[["reducedDims"]]@listData[["PCA"]][,2], col = d_sce$slingPseudotime_1) )
FeaturePlot(d_obj, reduction = "pca", features = "sim_time") 


require(scales)

# Create vector with levels of object@ident
identities <- levels(d_obj$RNA_snn_res.0.5)

# Create vector of default ggplot2 colors
my_color_palette <- hue_pal()(length(identities))

plot(reducedDims(d_sce)$PHATE, col = my_color_palette[d_sce$RNA_snn_res.0.5], pch = 16, cex = 0.5, bty='l')
lines(SlingshotDataSet(d_sce), col = "black", lwd = 2)

#c_sce$slingPseudotime_1 <- c_sce$sim_time
#d_sce$slingPseudotime_1 <- d_sce$sim_time

pseudo_end <- min(c(max(d_sce$slingPseudotime_1, c_sce$slingPseudotime_1, na.rm = T)))
window <- pseudo_end / 45

c_data <- as.matrix(c_sce@assays@data$counts)
d_data <- as.matrix(d_sce@assays@data$counts)

#c_data <- c_data[ which(row.names(c_data) %in% features),]
#d_data <- d_data[ which(row.names(d_data) %in% features),]

c_pseudo <- as.numeric(c_sce$slingPseudotime_1)
d_pseudo <- as.numeric(d_sce$slingPseudotime_1)

numPts = 50
interGlobalc = cellAlign::interWeights(expDataBatch = c_data, trajCond = c_pseudo,
                                       winSz = window, numPts = numPts)
interGlobald = cellAlign::interWeights(expDataBatch = d_data, trajCond = d_pseudo,
                                       winSz = window, numPts = numPts)

interScaledGlobalc = cellAlign::scaleInterpolate(interGlobalc)
interScaledGlobald = cellAlign::scaleInterpolate(interGlobald)

interScaledGlobalc$scaledData <- interScaledGlobalc$scaledData[ which( row.names(interScaledGlobalc$scaledData) %in% features), ]
interScaledGlobald$scaledData <- interScaledGlobald$scaledData[ which( row.names(interScaledGlobald$scaledData) %in% features), ]

interScaledGlobalc$scaledError <- interScaledGlobalc$scaledError[ which( row.names(interScaledGlobalc$scaledError) %in% features), ]
interScaledGlobald$scaledError <- interScaledGlobald$scaledError[ which( row.names(interScaledGlobald$scaledError) %in% features), ]



alignment = globalAlign(interScaledGlobalc$scaledData, interScaledGlobald$scaledData,
                        scores = list(query = interScaledGlobalc$traj, 
                                      ref = interScaledGlobald$traj),
                        sigCalc = F, numPerm = 20, dist.method = "Euclidean", normDist = F)
plotAlign(alignment)

alignment = globalAlign(interGlobalc$interpolatedVals[ which( row.names(interGlobalc$interpolatedVals) %in% features), ], interGlobald$interpolatedVals[ which( row.names(interGlobald$interpolatedVals) %in% features), ],
                        scores = list(query = interGlobalc$traj, 
                                      ref = interGlobald$traj),
                        sigCalc = F, numPerm = 20, dist.method = "correlation")
plotAlign(alignment)

