library(Seurat)
library(ggplot2)
library(RColorBrewer)
library(SingleCellExperiment)
library(slingshot)
library(tradeSeq)
library(ggborderline)

tryp_integrated <- readRDS("path/to/integrated/object/WT_ZC3H20_with_phate.rds")

DimPlot(tryp_integrated, reduction = "phate")

tryp_integrated <- AddMetaData(tryp_integrated, tryp_integrated@active.ident, "type")

#subset out stumpy
tryp_integrated <- subset(tryp_integrated, type != "SS A" & type != "SS B")

DimPlot(tryp_integrated, reduction = "phate", group.by = "cell_type")

tryp_integrated <- AddMetaData(tryp_integrated, tryp_integrated@active.ident, "cell_type")

tryp_sce <- as.SingleCellExperiment(tryp_integrated, assay = "RNA")
tryp_sce <- slingshot(tryp_sce, reducedDim = 'PHATE', clusterLabels = tryp_sce@colData@listData[["cell_type"]], start.clus = "LS A.1")
tryp_sling <- SlingshotDataSet(tryp_sce)

#tradeseq
set.seed(3)
icMat <- evaluateK(counts = as.matrix(assays(tryp_sce)$counts),
                   sds = tryp_sling,
                   nGenes = 300,
                   k = 4:10)

pseudotime <- slingPseudotime(tryp_sling, na = FALSE)
cellWeights <- slingCurveWeights(tryp_sling)

sce <- fitGAM(tryp_sce, conditions = factor(colData(tryp_sce)$line),
              nknots = 9, verbose = FALSE)

WT_sce <- subset(tryp_sce, , line == "WT")
KO_sce <- subset(tryp_sce, , line == "ZC3H20_KO")

#knots are based on the quantiles
quantile_limits <- quantile(sce$slingPseudotime_1, probs = seq(0,1,0.125))

#PCT

invalid_genes_list <- list()

windows <- list(c(1,2,3), c(3,4,5), c(5,6,7), c(7,8,9))

for(i in 1:length(windows)){
  knot_min <- min(windows[[i]])
  knot_max <- max(windows[[i]])
  
  invalid_genes <- c()
  
  current_WT <- subset(WT_sce, ,slingPseudotime_1 > quantile_limits[knot_min] & 
                         slingPseudotime_1 < quantile_limits[knot_max])
  
  current_KO <- subset(KO_sce, ,slingPseudotime_1 > quantile_limits[knot_min] & 
                         slingPseudotime_1 < quantile_limits[knot_max])
  
  count_mtx_WT <- 1 - (rowCounts(current_WT@assays@data@listData$counts, value = 0) /length(colnames(current_WT@assays@data@listData$counts)))
  names(count_mtx_WT) <- row.names(current_WT)
  
  count_mtx_KO <- 1 - (rowCounts(current_KO@assays@data@listData$counts, value = 0) /length(colnames(current_KO@assays@data@listData$counts)))
  names(count_mtx_KO) <- row.names(current_KO)
  
  invalid_WT <- count_mtx_WT[which(count_mtx_WT < 0.1)]
  invalid_KO <- count_mtx_KO[which(count_mtx_KO < 0.1)]
  
  invalid_genes <- intersect(names(invalid_KO), names(invalid_WT))
  
  invalid_genes_list[[i]] <- invalid_genes
}


cond_test_results <- list()

windows <- list(c(1,2,3), c(3,4,5), c(5,6,7), c(7,8,9))

for (i in 1:4){
  cond_current_results <- conditionTest(sce, knots = c(min(windows[[i]]), max(windows[[i]])), l2fc = 0.5)
  cond_test_results[[i]] <- cond_current_results
}

results_list_cp <- cond_test_results

cond_test_results <- results_list_cp

for(i in 1:4){
  current <- cond_test_results[[i]]
  current$padj <- p.adjust(current$pvalue, "bonferroni")
  current <- current[which(is.na(current$padj) == F),]
  current <- subset(current, padj < 0.05)
  current <- current[which( !(row.names(current) %in% invalid_genes) ),]
  cond_test_results[[i]] <- current
  
}

#Load in TrAGEDy aligned datasets
WT_sce <- readRDS("path/to/TrAGEDy/aligned/WT/WT_aligned.rds")
KO_sce <- readRDS("path/to/TrAGEDy/aligned/KO/KO_aligned.rds")

WT_sce_GAM <- fitGAM(WT_sce@assays@data$counts, pseudotime = WT_sce$newPseudotime, cellWeights = rep(1, length(WT_sce$newPseudotime)))

KO_sce_GAM <- fitGAM(KO_sce@assays@data$counts, pseudotime = KO_sce$newPseudotime, cellWeights = rep(1, length(KO_sce$newPseudotime)))

WT_sce$cell_id <- colnames(WT_sce)
KO_sce$cell_id <- colnames(KO_sce)

WT_sce <- subset(WT_sce, , Status == "align")
KO_sce <- subset(KO_sce, , Status == "align")

aligned_WT <- WT_sce$cell_id
aligned_KO <- KO_sce$cell_id

WT_sce_GAM$cell_id <- colnames(WT_sce_GAM)
KO_sce_GAM$cell_id <- colnames(KO_sce_GAM)

WT_sce_GAM <- subset(WT_sce_GAM, , cell_id %in% aligned_WT)
KO_sce_GAM <- subset(KO_sce_GAM, , cell_id %in% aligned_KO)

g <- paste0("Tbrucei---","Tb927.9.5900")

obj1_WT <- plotSmoothers(WT_sce_GAM, assays(WT_sce_GAM)$counts, gene = g, lwd = 4)
obj1_KO <- plotSmoothers(KO_sce_GAM, assays(KO_sce_GAM)$counts, gene = g, lwd = 4)

ggplot() + geom_point(aes(obj1_KO$data$time, log1p(obj1_KO$data$gene_count), col = "KO")) + geom_line(aes(obj1_KO$layers[[2]]$data$time, obj1_KO$layers[[2]]$data$gene_count, col = "KO"), lwd = 4) +
  geom_borderline(aes(obj1_KO$layers[[2]]$data$time, obj1_KO$layers[[2]]$data$gene_count, col = "KO"), size=2, bordercolour = "black") + 
  geom_point(aes(obj1_WT$data$time, log1p(obj1_WT$data$gene_count), col = "WT")) + geom_line(aes(obj1_WT$layers[[2]]$data$time, obj1_WT$layers[[2]]$data$gene_count, col = "WT"), lwd = 4) +
  geom_borderline(aes(obj1_WT$layers[[2]]$data$time, obj1_WT$layers[[2]]$data$gene_count, col = "WT"), size=2, bordercolour = "black")+
  xlab("Pseudotime") +ylab("Log(expression + 1)") + ggtitle("Glutamate dehydrogenase") + theme(text = element_text(size = 30)) 


#2-amino-3-ketobutyrate coenzyme A ligase (KBL)
g <- paste0("Tbrucei---","Tb927.8.6060")

obj1_WT <- plotSmoothers(WT_sce_GAM, assays(WT_sce_GAM)$counts, gene = g, lwd = 4)
obj1_KO <- plotSmoothers(KO_sce_GAM, assays(KO_sce_GAM)$counts, gene = g, lwd = 4)

ggplot() + geom_point(aes(obj1_KO$data$time, log1p(obj1_KO$data$gene_count), col = "KO")) + geom_line(aes(obj1_KO$layers[[2]]$data$time, obj1_KO$layers[[2]]$data$gene_count, col = "KO"), lwd = 4) +
  geom_borderline(aes(obj1_KO$layers[[2]]$data$time, obj1_KO$layers[[2]]$data$gene_count, col = "KO"), size=2, bordercolour = "black") + 
  geom_point(aes(obj1_WT$data$time, log1p(obj1_WT$data$gene_count), col = "WT")) + geom_line(aes(obj1_WT$layers[[2]]$data$time, obj1_WT$layers[[2]]$data$gene_count, col = "WT"), lwd = 4) +
  geom_borderline(aes(obj1_WT$layers[[2]]$data$time, obj1_WT$layers[[2]]$data$gene_count, col = "WT"), size=2, bordercolour = "black")+
  xlab("Pseudotime") +ylab("Log(expression + 1)") + ggtitle("AKCT") + theme(text = element_text(size = 30)) 


g <- paste0("Tbrucei---","Tb927.6.2790")

obj1_WT <- plotSmoothers(WT_sce_GAM, assays(WT_sce_GAM)$counts, gene = g, lwd = 4)
obj1_KO <- plotSmoothers(KO_sce_GAM, assays(KO_sce_GAM)$counts, gene = g, lwd = 4)

ggplot() + geom_point(aes(obj1_KO$data$time, log1p(obj1_KO$data$gene_count), col = "KO")) + geom_line(aes(obj1_KO$layers[[2]]$data$time, obj1_KO$layers[[2]]$data$gene_count, col = "KO"), lwd = 4) +
  geom_borderline(aes(obj1_KO$layers[[2]]$data$time, obj1_KO$layers[[2]]$data$gene_count, col = "KO"), size=2, bordercolour = "black") + 
  geom_point(aes(obj1_WT$data$time, log1p(obj1_WT$data$gene_count), col = "WT")) + geom_line(aes(obj1_WT$layers[[2]]$data$time, obj1_WT$layers[[2]]$data$gene_count, col = "WT"), lwd = 4) +
  geom_borderline(aes(obj1_WT$layers[[2]]$data$time, obj1_WT$layers[[2]]$data$gene_count, col = "WT"), size=2, bordercolour = "black")+
  xlab("Pseudotime") +ylab("Log(expression + 1)") + ggtitle("L-threonine 3-dehydrogenase") + theme(text = element_text(size = 30)) 

g <- paste0("Tbrucei---","Tb927.6.4440")

obj1_WT <- plotSmoothers(WT_sce_GAM, assays(WT_sce_GAM)$counts, gene = g, lwd = 4)
obj1_KO <- plotSmoothers(KO_sce_GAM, assays(KO_sce_GAM)$counts, gene = g, lwd = 4)

ggplot() + geom_point(aes(obj1_KO$data$time, log1p(obj1_KO$data$gene_count), col = "KO")) + geom_line(aes(obj1_KO$layers[[2]]$data$time, obj1_KO$layers[[2]]$data$gene_count, col = "KO"), lwd = 4) +
  geom_borderline(aes(obj1_KO$layers[[2]]$data$time, obj1_KO$layers[[2]]$data$gene_count, col = "KO"), size=2, bordercolour = "black") + 
  geom_point(aes(obj1_WT$data$time, log1p(obj1_WT$data$gene_count), col = "WT")) + geom_line(aes(obj1_WT$layers[[2]]$data$time, obj1_WT$layers[[2]]$data$gene_count, col = "WT"), lwd = 4) +
  geom_borderline(aes(obj1_WT$layers[[2]]$data$time, obj1_WT$layers[[2]]$data$gene_count, col = "WT"), size=2, bordercolour = "black")+
  xlab("Pseudotime") +ylab("Log(expression + 1)") + ggtitle("RBP42") + theme(text = element_text(size = 30)) 

g <- paste0("Tbrucei---","Tb927.10.12800")

obj1_WT <- plotSmoothers(WT_sce_GAM, assays(WT_sce_GAM)$counts, gene = g, lwd = 4)
obj1_KO <- plotSmoothers(KO_sce_GAM, assays(KO_sce_GAM)$counts, gene = g, lwd = 4)

ggplot() + geom_point(aes(obj1_KO$data$time, log1p(obj1_KO$data$gene_count), col = "KO")) + geom_line(aes(obj1_KO$layers[[2]]$data$time, obj1_KO$layers[[2]]$data$gene_count, col = "KO"), lwd = 4) +
  geom_borderline(aes(obj1_KO$layers[[2]]$data$time, obj1_KO$layers[[2]]$data$gene_count, col = "KO"), size=2, bordercolour = "black") + 
  geom_point(aes(obj1_WT$data$time, log1p(obj1_WT$data$gene_count), col = "WT")) + geom_line(aes(obj1_WT$layers[[2]]$data$time, obj1_WT$layers[[2]]$data$gene_count, col = "WT"), lwd = 4) +
  geom_borderline(aes(obj1_WT$layers[[2]]$data$time, obj1_WT$layers[[2]]$data$gene_count, col = "WT"), size=2, bordercolour = "black")+
  xlab("Pseudotime") +ylab("Log(expression + 1)") + ggtitle("ZC3H38") + theme(text = element_text(size = 30)) 


g <- paste0("Tbrucei---","Tb927.7.6860")

obj1_WT <- plotSmoothers(WT_sce_GAM, assays(WT_sce_GAM)$counts, gene = g, lwd = 4)
obj1_KO <- plotSmoothers(KO_sce_GAM, assays(KO_sce_GAM)$counts, gene = g, lwd = 4)

ggplot() + geom_point(aes(obj1_KO$data$time, log1p(obj1_KO$data$gene_count), col = "KO")) + geom_line(aes(obj1_KO$layers[[2]]$data$time, obj1_KO$layers[[2]]$data$gene_count, col = "KO"), lwd = 4) +
  geom_borderline(aes(obj1_KO$layers[[2]]$data$time, obj1_KO$layers[[2]]$data$gene_count, col = "KO"), size=2, bordercolour = "black") + 
  geom_point(aes(obj1_WT$data$time, log1p(obj1_WT$data$gene_count), col = "WT")) + geom_line(aes(obj1_WT$layers[[2]]$data$time, obj1_WT$layers[[2]]$data$gene_count, col = "WT"), lwd = 4) +
  geom_borderline(aes(obj1_WT$layers[[2]]$data$time, obj1_WT$layers[[2]]$data$gene_count, col = "WT"), size=2, bordercolour = "black")+
  xlab("Pseudotime") +ylab("Log(expression + 1)") + ggtitle("ESAG5") + theme(text = element_text(size = 30)) 






