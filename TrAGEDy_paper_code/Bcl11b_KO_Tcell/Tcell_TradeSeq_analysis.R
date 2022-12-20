library(Seurat)
library(ggplot2)
library(RColorBrewer)
library(SingleCellExperiment)
library(slingshot)
library(tradeSeq)

#Read in dataset with pseudotime values (not TrAGEDy aligned)
Tcell_sce <- readRDS("TrAGEDy_tests/zhou_Tcell/Tcell_combined_pseudotime.rds")

Tcell_sling <- SlingshotDataSet(Tcell_sce)

weights_subset <- Tcell_sling@curves$Lineage1$w[which(names(Tcell_sling@curves$Lineage1$w) %in% colnames(Tcell_sce))]

#tradeseq
set.seed(3)
icMat <- evaluateK(counts = as.matrix(assays(Tcell_sce)$counts),
                   pseudotime = Tcell_sce$slingPseudotime_1,
                   cellWeights = weights_subset,
                   conditions = factor(colData(Tcell_sce)$condition),
                   nGenes = 300,
                   k = 4:12)

Tcell_GAM <- fitGAM(Tcell_sce@assays@data$counts, pseudotime = Tcell_sce$slingPseudotime_1, cellWeights = weights_subset, conditions = as.factor(Tcell_sce$condition), nknots = 8)

Tcell_GAM$cell_id <- colnames(Tcell_GAM)
Tcell_sce$cell_id <- colnames(Tcell_sce)

WT_sce <- subset(Tcell_sce, , condition == "WT")
KO_sce <- subset(Tcell_sce, , condition == "Bcl11b_KO")

#knots are based on the quantiles
quantile_limits <- quantile(Tcell_sce$slingPseudotime_1, probs = seq(0,1,0.13)) 

#PCT
invalid_genes_list <- list()

windows <- list(c(1,2), c(2,3), c(3,4), c(4,5), c(5,6), c(6,7), c(7,8))

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

windows <- list(c(1,2), c(2,3), c(3,4), c(4,5), c(5,6), c(6,7), c(7,8))


for (i in 1:length(windows)){
  cond_current_results <- conditionTest(Tcell_GAM, knots = c(windows[[i]][1], windows[[i]][2]), l2fc = 0.5)
  cond_test_results[[i]] <- cond_current_results
}

results_list_cp <- cond_test_results

cond_test_results <- results_list_cp

for(i in 1:7){
  invalid_genes <- invalid_genes_list[[i]]
  current <- cond_test_results[[i]]
  current$padj <- p.adjust(current$pvalue, "bonferroni")
  current <- current[which(is.na(current$padj) == F),]
  current <- subset(current, padj < 0.05)
  current <- current[which( !(row.names(current) %in% invalid_genes) ),]
  cond_test_results[[i]] <- current
  
}


unique_gene <- c()

for(i in 1:7){
  unique_gene <- append(unique_gene, row.names(cond_test_results[[i]]) )
}

unique_gene <- unique(unique_gene)



for( i in 1:7){
  write.csv(cond_test_results[[i]], paste0(".../trade_seq_win",i,".csv") )
}

#Read in dataset with TrAGEDy aligned pseudotime
WT_sce <- readRDS("/datastore/Ross/zhou_Tcell/WT_aligned_allCells.rds")
KO_sce <- readRDS("/datastore/Ross/zhou_Tcell/KO_aligned_allCells.rds")

WT_sling <- SlingshotDataSet(WT_sce)

weights_subset <- WT_sling@curves$Lineage1$w[which(names(WT_sling@curves$Lineage1$w) %in% colnames(WT_sce))]


#tradeseq
set.seed(3)
icMat <- evaluateK(counts = as.matrix(assays(WT_sce)$counts),
                   pseudotime = WT_sce$newPseudotime,
                   cellWeights = rep(1, length(WT_sce$newPseudotime)),
                   nGenes = 300,
                   k = 4:12)

#Need to rerun this
WT_gam <- fitGAM(WT_sce@assays@data$counts, pseudotime = WT_sce$newPseudotime, cellWeights = rep(1, length(WT_sce$newPseudotime)), nknots = 8)

#tradeseq
set.seed(3)
icMat <- evaluateK(counts = as.matrix(assays(KO_sce)$counts),
                   pseudotime = KO_sce$newPseudotime,
                   cellWeights = rep(1, length(KO_sce$newPseudotime)),
                   nGenes = 300,
                   k = 4:12)

KO_gam <- fitGAM(KO_sce@assays@data$counts, pseudotime = KO_sce$newPseudotime, cellWeights = rep(1, length(KO_sce$newPseudotime)), nknots = 10)

library(ggborderline)

WT_sce$cell_id <- colnames(WT_sce)
KO_sce$cell_id <- colnames(KO_sce)

WT_sce <- subset(WT_sce, , Status == "align")
KO_sce <- subset(KO_sce, , Status == "align")

aligned_WT <- WT_sce$cell_id
aligned_KO <- KO_sce$cell_id

WT_gam$cell_id <- colnames(WT_gam)
KO_gam$cell_id <- colnames(KO_gam)

WT_gam <- subset(WT_gam, , cell_id %in% aligned_WT)
KO_gam <- subset(KO_gam, , cell_id %in% aligned_KO)

g <- "Lat"
obj1_WT <- plotSmoothers(WT_gam, assays(WT_gam)$counts, gene = g, lwd = 4)
obj1_KO <- plotSmoothers(KO_gam, assays(KO_gam)$counts, gene = g, lwd = 4)

ggplot() + geom_point(aes(obj1_KO$data$time, log1p(obj1_KO$data$gene_count), col = "KO")) + geom_line(aes(obj1_KO$layers[[2]]$data$time, obj1_KO$layers[[2]]$data$gene_count, col = "KO"), lwd = 4) +
  geom_borderline(aes(obj1_KO$layers[[2]]$data$time, obj1_KO$layers[[2]]$data$gene_count, col = "KO"), size=2, bordercolour = "black") + 
  geom_point(aes(obj1_WT$data$time, log1p(obj1_WT$data$gene_count), col = "WT")) + geom_line(aes(obj1_WT$layers[[2]]$data$time, obj1_WT$layers[[2]]$data$gene_count, col = "WT"), lwd = 4) +
  geom_borderline(aes(obj1_WT$layers[[2]]$data$time, obj1_WT$layers[[2]]$data$gene_count, col = "WT"), size=2, bordercolour = "black")+
  xlab("Pseudotime") +ylab("Log(expression + 1)") + ggtitle(g) + theme(text = element_text(size = 30)) 

g <- "Gsr"
obj1_WT <- plotSmoothers(WT_gam, assays(WT_gam)$counts, gene = g, lwd = 4)
obj1_KO <- plotSmoothers(KO_gam, assays(KO_gam)$counts, gene = g, lwd = 4)

ggplot() + geom_point(aes(obj1_KO$data$time, log1p(obj1_KO$data$gene_count), col = "KO")) + geom_line(aes(obj1_KO$layers[[2]]$data$time, obj1_KO$layers[[2]]$data$gene_count, col = "KO"), lwd = 4) +
  geom_borderline(aes(obj1_KO$layers[[2]]$data$time, obj1_KO$layers[[2]]$data$gene_count, col = "KO"), size=2, bordercolour = "black") + 
  geom_point(aes(obj1_WT$data$time, log1p(obj1_WT$data$gene_count), col = "WT")) + geom_line(aes(obj1_WT$layers[[2]]$data$time, obj1_WT$layers[[2]]$data$gene_count, col = "WT"), lwd = 4) +
  geom_borderline(aes(obj1_WT$layers[[2]]$data$time, obj1_WT$layers[[2]]$data$gene_count, col = "WT"), size=2, bordercolour = "black")+
  xlab("Pseudotime") +ylab("Log(expression + 1)") + ggtitle(g) + theme(text = element_text(size = 30)) 


g <- "Hivep3"
obj1_WT <- plotSmoothers(WT_gam, assays(WT_gam)$counts, gene = g, lwd = 4)
obj1_KO <- plotSmoothers(KO_gam, assays(KO_gam)$counts, gene = g, lwd = 4)

ggplot() + geom_point(aes(obj1_KO$data$time, log1p(obj1_KO$data$gene_count), col = "KO")) + geom_line(aes(obj1_KO$layers[[2]]$data$time, obj1_KO$layers[[2]]$data$gene_count, col = "KO"), lwd = 4) +
  geom_borderline(aes(obj1_KO$layers[[2]]$data$time, obj1_KO$layers[[2]]$data$gene_count, col = "KO"), size=2, bordercolour = "black") + 
  geom_point(aes(obj1_WT$data$time, log1p(obj1_WT$data$gene_count), col = "WT")) + geom_line(aes(obj1_WT$layers[[2]]$data$time, obj1_WT$layers[[2]]$data$gene_count, col = "WT"), lwd = 4) +
  geom_borderline(aes(obj1_WT$layers[[2]]$data$time, obj1_WT$layers[[2]]$data$gene_count, col = "WT"), size=2, bordercolour = "black")+
  xlab("Pseudotime") +ylab("Log(expression + 1)") + ggtitle(g) + theme(text = element_text(size = 30)) 

g <- "Gzmb"
obj1_WT <- plotSmoothers(WT_gam, assays(WT_gam)$counts, gene = g, lwd = 4)
obj1_KO <- plotSmoothers(KO_gam, assays(KO_gam)$counts, gene = g, lwd = 4)

ggplot() + geom_point(aes(obj1_KO$data$time, log1p(obj1_KO$data$gene_count), col = "KO")) + geom_line(aes(obj1_KO$layers[[2]]$data$time, obj1_KO$layers[[2]]$data$gene_count, col = "KO"), lwd = 4) +
  geom_borderline(aes(obj1_KO$layers[[2]]$data$time, obj1_KO$layers[[2]]$data$gene_count, col = "KO"), size=2, bordercolour = "black") + 
  geom_point(aes(obj1_WT$data$time, log1p(obj1_WT$data$gene_count), col = "WT")) + geom_line(aes(obj1_WT$layers[[2]]$data$time, obj1_WT$layers[[2]]$data$gene_count, col = "WT"), lwd = 4) +
  geom_borderline(aes(obj1_WT$layers[[2]]$data$time, obj1_WT$layers[[2]]$data$gene_count, col = "WT"), size=2, bordercolour = "black")+
  xlab("Pseudotime") +ylab("Log(expression + 1)") + ggtitle(g) + theme(text = element_text(size = 30)) 

g <- "Jak1"
obj1_WT <- plotSmoothers(WT_gam, assays(WT_gam)$counts, gene = g, lwd = 4)
obj1_KO <- plotSmoothers(KO_gam, assays(KO_gam)$counts, gene = g, lwd = 4)

ggplot() + geom_point(aes(obj1_KO$data$time, log1p(obj1_KO$data$gene_count), col = "KO")) + geom_line(aes(obj1_KO$layers[[2]]$data$time, obj1_KO$layers[[2]]$data$gene_count, col = "KO"), lwd = 4) +
  geom_borderline(aes(obj1_KO$layers[[2]]$data$time, obj1_KO$layers[[2]]$data$gene_count, col = "KO"), size=2, bordercolour = "black") + 
  geom_point(aes(obj1_WT$data$time, log1p(obj1_WT$data$gene_count), col = "WT")) + geom_line(aes(obj1_WT$layers[[2]]$data$time, obj1_WT$layers[[2]]$data$gene_count, col = "WT"), lwd = 4) +
  geom_borderline(aes(obj1_WT$layers[[2]]$data$time, obj1_WT$layers[[2]]$data$gene_count, col = "WT"), size=2, bordercolour = "black")+
  xlab("Pseudotime") +ylab("Log(expression + 1)") + ggtitle(g) + theme(text = element_text(size = 30)) 

g <- "Nfkb1"
obj1_WT <- plotSmoothers(WT_gam, assays(WT_gam)$counts, gene = g, lwd = 4)
obj1_KO <- plotSmoothers(KO_gam, assays(KO_gam)$counts, gene = g, lwd = 4)

ggplot() + geom_point(aes(obj1_KO$data$time, log1p(obj1_KO$data$gene_count), col = "KO")) + geom_line(aes(obj1_KO$layers[[2]]$data$time, obj1_KO$layers[[2]]$data$gene_count, col = "KO"), lwd = 4) +
  geom_borderline(aes(obj1_KO$layers[[2]]$data$time, obj1_KO$layers[[2]]$data$gene_count, col = "KO"), size=2, bordercolour = "black") + 
  geom_point(aes(obj1_WT$data$time, log1p(obj1_WT$data$gene_count), col = "WT")) + geom_line(aes(obj1_WT$layers[[2]]$data$time, obj1_WT$layers[[2]]$data$gene_count, col = "WT"), lwd = 4) +
  geom_borderline(aes(obj1_WT$layers[[2]]$data$time, obj1_WT$layers[[2]]$data$gene_count, col = "WT"), size=2, bordercolour = "black")+
  xlab("Pseudotime") +ylab("Log(expression + 1)") + ggtitle(g) + theme(text = element_text(size = 30)) 


