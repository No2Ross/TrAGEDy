#Need to rerun this again before using the objects



library(Seurat)
library(ggplot2)
library(RColorBrewer)
library(SingleCellExperiment)
library(slingshot)
library(tradeSeq)


library(unixtools)

mean.fxn <- function(x) {
  return(log(x = (rowSums(x = expm1(x = x)) + 1)/NCOL(x), base = 2))
  
}

setwd("/datastore/Ross")

#set.tempdir("pooling-scRNA/")
ulimit::memory_limit(70000)

Tcell_sce <- readRDS("/datastore/Ross/TrAGEDy_V2/Tcell/objects/phate_cell_cycle_UMI_regress_sce_integrated.rds")

Tcell_sling <- SlingshotDataSet(Tcell_sce)

weights_subset <- Tcell_sling@curves$Lineage1$w[which(names(Tcell_sling@curves$Lineage1$w) %in% colnames(Tcell_sce))]


#tradeseq
# set.seed(3)
# icMat <- evaluateK(counts = as.matrix(assays(Tcell_sce)$counts),
#                    pseudotime = Tcell_sce$slingPseudotime_1,
#                    cellWeights = weights_subset,
#                    conditions = factor(colData(Tcell_sce)$condition),
#                    nGenes = 300,
#                    k = 4:12)
# 
# rm(Tcell_sling)
# rm(icMat)
# gc()
start.time <- Sys.time()


rm(Tcell_sling)
gc()

Tcell_GAM <- fitGAM(Tcell_sce@assays@data$counts, pseudotime = Tcell_sce$slingPseudotime_1, cellWeights = weights_subset, conditions = as.factor(Tcell_sce$condition), nknots = 7, parallel = F)

# saveRDS(Tcell_GAM, "/datastore/Ross/TrAGEDy_V2/Tcell/objects/Tcell_FitGAM.rds")
# 
# Tcell_GAM <- readRDS("/datastore/Ross/TrAGEDy_V2/Tcell/objects/Tcell_FitGAM.rds")

Tcell_GAM$cell_id <- colnames(Tcell_GAM)
Tcell_sce$cell_id <- colnames(Tcell_sce)

WT_sce <- subset(Tcell_sce, , condition == "WT")
KO_sce <- subset(Tcell_sce, , condition == "Bcl11b_KO")

#knots are based on the quantiles
quantile_limits <- quantile(Tcell_sce$slingPseudotime_1, probs = seq(0,1,0.13)) 

#PCT
invalid_genes_list <- list()

windows <- list(c(1,2), c(2,3), c(3,4), c(4,5), c(5,6), c(6,7))

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
  
  invalid_genes_pct <- intersect(names(invalid_KO), names(invalid_WT))
  
  #Log2fc section
  
  WT_log2fc <- mean.fxn(current_WT@assays@data@listData$logcounts)
  KO_log2fc <- mean.fxn(current_KO@assays@data@listData$logcounts)
  
  current_log2fc <- WT_log2fc - KO_log2fc
  
  invalid_genes_log2fc <- current_log2fc[current_log2fc < 0.75 & current_log2fc > -0.75]
  
  invalid_genes_list[[i]] <- unique(c(invalid_genes_pct , names(invalid_genes_log2fc)))
}


cond_test_results <- list()

windows <- list(c(1,2), c(2,3), c(3,4), c(4,5), c(5,6), c(6,7))


for (i in 1:length(windows)){
  cond_current_results <- conditionTest(Tcell_GAM, knots = c(windows[[i]][1], windows[[i]][2]), l2fc = 0)
  cond_test_results[[i]] <- cond_current_results
}

results_list_cp <- cond_test_results



cond_test_results <- results_list_cp

#adjust p_value and subset

for(i in 1:6){
  invalid_genes <- invalid_genes_list[[i]]
  current <- cond_test_results[[i]]
  current$padj <- p.adjust(current$pvalue, "bonferroni")
  current <- current[which(is.na(current$padj) == F),]
  current <- subset(current, padj < 0.05)
  current <- current[which( !(row.names(current) %in% invalid_genes) ),]
  cond_test_results[[i]] <- current
  
}

#1223.178 minutes
end_time <- Sys.time()
time_taken <- end_time - start.time

unique_gene <- c()

for(i in 1:6){
  unique_gene <- append(unique_gene, row.names(cond_test_results[[i]]) )
}

unique_gene <- unique(unique_gene)



for( i in 1:6){
  write.csv(cond_test_results[[i]], paste0("/datastore/Ross/TrAGEDy_V2/Tcell/DEG_files/trade_seq_win",i,".csv") )
}


