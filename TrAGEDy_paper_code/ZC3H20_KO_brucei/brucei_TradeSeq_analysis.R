library(Seurat)
library(ggplot2)
library(RColorBrewer)
library(SingleCellExperiment)
library(slingshot)
library(tradeSeq)

library(BiocParallel)

library(unixtools)

mean.fxn <- function(x) {
  return(log(x = (rowSums(x = expm1(x = x)) + 1)/NCOL(x), base = 2))
  
}

setwd("/datastore/Ross/TrAGEDy_V2/Tbrucei")

#set.tempdir("pooling-scRNA/")
ulimit::memory_limit(60000)

tryp_sce <- readRDS("/datastore/Ross/TrAGEDy_V2/Tbrucei/objects/WT_ZC3H20_with_phate_sce.rds")
tryp_sling <- SlingshotDataSet(tryp_sce)

weights_subset <- tryp_sling@curves$Lineage1$w[which(names(tryp_sling@curves$Lineage1$w) %in% colnames(tryp_sce))]


#tryp_sce <- subset(tryp_sce, , is.na(slingPseudotime_1) == FALSE)

#tryp_sce$slingPseudotime_2 <- NULL

# set.seed(3)
# icMat <- evaluateK(counts = as.matrix(assays(tryp_sce)$counts),
#                    pseudotime = tryp_sce$slingPseudotime_1,
#                    cellWeights = weights_subset,
#                    conditions = factor(colData(tryp_sce)$line),
#                    nGenes = 300,
#                    k = 4:12)
# 
# rm(Tcell_sling)
# rm(icMat)
# gc()

start.time = Sys.time()

rm(tryp_integrated)
rm(tryp_sling)
gc()

sce_gam <- fitGAM(tryp_sce, conditions = factor(colData(tryp_sce)$line),
              nknots = 8, verbose = FALSE)

# saveRDS(sce_gam, "/datastore/Ross/TrAGEDy_V2/Tbrucei/objects/brucei_sce_GAM.rds")
# saveRDS(tryp_sce, "/datastore/Ross/TrAGEDy_V2/Tbrucei/objects/brucei_sce.rds")
# 
# 
# sce_gam <- readRDS("/datastore/Ross/TrAGEDy_V2/Tbrucei/objects/brucei_sce_GAM.rds")
# tryp_sce <- readRDS("/datastore/Ross/TrAGEDy_V2/Tbrucei/objects/brucei_sce.rds")

sce_gam$cell_id <- colnames(sce_gam)
tryp_sce$cell_id <- colnames(tryp_sce)

WT_sce <- subset(tryp_sce, , line == "WT")
KO_sce <- subset(tryp_sce, , line == "ZC3H20_KO")

#knots are based on the quantiles
quantile_limits <- quantile(tryp_sce$slingPseudotime_1, probs = seq(0,1,0.11)) 

#PCT
invalid_genes_list <- list()

windows <- list(c(1,2), c(3,4), c(5,6), c(7,8))

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
  
  invalid_genes_log2fc <- current_log2fc[current_log2fc < 0.5 & current_log2fc > -0.5]
  
  invalid_genes_list[[i]] <- unique(c(invalid_genes_pct , names(invalid_genes_log2fc)))

}




cond_test_results <- list()


for (i in 1:length(windows)){
  cond_current_results <- conditionTest(sce_gam, knots = c(windows[[i]][1], windows[[i]][2]), l2fc = 0)
  cond_test_results[[i]] <- cond_current_results
}

results_list_cp <- cond_test_results



cond_test_results <- results_list_cp

#adjust p_value and subset

for(i in 1:4){
  invalid_genes <- invalid_genes_list[[i]]
  current <- cond_test_results[[i]]
  current$padj <- p.adjust(current$pvalue, "bonferroni")
  current <- current[which(is.na(current$padj) == F),]
  current <- subset(current, padj < 0.05)
  current <- current[which( !(row.names(current) %in% invalid_genes) ),]
  cond_test_results[[i]] <- current
  
}

end_time <- Sys.time()
time_taken <- end_time - start.time

#Time difference of 177.716 hours
#Time difference of 177.4232 hours

unique_gene <- c()

for(i in 1:4){
  unique_gene <- append(unique_gene, row.names(cond_test_results[[i]]) )
}

unique_gene <- unique(unique_gene)



for( i in 1:4){
  write.csv(cond_test_results[[i]], paste0("/datastore/Ross/TrAGEDy_V2/Tbrucei/DEG_files/trade_seq_win",i,".csv") )
}


