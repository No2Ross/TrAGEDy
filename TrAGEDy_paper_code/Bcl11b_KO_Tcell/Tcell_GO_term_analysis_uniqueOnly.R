library(clusterProfiler)
library(Seurat)
library(matrixStats)
library(org.Mm.eg.db)
library(ggplot2)

combined_seurat <- readRDS("/Users/rosslaidlaw/R/TrAGEDy_V2/Tcell_Bcl11b_KO/objects/phate_cell_cycle_UMI_regress_seurat_integrated.rds")

#The gene background are those that are expressed in at least 5% of cells in the dataset

count_mtx <- 1 - (rowCounts(as.matrix(combined_seurat@assays$RNA@layers$counts), value = 0) /length(colnames(combined_seurat)))
names(count_mtx) <- row.names(combined_seurat)

bkg_genes <- count_mtx[which(count_mtx > 0.05)]
bkg_genes <- bitr(names(bkg_genes), fromType = "SYMBOL",
                  toType = c("ENTREZID", "SYMBOL"),
                  OrgDb = org.Mm.eg.db)$ENTREZID

length(bkg_genes)

tragedy_genes <- read.csv("/Users/rosslaidlaw/R/TrAGEDy_V2/Tcell_Bcl11b_KO/tables/TrAGEDy_noTradeSeq_noSeurat.csv", row.names = 1)$x
tradeseq_genes <- read.csv("/Users/rosslaidlaw/R/TrAGEDy_V2/Tcell_Bcl11b_KO/tables/tradeseq_noTrAGEDy_noSeurat.csv", row.names = 1)$x
seurat_genes <- read.csv("/Users/rosslaidlaw/R/TrAGEDy_V2/Tcell_Bcl11b_KO/tables/seurat_noTrAGEDy_noTradeseq.csv", row.names = 1)$x

#Map symbol to entrez
tradeseq_gene.df <- bitr(tradeseq_genes, fromType = "SYMBOL",
                         toType = c("ENTREZID", "SYMBOL"),
                         OrgDb = org.Mm.eg.db)
tradeseq_genes_entrez <- tradeseq_gene.df$ENTREZID

tragedy_gene.df <- bitr(tragedy_genes, fromType = "SYMBOL",
                        toType = c("ENTREZID", "SYMBOL"),
                        OrgDb = org.Mm.eg.db)
tragedy_genes_entrez <- tragedy_gene.df$ENTREZID

seurat_gene.df <- bitr(seurat_genes, fromType = "SYMBOL",
                       toType = c("ENTREZID", "SYMBOL"),
                       OrgDb = org.Mm.eg.db)
seurat_genes_entrez <- seurat_gene.df$ENTREZID


tragedy_genes_go <- enrichGO(tragedy_genes_entrez, universe = bkg_genes,
                             OrgDb = "org.Mm.eg.db",
                             ont = "BP",
                             pvalueCutoff = 0.01)

tradeseq_genes_go <- enrichGO(tradeseq_genes_entrez, universe = bkg_genes,
                              OrgDb = "org.Mm.eg.db",
                              ont = "BP",
                              pvalueCutoff = 0.01)

seurat_genes_go <- enrichGO(seurat_genes_entrez, universe = bkg_genes,
                              OrgDb = "org.Mm.eg.db",
                              ont = "BP",
                            pvalueCutoff = 0.01)

sig_go_tradeseq <- subset(tradeseq_genes_go@result, p.adjust < 0.01)$Description

sig_go_tragedy <-subset(tragedy_genes_go@result, p.adjust < 0.01)$Description

sig_go_seurat <-subset(seurat_genes_go@result, p.adjust < 0.01)$Description


plot_df <- data.frame("Method" = c("Seurat", "TrAGEDy", "TradeSeq"), 
                      "num_go_terms" = c( dim(seurat_genes_go)[1], dim(tragedy_genes_go)[1], dim(tradeseq_genes_go)[1] ))

pdf("/Users/rosslaidlaw/R/TrAGEDy_V2/Tcell_Bcl11b_KO/plots/GO_term_barplot.pdf")
ggplot(data=plot_df, aes(x=Method, y=num_go_terms, fill = Method)) +
  geom_bar(stat="identity") + xlab("") + ylab("Number of significant GO terms") + theme_minimal()
dev.off()



sig_go_tradeseq <- subset(tradeseq_genes_go@result, p.adjust < 0.01)$Description

sig_go_tragedy <-subset(tragedy_genes_go@result, p.adjust < 0.01)$Description

sig_go_seurat <-subset(seurat_genes_go@result, p.adjust < 0.01)$Description

#Save out dataframe of results
write.csv(subset(tradeseq_genes_go@result, p.adjust < 0.01),"/Users/rosslaidlaw/R/TrAGEDy_V2/Tcell_Bcl11b_KO/tables/tradeseq_GO_unique.csv")
write.csv(subset(tragedy_genes_go@result, p.adjust < 0.01),"/Users/rosslaidlaw/R/TrAGEDy_V2/Tcell_Bcl11b_KO/tables/tragedy_GO_unique.csv")
write.csv(subset(seurat_genes_go@result, p.adjust < 0.01),"/Users/rosslaidlaw/R/TrAGEDy_V2/Tcell_Bcl11b_KO/tables/seurat_GO_unique.csv")


