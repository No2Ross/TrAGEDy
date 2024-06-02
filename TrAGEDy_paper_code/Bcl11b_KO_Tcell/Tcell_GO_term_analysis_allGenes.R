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

tragedy_genes <- read.csv("/Users/rosslaidlaw/R/TrAGEDy_V2/Tcell_Bcl11b_KO/tables/TrAGEDy_all.csv", row.names = 1)$x
tradeseq_genes <- read.csv("/Users/rosslaidlaw/R/TrAGEDy_V2/Tcell_Bcl11b_KO/tables/tradeseq_all.csv", row.names = 1)$x
seurat_genes <- read.csv("/Users/rosslaidlaw/R/TrAGEDy_V2/Tcell_Bcl11b_KO/tables/seurat_all.csv", row.names = 1)$x

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


dim(tradeseq_genes_go)
dim(tragedy_genes_go)
dim(seurat_genes_go)

sig_go_tradeseq <- subset(tradeseq_genes_go@result, p.adjust < 0.01)$Description

sig_go_tragedy <-subset(tragedy_genes_go@result, p.adjust < 0.01)$Description

sig_go_seurat <-subset(seurat_genes_go@result, p.adjust < 0.01)$Description

#Save out dataframe of results
write.csv(subset(tradeseq_genes_go@result, p.adjust < 0.01),"/Users/rosslaidlaw/R/TrAGEDy_V2/Tcell_Bcl11b_KO/tables/tradeseq_GO_all.csv")
write.csv(subset(tragedy_genes_go@result, p.adjust < 0.01),"/Users/rosslaidlaw/R/TrAGEDy_V2/Tcell_Bcl11b_KO/tables/tragedy_GO_all.csv")
write.csv(subset(seurat_genes_go@result, p.adjust < 0.01),"/Users/rosslaidlaw/R/TrAGEDy_V2/Tcell_Bcl11b_KO/tables/seurat_GO_all.csv")

goplot(tradeseq_genes_go)
goplot(tragedy_genes_go)
goplot(seurat_genes_go)

intersect_first <- intersect(subset(tradeseq_genes_go@result, p.adjust < 0.01)$Description, 
             subset(seurat_genes_go@result, p.adjust < 0.01)$Description)

intersect_second <- intersect(intersect_first,
                             subset(tragedy_genes_go@result, p.adjust < 0.01)$Description)

intersect_second

plot_df <- data.frame("Method" = c("Seurat", "TrAGEDy", "TradeSeq"), 
                      "num_go_terms" = c( dim(seurat_genes_go)[1], dim(tragedy_genes_go)[1], dim(tradeseq_genes_go)[1] ))

pdf("/Users/rosslaidlaw/R/TrAGEDy_V2/Tcell_Bcl11b_KO/plots/GO_term_barplot_allGenes.pdf")
ggplot(data=plot_df, aes(x=Method, y=num_go_terms, fill = Method)) +
  geom_bar(stat="identity") + xlab("") + ylab("Number of significant GO terms") + theme_minimal()
dev.off()


library(VennDiagram)
library(RColorBrewer)
myCol <- c("#619cff", "#f8766d", "#00ba38")

# Chart
venn.diagram(
  x = list(sig_go_tragedy, sig_go_seurat, sig_go_tradeseq),
  category.names = c("TrAGEDy" , "Seurat" , "Tradeseq"),
  filename = "/Users/rosslaidlaw/R/TrAGEDy_V2/Tcell_Bcl11b_KO/plots/venn_diagram_GO_term_overlap.png",
  lwd = 2,
  lty = 'blank',
  fill = myCol,
)

tradeseq_unique <- setdiff(sig_go_tradeseq, c(sig_go_tragedy, sig_go_seurat))
seurat_unique <- setdiff(sig_go_seurat, c(sig_go_tragedy, sig_go_tradeseq))
tragedy_unique <- setdiff(sig_go_tragedy, c(sig_go_tradeseq, sig_go_seurat))

tragedy_genes_go_df <- tragedy_genes_go@result
tradeseq_genes_go_df <- tradeseq_genes_go@result
seurat_genes_go_df <- seurat_genes_go@result

#Check p-value for common GO terms and see who has the lowest overall adjusted p-value
tragedy_genes_go_df_intersect <- tragedy_genes_go_df[tragedy_genes_go_df$Description %in% intersect_second,]
tradeseq_genes_go_df_intersect <- tradeseq_genes_go_df[tradeseq_genes_go_df$Description %in% intersect_second, ]
seurat_genes_go_df_intersect <- seurat_genes_go_df[seurat_genes_go_df$Description %in% intersect_second, ]

tragedy_genes_go_df_intersect <- tragedy_genes_go_df_intersect[ order( match(tragedy_genes_go_df_intersect$Description, intersect_second) ),]
seurat_genes_go_df_intersect <- seurat_genes_go_df_intersect[ order( match(seurat_genes_go_df_intersect$Description, intersect_second) ), ]
tradeseq_genes_go_df_intersect <- tradeseq_genes_go_df_intersect[ order( match(tradeseq_genes_go_df_intersect$Description, intersect_second) ), ]


tragedy_go_plot <- data.frame("Description" = intersect_second,
                                            "Pvalue" = tragedy_genes_go_df_intersect$pvalue,
                                            "Method" = "TrAGEDy")

tradeseq_go_plot <- data.frame("Description" = intersect_second,
                                            "Pvalue" = tradeseq_genes_go_df_intersect$pvalue,
                                            "Method" = "tradeSeq")

seurat_go_plot <- data.frame("Description" = intersect_second,
                                            "Pvalue" = seurat_genes_go_df_intersect$pvalue,
                                           "Method" = "Seurat")

plot_go_df <- rbind(tragedy_go_plot, rbind(tradeseq_go_plot, seurat_go_plot))

pdf("/Users/rosslaidlaw/R/TrAGEDy_V2/Tcell_Bcl11b_KO/plots/GO_shared_pvalue.pdf")
ggplot(plot_go_df) + geom_point(aes(x=-log10(Pvalue), y= Description, col = Method), size=3)
dev.off()

tragedy_genes_go_df_filter <- tragedy_genes_go_df[tragedy_genes_go@result$Description %in% tragedy_unique,]
tradeseq_genes_go_df_filter <- tradeseq_genes_go_df[tradeseq_genes_go@result$Description %in% tradeseq_unique,]
seurat_genes_go_df_filter <- seurat_genes_go_df[seurat_genes_go@result$Description %in% seurat_unique,]

