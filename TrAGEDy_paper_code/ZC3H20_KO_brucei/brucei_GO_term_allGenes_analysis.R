library(stringr)
library(ggplot2)

setwd("/Users/rosslaidlaw/R/TrAGEDy_v2/Tbrucei_Zc3h20_KO/tables")

seurat_GO <- read.delim("seurat_GO_all.tsv", sep = "\t")
tragedy_GO <- read.delim("tragedy_GO_all.tsv", sep = "\t")
tradeseq_GO <- read.delim("tradeseq_GO_all.tsv", sep = "\t")

#Subset where BH is less than 0.05
seurat_GO <- subset(seurat_GO, Benjamini < 0.01)
tragedy_GO <- subset(tragedy_GO, Benjamini < 0.01)
tradeseq_GO <- subset(tradeseq_GO, Benjamini < 0.01)

write.csv(subset(tragedy_GO, Benjamini < 0.01),"/Users/rosslaidlaw/R/TrAGEDy_V2/Tbrucei_Zc3h20_KO/tables/tragedy_GO_all_filtered.csv")
write.csv(subset(tradeseq_GO, Benjamini < 0.01),"/Users/rosslaidlaw/R/TrAGEDy_V2/Tbrucei_Zc3h20_KO/tables/tradeseq_GO_all_filtered.csv")
write.csv(subset(seurat_GO, Benjamini < 0.01),"/Users/rosslaidlaw/R/TrAGEDy_V2/Tbrucei_Zc3h20_KO/tables/seurat_GO_all_filtered.csv")


sig_go_seurat <- seurat_GO$Name
sig_go_tragedy <- tragedy_GO$Name
sig_go_tradeseq <- tradeseq_GO$Name

plot_df <- data.frame("Method" = c("Seurat", "TrAGEDy", "TradeSeq"), 
                      "num_go_terms" = c( dim(seurat_GO)[1], dim(tragedy_GO)[1], dim(tradeseq_GO)[1] ))

pdf("/Users/rosslaidlaw/R/TrAGEDy_V2/Tbrucei_Zc3h20_KO/plots/GO_term_barplot_allGenes.pdf")
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
  filename = "/Users/rosslaidlaw/R/TrAGEDy_v2/Tbrucei_Zc3h20_KO/plots/venn_diagram_allGenes_GO_term_overlap.png",
  lwd = 2,
  lty = 'blank',
  fill = myCol,
)

tragedy_unique_go <- subset(tragedy_GO, ID %in% setdiff(tragedy_GO$ID, c(tradeseq_GO$ID, seurat_GO$ID)))
tradeseq_unique_go <-subset(tradeseq_GO, ID %in% setdiff(tradeseq_GO$ID, c(tragedy_GO$ID, seurat_GO$ID)))

intersect_first <- intersect(tragedy_GO$Name, seurat_GO$Name)
intersect_second <- intersect(intersect_first, tradeseq_GO$Name)

#Check p-value for common GO terms and see who has the lowest overall adjusted p-value
tragedy_genes_go_df_intersect <- tragedy_GO[tragedy_GO$Name %in% intersect_second,]
tradeseq_genes_go_df_intersect <- tradeseq_GO[tradeseq_GO$Name %in% intersect_second, ]
seurat_genes_go_df_intersect <- seurat_GO[seurat_GO$Name %in% intersect_second, ]

tragedy_genes_go_df_intersect <- tragedy_genes_go_df_intersect[ order( match(tragedy_genes_go_df_intersect$Name, intersect_second) ),]
seurat_genes_go_df_intersect <- seurat_genes_go_df_intersect[ order( match(seurat_genes_go_df_intersect$Name, intersect_second) ), ]
tradeseq_genes_go_df_intersect <- tradeseq_genes_go_df_intersect[ order( match(tradeseq_genes_go_df_intersect$Name, intersect_second) ), ]


tragedy_go_plot <- data.frame("Description" = intersect_second,
                              "Pvalue" = tragedy_genes_go_df_intersect$P.value,
                              "Method" = "TrAGEDy")

tradeseq_go_plot <- data.frame("Description" = intersect_second,
                               "Pvalue" = tradeseq_genes_go_df_intersect$P.value,
                               "Method" = "tradeSeq")

seurat_go_plot <- data.frame("Description" = intersect_second,
                             "Pvalue" = seurat_genes_go_df_intersect$P.value,
                             "Method" = "Seurat")

plot_go_df <- rbind(tragedy_go_plot, rbind(tradeseq_go_plot, seurat_go_plot))

pdf("/Users/rosslaidlaw/R/TrAGEDy_V2/Tbrucei_Zc3h20_KO/plots/GO_shared_pvalue.pdf")
ggplot(plot_go_df) + geom_point(aes(x=-log10(Pvalue), y= Description, col = Method), size=3)
dev.off()
