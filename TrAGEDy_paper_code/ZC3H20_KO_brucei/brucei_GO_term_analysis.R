library(stringr)
library(ggplot2)

setwd("/Users/rosslaidlaw/R/TrAGEDy_v2/Tbrucei_Zc3h20_KO/tables")

seurat_GO <- read.delim("seurat_GO_unique.tsv", sep = "\t")
tragedy_GO <- read.delim("tragedy_GO_unique.tsv", sep = "\t")
tradeseq_GO <- read.delim("tradeseq_GO_unique.tsv", sep = "\t")


#Subset where BH is less than 0.05
seurat_GO <- subset(seurat_GO, Benjamini < 0.01)
tragedy_GO <- subset(tragedy_GO, Benjamini < 0.01)
tradeseq_GO <- subset(tradeseq_GO, Benjamini < 0.01)

write.csv(subset(tragedy_GO, Benjamini < 0.01),"/Users/rosslaidlaw/R/TrAGEDy_V2/Tbrucei_Zc3h20_KO/tables/tragedy_GO_unique.csv")

sig_go_seurat <- seurat_GO$Name
sig_go_tragedy <- tragedy_GO$Name
sig_go_tradeseq <- tradeseq_GO$Name

plot_df <- data.frame("Method" = c("Seurat", "TrAGEDy", "TradeSeq"), 
                      "num_go_terms" = c( dim(seurat_GO)[1], dim(tragedy_GO)[1], dim(tradeseq_GO)[1] ))

pdf("/Users/rosslaidlaw/R/TrAGEDy_V2/Tbrucei_Zc3h20_KO/plots/GO_term_barplot_uniqueGenes.pdf")
ggplot(data=plot_df, aes(x=Method, y=num_go_terms, fill = Method)) +
  geom_bar(stat="identity") + xlab("") + ylab("Number of significant GO terms") + theme_minimal()
dev.off()

tragedy_unique_go <- subset(tragedy_GO, ID %in% setdiff(tragedy_GO$ID, c(tradeseq_GO$ID, seurat_GO$ID)))
tradeseq_unique_go <-subset(tradeseq_GO, ID %in% setdiff(tradeseq_GO$ID, c(tragedy_GO$ID, seurat_GO$ID)))
