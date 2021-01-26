
# Plot Wang. et al CRISPR screen data (Figure S2)

# load libraries
library(ggplot2)
library(cowplot)
library(reshape2)
library(matrixStats)
library(readxl)
library(dplyr)
library(ggrepel)
library(RColorBrewer)
library(ggrastr)
library(biomaRt)

# read data
data <- as.data.frame(read_excel("Cell 2017 Wang.xlsx"))
colnames(data) <- data[1,]
colnames(data) <- toupper(gsub("-|;|/", "", colnames(data)))
rownames(data) <- data[,1]
data <- data[-1,-c(1,2,9)] # remove annotation columns and NB4 replicate
data <- data.matrix(data)

colnames(data)[6] <- "NB4"


# obtain gene annotations
ensembl_hs_mart <- useMart(biomart="ensembl", dataset="hsapiens_gene_ensembl")
ensembl_df <- getBM(attributes=c("ensembl_gene_id", "hgnc_symbol", "description", "chromosome_name", "start_position"),
                    mart=ensembl_hs_mart)
gene_annot <- ensembl_df[match(rownames(data), ensembl_df$hgnc_symbol),] # match biomart data using Ensembl gene symbols


# calculate differential CRISPR scores between erythroid/megakaryocytic and other AMLs
erythroid_vs_other <- rowMeans(data[,!colnames(data)%in%c("HEL", "TF1")]) - rowMeans(data[,c("HEL", "TF1")])
erythroid_vs_other <- as.data.frame(erythroid_vs_other)
rownames(erythroid_vs_other) <- rownames(data)
erythroid_vs_other$gene <- rownames(data)
erythroid_vs_other$position_rank <- rank(gene_annot$start_position)

# order by differential CS
erythroid_vs_other <- erythroid_vs_other[order(erythroid_vs_other$erythroid_vs_other),]

# color selected genes
erythroid_vs_other$coloredgene <- erythroid_vs_other$gene
erythroid_vs_other$coloredgene[!erythroid_vs_other$coloredgene%in%c("BCL2L1", "BCL2", "MCL1")] <- "Other"


# plot scatter plot of genes ranked by differential CS and genomic position
topbottom10 <- erythroid_vs_other[c(1:10, nrow(erythroid_vs_other):(nrow(erythroid_vs_other)-10)),]
topbottom10 <- topbottom10[!topbottom10$gene%in%c("BCL2L1", "BCL2", "MCL1"),]
bcl <- erythroid_vs_other[erythroid_vs_other$gene%in%c("BCL2L1", "BCL2", "MCL1"),]

pdf("FigS2A_Wang_CRISPR_scatter.pdf", height = 5, width = 6)
ggplot(erythroid_vs_other,
       aes(x = position_rank,
           y = erythroid_vs_other,
           color = coloredgene)) +
  geom_point_rast(size = 0.5) +
  geom_point(data = topbottom10, aes(x = position_rank, y = erythroid_vs_other), size = 2, color = "grey50") +
  geom_point(data = bcl, aes(x = position_rank, y = erythroid_vs_other), size = 4) +
  theme_cowplot() +
  scale_color_manual(values = c("BCL2L1" = brewer.pal(11, "RdBu")[2], "BCL2" = brewer.pal(11, "RdBu")[10], "MCL1" = brewer.pal(11, "RdBu")[10], "Other" = "grey70")) +
  geom_text_repel(data = bcl, aes(label = as.character(gene)), point.padding = 0.2, color = "black") +
  geom_text_repel(data = topbottom10, aes(label = as.character(gene)), point.padding = 0.2, color = "grey50") +
  xlab("Genes") +
  ylab("Differential CRISPR score") +
  theme(legend.position="none",
        plot.title = element_text(hjust = 0, face = "plain")) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  ggtitle("Wang et al. CRISPR screen")
dev.off()


# plot individual genes across cell lines

# BCL2L1
data_BCL2L1 <- melt(data["BCL2L1",])
data_BCL2L1$variable <- rownames(data_BCL2L1)
data_BCL2L1$Subtype <- "Other AML"
data_BCL2L1$Subtype[rownames(data_BCL2L1)%in%c("TF1", "HEL")] <- "Erythroid AML"

pdf("FigS2B_Wang_CRISPR_boxplot_BCL2L1.pdf", height = 4, width = 5)
ggplot(data_BCL2L1, aes(y = reorder(variable, -value), x = value, fill = Subtype)) +
  geom_bar(stat = "identity") +
  theme_cowplot() +
  theme(axis.text.x = element_text(angle = 0, hjust = 1, vjust = 0.5),
        plot.title = element_text(hjust = 1, face = "plain")) +
  scale_fill_manual(values = c("Erythroid AML" = brewer.pal(11, "RdBu")[2], "Other AML" = "grey70")) +
  xlab("CRISPR score") +
  ylab("") +
  ggtitle("BCL2L1 (Wang et al.)") +
  labs(fill = "")
dev.off()


