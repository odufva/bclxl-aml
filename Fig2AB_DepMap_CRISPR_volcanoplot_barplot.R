
# Plot erythroid/megakaryoblastic vs other AML CRISPR gene effect scores (Figures 2A-B)

# load libraries
library(ggplot2)
library(cowplot)
library(reshape2)
library(matrixStats)
library(readxl)
library(dplyr)
library(ggrepel)
library(data.table)
library(RColorBrewer)
library(ggrastr)

# read sample info
annot <- fread("sample_info.csv")

# read DepMap CRISPR gene effect data
data <- fread("Achilles_gene_effect.csv")

annot <- annot[match(data$DepMap_ID, annot$DepMap_ID),]

# prepare data matrix
cell_lines <- annot$stripped_cell_line_name
data$DepMap_ID <- NULL
colnames(data) <- gsub("\\ \\(.*", "", colnames(data))
data <- data.matrix(t(data))
colnames(data) <- cell_lines

# check which erythroid/megakaryoblastic cell lines are available in Depmap
erythroid <- colnames(data)[colnames(data) %in% c("HEL", "TF1", "CMK", "OCIM1", "HEL9217", "M07", "M07E", "F36P")]
erythroid
other <- annot$stripped_cell_line_name[grepl("AML", annot$Subtype)]
other <- other[!other %in% erythroid]
other

# test for differential dependency between erythroid/megakaryoblastic and other cell lines
t_test <- function(GENE){
  p <- t.test(data[GENE,erythroid], data[GENE,other])$p.value
  mean1 <- mean(as.numeric(na.omit(data[GENE,erythroid])))
  mean2 <- mean(as.numeric(na.omit(data[GENE,other])))
  return(data.frame(p, mean1, mean2))
}

result <- do.call(rbind, lapply(rownames(data), t_test))
result$q <- p.adjust(result$p, method = "fdr") # adjust p values
result$dCS <- result$mean2 - result$mean1
result$log10p <- -log10(result$p)
result$log10q <- -log10(result$q)
result$gene <- rownames(data)
result_t <- result[,c(8, 1:7)]
result_t <- result_t %>% arrange(p)

# volcano plot (Figure 2A)
cutoff <- result_t$q < 0.25
cutoff2 <- result_t$p < 0.05 & abs(result_t$dCS) > 0.6
selected <- result_t$gene %in% c("BCL2L1", "BCL2", "MCL1")


p <- ggplot(result_t, aes(x = dCS, y = log10p)) +
  geom_point(color = "grey70", size = 1.25) +
  geom_point(data = result_t[selected|cutoff2,], size = 4,
             color = ifelse(result_t[selected|cutoff2,"dCS"] > 0, brewer.pal(11, "RdBu")[2], brewer.pal(11, "RdBu")[10])) +
  geom_text_repel(data = result_t[selected|cutoff2,], aes(label = as.character(gene)), point.padding = 0.15) +
  ylab("-log10(P value)") +
  xlab("Differential CRISPR gene effect") +
  geom_hline(yintercept = 3.243397, linetype = "dashed") +
  theme_cowplot()

pdf("Fig2A_DepMap_CRISPR_volcanoplot.pdf", height = 5, width = 5)
p
dev.off()

# print result table
write.table(result_t, "DepMap_CRISPR_erythroid_megakaryoblastic.txt", quote = F, row.names = F, sep = "\t")


# plot individual genes (BCL2L1) across cell lines as bar plot(Figure 2B)
plot_barplot <- function(GENE){
  
  data_gene <- reshape2::melt(data[GENE,c(erythroid,other)])
  data_gene$variable <- rownames(data_gene)
  data_gene$Subtype <- "Other AML"
  data_gene$Subtype[rownames(data_gene)%in%erythroid] <- "Erythroid AML"
  
  pdf(paste0("Fig2B_DepMap_CRISPR_barplot_", GENE, ".pdf"), height = 4, width = 5.5)
  p <- ggplot(data_gene, aes(y = reorder(variable, -value), x = value, fill = Subtype)) +
    geom_bar(stat = "identity") +
    theme_cowplot() +
    theme(axis.text.x = element_text(angle = 0, hjust = 1, vjust = 0.5),
          plot.title = element_text(hjust = 0, face = "italic")) +
    scale_fill_manual(values = c("Erythroid AML" = brewer.pal(11, "RdBu")[2], "Other AML" = "grey70")) +
    xlab("CRISPR gene effect") +
    ggtitle(GENE) +
    labs(fill = "") +
    ylab("")
  print(p)
  dev.off()
}

plot_barplot("BCL2L1")


