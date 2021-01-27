
# Plot erythroid/megakaryoblastic vs other AML RNAi dependency scores (Figures 2C-D)

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

# read data
data <- read.csv("portal-RNAi_merged-2018-06-21.csv", header = TRUE)
rownames(data) <- data[,1]
data[,1] <- NULL
data <- as.matrix(data)

# read annotations
annot <- read.csv(file = "CCLE_sample_info_file_2012-10-18_modified.csv", header = TRUE, sep = ";") # read in CCLE annotation file (2012 publication) with lineage data added
cell_line_names <- read.delim(file = "CCLE_cell_line_names_AUG2017.txt", sep = "\t") # read in cell line names from CCLE portal to get updated cell line names
colnames(data) <- gsub("^X", "", colnames(cell_line_names)[match(gsub("\\.", "", gsub("\\.2", "", gsub("_.*", "", gsub("NHA_HT_DD", "NHAHTDD", gsub("CACO.2", "CACO2", gsub("^X", "", toupper(colnames(data)))))))), gsub("_.*", "", gsub("^X", "", colnames(cell_line_names))))])

# subset annotation table to match data matrix
annot_rnai <- annot[match(colnames(data), annot$CCLE.name),]

# select AML and CML cell lines
rnai_aml <- data[,annot_rnai$Hist.Subtype1 %in% c("acute_myeloid_leukaemia", "blast_phase_chronic_myeloid_leukaemia")]
rnai_aml[rnai_aml == "NaN"] <- NA

# define erythroid/megakaryoblastic and other cell lines
erythroid <- colnames(rnai_aml)[grepl("CMK|CMK115|F36P|HEL|HEL9217|MOLM16|LAMA84|K562", colnames(rnai_aml))]
erythroid
other <- colnames(rnai_aml)[!colnames(rnai_aml) %in% erythroid]
other

# include only genes with more than one non-NA value per group
rnai_aml <- rnai_aml[rowSums(!is.na(rnai_aml[,erythroid]))>1 & rowSums(!is.na(rnai_aml[,other]))>1,]

# test for differential dependency between erythroid/megakaryoblastic and other cell lines
t_test <- function(GENE){
  data <- rnai_aml
  p <- t.test(as.numeric(data[GENE,erythroid]), as.numeric(data[GENE,other]))$p.value
  mean1 <- mean(as.numeric(na.omit(data[GENE,erythroid])))
  mean2 <- mean(as.numeric(na.omit(data[GENE,other])))
  return(data.frame(p, mean1, mean2))
}

result <- do.call(rbind, lapply(rownames(rnai_aml), t_test))
result$q <- p.adjust(result$p, method = "fdr") # adjust p values
result$dCS <- result$mean2 - result$mean1
result$log10p <- -log10(result$p)
result$log10q <- -log10(result$q)
result$gene <- rownames(rnai_aml)
result_t <- result[,c(8, 1:7)]
result_t <- result_t %>% arrange(p)

# volcano plot (Figure 2C)
cutoff <- result_t$q < 0.25
cutoff2 <- result_t$p < 0.05 & abs(result_t$dCS) > 0.5
selected <- result_t$gene %in% c("BCL2L1", "BCL2", "MCL1")

p <- ggplot(result_t, aes(x = dCS, y = log10p)) +
  geom_point(color = "grey70", size = 1.25) +
  geom_point(data = result_t[selected|cutoff2,], size = 4,
             color = ifelse(result_t[selected|cutoff2,"dCS"] > 0, brewer.pal(11, "RdBu")[2], brewer.pal(11, "RdBu")[10])) +
  geom_text_repel(data = result_t[selected|cutoff2,], aes(label = as.character(gene)), point.padding = 0.15) +
  ylab("-log10(P value)") +
  xlab("Differential RNAi dependency score") +
  geom_hline(yintercept = 2.956509, linetype = "dashed") +
  theme_cowplot()


pdf("Fig2C_RNAi_volcanoplot.pdf", height = 5, width = 5)
p
dev.off()

# print result table
write.table(result_t, "TableS4_RNAi_erythroid_megakaryoblastic.txt", quote = F, row.names = F, sep = "\t")


# plot individual genes across cell lines

# BCL2L1 (Figure 2D)
data_BCL2L1 <- melt(rnai_aml["BCL2L1",])
rownames(data_BCL2L1) <- gsub("_HAEMATOPOIETIC_AND_LYMPHOID_TISSUE", "", rownames(data_BCL2L1))
data_BCL2L1$variable <- rownames(data_BCL2L1)
data_BCL2L1$Subtype <- "Other AML"
data_BCL2L1$Subtype[rownames(data_BCL2L1)%in%c("F36P", "HEL", "HEL9217")] <- "Erythroid AML"
data_BCL2L1$Subtype[rownames(data_BCL2L1)%in%c("CMK", "CMK115", "MOLM16")] <- "Megakaryoblastic AML"
data_BCL2L1$Subtype[rownames(data_BCL2L1)%in%c("K562", "LAMA84")] <- "Erythroid CML"

pdf("Fig2D_RNAi_barplot_BCL2L1.pdf", height = 4, width = 6)
ggplot(data_BCL2L1, aes(y = reorder(variable, -value), x = value, fill = Subtype)) +
  geom_bar(stat = "identity") +
  theme_cowplot() +
  theme(axis.text.x = element_text(angle = 0, hjust = 1, vjust = 0.5),
        plot.title = element_text(hjust = 0, face = "plain")) +
  scale_fill_manual(values = c("Erythroid AML" = brewer.pal(11, "RdBu")[2], "Megakaryoblastic AML" = "#d35f5f", "Erythroid CML" = "#ff5555", "Other AML" = "grey50")) +
  xlab("Dependency score") +
  ylab("") +
  ggtitle("BCL2L1 RNAi") +
  labs(fill = "")
dev.off()



