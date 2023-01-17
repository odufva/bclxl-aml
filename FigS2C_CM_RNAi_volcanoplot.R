
# Plot erythroid CML RNAi dependency scores (Figure S2C for revision)

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
rnai_cml <- data[,annot_rnai$Hist.Subtype1 %in% c("blast_phase_chronic_myeloid_leukaemia")]
rnai_cml[rnai_cml == "NaN"] <- NA

# define erythroid/megakaryoblastic and other cell lines
erythroid <- colnames(rnai_cml)[grepl("CMK|CMK115|F36P|HEL|HEL9217|MOLM16|LAMA84|K562", colnames(rnai_cml))]
erythroid
other <- colnames(rnai_cml)[!colnames(rnai_cml) %in% erythroid]
other

# include only genes with more than one non-NA value per group
rnai_cml <- rnai_cml[rowSums(!is.na(rnai_cml[,erythroid]))>1,]


# plot individual genes across cell lines

# BCL2L1 (Figure S2C)
data_BCL2L1 <- melt(rnai_cml["BCL2L1",])
rownames(data_BCL2L1) <- gsub("_HAEMATOPOIETIC_AND_LYMPHOID_TISSUE", "", rownames(data_BCL2L1))
data_BCL2L1$variable <- rownames(data_BCL2L1)
data_BCL2L1$Subtype <- "Other AML"
data_BCL2L1$Subtype[rownames(data_BCL2L1)%in%c("F36P", "HEL", "HEL9217")] <- "Erythroid AML"
data_BCL2L1$Subtype[rownames(data_BCL2L1)%in%c("CMK", "CMK115", "MOLM16")] <- "Megakaryoblastic AML"
data_BCL2L1$Subtype[rownames(data_BCL2L1)%in%c("K562", "LAMA84")] <- "Erythroid CML"

pdf("FigS2D_RNAi_CML_barplot_BCL2L1.pdf", height = 1.5, width = 6)
ggplot(data_BCL2L1, aes(y = reorder(variable, -value), x = value, fill = Subtype)) +
  geom_bar(stat = "identity") +
  theme_cowplot() +
  theme(axis.text.x = element_text(angle = 0, hjust = 1, vjust = 0.5),
        plot.title = element_text(hjust = 0, face = "plain")) +
  scale_fill_manual(values = c("Erythroid CML" = "#ff5555")) +
  xlab("Dependency score") +
  ylab("") +
  ggtitle("BCL2L1 RNAi") +
  labs(fill = "")
dev.off()



