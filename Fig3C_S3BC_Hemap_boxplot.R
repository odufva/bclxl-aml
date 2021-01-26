
# Plot BCL2 family gene expression boxplots across Hemap subtypes (Figures 3C and S3B-C)

# load libraries
library(parallel)
library(gridExtra)
library(reshape2)
library(ggplot2)
library(ggrepel)
library(cowplot)

# load plotting functions
source("Rscript_plotting_functions_hemap_boxplot.R")

# load Hemap gene expression data
data=get(load("data9544_with_gene_symbols.RData"))

# load Hemap annotations
annot = get(load("Hemap_immunology_Annotations.Rdata"))

# subset data to annotations
data=data[annot[,1],]

# create new annotation vector for plotting boxplot

# add normal cell type annotations
annot$Category.specifying.subtype.boxplot <- annot$plotNormals

# add cancer type annotations
annot$Category.specifying.subtype.boxplot[annot$Sample.type%in%c("Cancer", "Prolif")] <- annot$Category.specifying.subtype[annot$Sample.type%in%c("Cancer", "Prolif")]

# fix annotations
annot$Category.specifying.subtype.boxplot[annot$Category.specifying.lineage.tumor.origin=="MDS"] <- "MDS"
annot$Category.specifying.subtype.boxplot[annot$Category.specifying.subtype=="AILT"] <- "AITL"
annot$Category.specifying.subtype.boxplot[annot$Category.specifying.subtype=="LC"] <- "LCH"
annot$Category.specifying.subtype.boxplot[annot$Category.specifying.subtype=="B-CLL"] <- "CLL"


## ---------------------------------------------------------------------------

genelist = c("BCL2L1", "BCL2", "MCL1")

# boxplots across Hemap subtypes
logicalVectors=GET_LOGICAL(annovector = list(annot$Category.specifying.subtype.boxplot), filterv = annot$Category.specifying.subtype.boxplot%in%c("M0", "M1", "M2", "M3", "M4", "M4eo", "M5", "M6", "M7", "CML", "MDS", "CLL", "pre-B-ALL", "MM", "T-ALL", "HCL", "DLBCL", "FL", "CHL", "NLPHL", "MCL", "MZL", "BL", "MALT", "PTCLNOS", "ALCL", "AITL", "ATL", "CTCL", "HSTCL", "ENKTL")&annot$Sample.type%in%c("Cancer", "Prolif"))
p.all=lapply(genelist, FUN_PLOT, logicalVectors, namesLV=names(logicalVectors), data=data, ORDER=T, col = "grey70")
ggsave(paste0("Fig3C_S3BC_", paste(genelist, collapse = "_"), "_Hemap_boxplot.pdf"), do.call(marrangeGrob, append(list(grobs=p.all, ncol=1, nrow = 1),list(top=NULL))), width = 200 , height = 100, units = "mm", dpi=250)

