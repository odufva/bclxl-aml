
# Plot BCL-2 family expression in TCGA AML data (Figure 3C)

# load libraries
library(parallel)
library(gridExtra)
library(reshape2)
library(ggplot2)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)

# functions
source("Rscript_plotting_functions_hemap_boxplot.R")

FUN_PLOT_BCL2L1=function(gene, logicalVectors, namesLV, data=NULL, matrix=NULL, col=NULL, ORDER=F, RANGE=NULL) {
  if(is.null(matrix)&is.null(data))stop("No data to plot, check data/matrix")
  
  GNAME=gsub("N:....:|:::::|DUFVA_", "", gene)
  GNAME=gsub("_", " ", GNAME)
  namesLV=gsub("Cancer_", " ", namesLV)
  
  if(is.null(col)){
    # http://tools.medialab.sciences-po.fr/iwanthue/
    col=c("#d7a85b",
          "#4d56b9",
          "#acb839",
          "#5e2883",
          "#42c87f",
          "#bf4da5",
          "#75b550",
          "#8479e6",
          "#cea632",
          "#5488e3",
          "#d38725",
          "#3e397f",
          "#a4a94e",
          "#be7cde",
          "#4d7122",
          "#8460b5",
          "#62ac6a",
          "#86275d",
          "#43c8ac",
          "#cf483f",
          "#748ed7",
          "#ca6e37",
          "#de91dc",
          "#926a26",
          "#94589d",
          "#822c17",
          "#d76092",
          "#d2745b",
          "#b24258",
          "#d35760")[seq(logicalVectors)]
  }
  
  
  if(!is.null(matrix)){
    gene2=ifelse(grepl("GEXP", gene), gene, paste0("'N:GEXP:", gene, ":::::'"))
    D=as.numeric(read.delim(pipe(paste0("grep -Fw ", gene2, " ", matrix)), row.names = 1, header=F))
  }
  
  if(!is.null(data)){
    D = data[,grepl(paste("GEXP:", paste(gene,  ":chr", sep = ""), sep = ""), colnames(data))]
  }
  
  bplot_list=lapply(logicalVectors, function(v){
    D[v]
  })
  names(bplot_list)=gsub("_", " ", namesLV)
  
  if(ORDER){
    ord=c(7,3,2,1,6,5,4,9,8)
    bplot_list=bplot_list[ord]
    
  }
  
  plots=FUNCTION_PLOT_LIST(bplot_list, gene, col, ORDER, RANGE)
  return(plots)
}


# load TCGA AML data
load("DUFVA_TCGA_AML_FM_meth.Rdata")

# select samples with gene expression data
data <- fm[,is.na(fm["C:SAMP:cancermap_cluster:::::",])==FALSE]

# transform
t.data <- as.data.frame(t(data))

# shorten M0
levels(t.data$`C:CLIN:leukemia_french_american_british_morphology_code:::::`)[levels(t.data$`C:CLIN:leukemia_french_american_british_morphology_code:::::`)=="M0_Undifferentiated"] <- "M0"

# define vcetor for plotting FAB subtypes
logicalVectors=GET_LOGICAL(annovector = list(t.data$`C:CLIN:leukemia_french_american_british_morphology_code:::::`)) 

# gene to plot
gene="BCL2L1"

pdf("Fig3B_TGCA_AML_FAB_BCL2L1_boxplot.pdf", height = 3, width = 5)
FUN_PLOT_BCL2L1(gene, logicalVectors, namesLV=names(logicalVectors), data=t.data, ORDER=T, col="grey50")
dev.off()

# heatmap of median expression of BCL2 family genes in AML FAB subtypes

#subset to BCL2 family genes
genelist <- c("BCL2L1", "BCL2", "MCL1", "BCL2L2", "BCL2A1")
genelist_gexp <- paste("GEXP:", paste(genelist,  ":chr", sep = ""), sep = "")
bcl <- t.data[,grepl(paste(genelist_gexp, collapse = "|"), colnames(t.data))]

bcl_matrix <- apply(as.matrix.noquote(bcl),2,as.numeric)
rownames(bcl_matrix) <- rownames(bcl)

# medians by FAB subtypes
bcl_fab <- aggregate(bcl_matrix, by = list(FAB = t.data$`C:CLIN:leukemia_french_american_british_morphology_code:::::`), FUN = median)

rownames(bcl_fab) <- bcl_fab$FAB
colnames(bcl_fab) <- gsub(":.*", "", gsub("N:GEXP:", "", colnames(bcl_fab)))
bcl_fab$FAB <- NULL
bcl_fab <- t(bcl_fab)

mat <- t(apply(bcl_fab, 1, scale))
rownames(mat) <- rownames(bcl_fab)
colnames(mat) <- colnames(bcl_fab)
mat <- mat[genelist,]

pdf("Fig3B_TCGA_AML_BCL_lineage.pdf", height = 1.5, width = 4)
Heatmap(mat, 
        name = "Z-score", 
        col = colorRamp2(seq(-2.5, 2.5, 0.5), rev(brewer.pal(11, "RdBu"))),
        row_names_side = "right", 
        row_names_gp = gpar(fontsize = 12, fontface = "italic"),
        show_column_names = FALSE,
        cluster_columns = FALSE,
        cluster_rows = FALSE
)
dev.off()       



