
# Plot volcano plot of Hemap erythroleukemia vs other hematological malignancies (Figure S3A)

# load libraries
library(parallel)
library(gridExtra)
library(reshape2)
library(ggplot2)
library(ggrepel)
library(cowplot)
library(dplyr)
library(RColorBrewer)
library(ggrastr)


source("Rscript_plotting_functions_hemap_boxplot.R")
source("statistics_wrappers.R")

# statistics function
FUN_STATISTICS <- function(genelist, data, logicalVectors, logicalVector_normals, ALTERNATIVE="greater", adj.method="BH", CORES=1){
  library(parallel)
  data=do.call(rbind, mclapply(genelist, TestGeneWilcox, data, logicalVectors, logicalVector_normals, ALTERNATIVE, mc.cores = CORES))
  data$adj.p=p.adjust(data$P.value, method = adj.method)
  
  data$adj.method=adj.method
  
  # signifcode
  data$signifCode=""
  data$signifCode[as.numeric(data$adj.p)<0.1]="."
  data$signifCode[as.numeric(data$adj.p)<0.05]="*"
  data$signifCode[as.numeric(data$adj.p)<0.01]="**"
  data$signifCode[as.numeric(data$adj.p)<0.001]="***"
  data$signifCode[as.numeric(data$adj.p)<0.0001]="****"
  
  if(adj.method=="BH"){
    colnames(data)[colnames(data)%in%"adj.p"]="FDR"
  }
  
  return(data)
}

# load Hemap gene expression data
data=get(load("data9544_with_gene_symbols.RData"))

# load Hemap annotations
annot = get(load("Hemap_immunology_Annotations.Rdata"))

df=t(data[annot[,1],])

# define groups for testing
logicalVectors_erythroid=GET_LOGICAL(annovector = list(annot$Category.specifying.subtype), filterv = annot$Category.specifying.lineage.tumor.origin=="AML"&annot$Category.specifying.subtype=="M6"&annot$Sample.type=="Cancer")
logicalVectors_nonerythroid=GET_LOGICAL(annovector = list(annot$Sample.type), filterv = annot$Category.specifying.lineage.tumor.origin!="M6"&annot$Sample.type=="Cancer")

# test differential expression
stat <- FUN_STATISTICS(rownames(df), df, logicalVectors_erythroid, logicalVectors_nonerythroid, ALTERNATIVE="two.sided", CORES=1)

# print table for supplement
stat_suppl <- stat[,c(1,2,3,6)]
colnames(stat_suppl) <- c("Gene", "Fold change (log2)", "p", "FDR")
stat_suppl <- stat_suppl %>% arrange(p)

write.table(stat_suppl, "M6_vs_other_differential_expression.txt", row.names = F, quote = F, sep = "\t")

# plot volcano plot
  stat_log <- stat
  stat_log$FC <- log2(stat_log$FC)
  stat_log$FDR <- -log10(stat_log$FDR)
  stat_log$genelabel <- ifelse(stat$gene=="BCL2L1", "BCL2L1", "other")
  
  pdf("FigS3A_M6_vs_other_volcanoplot.pdf", width = 6, height = 4)
  ggplot(stat_log, aes(x=FC, y=FDR, label = gene)) +
    geom_point_rast(color = "grey70", size = 0.5) +
    geom_point(data = stat_log[stat_log$gene == "BCL2L1",], color = brewer.pal(11, "RdBu")[2], size = 4) +
    geom_point(data = stat_log[stat_log$gene %in% c("BCL2", "MCL1"),], color = brewer.pal(11, "RdBu")[10], size = 4) +
    geom_text_repel(data = stat_log[stat_log$gene %in% c("BCL2L1", "BCL2", "MCL1"),], aes(label=ifelse(gene%in%c("BCL2L1", "BCL2", "MCL1"), as.character(gene),'')), point.padding = 0.2) +
    scale_color_manual(values = c("BCL2L1"= "red", "other"="grey70")) +
    ylab("-log10(FDR)") +
    xlab("Fold change erythroleukemia vs\nother hematological malignancies (log2)") +
    guides(color = F) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
    theme_cowplot()
  dev.off()
