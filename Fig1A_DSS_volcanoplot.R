
# Plot volcano plot of erythroid/megakaryoblastic vs other AML cell lines (Figure 1A)

# load libraries
library(readxl)
library(ggplot2)
library(ggrepel)
library(dplyr)
library(cowplot)
library(reshape2)
library(gridExtra)
library(RColorBrewer)
library(ggpubr)

# read erythroid/megakaryoblastic leukemia cell lines DSS
path = getwd()
filenames = list.files(path = path, pattern="Original_DSS2_2018", recursive = TRUE)
filenames

readdata <- function(filename) {
  df <- as.data.frame(read_excel(paste(path, filename, sep = "/")))
  vec <- df[, 3]
  names(vec) <- df[, 2]
  return(vec)
}

dss <- do.call(cbind, lapply(filenames[!grepl("js", filenames)], readdata))
colnames(dss) <- toupper(gsub("_.*", "", gsub("^.....", "", filenames[!grepl("js", filenames)])))

# read other AML cell lines DSS
dss_other <- data.frame(read_excel("MOLM13_MV411_OCIAML2_Original_DSS2_2019-05-21.xlsx"))
rownames(dss_other) <- dss_other$DRUG_NAME
dss_other <- dss_other[,-c(1,2)]

# merge erythroid and other DSS
dss <- merge(dss, dss_other, by.x = 0, by.y = 0)
rownames(dss) <- dss$Row.names
dss$Row.names <- NULL

# remove drugs with DSS = 0
dss <- dss[!rowSums(dss)==0,]

# test for differential drug sensitivity between erythroid/megakaryoblastic and other cell lines
t_test <- function(DRUG){
  p <- t.test(dss[DRUG,c(1:6)], dss[DRUG,c(7:9)])$p.value
  mean1 <- mean(as.numeric(dss[DRUG,c(1:6)]))
  mean2 <- mean(as.numeric(dss[DRUG,c(7:9)]))
  return(data.frame(p, mean1, mean2))
}

result <- do.call(rbind, lapply(rownames(dss), t_test))
result$q <- p.adjust(result$p, method = "fdr") # adjust p values
result$dDSS <- result$mean1 - result$mean2
result$log10p <- -log10(result$p)
result$log10q <- -log10(result$q)
result$drug <- rownames(dss)
result_t <- result[,c(8, 1:7)]
result_t <- result_t %>% arrange(p)

# volcano plot
cutoff <- result_t$q < 0.05
cutoff2 <- result_t$q < 0.05 & abs(result_t$dDSS) > 15
selected <- result_t$drug %in% c("A-1155463", "A-1331852", "Venetoclax", "Navitoclax", "S-63845")

p <- ggplot(result_t, aes(x = dDSS, y = log10p)) +
  geom_point(size = ifelse(cutoff, 1.5, 1),
             color = ifelse(cutoff, "grey50", "grey70")) +
  geom_point(data = result_t[selected|cutoff2,], size = 4,
             color = ifelse(result_t[selected|cutoff2,"dDSS"] > 0, brewer.pal(11, "RdBu")[2], brewer.pal(11, "RdBu")[10])) +
  geom_text_repel(aes(label = ifelse(selected|cutoff2, as.character(drug), ""))) +
  ylab("-log10(p value)") +
  xlab("Differential DSS") +
  theme_cowplot()

pdf("Fig1A_DSS_volcanoplot.pdf", height = 5, width = 5)
p
dev.off()


