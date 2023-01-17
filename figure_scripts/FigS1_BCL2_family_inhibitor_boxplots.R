
# Plot box plots of DSS in erythroid/megakaryoblastic and other leukemia cell lines (Figure S1)

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

# read data
dss_custom <- read_excel("BCL2_family_custom_plate_Original_DSS2_2018-06-08.xlsx", sheet = 1) %>% select(-ID)

# prepare data frame
df <- melt(dss_custom)
colnames(df) <- gsub("variable", "cell_line", gsub("value", "DSS", gsub("DRUG_NAME", "drug", colnames(df))))
df <- df %>% mutate(cancer_type = ifelse(cell_line %in% c("CMK", "F36P", "HEL", "M07", "TF1", "OCIM1"), "Erythroid/\nmegakaryoblastic", "Other")) %>%
  filter(cell_line != "UT7") %>% # remove, Western failed
  mutate(drug_mechanism = gsub("A-1331852", "A-1331852\n(BCL-XL)",
                               gsub("A-1155463", "A-1155463\n(BCL-XL)",
                                    gsub("WEHI-539", "WEHI-539\n(BCL-XL)",
                                         gsub("Venetoclax", "Venetoclax\n(BCL-2)",
                                              gsub("S-63845", "S-63845\n(MCL-1)",
                                                   gsub("Navitoclax", "Navitoclax\n(BCL-2, BCL-XL)",
                                                        gsub("Sabutoclax", "Sabutoclax\n(pan-BCL-2 family)",
                                                             gsub("A-1210477", "A-1210477\n(MCL-1)", drug)))))))))


# boxplot function
plot_boxplot <- function(DRUG, YMAX, YLABEL){
  ggplot(df[df$drug==DRUG,], aes(x = cancer_type, y = DSS, fill = cancer_type)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(width = 0.2) +
    xlab("") +
    scale_fill_manual(values = c("Erythroid/\nmegakaryoblastic" = brewer.pal(11, "RdBu")[2], "Other" = "grey50")) +
    theme_cowplot() +
    ggtitle(df$drug_mechanism[df$drug==DRUG]) +
    theme(plot.title = element_text(hjust = 0.5, face = "plain"),
          axis.text.x = element_text(angle = 45, hjust = 1),
          plot.margin = unit(c(0,0.1,-0.5,0.25),"cm")) +
    scale_y_continuous(limits = c(-0.5,YMAX)) +
    guides(fill = F) +
    stat_compare_means(method = "wilcox.test",
                       label = "p.format",
                       label.x = 1.1,
                       label.y = YLABEL)
}


p1 <- plot_boxplot("A-1331852", 43, 41)
p2 <- plot_boxplot("A-1155463", 43, 41)
p3 <- plot_boxplot("WEHI-539", 43, 41)
p4 <- plot_boxplot("Venetoclax", 43, 41)
p5 <- plot_boxplot("S-63845", 43, 41)
p6 <- plot_boxplot("Navitoclax", 43, 41)
p7 <- plot_boxplot("Sabutoclax", 43, 41)
p8 <- plot_boxplot("A-1210477", 43, 41)

pdf("FigS1_BCL2_family_inhibitor_boxplots.pdf", height = 8, width = 10)
grid.arrange(p1, p2, p3, p4, p5, p6, p7, p8, ncol = 4)
dev.off()

