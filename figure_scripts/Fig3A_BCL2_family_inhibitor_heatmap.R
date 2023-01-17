
# Plot heatmap of BCL-2 fmaily inhibitor custom plate DSS scores (Figure 3A, revised version)

# load libraries
library(readxl)
library(data.table)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(dplyr)
library(ggplot2)
library(cowplot)

# read data
dss <- read_excel("BCL2_family_custom_plate_Original_DSS2_2018-06-08.xlsx", sheet = 1)
ccle_annot <- fread("sample_info.csv", data.table = F)
ccle_gexp <- fread("revision_blood_2022/CCLE_expression.csv", data.table = F) # 22Q2 data (05/22)
ccle_mut_damaging <- fread("revision_blood_2022/CCLE_mutations_bool_damaging.csv", data.table = F) # 22Q2 data (05/22)
ccle_mut_hotspot <- fread("revision_blood_2022/CCLE_mutations_bool_hotspot.csv", data.table = F) # 22Q2 data (05/22)
ccle_mut <- fread("revision_blood_2022/CCLE_mutations.csv", data.table = F) # 22Q2 data (05/22) 


# create matrix for heatmap
mat <- as.matrix(dss[, c(3:ncol(dss))])
rownames(mat) <- dss$DRUG_NAME

# order according to Western blot
cell_lines <- c("F36P", "OCIM1", "M07", "HNT34", "TF1", "CMK", "HEL", "PL21", "SKM1", "KG1", "SH1", "THP1", "GDM1", "KASUMI1", "OCIAML3", "MV411", "SH2", "NOMO1", "HL60", "ML2", "MOLM13")
genes <- c("BCL2L1", "BCL2", "MCL1", "GATA1", "GFI1B", "TAL1", "KLF1", "NFE2")
mat <- mat[c("A-1331852", "A-1155463", "Navitoclax", "Venetoclax", "S-63845"), cell_lines]
rownames(mat) <- c("A-1331852 (BCL-XL)", "A-1155463 (BCL-XL)", "Navitoclax (BCL-XL/BCL-2)", "Venetoclax (BCL-2)", "S-63485 (MCL-1)")
ht_dss <- Heatmap(mat,
                  name = "DSS",
                  cluster_rows = F,
                  cluster_columns = F,
                  #col = colorRamp2(c(0,2.5,5,7.5,10,12.5,15,20,25,30,40), rev(brewer.pal(11, "PuOr"))),
                  col = colorRamp2(c(0,2.5,5,7.5,10,12.5,15,20,25,30,40), pals::ocean.tempo(11)),
                  rect_gp = gpar(col = "white", lty = 1, lwd = 1),
                  row_names_side = "left",
                  show_column_names = F,
                  row_names_gp = gpar(type = "Arial", fontsize = 8),
                  heatmap_legend_param = list(title = "DSS",
                                              title_gp = gpar(fontsize = 8),
                                              labels_gp = gpar(fontsize = 8),
                                              grid_height = unit(0.2, "cm"),
                                              grid_width = unit(2, "mm"))
)

# gene expression
annot <- ccle_annot[ccle_annot$stripped_cell_line_name %in% c("F36P", "OCIM1", "M07", "HNT34", "TF1", "CMK", "HEL", "PL21", "SKM1", "KG1", "SH1", "THP1", "GDM1", "KASUMI1", "OCIAML3", "MV411", "SH2", "NOMO1", "HL60", "ML2", "MOLM13"),]
gexp <- ccle_gexp[match(annot$DepMap_ID, ccle_gexp$V1),]
gexp <- gexp[!is.na(gexp$V1),]
rownames(gexp) <- annot$stripped_cell_line_name[match(gexp$V1, annot$DepMap_ID)]
colnames(gexp) <- gsub("\\ .*", "", colnames(gexp))
gexp <- gexp[cell_lines,genes]
mat_gexp <- t(gexp)
mat_scaled <- t(apply(mat_gexp, 1, scale))
colnames(mat_scaled) <- colnames(mat_gexp)

ht_gexp <- Heatmap(mat_scaled,
                   name = "Expression (Z-score)",
                   cluster_rows = F,
                   cluster_columns = F,
                   col = colorRamp2(seq(-2, 2, length.out = 11), rev(brewer.pal(11, "RdBu"))),
                   rect_gp = gpar(col = "white", lty = 1, lwd = 1),
                   row_names_side = "left",
                   show_column_names = F,
                   row_names_gp = gpar(type = "Arial", fontface = "italic", fontsize = 8),
                   heatmap_legend_param = list(title = "Expression \nZ-score",
                                               title_gp = gpar(fontsize = 8),
                                               labels_gp = gpar(fontsize = 8),
                                               grid_height = unit(0.2, "cm"),
                                               grid_width = unit(2, "mm"))
)


# TP53 mutations


mut_tp53 <- ccle_mut %>% filter(Hugo_Symbol == "TP53", DepMap_ID %in% annot$DepMap_ID)

mut_celllines <- mut_tp53 %>% select(Variant_Classification, DepMap_ID) %>% tidyr::pivot_wider(names_from = "DepMap_ID", values_from = "Variant_Classification")
colnames(mut_celllines) <- annot$stripped_cell_line_name[match(colnames(mut_celllines), annot$DepMap_ID)]

wt_celllines <- colnames(mat)[!colnames(mat) %in% colnames(mut_celllines)]

wt_celllines_vec <- rep("None", length(wt_celllines))
names(wt_celllines_vec) <- wt_celllines

mat_mut <- cbind(mut_celllines, t(as.data.frame(wt_celllines_vec)))

mat_mut$TF1 <- "Frame_Shift_Del_Splice_Site"
mat_mut$OCIM1 <- "Missense_Mutation_Splice_Site"

mat_mut <- t(unlist(mat_mut))

mat_mut <- gsub("_S", " s", gsub("_Sp", " + sp", gsub("_Sh", " sh", gsub("_Del", " deletion", gsub("_M", " m", mat_mut)))))

mat_mut <- t(as.matrix(mat_mut[,colnames(mat)]))
rownames(mat_mut) <- "TP53"


cols <- c(NineteenEightyR::sunset1(), "grey90")
names(cols) <- c("Missense mutation", "Frame shift deletion", "Splice site", "Frame shift deletion + splice site", "Missense mutation + splice site", "None")

ht_mut <- Heatmap(mat_mut,
                  name = "TP53 mutation",
                  cluster_rows = F,
                  cluster_columns = F,
                  col = cols,
                  rect_gp = gpar(col = "white", lty = 1, lwd = 1),
                  row_names_side = "left",
                  show_column_names = F,
                  row_names_gp = gpar(type = "Arial", fontface = "italic", fontsize = 8),
                  heatmap_legend_param = list(title = "TP53 mutation",
                                              title_gp = gpar(fontsize = 8),
                                              labels_gp = gpar(fontsize = 8),
                                              grid_height = unit(0.2, "cm"),
                                              grid_width = unit(2, "mm"))
)



# plot heatmap
pdf("Fig3A_BCL2_family_inhibitor_custom_heatmap_revision.pdf", height = 2.75, width = 9)
ht_gexp %v% ht_dss %v% ht_mut 
dev.off()

# drug sensitivity  vs TP53 mutations

df <- data.frame(cell_line = colnames(mat_mut),
           tp53_mut = as.character(ifelse(mat_mut=="None", "TP53 WT", "TP53-mutated")),
   cancer_type = c("M6/M7", "M6/M7", "M6/M7", "Other", "M6/M7", "M6/M7", "M6/M7", rep("Other", 14)))

df_dss <- mat  %>% as.data.frame %>% mutate(drug = rownames(mat))%>%  tidyr::pivot_longer(cols = 1:ncol(mat), names_to = "cell_line", values_to = "dss")

df <- df %>% left_join(df_dss)


# facet by cancer type and group by TP53 mtuation status

df_sd <- df  %>% 
  group_by(cancer_type, tp53_mut) %>% 
  summarize(mean = mean(dss), sd = sd(dss))

df_tp53 <- df %>% left_join(df_sd)


ggplot(df_tp53, aes(x = tp53_mut, y = dss, fill = tp53_mut)) +
  geom_bar(stat = "summary", color = "black") +
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), width = 0.3) +
  geom_jitter(width = 0.2) +
  facet_wrap(. ~ cancer_type + drug, ncol = 2, dir = "v") +
  theme_cowplot() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        strip.background = element_blank()) +
  scale_y_continuous(expand = c(0,0), limits = c(-0.25,32)) +
  scale_fill_manual(values = c("grey30", "grey70")) +
  ylab("Drug sensitivity score") +
  xlab("")
  

# facet by TP53 mtuation status and group by cancer type

df_sd <- df  %>% 
  group_by(tp53_mut, cancer_type, drug) %>% 
  summarize(mean = mean(dss), sd = sd(dss))

df_cancertype <- df %>% left_join(df_sd)

df_cancertype$drug <- factor(df_cancertype$drug, levels = c("A-1331852 (BCL-XL)", "A-1155463 (BCL-XL)", "Navitoclax (BCL-XL/BCL-2)", "Venetoclax (BCL-2)", "S-63485 (MCL-1)" ))

ggplot(df_cancertype, aes(x = cancer_type, y = dss, fill = cancer_type)) +
  geom_errorbar(aes(ymin = mean-1, ymax = mean + sd), width = 0.2) +
  geom_bar(stat = "summary", color = "black") +
  # geom_errorbar(aes(ymin = ifelse(mean - sd < 0, 0, mean - sd), ymax = mean + sd), width = 0.3) +
  geom_jitter(width = 0.2) +
  # facet_wrap(. ~ tp53_mut + drug, ncol = 2, dir = "v") +
  facet_grid(drug ~ tp53_mut, switch = "y") +
  theme_cowplot() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        strip.background = element_blank(),
        strip.placement = "outside") +
  scale_y_continuous(expand = c(0,0), limits = c(-1,45)) +
  scale_fill_manual(values = c("darkred", "grey70")) +
  ylab("Drug sensitivity score") +
  xlab("") +
  guides(fill = "none") +
  ggpubr::stat_compare_means(comparisons = list(c("M6/M7", "Other")), method = "wilcox.test")

ggsave("FigS3A_TP53mut_barplot.pdf", height = 13, width = 3.5)


ggplot(df_cancertype, aes(x = cancer_type, y = dss, fill = cancer_type)) +
  geom_boxplot(color = "black", outlier.shape = NA) +
  geom_jitter(width = 0.2) +
  # facet_wrap(. ~ tp53_mut + drug, ncol = 2, dir = "v") +
  facet_grid(drug ~ tp53_mut, switch = "y") +
  theme_cowplot() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        strip.background = element_blank(),
        strip.placement = "outside") +
  scale_y_continuous(expand = c(0,0), limits = c(-0.25,32)) +
  scale_fill_manual(values = c("darkred", "grey70")) +
  ylab("Drug sensitivity score") +
  xlab("") +
  guides(fill = "none")

ggsave("FigS3A_TP53mut_boxplot.pdf", height = 9, width = 3)


