
# Plot heatmap of BCL-2 fmaily inhibitor custom plate DSS scores (Figure 3A)

# load libraries
library(readxl)
library(data.table)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)

# read data
dss <- read_excel("BCL2_family_custom_plate_Original_DSS2_2018-06-08.xlsx", sheet = 1)
ccle_annot <- fread("sample_info.csv", data.table = F)
ccle_gexp <- fread("CCLE_expression.csv", data.table = F)

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
        col = colorRamp2(c(0,2.5,5,7.5,10,12.5,15,20,25,30,40), rev(brewer.pal(11, "PuOr"))),
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

# plot heatmap
pdf("Fig3A_BCL2_family_inhibitor_custom_heatmap.pdf", height = 2.5, width = 9)
ht_gexp %v% ht_dss
dev.off()

