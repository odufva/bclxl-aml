
# Analyze AML-1 patient flow DSRT data (Figure 5C-E, S6B)
# Live singlets gated in FlowJo

# load libraries
library(readxl)
library(dplyr)
library(CATALYST)
library(ggplot2)
library(cowplot)
library(gridExtra)
library(wesanderson)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)


md <- "FC_DSRT_metadata_AML1.xlsx"
md <- read_excel(md)
head(data.frame(md))

fs <- read.flowSet(path = "AML1_singlets_live/", pattern = ".fcs")

panel <- "FC_DSRT_panel.xlsx"
panel <- read_excel(panel)
head(data.frame(panel))

# spot check that all panel columns are in the flowSet object
all(panel$fcs_colname %in% colnames(fs))

# specify levels for conditions & sample IDs to assure desired ordering
md$condition <- factor(md$condition, levels = c("DMSO", "C1", "C2", "C3", "C4", "C5", "C6", "C7"))
md$sample_id <- factor(md$sample_id, 
                       levels = md$sample_id[order(md$condition)])

# construct SingleCellExperiment
sce <- prepData(fs, panel, md, features = panel$fcs_colname)

# plot heatmap of marker expression in metaclusters
set.seed(1234)
sce <- cluster(sce, features = "type",
               xdim = 10, ydim = 10, maxK = 20, seed = 1234)

p <- plotExprHeatmap(sce, features = "type", 
                     by = "cluster_id", k = "meta20", 
                     bars = TRUE, perc = TRUE)

pdf("FigS6B_FC_DSRT_median_expression_clusters_heatmap.pdf")
p
dev.off()


# run UMAP on at most 1000 cells per sample
set.seed(1234)
sce <- runDR(sce, "TSNE", cells = 1e3, features = "type")
sce <- runDR(sce, "UMAP", cells = 1e3, features = "type")


# plot metaclusters colored on UMAP
p <- plotDR(sce, "UMAP", color_by = "meta20")
lgd <- get_legend(p)
p <- p + theme(legend.position = "none")
plot_grid(p1, p2, lgd, nrow = 1, rel_widths = c(5, 5, 2))
ggsave("umap_clusters.pdf", height = 7, width = 5)


# cluster merging 
merging_table1 <- "FC_DSRT_cluster_merging_AML1.xlsx"
merging_table1 <- read_excel(merging_table1)
head(data.frame(merging_table1))

# convert to factor with merged clusters in desired order
merging_table1$new_cluster <- factor(merging_table1$new_cluster, 
                                     levels = c("Blasts", "Lymphocytes", "Erythrocytes",
                                                "Other"))

# apply manual merging
sce <- mergeClusters(sce, k = "meta20", 
                     table = merging_table1, id = "merging1")

# Plot figures
pal <- c(wes_palette("Cavalcanti1", 5), wes_palette("GrandBudapest1", 4), wes_palette("Royal2", 4), "grey70")
names(pal) <- c("CMP", "GMP", "HSC", "Monocytes", "MPP", "Neutrophils", "Blasts", "Megakaryocytes", "Erythrocytes",
                "Lymphocytes", "NK cells", "Plasma cells", "CLP", "Other")

# UMAPs (Figure 5E)
p1 <- plotDR(sce, "UMAP", color_by = "CD34") + theme_void()
p2 <- plotDR(sce, "UMAP", color_by = "CD117") + theme_void()
p3 <- plotDR(sce, "UMAP", color_by = "CD38") + theme_void()
p4 <- plotDR(sce, "UMAP", color_by = "CD45") + theme_void()

plot_grid(p1, p2, p3, p4, ncol = 2)
ggsave("Fig5C_umap_markers_selected.pdf", height = 4, width = 5)

plotDR(sce, "UMAP", color_by = "merging1") + theme_void() +  scale_color_manual(values = pal)
ggsave("Fig5C_umap_clusters.pdf", height = 3, width = 4)

plotDR(sce[,sce$sample_id%in%c("DMSO1", "a1331852_C6", "venetoclax_C6", "navitoclax_C6")],
       "UMAP", color_by = "merging1", facet_by = "sample_id") +
  theme_void() +
  scale_color_manual(values = pal)
ggsave("Fig5E_umap_concentrations_selected.pdf", height = 4, width = 4.5)


# Drug sensitivity heatmaps (Figure 5D)
p <- plotExprHeatmap(sce, features = "type", 
                     by = "cluster_id", k = "merging1", 
                     bars = TRUE, perc = TRUE)

counts <- data.frame(table(sample_ids(sce), cluster_ids(sce, k = "merging1"))) %>%
  tidyr::pivot_wider(names_from = "Var2", values_from = Freq)

norm_counts <- counts
norm_counts$Blasts <- as.numeric(norm_counts$Blasts/mean(norm_counts$Blasts[1:2]))
norm_counts$Lymphocytes <- as.numeric(norm_counts$Lymphocytes/mean(norm_counts$Lymphocytes[1:2]))
norm_counts$Erythrocytes <- as.numeric(norm_counts$Erythrocytes/mean(norm_counts$Erythrocytes[1:2]))
norm_counts$Other <- as.numeric(norm_counts$Other/mean(norm_counts$Other[1:2]))

mat1 <- t(norm_counts[grepl("DMSO1|a133", norm_counts$Var1),c("Blasts", "Lymphocytes", "Erythrocytes")])
colnames(mat1) <- c("DMSO", "0.375", "1.25", "3.75", "12.5", "37.5", "125", "375")
h1 <- Heatmap(mat1,
              cluster_rows = F,
              cluster_columns = F,
              column_names_rot = 0,
              column_names_centered = T,
              column_names_gp = gpar(fontsize(2)),
              col = colorRamp2(seq(0, 1, length.out = 11), rev(brewer.pal(11, "PuOr"))),
              rect_gp = gpar(col= "white", lwd = unit(0.4, "mm")),
              heatmap_legend_param = list(title = "Expression (log2)\nZ-score",
                                          title_gp = gpar(fontsize = 5),
                                          labels_gp = gpar(fontsize = 5),
                                          grid_height = unit(0.2, "cm"),
                                          grid_width = unit(2, "mm"))
)

mat2 <- t(norm_counts[grepl("DMSO1|navi", norm_counts$Var1),c("Blasts", "Lymphocytes", "Erythrocytes")])
colnames(mat2) <- c("DMSO", "1.25", "3.75", "12.5", "37.5", "125", "375", "1250")
h2 <- Heatmap(mat2,
              cluster_rows = F,
              cluster_columns = F,
              column_names_rot = 0,
              column_names_centered = T,
              column_names_gp = gpar(fontsize(2)),
              col = colorRamp2(seq(0, 1, length.out = 11), rev(brewer.pal(11, "PuOr"))),
              rect_gp = gpar(col= "white", lwd = unit(0.4, "mm")),
              heatmap_legend_param = list(title = "Expression (log2)\nZ-score",
                                          title_gp = gpar(fontsize = 5),
                                          labels_gp = gpar(fontsize = 5),
                                          grid_height = unit(0.2, "cm"),
                                          grid_width = unit(2, "mm"))
)

mat3 <- t(norm_counts[grepl("DMSO1|vene", norm_counts$Var1),c("Blasts", "Lymphocytes", "Erythrocytes")])
colnames(mat3) <- c("DMSO", "0.375", "1.25", "3.75", "12.5", "37.5", "125", "375")
h3 <- Heatmap(mat3,
              cluster_rows = F,
              cluster_columns = F,
              column_names_rot = 0,
              column_names_centered = T,
              column_names_gp = gpar(fontsize(2)),
              col = colorRamp2(seq(0, 1, length.out = 11), rev(brewer.pal(11, "PuOr"))),
              rect_gp = gpar(col= "white", lwd = unit(0.4, "mm")),
              heatmap_legend_param = list(title = "Expression (log2)\nZ-score",
                                          title_gp = gpar(fontsize = 5),
                                          labels_gp = gpar(fontsize = 5),
                                          grid_height = unit(0.2, "cm"),
                                          grid_width = unit(2, "mm"))
)


pdf("Fig5D_A1331852_heatmap.pdf", height = 1, width = 6)
h1
dev.off()

pdf("Fig5D_Navitoclax_heatmap.pdf", height = 1, width = 6)
h2
dev.off()

pdf("Fig5D_Venetoclax_heatmap.pdf", height = 1, width = 6)
h3
dev.off()

