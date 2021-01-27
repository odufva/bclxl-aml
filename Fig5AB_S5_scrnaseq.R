
# Analyze scRNA-seq data from patient AML-1 with megakaryocytic differentiation (Figures 5 and S5)

# load libraries
library(Seurat)
library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(wesanderson)
library(cowplot)
library(data.table)
library(viridis)

# load data
bcl_seurat <- readRDS("bcl_seurat.rds")

# Simplify SingleR annotations
prop.table(table(bcl_seurat$singler_blueprint_pred)) > 0.01 # list cell types present in less thatn 1% total cells
bcl_seurat$cell_type <- ifelse(grepl("T", bcl_seurat$singler_blueprint_pred), "T cells",
                               ifelse(grepl("CMP|GMP|MPP|CLP|rare|Plasma|Neutrophils", bcl_seurat$singler_blueprint_pred), "Other", bcl_seurat$singler_blueprint_pred))
bcl_seurat$cell_type <- factor(bcl_seurat$cell_type, levels = c("MEP", "T cells", "Megakaryocytes", "Erythrocytes", "NK cells", "Monocytes", "Other"))

# Cell type colors
pal <- c(wes_palette("Cavalcanti1", 5), wes_palette("GrandBudapest1", 4), wes_palette("Royal2", 4), "grey70")
names(pal) <- c("CMP", "GMP", "HSC", "Monocytes", "MPP", "Neutrophils", "MEP", "Megakaryocytes", "Erythrocytes",
                "T cells", "NK cells", "Plasma cells", "CLP", "Other")

# SingleR UMAP (Figure 5A)
DimPlot(bcl_seurat, reduction = "umap", group.by = "cell_type", label = T, repel = T) +
  labs(x = "UMAP 1", y = "UMAP 2") +
  scale_color_manual(values = pal) +
  ggtitle("Cell types") +
  theme(plot.title = element_text(face = "plain"))
ggsave("Fig5A_umap_singler.pdf", width = 5.75, height = 4)


# Cell type fractions bar plot (Figure 5B)
bcl_prop <- bcl_seurat@meta.data %>%
  group_by(cell_type) %>%
  summarise(n = n()) %>%
  mutate(prop = 100*(n / sum(n))) %>%
  mutate(sample = "sample")

bcl_prop %>% mutate(cell_type = factor(cell_type, levels = c("MEP", "Megakaryocytes", "Erythrocytes", "Monocytes", "T cells", "NK cells", "Other"))) %>%
  ggplot(aes(x = cell_type, y = prop, fill = cell_type, label = cell_type)) +
  geom_col() +
  theme_cowplot() +
  labs(x = "") +
  ylab("%") +
  guides(fill = F) +
  scale_y_continuous(expand = c(0,0)) +
  scale_fill_manual(values = pal) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank())
ggsave("Fig5B_barplot_celltypes.pdf", width = 3, height = 1.25)

# Dot plot (Figure 5B)
genes = c("ITGA2B", "GP1BA", "GP1BB", "ITGB3", "TFRC",
          "NFE2", "GATA1", "RUNX1", "TAL1", "FLI1",
          "CD34", "KIT", "CD38",
          "BCL2L1", "BCL2", "MCL1")

bcl_seurat$cell_type <- factor(bcl_seurat$cell_type, levels = c("MEP", "Megakaryocytes", "Erythrocytes", "Monocytes", "T cells", "NK cells", "Other"))
Idents(bcl_seurat) <- bcl_seurat$cell_type

DotPlot(bcl_seurat, features = rev(genes)) +
  coord_flip() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text.y = element_text(face = "italic"),
        axis.text = element_text(size = 10)) +
  scale_color_distiller(palette = "Reds", values = seq(0, 1, length.out = 9), direction = 1) +
  xlab("") +
  ylab("")

ggsave("Fig5B_dotplot.pdf", height = 4.5, width = 4.5)

# Clusters and DEG
markers <- fread("all_markers_1e3.txt", data.table = F)
topmarkers <- markers %>% group_by(cluster) %>% top_n(n = 30, wt = avg_logFC)
bcl_seurat$cluster <- factor(bcl_seurat$cluster, labels = unique(markers$cluster))

clusters_order <- c("0 MEP", "2 MEP", "4 MEP", "6 MEP",
                    "8 MEP", "10 MEP", "13 MEP", "16 MEP",
                    "9 Megakaryocytes", "12 MEP/megakaryocytes",
                    "3 Erythrocytes", "15 Erythrocytes", "14 Monocytes",
                    "1 CD4", "5 CD4", "7 CD8", "11 T/NK", "17 cycling")

bcl_seurat$cluster <- factor(bcl_seurat$cluster, levels = clusters_order)
Idents(bcl_seurat) <- "cluster"

colors_cluster <- c("#ED2891", "#B2509E", "#D49DC7", "#C1A72F", "#E8C51D",
                    "#F9ED32", "#104A7F", "#9EDDF9", "#007EB5", "#CACCDB",
                    "#6E7BA2", "#DAF1FC", "#00AEEF", "#F6B667", "#D97D25",
                    "#FBE3C7", "#F89420", "#97D1A9")

# Cluster UMAP
DimPlot(bcl_seurat, reduction = "umap", group.by = "cluster", label = T, repel = T) +
  labs(x = "UMAP 1", y = "UMAP 2") +
  scale_color_manual(values = colors_cluster) +
  ggtitle("Clusters") +
  theme(plot.title = element_text(face = "plain"))

ggsave("FigS5A_umap_clusters.pdf", width = 9, height = 7)

# Cluster DEG heatmap
markers$cluster <- factor(markers$cluster, levels = clusters_order)
topmarkers <- markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC) %>% arrange(cluster)

pdf("FigS5B_deg_heatmap.pdf", height = 17, width = 15)
DoHeatmap(subset(bcl_seurat, downsample = 50),
          features = topmarkers$gene,
          raster = F,
          angle = 90,
          group.colors = colors_cluster) +
  scale_fill_distiller(palette = "RdBu", values = seq(0, 1, length.out = 9), direction = -1)
dev.off()

# SingleR cell type DEG
Idents(bcl_seurat) <- "cell_type"
mep_mega_markers <- FindMarkers(bcl_seurat, test.use = "t", verbose = T, max.cells.per.ident = 1e3,
                                only.pos = F, return.thresh = 0.05, logfc.threshold = 0.1,
                                ident.1 = c("MEP", "Megakaryocytes"))
mep_mega_markers$gene <- rownames(mep_mega_markers)

fwrite(mep_mega_markers, "TableS5_mep_mega_markers_1e3.txt", sep = "\t", quote = F, row.names = F)

# Cluster and SingleR identities and UMAP coordinates
df <- data.frame(Barcode = gsub("-1", "", bcl_seurat$barcode),
                 Cluster = bcl_seurat$cluster,
                 SingleR = bcl_seurat$singler_blueprint_pred,
                 SingleR_simple = bcl_seurat$cell_type)
df <- merge(df,  bcl_seurat[["umap"]]@cell.embeddings, by = "row.names")
df <- df[,-1]

fwrite(df, "TableS5_metadata.txt", sep = "\t", quote = F, row.names = F)




