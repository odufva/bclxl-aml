
# Analyze scRNA-seq data from patient AML-1 with megakaryocytic differentiation (Figure 6)

# load libraries
library(Seurat)
library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(wesanderson)

# load data
blast_seurat <- readRDS("blastpre_seurat.rds")

# UMAP with SingleR cell type annotations colored
pal <- c(wes_palette("Cavalcanti1", 5), wes_palette("GrandBudapest1", 4), wes_palette("Royal2", 5))
names(pal) <- c("CMP", "GMP", "HSC", "Monocytes", "MPP", "Neutrophils", "MEP", "Megakaryocytes", "Erythrocytes",
                "T cells", "NK cells", "Plasma cells", "CLP", "Other")

DimPlot(blast_seurat, reduction = "umap", group.by = "singler_blueprint_pred", label = T, repel = T) +
  labs(x = "UMAP 1", y = "UMAP 2") +
  scale_color_manual(values = pal) +
  ggtitle("Cell types") +
  theme(plot.title = element_text(face = "plain"))
ggsave("Fig6A_umap_singler.pdf", width = 5.75, height = 4)

# UMAP with samples colored
pal <- wes_palette("Darjeeling1", 5)[c(5,4,3,1,1)]
DimPlot(blast_seurat, reduction = "umap", group.by = "orig.ident", label = T, repel = T) +
  labs(x = "UMAP 1", y = "UMAP 2") +
  scale_color_manual(values = pal) +
  ggtitle("Patients") +
  theme(plot.title = element_text(face = "plain"))
ggsave("Fig6B_umap_sample.pdf", width = 5.5, height = 4)

# UMAP with BCL-2 family gene expression colored
genes = c("BCL2L1", "BCL2", "MCL1")
FeaturePlot(blast_seurat, reduction = "umap", features = genes, order = T,
            ncol = 3, label = F, cols = c("gray90", "tomato"), max.cutoff = c(3,3,15)) +
  labs(x = "UMAP 1", y = "UMAP 2") &
  theme(plot.title = element_text(face = "italic")) &
  scale_color_distiller(palette = "Reds", values = seq(0, 1, length.out = 9), direction = 1)
ggsave("Fig6B_umap_features.pdf", width = 4.5*3, height = 4)

# Cluster and SingleR identities and UMAP coordinates
df <- data.frame(Barcode = gsub("-1", "", blast_seurat$barcode),
                 Cluster = blast_seurat$cluster,
                 SingleR = blast_seurat$singler_blueprint_pred)
df <- merge(df,  blast_seurat[["umap"]]@cell.embeddings, by = "row.names")
df <- df[,-1]

fwrite(df, "TableS5_metadata_blast.txt", sep = "\t", quote = F, row.names = F)


