
# Plot expression of BCL2L1 with GFI1B and GATA1 knockdown in data from Replogle et al. Cell 2022 (Figure S2D) 

library(anndata)
library(dplyr)
library(ggplot2)
library(patchwork)

data <- read_h5ad("/Users/odufva/Downloads/K562_gwps_normalized_bulk_01.h5ad")

df <- as.data.frame(data)

rownames(df)[grepl("GFI1B", rownames(df))]
rownames(df)[grepl("GATA1", rownames(df))]
rownames(df)[grepl("BCL2L1", rownames(df))]
rownames(df)[grepl("HBB", rownames(df))]

df["3343_GFI1B_P1P2_ENSG00000165702", "ENSG00000165702"]
df["3343_GFI1B_P1P2_ENSG00000165702", "ENSG00000171552"]

df_gfi1b <- as.data.frame(t(df["3343_GFI1B_P1P2_ENSG00000165702",]))

df_gfi1b <- df_gfi1b %>% arrange(`3343_GFI1B_P1P2_ENSG00000165702`)
df_gfi1b$gene <- rownames(df_gfi1b)

df_gfi1b$gene <- factor(df_gfi1b$gene, levels = rownames(df_gfi1b))


p1 <- ggplot(df_gfi1b, aes(x = gene, y = `3343_GFI1B_P1P2_ENSG00000165702`)) +
  geom_bar(stat = "identity", color = "grey70") +
  geom_bar(data = df_gfi1b["ENSG00000171552",], stat = "identity", color = "darkred") +
  geom_bar(data = df_gfi1b["ENSG00000170180",], stat = "identity", color = "darkred") +
  theme_cowplot() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.title = element_text(face = "plain", hjust = 0.5)) +
  xlab("Gene") +
  ylab("Expression Z-score") +
  ggtitle("GFI1B knockdown") +
  scale_x_discrete(expand = c(0.01,0.01)) +
  geom_text_repel(data = df_gfi1b["ENSG00000171552",], aes(label = ifelse(gene == "ENSG00000171552", "BCL2L1", "")), fontface = "italic") +
  geom_text_repel(data = df_gfi1b["ENSG00000170180",], aes(label = ifelse(gene == "ENSG00000170180", "GYPA", "")), fontface = "italic")



# GATA1

df_gata1 <- as.data.frame(t(df["3291_GATA1_P1P2_ENSG00000102145",]))

df_gata1 <- df_gata1 %>% arrange(`3291_GATA1_P1P2_ENSG00000102145`)
df_gata1$gene <- rownames(df_gata1)

df_gata1$gene <- factor(df_gata1$gene, levels = rownames(df_gata1))


p2 <- ggplot(df_gata1, aes(x = gene, y = `3291_GATA1_P1P2_ENSG00000102145`)) +
  geom_bar(stat = "identity", color = "grey70") +
  geom_bar(data = df_gata1["ENSG00000171552",], stat = "identity", color = "darkred") +
  geom_bar(data = df_gata1["ENSG00000170180",], stat = "identity", color = "darkred") +
  theme_cowplot() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.title = element_text(face = "plain", hjust = 0.5)) +
  xlab("Gene") +
  ylab("Expression Z-score") +
  ggtitle("GATA1 knockdown") +
  scale_x_discrete(expand = c(0.01,0.01)) +
  geom_text_repel(data = df_gata1["ENSG00000171552",], aes(label = ifelse(gene == "ENSG00000171552", "BCL2L1", "")), fontface = "italic") +
  geom_text_repel(data = df_gata1["ENSG00000170180",], aes(label = ifelse(gene == "ENSG00000170180", "GYPA", "")), fontface = "italic")


p1 + p2
ggsave("FigS2D_GFI1B_GATA1.pdf", height = 5, width = 12)
