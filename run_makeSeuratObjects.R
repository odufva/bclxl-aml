
## The index patient (9119) samples
folders        <- list.dirs("data/scRNAseq/9119/", recursive = T)[-c(1,17)]
scrnaseq_files <- lapply(folders, function(x){message(getSeuratName(x)); Read10X(data.dir = x) %>% CreateSeuratObject(project = getSeuratName(x), min.cells = 3, min.features = 200)})
bcl_seurat     <- scrnaseq_files[[1]]

## Basic QC
bcl_seurat <- bcl_seurat %>% getQC()

## Get SingleR predictions; omit predictions from cell types rare than 10 cells
bcl_seurat             <- bcl_seurat %>% getSingler()
relevant_hpca_clusters <- bcl_seurat@meta.data %>% group_by(singler_hpca_pred) %>% summarise(n=n()) %>% filter(n>=10) %>% pull(singler_hpca_pred)
relevant_blue_clusters <- bcl_seurat@meta.data %>% group_by(singler_blueprint_pred) %>% summarise(n=n()) %>% filter(n>=10) %>% pull(singler_blueprint_pred)

bcl_seurat$singler_hpca_pred      <- ifelse(bcl_seurat$singler_hpca_pred %in% relevant_hpca_clusters, bcl_seurat$singler_hpca_pred, "rare")
bcl_seurat$singler_blueprint_pred <- ifelse(bcl_seurat$singler_blueprint_pred %in% relevant_blue_clusters, bcl_seurat$singler_blueprint_pred, "rare")

## Get doublets
bcl_seurat <- bcl_seurat %>% getDoublets()
bcl_seurat <- subset(bcl_seurat, hybrid_doublet_score < 1.8)

## Get Seurat
clonality_genes <- getClonalityGenes(bcl_seurat)
unwanted_genes  <- getUnwantedGenes(bcl_seurat)

bcl_seurat <- bcl_seurat %>% preprocessSeurat(cells.to.use = colnames(bcl_seurat))
bcl_seurat <- bcl_seurat %>% getClustering()

## Make scVI input
dir.create("results/scvi/", showWarnings = F)
dir.create("results/scvi/input_files/", showWarnings = F)
bcl_seurat$orig.ident <- gsub("\\/", "\\_", bcl_seurat$orig.ident)
bcl_seurat %>% getScviInput(folder = "results/scvi/input_files/")

## Get scVI results
latents    <- fread("results/scvi/results/bcl_latent.csv")
bcl_seurat <- bcl_seurat %>% putLatentsSeurat(latent = latents)
bcl_seurat <- bcl_seurat %>% getLatentClustering() %>% fixSeurat()

## Decide on clustering
bcl_seurat %>% plotClustering()
Idents(bcl_seurat) <- bcl_seurat$RNA_snn_res.0.7 %>% as.character() %>% as.numeric() %>% as.factor()
saveRDS(bcl_seurat, "results/bcl_seurat.rds")

## Get DEGs
all_markers <- FindAllMarkers(bcl_seurat, test = "t")
fwrite(all_markers, "results/de.txt", sep = "\t", quote = F, row.names = F)











##### === integrate with previous hemap samples

folders        <- list.dirs("data/scRNAseq/", recursive = T)[-1]
scrnaseq_files <- lapply(folders, function(x){message(getSeuratName(x)); Read10X(data.dir = x) %>% CreateSeuratObject(project = getSeuratName(x), min.cells = 3, min.features = 200)})
hemap_seurat   <- merge(scrnaseq_files[[1]], scrnaseq_files[-1], add.cell.ids = getSeuratName(folders))

## Basic QC
hemap_seurat <- hemap_seurat %>% getQC()

## Get SingleR predictions; omit predictions from cell types rare than 10 cells
hemap_seurat           <- hemap_seurat %>% getSingler()
relevant_hpca_clusters <- hemap_seurat@meta.data %>% group_by(singler_hpca_pred) %>% summarise(n=n()) %>% filter(n>=10) %>% pull(singler_hpca_pred)
relevant_blue_clusters <- hemap_seurat@meta.data %>% group_by(singler_blueprint_pred) %>% summarise(n=n()) %>% filter(n>=10) %>% pull(singler_blueprint_pred)

hemap_seurat$singler_hpca_pred      <- ifelse(hemap_seurat$singler_hpca_pred %in% relevant_hpca_clusters, hemap_seurat$singler_hpca_pred, "rare")
hemap_seurat$singler_blueprint_pred <- ifelse(hemap_seurat$singler_blueprint_pred %in% relevant_blue_clusters, hemap_seurat$singler_blueprint_pred, "rare")

## Get doublets
hemap_seurat <- hemap_seurat %>% getDoublets()
hemap_seurat <- subset(hemap_seurat, hybrid_doublet_score < 1.8)

## Get Seurat
hemap_seurat <- hemap_seurat %>% preprocessSeurat(cells.to.use = colnames(hemap_seurat))

## Focus only on blasts
celltypes.to.keep <- c("CMP", "GMP", "HSC", "Megakaryocytes", "MEP", "MPP", "Monocytes")
cells.to.keep     <- hemap_seurat@meta.data %>% filter(singler_blueprint_pred %in% celltypes.to.keep) %>% filter(orig.ident != "FH_9119_R1") %>% pull(barcode)
blastpre_seurat   <- subset(hemap_seurat, cells = cells.to.keep)
blastpre_seurat   <- blastpre_seurat %>% getLatentUMAP()

# ## Get scVI input files
dir.create("results/scvi/input_files/hemap_blast_pre/", showWarnings = F)
blastpre_seurat %>% getScviInput(folder = "results/scvi/input_files/hemap_blast_pre/")

## Get scVI results
latents      <- fread("results/scvi/output/hemap_blast_pre_latent.csv")
blastpre_seurat <- blastpre_seurat %>% putLatentsSeurat(latent = latents)
blastpre_seurat <- blastpre_seurat %>% getLatentClustering() %>% fixSeurat()
blastpre_seurat %>% plotClustering()

Idents(blastpre_seurat) <- blastpre_seurat$singler_blueprint_pred
saveRDS(blastpre_seurat, "results/blastpre_seurat.rds")
