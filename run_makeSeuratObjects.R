
folders        <- list.dirs("data/scRNAseq/", recursive = T)[-c(1,17)]
scrnaseq_files <- lapply(folders, function(x){message(getSeuratName(x)); Read10X(data.dir = x) %>% CreateSeuratObject(project = getSeuratName(x), min.cells = 3, min.features = 200)})
bcl_seurat     <- scrnaseq_files[[1]]

## Basic QC
dir.create("results/qc/", showWarnings = F)
dir.create("results/qc/before_1/", showWarnings = F)
dir.create("results/qc/after_1/", showWarnings = F)

bcl_seurat  <- PercentageFeatureSet(bcl_seurat, pattern = "^MT-", col.name = "percent.mt")
bcl_seurat  <- PercentageFeatureSet(bcl_seurat, pattern = "^RP", col.name = "percent.ribo")
bcl_seurat  <- PercentageFeatureSet(bcl_seurat, features = cycle.genes, col.name = "percent.cycle")
bcl_seurat@meta.data$barcode   <- colnames(bcl_seurat)

bcl_seurat %>% plotQC(folder = "results/qc/before_1/")
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

## Get scVI
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
ggsave("results/qc/after_1/scatter_clustering.png", width = 5, height = 4)

Idents(bcl_seurat) <- bcl_seurat$RNA_snn_res.0.7 %>% as.character() %>% as.numeric() %>% as.factor()
bcl_seurat$cluster <- Idents(bcl_seurat)
saveRDS(bcl_seurat, "results/bcl_seurat.rds")

all_markers <- FindAllMarkers(bcl_seurat, test = "t")
fwrite(all_markers, "results/de.txt", sep = "\t", quote = F, row.names = F)







##### === integrate with previous hemap samples

folders        <- list.dirs("data/scRNAseq/", recursive = T)[-c(1,17)]
scrnaseq_files <- lapply(folders, function(x){message(getSeuratName(x)); Read10X(data.dir = x) %>% CreateSeuratObject(project = getSeuratName(x), min.cells = 3, min.features = 200)})
bcl_seurat     <- scrnaseq_files[[1]]

## Basic QC
dir.create("results/qc/", showWarnings = F)
dir.create("results/qc/before_1/", showWarnings = F)
dir.create("results/qc/after_1/", showWarnings = F)

bcl_seurat  <- PercentageFeatureSet(bcl_seurat, pattern = "^MT-", col.name = "percent.mt")
bcl_seurat  <- PercentageFeatureSet(bcl_seurat, pattern = "^RP", col.name = "percent.ribo")
bcl_seurat  <- PercentageFeatureSet(bcl_seurat, features = cycle.genes, col.name = "percent.cycle")
bcl_seurat@meta.data$barcode   <- colnames(bcl_seurat)

bcl_seurat %>% plotQC(folder = "results/qc/before_1/")
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

## Get scVI
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
ggsave("results/qc/after_1/scatter_clustering.png", width = 5, height = 4)

Idents(bcl_seurat) <- bcl_seurat$RNA_snn_res.0.7 %>% as.character() %>% as.numeric() %>% as.factor()
bcl_seurat$cluster <- Idents(bcl_seurat)
saveRDS(bcl_seurat, "results/bcl_seurat.rds")

all_markers <- FindAllMarkers(bcl_seurat, test = "t")
fwrite(all_markers, "results/de.txt", sep = "\t", quote = F, row.names = F)
