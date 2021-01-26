
me=system("whoami", intern = TRUE)
setwd(paste0("/Users/", me, "/Dropbox/bclxl/"))

library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(gridExtra)
library(data.table)
library(Seurat)

## Set global ggplot2-themes
theme_set(theme_classic(base_size = 12))

add_guide   <- guides(colour = guide_legend(override.aes = list(size=5)))

getPalette  <- colorRampPalette(brewer.pal(9, "Set1"))
getPalette2 <- colorRampPalette(brewer.pal(8, "Set2"))
getPalette3 <- colorRampPalette(brewer.pal(9, "Set3"))
getPalette4 <- colorRampPalette(brewer.pal(9, "Pastel1"))
getPalette5 <- colorRampPalette(brewer.pal(8, "Pastel2"))


getDoublets <- function(seurat_object){

  require(scds)

  # Annotate doublet using co-expression based doublet scoring:
  sce_object <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = seurat_object@assays$RNA@counts), colData = seurat_object@meta.data)
  sce_object <- cxds(sce_object)
  sce_object <- bcds(sce_object)
  sce_object <- cxds_bcds_hybrid(sce_object)

  ## Add into Seurat
  seurat_object$cxds_doublet_score   <- SingleCellExperiment::colData(sce_object)$cxds_score
  seurat_object$bcds_doublet_score   <- SingleCellExperiment::colData(sce_object)$bcds_score

  seurat_object$hybrid_doublet_score <- SingleCellExperiment::colData(sce_object)$hybrid_score
  seurat_object$cxds_doublet_score_norm <- c(SingleCellExperiment::colData(sce_object)$cxds_score - min(SingleCellExperiment::colData(sce_object)$cxds_score)) / max(SingleCellExperiment::colData(sce_object)$cxds_score)
  seurat_object$bcds_doublet_score_norm <- c(SingleCellExperiment::colData(sce_object)$bcds_score - min(SingleCellExperiment::colData(sce_object)$bcds_score)) / max(SingleCellExperiment::colData(sce_object)$bcds_score)
  return(seurat_object)

}

getClustering <- function(seurat_object){

  ## Clustering
  res           <- c(seq(0.1, 1, 0.1), seq(1.2, 2, 0.2), 2.5, 3)
  seurat_object <- FindNeighbors(seurat_object, reduction = "pca", dims = c(1:ncol(seurat_object@reductions$pca@cell.embeddings)))
  seurat_object <- FindClusters(object = seurat_object, resolution = res, verbose = F)
  return(seurat_object)

}

extractPair1 <- function(str1){
  # strsplit(str1, "[_]")[[1]][1]
  sub("\\_.*", "", str1)
}

extractPair2 <- function(str1){
  # strsplit(str1, "[_]")[[1]][1]
  sub(".*\\_", "", str1)
}

getNewClusters <- function(df, cluster){
  df %>% mutate(cluster = do.call(what = "c", extractClusterNumber(cluster)) %>% as.factor() %>% getClusterPhenotypesAANew()) %>% return()
}

getSingler <- function(seurat_object, cluster = NULL, method = NULL, sample = NULL){

  hpca.se   <- SingleR::HumanPrimaryCellAtlasData()
  blueprint <- SingleR::BlueprintEncodeData()
  ## @ params
  ## cluster = possible cluster vec, if not provided, tries to find in meta.data$cluster
  ## method = if "cluster", then performs preds based on clusters, not cells
  ## sample = to subsample or not

  if(!is.null(sample)){

    set.seed(123)
    seurat_object <- subset(seurat_object, cells = colnames(seurat_object)[sample(1:ncol(seurat_object), sample)])

  }

  sce       <- as.SingleCellExperiment(seurat_object)

  ## Predictions
  if(is.null(method)){
    pred.hca <- SingleR::SingleR(test = sce, ref = hpca.se, assay.type.test = 1,   labels = hpca.se$label.fine)
    pred.blu <- SingleR::SingleR(test = sce, ref = blueprint, assay.type.test = 1, labels = blueprint$label.fine)

    if(is.null(sample)){
      seurat_object$singler_hpca_pred      <- pred.hca$first.labels
      seurat_object$singler_blueprint_pred <- pred.blu$first.labels
      return(seurat_object)
    }

    else{
      df <- data.frame(barcode = rownames(pred.hca), cluster = seurat_object$cluster, singler_hpca_pred = pred.hca$labels, singler_blueprint_pred = pred.blu$labels)
      return(df)
    }

  }


  if(method == "cluster"){
    if(is.null(cluster)){
      cluster=seurat_object$cluster
    }
    pred.hca <- SingleR::SingleR(test = sce, ref = hpca.se, assay.type.test = 1,   labels = hpca.se$label.fine, method = "cluster", clusters = cluster)
    pred.blu <- SingleR::SingleR(test = sce, ref = blueprint, assay.type.test = 1, labels = blueprint$label.fine, method = "cluster", clusters = cluster)
    df <- data.frame(cluster = rownames(pred.hca), singler_hpca_pred = pred.hca$labels, singler_blueprint_pred = pred.blu$labels)
    return(df)
  }
}

plotClustering <- function(seurat_object){

  res <- c(seq(0.1, 1, 0.1), seq(1.2, 2, 0.2), 2.5, 3)
  clustering_columns <- grep("res", colnames(seurat_object@meta.data), value = T)
  clustering_columns <- clustering_columns[order(substr(clustering_columns, 13, nchar(clustering_columns)) %>% as.numeric())]

  q <- NULL; i <- 1

  for(clustering_column in clustering_columns){
    q[[i]] <- seurat_object@meta.data[,clustering_column] %>% levels %>% length
    i <- i + 1
  }

  data.frame(resolution = res, nClusters = do.call(q, what="c")) %>%
    ggplot(aes((resolution),nClusters), label = nClusters) + geom_point(shape = 21) + theme_bw()

}

fixSeurat <- function(seurat_object){

  ## Fix meta data if it brokes

  meta.data           <- seurat_object@meta.data
  count.data          <- seurat_object@assays$RNA@counts
  scale.data          <- seurat_object@assays$RNA@scale.data
  # hvg                 <- VariableFeatures(seurat_object)

  # pca_dimred          <- seurat_object[["pca"]]
  # umap_dimred         <- seurat_object[["umap"]]
  latent_dimred       <- seurat_object[["latent"]]
  latent_umap_dimred  <- seurat_object[["latent_umap"]]

  rownames(meta.data) <- meta.data$barcode

  old_idents <- Idents(seurat_object)
  new_seurat <- CreateSeuratObject(counts = count.data)

  new_seurat@meta.data             <- meta.data
  new_seurat@assays$RNA@counts     <- count.data
  new_seurat@assays$RNA@scale.data <- scale.data
  # VariableFeatures(seurat_object)  <- hvg

  # new_seurat[["pca"]]              <- pca_dimred
  # new_seurat[["umap"]]             <- umap_dimred
  new_seurat[["latent"]]           <- latent_dimred
  new_seurat[["latent_umap"]]      <- latent_umap_dimred
  Idents(new_seurat) <- old_idents
  return(new_seurat)

}

getLatentClustering <- function(seurat_object){

  ## Clustering
  res        <- c(seq(0.1, 1, 0.1), seq(1.2, 2, 0.2), 2.5, 3)
  seurat_object <- FindNeighbors(seurat_object, reduction = "latent", dims = c(1:ncol(seurat_object@reductions$latent@cell.embeddings)))
  seurat_object <- FindClusters(object = seurat_object, resolution = res, verbose = F)
  return(seurat_object)

}

putLatentsSeurat <- function(seurat_object, latent){

  latent_umap <- uwot::umap(latent) %>% as.data.frame() %>% dplyr::rename(UMAP1=V1, UMAP2=V2)

  latent      <- as.matrix(latent)
  latent_umap <- as.matrix(latent_umap)

  rownames(latent)      <- colnames(seurat_object)
  rownames(latent_umap) <- colnames(seurat_object)

  latent_dim_red            <- CreateDimReducObject(key = "latent", embeddings = as.matrix(x = latent))
  latent_umap_dim_red       <- CreateDimReducObject(key = "latent_umap", embeddings = as.matrix(x = latent_umap))

  seurat_object[['latent']]      <- latent_dim_red
  seurat_object[['latent_umap']] <- latent_umap_dim_red
  return(seurat_object)
}

getScviInput <- function(seurat_object, folder){

  dir.create(folder, showWarnings = F)
  genes_to_keep <- list()
  i <- 1
  Idents(seurat_object) <- seurat_object$orig.ident

  for(patient in unique(seurat_object$orig.ident)){

    message(patient)
    seurat_temp <- subset(seurat_object, idents = patient)
    counts_temp <- seurat_temp@assays$RNA@counts %>% as.data.frame

    genes_to_keep[[i]] <- rownames(counts_temp)
    i <- i + 1

    counts_temp <- seurat_object@assays$RNA@counts[ ,seurat_object$orig.ident == patient] %>% as.data.frame
    counts_temp <- counts_temp[!rownames(counts_temp) %in% clonality_genes, ]
    data.table::fwrite(counts_temp, paste0(folder, patient, ".csv"), sep = ",", quote = F, row.names = T, col.names = T)

  }
}

preprocessSeurat <- function(orig_object, cells.to.use){

  ## Subset object
  object <- subset(orig_object, cells = cells.to.use)

  # orig_object@meta.data$barcode
  temp_meta <- orig_object@meta.data[as.character(orig_object@meta.data$barcode) %in% cells.to.use, ]
  temp_meta <- temp_meta[match(colnames(object), temp_meta$barcode), ]
  temp_meta$barcode == colnames(object)
  object@meta.data <- temp_meta

  ## Normalize and find HVGs
  object  <- NormalizeData(object, normalization.method = "LogNormalize", scale.factor = 10000)
  object  <- FindVariableFeatures(object, selection.method = "vst", nfeatures = 2000, clip.max = 10)

  ## Remove clonality genes
  hvg     <- VariableFeatures(object)
  too_hvg <- HVFInfo(object = object) %>% add_rownames(var = "gene") %>% filter(variance.standardized > 10) %>% pull("gene") %>% as.character()
  hvg     <- hvg[!hvg %in% too_hvg]
  hvg     <- hvg[!hvg %in% clonality_genes]
  hvg     <- hvg[!hvg %in% unwanted_genes]

  VariableFeatures(object) <- hvg
  # plotHVG(object, 30) #+ ylim(values = c(0,10))

  ## Scale data
  object <- ScaleData(object, features = hvg)

  ## PCA data
  object <- RunPCA(object, features = hvg, npcs = 50)
  nPCs   <- sum(object[["pca"]]@stdev > 2)
  print(paste("nPCs:", nPCs))

  ## RunUMAP does not work
  object <- RunUMAP(object, dims = 1:nPCs, learning.rate = 1)

  # Meanwhile try something hacky-ish
  # umap_df <- object[["pca"]]@cell.embeddings[,1:nPCs] %>% umapr::umap() %>% select(UMAP1:UMAP2)
  # umap_df <- CreateDimReducObject(key = "umap", embeddings = as.matrix(x = umap_df))
  # object[["umap"]] <- umap_df

  return(object)

}

getQC <- function(seurat_object){

  ###################

  min_mito     <- 0
  max_mito     <- 15

  min_ribo     <- 5
  max_ribo     <- 50

  min_features <- 300
  max_features <- 5e3

  min_counts   <- 1e3
  max_counts   <- 30e3


  ###################

  seurat_object@meta.data$barcode <- colnames(seurat_object)

  ## In total, we remove with the following conditions:
  qc_df <- seurat_object@meta.data %>% as.data.frame()

  percent_mito_outlier <- qc_df %>% dplyr::filter(percent.mt   > max_mito     | percent.mt   < min_mito)     %>% pull(barcode) %>% as.character()
  percent_ribo_outlier <- qc_df %>% dplyr::filter(percent.ribo > max_ribo     | percent.ribo < min_ribo)     %>% pull(barcode) %>% as.character()
  features_outlier     <- qc_df %>% dplyr::filter(nFeature_RNA < min_features | nFeature_RNA > max_features) %>% pull(barcode) %>% as.character()
  umis_outlier         <- qc_df %>% dplyr::filter(nCount_RNA   > max_counts   | nCount_RNA   < min_counts)   %>% pull(barcode) %>% as.character()

  outlier_cells        <- c(percent_mito_outlier,
                            percent_ribo_outlier,
                            features_outlier,
                            umis_outlier)

  reason               <- c(rep("percent_mito_outlier", length(percent_mito_outlier)),
                            rep("percent_ribo_outlier", length(percent_ribo_outlier)),
                            rep("features_outlier",     length(features_outlier)),
                            rep("umis_outlier",         length(umis_outlier)))

  outlier_df <- data.frame(barcode = outlier_cells, reason = reason) %>% dplyr::mutate(from = extractName(barcode)) #, 1, 10))

  ## Remove the cells from Seurat-object and save a new seurat-object
  cells.to.use  <- colnames(seurat_object)[!colnames(seurat_object) %in% outlier_df$barcode]
  seurat_object <- subset(seurat_object, cells = cells.to.use)
  return(seurat_object)

}

getLatentClustering <- function(seurat_object){

  ## Clustering
  res        <- c(seq(0.1, 1, 0.1), seq(1.2, 2, 0.2), 2.5, 3)
  seurat_object <- FindNeighbors(seurat_object, reduction = "latent", dims = c(1:ncol(seurat_object@reductions$latent@cell.embeddings)))
  seurat_object <- FindClusters(object = seurat_object, resolution = res, verbose = F)
  return(seurat_object)

}

extractClusterNumber <- function(strs){

  p <- NULL
  i <- 1
  for(str1 in strs){
    p[[i]] <- strsplit(str1, "[ ]")[[1]][1]
    i <- i + 1
  }

  return(p)

}

facets_nice <- theme(strip.background = element_rect(fill="grey96"), strip.text = element_text(colour = 'black'))

extractName = function(str1){
  # strsplit(str1, "[_]")[[1]][1]
  sub("\\_.*", "", str1)
}

extractFileName = function(str1){
  # strsplit(str1, "[_]")[[1]][1]
  sub(".*\\/", "", str1)
}

extractSeuratName <- function(str1){

  str1 <- substr(str1, 1, nchar(str1) - 27)
  extractFileName(str1)
}

extractTimepoint <- function(strs){

  strs2 <-NULL
  i <- 1
  for(str1 in strs){
    strs2[[i]] <- strsplit(str1, "[_]")[[1]][2]
    i <- i + 1
  }

  # return(strs2)
  return(factor(strs2, levels = c("311019", "141119", "021219", "301219", "130120", "280120")))

}

getClonalityGenes <- function(object){

  clonality_genes <- c(grep("^TRAV", rownames(object), value = T), grep("^TRBV", rownames(object), value = T),
                       grep("^TRGV", rownames(object), value = T), grep("^TRDV", rownames(object), value = T),
                       grep("^IGLV", rownames(object), value = T), grep("^IGLC", rownames(object), value = T),
                       grep("^IGLL", rownames(object), value = T), grep("^IGKV", rownames(object), value = T),
                       grep("^IGHV", rownames(object), value = T), grep("^IGKC", rownames(object), value = T),
                       grep("^IGH", rownames(object), value = T),  grep("^IGK", rownames(object), value = T))

}

getUnwantedGenes <- function(object){

  unwanted_variation <- c(grep("^LINC", rownames(object), value = T), grep("^AC", rownames(object), value = T),
                          grep("^AL", rownames(object), value = T),
                          grep("^MT-", rownames(object), value = T), grep("^RP", rownames(object), value = T))

}

fixSeurat <- function(seurat_object){

  ## Fix meta data if it brokes

  meta.data           <- seurat_object@meta.data
  count.data          <- seurat_object@assays$RNA@counts
  scale.data          <- seurat_object@assays$RNA@scale.data
  # hvg                 <- VariableFeatures(seurat_object)

  # pca_dimred          <- seurat_object[["pca"]]
  # umap_dimred         <- seurat_object[["umap"]]
  latent_dimred       <- seurat_object[["latent"]]
  latent_umap_dimred  <- seurat_object[["latent_umap"]]

  rownames(meta.data) <- meta.data$barcode

  old_idents <- Idents(seurat_object)
  new_seurat <- CreateSeuratObject(counts = count.data)

  new_seurat@meta.data             <- meta.data
  new_seurat@assays$RNA@counts     <- count.data
  new_seurat@assays$RNA@scale.data <- scale.data
  # VariableFeatures(seurat_object)  <- hvg

  # new_seurat[["pca"]]              <- pca_dimred
  # new_seurat[["umap"]]             <- umap_dimred
  new_seurat[["latent"]]           <- latent_dimred
  new_seurat[["latent_umap"]]      <- latent_umap_dimred
  Idents(new_seurat) <- old_idents
  return(new_seurat)

}

extractCoarsePhenotype <- function(strs){

  p <- NULL
  i <- 1
  for(str1 in strs){
    p[[i]] <- strsplit(str1, "[ ]")[[1]][2]
    i <- i + 1
  }

  return(p)

}

getLatentUMAP <- function(seurat_object){

  umap_df           <- seurat_object[["latent"]]@cell.embeddings %>% uwot::umap()
  colnames(umap_df) <- c("latent_umap1", "latent_umap2")
  rownames(umap_df) <- colnames(seurat_object)
  umap_df           <- CreateDimReducObject(key = "latent_umap", embeddings = as.matrix(x = umap_df))
  seurat_object[['latent_umap']] <- umap_df

  return(seurat_object)

}
