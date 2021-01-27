#' @param seurat.object Seurat class object or data matrix that can be converted to Seurat object. Can also be a list of data matrices (RNA+citeseq etc.), Must be named by assay type with four uppercase letters (PROT, etc) (also found from seurat object with this name). 
#' @param regress.cell.label Character vector with as many columns as seurat.object Contains batch information for each cell. Uses fastMNN/SCTtransform to batch correct the data (usually different samples or experiments).
#' @param check.pcs Plot jackstraw and elbowplot to optimize the number of PCA componens
#' @param plot.umap Umap plot for clusters and for sample ID, if batch corrected also batch clustering
#' @param resolution Resolution parameter for clustering
#' @param nFeature.min Filter cells with < nFeature.min
#' @param nFeature.max Filter cells with > nFeature.max
#' @param percent.mitoDNA Filter cells with mitochondrial genes > percent.mitoDNA
#' @param singleR Annotate each cell using singleR built in annotations. Mainly uses Encode+blueprint annotations.
#' @param cores Number of cores to use for Seurat/singleR/MNNcorrect
#' @param add.gm.score Adds geneset score (using geometric mean) to fm, must be a named list of genes. Geneset name is added to fm and seurat object.
#' @param ... Pass parameters to Seurat tools
sc.data.analysis=function(scmat, name="scRNA_object", nr.pcs=30, regress.cell.label=NULL, batch.correction.method=c("SCTtransform","MNNcorrect"), auto.order=F, normalize.input=T, check.pcs=T, plot.umap=T, resolution=0.5, nFeature.min=0, nFeature.max=10000, percent.mitoDNA=25, singleR=T, DE.analysis=T, singleR.reference=c("blueprint_encode", "HPCA"), cores=6, add.gm.score=NULL, ...){
  
  # determine computational resources:
  future::plan("multiprocess", workers = cores)
  options(future.globals.maxSize= 5e+09)
  
  # set method for batch correction
  batch.correction.method=match.arg(batch.correction.method)
  
  # list of data matrix, can be citeseq data:
  if(class(scmat)=="list"){
    cat("Input is list, splitting to data matrix by type and adding to Seurat", sep="\n\n")
    
    list.additional=lapply(2:length(scmat), function(i)CreateAssayObject(scmat[[i]]))
    names(list.additional)=names(scmat)[2:length(scmat)]
    scmat <- CreateSeuratObject(scmat[[1]])
  }else{
    list.additional=NULL
    # if class is not seurat object
    if(!class(scmat)=="Seurat"){
      cat("Input is data matrix, converting to Seurat object", sep="\n\n")
      scmat <- CreateSeuratObject(scmat)
    }else{
      cat("Input is Seurat object", sep="\n\n")
    }
  }
  DefaultAssay(object = scmat) <- "RNA"
  
  # is it a vector of pcs or single value?
  if(length(nr.pcs)==1)nr.pcs=seq(nr.pcs)
  
  # QC
  scmat=PercentageFeatureSet(scmat, pattern = "^MT-", col.name = "percent.mt")
  r1=VlnPlot(scmat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
  ggplot2::ggsave(r1, filename = paste0(name, "_QC.pdf"), width = 12)
  
  # Filtering
  nFeature_RNA <- FetchData(object = scmat, vars = "nFeature_RNA")
  pct.mt <- FetchData(object = scmat, vars = "percent.mt")
  filt=which(x = nFeature_RNA > nFeature.min & nFeature_RNA < nFeature.max & pct.mt < percent.mitoDNA)
  cat(paste("Dimension before/after filtering:", paste(paste(dim(scmat), collapse="x"), paste(dim(scmat[, filt]), collapse="x"), collapse = "/")), sep="\n\n")
  scmat=scmat[, filt]
  
  if(normalize.input){
    assay="SCT"
    scmat=SCTransform(scmat, variable.features.n = 3000) #,...)
    
    # seurat standard processing
    scmat=RunPCA(scmat, npcs = 100)
    reduction="pca"
  }else{
    assay="RNA"
    
    # assume that the data is normalized
    scmat <- FindVariableFeatures(object = scmat)
    scmat <- ScaleData(object = scmat)
    scmat=RunPCA(scmat, npcs = 100)
    reduction="pca"
  }
  
  # set default assay to use:
  DefaultAssay(object = scmat) <- assay
  
  # use SCTransform to remove batch effect, used in umap and clustering
  # else use SCT to normalize data
  if(!is.null(regress.cell.label)&batch.correction.method=="SCTtransform"){
    batch.df=data.frame("batch"=regress.cell.label[filt])
    rownames(batch.df)=colnames(scmat)
    scmat[["batch"]] = batch.df
    scmat=SCTransform(scmat,vars.to.regress="batch", variable.features.n = 3000)
    
    # seurat standard processing
    scmat=RunPCA(scmat, npcs = 100)
    
    name=paste0(name, "_", batch.correction.method)
    reduction="pca"
  }
  
  # optimize pcs number:
  if(check.pcs){
    
    # setting here original RNA to define PCAs, SCT not working yet, wait for a fix
    # also the issue here how to deal with the mnnCorrect?
    # DefaultAssay(object = scmat) <- "RNA"
    # scmat <- NormalizeData(scmat,display.progress = FALSE)
    # scmat <- FindVariableFeatures(scmat,do.plot = F,display.progress = FALSE)
    # scmat <- ScaleData(scmat)
    # scmat=RunPCA(scmat, npcs = 100)
    
    # determine how many pcs should be used:
    r3=ElbowPlot(scmat, ndims = 100, reduction = "pca")
    ggplot2::ggsave(r3, filename = paste0(name, "_elbow.pdf"))
    
    # find optimal number of PCs to use:
    scmat <- JackStraw(scmat, num.replicate = 100, dims = 100, reduction = "pca")
    scmat <- ScoreJackStraw(scmat, dims = 1:100, reduction = "pca")
    r4=JackStrawPlot(scmat, dims = 1:100)
    ggplot2::ggsave(r4, filename = paste0(name, "_Jackstraw.pdf"), width = 12)
    
    cat("PCs checking done", sep="\n\n")
  }
  
  # use MNN to remove batch effect, used in umap and clustering:
  if(!is.null(regress.cell.label)&batch.correction.method=="MNNcorrect"){
    library(BiocParallel)
    library(batchelor)
    
    cat("FastMNN batch correction...", sep="\n\n")
    
    if(is.factor(regress.cell.label)&!auto.order){
      order=levels(regress.cell.label)
      cat("FastMNN batch correction order:", sep="\n\n")
      cat(order, sep="-->")
      cat("", sep="\n\n")
    }
    
    if(!is.factor(regress.cell.label)&!auto.order){
      order=unique(regress.cell.label)
      cat("FastMNN batch correction order:", sep="\n\n")
      cat(order, sep="-->")
      cat("", sep="\n\n")
    }
    
    if(auto.order){
      order=unique(regress.cell.label)
      cat("FastMNN batch correction order automatic", sep="\n\n")
    }
    
    batch.df=data.frame("batch"=regress.cell.label[filt])
    rownames(batch.df)=colnames(scmat)
    scmat[["batch"]] = batch.df
    
    if(normalize.input){
    # http://htmlpreview.github.io/?https://github.com/satijalab/seurat-wrappers/blob/master/docs/fast_mnn.html
    scmat <- NormalizeData(scmat)
    scmat <- FindVariableFeatures(scmat)
    }
    
    object.list = SplitObject(scmat, split.by = "batch")[order(order)]
    
    scmat <- SeuratWrappers::RunFastMNN(object.list, reduction.name = "mnnCorrect", auto.order=auto.order, BPPARAM=MulticoreParam(cores))
    
    # set defaults 
    reduction="mnnCorrect"
    cat("FastMNN batch correction... done", sep="\n\n")
    
  }
  
  # Basic seurat processing with umap:
  scmat=RunUMAP(scmat, dims = nr.pcs, reduction = reduction)
  scmat=FindNeighbors(scmat, dims = nr.pcs, reduction = reduction)
  scmat=FindClusters(scmat, resolution = resolution)
  
  # add celltype automatically
  if(singleR){
    cat("Annotating with singleR", sep="\n\n")
    
    if(!file.exists(paste0(name, "_singler_object.Rdata"))){
      library("SingleR")
      
      counts <- GetAssayData(scmat, assay = assay, slot = "counts")
      
      if(singleR.reference=="blueprint_encode"){
        if(dim(counts)[2]<100000){
          singler = CreateSinglerObject(counts, annot = NULL, clusters=Idents(scmat), ref.list = list(blueprint_encode), project.name=name, fine.tune = T,do.signatures = T, numCores=cores)
        }else{
          cat("SingleR in batches", sep="\n\n")
          singler = CreateBigSingleRObject.custom(counts, annot = NULL, clusters=Idents(scmat),  N=30000, ref.list = list(blueprint_encode), xy = Embeddings(scmat[["umap"]]), project.name=name, fine.tune = T, do.signatures = T, numCores=cores)
        }
      }
      
      if(singleR.reference=="HPCA"){
        if(dim(counts)[2]<100000){
          singler = CreateSinglerObject(counts, annot = NULL, clusters=Idents(scmat), ref.list = list(hpca), project.name=name, fine.tune = T,do.signatures = T, numCores=cores)
        }else{
          cat("SingleR in batches", sep="\n\n")
          singler = CreateBigSingleRObject.custom(counts, annot = NULL, clusters=Idents(scmat), N=30000, ref.list = list(hpca), project.name=name, xy = Embeddings(scmat[["umap"]]), fine.tune = T,do.signatures = T, numCores=cores)
        }
      }
      
      # add original identifiers
      singler$meta.data$orig.ident = scmat@meta.data$orig.ident # the original identities, if not supplied in 'annot'
      
      ## if using Seurat v3.0 and over use:
      singler$meta.data$xy = scmat@reductions$umap@cell.embeddings # the tSNE coordinates
      singler$meta.data$clusters = scmat@active.ident # the Seurat clusters (if 'clusters' not provided)
      
      # add annot from hpca
      # singler2 = CreateSinglerObject(counts, annot = NULL, clusters=Idents(scmat), ref.list = list(hpca), project.name=name, fine.tune = T,do.signatures = T, numCores=cores)
      # singler$singler[[1]]$SingleR.single$labels[singler2$singler[[1]]$SingleR.single$labels%in%"Pre-B_cell_CD34-"]="pre-B-cells"
      # singler$singler[[1]]$SingleR.single$labels[singler2$singler[[1]]$SingleR.single$labels%in%"Pro-B_cell_CD34+"]="pro-B-cells"
      # singler$singler[[1]]$SingleR.clusters$labels[singler2$singler[[1]]$SingleR.clusters$labels%in%"Pre-B_cell_CD34-"]="pre-B-cells"
      # singler$singler[[1]]$SingleR.clusters$labels[singler2$singler[[1]]$SingleR.clusters$labels%in%"Pro-B_cell_CD34+"]="pro-B-cells"
      
      # these names are not intuitive, replace them, mentioned that these are actually naive http://comphealth.ucsf.edu/SingleR/SupplementaryInformation2.html:
      singler$singler[[1]]$SingleR.single$labels[singler$singler[[1]]$SingleR.single$labels%in%"CD4+ T-cells"]="CD4+ Tn"
      singler$singler[[1]]$SingleR.single$labels[singler$singler[[1]]$SingleR.single$labels%in%"CD8+ T-cells"]="CD8+ Tn"
      singler$singler[[1]]$SingleR.clusters$labels[singler$singler[[1]]$SingleR.clusters$labels%in%"CD4+ T-cells"]="CD4+ Tn"
      singler$singler[[1]]$SingleR.clusters$labels[singler$singler[[1]]$SingleR.clusters$labels%in%"CD8+ T-cells"]="CD8+ Tn"
      
      save(singler, file=paste0(name, "_singler_object.Rdata"))
      
    }else{
      library("SingleR")
      load(paste0(name, "_singler_object.Rdata"))
      cat(paste0("Loaded existing SingleR object (", paste0(name, "_singler_object.Rdata"), "), remove/rename it if you want to re-compute."), sep="\n\n")
      
    }
    # can be added as cluster id
    cluster=singler$singler[[1]]$SingleR.clusters$labels
    cluster.main=singler$singler[[1]]$SingleR.clusters.main$labels
    
    # order based on cell type and add celltype to cluster:
    lineage=c("HSC","MPP","CLP","CMP","GMP","MEP","Monocytes","DC","Macrophages","Macrophages M1","Macrophages M2","CD4+ Tn","CD4+ T-cells", "CD4+ Tcm", "CD4+ Tem","CD8+ Tn","CD8+ T-cells","CD8+ Tcm","CD8+ Tem","NK cells","Tregs","naive B-cells","Memory B-cells","Class-switched memory B-cells","Plasma cells","Endothelial cells","Neutrophils","Eosinophils","Fibroblasts","Smooth muscle","Erythrocytes","Megakaryocytes")
    lineage=unique(c(lineage,blueprint_encode$types, hpca$types))
    lineage=lineage[lineage%in%singler$singler[[1]]$SingleR.single$labels]
    
    scmat[["SingleR.label.main"]]=singler$singler[[1]]$SingleR.single.main$labels
    scmat[["SingleR.label.cluster"]]=paste(singler$singler[[1]]$SingleR.single$labels, Idents(scmat))
    scmat[["SingleR.label"]]=factor(singler$singler[[1]]$SingleR.single$labels, levels=lineage[lineage%in%unique(singler$singler[[1]]$SingleR.single$labels)])
    
    identityVector.samples=as.character(Idents(scmat))
    clusters.samples=Idents(scmat)
    
    for(j in seq(cluster)){
      identityVector.samples[clusters.samples%in%rownames(cluster)[j]]=paste0(cluster[j], ":", rownames(cluster)[j])
    }
    
    # order and cluster identity;
    cluster=cluster[order(match(cluster[,1],lineage)),,drop=F]
    cat(paste(levels(Idents(scmat)), collapse=","), paste(cluster, collapse=","), sep="\t\t")   
    
    scmat[["SingleR.cluster"]]=gsub(":.*.", "", Idents(scmat)) # just the "cluster label" per cell
    
    # order seurat clusters too:
    Idents(scmat)=factor(identityVector.samples, levels=paste0(cluster, ":", rownames(cluster)))
    
    r4=DimPlot(object = scmat, group.by = "SingleR.label", reduction = 'umap' , label = TRUE,  pt.size = 0.5, repel = T)  + NoLegend() + NoAxes()
    ggplot2::ggsave(r4, filename = paste0(name, "_UMAP_singleR.pdf"), width = 6, height = 6)
    
    cat("\nSingleR done", sep="\n\n")
    
    
    # new singleR:
    # counts <- GetAssayData(scmat, assay = "RNA", slot = "counts")
    # 
    # library(SingleR)
    # library(scater)
    # 
    # ref=NovershternHematopoieticData()	
    # 
    # common <- intersect(rownames(counts), rownames(ref))
    # ref <- ref[common,]
    # counts <- counts[common,]
    # counts.log <- logNormCounts(counts)
    # 
    # singler <- SingleR(test = counts.log, ref = ref, labels = ref$label.main, numCores=cores)
    # 
    # # order based on cell type and add celltype to cluster:
    # lineage=c("HSC","MPP","CLP","CMP","GMP","MEP","Monocytes","DC","Macrophages","Macrophages M1","Macrophages M2","CD4+ Tn","CD4+ T-cells", "CD4+ Tcm", "CD4+ Tem","CD8+ Tn","CD8+ T-cells","CD8+ Tcm","CD8+ Tem","NK cells","Tregs","naive B-cells","Memory B-cells","Class-switched memory B-cells","Plasma cells","Endothelial cells","Neutrophils","Eosinophils","Fibroblasts","Smooth muscle","Erythrocytes","Megakaryocytes")
    # lineage=unique(c(lineage,ref$label.main))
    # 
    # scmat[["SingleR.label.main"]]=singler$first.labels
    # scmat[["SingleR.label"]]=factor(singler$labels, levels=lineage[lineage%in%unique(singler$labels)])
    # 
    # r4=DimPlot(object = scmat, group.by = "SingleR.label", reduction = 'umap' , label = TRUE,  pt.size = 0.5, repel = T)  + NoLegend() + NoAxes()
    # ggplot2::ggsave(r4, filename = paste0(name, "_UMAP_singleR.pdf"), width = 6, height = 6)
    # 
    # cat("\nSingleR done", sep="\n\n")
    
  }
  
  # plot clusters and samples
  if(plot.umap){
    if(!is.null(regress.cell.label)){
      r6=DimPlot(scmat, group.by = c("batch"), ncol = 2, label = T, repel = T)  + NoLegend() + NoAxes()
      ggplot2::ggsave(r6, filename = paste0(name, "_UMAP_batch.pdf"), width = 6, height = 6)
      r6=DimPlot(scmat, group.by = c("ident"), ncol = 2, label = T, repel = T)  + NoLegend() + NoAxes()
      ggplot2::ggsave(r6, filename = paste0(name, "_UMAP_clusters.pdf"), width = 6, height = 6)
    }else{
      r6=DimPlot(scmat, group.by = c("ident"), ncol = 1, label = T, repel = T)  + NoLegend() + NoAxes()
      ggplot2::ggsave(r6, filename = paste0(name, "_UMAP_clusters.pdf"), width = 6, height = 6)
    }
  }
  
  # add immunoscores to FM format data matrix
  fm.f <- GetAssayData(scmat)
  
  # add gm based score:
  add.scores=list(HLAIScore=c("B2M", "HLA-A", "HLA-B","HLA-C"), HLAIIScore=c("HLA-DMA","HLA-DMB","HLA-DPA1","HLA-DPB1","HLA-DRA","HLA-DRB1"), CytolyticScore=c("GZMA", "PRF1", "GNLY", "GZMH", "GZMM"))
  if(!is.null(add.gm.score))add.scores=append(add.scores, add.gm.score)
  
  gm.objects=do.call(rbind, lapply(seq(add.scores), function(i){
    dat3=fm.f[rownames(fm.f)%in%add.scores[[i]],]
    gm=t(apply(dat3, 2, gm_mean)) # done to normalized values
    rownames(gm)=names(add.scores)[i]
    return(gm)
  }))
  
  # also add to seurat object:
  for(i in seq(add.scores)){
    scmat[[names(add.scores)[i]]] <- gm.objects[i,]
  }
  
  # gene expression and scores to fm:
  rownames(fm.f)=paste0("N:GEXP:", rownames(fm.f))
  rownames(gm.objects)=paste0("N:SAMP:", rownames(gm.objects))
  fm=rbind(gm.objects, fm.f)
  
  # go through each data matrix, add to seurat and fm if more matrices
  if(!is.null(list.additional)){
    
    for(i in seq(list.additional)){
      mat.add=list.additional[[i]]
      mat.add=subset(mat.add,cells = colnames(scmat))
      scmat[[names(list.additional)[i]]] <- mat.add
      scmat <- NormalizeData(scmat, assay = names(list.additional)[i], normalization.method = "CLR") # test scttransform too
      scmat <- ScaleData(scmat, assay = names(list.additional)[i])
      
      add.data.normalized <- data.matrix(GetAssayData(scmat, assay = names(list.additional)[i], slot = "data"))
      rownames(add.data.normalized)=paste0("N:",names(list.additional)[i], ":", rownames(add.data.normalized))
      fm=rbind(fm, add.data.normalized)
    }
  }
  
  # DE analysis:
  if(DE.analysis){
    markers.all=FindAllMarkers(scmat, only.pos = T)
    save(list=c("fm", "scmat", "markers.all"), file=paste0(name, "_scRNA.Rdata"))
  }else{
    save(list=c("fm", "scmat"), file=paste0(name, "_scRNA.Rdata"))
  }
  
  return(paste0(name, "_scRNA.Rdata"))
}

# make sure  data is on a linear non-log scale
# check 0 values, if a lot between 0-1 better to use add.one. If counts, remove works too quite well
# zero fix methods: https://arxiv.org/pdf/1806.06403.pdf
gm_mean=function(x, na.rm=TRUE, zero.handle=c("remove", "add.one", "zero.propagate")){
  zero.handle=match.arg(zero.handle)
  if(any(x < 0, na.rm = TRUE)){
    warning("Negative values produced NaN - is the data on linear - non-log scale?")
    return(NaN)
  }
  
  if(zero.handle=="remove"){
    return(exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x)))
  }
  
  if(zero.handle=="add.one"){
    return(exp(mean(log(x+1), na.rm = na.rm))-1)
  }
  
  if(zero.handle=="zero.propagate"){
    if(any(x == 0, na.rm = TRUE)){
      return(0)
    }
    return(exp(mean(log(x), na.rm = na.rm)))
  }
  
}

get.cluster.label.singleR=function(singler){
  
  # how to add this?
  cluster=singler$singler[[1]]$SingleR.clusters$labels
  cluster.main=singler$singler[[1]]$SingleR.clusters.main$labels
  
  # order based on cell type:
  # cat(paste0("'", colors.group[,1], "'"), sep=",")
  lineage=c('HSC','MPP','CLP','CMP','GMP','MEP','Monocytes','DC','Macrophages','Macrophages M1','Macrophages M2','CD4+ T-cells','CD8+ T-cells','Tregs','T_cell:gamma-delta','T_cell:CD4+_Naive','T_cell:CD8+_naive','T_cell:Treg:Naive','CD4+ Tcm','CD4+ Tem','CD8+ Tcm','CD8+ Tem','NK cells','naive B-cells','Memory B-cells','B-cells','Class-switched memory B-cells','Plasma cells','Endothelial cells','Neutrophils','Eosinophils','Fibroblasts','Smooth muscle','Erythroblast','Erythrocytes','Megakaryocytes')
  
  cluster=cluster[order(match(cluster[,1],lineage)),,drop=F]
  
  # order and cluster identity;
  cat(paste(rownames(cluster), collapse=","), paste(cluster, collapse=","), sep="\t")
  cat("", sep="\n")
  return(cluster)
}

# found bug in reference list, was missing, also min genes was 200 not 0 causing errors.
CreateBigSingleRObject.custom=function (counts, annot = NULL, project.name, xy, clusters, N = 10000,
                                        min.genes = 0, technology = "10X", species = "Human", citation = "",
                                        ref.list = list(), normalize.gene.length = F, variable.genes = "de",
                                        fine.tune = T, reduce.file.size = T, do.signatures = F, do.main.types = T,
                                        temp.dir = getwd(), numCores = SingleR.numCores){
  n = ncol(counts)
  s = seq(1, n, by = N)
  dir.create(paste0(temp.dir, "/singler.temp/"), showWarnings = FALSE)
  for (i in s) {
    print(i)
    A = seq(i, min(i + N - 1, n))
    
    singler = CreateSinglerObject(counts[, A], annot = annot[A], ref.list = ref.list,
                                  project.name = project.name, min.genes = min.genes,
                                  technology = technology, species = species, citation = citation,
                                  do.signatures = do.signatures, clusters = NULL, numCores = numCores)
    
    
    save(singler, file = paste0(temp.dir, "/singler.temp/",
                                project.name, ".", i, ".RData"))
  }
  singler.objects.file <- list.files(paste0(temp.dir, "/singler.temp/"),
                                     pattern = "RData", full.names = T)
  singler.objects = list()
  for (i in 1:length(singler.objects.file)) {
    load(singler.objects.file[[i]])
    singler.objects[[i]] = singler
  }
  singler = SingleR.Combine.custom(singler.objects, order = colnames(counts),expr = counts, clusters = clusters, xy = xy)
  singler
}

# bug also in this.
SingleR.Combine.custom=function (singler.list, order = NULL, clusters = NULL, expr = NULL, xy = NULL) 
{
  singler = c()
  singler$singler = singler.list[[1]]$singler
  for (j in 1:length(singler.list[[1]]$singler)) {
    singler$singler[[j]]$SingleR.cluster = c()
    singler$singler[[j]]$SingleR.cluster.main = c()
    singler$singler[[j]]$SingleR.single$clusters = c()
  }
  singler$meta.data = singler.list[[1]]$meta.data
  singler$meta.data$clusters = c()
  singler$meta.data$xy = c()
  singler$meta.data$data.sets = rep(singler$meta.data$project.name, 
                                    length(singler$meta.data$orig.ident))
  for (i in 2:length(singler.list)) {
    for (j in 1:length(singler$singler)) {
      if (singler.list[[i]]$singler[[j]]$about$RefData != 
          singler.list[[1]]$singler[[j]]$about$RefData) {
        stop("The objects are not ordered by the same reference data.")
      }
      singler$singler[[j]]$about$Organism = c(singler$singler[[j]]$about$Organism, 
                                              singler.list[[i]]$singler[[j]]$about$Organism)
      singler$singler[[j]]$about$Citation = c(singler$singler[[j]]$about$Citation, 
                                              singler.list[[i]]$singler[[j]]$about$Citation)
      singler$singler[[j]]$about$Technology = c(singler$singler[[j]]$about$Technology, 
                                                singler.list[[i]]$singler[[j]]$about$Technology)
      singler$singler[[j]]$SingleR.single$labels = rbind(singler$singler[[j]]$SingleR.single$labels, 
                                                         singler.list[[i]]$singler[[j]]$SingleR.single$labels)
      if (!is.null(singler$singler[[j]]$SingleR.single$labels1)) {
        singler$singler[[j]]$SingleR.single$labels1 = rbind(singler$singler[[j]]$SingleR.single$labels1, 
                                                            singler.list[[i]]$singler[[j]]$SingleR.single$labels1)
      }
      singler$singler[[j]]$SingleR.single$scores = rbind(singler$singler[[j]]$SingleR.single$scores, 
                                                         singler.list[[i]]$singler[[j]]$SingleR.single$scores)
      singler$singler[[j]]$SingleR.single.main$labels = rbind(singler$singler[[j]]$SingleR.single.main$labels, 
                                                              singler.list[[i]]$singler[[j]]$SingleR.single.main$labels)
      if (!is.null(singler$singler[[j]]$SingleR.single.main$labels1)) {
        singler$singler[[j]]$SingleR.single.main$labels1 = rbind(singler$singler[[j]]$SingleR.single.main$labels1, 
                                                                 singler.list[[i]]$singler[[j]]$SingleR.single.main$labels1)
      }
      singler$singler[[j]]$SingleR.single.main$scores = rbind(singler$singler[[j]]$SingleR.single.main$scores, 
                                                              singler.list[[i]]$singler[[j]]$SingleR.single.main$scores)
      singler$singler[[j]]$SingleR.single$cell.names = c(singler$singler[[j]]$SingleR.single$cell.names, 
                                                         singler.list[[i]]$singler[[j]]$SingleR.single$cell.names)
      singler$singler[[j]]$SingleR.single.main$cell.names = c(singler$singler[[j]]$SingleR.single.main$cell.names, 
                                                              singler.list[[i]]$singler[[j]]$SingleR.single.main$cell.names)
      if (!is.null(singler$singler[[j]]$SingleR.single.main$pval)) {
        singler$singler[[j]]$SingleR.single.main$pval = c(singler$singler[[j]]$SingleR.single.main$pval, 
                                                          singler.list[[i]]$singler[[j]]$SingleR.single.main$pval)
      }
      if (!is.null(singler$singler[[j]]$SingleR.single$pval)) {
        singler$singler[[j]]$SingleR.single$pval = c(singler$singler[[j]]$SingleR.single$pval, 
                                                     singler.list[[i]]$singler[[j]]$SingleR.single$pval)
      }
    }
    singler$meta.data$project.name = paste(singler$meta.data$project.name, 
                                           singler.list[[i]]$meta.data$project.name, sep = "+")
    singler$meta.data$orig.ident = c(singler$meta.data$orig.ident, 
                                     singler.list[[i]]$meta.data$orig.ident)
    singler$meta.data$data.sets = c(singler$meta.data$data.sets, 
                                    rep(singler.list[[i]]$meta.data$project.name, length(singler.list[[i]]$meta.data$orig.ident)))
  }
  for (j in 1:length(singler$singler)) {
    if (!is.null(order)) {
      singler$singler[[j]]$SingleR.single$labels = singler$singler[[j]]$SingleR.single$labels[order, 
                                                                                              ]
      if (!is.null(singler$singler[[j]]$SingleR.single$labels1)) {
        singler$singler[[j]]$SingleR.single$labels1 = singler$singler[[j]]$SingleR.single$labels1[order, 
                                                                                                  ]
      }
      singler$singler[[j]]$SingleR.single$scores = singler$singler[[j]]$SingleR.single$scores[order, 
                                                                                              ]
      singler$singler[[j]]$SingleR.single$cell.names = singler$singler[[j]]$SingleR.single$cell.names[order]
      singler$singler[[j]]$SingleR.single.main$labels = singler$singler[[j]]$SingleR.single.main$labels[order, 
                                                                                                        ]
      if (!is.null(singler$singler[[j]]$SingleR.single.main$labels1)) {
        singler$singler[[j]]$SingleR.single.main$labels1 = singler$singler[[j]]$SingleR.single.main$labels1[order, 
                                                                                                            ]
      }
      singler$singler[[j]]$SingleR.single.main$scores = singler$singler[[j]]$SingleR.single.main$scores[order, 
                                                                                                        ]
      singler$singler[[j]]$SingleR.single.main$cell.names = singler$singler[[j]]$SingleR.single.main$cell.names[order]
      if (!is.null(singler$singler[[j]]$SingleR.single$pval)) {
        singler$singler[[j]]$SingleR.single$pval = singler$singler[[j]]$SingleR.single$pval[order]
        singler$singler[[j]]$SingleR.single.main$pval = singler$singler[[j]]$SingleR.single.main$pval[order]
      }
    }
  }
  if (!is.null(clusters) && !is.null(expr)) {
    for (j in 1:length(singler$singler)) {
      if (is.character(singler$singler[[j]]$about$RefData)) {
        ref = get(tolower(singler$singler[[j]]$about$RefData))
      }
      else {
        ref = singler$singler[[j]]$about$RefData
      }
      singler$singler[[j]]$SingleR.clusters = SingleR("cluster", 
                                                      expr, ref$data, types = ref$types, clusters = factor(clusters), 
                                                      sd.thres = ref$sd.thres, genes = "de", fine.tune = T)
      singler$singler[[j]]$SingleR.clusters.main = SingleR("cluster", 
                                                           expr, ref$data, types = ref$main_types, clusters = factor(clusters), 
                                                           sd.thres = ref$sd.thres, genes = "de", fine.tune = T)
    }
    singler$meta.data$clusters = clusters
    if (!is.null(xy)) {
      singler$meta.data$xy = xy
    }
  }
  singler
}