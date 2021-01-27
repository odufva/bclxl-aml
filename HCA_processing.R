# # HumanCellAtlas data
# file="/research/groups/biowhat_share/public_data/scRNAseq/Human_Cell_Atlas_Preview_Datasets/scanpy/hg19/new_2019/cellranger3/HCA.h5ad"
# library(Seurat)
# hca <- ReadH5AD(file)
# Idents(hca) <- "louvain"
# save(hca, file="/research/groups/biowhat_share/public_data/scRNAseq/Human_Cell_Atlas_Preview_Datasets/scanpy/hg19/new_2019/cellranger3/HumanCellAtlas_scRNA.Rdata")

# file="/research/groups/biowhat_share/public_data/scRNAseq/Human_Cell_Atlas_Preview_Datasets/scanpy/hg19/new_2019/cellranger3/HCA_raw.h5ad"
# counts=t(data.table::fread("/research/groups/biowhat_share/public_data/scRNAseq/Human_Cell_Atlas_Preview_Datasets/scanpy/hg19/new_2019/cellranger3/data/raw/X_allgenes.csv.gz", data.table = F))
# rownames(counts)=scan("/research/groups/biowhat_share/public_data/scRNAseq/Human_Cell_Atlas_Preview_Datasets/scanpy/hg19/new_2019/cellranger3/data/raw/var_allgenes.txt","genes")
# colnames(counts)=scan("/research/groups/biowhat_share/public_data/scRNAseq/Human_Cell_Atlas_Preview_Datasets/scanpy/hg19/new_2019/cellranger3/data/raw/obs_allgenes.txt","sample")
# save(counts, file="/research/groups/biowhat_share/public_data/scRNAseq/Human_Cell_Atlas_Preview_Datasets/scanpy/hg19/new_2019/cellranger3/HumanCellAtlas_scRNA_rawCounts.Rdata")

#******************************************* all cell types *********************************************

library(Matrix)
library(Seurat)
library(dplyr)
source("functions.scRNA.analysis.R")

load("HumanCellAtlas_scRNA_rawCounts.Rdata")
load("HumanCellAtlas_scRNA.Rdata")

name="HumanCellAtlas"

library("HGNChelper")
update=checkGeneSymbols(rownames(counts))
genelist1=update[,1]
genelist1[!is.na(update[,3])]=update[!is.na(update[,3]),3]
rownames(counts)=genelist1

# batch
batch=gsub("_HiSeq_.*.", "", colnames(counts))

# similar name between filtered HCA.h5ad matrix and non-filtered count matrix
name=paste0(colnames(counts), "-", gsub("_HiSeq_.*.", "", colnames(counts)))

test1=sc.data.analysis(scmat = counts[,name%in%colnames(hca)], regress.cell.label = batch, batch.correction.method = "MNNcorrect", name="HCA", nr.pcs = 50, check.pcs=F, plot.umap = T, cores = 10, percent.mitoDNA = 10, nFeature.min = 0, nFeature.max = 6000)
