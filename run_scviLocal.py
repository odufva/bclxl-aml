import scvi; import scanpy as sc; import os; import numpy as np
os.chdir("/scratch/cs/csb/projects/bclxl/")

## Load and preprocess
seurat_data = scvi.data.read_h5ad("results/scvi/input_files/hemap_blast_pre/blastpre_seurat.h5ad")
seurat_data.layers["counts"] = seurat_data.X.copy() # preserve counts
sc.pp.normalize_total(seurat_data, target_sum=1e4)
sc.pp.log1p(seurat_data)
seurat_data.raw = seurat_data # freeze the state in `.raw`

## Train
scvi.data.setup_anndata(seurat_data, layer = "counts", batch_key = "orig.ident")
seurat_data_model=scvi.model.SCVI(seurat_data)
seurat_data_model.train()

## Save it
seurat_data_model.save("results/scvi/output/blast_seurat_diet/")
latent = seurat_data_model.get_latent_representation()
np.savetxt("results/scvi/output/hemap_blast_pre_latent.csv", latent, delimiter=",")
