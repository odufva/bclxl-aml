## Init
import os
import numpy as np
import numpy.random as random
import pandas as pd

from scvi.dataset.dataset import GeneExpressionDataset
from scvi.dataset.csv import CsvDataset
from scvi.inference import UnsupervisedTrainer
from scvi.models import SCANVI, VAE
from scvi.inference.autotune import auto_tune_scvi_model

from umap import UMAP

import torch
import scanpy as sc
import louvain

import logging
import pickle
from hyperopt import hp

use_cuda = True
n_epochs_all = None
save_path = ''
show_plot = True

os.chdir("/Users/janihuuh/Dropbox/bclxl/")


###### Only the two samples
## Download data
FH_9119_3   = CsvDataset(filename='results/scvi/input_files/bcl2/FH_9119_3.csv', save_path='', sep=',', new_n_genes=False)
FH_9119_3R1 = CsvDataset(filename='results/scvi/input_files/bcl2/FH_9119_R1.csv', save_path='', sep=',', new_n_genes=False)

## Combine
all_dataset = GeneExpressionDataset()
all_dataset.populate_from_per_batch_list(Xs = [FH_9119_3.X, FH_9119_3R1.X])

## Train
vae      = VAE(all_dataset.nb_genes, n_batch=all_dataset.n_batches, n_labels=all_dataset.n_labels, n_hidden=128, n_latent=30, n_layers=2, dispersion='gene')
trainer  = UnsupervisedTrainer(vae, all_dataset, train_size=1.0)
trainer.train(n_epochs=50)

## Sample posterior to get latent representation and save those embeddings
full = trainer.create_posterior(trainer.model, all_dataset, indices=np.arange(len(all_dataset)))
latent, batch_indices, labels = full.sequential().get_latent()
batch_indices = batch_indices.ravel()

## Save
torch.save(trainer.model.state_dict(), 'results/scvi/output/bcl2.pkl')
np.savetxt("results/scvi/output/bcl2_latent.csv", latent, delimiter=",")
np.savetxt("results/scvi/output/bcl2_indices.csv", batch_indices, delimiter=",")

###### only pre sample + three hemap but only blasts
## Download data
t3667       = CsvDataset(filename='results/scvi/input_files/hemap_blast_pre/3667.csv', save_path='', sep=',', new_n_genes=False)
t5897       = CsvDataset(filename='results/scvi/input_files/hemap_blast_pre/5897.csv', save_path='', sep=',', new_n_genes=False)
t9119       = CsvDataset(filename='results/scvi/input_files/hemap_blast_pre/6386.csv', save_path='', sep=',', new_n_genes=False)
FH_9119_3   = CsvDataset(filename='results/scvi/input_files/hemap_blast_pre/FH_9119_3.csv', save_path='', sep=',', new_n_genes=False)

## Combine
all_dataset = GeneExpressionDataset()
all_dataset.populate_from_per_batch_list(Xs = [t3667.X, t5897.X, t9119.X, FH_9119_3.X])

## Train
vaehemap = VAE(all_dataset.nb_genes, n_batch=all_dataset.n_batches, n_labels=all_dataset.n_labels, n_hidden=128, n_latent=30, n_layers=2, dispersion='gene')
trainerhemap = UnsupervisedTrainer(vaehemap, all_dataset, train_size=1.0)
trainerhemap.train(n_epochs=50)

## Sample posterior to get latent representation and save those embeddings
full = trainerhemap.create_posterior(trainerhemap.model, all_dataset, indices=np.arange(len(all_dataset)))
latent, batch_indices, labels = full.sequential().get_latent()
batch_indices = batch_indices.ravel()

## Save
torch.save(trainerhemap.model.state_dict(), 'results/scvi/output/hemap_blast_pre.pkl')
np.savetxt("results/scvi/output/hemap_blast_pre_latent.csv", latent, delimiter=",")
np.savetxt("results/scvi/output/hemap_blast_pre_index.csv", batch_indices, delimiter=",")
