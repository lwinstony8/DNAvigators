# This VAE trains a separate model for each cell type, creating separate latent spaces
# These outputs may potentially be combined for a representation using concat, which could serve as a valuable input back to multigrate

import scarches as sca
import scanpy as sc
import anndata as ad
import numpy as np
import muon
import gdown
import json

import pickle

import warnings
warnings.filterwarnings("ignore")

sc.set_figure_params(figsize=(4, 4), fontsize=8)

cell_types = ['B', 'Myeloid', 'T', 'Other Lympoid', 'Erythroid', 'NK', 'HSC', 'Plasma']
sample_names = ['s1d1', 's1d2', 's1d3', 's2d1', 's2d4', 's2d5', 's3d3', 's3d6', 's3d7', 's3d10', 's4d1', 's4d8', 's4d9']

for cell_type in cell_types:
    rna_multiome = sc.read(f'./data/multigrate_processed/l1_cell_types/{cell_type}/rna_multiome.h5ad')
    atac = sc.read(f'./data/multigrate_processed/l1_cell_types/{cell_type}/atac.h5ad')
    rna_cite = sc.read(f'./data/multigrate_processed/l1_cell_types/{cell_type}/rna_cite.h5ad')
    adt = sc.read(f'./data/multigrate_processed/l1_cell_types/{cell_type}/adt.h5ad')


    # define the reference and the query batches
    cite_reference_batches = sample_names
    multiome_reference_batches = sample_names
    # reference
    rna_multiome = rna_multiome[rna_multiome.obs['batch'].isin(multiome_reference_batches)].copy()
    atac = atac[atac.obs['batch'].isin(multiome_reference_batches)].copy()
    rna_cite = rna_cite[rna_cite.obs['batch'].isin(cite_reference_batches)].copy()
    adt = adt[adt.obs['batch'].isin(cite_reference_batches)].copy()

    adata = sca.models.organize_multiome_anndatas(
        adatas = [[rna_cite, rna_multiome], [None, atac], [adt, None]],    # a list of anndata objects per modality, RNA-seq always goes first
        layers = [['counts', 'counts'], [None, 'log-norm'], ['clr', None]], # if need to use data from .layers, if None use .X
    )

    sca.models.MultiVAE.setup_anndata(
        adata,
        categorical_covariate_keys=['Modality', 'Samplename'],
        rna_indices_end=4000,
    )

    model = sca.models.MultiVAE(
        adata,
        losses=['nb', 'mse', 'mse'],
        loss_coefs={'kl': 1e-1,
                'integ': 3000,
                },
        integrate_on='Modality',
        mmd='marginal',
    )

    model.train(max_epochs=25,
                early_stopping = False)

    model.get_latent_representation()
    adata.obsm['latent_ref'] = adata.obsm['latent'].copy()

    path = f'./data/multigrate_processed/l1_cell_types/{cell_type}/reference_latent_all_samples.pkl'

    with open(path, 'wb') as file:
        pickle.dump(adata, file)

    # adata.write(path)

    sc.pp.neighbors(adata, use_rep='latent')
    sc.tl.umap(adata)
    sc.pl.umap(adata, color=['l1_cell_type', 'l2_cell_type', 'Modality', 'Samplename'], frameon=False, ncols=1, save=f'multivae_reference_{cell_type}_all_samples.png')
