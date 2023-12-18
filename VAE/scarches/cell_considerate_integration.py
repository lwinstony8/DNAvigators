# This VAE encodes cell type but does not train separately

import scarches as sca
import scanpy as sc
import anndata as ad
import numpy as np
import muon
import gdown
import json

# import multigrate 

import warnings
warnings.filterwarnings("ignore")

sc.set_figure_params(figsize=(4, 4), fontsize=8)

rna_multiome = sc.read(f'./data/multigrate_processed/rna_multiome.h5ad')
atac = sc.read(f'./data/multigrate_processed/atac_multiome.h5ad')
rna_cite = sc.read(f'./data/multigrate_processed/rna_cite.h5ad')
adt = sc.read(f'./data/multigrate_processed/adt_cite.h5ad')


# define the reference and the query batches
cite_reference_batches = ['s1d1', 's1d2', 's1d3']
multiome_reference_batches = ['s1d1', 's1d2', 's1d3']
cite_query_batches = ['s2d1', 's2d4']
multiome_query_batches = ['s2d1', 's2d4']
# query
rna_multiome_query = rna_multiome[rna_multiome.obs['batch'].isin(multiome_query_batches)].copy()
atac_query = atac[atac.obs['batch'].isin(multiome_query_batches)].copy()
rna_cite_query = rna_cite[rna_cite.obs['batch'].isin(cite_query_batches)].copy()
adt_query = adt[adt.obs['batch'].isin(cite_query_batches)].copy()
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
    categorical_covariate_keys=['Modality', 'Samplename', 'l2_cell_type'],
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

model.train(max_epochs=10)

model.get_latent_representation()
adata.obsm['latent_ref'] = adata.obsm['latent'].copy()
# adata.write(f'./data/multigrate_processed/reference_only.h5ad')

sc.pp.neighbors(adata, use_rep='latent')
sc.tl.umap(adata)
sc.pl.umap(adata, color=['l1_cell_type', 'l2_cell_type', 'Modality', 'Samplename'], frameon=False, ncols=1, save='multivae_reference.png')


query = sca.models.organize_multiome_anndatas(
    adatas = [[rna_cite_query, rna_multiome_query], [None, atac_query], [adt_query, None]],
    layers = [['counts', 'counts'], [None, 'log-norm'], ['clr', None]],
)

sca.models.MultiVAE.setup_anndata(
    query,
    categorical_covariate_keys=['Modality', 'Samplename', 'l2_cell_type'],
    rna_indices_end=4000,
)

idx_atac_query = query.obs['Samplename'] == 'site2_donor4_multiome'
idx_scrna_query = query.obs['Samplename'] == 'site2_donor1_cite'

idx_mutiome_query = query.obs['Samplename'] == 'site2_donor1_multiome'
idx_cite_query = query.obs['Samplename'] == 'site2_donor4_cite'

query[idx_atac_query, :4000].X = 0
query[idx_scrna_query, 4000:].X = 0

q_model = sca.models.MultiVAE.load_query_data(query, model)
q_model.train(weight_decay=0, max_epochs=10)

q_model.get_latent_representation(adata=query)
q_model.get_latent_representation(adata=adata)

adata.obs['reference'] = 'reference'
query.obs['reference'] = 'query'

adata.obs['type_of_query'] = 'reference'
query.obs.loc[idx_atac_query, 'type_of_query'] = 'ATAC query'
query.obs.loc[idx_scrna_query, 'type_of_query'] = 'scRNA query'
query.obs.loc[idx_mutiome_query, 'type_of_query'] = 'multiome query'
query.obs.loc[idx_cite_query, 'type_of_query'] = 'CITE-seq query'

adata_both = ad.concat([adata, query])
sc.pp.neighbors(adata_both, use_rep='latent')
sc.tl.umap(adata_both)
sc.pl.umap(adata_both, color=['l1_cell_type', 'l2_cell_type', 'reference', 'Modality', 'Samplename'], ncols=1, frameon=False, save='multivae_both.png')

sc.pl.umap(
    adata_both,
    color='type_of_query',
    ncols=1,
    frameon=False,
    groups=['CITE-seq query'],
    save='multivae_cite_query.png'
)

sc.pl.umap(
    adata_both,
    color='type_of_query',
    ncols=1,
    frameon=False,
    groups=['multiome query'],
    save='multivae_multiome_query.png'
)

sc.pl.umap(
    adata_both,
    color='type_of_query',
    ncols=1,
    frameon=False,
    groups=['scRNA query'],
    save='multivae_scRNA_query.png'
)

sc.pl.umap(
    adata_both,
    color='type_of_query',
    ncols=1,
    frameon=False,
    groups=['ATAC query'],
    save='multivae_atac_query.png'
)
