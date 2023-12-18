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

#adata_dict = {}
adatas = []
for cell_type in cell_types:
    with open(f'./data/multigrate_processed/l1_cell_types/{cell_type}/reference_latent_all_samples.pkl', 'rb') as file:
        adata = pickle.load(file)
    adatas.append(adata)

adata_concat = ad.concat(adatas,merge=None, axis=0)

#(149989, 24134)
print(adata_concat.shape)
#(149989, 15)
print(adata_concat.obsm['latent'].shape)

new_adata = 

# Potential idea-- merge latent spaces of each celltype and run the model on the concatenated spaces.

# adata_concat = [adata_dict[cell_type] for cell_type in adata_dict.keys()]

# adata_concat = sca.models.organize_multiome_anndatas(
#     adatas = [[adata_concat]],    # a list of anndata objects per modality, RNA-seq always goes first
#     layers = [['log-norm']], # if need to use data from .layers, if None use .X
# )

# adata_concat = sca.models.organize_multiome_anndatas(
#     adatas = [[adata_concat]],    # a list of anndata objects per modality, RNA-seq always goes first
#     layers = [['log-norm']*len(adata)], # if need to use data from .layers, if None use .X
# )

# print(adata_concat.shape)

# sca.models.MultiVAE.setup_anndata(
#     adata_concat,
#     categorical_covariate_keys=['Modality', 'Samplename', 'l1_cell_type'],
#     rna_indices_end=4000,
# )

# print(adata_concat.shape)
# print(adata.uns)

# model = sca.models.MultiVAE(
#     adata_concat,
#     losses=['nb', 'mse', 'mse'],
#     loss_coefs={'kl': 1e-1,
#                'integ': 3000,
#                },
#     integrate_on='l1_cell_type',
#     mmd='marginal',
# )

# model.train(max_epochs=10)

# model.get_latent_representation()
# adata_concat.obsm['latent_ref'] = adata_concat.obsm['latent'].copy()

# sc.pp.neighbors(adata_concat, use_rep='latent')
# sc.tl.umap(adata_concat)
# sc.pl.umap(adata_concat, color=['l1_cell_type', 'l2_cell_type', 'Modality', 'Samplename'], frameon=False, ncols=1, save='multivae_reference.png')