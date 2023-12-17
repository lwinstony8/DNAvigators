import numpy as np
import scanpy as sc
from scipy import sparse
import json

from sklearn.cluster import KMeans
from sklearn.metrics import mean_squared_error, mean_absolute_error
from sklearn.metrics.cluster import adjusted_rand_score, normalized_mutual_info_score
from scipy.stats import pearsonr
from sklearn.metrics.pairwise import cosine_similarity

# from training.api import simpleGNN
from training.model import run_model
import os
import pandas as pd

import sys

# adata = sc.read_h5ad(f'./data/haniffa21.processed.h5ad')
multiome = sc.read(f'./data/GSE194122_openproblems_neurips2021_multiome_BMMC_processed.h5ad')

multiome = multiome[:, multiome.var['feature_types']=='GEX']
multiome.X = multiome.layers['counts']
# atac = multiome[:, multiome.var['feature_types']=='ATAC'].copy()

# with open(f'./data/celltype_harmonize.json', 'r') as f:
#     harmonized_celltypes = json.load(f)

# rna_multiome.obs['l1_cell_type'] = rna_multiome.obs['cell_type'].map(harmonized_celltypes['multi_ct_l1_map'])
# rna_multiome.obs['l2_cell_type'] = rna_multiome.obs['cell_type'].map(harmonized_celltypes['multi_ct_l2_map'])

# atac.obs['l1_cell_type'] = atac.obs['cell_type'].map(harmonized_celltypes['multi_ct_l1_map'])
# atac.obs['l2_cell_type'] = atac.obs['cell_type'].map(harmonized_celltypes['multi_ct_l2_map'])

print('after reading')

multiome = run_model(adata=multiome,
                  layer_type='GATConv',
                  epochs=3000,
                  lr=0.001,
                  patience=200,
                  heads=3)

folder = os.path.exists('./results')

if not folder:
    os.makedirs('./results')

multiome.write(f'./results/predicted_rna_multiome.h5ad')