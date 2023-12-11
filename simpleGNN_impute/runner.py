import numpy as np
import scanpy as sc
from scipy import sparse

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

adata = sc.read_h5ad(f'data/haniffa21.processed.h5ad')

print('after reading')

adata = run_model(adata=adata,
                  layer_type='GATConv',
                  epochs=3000,
                  lr=0.001,
                  patience=200,
                  heads=3)

folder = os.path.exists('./results')

if not folder:
    os.makedirs('./results')

adata.write(f'./results/predicted_adata.h5ad')