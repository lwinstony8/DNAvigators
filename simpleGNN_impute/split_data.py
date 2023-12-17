import scarches as sca
import scanpy as sc
import anndata as ad
import numpy as np
import muon
import gdown
import json
import os

# import multigrate 

import warnings
warnings.filterwarnings("ignore")

sc.set_figure_params(figsize=(4, 4), fontsize=8)

rna_multiome = sc.read(f'./data/multigrate_processed/rna_multiome.h5ad')
atac = sc.read(f'./data/multigrate_processed/atac_multiome.h5ad')
rna_cite = sc.read(f'./data/multigrate_processed/rna_cite.h5ad')
adt = sc.read(f'./data/multigrate_processed/adt_cite.h5ad')

# print(rna_multiome.obs['l2_cell_type'])

cell_types = rna_multiome.obs['l2_cell_type'].unique()

for cell_type in cell_types:
    folder = f'data/multigrate_processed/l2_cell_types/{cell_type}'
    os.makedirs(folder)
    subset_rna_multiome = rna_multiome[rna_multiome.obs['l2_cell_type'] == cell_type].copy()
    subset_rna_cite = rna_cite[rna_cite.obs['l2_cell_type'] == cell_type].copy()
    subset_atac = atac[atac.obs['l2_cell_type'] == cell_type].copy()
    subset_adt = adt[adt.obs['l2_cell_type'] == cell_type].copy()

    subset_rna_multiome.write(f'./data/multigrate_processed/l2_cell_types/{cell_type}/rna_multiome.h5ad')
    subset_rna_cite.write(f'./data/multigrate_processed/l2_cell_types/{cell_type}/rna_cite.h5ad')
    subset_atac.write(f'./data/multigrate_processed/l2_cell_types/{cell_type}/atac.h5ad')
    subset_adt.write(f'./data/multigrate_processed/l2_cell_types/{cell_type}/adt.h5ad')

