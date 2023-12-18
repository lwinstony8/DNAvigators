import numpy as np
import scanpy as sc
import anndata as ad

import pandas as pd

import sys
import json
import muon


cite = sc.read(f'./data/GSE194122_openproblems_neurips2021_cite_BMMC_processed.h5ad')
multiome = sc.read(f'./data/GSE194122_openproblems_neurips2021_multiome_BMMC_processed.h5ad')

rna_cite = cite[:, cite.var['feature_types'] == 'GEX'].copy()
adt = cite[:, cite.var['feature_types'] == 'ADT'].copy()

rna_multiome = multiome[:, multiome.var['feature_types'] == 'GEX'].copy()
atac = multiome[:, multiome.var['feature_types'] == 'ATAC'].copy()

# concat
rna = ad.concat([rna_cite, rna_multiome])
# normalize
rna.X = rna.layers['counts'].copy()
sc.pp.normalize_total(rna, target_sum=1e4)
sc.pp.log1p(rna)
# subset to hvg
sc.pp.highly_variable_genes(rna, n_top_genes=4000, batch_key='Samplename')
rna = rna[:, rna.var.highly_variable].copy()
# split again
rna_cite = rna[rna.obs['Modality'] == 'cite'].copy()
rna_multiome = rna[rna.obs['Modality'] == 'multiome'].copy()

adt.X = adt.layers['counts'].copy()
muon.prot.pp.clr(adt)
adt.layers['clr'] = adt.X.copy()

atac.X = atac.layers['counts'].copy()
sc.pp.normalize_total(atac, target_sum=1e4)
sc.pp.log1p(atac)
atac.layers['log-norm'] = atac.X.copy()
sc.pp.highly_variable_genes(atac, n_top_genes=20000, batch_key='batch')
atac = atac[:, atac.var.highly_variable].copy()

with open('./data/celltype_harmonize.json', 'r') as f:
    harmonized_celltypes = json.load(f)

rna_multiome.obs['l1_cell_type'] = rna_multiome.obs['cell_type'].map(harmonized_celltypes['multi_ct_l1_map'])
rna_multiome.obs['l2_cell_type'] = rna_multiome.obs['cell_type'].map(harmonized_celltypes['multi_ct_l2_map'])

atac.obs['l1_cell_type'] = atac.obs['cell_type'].map(harmonized_celltypes['multi_ct_l1_map'])
atac.obs['l2_cell_type'] = atac.obs['cell_type'].map(harmonized_celltypes['multi_ct_l2_map'])

rna_cite.obs['l1_cell_type'] = rna_cite.obs['cell_type'].map(harmonized_celltypes['cite_ct_l1_map'])
rna_cite.obs['l2_cell_type'] = rna_cite.obs['cell_type'].map(harmonized_celltypes['cite_ct_l2_map'])

adt.obs['l1_cell_type'] = adt.obs['cell_type'].map(harmonized_celltypes['cite_ct_l1_map'])
adt.obs['l2_cell_type'] = adt.obs['cell_type'].map(harmonized_celltypes['cite_ct_l2_map'])

rna_multiome.write(f'./data/multigrate_processed/rna_multiome.h5ad')
atac.write(f'./data/multigrate_processed/atac_multiome.h5ad')
rna_cite.write(f'./data/multigrate_processed/rna_cite.h5ad')
adt.write(f'./data/multigrate_processed/adt_cite.h5ad')