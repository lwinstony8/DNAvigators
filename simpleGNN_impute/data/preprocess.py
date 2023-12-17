import os
import scanpy as sc
from scipy import sparse

adata = sc.read('./data/haniffa21.processed.h5ad')
print('finished readingmini')
adata.X = sparse.csr_matrix(adata.X)

sc.pp.calculate_qc_metrics(adata, percent_top=None, log1p=False, inplace=True)

adata = adata[adata.obs.total_counts < 75000, :]

sc.pl.scatter(adata, x='total_counts', y='n_genes_by_counts', save='Klein_true_scatter')
sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts'], jitter=False, multi_panel=True, save='Klein_true_violin')

sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)

adata.raw = adata


folder = os.path.exists('./data/processed')

if not folder:
    os.makedirs('./data/processed')

adata.write('./data/processed/haniffa.h5ad')

