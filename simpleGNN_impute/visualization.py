import scanpy as sc
import sys

args = sys.argv
masked_prob = str(float(args[1])).replace('.','')
true_adata = sc.read('data/processed/Klein.h5ad')
masked_adata = sc.read(f'data/masked/Klein_{masked_prob}.h5ad')
predicted_adata = sc.read(f'results/predicted_adata_{masked_prob}.h5ad')

sc.pp.neighbors(true_adata)
sc.tl.umap(true_adata)
sc.pl.umap(true_adata, color='cluster', save='Klein_true_umap.png')

sc.pp.neighbors(masked_adata)
sc.tl.umap(masked_adata)
sc.pl.umap(masked_adata, color='cluster', save=f'Klein_masked_umap_{masked_prob}.png')

sc.pp.neighbors(predicted_adata)
sc.tl.umap(predicted_adata)
sc.pl.umap(predicted_adata, color='cluster', save=f'Klein_predicted_umap_{masked_prob}.png')

