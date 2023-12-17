import scanpy as sc

rna_multiome = sc.read(f'./data/multigrate_processed/rna_multiome.h5ad')
sample_names = rna_multiome.obs['batch'].unique().tolist()
print(sample_names)
# cell_types = rna_multiome.obs['l1_cell_type'].unique().tolist()
# print(cell_types)