import pickle

cell_types = ['B', 'Myeloid', 'T', 'Other Lympoid', 'Erythroid', 'NK', 'HSC', 'Plasma']

adata_dict = {}
for cell_type in cell_types:
    with open(f'./data/multigrate_processed/l1_cell_types/{cell_type}/reference_latent_all_samples.pkl', 'rb') as file:
        adata = pickle.load(file)
    adata_dict[cell_type] = adata

adata_concat = sc.anndata.concat([adata_dict[cell_type] for cell_type in adata_dict.values()],merge=None)

