import scanpy as sc

rna_multiome = sc.read(f'./data/multigrate_processed/rna_multiome.h5ad')
sample_names = rna_multiome.obs['batch'].unique().tolist()
print(sample_names)

def explore_obs(data, tag):
    values, counts = np.unique(data.obs[tag], return_counts=True)
    print(f"{len(values)} unique {tag} elements")
    print(f"--------------------------------------------")
    
    max_value_len = max([len(str(value)) for value in values])
    max_count_len = max([len(str(count)) for count in counts])
    
    for i, (value, count) in enumerate(reversed(sorted(list(zip(values, counts)), key = lambda x: x[1]))):
        print(f"{value:<{max_value_len+3}} {count:>{max_count_len}}")

# cell_types = rna_multiome.obs['l1_cell_type'].unique().tolist()
# print(cell_types)