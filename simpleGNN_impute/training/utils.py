import torch
import numpy as np
import networkx as nx
import scanpy as sc

from torch_geometric.data import Data
from sklearn.decomposition import PCA
from sklearn.neighbors import kneighbors_graph

def normalize(adata):
    # sc.pp.filter_genes(adata, min_counts=1)
    # sc.pp.filter_cells(adata, min_counts=1)

    sc.pp.normalize_per_cell(adata)
    adata.obs['size_factors'] = adata.obs.n_counts / np.median(adata.obs.n_counts)

    sc.pp.log1p(adata)
    sc.pp.scale(adata)

    return adata

def train_val_split(adata, train_size=0.6, val_size=0.2, test_size=0.2):
    assert train_size + val_size + test_size == 1
    adata = adata.copy()
    cell_nums = adata.n_obs
    test_val = np.random.choice(cell_nums, 
                                int(cell_nums * (val_size + test_size)), 
                                replace=False)
    train_index = [i for i in list(range(cell_nums)) if i not in test_val]
    test_index = np.random.choice(test_val, int(len(test_val) * (test_size / (val_size + test_size))), replace=False)
    val_index = [i for i in test_val if i not in test_index]
    
    tmp = np.zeros(cell_nums, dtype=bool)
    tmp[train_index] = True
    adata.obs['train_index'] = tmp

    tmp = np.zeros(cell_nums, dtype=bool)
    tmp[val_index] = True
    adata.obs['val_index'] = tmp

    tmp = np.zeros(cell_nums, dtype=bool)
    tmp[test_index] = True
    adata.obs['test_index'] = tmp

    return adata


def kneighbor(adata, n_components=50, k=5):
    pca = PCA(n_components=n_components)
    data_pca = pca.fit_transform(adata.X)

    A = kneighbors_graph(data_pca, k)
    G = nx.from_numpy_array(A.todense())

    edges = []
    for (u, v) in G.edges():
        edges.append([u, v])
        edges.append([v, u])

    edges = np.array(edges).T

    return edges


def adata2gdata(adata):
    edges = kneighbor(adata)

    edges = torch.tensor(edges, dtype=torch.long)
    features = torch.tensor(adata.X, dtype=torch.float)
    size_factors = torch.tensor(adata.obs.size_factors, dtype=torch.float).reshape(-1, 1)
    
    labels = torch.tensor(adata.raw.X.A, dtype=torch.float)

    gdata = Data(x=features, y=labels, edge_index=edges, size_factors=size_factors)
    gdata.train_mask = torch.tensor(adata.obs.train_index, dtype=torch.bool)
    gdata.val_mask = torch.tensor(adata.obs.val_index, dtype=torch.bool)

    return gdata