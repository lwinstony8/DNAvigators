import torch
import torch.nn.functional as F
from torch_geometric.nn import GCNConv, GATConv

from .train import train
from .utils import adata2gdata, train_val_split, normalize

def layer(layer_type, **kwargs):
    if layer_type == 'GATConv':
        return GATConv(in_channels=kwargs['in_channels'],
                       out_channels=kwargs['out_channels'],
                       heads=kwargs['heads'], 
                       concat=kwargs['concat'],
                       dropout=kwargs['dropout'])
    if layer_type == 'GCNConv':
        return GCNConv(in_channels=kwargs['in_channels'],
                       out_channels=kwargs['out_channels'])

    
class simpleGNN(torch.nn.Module):
    def __init__(self, input_dim, h_dim=512, layer_type='GATConv', heads=3):
        super(simpleGNN, self).__init__()

        self.gnn_conv1 = layer(layer_type=layer_type, 
                               in_channels=input_dim,
                               out_channels=h_dim,
                               heads=heads,
                               concat=False,
                               dropout=0.6)
        self.gnn_norm1 = torch.nn.BatchNorm1d(h_dim)

        self.gnn_conv2 = layer(layer_type=layer_type, 
                               in_channels=h_dim,
                               out_channels=h_dim,
                               heads=heads,
                               concat=False,
                               dropout=0.6)
        self.gnn_norm2 = torch.nn.BatchNorm1d(h_dim)

        self.dense_linear1 = torch.nn.Linear(h_dim, h_dim)
        self.dense_norm1 = torch.nn.BatchNorm1d(h_dim)
        self.dense_linear2 = torch.nn.Linear(h_dim, input_dim)

    def gnn(self, x, edge_index):
        x = F.relu(self.gnn_norm1(self.gnn_conv1(x, edge_index)))
        x = F.dropout(x, p=0.5, training=self.training)
        x = F.relu(self.gnn_norm2(self.gnn_conv2(x, edge_index)))
        x = F.dropout(x, p=0.5, training=self.training)
        return x
    
    def dense(self, x):
        x = F.relu(self.dense_norm1(self.dense_linear1(x)))
        x = F.relu(self.dense_linear2(x))
        return x
    
    def forward(self, x, edge_index, size_factors):
        z = self.gnn(x, edge_index)
        x = self.dense(z)
        x = x * size_factors
        return x
    
def run_model(adata,
              layer_type='GATConv',
              epochs=3000,
              lr=0.001,
              patience=200,
              heads=3):
    
    input_dim = adata.n_vars

    model = simpleGNN(input_dim=input_dim,
                      layer_type=layer_type,
                      heads=heads,
                      h_dim=512)
    adata = normalize(adata)
    adata = train_val_split(adata)
    gdata = adata2gdata(adata)

    train(gdata=gdata,
          model=model,
          epochs=epochs,
          lr=lr,
          patience=patience)
    
    pred = model(gdata.x, gdata.edge_index, gdata.size_factors)
    adata.X = pred.detach().cpu()
    return adata