from .model import simpleGNN as Model
from .train import train
from .utils import adata2gdata, train_val_split, normalize


def simpleGNN(adata,
              layer_type='GATConv',
              epochs=3000,
              lr=0.001,
              patience=200,
              heads=3,):
    input_dim = adata.n_vars
    model = Model(input_dim=input_dim,
                  h_dim=512,
                  layer_type=layer_type,
                  heads=heads)
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