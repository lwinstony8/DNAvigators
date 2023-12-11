import os
import time
import glob
import torch

def train(gdata,
          model,
          epochs=3000,
          lr=0.001,
          patience=200):
    
    device = torch.device('cuda')
    model = model.to(device)
    optimizer = torch.optim.Adam(model.parameters(), lr=lr)
    
    lossFunc = torch.nn.MSELoss(reduction='mean')
    gdata = gdata.to(device)
    loss_values = []
    time_start = time.time()
    best_epoch = 0
    best_loss_val = float('inf')
    current_patience = 0

    def train_wrapper(epoch):
        model.train()
        optimizer.zero_grad()

        pred = model(gdata.x, gdata.edge_index, gdata.size_factors)

        dropout_pred = pred[gdata.train_mask]
        dropout_true = gdata.y[gdata.train_mask]

        loss_train = lossFunc(dropout_pred, dropout_true)

        loss_train.backward()
        optimizer.step()

        # Validation
        model.eval()
        pred = model(gdata.x, gdata.edge_index, gdata.size_factors)

        dropout_pred = pred[gdata.val_mask]
        dropout_true = gdata.y[gdata.val_mask]

        loss_val = lossFunc(dropout_pred, dropout_true)

        current_step = epoch+1

        if current_step % 10 == 0:
            print(f'Epoch: {current_step}',
                  f'loss_train: {loss_train.data.item():.4f}',
                  f'loss_val: {loss_val.data.item():.4f}')
        
        return loss_val.data.item()
    
    folder = os.path.exists('./models')
    if not folder:
        os.makedirs('./models')

    for epoch in range(epochs):
        current_loss = train_wrapper(epoch)
        loss_values.append(current_loss)

        current_step = epoch + 1

        if current_loss < best_loss_val:
            torch.save(model.state_dict(), f'./models/{epoch}.pkl')
            best_epoch = epoch
            best_loss_val = current_loss
            current_patience = 0
        else:
            current_patience += 1

        if current_patience >= patience:
            print(f'Out of patience at epoch {current_step}')
            break

    model_files = glob.glob('./models/*.pkl')
    for model_file in model_files:
        if int(model_file.split('/')[-1].split('.')[0]) != best_epoch:
            os.remove(model_file)
    
    print(f'Time elapsed: {(time.time()-time_start):.4f}s')

    model.load_state_dict(torch.load(f'./models/{best_epoch}.pkl'))
    

