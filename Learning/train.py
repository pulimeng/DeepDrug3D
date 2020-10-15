import argparse

import os
import random
import pandas as pd
import numpy as np
import time

import torch

from load_data import VoxelDataset
from torch.utils.data import DataLoader, Subset

from model import DeepDrug3D

from sklearn.metrics import confusion_matrix
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import accuracy_score

seed = 12306
random.seed(seed)
torch.manual_seed(seed)
if torch.cuda.is_available():
    torch.cuda.manual_seed_all(seed)
    
device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
print('Current device: ' + str(device))

def main(opt):
    in_channel = 14
    model = DeepDrug3D(in_channel)
    print(model)
    model = model.to(device)
    criterion = torch.nn.CrossEntropyLoss()
    
    if opt.opath is None:
        os.mkdir('./logs')
        opt.opath = './logs'
    
    labels = pd.read_csv(opt.lpath, names=['id','class'])
    xid = labels['id'].tolist()
    ys = labels['class'].tolist()
    dataset = VoxelDataset(label_file=opt.lpath, root_dir=opt.path)
    kfold = StratifiedKFold(n_splits=5, shuffle=True, random_state=seed)
    bs = opt.bs
    f_cnt = 0
    for train_id, val_id in kfold.split(xid, ys):
        train_set = Subset(dataset, train_id)
        train_loader = DataLoader(train_set, batch_size=bs, shuffle=True)
        val_set = Subset(dataset, val_id)
        val_loader = DataLoader(val_set, batch_size=bs, shuffle=True)
        
        tr_losses = np.zeros((opt.epoch,))
        tr_accs = np.zeros((opt.epoch,))
        val_losses = np.zeros((opt.epoch,))
        val_accs = np.zeros((opt.epoch,))
        
        model.reset_parameters()
        optimizer = torch.optim.Adam(model.parameters(), lr=opt.lr)
        
        best_val_loss = 1e6
        
        print('===================Fold {} starts==================='.format(f_cnt+1))
        for epoch in range(50):
            s = time.time()
            
            model.train()
            losses = 0
            acc = 0
            
            for i, sampled_batch in enumerate(train_loader):
                data = sampled_batch['voxel']
                y = sampled_batch['label'].squeeze()
                data = data.type(torch.FloatTensor)
                if in_channel == 1:
                   data = torch.unsqueeze(data,1)
                y = y.to(device)
                data = data.to(device)
                optimizer.zero_grad()
                output = model(data)
                loss = criterion(output, y)
                loss.backward()
                optimizer.step()
                
                y_true = y.cpu().numpy()
                y_pred = output.data.cpu().numpy().argmax(axis=1)
                acc += accuracy_score(y_true, y_pred)*100
                losses += loss.data.cpu().numpy()

            tr_losses[epoch] = losses/(i+1)
            tr_accs[epoch] = acc/(i+1)
            
            model.eval()
            v_losses = 0
            v_acc = 0
            y_preds = []
            y_trues = []
            
            for j, sampled_batch in enumerate(val_loader):
                data = sampled_batch['voxel']
                y = sampled_batch['label'].squeeze()
                data = data.type(torch.FloatTensor)
                if in_channel == 1:
                   data = torch.unsqueeze(data,1)
                y = y.to(device)
                data = data.to(device)
                with torch.no_grad():
                    output = model(data)
                    loss = criterion(output, y)
                
                y_pred = output.data.cpu().numpy().argmax(axis=1)
                y_true = y.cpu().numpy()
                y_trues += y_true.tolist()
                y_preds += y_pred.tolist()
                v_acc += accuracy_score(y_true, y_pred)*100
                v_losses += loss.data.cpu().numpy()
                
            cnf = confusion_matrix(y_trues, y_preds)        
            val_losses[epoch] = v_losses/(j+1)
            val_accs[epoch] = v_acc/(j+1)
            
            current_val_loss = v_losses/(j+1)
            if current_val_loss < best_val_loss:
                best_val_loss = current_val_loss
                torch.save(model.state_dict(), os.path.join(opt.opath, 'best_model_fold_{}.ckpt'.format(f_cnt+1)))
            
            print('Epoch: {:03d} | time: {:.4f} seconds\n'
                  'Train Loss: {:.4f} | Train accuracy {:.4f}\n'
                  'Validation Loss: {:.4f} | Validation accuracy {:.4f} | Best {:.4f}'.format(epoch+1, time.time()-s, losses/(i+1),
                                    acc/(i+1), v_losses/(j+1), v_acc/(j+1), best_val_loss))
            print('Validation confusion matrix:')
            print(cnf)

        print('===================Fold {} ends==================='.format(f_cnt+1))
        np.save(os.path.join(opt.opath, 'train_loss_{}.npy'.format(f_cnt+1)), tr_losses)
        np.save(os.path.join(opt.opath, 'train_acc_{}.npy'.format(f_cnt+1)), tr_accs)
        np.save(os.path.join(opt.opath, 'val_loss_{}.npy'.format(f_cnt+1)), val_losses)
        np.save(os.path.join(opt.opath, 'val_acc_{}.npy'.format(f_cnt+1)), val_accs)
    
        f_cnt += 1
        
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--path', type=str, required=True, help='path to data folder')
    parser.add_argument('--lpath', type=str, required=True, help='path to label file')
    parser.add_argument('--opath', type=str, required=False, help='output folder name')
    parser.add_argument('--bs', type=int, required=True, help='batch size')
    parser.add_argument('--lr', type=float, required=True, help='learning rate')
    parser.add_argument('--epoch', type=int, action='store_true', help='number of epochs to train for')
    
    opt = parser.parse_args()
    main(opt)
