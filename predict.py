import argparse
import h5py

import torch
import torch.nn.functional as F
from model import DeepDrug3D

def load_data(path):
    with h5py.File(path, 'r') as f:
        X = f['X'][()]
    X = torch.FloatTensor(X)
    return X

def predict(path, model_path):
    device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")
    print('Current device: ' + str(device))
    
    net = DeepDrug3D(14)
    net.load_state_dict(torch.load(model_path))
    net.to(device)
    data = load_data(path)
    data = data.to(device)
    output = net(data)
    proba = F.softmax(output, dim=1)
    score = proba.data.cpu().numpy()
    print('The probability of pocket provided binds with ATP ligands: {:.4f}'.format(score[0]))
    print('The probability of pocket provided binds with Heme ligands: {:.4f}'.format(score[1]))
    print('The probability of pocket provided binds with other ligands: {:.4f}'.format(score[2]))

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--f', type=str, required=True, help='Input h5 file name')
    parser.add_argument('--m', type=str, required=True, help='Path to the trained model weights')
    opt = parser.parse_args()

    predict(opt.f, opt.m)