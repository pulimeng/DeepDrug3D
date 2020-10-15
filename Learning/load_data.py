import os
import pandas as pd
import h5py

import torch
from torch.utils.data import Dataset

def load_feature(path):
    with h5py.File(path, 'r') as f:
        X = f['X'][()]
    X = torch.FloatTensor(X)
    return X

class VoxelDataset(Dataset):

    def __init__(self, label_file, root_dir):
        """
        Args:
            label_file (string): Path to the csv file with labels.
            root_dir (string): Directory with all the voxel data.
            transform (callable, optional): Optional transform to be applied
                on a sample.
        """
        self.labels = pd.read_csv(label_file)
        self.root_dir = root_dir

    def __len__(self):
        return len(self.labels)

    def __getitem__(self, idx):
        if torch.is_tensor(idx):
            idx = idx.tolist()

        voxel_name = os.path.join(self.root_dir,
                                self.labels.iloc[idx, 0]) + '.h5'
        voxel = load_feature(voxel_name)
        label = self.labels.iloc[idx, 1]
        sample = {'voxel': voxel, 'label': label}

        return sample
