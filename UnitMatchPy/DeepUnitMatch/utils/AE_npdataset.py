import os
import random
import numpy as np
import h5py
from torch.utils.data import Dataset


if __name__ == '__main__':
    from train_utils import read_good_ids
else:
    from utils.train_utils import read_good_ids

    
class AE_NeuropixelsDataset(Dataset):
    """
    Dataset class for Autoencoder training on Neuropixels data.
    Each item is a waveform of shape (60, 30) corresponding to a single waveform
    from a good unit. The dataset randomly selects between the first and second
    half of the waveform data if available.

    Args:
        root (str): Root directory containing mouse subdirectories with .npy files.
        batch_size (int): Number of units to load per batch.
        mice (list, optional): List of mouse names to include. If None, all mice in root are used.
    """

    def __init__(self, root: str, batch_size = 32, mice = None):
        self.root = root
        if mice is not None:
            self.mouse_names = mice
        else:
            self.mouse_names = os.listdir(self.root)
        self.batch_size = batch_size

        self.np_file_names = read_good_ids(self.root, self.batch_size, self.mouse_names, finetune=False)
        self.n_neurons = len(self.np_file_names)
        if self.n_neurons < 1:
            print("No data! Try reducing batch size?")

    def __len__(self):
        return self.n_neurons

    def __getitem__(self, i):
        file_name = self.np_file_names[i]
        # Randomly pick 0 or 1 to choose the first half or second half of the data
        half = random.randint(0, 1)
        with h5py.File(file_name, 'r') as f:
            waveform = f['waveform'][()] 
            # MaxSitepos = f['MaxSitepos'][()]
        if half == 0:
            data = waveform[..., 0]  # First half
        else:
            try:
                data = waveform[..., 1]  # Second half
            except:
                data = waveform[..., 0]  # First half
        # Handle data being the wrong shape
        if data.shape != (60,30):
            print(f"File: {self.np_file_names[i]} was wrong shape: {data.shape}")
            data = np.zeros((60,30))
        return data

