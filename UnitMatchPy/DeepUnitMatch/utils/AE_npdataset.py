import os
import random
import numpy as np
import h5py
from torch.utils.data import Dataset
from DeepUnitMatch.utils.train_utils import read_good_ids


class AE_NeuropixelsDataset(Dataset):
    """
    Dataset class for Autoencoder training on Neuropixels data.
    Each item is a waveform of shape (60, 30) corresponding to a single waveform
    from a good unit. The dataset randomly selects between the first and second
    half of the waveform data if available.

    Args:
        - save_path (str): Path where processed waveform files (snippets) are saved.
        - batch_size (int): Number of units to process in each batch.
    """

    def __init__(self, save_path: str, batch_size = 32):
        self.save_path = os.path.join(save_path, "processed_waveforms")
        self.batch_size = batch_size

        processed_waveform_sessions = os.listdir(self.save_path)
        
        self.np_file_names = []
        for directory in processed_waveform_sessions:
            snippet_path = os.path.join(self.save_path, directory)
            self.np_file_names.extend([os.path.join(snippet_path, filename) for filename in os.listdir(snippet_path)])

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

class AE_NeuropixelsDataset_cortexlab(Dataset):
    """
    The internal version of the dataset class. This works with our specific data structure:
    root/mouse/probe/location/experiment/processed_waveforms/*.npy
    root/mouse/probe/location/experiment/metadata.json
    The metadata.json file contains the good unit indices.
    
    This is left here in case user's data is structured this way, as it allows you to use the train_AE.py script's
    run() function from a command line/bash script. 
    """

    def __init__(self, root: str, batch_size = 32, mice = None):
        self.root = root
        if mice == 10:
            self.mouse_names = ['AL031', 'AL032', 'AL036', 'AV008', 'CB015', 'CB016', 'CB017', 'CB018', 'CB020', 'EB019']
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
        # file_name, half = self.np_file_names[i]
        # data = np.load(file_name)
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
    
