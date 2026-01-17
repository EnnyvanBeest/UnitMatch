import os
from pathlib import Path
import random
import numpy as np
import h5py
from torch.utils.data import Dataset, Sampler
from utils.helpers import get_unit_id

class NeuropixelsDataset(Dataset):
    def __init__(self, data_dir:str, batch_size=1, mode='val'):
        """
        Initialises a dataset for testing or training.

        Args:
            data_dir: the root (absolute) directory under which the processed data can be found.
            batch_size: the min. number of units in a training/testing batch.
            mode: 'train' or 'val'
        """
        self.data_dir = Path(data_dir).resolve()

        self.batch_size = batch_size
        self.mode = mode
        self.experiment_unit_map = {}  # Maps experiment to its units

        print(self.data_dir, "is the data directory")

        for id, session in enumerate(os.listdir(self.data_dir)):
            self.experiment_unit_map[id] = self.select_good_units_files(os.path.join(self.data_dir, session), load_pre_merge=False)

        self.all_files = [(exp, file) for exp, files in self.experiment_unit_map.items() for file in files]
                    
        if len(self.all_files) < 1:
            print("No data in test dataset! Try a smaller batch size?")
        else:
            print(f"Initialised with {len(self.all_files)} files in the dataset.")

    def __len__(self):
        return len(self.all_files)

    def _normalize_waveform(self, waveform):
        # max-min normalization
        max_val = np.max(waveform)
        min_val = np.min(waveform)
        return (waveform - min_val) / (max_val - min_val)
        
    def __getitem__(self, i):
        experiment_path, neuron_file = self.all_files[i]
        with h5py.File(neuron_file, 'r') as f:
            waveform = f['waveform'][()] # waveform [T,C,2]
            MaxSitepos = f['MaxSitepos'][()]
        if waveform.shape != (60,30,2):
            waveform = np.zeros((60,30,2))
            # assert False, f"Waveform shape is not (60,30,2) but {waveform.shape}"
        ## do data augmentation 
        if self.mode == 'train':
            waveform_fh = self._augment_original(waveform[..., 0])
            waveform_sh = self._augment_original(waveform[..., 1])
        else:
            waveform_fh = waveform[..., 0]
            waveform_sh = waveform[..., 1]
            
        return waveform_fh, waveform_sh, MaxSitepos, experiment_path, neuron_file
    
    def select_good_units_files(self, directory, load_pre_merge:bool=True):
        """
        Selects the filenames of the good units based on the good_units_value array.
        Args:
        - directory (str): The directory containing the unit files.
        - good_units_value (list or numpy.ndarray): An array where a value of 1 indicates a good unit.
        - load_pre_merge (bool): Whether to load the pre-merge data.
        Returns:
        - list: A list of filenames corresponding to the good units.
        """
        files = os.listdir(directory)
        merges = {}
        removes = []
        indices = []
        for file in files:
            if load_pre_merge:
                if '+' in file:
                    removes.append(file)
            else:
                if '+' in file:
                    f = file.replace("Unit",'')
                    f = f.replace(".npy", '')
                    id1 = int(f[:f.find('+')])
                    id2 = int(f[f.find('+')+1:])
                    merges[id1] = id2
                if '#' in file:
                    f = file.replace("Unit",'')
                    f = f.replace(".npy", '')
                    id1 = int(f[:f.find('#')])
                    removes.append(id1)
            indices.append(get_unit_id(file))
        indices = list(set(indices))  # Remove duplicates
        good_units_files = []
        for index in indices:
            if index in merges.keys():
                filename = f"Unit{index}+{merges[index]}.npy"
            elif index in merges.values() or index in removes:
                # don't load a unit if we already loaded the unit it merged with
                # or if it's a unit we wanted to remove
                continue
            else:   # load the unit, ignoring the # if it is there (for pre-merge data)
                filename = f"Unit{index}.npy"
                withhash = f"Unit{index}#.npy"
            filepath = os.path.join(directory, filename)
            if os.path.exists(filepath):  # Check if file exists before adding
                good_units_files.append(filepath)
            elif load_pre_merge and os.path.exists(os.path.join(directory, withhash)):
                good_units_files.append(os.path.join(directory, withhash))
            else:
                print(f"Warning: Expected file {filepath} does not exist.")
        return good_units_files


class TrainExperimentBatchSampler(Sampler):
    def __init__(self, data_source, batch_size, shuffle=False):
        self.data_source = data_source
        self.batch_size = batch_size
        self.shuffle = shuffle
        self.experiment_batches = self._create_batches()

    def _create_batches(self):
        batches = []
        for experiment, unit_paths in self.data_source.experiment_unit_map.items():
            file_to_idx = {file: idx for idx, (exp, file) in enumerate(self.data_source.all_files) if exp == experiment}
            experiment_indices = [file_to_idx[file] for file in unit_paths]
            batches.append(experiment_indices)
        return batches

    def __iter__(self):
        iter_batches = []
        for experiment_indices in self.experiment_batches:
            if self.shuffle:
                random.shuffle(experiment_indices)
            for i in range(0, len(experiment_indices), self.batch_size):
                batch = experiment_indices[i:i + self.batch_size]
                # Check if the last batch is smaller than batch_size
                if len(batch) < self.batch_size:
                    # Resample additional items from the experiment_indices to fill the batch
                    shortfall = self.batch_size - len(batch)
                    additional_samples = random.choices(experiment_indices, k=shortfall)
                    batch.extend(additional_samples)
                iter_batches.append(batch)
        if self.shuffle:
            random.shuffle(iter_batches)
        return iter(iter_batches)

    def __len__(self):
        total_batches = sum((len(exp_indices) + self.batch_size - 1) // self.batch_size for exp_indices in self.experiment_batches)
        return total_batches

class ValidationExperimentBatchSampler(Sampler):
    """
    Creates one batch per experiment with all data points for validation.
    Optionally shuffles data within each experiment batch in each iteration.
    """
    def __init__(self, data_source, shuffle=False):
        self.data_source = data_source
        self.shuffle = shuffle
        self.experiment_batches = self._create_batches()
        print(f"No. of experiment batches: {len(self.experiment_batches)}")

    def _create_batches(self):
        batches = []
        for experiment, unit_paths in self.data_source.experiment_unit_map.items():
            # Create a mapping from file paths to indices
            file_to_idx = {file: idx for idx, (exp, file) in enumerate(self.data_source.all_files) if exp == experiment}
            experiment_indices = [file_to_idx[file] for file in unit_paths]
            # Each experiment is a single batch with all its units
            batches.append(experiment_indices)
        return batches

    def __iter__(self):
        iter_batches = []
        for experiment_indices in self.experiment_batches:
            # Shuffle the indices within each experiment if required
            if self.shuffle:
                random.shuffle(experiment_indices)
            iter_batches.append(experiment_indices)
        return iter(iter_batches)

    def __len__(self):
        return len(self.experiment_batches)
