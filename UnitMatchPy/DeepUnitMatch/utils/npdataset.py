import os
from pathlib import Path
import random
import numpy as np
import h5py
from torch.utils.data import Dataset, Sampler
from utils.helpers import get_unit_id


def _load_good_unit_ids_from_labels(session_dir: str):
    """
    Attempt to mirror UnitMatchPy.utils.load_good_waveforms() ordering:
    - Prefer BombCell labels (cluster_bc_unitType.tsv): keep 'GOOD' and 'NON-SOMA GOOD'
    - Else fall back to KiloSort labels (cluster_group.tsv): keep 'good'
    - Preserve file order (no sorting of IDs)
    Returns list[int] of good unit IDs, or None if no label file found.
    """
    label_candidates = [
        "cluster_bc_unitType.tsv",
        "cluster_group.tsv",
    ]
    label_path = None
    for name in label_candidates:
        candidate = os.path.join(session_dir, name)
        if os.path.exists(candidate):
            label_path = candidate
            break
    if label_path is None:
        return None

    good_ids = []
    with open(label_path, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) < 2:
                continue

            # Skip header rows like "cluster_id\tgroup"
            try:
                unit_id = int(parts[0])
            except ValueError:
                continue

            label = str(parts[1]).strip().lower()
            if os.path.basename(label_path) == "cluster_bc_unitType.tsv":
                if label in {"good", "non-soma good"}:
                    good_ids.append(unit_id)
            else:
                if label == "good":
                    good_ids.append(unit_id)
    return good_ids


def _unit_id_to_filepath(session_dir: str, unit_id: int):
    """
    Find a per-unit file for a given unit_id.
    Prefer the exact UnitMatch-style filename, but fall back to any matching prefix.
    """
    exact = os.path.join(session_dir, f"Unit{unit_id}_RawSpikes.npy")
    if os.path.exists(exact):
        return exact

    prefix = f"Unit{unit_id}"
    matches = sorted(
        f
        for f in os.listdir(session_dir)
        if (
            f.startswith(prefix)
            and f.endswith("_RawSpikes.npy")
            and (not os.path.isdir(os.path.join(session_dir, f)))
        )
    )
    if matches:
        # Prefer non-removed/non-merged filenames if multiple exist.
        matches = sorted(matches, key=lambda name: ("#" in name, "+" in name, name))
        return os.path.join(session_dir, matches[0])
    return None


class NeuropixelsDataset(Dataset):
    def __init__(self, save_path: str, batch_size = 32, mode='val'):
        """
        Initialises a dataset for testing or training.

        Args:
            save_path: the directory under which the processed data can be found.
            batch_size: the min. number of units in a training/testing batch.
            mode: 'train' or 'val'
        """

        self.save_path = os.path.join(save_path, "processed_waveforms")
        self.batch_size = batch_size
        self.mode = mode
        
        self.experiment_unit_map = {}

        for id, session in enumerate(os.listdir(self.save_path)):
            full_session_path = os.path.join(self.save_path, session)
            self.experiment_unit_map[id] = [os.path.join(full_session_path, file) for file in os.listdir(full_session_path)]

        self.all_files = [(exp, file) for exp, files in self.experiment_unit_map.items() for file in files]

        if len(self.all_files) < 1:
            print("No data in test dataset! Try a smaller batch size?")
        else:
            print(f"Initialised with {len(self.all_files)} files in the dataset.")

    def __len__(self):
        return len(self.all_files)

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

    def _augment_original(self, data):
        # Apply random augmentations to data, shape [T,C]
        roll_choice = random.choice(["roll_up", "roll_down", "none"])

        if roll_choice == "roll_up":
            # implement roll_up augmentation
            C = data.shape[1]  # Number of channels
            # Indices for odd channels, excluding the last one if C is odd
            odd_indices = np.arange(0, C - 1, 2)
            # Indices for even channels, excluding the last one
            even_indices = np.arange(1, C - 1, 2)
            # Shift odd channels up, excluding the last odd channel
            if len(odd_indices) > 1:  # Check if there are at least 2 odd channels to roll
                data[:, odd_indices[:-1]] = data[:, odd_indices[1:]]
            # Shift even channels up, excluding the last even channel
            if len(even_indices) > 1:  # Check if there are at least 2 even channels to roll
                data[:, even_indices[:-1]] = data[:, even_indices[1:]]
        elif roll_choice == "roll_down":
            # implement roll_down augmentation
            C = data.shape[1]  # Number of channels
            odd_indices = np.arange(2, C, 2)
            even_indices = np.arange(3, C, 2)
            if len(odd_indices) > 0:  # Check if there are odd channels to roll
                data[:, odd_indices] = data[:, odd_indices - 2]
            if len(even_indices) > 0:  # Check if there are even channels to roll
                data[:, even_indices] = data[:, even_indices - 2]

        return data
    

class NeuropixelsDataset_cortexlab(Dataset):
    def __init__(self, data_dir: str, batch_size=1, mode="val", unit_order: str = "filesystem"):
        """
        Initialises a dataset for testing or training.

        Args:
            data_dir: the root (absolute) directory under which the processed data can be found.
            batch_size: the min. number of units in a training/testing batch.
            mode: 'train' or 'val'
            unit_order: 'filesystem' (default) or 'unitmatch' to mirror UnitMatch TSV order.
        """
        self.data_dir = Path(data_dir).resolve()

        self.batch_size = batch_size
        self.mode = mode
        self.unit_order = unit_order
        self.experiment_unit_map = {}  # Maps experiment to its units

        print(self.data_dir, "is the data directory")

        sessions = list(os.listdir(self.data_dir))
        # Deterministic session ordering. Prefer numeric ordering when folder names are integers.
        sessions = sorted(sessions, key=lambda s: int(s) if str(s).isdigit() else str(s))
        for id, session in enumerate(sessions):
            session_dir = os.path.join(self.data_dir, session)
            if self.unit_order == "unitmatch":
                unit_ids = _load_good_unit_ids_from_labels(session_dir)
                if unit_ids is None:
                    # Fallback: keep existing behavior if label files aren't present
                    self.experiment_unit_map[id] = self.select_good_units_files(session_dir, load_pre_merge=False)
                else:
                    ordered_files = []
                    for unit_id in unit_ids:
                        fp = _unit_id_to_filepath(session_dir, unit_id)
                        if fp is not None:
                            ordered_files.append(fp)
                    self.experiment_unit_map[id] = ordered_files
            elif isinstance(self.unit_order, (list, tuple)):
                # Explicit ordering: unit_order is a list (per experiment) of unit IDs in the desired order.
                try:
                    unit_ids = list(self.unit_order[id])
                except Exception:
                    raise ValueError("When unit_order is a list/tuple, it must provide unit IDs per session in order.")
                ordered_files = []
                for unit_id in unit_ids:
                    fp = _unit_id_to_filepath(session_dir, int(unit_id))
                    if fp is not None:
                        ordered_files.append(fp)
                self.experiment_unit_map[id] = ordered_files
            else:
                self.experiment_unit_map[id] = self.select_good_units_files(session_dir, load_pre_merge=False)

        self.all_files = [(exp, file) for exp, files in self.experiment_unit_map.items() for file in files]
                     
        if len(self.all_files) < 1:
            print("No data in test dataset! Try a smaller batch size?")
        else:
            print(f"Initialised with {len(self.all_files)} files in the dataset.")

    def __len__(self):
        return len(self.all_files)

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
    
    def _augment_original(self, data):
        # Apply random augmentations to data, shape [T,C]
        roll_choice = random.choice(["roll_up", "roll_down", "none"])

        if roll_choice == "roll_up":
            # implement roll_up augmentation
            C = data.shape[1]  # Number of channels
            # Indices for odd channels, excluding the last one if C is odd
            odd_indices = np.arange(0, C - 1, 2)
            # Indices for even channels, excluding the last one
            even_indices = np.arange(1, C - 1, 2)
            # Shift odd channels up, excluding the last odd channel
            if len(odd_indices) > 1:  # Check if there are at least 2 odd channels to roll
                data[:, odd_indices[:-1]] = data[:, odd_indices[1:]]
            # Shift even channels up, excluding the last even channel
            if len(even_indices) > 1:  # Check if there are at least 2 even channels to roll
                data[:, even_indices[:-1]] = data[:, even_indices[1:]]
        elif roll_choice == "roll_down":
            # implement roll_down augmentation
            C = data.shape[1]  # Number of channels
            odd_indices = np.arange(2, C, 2)
            even_indices = np.arange(3, C, 2)
            if len(odd_indices) > 0:  # Check if there are odd channels to roll
                data[:, odd_indices] = data[:, odd_indices - 2]
            if len(even_indices) > 0:  # Check if there are even channels to roll
                data[:, even_indices] = data[:, even_indices - 2]

        return data
    
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
        files = sorted(os.listdir(directory))
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
                    f = f.replace("_RawSpikes.npy", '')
                    id1 = int(f[:f.find('+')])
                    id2 = int(f[f.find('+')+1:])
                    merges[id1] = id2
                if '#' in file:
                    f = file.replace("Unit",'')
                    f = f.replace("_RawSpikes.npy", '')
                    id1 = int(f[:f.find('#')])
                    removes.append(id1)
            indices.append(get_unit_id(file))
        indices = sorted(set(indices))  # Remove duplicates, deterministic order
        good_units_files = []
        for index in indices:
            if index in merges.keys():
                filename = f"Unit{index}+{merges[index]}_RawSpikes.npy"
            elif index in merges.values() or index in removes:
                # don't load a unit if we already loaded the unit it merged with
                # or if it's a unit we wanted to remove
                continue
            else:   # load the unit, ignoring the # if it is there (for pre-merge data)
                filename = f"Unit{index}_RawSpikes.npy"
                withhash = f"Unit{index}#_RawSpikes.npy"
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
