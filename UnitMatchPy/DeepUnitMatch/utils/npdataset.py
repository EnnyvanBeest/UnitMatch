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


def roll_one_row(data, choice=None):
    """
    Training augmentation of one waveform snippet (shape [T, C]): shift it by
    one electrode row up or down, or leave it unchanged (each with p = 1/3).

    The snippet's channels interleave the probe's two columns (slot 2k and
    2k+1 are the two sites of row k, ordered by depth), so moving every
    channel by 2 slots moves the whole footprint one row along the probe;
    the row at the edge that has no neighbour keeps its own values. It is
    applied independently to the two halves of a neuron, so the two copies
    of a neuron can end up as much as two rows apart. This mimics the peak
    channel (on which the snippet is centred) being picked one row off.

    Fixed: the "up" shift used to skip the last channel of the second column
    (slot C-2 never received slot C's values), so that column's top row was
    not shifted.
    """
    if choice is None:
        choice = random.choice(["roll_up", "roll_down", "none"])
    if choice == "roll_up":
        data[:, :-2] = data[:, 2:].copy()  # slot i <- slot i+2; last row keeps its values
    elif choice == "roll_down":
        data[:, 2:] = data[:, :-2].copy()  # slot i <- slot i-2; first row keeps its values
    return data


class NeuropixelsDataset(Dataset):
    def __init__(self, save_path: str, batch_size=32, mode="val"):
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
            self.experiment_unit_map[id] = [
                os.path.join(full_session_path, file)
                for file in os.listdir(full_session_path)
            ]

        self.all_files = [
            (exp, file)
            for exp, files in self.experiment_unit_map.items()
            for file in files
        ]

        if len(self.all_files) < 1:
            print("No data in test dataset! Try a smaller batch size?")
        else:
            print(f"Initialised with {len(self.all_files)} files in the dataset.")

    def __len__(self):
        return len(self.all_files)

    def __getitem__(self, i):
        experiment_path, neuron_file = self.all_files[i]
        with h5py.File(neuron_file, "r") as f:
            waveform = f["waveform"][()]  # waveform [T,C,2]
            MaxSitepos = f["MaxSitepos"][()]
        if waveform.shape != (60, 30, 2):
            waveform = np.zeros((60, 30, 2))
            # assert False, f"Waveform shape is not (60,30,2) but {waveform.shape}"
        ## do data augmentation
        if self.mode == "train":
            waveform_fh = self._augment_original(waveform[..., 0])
            waveform_sh = self._augment_original(waveform[..., 1])
        else:
            waveform_fh = waveform[..., 0]
            waveform_sh = waveform[..., 1]

        return waveform_fh, waveform_sh, MaxSitepos, experiment_path, neuron_file

    def _augment_original(self, data):
        return roll_one_row(data)


class NeuropixelsDataset_cortexlab(Dataset):
    def __init__(
        self, data_dir: str, batch_size=1, mode="val", unit_order: str = "filesystem"
    ):
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
        sessions = sorted(
            sessions, key=lambda s: int(s) if str(s).isdigit() else str(s)
        )
        for id, session in enumerate(sessions):
            session_dir = os.path.join(self.data_dir, session)
            if self.unit_order == "unitmatch":
                unit_ids = _load_good_unit_ids_from_labels(session_dir)
                if unit_ids is None:
                    # Fallback: keep existing behavior if label files aren't present
                    self.experiment_unit_map[id] = self.select_good_units_files(
                        session_dir, load_pre_merge=False
                    )
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
                    raise ValueError(
                        "When unit_order is a list/tuple, it must provide unit IDs per session in order."
                    )
                ordered_files = []
                for unit_id in unit_ids:
                    fp = _unit_id_to_filepath(session_dir, int(unit_id))
                    if fp is not None:
                        ordered_files.append(fp)
                self.experiment_unit_map[id] = ordered_files
            else:
                self.experiment_unit_map[id] = self.select_good_units_files(
                    session_dir, load_pre_merge=False
                )

        self.all_files = [
            (exp, file)
            for exp, files in self.experiment_unit_map.items()
            for file in files
        ]

        if len(self.all_files) < 1:
            print("No data in test dataset! Try a smaller batch size?")
        else:
            print(f"Initialised with {len(self.all_files)} files in the dataset.")

    def __len__(self):
        return len(self.all_files)

    def __getitem__(self, i):
        experiment_path, neuron_file = self.all_files[i]
        with h5py.File(neuron_file, "r") as f:
            waveform = f["waveform"][()]  # waveform [T,C,2]
            MaxSitepos = f["MaxSitepos"][()]
        if waveform.shape != (60, 30, 2):
            waveform = np.zeros((60, 30, 2))
            # assert False, f"Waveform shape is not (60,30,2) but {waveform.shape}"
        ## do data augmentation
        if self.mode == "train":
            waveform_fh = self._augment_original(waveform[..., 0])
            waveform_sh = self._augment_original(waveform[..., 1])
        else:
            waveform_fh = waveform[..., 0]
            waveform_sh = waveform[..., 1]

        return waveform_fh, waveform_sh, MaxSitepos, experiment_path, neuron_file

    def _augment_original(self, data):
        return roll_one_row(data)

    def select_good_units_files(self, directory, load_pre_merge: bool = True):
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
                if "+" in file:
                    removes.append(file)
            else:
                if "+" in file:
                    f = file.replace("Unit", "")
                    f = f.replace("_RawSpikes.npy", "")
                    id1 = int(f[: f.find("+")])
                    id2 = int(f[f.find("+") + 1 :])
                    merges[id1] = id2
                if "#" in file:
                    f = file.replace("Unit", "")
                    f = f.replace("_RawSpikes.npy", "")
                    id1 = int(f[: f.find("#")])
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
            else:  # load the unit, ignoring the # if it is there (for pre-merge data)
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


def _session_indices(data_source, experiments=None):
    """{session: [dataset indices of its units]} (all sessions, or only `experiments`)."""
    wanted = None if experiments is None else set(experiments)
    file_to_idx = {
        (exp, file): idx for idx, (exp, file) in enumerate(data_source.all_files)
    }
    return {
        experiment: [file_to_idx[(experiment, file)] for file in unit_paths]
        for experiment, unit_paths in data_source.experiment_unit_map.items()
        if wanted is None or experiment in wanted
    }


class TrainExperimentBatchSampler(Sampler):
    """
    Training batches for the contrastive loss: every batch holds units of one
    recording session only (the loss contrasts neurons recorded together).

    Each session is split into the smallest number of batches of at most
    batch_size units, all of (nearly) equal size, e.g. 52 units -> 26 + 26,
    25 units -> one batch of 25. Every unit is used exactly once per epoch and
    never appears twice in a batch. (Previously the last batch of a session
    was padded to batch_size by drawing units of the same session with
    replacement, so a neuron could be its own negative.) Sessions with fewer
    than 2 units are skipped: they have no negatives.

    experiments: optional subset of sessions (e.g. the training split).
    """

    def __init__(self, data_source, batch_size, shuffle=False, experiments=None):
        self.data_source = data_source
        self.batch_size = batch_size
        self.shuffle = shuffle
        self.experiment_batches = [
            idx for idx in _session_indices(data_source, experiments).values() if len(idx) >= 2
        ]

    def _n_batches(self, n_units):
        return -(-n_units // self.batch_size)  # ceil

    def __iter__(self):
        iter_batches = []
        for experiment_indices in self.experiment_batches:
            indices = list(experiment_indices)
            if self.shuffle:
                random.shuffle(indices)
            for chunk in np.array_split(np.array(indices), self._n_batches(len(indices))):
                iter_batches.append(chunk.tolist())
        if self.shuffle:
            random.shuffle(iter_batches)
        return iter(iter_batches)

    def __len__(self):
        return sum(self._n_batches(len(idx)) for idx in self.experiment_batches)


class ValidationExperimentBatchSampler(Sampler):
    """
    Creates one batch per experiment with all data points for validation.
    Optionally shuffles data within each experiment batch in each iteration.

    experiments: optional subset of sessions (e.g. the held-out validation split).
    """

    def __init__(self, data_source, shuffle=False, experiments=None):
        self.data_source = data_source
        self.shuffle = shuffle
        self.experiment_batches = list(_session_indices(data_source, experiments).values())
        print(f"No. of experiment batches: {len(self.experiment_batches)}")

    def __iter__(self):
        iter_batches = []
        for experiment_indices in self.experiment_batches:
            indices = list(experiment_indices)
            # Shuffle the indices within each experiment if required
            if self.shuffle:
                random.shuffle(indices)
            iter_batches.append(indices)
        return iter(iter_batches)

    def __len__(self):
        return len(self.experiment_batches)


def split_sessions(data_source, val_fraction=0.05, seed=0, min_sessions=20):
    """
    Hold out whole recording sessions for validating the contrastive
    fine-tuning (batches are built per session, so units of one session must
    not be split between training and validation).

    Returns (train_sessions, val_sessions) as lists of experiment keys.
    A fixed seed makes the split reproducible. With fewer than min_sessions
    sessions nothing is held out (val_sessions empty).
    """
    sessions = sorted(data_source.experiment_unit_map)
    if len(sessions) < min_sessions or val_fraction <= 0:
        return sessions, []
    n_val = max(1, int(round(val_fraction * len(sessions))))
    val = sorted(random.Random(seed).sample(sessions, n_val))
    val_set = set(val)
    return [s for s in sessions if s not in val_set], val
