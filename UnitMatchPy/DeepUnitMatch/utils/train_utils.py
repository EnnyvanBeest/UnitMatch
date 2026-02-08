import os
import json


def read_good_ids(root, batch_size, mouse_names, finetune:bool):
    """
    Output depends on whether you want to pre-train the encoder (finetune=False)
    or finetune it via contrastive learning (finetune=True). This must be specified.
    """
    np_file_names = []
    experiment_unit_map = {}

    for name in mouse_names:
        name_path = os.path.join(root, name)
        probes = os.listdir(name_path)
        for probe in probes:
            probe_path = os.path.join(name_path, probe)
            locs = os.listdir(probe_path)
            for loc in locs:
                loc_path = os.path.join(probe_path, loc)
                experiments = os.listdir(loc_path)
                for experiment in experiments:
                    experiment_path = os.path.join(loc_path, experiment)
                    good_units_files = read_good_files(experiment_path, batch_size)
                    if good_units_files is None:
                        continue
                    if finetune:
                        experiment_unit_map[experiment_path] = good_units_files
                    else:
                        for file in good_units_files:
                            np_file_names.append(file)
    
    if finetune:
        return experiment_unit_map
    else:
        return np_file_names

def read_good_files(experiment_path, batch_size, load_pre_merge:bool=True):
    if not os.path.isdir(experiment_path):
        return None
    try:
        metadata_file = os.path.join(experiment_path, "metadata.json")
        metadata = json.load(open(metadata_file))
    except:
        print(experiment_path)
        raise ValueError("Did not find metadata.json file for this experiment")
    good_units_index = metadata["good_ids"]
    len_good_units = sum(good_units_index)
    if len_good_units < batch_size:
        return None
    good_units_files = select_good_units_files(os.path.join(experiment_path, 'processed_waveforms'), 
                                               good_units_index, load_pre_merge)
    return good_units_files

def select_good_units_files(directory, good_units_value, load_pre_merge:bool=True):
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
    good_units_files = []
    for index, is_good in enumerate(good_units_value):
        if is_good:
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

