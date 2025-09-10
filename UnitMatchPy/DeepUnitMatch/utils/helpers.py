import pandas as pd
import numpy as np
import os
import h5py


def create_dataframe(good_units, prob_matrix, session_list=None):
    # Create dataframe for all session pairs

    # Get number of sessions and units per session
    if session_list is None:
        n_sessions = len(good_units)
        session_list = list(range(n_sessions))
    session_switch = np.cumsum([len(units.squeeze()) for units in good_units])
    session_switch = np.insert(session_switch, 0, 0)    

    d = {
        "RecSes1": [],
        "RecSes2": [],
        "ID1": [],
        "ID2": [],
        "Prob": [],
    }
    
    # Generate all pairs of sessions (including within-session pairs)
    for i, ses1 in enumerate(session_list):
        for j, ses2 in enumerate(session_list):
            n_units_ses1 = len(good_units[i].squeeze())
            n_units_ses2 = len(good_units[j].squeeze())
            
            # Get unit IDs for this session
            units_ses1 = np.arange(n_units_ses1)
            units_ses2 = np.arange(n_units_ses2)
            
            # Extract probability matrix block for this session pair
            prob_block = prob_matrix[session_switch[i]:session_switch[i+1], 
                                   session_switch[j]:session_switch[j+1]]
            
            # Add all unit pairs from this session pair
            d["RecSes1"].extend([ses1] * (n_units_ses1 * n_units_ses2))
            d["RecSes2"].extend([ses2] * (n_units_ses1 * n_units_ses2))
            d["ID1"].extend(np.repeat(units_ses1, n_units_ses2).tolist())
            d["ID2"].extend(np.tile(units_ses2, n_units_ses1).tolist())
            d["Prob"].extend(prob_block.ravel().tolist())

    df = pd.DataFrame(d)
    return df
    

def get_unit_ids(filepaths: list) -> list:
    # Convert to pandas Series for vectorized operations
    fps = pd.Series([os.path.basename(fp) for fp in filepaths])
    
    # Check format validity
    valid_mask = (fps.str[:4] == "Unit") & (fps.str[-4:] == ".npy")
    
    if not valid_mask.all():
        invalid_files = [filepaths[i] for i in valid_mask[~valid_mask].index]
        raise ValueError(f"Invalid filepath format for files: {invalid_files}")
    
    # Vectorized string operations
    ids = (fps.str.replace("Unit", "", regex=False)
              .str.replace(".npy", "", regex=False)
              .str.split('+').str[0]
              .str.replace("#", "", regex=False))
    
    # Convert to integers with error handling
    try:
        return ids.astype(int).tolist()
    except ValueError as e:
        # Find problematic entries
        for i, id_str in enumerate(ids):
            try:
                int(id_str)
            except ValueError:
                print(f"Problematic ID: {id_str}")
                raise ValueError(f"Invalid filepath format for waveform: {filepaths[i]}")
            
def get_unit_id(filepath:str):
    fp = os.path.basename(filepath)
    if fp[:4] == "Unit" and fp[-4:] == ".npy":
        fp = fp.replace("Unit", "")
        id = fp.replace(".npy", "")
        if '+' in id:
            id = id[:id.find('+')]
        if '#' in id:
            id = id.replace("#", "")
        try:
            return int(id)
        except:
            print(id)
            raise ValueError(f"Invalid filepath format for this waveform: {filepath}")
    else:
        raise ValueError(f"Invalid filepath format for this waveform: {filepath}")

def read_pos(path, skip_removed_units=True):
    files = os.listdir(path)
    x = []
    y = []
    filenames = []

    for file in files:
        if skip_removed_units:
            if '+' in file or '#' in file:
                continue
        fp = os.path.join(path, file)
        with h5py.File(fp, 'r') as f:
            # waveform = f['waveform'][()] 
            MaxSitepos = f['MaxSitepos'][()]
        x.append(MaxSitepos[0])
        y.append(MaxSitepos[1])
        filenames.append(get_unit_id(file))
    output = {"ID":filenames, "x":x, "y":y}

    return pd.DataFrame(output)


