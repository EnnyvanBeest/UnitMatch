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
    session_switch = np.cumsum([np.asarray(units).squeeze().shape[0] for units in good_units])
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
            units_ses1 = np.asarray(good_units[i]).squeeze().astype(int)
            units_ses2 = np.asarray(good_units[j]).squeeze().astype(int)
            n_units_ses1 = units_ses1.shape[0]
            n_units_ses2 = units_ses2.shape[0]
             
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
            
def get_unit_id(filepath:str):
    fp = os.path.basename(filepath)
    if fp[:4] == "Unit" and fp[-14:] == "_RawSpikes.npy":
        fp = fp.replace("Unit", "")
        id = fp.replace("_RawSpikes.npy", "")
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
    output = {"ID": filenames, "x": x, "y": y}
    return pd.DataFrame(output).sort_values("ID").reset_index(drop=True)
