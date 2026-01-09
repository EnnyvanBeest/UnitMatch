import numpy as np
import pandas as pd
import os, sys
import h5py
sys.path.insert(0, os.getcwd())
sys.path.insert(0, os.path.dirname(os.getcwd()))
sys.path.insert(0, os.path.join(os.path.dirname(os.getcwd()), "DeepUnitMatch"))
from testing.test import inference, get_threshold, directional_filter_df
from utils.helpers import create_dataframe


def merge_split_unit(row, spk_df:pd.DataFrame, data_dir):
    """
    This function will do the following to a pair of units that we decide are a split:
        - Combine their spike trains in spk_df.
        - Compute a weighted average of their waveforms and save it.
    """
    # merge the spike trains in spk_df
    clus = spk_df.loc[spk_df["recses"]==row["RecSes1"], "clusters"].item()
    times = spk_df.loc[spk_df["recses"]==row["RecSes1"], "times"].item()
    merged = times[np.where((clus==row['ID1']) | (clus==row['ID2']))]
    clus[np.where(clus==row["ID2"])] = row['ID1']               # reassign Unit 2 spikes to Unit 1

    # save a new waveform which is the weighted average of the two waveforms we are merging
    weight1 = len(clus[np.where((clus==row['ID1']))]) / len(merged)
    weight2 = len(clus[np.where((clus==row['ID2']))]) / len(merged)
    
    exp_id = spk_df.loc[spk_df['recses']==row['RecSes1'], 'recses'].item()
    waveforms_dir = os.path.join(data_dir, str(exp_id-1))
    wave1 = os.path.join(waveforms_dir, f"Unit{str(int(row['ID1']))}.npy")
    wave2 = os.path.join(waveforms_dir, f"Unit{str(int(row['ID2']))}.npy")
    with h5py.File(wave1, 'r') as f1:
        waveform1 = f1['waveform'][()]
        msp1 = f1['MaxSitepos'][()]
    with h5py.File(wave2, 'r') as f2:
        waveform2 = f2['waveform'][()]
    new_waveform = weight1 * waveform1 + weight2 * waveform2
    save_path = os.path.join(waveforms_dir, f"Unit{str(int(row['ID1']))}+{str(int(row['ID2']))}.npy")
    new_data = {
        "waveform": new_waveform,          # (60,30,2) 
        "MaxSitepos": msp1
    }
    with h5py.File(save_path, 'w') as f:
        for key, value in new_data.items():
            f.create_dataset(key, data=value)
    return spk_df

def load_spk_data(KSdirs):
    # Load the spike trains for all sessions

    data = {
        "recses" : [],
        "clusters" : [],
        "times" : [],
        "paths" : []
    }
    for idx, ses in enumerate(KSdirs):
        clusters = np.load(os.path.join(ses, "spike_clusters.npy"))
        times = np.load(os.path.join(ses, "spike_times.npy")) / 3e4  # convert to seconds

        data['recses'].append(idx)
        data['clusters'].append(clusters)
        data['times'].append(times)
        data['paths'].append(ses)
    return pd.DataFrame(data)

def estimate_C(spike_train, t_r, t_c):
    n_v = 0
    for i in range(spike_train.shape[0]):
        spike1 = spike_train[i]
        for spike2 in spike_train[i+1:]:
            diff = spike2-spike1
            if diff < (t_r - t_c):
                if diff > t_c:
                    n_v += 1
            else:
                break
    N = len(spike_train)
    T = max(spike_train) - min(spike_train) - 2*N*t_c
    if (1 - (2 * n_v * T)/(N**2 * t_r)) > 0:
        C = 0.5 * (1 - np.sqrt(1 - (2 * n_v * T)/(N**2 * t_r)))
        if C == 0:
            C = 0.01
    else:
        C = 1
    return C

def C_ratios(within, threshold_dnn, spk_df, null=False, t_r=0.002, t_c=0.0001):

    # Get the putative set of split units as off-diagonal matches that pass the thresholds
    matches_within = within.loc[within["Prob"]>=threshold_dnn]
    l = len(matches_within)
    if null:
        matches_within = within.loc[within["Prob"]<threshold_dnn].nsmallest(3*l)
        matches_within = directional_filter_df(matches_within)
        matches_within = matches_within.nlargest(l, ['dist'])
    else:
        matches_within = directional_filter_df(matches_within)
    putative_splits = matches_within.loc[matches_within["ID1"]!=matches_within["ID2"]]
    putative_splits_1dir = putative_splits.loc[putative_splits["ID1"]<putative_splits["ID2"]]

    ratios = {}
    for idx, row in putative_splits_1dir.iterrows():
        clus = spk_df.loc[spk_df["recses"]==row["RecSes1"], "clusters"].item()
        times = spk_df.loc[spk_df["recses"]==row["RecSes1"], "times"].item()
        if len(times[np.where((clus==row['ID1']))]) > len(times[np.where((clus==row['ID2']))]):
            non_merged = times[np.where((clus==row['ID1']))]
        else:
            non_merged = times[np.where((clus==row['ID2']))]
        merged = times[np.where((clus==row['ID1']) | (clus==row['ID2']))]
        C = estimate_C(merged, t_r, t_c)
        C_pm = estimate_C(non_merged, t_r, t_c)
        ratios[idx] = C / C_pm
    return ratios

def merge_and_remove_splits(param, prob_matrix, session_id, model, data_dir, max_C_ratio=1):
    """
    This function performs the entire pipeline for merging/removing split units.
    Requires that inference has already been done pre-merging.
    Output is the following:
        - Saves new waveforms for merged units (weighted average of the two)
        - Saves new spike data in merged_spikes.npy. Does not overwrite spikes.npy.
        - Does NN inference and functional score computation again taking merges into account.
    """

    spk_df = load_spk_data(param['KS_dirs'])
    to_merge = set()                        # match table indices corresponding to rows to merge
    to_remove = set()                       # unique unit identifiers to remove

    n_sessions = len(param['good_units'])
    session_switch = np.cumsum([len(units.squeeze()) for units in param['good_units']])
    session_switch = np.insert(session_switch, 0, 0)    

    df = create_dataframe(param['good_units'], prob_matrix)

    df['uid1'] = df['RecSes1'] * 1e6 + df['ID1']
    df['uid2'] = df['RecSes2'] * 1e6 + df['ID2']

    thresh = get_threshold(prob_matrix, session_id, True)

    within = df.loc[(df["RecSes1"]==df["RecSes2"])]
    across = df.loc[(df["RecSes1"]!=df["RecSes2"])]
    # Correct for different median similarities between within- and across-day sets.
    diff = np.median(within["Prob"]) - np.median(across["Prob"])
    thresh = thresh - diff

    match_C_ratios = C_ratios(within, thresh, spk_df, null=False)
    merges = [(idx, df.loc[idx,"uid1"].item(), df.loc[idx,"uid2"].item(), ratio)
                for idx,ratio in match_C_ratios.items() if ratio < max_C_ratio]
    merges = sorted(merges, key=lambda x:x[-1])     # sort merges by C ratio
    
    # Now we have a list of potential merges, we need to make sure we don't merge a unit more than once.
    ids_being_merged = set()
    for merge in merges:
        if merge[2] in ids_being_merged or merge[1] in ids_being_merged:
            # do not merge this pair as at least one of them is already being merged.
            continue
        else:
            ids_being_merged.update([merge[2], merge[1]])
            to_remove.add(merge[2])                     # remove the larger ID from this pair
            to_merge.add(merge[0])                      # merge this pair

    print(f'Units to remove: {len(to_remove)}')
    print(f'Unit to merge: {len(to_merge)}')

    if len(to_remove) == 0 and len(to_merge) == 0:
        print("No units to merge or remove. Exiting.")
        return prob_matrix, param, session_id, session_switch
    
    df_to_merge = df.loc[list(to_merge),:]
    for __, row in df_to_merge.iterrows():
        spk_df = merge_split_unit(row, spk_df, data_dir)
    for __, row in spk_df.iterrows():
        data = {
            "spkclus": spk_df.loc[spk_df['recses']==row['recses'], "clusters"].item(),
            "spktimes": spk_df.loc[spk_df['recses']==row['recses'], "times"].item()
        }
        old_spikes_path = spk_df.loc[spk_df['recses']==row['recses'], "paths"].item()
        new_path = os.path.join(os.path.dirname(old_spikes_path), "merged_spikes.npy")
        with h5py.File(new_path, 'w') as f:
            for key, value in data.items():
                f.create_dataset(key, data=value)
                
    # Handle removals for multiple sessions
    sessions = [int(i//1e6) for i in to_remove]
    ids = [int(i%1e6) for i in to_remove]
    removals = {i: [] for i in range(n_sessions)}  # Initialize removals dict for all sessions
    
    for s, id in zip(sessions, ids):
        dir = os.path.join(data_dir, str(s-1))
        file = os.path.join(dir, f"Unit{id}_RawSpikes.npy")
        newfile = os.path.join(dir, f"Unit{id}#_RawSpikes.npy")
        os.rename(file, newfile)
        removals[s].append(id)

    # Update good_units for all sessions
    if len(to_remove) > 0:
        for ses in range(n_sessions):
            if removals[ses]:  # If there are removals in this session
                # For sessions other than the first, we need to adjust the indices
                if ses == 0:
                    param['good_units'][ses] = np.delete(param['good_units'][ses].squeeze(), removals[ses])
                else:
                    # Calculate offset for this session
                    offset = session_switch[ses]
                    adjusted_removals = [r - offset for r in removals[ses] if r >= offset]
                    if adjusted_removals:
                        param['good_units'][ses] = np.delete(param['good_units'][ses].squeeze(), adjusted_removals)
                    else:
                        param['good_units'][ses] = param['good_units'][ses].squeeze()
            else:
                param['good_units'][ses] = param['good_units'][ses].squeeze()
    else:
        for ses in range(n_sessions):
            param['good_units'][ses] = param['good_units'][ses].squeeze()

    # Recalculate session information
    param['n_units'] = sum(len(param['good_units'][ses]) for ses in range(n_sessions))
    session_id = np.concatenate([np.full(len(param['good_units'][ses]), ses) for ses in range(n_sessions)]).astype(int)
    new_dnnsim = inference(model, data_dir)
    
    # Update session_switch based on new good_units lengths
    new_session_switch = np.cumsum([len(param['good_units'][ses]) for ses in range(n_sessions)])
    new_session_switch = np.insert(new_session_switch, 0, 0)
    
    return new_dnnsim, param, session_id, new_session_switch

def undo_merge(data_dir):
    """
    Undoes all changes to waveform file names made by the merge_and_remove_splits function.
    Also resets the param object.
    """

    for exp in os.listdir(data_dir):
        waveforms_path = os.path.join(data_dir, exp)
        for file in os.listdir(waveforms_path):
            if "#" in file:
                new_filename = file.replace("#", "")
                os.rename(os.path.join(waveforms_path, file), 
                            os.path.join(waveforms_path, new_filename))
            if "+" in file:
                os.remove(os.path.join(waveforms_path, file))


if __name__ == "__main__":

    pass