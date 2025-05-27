import numpy as np
import pandas as pd
import os, sys
import h5py, tqdm, mat73
sys.path.insert(0, os.getcwd())
sys.path.insert(0, os.path.join(os.getcwd(), "testing"))
from utils.myutil import get_exp_id, mtpath_to_expids, get_threshold
from testing.test import inference, directional_filter, spatial_filter
from testing.compute_func_scores import post_merge_ISI, compute_rpc


def merge_split_unit(mt_path, row, spk_df:pd.DataFrame):
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
    exp_id = spk_df.loc[spk_df['recses']==row['RecSes1'], 'exp_id'].item()
    waveforms_dir = os.path.join(os.path.dirname(mt_path), exp_id, 'processed_waveforms')
    wave1 = os.path.join(waveforms_dir, f"Unit{str(int(row['ID1']))}_RawSpikes.npy")
    wave2 = os.path.join(waveforms_dir, f"Unit{str(int(row['ID2']))}_RawSpikes.npy")
    with h5py.File(wave1, 'r') as f1:
        waveform1 = f1['waveform'][()]
        msp1 = f1['MaxSitepos'][()]
    with h5py.File(wave2, 'r') as f2:
        waveform2 = f2['waveform'][()]
    new_waveform = weight1 * waveform1 + weight2 * waveform2
    save_path = os.path.join(waveforms_dir, f"Unit{str(int(row['ID1']))}+{str(int(row['ID2']))}_RawSpikes.npy")
    new_data = {
        "waveform": new_waveform,          # (60,30,2) 
        "MaxSitepos": msp1
    }
    with h5py.File(save_path, 'w') as f:
        for key, value in new_data.items():
            f.create_dataset(key, data=value)
    return spk_df

def load_spk_data(mt_path, mt):
    # Load the spike trains for all the sessions in mt
    exp_ids, metadata = mtpath_to_expids(mt_path, mt)
    sessions = mt["RecSes1"].unique()
    data = {
        "recses" : [],
        "exp_id" : [],
        "clusters" : [],
        "times" : [],
        "paths" : []
    }
    for ses in sessions:
        spk_path = os.path.join(os.path.dirname(mt_path), exp_ids[ses], "spikes.npy")
        with h5py.File(spk_path, 'r') as f:
            clusters = f['spkclus'][()]
            times = f['spktimes'][()]
        data['recses'].append(ses)
        data['exp_id'].append(exp_ids[ses])
        data['clusters'].append(clusters)
        data['times'].append(times)
        data['paths'].append(spk_path)
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

def C_ratios(within, threshold_dnn, threshold_um, mt_path, spk_df, null=False, t_r=0.002, t_c=0.0001):

    # Get the putative set of split units as off-diagonal matches that pass the thresholds
    matches_within = within.loc[within["DNNSim"]>=threshold_dnn, ["ID1", "ID2", "RecSes1", "RecSes2", "DNNSim", "MatchProb"]]
    matches_within = matches_within.loc[matches_within["MatchProb"]>=threshold_um, ["ID1", "ID2", "RecSes1", "RecSes2", "DNNSim", "MatchProb"]]
    l = len(matches_within)
    if null:
        matches_within = within.loc[within["DNNSim"]<threshold_dnn, ["ID1", "ID2", "RecSes1", "RecSes2", "DNNSim", "MatchProb"]].nsmallest(3*l, ["DNNSim"])
        matches_within = matches_within.loc[matches_within["MatchProb"] < 0.5]
        matches_within = directional_filter(matches_within)
        matches_within = spatial_filter(mt_path, matches_within, 20000)     # just to get the distance column in the df
        matches_within = matches_within.nlargest(l, ['dist'])
    else:
        matches_within = directional_filter(matches_within)
        matches_within = spatial_filter(mt_path, matches_within, 20)
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

def merge_and_remove_splits(test_data_root, max_C_ratio=1, mt_name="new_matchtable.csv", 
                            dnn_metric="DNNSim", um_metric="MatchProb"):
    """
    This function performs the entire pipeline for merging/removing split units across all mice.
    Requires that inference has already been done pre-merging (so that DNNSim values are in the match tables).
    Output is the following:
        - Saves new waveforms for merged units (weighted average of the two)
        - Saves new spike data in merged_spikes.npy. Does not overwrite spikes.npy.
        - Saves a new match table with merges taken into account called merged_mt.csv.
    
    Args:
        - test_data_root: root directory containing all the test data with the split units we want to merge.
        - max_C_ratio: Defaults to 1. Increasing this means more splits are merged. Reducing it means fewer splits are merged.
        - mt_name: the name of the match table containing the split units in each mouse/probe/loc directory.
        - dnn_metric: usually DNNSim, but can also be DNNProb.
        - um_metric: usually MatchProb, but can also be some other metric from UnitMatch.
    """
    server_root = r"\\znas\Lab\Share\UNITMATCHTABLES_ENNY_CELIAN_JULIE\FullAnimal_KSChanMap"

    for mouse in tqdm(os.listdir(test_data_root)):
        mouse_path = os.path.join(test_data_root, mouse)
        probes = os.listdir(mouse_path)
        for probe in probes:
            probe_path = os.path.join(mouse_path, probe)
            locs = os.listdir(probe_path)
            for loc in locs:
                mt_path = os.path.join(test_data_root, mouse, probe, loc, mt_name)
                if not os.path.exists(mt_path):
                    continue
                # if os.path.exists(os.path.join(test_data_root, mouse, probe, loc, "merged_mt.csv")):
                #     continue
                um_path = os.path.join(server_root, mouse, probe, loc, "UnitMatch", "UnitMatch.mat")
                um = mat73.loadmat(um_path)
                path_list = um["UMparam"]["KSDir"]
                path_dict = {}
                for idx, path in enumerate(path_list):
                    p = get_exp_id(path, mouse)
                    path_dict[int(idx+1)] = p
                try:
                    mt = pd.read_csv(mt_path)
                except:
                    mt = redo_inference(test_data_root, mouse, probe, loc)
                if len(mt) < 400:
                    # Need at least 20 functional matches, so 20^2 is min mt length
                    continue
                sessions = mt["RecSes1"].unique()
                mt['uid1'] = mt['RecSes1'] * 1e6 + mt['ID1']
                mt['uid2'] = mt['RecSes2'] * 1e6 + mt['ID2']
                spk_df = load_spk_data(mt_path, mt)
                to_merge = set()                        # match table indices corresponding to rows to merge
                to_remove = set()                       # unique (to this mt) unit identifiers to remove

                # We need to loop through pairs of sessions in order to define a threshold for matches
                for r1 in sessions:
                    for r2 in sessions:
                        if r1 >= r2 or abs(r2-r1)>1:
                            continue
                        df = mt.loc[(mt["RecSes1"].isin([r1,r2])) & (mt["RecSes2"].isin([r1,r2])),:]
                        if len(df) < 10:
                            continue
                        try:
                            thresh = get_threshold(df, metric=dnn_metric, vis=False)
                        except:
                            # some match tables contain NaNs in the DNNSim column. We need to redo inference for these.
                            mt = redo_inference(test_data_root, mouse, probe, loc)
                            df = mt.loc[(mt["RecSes1"].isin([r1,r2])) & (mt["RecSes2"].isin([r1,r2])),:]
                            thresh = get_threshold(df, metric=dnn_metric, vis=False)
                        within = df.loc[(df["RecSes1"]==df["RecSes2"]), [dnn_metric, um_metric, "ID1", "ID2", "RecSes1", "RecSes2", "uid1", "uid2"]]
                        across = df.loc[(df["RecSes1"]!=df["RecSes2"]), [dnn_metric, um_metric, "ID1", "ID2", "RecSes1", "RecSes2", "uid1", "uid2"]]
                        # Correct for different median similarities between within- and across-day sets.
                        diff = np.median(within[dnn_metric]) - np.median(across[dnn_metric])
                        thresh = thresh - diff
                        diff_um = np.median(within[um_metric]) - np.median(across[um_metric])
                        thresh_um = 0.5 - diff_um
                        # Now we have thresholds, we discard the units from day 2 otherwise things get complicated
                        # We want to attempt merging in each session exactly once.
                        if r2 != sessions.max():
                            match_C_ratios = C_ratios(within.loc[within['RecSes1']==r1], thresh, 
                                                  thresh_um, mt_path, spk_df, null=False)
                        else:
                            # The exception is for the last pair, where we want to do both days
                            # Otherwise, we would skip the last session altogether
                            match_C_ratios = C_ratios(within, thresh, thresh_um, mt_path, spk_df, null=False)
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
                
                df_to_merge = mt.loc[list(to_merge),:]
                for __, row in df_to_merge.iterrows():
                    spk_df = merge_split_unit(mt_path, row, spk_df)
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
                sessions = [int(i//1e6) for i in to_remove]
                ids = [int(i%1e6) for i in to_remove]
                for s, id in zip(sessions, ids):
                    exp_id = path_dict[s]
                    dir = os.path.join(os.path.dirname(mt_path), exp_id, "processed_waveforms")
                    file = os.path.join(dir, f"Unit{id}_RawSpikes.npy")
                    newfile = os.path.join(dir, f"Unit{id}#_RawSpikes.npy")
                    os.rename(file, newfile)

                merged_mt = mt[~mt[["uid1","uid2"]].apply(lambda x: x.isin(list(to_remove)), axis=1).any(axis=1)]
                merged_mt = merged_mt.drop(["uid1","uid2"], axis=1)
                merged_mt = merged_mt.loc[:, ~merged_mt.columns.str.contains('Unnamed')]
                merged_mt = update_mt(test_data_root, mouse, probe, loc, merged_mt)
                mt_save_path = os.path.join(os.path.dirname(mt_path), "merged_mt.csv")
                merged_mt.to_csv(mt_save_path, index=False)

def update_mt(test_data_root:str, mouse:str, probe:str, loc:str, merged_mt:pd.DataFrame):
    """
    Input is the merged match table. This just has the removed units removed.
    Output is the merged match table (same size as input) with updated DNNSim, refPopCorr and newISI columns.
    These need to be updated after merging. 
    DNNSim and newISI are updated in-place, whereas newRPC is added as a new column.
    """

    # recompute DNNSim
    merged_mt = inference(test_data_root, mouse, probe, loc, merged=True, merged_mt=merged_mt, model_name="incl_AV008")

    # recompute ISI
    mt_path = os.path.join(test_data_root, mouse, probe, loc, "merged_mt.csv")
    merged_mt = post_merge_ISI(merged_mt, mt_path)

    # recompute refPopCorr
    merged_mt = compute_rpc(merged_mt, mt_path)

    return merged_mt

def redo_inference(test_data_root, mouse, probe, loc):
    inference(test_data_root, mouse, probe, loc, "incl_AV008", overwrite=True)
    mt = pd.read_csv(os.path.join(test_data_root, mouse, probe, loc, "new_matchtable.csv"))
    mt = mt.drop(labels=['DNNProb', 'DNNSim'], axis=1) if 'DNNProb' in mt.columns else mt
    if "DNNProbincl_AV008" in mt.columns:
        mt = mt.rename(columns={"DNNProbincl_AV008":"DNNProb", "DNNSimincl_AV008":"DNNSim"})
    mt.to_csv(os.path.join(test_data_root, mouse, probe, loc, "new_matchtable.csv"), index=False)
    mt['uid1'] = mt['RecSes1'] * 1e6 + mt['ID1']
    mt['uid2'] = mt['RecSes2'] * 1e6 + mt['ID2']
    return mt

def undo_merge(test_data_root):
    """
    Undoes all changes to waveform file names made by the merge_and_remove_splits function.
    """

    for mouse in os.listdir(test_data_root):
        mouse_path = os.path.join(test_data_root, mouse)
        probes = os.listdir(mouse_path)
        for probe in probes:
            probe_path = os.path.join(mouse_path, probe)
            locs = os.listdir(probe_path)
            for loc in locs:
                loc_path = os.path.join(probe_path, loc)
                exps = os.listdir(loc_path)
                for exp in exps:
                    exp_path = os.path.join(loc_path, exp)
                    if os.path.isdir(exp_path):
                        waveforms_path = os.path.join(exp_path, "processed_waveforms")
                        for file in os.listdir(waveforms_path):
                            if "#" in file:
                                new_filename = file.replace("#", "")
                                os.rename(os.path.join(waveforms_path, file), 
                                          os.path.join(waveforms_path, new_filename))
                            if "+" in file:
                                os.remove(os.path.join(waveforms_path, file))


if __name__ == "__main__":

    # C_values("AL036", "19011116882", "3")
    # C_values("AL031", "19011116684", "1")
    # C_values("AL032", "19011111882", "2")
    # C_values("CB017", '19011110803', '2')
    # C_values("AV015", 'Probe0', 'IMRO_3')
    test_data_root = r"C:\Users\suyas\ALL_DATA"
    undo_merge(test_data_root)
    merge_and_remove_splits(test_data_root)