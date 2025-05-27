import numpy as np
import os, sys
sys.path.insert(0, os.getcwd())
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import zscore
from scipy.stats import fisher_exact
from scipy.ndimage import gaussian_filter1d
import mat73, h5py
import time
from tqdm import tqdm
from utils.myutil import mtpath_to_expids


def pairwise_histogram_correlation(A, B):
    # Initialize the output correlation matrix
    num_histograms = A.shape[0]
    correlation_matrix = np.zeros((num_histograms, num_histograms))

    # Compute pairwise correlations
    for i in range(num_histograms):
        for j in range(num_histograms):
            # Correlate the i-th histogram in A with the j-th histogram in B
            correlation_matrix[i, j] = np.corrcoef(A[i], B[j])[0, 1]
    
    return correlation_matrix

def times_to_spktrain(times):
    duration = times.max() - times.min()
    binsize = 0.001
    num_bins = int(duration / binsize)
    time_bins = np.linspace(0, duration, num_bins + 1)
    binned_spike_train, _ = np.histogram(times, bins=time_bins)
    spk_train = gaussian_filter1d(binned_spike_train, sigma=2)
    return spk_train

def get_spkinfo(server_root, data_root, mouse, probe, loc):
    mt_path = os.path.join(data_root, mouse, probe, loc, "new_matchtable.csv")
    mt = pd.read_csv(mt_path)
    um_path = os.path.join(server_root, mouse, probe, loc, "UnitMatch", "UnitMatch.mat")
    um = mat73.loadmat(um_path, verbose=False)
    recsesAll = um["UniqueIDConversion"]["recsesAll"].astype(int)
    goodid = um["UniqueIDConversion"]['GoodID'].astype(bool)
    recses = recsesAll[goodid]
    expids, metadata = mtpath_to_expids(mt_path, mt)
    OriID = um["UniqueIDConversion"]["OriginalClusID"].astype(int)
    # d = df.loc[(df["mouse"]==mouse)&(df["probe"]==probe)&(df["loc"]==loc)]
    # recs = d['recordings'].item()
    # exp_dict = {get_exp_id(recs[i], mouse):recs[i] for i in range(len(recs))}
    nclus = len(recses)
    if nclus**2 != len(mt):
        print("Warning - number of good units is inconsistent with length of match table.")
    return nclus, recses, expids, OriID, goodid, mt, mt_path

def get_ISI_histograms(data_root, server_root, mouse, probe, loc, ISIbins):

    nclus, recses, expids, OriID, goodid, mt, mt_path = get_spkinfo(server_root, data_root, mouse, probe, loc)
    ISIMat = np.zeros((len(ISIbins) - 1, 2, nclus))
    for clusid in tqdm(range(nclus)):
        session = recses[clusid]
        exp = expids[session]
        spikes_path = os.path.join(os.path.dirname(mt_path), exp, "spikes.npy")
        with h5py.File(spikes_path, 'r') as f:
            clusters = f['spkclus'][()]
            times = f['spktimes'][()]
        for cv in range(2):
            idx1 = np.where(clusters == OriID[goodid][clusid])[0]
            if idx1.size > 0:
                if idx1.size < 50 and cv == 0:
                    print(f"Warning: Fewer than 50 spikes for neuron {clusid}, please check your inclusion criteria")
                # Split idx1 into two halves
                if cv == 0:
                    idx1 = idx1[:len(idx1) // 2]
                else:
                    idx1 = idx1[len(idx1) // 2:]
                # nspkspersec, _ = np.histogram(times[idx1], bins=np.arange(min(times[idx1]), max(times[idx1]) + 1))
                ISIMat[:, cv, clusid], _ = np.histogram(np.diff(times[idx1].astype(float)), bins=ISIbins)
    return ISIMat, mt

def compute_values(mouse, probe, loc, new_col):
    test_data_root = os.path.join(os.path.dirname(os.getcwd()), "ALL_DATA")
    server_root = r"\\znas\Lab\Share\UNITMATCHTABLES_ENNY_CELIAN_JULIE\FullAnimal_KSChanMap"
    mice = os.listdir(test_data_root)
    # df = pd.DataFrame(read_datapaths(mice))

    # Define ISI bins
    if new_col=="newISI":
        ISIbins = np.concatenate(([0], 5 * 10 ** np.arange(-4, 0.1, 0.1)))
        for mouse in tqdm(mice):
            name_path = os.path.join(test_data_root, mouse)
            probes = os.listdir(name_path)
            for probe in probes:
                name_probe = os.path.join(name_path, probe)
                locations = os.listdir(name_probe)
                for loc in locations:
                    name_probe_location = os.path.join(name_probe, loc)
                    if not os.path.isdir(name_probe_location):
                        continue
        ISIMat, mt = get_ISI_histograms(test_data_root, server_root, mouse, probe, loc, ISIbins)
        A, B = ISIMat[:,0,:].T, ISIMat[:,1,:].T
        # now compute the value of choice
        num_histograms = A.shape[0]
        value_matrix = np.zeros((num_histograms, num_histograms))
        for i in range(num_histograms):
            for j in range(num_histograms):
                value_matrix[i,j] = np.sum(A[i] + B[j]) / np.sum((A[i] - B[j])**2)      # inverse chi-squared distance
        mt.insert(len(mt.columns), new_col, value_matrix.ravel())

    return mt
                
def write_newISI_to_all_matchtables():
    """
    Shouldn't need to call this function again - just here for posterity.
    """
    test_data_root = os.path.join(os.path.dirname(os.getcwd()), "ALL_DATA")
    server_root = r"\\znas\Lab\Share\UNITMATCHTABLES_ENNY_CELIAN_JULIE\FullAnimal_KSChanMap"
    mice = os.listdir(test_data_root)
    # df = pd.DataFrame(read_datapaths(mice))

    # Define ISI bins
    ISIbins = np.concatenate(([0], 5 * 10 ** np.arange(-4, 0.1, 0.1)))
    for mouse in tqdm(mice):
        name_path = os.path.join(test_data_root, mouse)
        probes = os.listdir(name_path)
        for probe in probes:
            name_probe = os.path.join(name_path, probe)
            locations = os.listdir(name_probe)
            for loc in locations:
                name_probe_location = os.path.join(name_probe, loc)
                if not os.path.isdir(name_probe_location):
                    continue
                mt_path = os.path.join(name_probe_location, "new_matchtable.csv")
                if not os.path.exists(mt_path):
                    continue
                mt = pd.read_csv(mt_path)
                if 'newISI' in mt.columns:
                    continue
                um_path = os.path.join(server_root, mouse, probe, loc, "UnitMatch", "UnitMatch.mat")
                um = mat73.loadmat(um_path, verbose=False)
                recsesAll = um["UniqueIDConversion"]["recsesAll"].astype(int)
                goodid = um["UniqueIDConversion"]['GoodID'].astype(bool)
                recses = recsesAll[goodid]
                expids, metadata = mtpath_to_expids(mt_path, mt)
                OriID = um["UniqueIDConversion"]["OriginalClusID"].astype(int)
                # d = df.loc[(df["mouse"]==mouse)&(df["probe"]==probe)&(df["loc"]==loc)]
                # recs = d['recordings'].item()
                # exp_dict = {get_exp_id(recs[i], mouse):recs[i] for i in range(len(recs))}
                # recompute = True
                nclus = len(recses)
                if nclus**2 != len(mt):
                    print("Warning - number of good units is inconsistent with length of match table.")
                ISIMat = np.zeros((len(ISIbins) - 1, 2, nclus))
                for clusid in tqdm(range(nclus)):
                    session = recses[clusid]
                    exp = expids[session]
                    spikes_path = os.path.join(os.path.dirname(mt_path), exp, "spikes.npy")
                    with h5py.File(spikes_path, 'r') as f:
                        clusters = f['spkclus'][()]
                        times = f['spktimes'][()]
                    for cv in range(2):
                        idx1 = np.where(clusters == OriID[goodid][clusid])[0]
                        if idx1.size > 0:
                            if idx1.size < 50 and cv == 0:
                                print(f"Warning: Fewer than 50 spikes for neuron {clusid}, please check your inclusion criteria")
                            # Split idx1 into two halves
                            if cv == 0:
                                idx1 = idx1[:len(idx1) // 2]
                            else:
                                idx1 = idx1[len(idx1) // 2:]
                            nspkspersec, _ = np.histogram(times[idx1], bins=np.arange(min(times[idx1]), max(times[idx1]) + 1))
                            ISIMat[:, cv, clusid], _ = np.histogram(np.diff(times[idx1].astype(float)), bins=ISIbins)
                correlation_matrix = pairwise_histogram_correlation(ISIMat[:, 0, :].T, ISIMat[:, 1, :].T)
                correlation_matrix = np.tanh(0.5*np.arctanh(correlation_matrix) + 0.5*np.arctanh(correlation_matrix.T))

                # Saving results in the DataFrame
                mt.insert(len(mt.columns), 'newISI', correlation_matrix.ravel())
                mt.to_csv(mt_path, index=False)              # write to the matchtable we opened already

                # write the same values to all other match tables.
                l = os.listdir(name_probe_location)
                for item in l:
                    if item[:14]!="new_matchtable" or item=="new_matchtable.csv":
                        continue
                    mt_path = os.path.join(name_probe_location, item)
                    mt = pd.read_csv(mt_path)
                    mt.insert(len(mt.columns), 'newISI', correlation_matrix.ravel())
                    mt.to_csv(mt_path, index=False)

def post_merge_ISI(mt, mt_path):
    server_root = r"\\znas\Lab\Share\UNITMATCHTABLES_ENNY_CELIAN_JULIE\FullAnimal_KSChanMap"
    ISIbins = np.concatenate(([0], 5 * 10 ** np.arange(-4, 0.1, 0.1)))
    
    expids, metadata = mtpath_to_expids(mt_path, mt)
    mouse, probe, loc = metadata['mouse'], metadata['probe'], metadata['loc']
    um_path = os.path.join(server_root, mouse, probe, loc, "UnitMatch", "UnitMatch.mat")
    um = mat73.loadmat(um_path, verbose=False)
    recsesAll = um["UniqueIDConversion"]["recsesAll"].astype(int)
    goodid = um["UniqueIDConversion"]['GoodID'].astype(bool)
    recses = recsesAll[goodid]
    OriID = um["UniqueIDConversion"]["OriginalClusID"].astype(int)
    nclus = len(recses)
    nclus_after_removal = int(np.sqrt(len(mt)))
    ISIMat = np.zeros((len(ISIbins) - 1, 2, nclus_after_removal))
    skips = 0
    for clusid in tqdm(range(nclus)):
        session = recses[clusid]
        exp = expids[session]
        spikes_path = os.path.join(os.path.dirname(mt_path), exp, "merged_spikes.npy")
        waveforms_path = os.path.join(os.path.dirname(mt_path), exp, "processed_waveforms")
        waves = os.listdir(waveforms_path)
        id = OriID[goodid][clusid]
        with h5py.File(spikes_path, 'r') as f:
            clusters = f['spkclus'][()]
            times = f['spktimes'][()]
        idx1 = np.where(clusters == id)[0]
        if len(idx1) == 0 or f"Unit{id}#_RawSpikes.npy" in waves:
            skips += 1
            continue
        for cv in range(2):
            if idx1.size > 0:
                if idx1.size < 50 and cv == 0:
                    print(f"Warning: Fewer than 50 spikes for neuron {clusid}, please check your inclusion criteria")
                # Split idx1 into two halves
                if cv == 0:
                    idx1 = idx1[:len(idx1) // 2]
                else:
                    idx1 = idx1[len(idx1) // 2:]
                ISIMat[:, cv, clusid-skips], _ = np.histogram(np.diff(times[idx1].astype(float)), bins=ISIbins)

    correlation_matrix = pairwise_histogram_correlation(ISIMat[:, 0, :].T, ISIMat[:, 1, :].T)
    correlation_matrix = np.tanh(0.5*np.arctanh(correlation_matrix) + 0.5*np.arctanh(correlation_matrix.T))

    # Saving results in the DataFrame
    if 'newISI' in mt.columns:
        mt['newISI'] = correlation_matrix.ravel()
    else:
        mt.insert(len(mt.columns), 'newISI', correlation_matrix.ravel())
    return mt

def cross_correlation_fingerprint(session_correlations_all, pairs, unit_to_take, recses_good, plt=True):
    """
    Computes cross-correlation fingerprints across sessions.

    Parameters:
    - session_correlations_all: List of dictionaries with 'fold1' and 'fold2' matrices for each session.
    - pairs: Array of neuron index pairs (Nx2) that likely match across sessions.
    - unit_to_take: Array of unit IDs.
    - recses_good: Array of session IDs for each unit.
    - plt: Boolean to control plotting (default is True).

    Returns:
    - FingerprintRAll: Cross-correlation fingerprint matrix.
    - SigMask: Significance mask (binary matrix for illustration, if needed).
    - AllSessionCorrelationsFingerprints: Nested list with fold correlations for each session pair.
    """

    # Parameters
    nclus = len(unit_to_take)
    ndays = len(session_correlations_all)
    rec_opt = np.unique(recses_good)

    # Session switch to mimic MATLAB cumulative indexing
    session_switch = [0] + np.cumsum([session['fold1'].shape[0] for session in session_correlations_all]).tolist()
    all_session_correlations_fingerprints = [[None for _ in range(ndays)] for _ in range(ndays)]
    fingerprint_r_all = np.full((nclus, nclus), np.nan)
    sig_mask = np.zeros((nclus, nclus), dtype=int)

    # Step 1: Average fold correlations per session
    session_correlations_per_day = []
    for session in session_correlations_all:
        avg_corr = np.tanh(np.nanmean(np.arctanh(np.stack([session['fold1'], session['fold2']])), axis=0))
        session_correlations_per_day.append(avg_corr)

    # Step 2: Compute fingerprints within and across days
    for did1 in range(ndays):
        for did2 in range(ndays):
            clus_idx_d1_all = range(session_switch[did1], session_switch[did1 + 1])
            clus_idx_d2_all = range(session_switch[did2], session_switch[did2 + 1])

            if did1 == did2:  # Within-session comparison
                fold1 = session_correlations_all[did1]['fold1']
                fold2 = session_correlations_all[did1]['fold2']

                fingerprint_r = np.corrcoef(fold1, fold2, rowvar=False)
                np.fill_diagonal(fingerprint_r, np.nan)

                all_session_correlations_fingerprints[did1][did2] = np.hstack([fold1, fold2])

            else:  # Across-session comparison
                dayopt = [did1, did2]
                session_corr_pair = []

                for idx, day in enumerate(dayopt):
                    pair_idx = np.logical_and(
                        recses_good[pairs[:, 0]] == rec_opt[day],
                        recses_good[pairs[:, 1]] == rec_opt[dayopt[1 - idx]]
                    )
                    pairs_tmp = pairs[pair_idx]

                    # Keep unique matches only
                    pairs_tmp = np.unique(pairs_tmp[:, idx], return_index=True)[1]
                    units_to_take_idx = pairs_tmp

                    # Extract correlation matrix for units
                    unit_indices = np.where(np.isin(unit_to_take, unit_to_take[units_to_take_idx]))[0]
                    session_corr_matrix = session_correlations_per_day[day][unit_indices]

                    # Remove diagonal entries (self-correlations)
                    for i, unit_idx in enumerate(units_to_take_idx):
                        session_corr_matrix[i, unit_indices == unit_idx] = np.nan

                    session_corr_pair.append(session_corr_matrix)

                # Correlation between two sessions
                if session_corr_pair[0].size == 0 or session_corr_pair[1].size == 0:
                    fingerprint_r = np.zeros((len(clus_idx_d1_all), len(clus_idx_d2_all)))
                else:
                    fingerprint_r = np.corrcoef(session_corr_pair[0], session_corr_pair[1], rowvar=False)

                all_session_correlations_fingerprints[did1][did2] = np.hstack(session_corr_pair)

            # Save fingerprint correlations
            fingerprint_r_all[np.ix_(clus_idx_d1_all, clus_idx_d2_all)] = fingerprint_r

    return fingerprint_r_all, sig_mask, all_session_correlations_fingerprints

def compute_rpc(mt:pd.DataFrame, mt_path):
    """
    Computes cross-correlation fingerprints for multiple sessions.
    Input is the full merged match table. Output has an extra column 'newRPC' with recomputed RPC values.
    """
    server_root = r"\\znas\Lab\Share\UNITMATCHTABLES_ENNY_CELIAN_JULIE\FullAnimal_KSChanMap"
    expids, metadata = mtpath_to_expids(mt_path, mt)
    mouse, probe, loc = metadata['mouse'], metadata['probe'], metadata['loc']
    um_path = os.path.join(server_root, mouse, probe, loc, "UnitMatch", "UnitMatch.mat")
    um = mat73.loadmat(um_path, verbose=False)
    session_correlations = []
    bin_size = um['UMparam']['binsz']

    for ses in mt["RecSes1"].unique():

        # Load the spike data
        spk_path = os.path.join(os.path.dirname(mt_path), expids[ses], 'spikes.npy')
        with h5py.File(spk_path, 'r') as f:
            clusters = f['spkclus'][()]
            times = f['spktimes'][()]
        good_neuron_ids = np.unique(mt.loc[mt['RecSes1']==ses, "ID1"].values)
        if len(good_neuron_ids) < 2:
            continue

        # Define edges for binning
        edges = np.arange(np.floor(np.min(times)) - bin_size/2, np.ceil(np.max(times)) + bin_size/2, bin_size)

        # Create histogram matrix for good neurons
        num_good_neurons = len(good_neuron_ids)
        num_bins = len(edges) - 1
        hist_matrix = np.zeros((num_good_neurons, num_bins))
        for idx, neuron_id in enumerate(good_neuron_ids):
            neuron_spike_times = times[np.where(clusters == neuron_id)]
            hist_matrix[idx, :] = np.histogram(neuron_spike_times, bins=edges)[0]

        # Compute cross-correlation matrices for each half
        fold1_corr = np.corrcoef(hist_matrix[:, :num_bins // 2])
        fold2_corr = np.corrcoef(hist_matrix[:, num_bins // 2:])
        np.fill_diagonal(fold1_corr, np.nan)
        np.fill_diagonal(fold2_corr, np.nan)

        session_correlations.append({'fold1': fold1_corr, 'fold2': fold2_corr, 'session':ses, 'ids':good_neuron_ids})

    mt.insert(len(mt.columns), "newRPC", '')
    for corrdict1 in session_correlations:
        for corrdict2 in session_correlations:
            session1, session2 = corrdict1['session'], corrdict2['session']
            df = mt.loc[(mt['RecSes1']==session1) & (mt['RecSes2']==session2),:]
            assert len(df) == corrdict1['fold1'].shape[0] * corrdict2['fold2'].shape[0]
            assert len(df) == len(corrdict1['ids']) * len(corrdict2['ids'])
            refpop,_ = remove_conflicts(df.loc[df['MatchProb']>0.5,:], 'MatchProb')
            N1, N2 = corrdict1['fold1'].shape[0], corrdict2['fold2'].shape[0]
            mat1_idx, mat2_idx = [], []
            for idx, row in refpop.iterrows():
                id1=row['ID1']
                id2=row['ID2']
                mat1_idx.append(np.where(corrdict1['ids']==id1)[0][0])
                mat2_idx.append(np.where(corrdict2['ids']==id2)[0][0])
            reduced1 = corrdict1['fold1'][:,mat1_idx]
            reduced2 = corrdict2['fold2'][:,mat2_idx]
            if len(refpop)==0:
                # no matches found -> set RPC to 0 (consistent with UnitMatch implementation)
                rpc = np.zeros((N1, N2))
            else:
                rpc = np.empty((N1, N2))
                for i,row1 in enumerate(reduced1):
                    for j,row2 in enumerate(reduced2):
                        nan_indices = np.argwhere(np.isnan(row1) | np.isnan(row2))
                        r1 = [item for index,item in enumerate(row1) if index not in nan_indices]
                        r2 = [item for index,item in enumerate(row2) if index not in nan_indices]
                        rpc[i,j] = np.corrcoef(r1,r2)[1,0]
            mt.loc[(mt['RecSes1']==session1) & (mt['RecSes2']==session2),"newRPC"] = rpc.ravel()

    return mt

if __name__=="__main__":
    # write_newISI_to_all_matchtables()            # recomputes the ISICorr values and stores in newISI column.
    data_root = os.path.join(os.path.dirname(os.getcwd()), "ALL_DATA")
    mouse, probe, loc = "AL031", "19011116684", "1"
    # mt = compute_values(mouse, probe, loc, 'wavelet')
    mt_path = os.path.join(data_root, mouse, probe, loc, "new_matchtable.csv")
    # print(test_metric(mt, 'wavelet', vis=True))
    # save_diagrams("AL036_chi", mouse, loc, 'chi', mt, mt_path)

    mt = pd.read_csv(mt_path)
    rpc = compute_rpc(mt, mt_path)