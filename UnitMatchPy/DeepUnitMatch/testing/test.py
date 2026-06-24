import os, sys
from pathlib import Path

# Ensure the DeepUnitMatch package root is on the path so `utils` is
# importable regardless of the caller's working directory.
sys.path.insert(0, str(Path(__file__).parent.parent))

from utils.losses import clip_sim, CustomClipLoss, Projector
from utils.npdataset import NeuropixelsDataset_cortexlab, ValidationExperimentBatchSampler
from utils.helpers import read_pos
import numpy as np
import matplotlib.pyplot as plt
from utils.mymodel import SpatioTemporalCNN_V2
import torch
from torch.utils.data import DataLoader
from tqdm import tqdm
import pandas as pd
from sklearn.neighbors import KernelDensity
from scipy.optimize import linear_sum_assignment
import importlib
from utils import helpers
importlib.reload(helpers)
from utils.helpers import *


def _parse_unitmatch_good_units(unit_label_paths):
    """
    Mirror UnitMatchPy.utils.load_good_waveforms() selection and ordering from TSVs.
    Returns list[list[int]] good unit IDs per session, preserving TSV row order.
    """
    if unit_label_paths is None:
        return None

    good_units = []
    is_bombcell = os.path.basename(unit_label_paths[0]) == "cluster_bc_unitType.tsv"

    for path in unit_label_paths:
        session_good = []
        with open(path, "r", encoding="utf-8", errors="ignore") as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                parts = line.split("\t")
                if len(parts) < 2:
                    continue
                try:
                    unit_id = int(parts[0])
                except ValueError:
                    # header row
                    continue
                label = str(parts[1]).strip().lower()
                if is_bombcell:
                    if label in {"good", "non-soma good"}:
                        session_good.append(unit_id)
                else:
                    if label == "good":
                        session_good.append(unit_id)
        good_units.append(session_good)

    return good_units



def load_trained_model(device="cpu", read_path=None):

    model = SpatioTemporalCNN_V2(n_channel=30,n_time=60,n_output=256).to(device)
    model = model.double()
    if read_path is None:
        current_dir = Path(__file__).parent.parent
        read_path = current_dir / "utils" / "model"
    checkpoint = torch.load(read_path)
    model.load_state_dict(checkpoint['model'])
    clip_loss = CustomClipLoss().to(device)
    clip_loss.load_state_dict(checkpoint['clip_loss'])
    model.eval()
    clip_loss.eval()

    # Load projector
    projector = Projector(input_dim=256, output_dim=128, hidden_dim=128, n_hidden_layers=1, dropout=0.1).to(device)
    projector = projector.double()

    # Can also return projector if needed
    return model

def reorder_by_depth(matrix:np.ndarray, pos1:np.ndarray, pos2) -> np.ndarray:
    """
    Matrix should compare just one recording session against another.
    """

    depths1 = pos1[:,1]
    depths2 = pos2[:,1]

    if matrix.shape[0] == depths1.shape[0]:
        sort_indices = np.argsort(depths1)
        sorted_matrix = matrix[sort_indices]
    else:
        sort_indices = np.argsort(depths2)
        sorted_matrix = matrix[sort_indices]
    if matrix.shape[1] == depths1.shape[0]:
        sort_indices = np.argsort(depths1)
        sorted_matrix = sorted_matrix[:, sort_indices]
    else:
        sort_indices = np.argsort(depths2)
        sorted_matrix = sorted_matrix[:, sort_indices]
    
    return sorted_matrix

def inference(model, data_dir, unit_label_paths=None):

    # If you pass `unit_label_paths` (the same TSV paths used by UnitMatchPy.utils.load_good_waveforms),
    # this will build the dataset in the exact same per-session unit order.
    good_units = _parse_unitmatch_good_units(unit_label_paths)
    unit_order = good_units if good_units is not None else "unitmatch"
    test_dataset = NeuropixelsDataset_cortexlab(data_dir, unit_order=unit_order)
    test_sampler = ValidationExperimentBatchSampler(test_dataset, shuffle = False)
    test_loader = DataLoader(test_dataset, batch_sampler=test_sampler)

    submatrices = []
    n_batches = len(test_loader)

    for estimates_i, _, positions_i, exp_ids_i, filepaths_i in tqdm(test_loader):
        # Forward pass
        enc_estimates_i = model(estimates_i)        # shape [bsz, 256]

        for _, candidates_j, positions_j, exp_ids_j, filepaths_j in tqdm(test_loader):
            enc_candidates_j = model(candidates_j)
            s = clip_sim(enc_estimates_i, enc_candidates_j)
            submatrices.append(s.detach().cpu().numpy())
    
    result_rows = []
    for i in range(n_batches):
        row_matrices = []
        for j in range(n_batches):
            matrix_idx = i * n_batches + j
            row_matrices.append(submatrices[matrix_idx])
        row = np.hstack(row_matrices)
        result_rows.append(row)
    
    result = np.vstack(result_rows)

    return result

def get_threshold(prob_matrix:np.ndarray, session_id, MAP=False):

    n = len(session_id)
    within_session = (session_id[:, None] == session_id).astype(int)
    pmat = prob_matrix[within_session==True].reshape(-1, 1)
    labels = np.eye(n)[within_session==True].reshape(-1, 1)

    # On-diagonal means same neuron. Off-diagonal means different neurons.
    on_diag = pmat[labels==1]
    off_diag = pmat[labels==0]

    # Kernel density estimation (distributions are more useful than histograms)
    kde_on = KernelDensity(kernel='gaussian', bandwidth=0.01).fit(on_diag.reshape(-1, 1))
    kde_off = KernelDensity(kernel='gaussian', bandwidth=0.01).fit(off_diag.reshape(-1, 1))
    x = np.linspace(min(off_diag), max(on_diag), 1000).reshape(-1, 1)
    y_on = np.exp(kde_on.score_samples(x))
    y_off = np.exp(kde_off.score_samples(x))

    # Find the threshold where the distributions intersect
    if MAP:
        thresh = np.argwhere(np.diff(np.sign(y_off*(n-1) - y_on)))
    else:
        thresh=np.argwhere(np.diff(np.sign(y_off - y_on)))
    if len(thresh) == 0:
        thresh = len(x) - 1
    elif len(thresh) > 1:
        thresh = thresh[-1]

    across = prob_matrix[within_session==False].reshape(-1, 1)
    diff = np.median(pmat) - np.median(across)

    return x[thresh].item() - diff

def directional_filter_df(matches: pd.DataFrame):
    filtered_matches = matches.copy()
    for idx, row in matches.iterrows():
        i1 = row["ID1"]
        i2 = row["ID2"]
        r1 = row["RecSes1"]
        r2 = row["RecSes2"]

        # Check for the reverse match
        reverse_match = matches.loc[
            (matches["RecSes1"] == r2) & (matches["ID1"] == i2) & 
            (matches["RecSes2"] == r1) & (matches["ID2"] == i1)
        ]

        if reverse_match.empty:
            filtered_matches = filtered_matches.drop(idx)  # Drop if no reverse match is found
    return filtered_matches

def directional_filter(sim_matrix, session_id, threshold):
    """
    Matrix version of directional_filter_df.
    Returns a boolean matrix where (i,j) is True only if both (i,j) and (j,i)
    exceed the threshold and belong to different sessions.
    """
    above_thresh = sim_matrix > threshold
    across_session = session_id[:, None] != session_id[None, :]
    candidates = above_thresh & across_session
    return candidates & candidates.T

def remove_conflicts(matches: pd.DataFrame, metric: str):
    
    # Find the best match for each neuron1 (RecSes1, ID1 combination)
    # For each group, get the index of the row with maximum metric value
    best_neuron1_indices = matches.groupby(['RecSes1', 'ID1'])[metric].idxmax()
    
    # Find the best match for each neuron2 (RecSes2, ID2 combination)  
    best_neuron2_indices = matches.groupby(['RecSes2', 'ID2'])[metric].idxmax()
    
    # A match is valid only if it's the best match for BOTH neurons involved
    # Take intersection of the two sets of indices
    valid_indices = set(best_neuron1_indices) & set(best_neuron2_indices)
    
    # Count conflicts (total matches minus valid matches)
    num_conflicts = len(matches) - len(valid_indices)
    
    # Filter to keep only valid matches
    filtered_matches = matches.loc[list(valid_indices)]

    return filtered_matches, num_conflicts

def remove_conflicts_matrix(sim_matrix, candidates):
    """
    Matrix version of remove_conflicts using the Hungarian algorithm.
    Finds the globally optimal 1-to-1 assignment that maximises total similarity.

    Parameters
    ----------
    sim_matrix : np.ndarray (n, n) — similarity scores
    candidates : np.ndarray (n, n) bool — mask of eligible matches

    Returns
    -------
    matches : np.ndarray (n, n) bool
    num_conflicts : int
    """
    active_rows = np.where(candidates.any(axis=1))[0]
    active_cols = np.where(candidates.any(axis=0))[0]

    if len(active_rows) == 0 or len(active_cols) == 0:
        return np.zeros_like(candidates, dtype=bool), int(candidates.sum())

    # Extract submatrix of active neurons only
    sub_sim = sim_matrix[np.ix_(active_rows, active_cols)]
    sub_cand = candidates[np.ix_(active_rows, active_cols)]

    # Set non-candidate entries to large negative so they're never assigned
    sub_cost = np.where(sub_cand, sub_sim, -1e9)

    # Hungarian algorithm — globally optimal 1-to-1 assignment
    ri, ci = linear_sum_assignment(sub_cost, maximize=True)

    # Keep only assignments that were actual candidates
    valid = sub_cand[ri, ci]
    ri, ci = ri[valid], ci[valid]

    # Map back to original indices
    matches = np.zeros_like(candidates, dtype=bool)
    matches[active_rows[ri], active_cols[ci]] = True

    num_conflicts = int(candidates.sum() - matches.sum())
    return matches, num_conflicts

def get_matches(df, sim_matrix, session_id, data_dir, dist_thresh):
    """
    Process the output probability matrix to get final set of matches across sessions.
    Output is the probability matrix of matches after thresholding and filtering and a boolean matrix indicating final matches.
    """

    sim_thresh = get_threshold(sim_matrix, session_id)

    matches = df.loc[df["Prob"] > sim_thresh].copy()

    matches = matches.loc[matches["RecSes1"] != matches["RecSes2"]]   # only keep across-session matches
    matches = directional_filter_df(matches)

    # spatial filtering
    sessions = np.unique(session_id)
    positions = {}
    for session in sessions:
        session_path = os.path.join(data_dir, str(session))
        positions[session] = read_pos(session_path)
    corrections = get_corrections(matches, positions)
    matches_with_dist = vectorised_drift_corrected_dist(corrections, positions, matches)
    matches = matches_with_dist.loc[matches_with_dist["dist"]<dist_thresh]

    matches = matches[[c for c in df.columns]]
    matches, confs = remove_conflicts(matches, "Prob")

    matches.insert(len(matches.columns), "match", True)
    df.insert(len(df.columns), "match", False)
    df.loc[matches.index, "match"] = True
    df.loc[(df['RecSes1'] == df['RecSes2']) & (df['ID1'] == df['ID2']), "match"] = True
    return df

def vectorised_drift_corrected_dist(corrections, positions, matches:pd.DataFrame, nocorr=False):
    """
    Vectorised drift-corrected distance computation.
    Returns matches dataframe with extra column for drift-corrected distance.
    """

    sessions = list(positions.keys())

    # Extract the x, y positions
    idx = matches.index
    coords1 = pd.concat([positions[session].loc[:, ['ID', 'x', 'y']]
                        .rename(columns={'ID': 'ID1', 'x': 'x1', 'y': 'y1'})
                        .assign(RecSes1=session)
                        for session in [sessions[0], sessions[1]]])

    coords2 = pd.concat([positions[session].loc[:, ['ID', 'x', 'y']]
                        .rename(columns={'ID': 'ID2', 'x': 'x2', 'y': 'y2'})
                        .assign(RecSes2=session)
                        for session in [sessions[1], sessions[0]]])

    # Merge with matches DataFrame
    newmatches = matches.merge(coords1, on=['RecSes1', 'ID1'])
    newmatches = newmatches.merge(coords2, on=['RecSes2', 'ID2'])
    newmatches = newmatches.drop_duplicates(subset=['ID1', 'ID2', 'RecSes1', 'RecSes2'])
    newmatches.index = idx
    matches = newmatches

    # Calculate shanks (rounded x coordinates)
    matches['shank1'] = matches['x1'].round(-2).astype(int)
    matches['shank2'] = matches['x2'].round(-2).astype(int)
    # Apply corrections where necessary
    if not nocorr:
        matches["idx"] = matches.index
        matches = matches.merge(
            corrections[['rec1', 'rec2', 'shank', 'ydiff']], 
            how='left', 
            left_on=['RecSes1', 'RecSes2', 'shank1'], 
            right_on=['rec1', 'rec2', 'shank']
        )
        matches.index = matches["idx"]
        matches.drop(['rec1', 'rec2', 'shank', 'idx'], axis=1, inplace=True)
        # Apply the drift correction conditionally
        matches['y2_corr'] = np.where(
            matches['shank1'] != matches['shank2'], 
            1000, 
            matches['y2'] - matches['ydiff'].fillna(0)
        )
    else:
        matches['y2_corr'] = matches['y2']
    # Calculate the Euclidean distance
    matches['dist'] = np.sqrt((matches['x1'] - matches['x2'])**2 + (matches['y1'] - matches['y2_corr'])**2)
    return matches

def get_corrections(matches, positions):

    pos1 = pd.concat([positions[key].assign(RecSes1=key) for key in positions.keys()], ignore_index=True)
    pos2 = pos1.rename(columns={"RecSes1": "RecSes2", "x": "x2", "y": "y2", "ID": "ID2"})

    # Merge matches with positions data to get all necessary columns in one DataFrame
    matches = (matches
               .merge(pos1, left_on=["RecSes1", "ID1"], right_on=["RecSes1", "ID"], how="left")
               .merge(pos2, left_on=["RecSes2", "ID2"], right_on=["RecSes2", "ID2"], how="left"))
    
    # Drop rows with missing positions data
    matches = matches.dropna(subset=["x", "y", "x2", "y2"])

    # Calculate shank and ydiff in a vectorized way
    matches['shank1'] = matches['x'].round(-2)
    matches['shank2'] = matches['x2'].round(-2)
    matches['ydiff'] = matches['y2'] - matches['y']

    # Filter only rows with matching shanks
    matches = matches[matches['shank1'] == matches['shank2']]

    # Group by session pairs and shank, then calculate median ydiff
    output = (matches.groupby(['RecSes1', 'RecSes2', 'shank1'])['ydiff']
              .median()
              .reset_index()
              .rename(columns={'shank1': 'shank','RecSes1':'rec1', 'RecSes2':'rec2'}))
    return output

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

def get_ISI_histograms(param, ISIbins):

    fs = 3e4
    nclus = param['n_units']
    KSdirs = param['KS_dirs']
    good_units = param['good_units']

    ISIMat = np.zeros((len(ISIbins) - 1, 2, nclus))
    index = 0
    for session in range(len(good_units)):

        times = np.load(os.path.join(KSdirs[session], "spike_times.npy")) / fs
        clusters = np.load(os.path.join(KSdirs[session], "spike_clusters.npy"))

        for clusid in tqdm(good_units[session].squeeze()):
            idx1 = np.where(clusters == clusid)[0]

            for cv in range(2):    
                if idx1.size > 0:
                    if idx1.size < 50 and cv == 0:
                        print(f"Warning: Fewer than 50 spikes for neuron {clusid}, please check your inclusion criteria")
                    # Split idx1 into two halves
                    if cv == 0:
                        idx1 = idx1[:len(idx1) // 2]
                    else:
                        idx1 = idx1[len(idx1) // 2:]
                    ISIMat[:, cv, index], _ = np.histogram(np.diff(times[idx1].astype(float)), bins=ISIbins)

            index += 1
    return ISIMat

def ISI_correlations(param):

    # Define ISI bins
    ISIbins = np.concatenate(([0], 5 * 10 ** np.arange(-4, 0.1, 0.1)))
    ISIMat = get_ISI_histograms(param, ISIbins)
    A, B = ISIMat[:,0,:].T, ISIMat[:,1,:].T                                                         # pull out each cross-validation fold

    return pairwise_histogram_correlation(A, B)                                                      # compute pairwise correlations

def get_binned_psth(param, bin_size=0.01):
    """
    Bin spike times into a population activity matrix per session.
    Returns list of (n_units_session x n_time_bins) arrays, one per session.
    bin_size is in seconds (default 10 ms, matching MATLAB UMparam.binsz).
    """
    fs = 3e4
    KSdirs = param['KS_dirs']
    good_units = param['good_units']

    session_psthas = []
    for session in range(len(good_units)):
        times = np.load(os.path.join(KSdirs[session], "spike_times.npy")) / fs
        clusters = np.load(os.path.join(KSdirs[session], "spike_clusters.npy"))

        edges = np.arange(times.min() - bin_size / 2, times.max() + bin_size, bin_size)

        session_units = good_units[session]
        if hasattr(session_units, 'squeeze'):
            session_units = session_units.squeeze()

        sr = np.zeros((len(session_units), len(edges) - 1))
        for uid, clusid in enumerate(session_units):
            idx = np.where(clusters == clusid)[0]
            if idx.size > 0:
                sr[uid, :], _ = np.histogram(times[idx], bins=edges)

        session_psthas.append(sr)

    return session_psthas


def _normalize_fingerprints(C):
    """
    Replace diagonal NaN with row mean, then standardize each row to zero mean / unit std.
    Units with zero variance (no spikes) get an all-zero fingerprint.
    """
    n = C.shape[0]
    C_f = C.copy()
    for k in range(n):
        C_f[k, k] = np.nanmean(C[k])
    C_f = np.nan_to_num(C_f, nan=0.0)
    m = C_f.mean(axis=1, keepdims=True)
    s = C_f.std(axis=1, keepdims=True)
    s[s == 0] = 1.0
    return (C_f - m) / s


def refpop_correlations(param, bin_size=0.01):
    """
    Cross-correlation population fingerprint (refPopCorr), analogous to the
    MATLAB ComputeFunctionalScores refPopCorr.

    For each unit the fingerprint is its row in the within-session pairwise
    correlation matrix (one fold per cross-validation half).
    refPopCorr[i, j] = Pearson r between the fold-1 fingerprint of unit i and
    the fold-2 fingerprint of unit j.  Higher values indicate more similar
    population coupling (likely the same unit).

    For cross-session pairs whose sessions have different unit counts, the
    fingerprints are truncated to the shorter length (reasonable when the
    same population is recorded across days).
    """
    session_psthas = get_binned_psth(param, bin_size)
    good_units = param['good_units']

    session_C1_norm, session_C2_norm = [], []
    for sr in session_psthas:
        n_fold = sr.shape[1] // 2
        C1 = np.corrcoef(sr[:, :n_fold])
        C2 = np.corrcoef(sr[:, n_fold:2 * n_fold])
        # Zero out NaN from units with no spikes before setting diagonal
        C1 = np.nan_to_num(C1, nan=0.0)
        C2 = np.nan_to_num(C2, nan=0.0)
        np.fill_diagonal(C1, np.nan)
        np.fill_diagonal(C2, np.nan)
        session_C1_norm.append(_normalize_fingerprints(C1))
        session_C2_norm.append(_normalize_fingerprints(C2))

    # Map each global unit index to its session and within-session position
    unit_session = []
    for session, units in enumerate(good_units):
        if hasattr(units, 'squeeze'):
            units = units.squeeze()
        unit_session.extend([session] * len(units))
    unit_session = np.array(unit_session)

    nclus = param['n_units']
    refPopCorr = np.zeros((nclus, nclus))

    for si in range(len(good_units)):
        for sj in range(len(good_units)):
            rows_i = np.where(unit_session == si)[0]
            rows_j = np.where(unit_session == sj)[0]
            C1n = session_C1_norm[si]   # (n_si, n_si)
            C2n = session_C2_norm[sj]   # (n_sj, n_sj)
            min_n = min(C1n.shape[1], C2n.shape[1])
            # Vectorised pairwise correlation via normalised dot product
            block = (C1n[:, :min_n] @ C2n[:, :min_n].T) / min_n
            refPopCorr[np.ix_(rows_i, rows_j)] = np.clip(block, -1, 1)

    return refPopCorr


def get_FR(param):
    """
    Compute mean firing rate (spikes/s) per unit for two cross-validation halves.
    Returns array of shape (2, n_units): FR[fold, unit].
    """
    fs = 3e4
    nclus = param['n_units']
    KSdirs = param['KS_dirs']
    good_units = param['good_units']

    FR = np.zeros((2, nclus))
    index = 0
    for session in range(len(good_units)):
        times = np.load(os.path.join(KSdirs[session], "spike_times.npy")) / fs
        clusters = np.load(os.path.join(KSdirs[session], "spike_clusters.npy"))

        session_units = good_units[session]
        if hasattr(session_units, 'squeeze'):
            session_units = session_units.squeeze()

        for clusid in session_units:
            idx = np.where(clusters == clusid)[0]
            if idx.size > 0:
                n = len(idx)
                for cv, cv_idx in enumerate([idx[:n // 2], idx[n // 2:]]):
                    t_cv = times[cv_idx]
                    if len(t_cv) > 1:
                        bins = np.arange(int(t_cv.min()), int(t_cv.max()) + 2)
                        if len(bins) > 1:
                            FR[cv, index] = np.nanmean(np.histogram(t_cv, bins=bins)[0])
            index += 1

    return FR


def FR_diff(param):
    """
    Pairwise absolute firing rate difference (cross-validated), analogous to
    the MATLAB ComputeFunctionalScores FRDiff.

    FRDiff[i, j] = |FR_fold2[i] - FR_fold1[j]|.
    Lower values indicate more similar firing rates (likely same unit).

    Note: the AUC function expects higher = better match, so pass
    -FR_diff(param) when calling AUC for this metric.
    """
    FR = get_FR(param)   # (2, n_units)
    return np.abs(FR[1, :, np.newaxis] - FR[0, np.newaxis, :])


def get_ISI_CV(param):
    """
    Compute the coefficient of variation (CV = std/mean) of interspike intervals
    per unit for two cross-validation halves.
    Returns array of shape (2, n_units): CV[fold, unit].
    Units with fewer than 2 spikes in a fold get CV = 0.
    """
    fs = 3e4
    nclus = param['n_units']
    KSdirs = param['KS_dirs']
    good_units = param['good_units']

    CV = np.zeros((2, nclus))
    index = 0
    for session in range(len(good_units)):
        times = np.load(os.path.join(KSdirs[session], "spike_times.npy")) / fs
        clusters = np.load(os.path.join(KSdirs[session], "spike_clusters.npy"))

        session_units = good_units[session]
        if hasattr(session_units, 'squeeze'):
            session_units = session_units.squeeze()

        for clusid in session_units:
            idx = np.where(clusters == clusid)[0]
            if idx.size > 1:
                n = len(idx)
                for cv, cv_idx in enumerate([idx[:n // 2], idx[n // 2:]]):
                    isis = np.diff(times[cv_idx].astype(float))
                    if len(isis) > 0 and isis.mean() > 0:
                        CV[cv, index] = isis.std() / isis.mean()
            index += 1

    return CV


def ISI_CV_diff(param):
    """
    Pairwise absolute difference in ISI coefficient of variation (cross-validated).
    CVDiff[i, j] = |CV_fold2[i] - CV_fold1[j]|.
    Lower values indicate more similar firing regularity (likely same unit).

    Note: the AUC function expects higher = better match, so pass
    -ISI_CV_diff(param) when calling AUC for this metric.
    """
    CV = get_ISI_CV(param)   # (2, n_units)
    return np.abs(CV[1, :, np.newaxis] - CV[0, np.newaxis, :])

def get_natim_responses(param):
    """
    Bin spike times into stimulus-triggered PSTH per session.
    Returns list of arrays with shape (n_time_bins, n_stimuli, n_units, max_repeats).
    Organized as: time × stimulus × neurons × repeats
    Time dimension: onset bins followed by offset bins (concatenated).
    """
    fs = 3e4
    nclus = param['n_units']
    KSdirs = param['KS_dirs']
    good_units = param['good_units']

    session_timecourses = []
    session_stimulus_responses = []
    for session in range(len(good_units)):
        times = np.load(os.path.join(KSdirs[session], "spike_times.npy")) / fs
        clusters = np.load(os.path.join(KSdirs[session], "spike_clusters.npy"))

        session_units = good_units[session]
        if hasattr(session_units, 'squeeze'):
            session_units = session_units.squeeze()
        n_units = len(session_units)

        # Try to load trial files; if they don't exist, use NaN arrays
        trial_path = os.path.dirname(os.path.dirname(KSdirs[session]))
        try:
            trials_IDs = np.load(os.path.join(trial_path, "trial.imageIDs.npy"))
            trials_offsetTimes = np.load(os.path.join(trial_path, "trial.offsetTimes.npy"))
            trials_onsetTimes = np.load(os.path.join(trial_path, "trial.onsetTimes.npy"))
        except FileNotFoundError:
            # Trial files don't exist; return NaN arrays
            bin_size = 0.005
            onset_window = [-0.3, 0.5]
            offset_window = [0.0, 0.5]
            onset_bins = np.arange(onset_window[0], onset_window[1] + bin_size, bin_size)
            offset_bins = np.arange(offset_window[0], offset_window[1] + bin_size, bin_size)
            n_time_bins_total = (len(onset_bins) - 1) + (len(offset_bins) - 1)
            
            timecourse = np.full((2, n_time_bins_total, n_units), np.nan)
            stimulus_response = np.full((2, 112, n_units), np.nan)
            
            session_timecourses.append(timecourse)
            session_stimulus_responses.append(stimulus_response)
            continue

        bin_size = 0.005
        onset_window = [-0.3, 0.5]
        offset_window = [0.0, 0.5]
        psth, _ = get_stimulus_triggered_psth(times, clusters, trials_onsetTimes, trials_offsetTimes, 
                                 trials_IDs, session_units, 
                                 onset_window=onset_window, offset_window=offset_window, 
                                 bin_size=bin_size)
        
        # psth shape: (time, stimulus, neurons, repeats)
        # Average timecourse over stimuli and repeats: (fold, time, neurons)
        timecourse1 = np.nanmean(psth[:,:,:,1::2], axis=(1, 3))
        timecourse2 = np.nanmean(psth[:,:,:,::2], axis=(1, 3))
        timecourse = np.stack([timecourse1, timecourse2], axis=0)  # shape: (fold, time, neurons)
        
        # Average response across repeats for each stimulus, using only first 200ms after onset
        # onset_window = [-0.3, 0.5], we want time 0 to 0.2s
        ms_after_onset = 0.2
        bins_in_window = int(ms_after_onset / bin_size)  # 40 bins
        onset_start_bin = int((0 - onset_window[0]) / bin_size)  # skip to time 0 (60 bins)
        onset_end_bin = onset_start_bin + bins_in_window  # 100
        stimulus_response1 = np.nanmean(psth[onset_start_bin:onset_end_bin, :, :, 1::2], axis=(0, 3))
        stimulus_response2 = np.nanmean(psth[onset_start_bin:onset_end_bin, :, :, ::2], axis=(0, 3))
        stimulus_response = np.stack([stimulus_response1, stimulus_response2], axis=0)  # shape: (fold, stimulus, neurons)

        session_timecourses.append(timecourse)
        session_stimulus_responses.append(stimulus_response)

    return session_stimulus_responses, session_timecourses

def get_stimulus_triggered_psth(times, clusters, trials_onsetTimes, trials_offsetTimes, 
                                 trials_IDs, good_units, 
                                 onset_window=[-0.3, 0.5], offset_window=[0, 0.5], 
                                 bin_size=0.005):
    """
    Extract PSTH responses around stimulus onset and offset for each neuron.
    Groups trials by stimulus ID, organizing responses by stimulus identity.
    Onset and offset responses are concatenated along the time axis.
    
    Parameters
    ----------
    times : np.ndarray
        Spike times in seconds
    clusters : np.ndarray
        Cluster IDs for each spike
    trials_onsetTimes : np.ndarray
        Stimulus onset times for each trial
    trials_offsetTimes : np.ndarray
        Stimulus offset times for each trial
    trials_IDs : np.ndarray
        Stimulus IDs for each trial (groups stimulus presentations)
    good_units : array-like
        Unit IDs to include
    onset_window : list
        Time window around stimulus onset [start, end] in seconds
    offset_window : list
        Time window around stimulus offset [start, end] in seconds
    bin_size : float
        Bin size in seconds
    
    Returns
    -------
    psth : np.ndarray
        Shape (n_onset_bins + n_offset_bins, n_stimulus_types, n_units, max_repeats)
        Dimensions: time × stimulus × neurons × repeats
        Time dimension: onset bins [0:n_onset_bins] then offset bins [n_onset_bins:]
    unique_stimulus_ids : np.ndarray
        Unique stimulus IDs in order
    """
    
    if hasattr(good_units, 'squeeze'):
        good_units = good_units.squeeze()
    
    good_units = np.array(good_units)
    n_units = len(good_units)
    
    # Get unique stimulus IDs and bin edges
    unique_stimulus_ids = np.unique(trials_IDs)
    n_stimulus_types = len(unique_stimulus_ids)
    
    onset_bins = np.arange(onset_window[0], onset_window[1] + bin_size, bin_size)
    offset_bins = np.arange(offset_window[0], offset_window[1] + bin_size, bin_size)
    n_onset_bins = len(onset_bins) - 1
    n_offset_bins = len(offset_bins) - 1
    n_time_bins_total = n_onset_bins + n_offset_bins
    
    # Count max repeats per stimulus to determine output size
    max_repeats = 0
    for stim_id in unique_stimulus_ids:
        n_repeats = np.sum(trials_IDs == stim_id)
        max_repeats = max(max_repeats, n_repeats)
    
    # Initialize output: time (onset+offset concatenated) × stimulus × neurons × repeats
    psth = np.zeros((n_time_bins_total, n_stimulus_types, n_units, max_repeats))
    
    # Extract PSTH for each stimulus type
    for stim_idx, stim_id in enumerate(unique_stimulus_ids):
        # Get all trials with this stimulus ID
        stim_trials = np.where(trials_IDs == stim_id)[0]
        
        for repeat_idx, trial_idx in enumerate(stim_trials):
            trial_onset = trials_onsetTimes[trial_idx]
            trial_offset = trials_offsetTimes[trial_idx]
            
            # Extract spikes around onset
            onset_start = trial_onset + onset_window[0]
            onset_end = trial_onset + onset_window[1]
            onset_spike_mask = (times >= onset_start) & (times <= onset_end)
            onset_spike_times = times[onset_spike_mask] - trial_onset
            onset_spike_clusters = clusters[onset_spike_mask]
            
            # Extract spikes around offset
            offset_start = trial_offset + offset_window[0]
            offset_end = trial_offset + offset_window[1]
            offset_spike_mask = (times >= offset_start) & (times <= offset_end)
            offset_spike_times = times[offset_spike_mask] - trial_offset
            offset_spike_clusters = clusters[offset_spike_mask]
            
            # Bin spikes for each unit
            for unit_idx, unit in enumerate(good_units):
                # Onset PSTH (first n_onset_bins of time axis)
                unit_onset_spikes = onset_spike_times[onset_spike_clusters == unit]
                if len(unit_onset_spikes) > 0:
                    hist_onset, _ = np.histogram(unit_onset_spikes, bins=onset_bins)
                    psth[:n_onset_bins, stim_idx, unit_idx, repeat_idx] = hist_onset
                
                # Offset PSTH (remaining n_offset_bins of time axis)
                unit_offset_spikes = offset_spike_times[offset_spike_clusters == unit]
                if len(unit_offset_spikes) > 0:
                    hist_offset, _ = np.histogram(unit_offset_spikes, bins=offset_bins)
                    psth[n_onset_bins:n_onset_bins + n_offset_bins, stim_idx, unit_idx, repeat_idx] = hist_offset
    
    return psth, unique_stimulus_ids

def natim_correlations(param):
    """
    Compute pairwise correlations of natural image stimulus fingerprints across neurons.
    Fingerprints combine average timecourse over stimuli + average response across stimuli.
    """
    session_stimulus_responses, session_timecourses = get_natim_responses(param)
    
    # Concatenate fingerprints from all sessions
    all_timecourses = np.concatenate(session_timecourses, axis=2)  # (fold, time, n_all_units)
    all_stimulus_responses = np.concatenate(session_stimulus_responses, axis=2)  # (fold, stimulus, n_all_units)

    # Compute pairwise correlations for timecourses and stimulus responses
    timecourse_corr = pairwise_histogram_correlation(all_timecourses[0,:,:].T, all_timecourses[1,:,:].T)
    stimulus_corr = pairwise_histogram_correlation(all_stimulus_responses[0,:,:].T, all_stimulus_responses[1,:,:].T)
    
    # Fisher z-transform
    z_timecourse = 0.5 * np.arctanh(timecourse_corr)
    z_stimulus = 0.5 * np.arctanh(stimulus_corr)

    # Average z-scores
    z_fingerprint = (z_timecourse + z_stimulus) / 2.0

    # Inverse transform with tanh
    fingerprint_corr = np.tanh(z_fingerprint)
    return fingerprint_corr

def AUC(matches:np.ndarray, func_metric:np.ndarray, session_id):
    """
    The AUC depends on a functional metric which is considered as ground truth. This is passed in via func_metric.
    """
    
    # Filter out NaN values
    nan_mask = np.all(~np.isfinite(func_metric),axis=1)
    func_metric = func_metric[~nan_mask, :][:, ~nan_mask]
    matches = matches[~nan_mask, :][:, ~nan_mask]
    session_id = session_id[~nan_mask]

    P = np.sum(matches)
    if P < 1:
        raise ValueError("No matches found - can't compute AUC")

    # Calculate AUCs using final sets of matches
    within_session = (session_id[:, None] == session_id).astype(bool)
    func_across = func_metric[~within_session]
    matches_across = matches[~within_session]


    sorted_indices = np.argsort(func_across)[::-1]

    tp, fp = 0,0
    N = len(func_across) - P
    recall, fpr = [], []

    for idx in sorted_indices:
        if matches_across[idx]:
            tp+=1
        else:
            fp+=1
        recall.append(tp/P)
        fpr.append(fp/N)

    # NumPy compatibility: `np.trapezoid` exists in newer NumPy versions;
    # `np.trapz` is the older equivalent.
    if hasattr(np, "trapezoid"):
        auc = np.trapezoid(recall, fpr)
    else:
        auc = np.trapz(recall, fpr)
    return auc


if __name__ == '__main__':

    pass
