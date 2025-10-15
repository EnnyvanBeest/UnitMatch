import os, sys

sys.path.insert(0, os.getcwd())
sys.path.insert(0, os.path.join(os.getcwd(), os.pardir))

from utils.losses import clip_sim, CustomClipLoss, Projector
from utils.npdataset import NeuropixelsDataset, ValidationExperimentBatchSampler
from utils.helpers import read_pos
import numpy as np
import matplotlib.pyplot as plt
from utils.mymodel import SpatioTemporalCNN_V2
import torch
from torch.utils.data import DataLoader
from tqdm import tqdm
import pandas as pd
from pathlib import Path
from sklearn.neighbors import KernelDensity
import importlib
from utils import helpers
importlib.reload(helpers)
from utils.helpers import *



def load_trained_model(device="cpu"):

    model = SpatioTemporalCNN_V2(n_channel=30,n_time=60,n_output=256).to(device)
    model = model.double()
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

def inference(model, data_dir):

    test_dataset = NeuropixelsDataset(data_dir)
    test_sampler = ValidationExperimentBatchSampler(test_dataset, shuffle = False)
    test_loader = DataLoader(test_dataset, batch_sampler=test_sampler)

    submatrices = []
    positions = []
    n_batches = len(test_loader)

    for estimates_i, _, positions_i, exp_ids_i, filepaths_i in tqdm(test_loader):
        # Forward pass
        enc_estimates_i = model(estimates_i)        # shape [bsz, 256]
        positions.append(positions_i)

        for _, candidates_j, positions_j, exp_ids_j, filepaths_j in tqdm(test_loader):
            enc_candidates_j = model(candidates_j)
            s = clip_sim(enc_estimates_i, enc_candidates_j)
            submatrices.append(reorder_by_depth(s.detach().cpu().numpy(), positions_i, positions_j))
    
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

def get_matches(df, sim_matrix, session_id, data_dir, pos_array, dist_thresh):
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
        positions[session]['ID'] = positions[session]['ID'].values - positions[session]['ID'].min()
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

def AUC(matches:np.ndarray, func_metric:np.ndarray, session_id):
    """
    The AUC depends on a functional metric which is considered as ground truth. This is passed in via func_metric.
    """

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
    auc = np.trapezoid(recall, fpr)
    return auc


if __name__ == '__main__':

    pass
