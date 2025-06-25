import os, sys

sys.path.insert(0, os.getcwd())
sys.path.insert(0, os.path.join(os.getcwd(), os.pardir))

from utils.losses import clip_sim, CustomClipLoss, Projector
from utils.npdataset import NeuropixelsDataset, ValidationExperimentBatchSampler
import numpy as np
import matplotlib.pyplot as plt
from utils.mymodel import SpatioTemporalCNN_V2
import torch
from torch.utils.data import DataLoader
from tqdm import tqdm
import pandas as pd
from pathlib import Path
from sklearn.neighbors import KernelDensity
import scipy

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

    for estimates_i, _, positions_i, exp_ids_i, filepaths_i in tqdm(test_loader):
        # Forward pass
        enc_estimates_i = model(estimates_i)        # shape [bsz, 256]
        positions.append(positions_i)

        for _, candidates_j, positions_j, exp_ids_j, filepaths_j in tqdm(test_loader):
            enc_candidates_j = model(candidates_j)
            s = clip_sim(enc_estimates_i, enc_candidates_j)
            submatrices.append(reorder_by_depth(s.detach().cpu().numpy(), positions_i, positions_j))
    
    result = np.vstack((np.hstack((submatrices[0], submatrices[1])), np.hstack((submatrices[2], submatrices[3]))))

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

    return x[thresh].item()

def directional_filter(matrix: np.ndarray):
    """
    Filter a 2D numpy array to keep only bidirectional relationships.
    If matrix[i][j] is non-zero, matrix[j][i] must also be non-zero, 
    otherwise both are set to zero.
    """
    nonzero_mask = matrix != 0
    nonzero_transpose_mask = matrix.T != 0
    
    # Keep only elements where both directions are non-zero
    bidirectional_mask = nonzero_mask & nonzero_transpose_mask
    
    # return the bidirectional mask
    return bidirectional_mask

def remove_conflicts(matrix, fill_value=0):
    """Optimal conflict removal using Hungarian algorithm."""
    
    m = np.array(matrix)
    cost_matrix = np.where(m != fill_value, -1*m, 0)  # negative for maximisation
    rows, cols = scipy.optimize.linear_sum_assignment(cost_matrix)
    
    result = np.full_like(m, fill_value)
    valid = m[rows, cols] != fill_value
    result[rows[valid], cols[valid]] = m[rows[valid], cols[valid]]
    
    return result

def get_final_matches(output_prob_matrix, session_id, matching_threshold=0.5):
    """
    Process the output probability matrix to get final set of matches across sessions.
    Output is the probability matrix of matches after thresholding and filtering and a boolean matrix indicating final matches.
    """
    within_session = (session_id[:, None] == session_id).astype(int)
    probs = output_prob_matrix.copy()
    probs[within_session == True] = 0                               # don't include within-session pairs
    matches = probs > matching_threshold
    bidirectional_mask = directional_filter(matches)
    probs[~bidirectional_mask] = 0                                  # make the matching critera more stringent
    probs = remove_conflicts(probs)                    # ensure each neuron is only matched once
    final_matches = probs > matching_threshold
    return probs, final_matches

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

def get_ISI_histograms(param, good_units, ISIbins):

    fs = 3e4
    nclus = param['n_units']
    KSdirs = param['KS_dirs']
    spktimes1 = np.load(os.path.join(KSdirs[0], "spike_times.npy")) / fs
    spktimes2 = np.load(os.path.join(KSdirs[1], "spike_times.npy")) / fs
    spktimes = [spktimes1, spktimes2]
    spkclus1 = np.load(os.path.join(KSdirs[0], "spike_clusters.npy"))
    spkclus2 = np.load(os.path.join(KSdirs[1], "spike_clusters.npy"))
    spkclus = [spkclus1, spkclus2]

    ISIMat = np.zeros((len(ISIbins) - 1, 2, nclus))
    index = 0
    for session in range(2):
        times = spktimes[session]
        clusters = spkclus[session]

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

def ISI_correlations(param, good_units):

    # Define ISI bins
    ISIbins = np.concatenate(([0], 5 * 10 ** np.arange(-4, 0.1, 0.1)))
    ISIMat = get_ISI_histograms(param, good_units, ISIbins)
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
    auc = np.trapz(recall, fpr)
    return auc


if __name__ == '__main__':

    pass
