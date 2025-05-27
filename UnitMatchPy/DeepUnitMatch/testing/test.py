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
# from testing.isi_corr import plot_edf5, get_corrections, vectorized_drift_corrected_dist, remove_conflicts, avg_across_directions
# from testing.func_match import func_matches, get_matches, AUC

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

# def spatial_filter(mt_path:str, matches:pd.DataFrame, dist_thresh=None, drift_corr=True, plot_drift=False):
#     """
#     Input is a dataframe of potential matches (according to some threshold, e.g. DNNSim)
#     Output is a reduced dataframe after filtering out matches that are spatially distant.
#     dist_thresh can take a numerical value, or set it to None to just reject the 50% of matches 
#     with greatest Euclidean distance.

#     Args:
#         mt_path (str): Path to the matchtable.
#         matches (pd.DataFrame): DataFrame containing potential matches.
#         dist_thresh (Optional[float]): Maximum allowed Euclidean distance. If None, filters out the 50% most distant matches.
#         drift_corr (bool): Whether to apply drift correction. Defaults to True.
#         plot_drift (bool): Whether to plot drift correction visualization. Defaults to False.
#     """
#     if len(matches) < 1:
#         return matches
#     exp_ids, metadata = mtpath_to_expids(mt_path, matches)
#     test_data_root = mt_path[:mt_path.find(metadata["mouse"])]
#     positions = {}
#     for recses, exp_id in exp_ids.items():
#         fp = os.path.join(test_data_root, metadata["mouse"], metadata["probe"], 
#                           metadata["loc"], exp_id, "processed_waveforms")
#         pos_dict = read_pos(fp)
#         positions[recses] = pd.DataFrame(pos_dict)
#     if drift_corr:
#         corrections = get_corrections(matches, positions)
#     # plot_distances(matches, positions, corrections=corrections)
#     matches_with_dist = vectorized_drift_corrected_dist(corrections, positions, matches)
#     if not dist_thresh:
#         matches_with_dist.sort_values(by = "dist", inplace=True)
#         return matches_with_dist.head(len(matches)//2)
#     else:
#         return matches_with_dist.loc[matches_with_dist["dist"]<dist_thresh]

def directional_filter(matches: pd.DataFrame):
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


if __name__ == '__main__':

    pass
