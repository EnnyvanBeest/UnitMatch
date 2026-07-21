import os
import sys
from pathlib import Path

# Ensure the DeepUnitMatch package root is on the path so `utils` is
# importable regardless of the caller's working directory.
sys.path.insert(0, str(Path(__file__).parent.parent))

from utils.losses import clip_sim, CustomClipLoss, Projector
from utils.npdataset import (
    NeuropixelsDataset_cortexlab,
    ValidationExperimentBatchSampler,
)
from utils.helpers import read_pos
import numpy as np
from utils.mymodel import SpatioTemporalCNN_V2
import torch
from torch.utils.data import DataLoader
from tqdm import tqdm
import pandas as pd
from sklearn.neighbors import KernelDensity
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


def load_trained_model(device="cpu", read_path=None, n_output=256):

    model = SpatioTemporalCNN_V2(n_channel=30, n_time=60, n_output=n_output).to(device)
    model = model.double()
    if read_path is None:
        current_dir = Path(__file__).parent.parent
        read_path = current_dir / "utils" / "model"
    checkpoint = torch.load(read_path, map_location=device)

    if "n_output" in checkpoint and checkpoint["n_output"] != n_output:
        raise ValueError(
            f"Checkpoint {read_path} was trained with n_output={checkpoint['n_output']}, "
            f"but n_output={n_output} was requested."
        )

    if "clip_loss" in checkpoint:
        # Fine-tuned (clip-loss) checkpoint: checkpoint["model"] is the encoder alone.
        model.load_state_dict(checkpoint["model"])
        clip_loss = CustomClipLoss().to(device)
        clip_loss.load_state_dict(checkpoint["clip_loss"])
        clip_loss.eval()
    else:
        # Autoencoder-only checkpoint: checkpoint["model"] holds encoder.*/decoder.*
        # keys for the full SpatioTemporalAutoEncoder_V2; checkpoint["encoder"] is
        # the encoder-only state dict we actually need here.
        model.load_state_dict(checkpoint["encoder"])
    model.eval()

    # Load projector
    projector = Projector(
        input_dim=256, output_dim=128, hidden_dim=128, n_hidden_layers=1, dropout=0.1
    ).to(device)
    projector = projector.double()

    # Can also return projector if needed
    return model


def reorder_by_depth(matrix: np.ndarray, pos1: np.ndarray, pos2) -> np.ndarray:
    """
    Matrix should compare just one recording session against another.
    """

    depths1 = pos1[:, 1]
    depths2 = pos2[:, 1]

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
    test_sampler = ValidationExperimentBatchSampler(test_dataset, shuffle=False)
    test_loader = DataLoader(test_dataset, batch_sampler=test_sampler)

    submatrices = []
    n_batches = len(test_loader)
    device = next(model.parameters()).device

    for estimates_i, _, positions_i, exp_ids_i, filepaths_i in tqdm(test_loader):
        # Forward pass
        enc_estimates_i = model(estimates_i.to(device))  # shape [bsz, 256]

        for _, candidates_j, positions_j, exp_ids_j, filepaths_j in tqdm(test_loader):
            enc_candidates_j = model(candidates_j.to(device))
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


def get_threshold(prob_matrix: np.ndarray, session_id, MAP=False):

    n = len(session_id)
    within_session = (session_id[:, None] == session_id).astype(int)
    pmat = prob_matrix[within_session == True].reshape(-1, 1)
    labels = np.eye(n)[within_session == True].reshape(-1, 1)

    # On-diagonal means same neuron. Off-diagonal means different neurons.
    on_diag = pmat[labels == 1]
    off_diag = pmat[labels == 0]

    # Kernel density estimation (distributions are more useful than histograms)
    kde_on = KernelDensity(kernel="gaussian", bandwidth=0.01).fit(
        on_diag.reshape(-1, 1)
    )
    kde_off = KernelDensity(kernel="gaussian", bandwidth=0.01).fit(
        off_diag.reshape(-1, 1)
    )
    x = np.linspace(min(off_diag), max(on_diag), 1000).reshape(-1, 1)
    y_on = np.exp(kde_on.score_samples(x))
    y_off = np.exp(kde_off.score_samples(x))

    # Find the threshold where the distributions intersect
    if MAP:
        thresh = np.argwhere(np.diff(np.sign(y_off * (n - 1) - y_on)))
    else:
        thresh = np.argwhere(np.diff(np.sign(y_off - y_on)))
    if len(thresh) == 0:
        thresh = len(x) - 1
    elif len(thresh) > 1:
        thresh = thresh[-1]

    across = prob_matrix[within_session == False].reshape(-1, 1)
    diff = np.median(pmat) - np.median(across)

    return x[thresh].item() - diff


def directional_filter(matches: pd.DataFrame):
    if len(matches) == 0:
        return matches

    # Create a set of all match tuples for O(1) lookup
    match_tuples = set(
        zip(matches["RecSes1"], matches["ID1"], matches["RecSes2"], matches["ID2"])
    )

    # For each match, check if its reverse exists in the set
    # Keep only matches where the reverse direction also exists
    valid_mask = [
        (row["RecSes2"], row["ID2"], row["RecSes1"], row["ID1"]) in match_tuples
        for _, row in matches.iterrows()
    ]

    return matches[valid_mask]


def directional_filter_matrix(sim_matrix, session_id, threshold):
    """
    Matrix version of directional_filter.
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
    best_neuron1_indices = matches.groupby(["RecSes1", "ID1"])[metric].idxmax()

    # Find the best match for each neuron2 (RecSes2, ID2 combination)
    best_neuron2_indices = matches.groupby(["RecSes2", "ID2"])[metric].idxmax()

    # A match is valid only if it's the best match for BOTH neurons involved
    # Take intersection of the two sets of indices
    valid_indices = set(best_neuron1_indices) & set(best_neuron2_indices)

    # Count conflicts (total matches minus valid matches)
    num_conflicts = len(matches) - len(valid_indices)

    # Filter to keep only valid matches
    filtered_matches = matches.loc[list(valid_indices)]

    return filtered_matches, num_conflicts


def get_matches(df, sim_matrix, session_id, data_dir, dist_thresh):
    """
    Process the output probability matrix to get final set of matches across sessions.
    Output is the probability matrix of matches after thresholding and filtering and a boolean matrix indicating final matches.
    """

    sim_thresh = get_threshold(sim_matrix, session_id)

    matches = df.loc[df["Prob"] > sim_thresh].copy()

    matches = matches.loc[
        matches["RecSes1"] != matches["RecSes2"]
    ]  # only keep across-session matches
    matches = directional_filter(matches)

    # spatial filtering
    sessions = np.unique(session_id)
    positions = {}
    for session in sessions:
        session_path = os.path.join(data_dir, str(session))
        positions[session] = read_pos(session_path)
    corrections = get_corrections(matches, positions)
    matches_with_dist = vectorised_drift_corrected_dist(corrections, positions, matches)
    matches = matches_with_dist.loc[matches_with_dist["dist"] < dist_thresh]

    matches = matches[[c for c in df.columns]]
    matches, confs = remove_conflicts(matches, "Prob")

    matches.insert(len(matches.columns), "match", True)
    df.insert(len(df.columns), "match", False)
    df.loc[matches.index, "match"] = True
    df.loc[(df["RecSes1"] == df["RecSes2"]) & (df["ID1"] == df["ID2"]), "match"] = True
    return df


def vectorised_drift_corrected_dist(
    corrections, positions, matches: pd.DataFrame, nocorr=False
):
    """
    Vectorised drift-corrected distance computation.
    Returns matches dataframe with extra column for drift-corrected distance.
    """

    sessions = list(positions.keys())

    # Extract the x, y positions
    idx = matches.index
    coords1 = pd.concat(
        [
            positions[session]
            .loc[:, ["ID", "x", "y"]]
            .rename(columns={"ID": "ID1", "x": "x1", "y": "y1"})
            .assign(RecSes1=session)
            for session in [sessions[0], sessions[1]]
        ]
    )

    coords2 = pd.concat(
        [
            positions[session]
            .loc[:, ["ID", "x", "y"]]
            .rename(columns={"ID": "ID2", "x": "x2", "y": "y2"})
            .assign(RecSes2=session)
            for session in [sessions[1], sessions[0]]
        ]
    )

    # Merge with matches DataFrame
    newmatches = matches.merge(coords1, on=["RecSes1", "ID1"])
    newmatches = newmatches.merge(coords2, on=["RecSes2", "ID2"])
    newmatches = newmatches.drop_duplicates(subset=["ID1", "ID2", "RecSes1", "RecSes2"])
    newmatches.index = idx
    matches = newmatches

    # Calculate shanks (rounded x coordinates)
    matches["shank1"] = matches["x1"].round(-2).astype(int)
    matches["shank2"] = matches["x2"].round(-2).astype(int)
    # Apply corrections where necessary
    if not nocorr:
        matches["idx"] = matches.index
        matches = matches.merge(
            corrections[["rec1", "rec2", "shank", "ydiff"]],
            how="left",
            left_on=["RecSes1", "RecSes2", "shank1"],
            right_on=["rec1", "rec2", "shank"],
        )
        matches.index = matches["idx"]
        matches.drop(["rec1", "rec2", "shank", "idx"], axis=1, inplace=True)
        # Apply the drift correction conditionally
        matches["y2_corr"] = np.where(
            matches["shank1"] != matches["shank2"],
            1000,
            matches["y2"] - matches["ydiff"].fillna(0),
        )
    else:
        matches["y2_corr"] = matches["y2"]
    # Calculate the Euclidean distance
    matches["dist"] = np.sqrt(
        (matches["x1"] - matches["x2"]) ** 2 + (matches["y1"] - matches["y2_corr"]) ** 2
    )
    return matches


def get_corrections(matches, positions):

    pos1 = pd.concat(
        [positions[key].assign(RecSes1=key) for key in positions.keys()],
        ignore_index=True,
    )
    pos2 = pos1.rename(
        columns={"RecSes1": "RecSes2", "x": "x2", "y": "y2", "ID": "ID2"}
    )

    # Merge matches with positions data to get all necessary columns in one DataFrame
    matches = matches.merge(
        pos1, left_on=["RecSes1", "ID1"], right_on=["RecSes1", "ID"], how="left"
    ).merge(pos2, left_on=["RecSes2", "ID2"], right_on=["RecSes2", "ID2"], how="left")

    # Drop rows with missing positions data
    matches = matches.dropna(subset=["x", "y", "x2", "y2"])

    # Calculate shank and ydiff in a vectorized way
    matches["shank1"] = matches["x"].round(-2)
    matches["shank2"] = matches["x2"].round(-2)
    matches["ydiff"] = matches["y2"] - matches["y"]

    # Filter only rows with matching shanks
    matches = matches[matches["shank1"] == matches["shank2"]]

    # Group by session pairs and shank, then calculate median ydiff
    output = (
        matches.groupby(["RecSes1", "RecSes2", "shank1"])["ydiff"]
        .median()
        .reset_index()
        .rename(columns={"shank1": "shank", "RecSes1": "rec1", "RecSes2": "rec2"})
    )
    return output


def pairwise_histogram_correlation(A, B):
    A = A.astype(float)
    B = B.astype(float)
    A -= np.nanmean(A, axis=1, keepdims=True)
    B -= np.nanmean(B, axis=1, keepdims=True)
    # Track rows that are undefined (any NaN remaining, or zero variance)
    undef_A = np.any(np.isnan(A), axis=1) | (
        np.linalg.norm(np.nan_to_num(A), axis=1) == 0
    )
    undef_B = np.any(np.isnan(B), axis=1) | (
        np.linalg.norm(np.nan_to_num(B), axis=1) == 0
    )
    A = np.nan_to_num(A)
    B = np.nan_to_num(B)
    norm_A = np.linalg.norm(A, axis=1, keepdims=True)
    norm_B = np.linalg.norm(B, axis=1, keepdims=True)
    norm_A[norm_A == 0] = 1
    norm_B[norm_B == 0] = 1
    result = (A / norm_A) @ (B / norm_B).T
    result[undef_A, :] = np.nan
    result[:, undef_B] = np.nan
    return result


def get_ISI_histograms(param, ISIbins):

    fs = 3e4
    nclus = param["n_units"]
    KSdirs = param["KS_dirs"]
    good_units = param["good_units"]

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
                        print(
                            f"Warning: Fewer than 50 spikes for neuron {clusid}, please check your inclusion criteria"
                        )
                    # Split idx1 into two halves
                    if cv == 0:
                        idx1 = idx1[: len(idx1) // 2]
                    else:
                        idx1 = idx1[len(idx1) // 2 :]
                    ISIMat[:, cv, index], _ = np.histogram(
                        np.diff(times[idx1].astype(float)), bins=ISIbins
                    )

            index += 1
    return ISIMat


def ISI_correlations(param):

    # Define ISI bins
    ISIbins = np.concatenate(([0], 5 * 10 ** np.arange(-4, 0.1, 0.1)))
    ISIMat = get_ISI_histograms(param, ISIbins)
    A, B = ISIMat[:, 0, :].T, ISIMat[:, 1, :].T  # pull out each cross-validation fold

    return pairwise_histogram_correlation(A, B)  # compute pairwise correlations


def get_binned_psth(param, bin_size=0.01):
    """
    Bin spike times into a population activity matrix per session.
    Returns list of (n_units_session x n_time_bins) arrays, one per session.
    bin_size is in seconds (default 10 ms, matching MATLAB UMparam.binsz).
    """
    fs = 3e4
    KSdirs = param["KS_dirs"]
    good_units = param["good_units"]

    session_psthas = []
    for session in range(len(good_units)):
        times = np.load(os.path.join(KSdirs[session], "spike_times.npy")) / fs
        clusters = np.load(os.path.join(KSdirs[session], "spike_clusters.npy"))

        edges = np.arange(times.min() - bin_size / 2, times.max() + bin_size, bin_size)

        session_units = good_units[session]
        if hasattr(session_units, "squeeze"):
            session_units = session_units.squeeze()

        sr = np.zeros((len(session_units), len(edges) - 1))
        for uid, clusid in enumerate(session_units):
            idx = np.where(clusters == clusid)[0]
            if idx.size > 0:
                sr[uid, :], _ = np.histogram(times[idx], bins=edges)

        session_psthas.append(sr)

    return session_psthas


def _pairwise_corr_cols(A, B):
    """
    Pairwise Pearson r between all column pairs of A (n×p) and B (n×q),
    equivalent to MATLAB corr(A, B, 'rows', 'pairwise').

    NaN entries (at most one per column in our use case) are replaced by the
    column mean before correlating — a negligible approximation when NaNs are
    sparse.  Columns that are entirely NaN or have zero variance yield 0.
    Returns (p, q).
    """
    A = A.copy().astype(float)
    B = B.copy().astype(float)
    for M in (A, B):
        for j in range(M.shape[1]):
            nan_mask = np.isnan(M[:, j])
            if nan_mask.all():
                M[:, j] = 0.0  # all-NaN column → constant zero
            elif nan_mask.any():
                M[nan_mask, j] = np.nanmean(M[:, j])
    p = A.shape[1]
    combined = np.vstack([A.T, B.T])  # (p+q) × n
    with np.errstate(invalid="ignore", divide="ignore"):
        C = np.corrcoef(combined)  # (p+q) × (p+q)
    return np.nan_to_num(C[:p, p:], nan=0.0)


def refpop_correlations(param, matches=None, bin_size=0.01):
    """
    Cross-correlation population fingerprint (refPopCorr), following
    MATLAB ComputeFunctionalScores / CrossCorrelationFingerPrint.

    Within-session: fingerprint is the cross-validated pairwise correlation
    (fold1 vs fold2), i.e. MATLAB corr(fold1, fold2, 'rows', 'pairwise').

    Cross-session: the reference population is restricted to matched units
    (equal-length by definition of a match). Each unit's fingerprint is its
    correlation profile with those matched reference units in its own session;
    these equal-length vectors are then compared across sessions.

    Parameters
    ----------
    param : dict
        Must contain 'KS_dirs', 'good_units', 'n_units'.
    matches : np.ndarray (n_units, n_units) bool, optional
        Final match matrix. Required for cross-session blocks; cross-session
        blocks are left as zero when not provided.
    bin_size : float
        PSTH bin size in seconds (default 10 ms, matching UMparam.binsz).
    """
    session_psthas = get_binned_psth(param, bin_size)
    good_units = param["good_units"]
    n_sessions = len(good_units)
    eps = 1e-7

    session_fold1, session_fold2, session_avg = [], [], []
    for sr in session_psthas:
        n_fold = sr.shape[1] // 2
        C1 = np.corrcoef(sr[:, :n_fold])
        C2 = np.corrcoef(sr[:, n_fold : 2 * n_fold])
        # Zero out NaN from units with no spikes before setting diagonal
        C1 = np.nan_to_num(C1, nan=0.0)
        C2 = np.nan_to_num(C2, nan=0.0)
        # z-transform average of the two folds (MATLAB: tanh(nanmean(atanh(...))))
        avg = np.tanh(
            0.5
            * (
                np.arctanh(np.clip(C1, -1 + eps, 1 - eps))
                + np.arctanh(np.clip(C2, -1 + eps, 1 - eps))
            )
        )
        # NaN the diagonal — self-correlation is excluded from fingerprints
        np.fill_diagonal(C1, np.nan)
        np.fill_diagonal(C2, np.nan)
        np.fill_diagonal(avg, np.nan)
        session_fold1.append(C1)
        session_fold2.append(C2)
        session_avg.append(avg)

    unit_session = []
    for si, units in enumerate(good_units):
        if hasattr(units, "squeeze"):
            units = units.squeeze()
        unit_session.extend([si] * len(units))
    unit_session = np.array(unit_session)

    nclus = param["n_units"]
    refPopCorr = np.zeros((nclus, nclus))

    for si in range(n_sessions):
        rows_i = np.where(unit_session == si)[0]
        n_si = len(rows_i)

        for sj in range(n_sessions):
            rows_j = np.where(unit_session == sj)[0]
            n_sj = len(rows_j)

            if si == sj:
                # Cross-validated within-session fingerprint:
                # corr(fold1, fold2, 'rows', 'pairwise') — n_si observations
                block = _pairwise_corr_cols(session_fold1[si], session_fold2[si])

            else:
                if matches is None:
                    block = np.zeros((n_si, n_sj))
                else:
                    # Restrict reference population to matched units only.
                    # Both sessions contribute exactly the same number of matched
                    # units (one per pair), so fingerprint vectors are equal-length.
                    match_block = matches[np.ix_(rows_i, rows_j)]  # (n_si × n_sj)
                    match_pos = np.argwhere(match_block)  # (n_matched × 2)

                    if len(match_pos) == 0:
                        block = np.zeros((n_si, n_sj))
                    else:
                        local_si = match_pos[
                            :, 0
                        ]  # local indices of matched units in si
                        local_sj = match_pos[
                            :, 1
                        ]  # local indices of matched units in sj

                        # Each row is one matched unit's correlation with all session units.
                        # Diagonal NaN (self-correlation) already set via fill_diagonal above.
                        SC1 = session_avg[si][local_si, :]  # (n_matched × n_si)
                        SC2 = session_avg[sj][local_sj, :]  # (n_matched × n_sj)

                        # corr(SC1, SC2, 'rows', 'pairwise') → (n_si × n_sj)
                        block = _pairwise_corr_cols(SC1, SC2)

            refPopCorr[np.ix_(rows_i, rows_j)] = block

    return refPopCorr


def get_FR(param):
    """
    Compute mean firing rate (spikes/s) per unit for two cross-validation halves.
    Returns array of shape (2, n_units): FR[fold, unit].
    """
    fs = 3e4
    nclus = param["n_units"]
    KSdirs = param["KS_dirs"]
    good_units = param["good_units"]

    FR = np.zeros((2, nclus))
    index = 0
    for session in range(len(good_units)):
        times = np.load(os.path.join(KSdirs[session], "spike_times.npy")) / fs
        clusters = np.load(os.path.join(KSdirs[session], "spike_clusters.npy"))

        session_units = good_units[session]
        if hasattr(session_units, "squeeze"):
            session_units = session_units.squeeze()

        for clusid in session_units:
            idx = np.where(clusters == clusid)[0]
            if idx.size > 0:
                n = len(idx)
                for cv, cv_idx in enumerate([idx[: n // 2], idx[n // 2 :]]):
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
    FR = get_FR(param)  # (2, n_units)
    return np.abs(FR[1, :, np.newaxis] - FR[0, np.newaxis, :])


def get_ISI_CV(param):
    """
    Compute the coefficient of variation (CV = std/mean) of interspike intervals
    per unit for two cross-validation halves.
    Returns array of shape (2, n_units): CV[fold, unit].
    Units with fewer than 2 spikes in a fold get CV = 0.
    """
    fs = 3e4
    nclus = param["n_units"]
    KSdirs = param["KS_dirs"]
    good_units = param["good_units"]

    CV = np.zeros((2, nclus))
    index = 0
    for session in range(len(good_units)):
        times = np.load(os.path.join(KSdirs[session], "spike_times.npy")) / fs
        clusters = np.load(os.path.join(KSdirs[session], "spike_clusters.npy"))

        session_units = good_units[session]
        if hasattr(session_units, "squeeze"):
            session_units = session_units.squeeze()

        for clusid in session_units:
            idx = np.where(clusters == clusid)[0]
            if idx.size > 1:
                n = len(idx)
                for cv, cv_idx in enumerate([idx[: n // 2], idx[n // 2 :]]):
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
    CV = get_ISI_CV(param)  # (2, n_units)
    return np.abs(CV[1, :, np.newaxis] - CV[0, np.newaxis, :])


def get_natim_responses(param):
    """
    Bin spike times into stimulus-triggered PSTH per session.
    Returns list of arrays with shape (n_time_bins, n_stimuli, n_units, max_repeats).
    Organized as: time × stimulus × neurons × repeats
    Time dimension: onset bins followed by offset bins (concatenated).
    """
    fs = 3e4
    nclus = param["n_units"]
    KSdirs = param["KS_dirs"]
    good_units = param["good_units"]

    session_timecourses = []
    session_stimulus_responses = []
    for session in range(len(good_units)):
        times = np.load(os.path.join(KSdirs[session], "spike_times.npy")) / fs
        clusters = np.load(os.path.join(KSdirs[session], "spike_clusters.npy"))

        session_units = good_units[session]
        if hasattr(session_units, "squeeze"):
            session_units = session_units.squeeze()
        n_units = len(session_units)

        # Try to load trial files; if they don't exist, use NaN arrays
        trial_path = os.path.dirname(os.path.dirname(KSdirs[session]))
        try:
            trials_IDs = np.load(os.path.join(trial_path, "trial.imageIDs.npy"))
            trials_offsetTimes = np.load(
                os.path.join(trial_path, "trial.offsetTimes.npy")
            )
            trials_onsetTimes = np.load(
                os.path.join(trial_path, "trial.onsetTimes.npy")
            )
        except FileNotFoundError:
            # Trial files don't exist; return NaN arrays
            bin_size = 0.005
            onset_window = [-0.3, 0.5]
            offset_window = [0.0, 0.5]
            onset_bins = np.arange(
                onset_window[0], onset_window[1] + bin_size, bin_size
            )
            offset_bins = np.arange(
                offset_window[0], offset_window[1] + bin_size, bin_size
            )
            n_time_bins_total = (len(onset_bins) - 1) + (len(offset_bins) - 1)

            timecourse = np.full((2, n_time_bins_total, n_units), np.nan)
            stimulus_response = np.full((2, 112, n_units), np.nan)

            session_timecourses.append(timecourse)
            session_stimulus_responses.append(stimulus_response)
            continue

        bin_size = 0.005
        onset_window = [-0.3, 0.5]
        offset_window = [0.0, 0.5]
        psth, _ = get_stimulus_triggered_psth(
            times,
            clusters,
            trials_onsetTimes,
            trials_offsetTimes,
            trials_IDs,
            session_units,
            onset_window=onset_window,
            offset_window=offset_window,
            bin_size=bin_size,
        )

        # psth shape: (time, stimulus, neurons, repeats)
        # Average timecourse over stimuli and repeats: (fold, time, neurons)
        timecourse1 = np.nanmean(psth[:, :, :, 1::2], axis=(1, 3))
        timecourse2 = np.nanmean(psth[:, :, :, ::2], axis=(1, 3))
        timecourse = np.stack(
            [timecourse1, timecourse2], axis=0
        )  # shape: (fold, time, neurons)

        # Average response across repeats for each stimulus, using only first 200ms after onset
        # onset_window = [-0.3, 0.5], we want time 0 to 0.2s
        ms_after_onset = 0.2
        bins_in_window = int(ms_after_onset / bin_size)  # 40 bins
        onset_start_bin = int(
            (0 - onset_window[0]) / bin_size
        )  # skip to time 0 (60 bins)
        onset_end_bin = onset_start_bin + bins_in_window  # 100
        stimulus_response1 = np.nanmean(
            psth[onset_start_bin:onset_end_bin, :, :, 1::2], axis=(0, 3)
        )
        stimulus_response2 = np.nanmean(
            psth[onset_start_bin:onset_end_bin, :, :, ::2], axis=(0, 3)
        )
        stimulus_response = np.stack(
            [stimulus_response1, stimulus_response2], axis=0
        )  # shape: (fold, stimulus, neurons)

        session_timecourses.append(timecourse)
        session_stimulus_responses.append(stimulus_response)

    return session_stimulus_responses, session_timecourses


def get_stimulus_triggered_psth(
    times,
    clusters,
    trials_onsetTimes,
    trials_offsetTimes,
    trials_IDs,
    good_units,
    onset_window=[-0.3, 0.5],
    offset_window=[0, 0.5],
    bin_size=0.005,
):
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

    if hasattr(good_units, "squeeze"):
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
                    psth[
                        n_onset_bins : n_onset_bins + n_offset_bins,
                        stim_idx,
                        unit_idx,
                        repeat_idx,
                    ] = hist_offset

    return psth, unique_stimulus_ids


def natim_correlations(param):
    """
    Compute pairwise correlations of natural image stimulus fingerprints across neurons.
    Fingerprints combine average timecourse over stimuli + average response across stimuli.
    """
    session_stimulus_responses, session_timecourses = get_natim_responses(param)

    # Concatenate fingerprints from all sessions
    all_timecourses = np.concatenate(
        session_timecourses, axis=2
    )  # (fold, time, n_all_units)
    all_stimulus_responses = np.concatenate(
        session_stimulus_responses, axis=2
    )  # (fold, stimulus, n_all_units)

    # Compute pairwise correlations for timecourses and stimulus responses
    timecourse_corr = pairwise_histogram_correlation(
        all_timecourses[0, :, :].T, all_timecourses[1, :, :].T
    )
    stimulus_corr = pairwise_histogram_correlation(
        all_stimulus_responses[0, :, :].T, all_stimulus_responses[1, :, :].T
    )

    # Fisher z-transform
    z_timecourse = 0.5 * np.arctanh(timecourse_corr)
    z_stimulus = 0.5 * np.arctanh(stimulus_corr)

    # Average z-scores
    z_fingerprint = (z_timecourse + z_stimulus) / 2.0

    # Inverse transform with tanh
    fingerprint_corr = np.tanh(z_fingerprint)
    return fingerprint_corr


def AUC(matches: np.ndarray, func_metric: np.ndarray, session_id):
    """
    The AUC depends on a functional metric which is considered as ground truth. This is passed in via func_metric.
    """

    # Filter out NaN values
    nan_mask = np.all(~np.isfinite(func_metric), axis=1)
    func_metric = func_metric[~nan_mask, :][:, ~nan_mask]
    matches = matches[~nan_mask, :][:, ~nan_mask]
    session_id = session_id[~nan_mask]

    # Only across-session pairs are meaningful for cross-session AUC.
    # P must be computed from across-session matches so that recall and FPR
    # are correctly normalised (DeepUnitMatch already filters to across-session
    # before calling this; UMPy's raw threshold matrix includes within-session
    # self-matches that would otherwise inflate P and deflate recall).
    within_session = (session_id[:, None] == session_id).astype(bool)
    func_across = func_metric[~within_session]
    matches_across = matches[~within_session]

    P = np.sum(matches_across)
    if P < 1:
        raise ValueError("No matches found - can't compute AUC")

    sorted_indices = np.argsort(func_across)[::-1]

    tp, fp = 0, 0
    N = len(func_across) - P
    recall, fpr = [], []

    for idx in sorted_indices:
        if matches_across[idx]:
            tp += 1
        else:
            fp += 1
        recall.append(tp / P)
        fpr.append(fp / N)

    # NumPy compatibility: `np.trapezoid` exists in newer NumPy versions;
    # `np.trapz` is the older equivalent.
    if hasattr(np, "trapezoid"):
        auc = np.trapezoid(recall, fpr)
    else:
        auc = np.trapz(recall, fpr)
    return auc


if __name__ == "__main__":
    pass
