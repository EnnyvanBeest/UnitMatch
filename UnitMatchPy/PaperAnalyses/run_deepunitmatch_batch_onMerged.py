# Batch wrapper: runs DeepUnitMatch and UMPy on the *merged* dataset produced by
# generate_merged_dataset.py.
#
# generate_merged_dataset.py builds a fully self-contained merged tree: for every
# original (pre-merge) DeepUnitMatch output folder
#    DUM_NONMERGED_DATAPATH/<X>/DeepUnitMatch
# (one UMparam.pickle-holding folder per probe/depth group) it mirrors *every*
# entry in that group's KS_dirs list — not just the ones that had good units in
# the original comparison — into
#    BASE_INPUT/<X>/DeepUnitMatch/<0, 1, 2, ...>/
# each numbered subfolder getting that original session's KS-style files
# (spike_times.npy, channel_positions.npy, cluster_info.tsv, and a
# cluster_bc_unitType.tsv that generate_merged_dataset.py always *derives* from
# that group's UnitMatch.mat GoodID rather than copying the session's own raw
# bombcell output — see write_synthetic_bc_unit_type_tsv() there for why) plus
# the merged RawSpikes waveforms (qMetrics/RawWaveforms) and a merged
# spike_clusters.npy (units merged within a session are folded together and the
# losing unit's RawSpikes file is deleted).
#
# We therefore treat each BASE_INPUT/<X>/DeepUnitMatch folder as a complete,
# self-contained KS-style dataset in its own right: every numbered subfolder is
# a session, and good units for that session come straight from its
# cluster_bc_unitType.tsv as generated above. This script itself never consults
# UnitMatch.mat/UMparam.pickle directly — that happened once, upstream, at
# generation time.
# A unit whose TSV row still says GOOD but whose RawSpikes file is gone (it lost
# a within-session merge) is simply skipped when loading waveforms.
#
# Results are mirrored to:
#    \\znas.cortexlab.net\Lab\Share\UNITMATCHTABLES_ENNY_CELIAN_JULIE\DeepUM_NatMeth2026_V3_OnMergedData
# with the same subfolder structure, split into DeepUnitMatch/ and UMPy/ subfolders.
# Waveforms are loaded once per group and shared between both pipelines.

import os
import sys
import copy
import datetime
import traceback
import numpy as np
import matplotlib

matplotlib.use("Agg")  # non-interactive backend for batch runs
import matplotlib.pyplot as plt

# ── project paths ───────────────────────────────────────────────────────────
_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _HERE)
sys.path.insert(0, os.path.dirname(_HERE))
sys.path.insert(0, os.path.join(_HERE, "DeepUnitMatch"))

import argparse

import batch_lock
import UnitMatchPy.default_params as default_params
import UnitMatchPy.utils as util
import UnitMatchPy.overlord as ov
import UnitMatchPy.bayes_functions as bf
import UnitMatchPy.assign_unique_id as aid
import UnitMatchPy.save_utils as su
import UnitMatchPy.metric_functions as mf
from DeepUnitMatch.utils import param_fun
from DeepUnitMatch.testing import test
from DeepUnitMatch.utils import helpers

try:
    from convert_python_batch_output_to_matlab import convert_python_output_to_matlab
except Exception:
    convert_python_output_to_matlab = None

# ── user settings ────────────────────────────────────────────────────────────
BASE_INPUT = r"\\znas.cortexlab.net\Lab\Share\UNITMATCHTABLES_ENNY_CELIAN_JULIE\DeepUM_NatMeth2026V2_merged\merged_data_v2"
BASE_OUTPUT = r"\\znas.cortexlab.net\Lab\Share\UNITMATCHTABLES_ENNY_CELIAN_JULIE\DeepUM_NatMeth2026_V3_OnMergedData"


DEVICE = "cuda" if test.torch.cuda.is_available() else "cpu"
print(f"Device: {DEVICE}")
THRESH = 0.5
# Sentinel-based skip check: a group is treated as done when its
# MatchingOverview.png exists AND (REDO_FROM_DATE is None, or the file is at
# least that new). Set REDO_FROM_DATE to the date a fix landed (e.g. the DUM
# adaptive-prior/threshold fix) so only output computed before that date gets
# redone -- results a concurrent/earlier run already redid since then are
# correctly recognised as up to date instead of being reprocessed forever.
# Set to None to fall back to plain "skip if present". Set to a far-future
# date for the old unconditional REDO=True behaviour.
REDO_FROM_DATE = datetime.datetime(2026, 7, 22, 19, 0, 0)
WRITE_MATLAB_COMPAT = False


# ── good-unit helpers ─────────────────────────────────────────────────────────


def build_good_units_from_labels(unit_label_paths):
    """
    Determine good units per session directly from each session's own unit-label
    TSV (as copied into the merged tree by generate_merged_dataset.py), using the
    same GOOD / NON-SOMA GOOD (bombcell) or 'good' (cluster_group.tsv) criteria as
    util.load_good_waveforms.

    Returns a list (one entry per session) of (N, 1) float arrays of cluster IDs.
    """
    good_units = []
    for path in unit_label_paths:
        labels = util.load_tsv(path)
        if os.path.basename(path) == "cluster_bc_unitType.tsv":
            idx = np.isin(labels[:, 1], ["GOOD", "NON-SOMA GOOD"])
        else:
            idx = labels[:, 1] == "good"
        ids = labels[idx, 0].astype(float)
        good_units.append(ids.reshape(-1, 1))
    return good_units


# ── merged-group discovery ───────────────────────────────────────────────────


def find_merged_groups():
    """
    Find merged-data groups: every 'DeepUnitMatch' folder under BASE_INPUT whose
    immediate children are numbered session folders (0, 1, 2, ...), as written by
    generate_merged_dataset.py. Each such folder is treated as a complete,
    self-contained KS-style dataset.

    Returns a list of merged_dir paths.
    """
    groups = []
    for root, dirs, _ in os.walk(BASE_INPUT):
        if os.path.basename(root) == "DeepUnitMatch" and any(d.isdigit() for d in dirs):
            groups.append(root)
    return groups


def load_waveforms_for_good_units(wave_paths, good_units_per_session, param):
    """
    Directly loads RawSpikes waveforms for the specified good units, bypassing
    the TSV-based unit-label files used by util.load_good_waveforms.

    Units with no RawSpikes file (e.g. merged away by generate_merged_dataset.py,
    which deletes the losing unit's waveform file after folding it into the
    surviving unit) are silently skipped rather than assumed present — the
    surviving unit keeps its original ID with the merged waveform.

    Returns the same tuple as util.load_good_waveforms.
    """
    n_sessions = len(wave_paths)
    waveforms = []
    actual_good_units = []
    successful_sessions = []

    for sess_idx in range(n_sessions):
        wave_path = wave_paths[sess_idx]
        g_units = good_units_per_session[sess_idx].flatten()

        if len(g_units) == 0:
            print(f"  Session {sess_idx}: no good units, skipping.")
            continue

        try:
            buf = None
            kept_ids = []
            kept_idx = []
            for j, uid in enumerate(g_units):
                p = os.path.join(wave_path, f"Unit{int(uid)}_RawSpikes.npy")
                if os.path.exists(p):
                    w = np.load(p)
                    if buf is None:
                        buf = np.zeros((len(g_units), *w.shape), dtype=w.dtype)
                    buf[j] = w
                    kept_ids.append(uid)
                    kept_idx.append(j)
                else:
                    print(f"  Note: {p} missing (unit likely merged away), skipping.")

            if not kept_ids:
                print(f"  Session {sess_idx}: no waveform files found, skipping.")
                continue

            buf = buf[kept_idx]
            waveforms.append(buf)
            actual_good_units.append(np.array(kept_ids, dtype=float).reshape(-1, 1))
            successful_sessions.append(sess_idx)

        except Exception as e:
            print(f"  Error loading session {sess_idx}: {e}")
        finally:
            try:
                del buf
            except NameError:
                pass

    if not waveforms:
        raise RuntimeError("No sessions loaded successfully.")

    if len(successful_sessions) < n_sessions:
        failed = [i for i in range(n_sessions) if i not in successful_sessions]
        print(
            f"  Warning: skipped {len(failed)} session(s) with no loadable waveforms: {failed}"
        )

    waveform = np.concatenate(waveforms, axis=0)
    n_units_per_session = np.array([w.shape[0] for w in waveforms], dtype=int)

    param["n_units"], session_id, session_switch, param["n_sessions"] = (
        util.get_session_data(n_units_per_session)
    )
    within_session = util.get_within_session(session_id, param)

    param["n_channels"] = waveform.shape[2]
    param["n_units_per_session"] = [len(g) for g in actual_good_units]

    actual_width = waveform.shape[1]
    param["spike_width"] = actual_width
    param["peak_loc"] = int(np.floor(actual_width / 2))
    param["waveidx"] = np.arange(
        param["peak_loc"] - 8, param["peak_loc"] + 15, dtype=int
    )

    return (
        waveform,
        session_id,
        session_switch,
        within_session,
        actual_good_units,
        param,
    )


# ── path helpers ─────────────────────────────────────────────────────────────


def get_save_dir(merged_dir):
    """Return the DeepUnitMatch output directory for a given merged-data group dir."""
    subfolder = os.path.relpath(os.path.dirname(merged_dir), BASE_INPUT)
    return os.path.join(BASE_OUTPUT, subfolder, "DeepUnitMatch")


def get_umpy_save_dir(merged_dir):
    """Return the UMPy output directory for a given merged-data group dir."""
    subfolder = os.path.relpath(os.path.dirname(merged_dir), BASE_INPUT)
    return os.path.join(BASE_OUTPUT, subfolder, "UMPy")


def results_exist(merged_dir):
    """Return True when the DeepUnitMatch sentinel output file is present and fresh (see REDO_FROM_DATE)."""
    sentinel = os.path.join(get_save_dir(merged_dir), "MatchingOverview.png")
    return batch_lock.sentinel_is_fresh(sentinel, REDO_FROM_DATE)


def ensure_matlab_compatible_output(save_dir):
    """Create UnitMatch.mat in save_dir from Python outputs when explicitly requested."""
    if not WRITE_MATLAB_COMPAT:
        return None

    out_path = os.path.join(save_dir, "UnitMatch.mat")
    if os.path.exists(out_path):
        return out_path

    if convert_python_output_to_matlab is None:
        print(
            "  WARNING: MATLAB-compatible conversion requested but the converter is unavailable."
        )
        return None

    required = [
        os.path.join(save_dir, "UMparam.pickle"),
        os.path.join(save_dir, "ClusInfo.pickle"),
        os.path.join(save_dir, "MatchProb.npy"),
        os.path.join(save_dir, "MatchTable.csv"),
        os.path.join(save_dir, "WaveformInfo.npz"),
    ]
    if not all(os.path.exists(p) for p in required):
        return None

    try:
        return convert_python_output_to_matlab(save_dir, output_mat_path=out_path)
    except Exception as e:
        print(f"  WARNING: could not create MATLAB-compatible output: {e}")
        return None


def umpy_results_exist(merged_dir):
    """Return True when the UMPy sentinel output file is present and fresh (see REDO_FROM_DATE)."""
    sentinel = os.path.join(get_umpy_save_dir(merged_dir), "MatchingOverview.png")
    return batch_lock.sentinel_is_fresh(sentinel, REDO_FROM_DATE)


def get_group_lock_path(merged_dir):
    """
    Lock file marking 'a run is currently processing this group' (DeepUnitMatch
    + UMPy together), so multiple machines pointed at the same BASE_INPUT/
    BASE_OUTPUT can split work across groups without double-processing one.
    See batch_lock.py. Named distinctly from the extramodels script's lock so
    the two scripts can run on the same group concurrently without blocking
    each other.
    """
    subfolder = os.path.relpath(os.path.dirname(merged_dir), BASE_INPUT)
    return os.path.join(BASE_OUTPUT, subfolder, ".processing.lock")


# ── shared session loader ─────────────────────────────────────────────────────


def _prepare_session(merged_dir):
    """
    Load and validate everything shared by both pipelines:
      merged KS dirs → params → probe check → merged waveforms.

    merged_dir is treated as a complete, self-contained KS-style dataset: every
    numbered subfolder (0, 1, 2, ...) is a session, discovered directly on disk
    (see module docstring for why we don't consult the pre-merge UnitMatch.mat /
    UMparam.pickle). Good units per session come from that session's own
    cluster_bc_unitType.tsv; a unit whose TSV row still says GOOD but lost a
    within-session merge (RawSpikes file gone) is dropped by
    load_waveforms_for_good_units.

    Returns a dict with all session data, or None on failure.
    Both run functions receive this dict and work on an independent copy of param
    so that mutations in one pipeline do not affect the other.
    """
    session_names = sorted(
        (
            d
            for d in os.listdir(merged_dir)
            if d.isdigit() and os.path.isdir(os.path.join(merged_dir, d))
        ),
        key=int,
    )
    if not session_names:
        print(f"  ERROR: no numbered session folders found under {merged_dir}")
        return None

    ks_dirs = [os.path.join(merged_dir, d) for d in session_names]
    print(f"  {len(ks_dirs)} session(s): {session_names}")

    try:
        wave_paths, unit_label_paths, channel_pos = util.paths_from_KS(ks_dirs)
    except Exception as e:
        print(f"  ERROR in paths_from_KS: {e}")
        traceback.print_exc()
        return None

    good_units_per_session = build_good_units_from_labels(unit_label_paths)
    for i, d in enumerate(ks_dirs):
        print(f"    [{i}] {d}  ({len(good_units_per_session[i])} good units)")

    param = {"KS_dirs": ks_dirs}
    param = default_params.get_default_param(param=param)
    param = util.get_probe_geometry(channel_pos[0], param)

    # ── probe compatibility check ────────────────────────────────────────────
    # channel_pos entries are (nChan, 2) or (nChan, 3); when 3-col the first
    # column is a shank/depth offset so x is column index 1, otherwise 0.
    cp = channel_pos[0]
    x_col = cp[:, 1] if cp.shape[1] == 3 else cp[:, 0]
    actual_n_xchannelpos = int(len(np.unique(x_col)))
    unique_x = np.unique(x_col)
    x_gaps = np.diff(np.sort(unique_x))
    n_shanks = int(np.sum(x_gaps > 50)) + 1
    if actual_n_xchannelpos != param["n_xchannelpos"] * n_shanks:
        print(
            f"  SKIPPING: probe has {actual_n_xchannelpos} x-column position(s) across "
            f"{n_shanks} shank(s), expected a multiple of {param['n_xchannelpos']} "
            f"({param['n_xchannelpos'] * n_shanks}). "
            f"Set param['n_xchannelpos'] = {actual_n_xchannelpos // n_shanks} to process this probe type."
        )
        return None

    print("Loading merged waveforms …")
    try:
        waveform, session_id, session_switch, within_session, good_units, param = (
            load_waveforms_for_good_units(wave_paths, good_units_per_session, param)
        )
    except Exception as e:
        print(f"  ERROR loading waveforms: {e}")
        traceback.print_exc()
        return None

    param["good_units"] = good_units
    print(
        f"  {waveform.shape[0]} units across {param['n_sessions']} session(s) after merging"
    )

    return {
        "merged_dir": merged_dir,
        "channel_pos": channel_pos,
        "waveform": waveform,
        "session_id": session_id,
        "session_switch": session_switch,
        "within_session": within_session,
        "good_units": good_units,
        "param": param,
    }


# ── DeepUnitMatch pipeline ────────────────────────────────────────────────────


def run_deep_unit_match(sess):
    """Run the full DeepUnitMatch pipeline (default trained model) for one pre-loaded session."""
    merged_dir = sess["merged_dir"]
    save_dir = get_save_dir(merged_dir)

    print("Loading DeepUnitMatch model …")
    model = test.load_trained_model(device=DEVICE)

    run_deep_unit_match_core(sess, save_dir, model, label="DeepUnitMatch")


def run_deep_unit_match_core(sess, save_dir, model, label="DeepUnitMatch"):
    """
    Run the full DeepUnitMatch pipeline for one pre-loaded session, given an
    already-loaded model and an explicit output directory.

    Factored out of run_deep_unit_match() so alternative trained models (see
    run_deepunitmatch_batch_onMerged_extramodels.py) can reuse the exact same
    inference/matching/saving logic against a different checkpoint and save_dir.
    """
    merged_dir = sess["merged_dir"]
    print(f"\n--- {label}: {merged_dir}")

    tmp_path = os.path.join(save_dir, "tmp_waveforms")
    print(f"Save dir : {save_dir}")
    os.makedirs(save_dir, exist_ok=True)
    os.makedirs(tmp_path, exist_ok=True)

    # work on independent copies so UMPy (running after) sees clean state
    param = copy.deepcopy(sess["param"])
    waveform = sess["waveform"]
    session_id = sess["session_id"]
    session_switch = sess["session_switch"]
    good_units = sess["good_units"]
    channel_pos = sess["channel_pos"]
    within_session = sess["within_session"]

    # ── preprocess with DeepUnitMatch (get_snippets → HDF5) ─────────────────
    print("Preprocessing waveforms (get_snippets) …")
    unit_ids = np.concatenate(param["good_units"]).squeeze()
    try:
        _, _, kept_idx = param_fun.get_snippets(
            waveform,
            channel_pos,
            session_id,
            save_path=tmp_path,
            unit_ids=unit_ids,
            param=param,
        )
    except Exception as e:
        print(f"  ERROR in get_snippets: {e}")
        traceback.print_exc()
        return

    # re-sync arrays if any units were rejected by get_snippets
    if len(kept_idx) < len(waveform):
        waveform, session_id, session_switch, within_session, good_units, param = (
            util.filter_units_by_index(
                waveform, session_id, session_switch, good_units, kept_idx, param
            )
        )

    # ── neural-net inference ─────────────────────────────────────────────────
    print("Running DeepUnitMatch inference …")
    data_dir = os.path.join(tmp_path, "processed_waveforms")
    try:
        sim_matrix = test.inference(model, data_dir)
    except Exception as e:
        print(f"  ERROR in inference: {e}")
        traceback.print_exc()
        return

    # ── Naive Bayes matching ─────────────────────────────────────────────────
    print("Running Naive Bayes matching …")
    clus_info = {
        "good_units": param["good_units"],
        "session_switch": session_switch,
        "session_id": session_id,
        "original_ids": np.concatenate(param["good_units"]),
    }
    extracted_wave_properties = ov.extract_parameters(
        waveform, channel_pos, clus_info, param
    )
    sessions = np.unique(session_id)

    # ── pre-pass: collect DNN labels for all session pairs → shared drift correction ──
    # Build a global labels matrix from neural-net matches across every pair, then
    # call mf.drift_n_sessions — the same function UMPy uses — so both pipelines
    # share a single, consistent drift-correction mechanism.
    n_total = waveform.shape[0]
    labels_full = np.eye(n_total)
    pair_matches_cache = {}

    print("  Pre-pass: collecting neural-net labels for drift correction ...")
    for r1 in sessions:
        for r2 in sessions:
            if r1 >= r2:
                continue
            mask = np.isin(session_id, [r1, r2])
            sim_mat = sim_matrix[mask][:, mask]
            indices = np.where(mask)[0]
            n = int(np.sum(mask))

            df = helpers.create_dataframe(
                [param["good_units"][r1], param["good_units"][r2]],
                sim_mat,
                session_list=[r1, r2],
            )
            matches = test.get_matches(
                df, sim_mat, session_id[indices], data_dir, dist_thresh=param["max_dist"]
            )
            pair_matches_cache[(r1, r2)] = matches

            subsessionid = np.array(
                [r1] * len(param["good_units"][r1])
                + [r2] * len(param["good_units"][r2])
            )
            labels_pair = np.eye(n)
            for (recses1, recses2), group in matches.groupby(by=["RecSes1", "RecSes2"]):
                asmatrix = (
                    group["match"]
                    .values.reshape(
                        len(param["good_units"][recses1]),
                        len(param["good_units"][recses2]),
                    )
                    .astype(int)
                )
                labels_pair[
                    np.ix_(subsessionid == recses1, subsessionid == recses2)
                ] = asmatrix
            labels_full[np.ix_(indices, indices)] = labels_pair

    # Apply drift correction on full arrays — identical mechanism to UMPy.
    # sim_matrix serves as total_score so get_good_matches can de-duplicate pairs.
    avg_centroid = extracted_wave_properties["avg_centroid"].copy()
    avg_waveform_per_tp = extracted_wave_properties["avg_waveform_per_tp"].copy()
    _, avg_centroid, avg_waveform_per_tp = mf.drift_n_sessions(
        labels_full.astype(bool),
        session_switch,
        avg_centroid,
        avg_waveform_per_tp,
        sim_matrix,
        param,
    )

    # ── Bayes loop: use drift-corrected arrays ────────────────────────────────
    probs = np.zeros(sim_matrix.shape)
    distance_matrix = np.zeros(sim_matrix.shape)

    for r1 in sessions:
        for r2 in sessions:
            if r1 >= r2:
                continue

            mask = np.isin(session_id, [r1, r2])
            sim_mat = sim_matrix[mask][:, mask]
            n = int(np.sum(mask))
            indices = np.where(mask)[0]

            matches = pair_matches_cache[(r1, r2)]

            labels = np.eye(sim_mat.shape[0])
            subsessionid = np.array(
                [r1] * len(param["good_units"][r1])
                + [r2] * len(param["good_units"][r2])
            )
            for (recses1, recses2), group in matches.groupby(by=["RecSes1", "RecSes2"]):
                asmatrix = (
                    group["match"]
                    .values.reshape(
                        len(param["good_units"][recses1]),
                        len(param["good_units"][recses2]),
                    )
                    .astype(int)
                )
                labels[np.ix_(subsessionid == recses1, subsessionid == recses2)] = (
                    asmatrix
                )

            # use drift-corrected waveform (masked to this session pair)
            avg_waveform_per_tp_pair = avg_waveform_per_tp[:, mask, :, :]
            avg_waveform_per_tp_flip = mf.flip_dim(avg_waveform_per_tp_pair, param, n)
            euclid_dist = mf.get_Euclidean_dist(avg_waveform_per_tp_flip, param, n)
            centroid_dist, _ = mf.centroid_metrics(euclid_dist, param)

            scores_to_incl = {"similarity": sim_mat, "distance": centroid_dist}
            n_units = int(np.sqrt(len(matches)))

            # Adaptive match-count prior, mirroring UMPy's post-drift-correction
            # n_expected_matches estimate (overlord.extract_metric_scores,
            # is_first_pass=False) instead of a fixed "2 expected matches per
            # unit" heuristic. raw_dist is the same peak-timepoint,
            # min-over-flips physical distance (um) that centroid_metrics()
            # above rescales into centroid_dist -- reused here unscaled.
            waveidx_arr = np.asarray(param["waveidx"])
            peak_idx = int(np.flatnonzero(waveidx_arr == param["peak_loc"])[0])
            raw_dist = np.nanmin(euclid_dist[:, peak_idx, :, :], axis=1)
            within_session_pair = within_session[np.ix_(mask, mask)]
            include_these_pairs_idx = raw_dist < param["max_dist"]
            param_pair = dict(param, n_units=n_units)
            thrs_opt = mf.get_threshold(
                sim_mat, within_session_pair, raw_dist, param_pair, is_first_pass=False
            )
            n_expected_matches = int(
                np.sum((sim_mat > thrs_opt) & include_these_pairs_idx)
            )
            n_candidate_pairs = max(int(np.sum(include_these_pairs_idx)), 1)
            prior_match = 1 - (n_expected_matches / n_candidate_pairs)
            priors = np.array([prior_match, 1 - prior_match])

            parameter_kernels = bf.get_parameter_kernels(
                scores_to_incl, labels, np.unique(labels), param
            )
            predictors = np.stack(list(scores_to_incl.values()), axis=2)
            probability = bf.apply_naive_bayes(
                parameter_kernels, priors, predictors, param, np.unique(labels)
            )
            prob_matrix = probability[:, 1].reshape(n_units, n_units)

            probs[np.ix_(mask, mask)] = prob_matrix
            distance_matrix[np.ix_(mask, mask)] = centroid_dist

    # ── final matches ────────────────────────────────────────────────────────
    match_threshold = param.get("match_threshold", THRESH)
    final_matches = test.directional_filter_matrix(probs, session_id, match_threshold)
    n_matches = int(np.sum(final_matches)) // 2
    print(f"  {n_matches} matches found (threshold={match_threshold})")

    # ── assign unique IDs ────────────────────────────────────────────────────
    UIDs = aid.assign_unique_id(probs, param, clus_info)

    # ── performance metrics (AUC against functional scores) ──────────────────
    functional_scores = {}
    try:
        isicorr = test.ISI_correlations(param)
        auc_isi = test.AUC(final_matches, isicorr, session_id)
        print(f"AUC (ISI correlations):            {auc_isi:.3f}")
        functional_scores["ISI_correlations"] = isicorr

        isikl = test.ISI_KL_divergence(param)
        auc_isikl = test.AUC(final_matches, -isikl, session_id)
        print(f"AUC (ISI KL divergence):           {auc_isikl:.3f}")
        functional_scores["ISI_KL_divergence"] = isikl

        isiwass = test.ISI_wasserstein_distance(param)
        auc_isiwass = test.AUC(final_matches, -isiwass, session_id)
        print(f"AUC (ISI Wasserstein distance):    {auc_isiwass:.3f}")
        functional_scores["ISI_wasserstein_distance"] = isiwass

        refpopcorr = test.refpop_correlations(param, matches=final_matches)
        auc_refpop = test.AUC(final_matches, refpopcorr, session_id)
        print(f"AUC (ref. pop. correlation):       {auc_refpop:.3f}")
        functional_scores["refpop_correlations"] = refpopcorr

        frdiff = test.FR_diff(param)
        auc_fr = test.AUC(final_matches, -frdiff, session_id)
        print(f"AUC (firing rate difference):      {auc_fr:.3f}")
        functional_scores["FR_diff"] = frdiff

        cvdiff = test.ISI_CV_diff(param)
        auc_cv = test.AUC(final_matches, -cvdiff, session_id)
        print(f"AUC (ISI CV difference):           {auc_cv:.3f}")
        functional_scores["ISI_CV_diff"] = cvdiff

        try:
            natimcorr = test.natim_correlations(param, merged_architecture=True)
            auc_natim = test.AUC(final_matches, natimcorr, session_id)
            print(f"AUC (nat. image correlations):     {auc_natim:.3f}")
            functional_scores["natim_correlations"] = natimcorr
        except Exception:
            pass  # natim data may not be available for every session
    except Exception as e:
        print(f"  WARNING: functional score computation failed: {e}")
        functional_scores = {}

    # ── save ─────────────────────────────────────────────────────────────────
    su.save_to_output(
        save_dir,
        {"distance": distance_matrix},
        np.argwhere(final_matches),
        probs,
        extracted_wave_properties["avg_centroid"],
        extracted_wave_properties["avg_waveform"],
        extracted_wave_properties["avg_waveform_per_tp"],
        extracted_wave_properties["max_site"],
        distance_matrix,
        final_matches,
        clus_info,
        param,
        UIDs=UIDs,
        matches_curated=None,
        save_match_table=True,
        functional_scores=functional_scores if functional_scores else None,
    )
    ensure_matlab_compatible_output(save_dir)
    su.save_auc_summary(
        save_dir, test.auc_summary_from_functional_scores(functional_scores, final_matches, session_id)
    )

    # ── save diagnostic figures ───────────────────────────────────────────────
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))

    im = axes[0].imshow(sim_matrix, cmap="viridis", aspect="auto")
    axes[0].set_title("Similarity matrix")
    axes[0].set_xlabel("Unit")
    axes[0].set_ylabel("Unit")
    fig.colorbar(im, ax=axes[0])

    im = axes[1].imshow(probs, cmap="viridis", aspect="auto")
    axes[1].set_title("Match probability")
    axes[1].set_xlabel("Unit")
    axes[1].set_ylabel("Unit")
    fig.colorbar(im, ax=axes[1])

    im = axes[2].imshow(final_matches, cmap="viridis", aspect="auto")
    axes[2].set_title(f"Final matches (n={n_matches})")
    axes[2].set_xlabel("Unit")
    axes[2].set_ylabel("Unit")
    fig.colorbar(im, ax=axes[2])

    fig.tight_layout()
    fig.savefig(os.path.join(save_dir, "MatchingOverview.png"), dpi=150)
    plt.close(fig)

    if functional_scores:
        score_meta = {
            "ISI_correlations": ("ISI correlations", "viridis", None),
            "ISI_KL_divergence": ("ISI KL divergence", "magma", None),
            "ISI_wasserstein_distance": ("ISI Wasserstein distance", "magma", None),
            "refpop_correlations": ("Ref. pop. correlations", "viridis", None),
            "FR_diff": ("Firing rate difference", "magma", None),
            "ISI_CV_diff": ("ISI CV difference", "magma", None),
            "natim_correlations": ("Nat. image correlations", "viridis", None),
        }
        keys = [k for k in score_meta if k in functional_scores]
        fig, axes = plt.subplots(1, len(keys), figsize=(5 * len(keys), 5))
        if len(keys) == 1:
            axes = [axes]
        for ax, key in zip(axes, keys):
            title, cmap, _ = score_meta[key]
            im = ax.imshow(functional_scores[key], cmap=cmap, aspect="auto")
            ax.set_title(title)
            ax.set_xlabel("Unit")
            ax.set_ylabel("Unit")
            fig.colorbar(im, ax=ax)
        fig.tight_layout()
        fig.savefig(os.path.join(save_dir, "FunctionalScores.png"), dpi=150)
        plt.close(fig)

    print(f"  Results saved to: {save_dir}")


# ── UMPy pipeline ─────────────────────────────────────────────────────────────


def run_umpy(sess):
    """Run the full UMPy pipeline for one pre-loaded session.

    Unlike DeepUnitMatch, the Naive Bayes here uses the extracted waveform
    metric scores (amplitude, spatial decay, waveform similarity, etc.) rather
    than a neural-net similarity score.
    """
    merged_dir = sess["merged_dir"]
    save_dir = get_umpy_save_dir(merged_dir)
    run_umpy_core(sess, save_dir, label="UMPy")


def run_umpy_core(sess, save_dir, label="UMPy"):
    """
    Run the full UMPy pipeline for one pre-loaded session, given an explicit
    output directory.

    Factored out of run_umpy() so alternative parameter sweeps (see
    run_deepunitmatch_batch_onMerged_maxdist_sweep.py) can reuse the exact same
    extraction/matching/saving logic against a modified param dict (e.g. a
    different max_dist) and a different save_dir -- mirrors how
    run_deep_unit_match_core was factored out of run_deep_unit_match() for the
    extramodels script.
    """
    merged_dir = sess["merged_dir"]
    print(f"\n--- {label}: {merged_dir}")

    print(f"Save dir : {save_dir}")
    os.makedirs(save_dir, exist_ok=True)

    # work on independent copies so param mutations don't bleed between pipelines
    param = copy.deepcopy(sess["param"])
    waveform = sess["waveform"]
    session_id = sess["session_id"]
    session_switch = sess["session_switch"]
    within_session = sess["within_session"]
    channel_pos = sess["channel_pos"]

    clus_info = {
        "good_units": param["good_units"],
        "session_switch": session_switch,
        "session_id": session_id,
        "original_ids": np.concatenate(param["good_units"]),
    }

    # ── extract waveform properties ──────────────────────────────────────────
    print("Extracting waveform properties …")
    extracted_wave_properties = ov.extract_parameters(
        waveform, channel_pos, clus_info, param
    )

    # ── extract metric scores (these go directly into Naive Bayes) ───────────
    print("Extracting metric scores …")
    try:
        total_score, candidate_pairs, scores_to_include, predictors = (
            ov.extract_metric_scores(
                extracted_wave_properties,
                session_switch,
                within_session,
                param,
                niter=2,
            )
        )
    except Exception as e:
        print(f"  ERROR in extract_metric_scores: {e}")
        traceback.print_exc()
        return

    # ── Naive Bayes matching ─────────────────────────────────────────────────
    print("Running Naive Bayes matching …")
    prior_match = 1 - (param["n_expected_matches"] / param["n_units"] ** 2)
    priors = np.array([prior_match, 1 - prior_match])
    labels = candidate_pairs.astype(int)
    cond = np.unique(labels)
    parameter_kernels = bf.get_parameter_kernels(
        scores_to_include, labels, cond, param, add_one=1
    )
    probability = bf.apply_naive_bayes(
        parameter_kernels, priors, predictors, param, cond
    )
    output_prob_matrix = probability[:, 1].reshape(param["n_units"], param["n_units"])

    match_threshold = param.get("match_threshold", THRESH)
    output_threshold = np.zeros_like(output_prob_matrix)
    output_threshold[output_prob_matrix > match_threshold] = 1
    matches = np.argwhere(output_threshold == 1)
    n_matches = len(matches) // 2
    print(f"  {n_matches} matches found (threshold={match_threshold})")

    # ── assign unique IDs ────────────────────────────────────────────────────
    UIDs = aid.assign_unique_id(output_prob_matrix, param, clus_info)

    # ── performance metrics (AUC against functional scores) ──────────────────
    functional_scores = {}
    final_matches_bool = output_threshold.astype(bool)
    try:
        isicorr = test.ISI_correlations(param)
        auc_isi = test.AUC(final_matches_bool, isicorr, session_id)
        print(f"AUC (ISI correlations):            {auc_isi:.3f}")
        functional_scores["ISI_correlations"] = isicorr

        isikl = test.ISI_KL_divergence(param)
        auc_isikl = test.AUC(final_matches_bool, -isikl, session_id)
        print(f"AUC (ISI KL divergence):           {auc_isikl:.3f}")
        functional_scores["ISI_KL_divergence"] = isikl

        isiwass = test.ISI_wasserstein_distance(param)
        auc_isiwass = test.AUC(final_matches_bool, -isiwass, session_id)
        print(f"AUC (ISI Wasserstein distance):    {auc_isiwass:.3f}")
        functional_scores["ISI_wasserstein_distance"] = isiwass

        refpopcorr = test.refpop_correlations(param, matches=final_matches_bool)
        auc_refpop = test.AUC(final_matches_bool, refpopcorr, session_id)
        print(f"AUC (ref. pop. correlation):       {auc_refpop:.3f}")
        functional_scores["refpop_correlations"] = refpopcorr

        frdiff = test.FR_diff(param)
        auc_fr = test.AUC(final_matches_bool, -frdiff, session_id)
        print(f"AUC (firing rate difference):      {auc_fr:.3f}")
        functional_scores["FR_diff"] = frdiff

        cvdiff = test.ISI_CV_diff(param)
        auc_cv = test.AUC(final_matches_bool, -cvdiff, session_id)
        print(f"AUC (ISI CV difference):           {auc_cv:.3f}")
        functional_scores["ISI_CV_diff"] = cvdiff

        try:
            natimcorr = test.natim_correlations(param)
            auc_natim = test.AUC(final_matches_bool, natimcorr, session_id)
            print(f"AUC (nat. image correlations):     {auc_natim:.3f}")
            functional_scores["natim_correlations"] = natimcorr
        except Exception:
            pass
    except Exception as e:
        print(f"  WARNING: functional score computation failed: {e}")
        functional_scores = {}

    # ── save ─────────────────────────────────────────────────────────────────
    avg_centroid = extracted_wave_properties["avg_centroid"]
    avg_waveform = extracted_wave_properties["avg_waveform"]
    avg_waveform_per_tp = extracted_wave_properties["avg_waveform_per_tp"]
    max_site = extracted_wave_properties["max_site"]

    su.save_to_output(
        save_dir,
        scores_to_include,
        matches,
        output_prob_matrix,
        avg_centroid,
        avg_waveform,
        avg_waveform_per_tp,
        max_site,
        total_score,
        output_threshold,
        clus_info,
        param,
        UIDs=UIDs,
        matches_curated=None,
        save_match_table=True,
        functional_scores=functional_scores if functional_scores else None,
    )
    ensure_matlab_compatible_output(save_dir)
    su.save_auc_summary(
        save_dir,
        test.auc_summary_from_functional_scores(functional_scores, final_matches_bool, session_id),
    )

    # ── save diagnostic figures ───────────────────────────────────────────────
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))

    im = axes[0].imshow(total_score, cmap="viridis", aspect="auto")
    axes[0].set_title("Total score")
    axes[0].set_xlabel("Unit")
    axes[0].set_ylabel("Unit")
    fig.colorbar(im, ax=axes[0])

    im = axes[1].imshow(output_prob_matrix, cmap="viridis", aspect="auto")
    axes[1].set_title("Match probability")
    axes[1].set_xlabel("Unit")
    axes[1].set_ylabel("Unit")
    fig.colorbar(im, ax=axes[1])

    im = axes[2].imshow(output_threshold, cmap="viridis", aspect="auto")
    axes[2].set_title(f"Final matches (n={n_matches})")
    axes[2].set_xlabel("Unit")
    axes[2].set_ylabel("Unit")
    fig.colorbar(im, ax=axes[2])

    fig.tight_layout()
    fig.savefig(os.path.join(save_dir, "MatchingOverview.png"), dpi=150)
    plt.close(fig)

    if functional_scores:
        score_meta = {
            "ISI_correlations": ("ISI correlations", "viridis", None),
            "ISI_KL_divergence": ("ISI KL divergence", "magma", None),
            "ISI_wasserstein_distance": ("ISI Wasserstein distance", "magma", None),
            "refpop_correlations": ("Ref. pop. correlations", "viridis", None),
            "FR_diff": ("Firing rate difference", "magma", None),
            "ISI_CV_diff": ("ISI CV difference", "magma", None),
            "natim_correlations": ("Nat. image correlations", "viridis", None),
        }
        keys = [k for k in score_meta if k in functional_scores]
        fig, axes = plt.subplots(1, len(keys), figsize=(5 * len(keys), 5))
        if len(keys) == 1:
            axes = [axes]
        for ax, key in zip(axes, keys):
            title, cmap, _ = score_meta[key]
            im = ax.imshow(functional_scores[key], cmap=cmap, aspect="auto")
            ax.set_title(title)
            ax.set_xlabel("Unit")
            ax.set_ylabel("Unit")
            fig.colorbar(im, ax=ax)
        fig.tight_layout()
        fig.savefig(os.path.join(save_dir, "FunctionalScores.png"), dpi=150)
        plt.close(fig)

    print(f"  Results saved to: {save_dir}")


# ── entry point ───────────────────────────────────────────────────────────────


def parse_args():
    parser = argparse.ArgumentParser(
        description="Run DeepUnitMatch and UMPy on the merged dataset."
    )
    parser.add_argument(
        "--write-matlab-compat",
        action="store_true",
        help="Also write a MATLAB-compatible UnitMatch.mat from the Python outputs.",
    )
    return parser.parse_args()


def main():
    args = parse_args()
    global WRITE_MATLAB_COMPAT
    WRITE_MATLAB_COMPAT = args.write_matlab_compat

    print(f"Scanning for merged-data groups under:\n  {BASE_INPUT}\n")

    groups = find_merged_groups()

    if not groups:
        print("No merged-data groups found.")
        return

    print(f"Found {len(groups)} group(s).\n")

    for i, merged_dir in enumerate(groups):
        print(f"\n[{i + 1}/{len(groups)}] {merged_dir}")

        run_deep = not results_exist(merged_dir)
        run_ump = not umpy_results_exist(merged_dir)

        if not run_deep and not run_ump:
            print("  Skipping both pipelines (results exist and are fresh).")
            continue

        lock_path = get_group_lock_path(merged_dir)
        with batch_lock.try_lock(lock_path) as acquired:
            if not acquired:
                print(f"  Skipping (already being processed by another run): {lock_path}")
                continue

            # re-check now that we hold the lock: another machine may have
            # finished this group while we were scanning/waiting for the lock
            run_deep = not results_exist(merged_dir)
            run_ump = not umpy_results_exist(merged_dir)
            if not run_deep and not run_ump:
                print("  Skipping both pipelines (completed by another run).")
                continue

            sess = _prepare_session(merged_dir)
            if sess is None:
                continue

            if run_deep:
                try:
                    run_deep_unit_match(sess)
                except Exception as e:
                    print(f"  DeepUnitMatch FAILED: {e}")
                    traceback.print_exc()
            else:
                print(
                    f"  Skipping DeepUnitMatch (results exist and are fresh): {get_save_dir(merged_dir)}"
                )

            if run_ump:
                try:
                    run_umpy(sess)
                except Exception as e:
                    print(f"  UMPy FAILED: {e}")
                    traceback.print_exc()
            else:
                print(
                    f"  Skipping UMPy (results exist and are fresh): {get_umpy_save_dir(merged_dir)}"
                )

    print("\nAll done.")


if __name__ == "__main__":
    main()
