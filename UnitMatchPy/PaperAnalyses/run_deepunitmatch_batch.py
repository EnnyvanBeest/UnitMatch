# Batch wrapper: runs DeepUnitMatch and UMPy on every UnitMatch.mat found under
#    \\znas.cortexlab.net\Lab\Share\UNITMATCHTABLES_ENNY_CELIAN_JULIE\FullAnimal_KSChanMap

# KS directories are read from UMparam.KSDir; good units are taken from
# UniqueIDConversion (OriginalClusID[GoodID], indexed per session via recsesAll).

# Results are mirrored to:
#    \\znas.cortexlab.net\Lab\Share\UNITMATCHTABLES_ENNY_CELIAN_JULIE\DeepUM_NatMeth2026V2
# with the same subfolder structure, split into DeepUnitMatch/ and UMPy/ subfolders.
# Waveforms are loaded once per mat file and shared between both pipelines.

import os
import sys
import copy
import traceback
import numpy as np
import scipy.io
import h5py
import matplotlib

matplotlib.use("Agg")  # non-interactive backend for batch runs
import matplotlib.pyplot as plt

# ── project paths ───────────────────────────────────────────────────────────
_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _HERE)
sys.path.insert(0, os.path.dirname(_HERE))
sys.path.insert(0, os.path.join(_HERE, "DeepUnitMatch"))

import argparse

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
BASE_INPUT = r"\\znas.cortexlab.net\Lab\Share\UNITMATCHTABLES_ENNY_CELIAN_JULIE\FullAnimal_KSChanMap"
BASE_OUTPUT = r"\\znas.cortexlab.net\Lab\Share\UNITMATCHTABLES_ENNY_CELIAN_JULIE\DeepUM_NatMeth2026V2"

DEVICE = "cuda" if test.torch.cuda.is_available() else "cpu"
print(f"Device: {DEVICE}")
THRESH = 0.5
REDO = True  # if True, rerun even when output already exists
WRITE_MATLAB_COMPAT = False


# ── MATLAB file loading ───────────────────────────────────────────────────────


def _decode_hdf5_str(f, ref_or_ds):
    """Decode a MATLAB char array stored as uint16 in an HDF5 file."""
    if isinstance(ref_or_ds, h5py.Reference):
        obj = f[ref_or_ds]
    else:
        obj = ref_or_ds
    chars = obj[()].flatten()
    return "".join(chr(int(c)) for c in chars)


def _hdf5_cellstr(f, dataset):
    """Read a MATLAB cell array of strings from an h5py dataset."""
    data = dataset[()]
    flat = data.flatten()
    return [_decode_hdf5_str(f, ref) for ref in flat]


def _load_mat_scipy(mat_path):
    """Load via scipy (MATLAB < v7.3). Returns (ks_dirs, orig_clus_id, recsesAll, good_id)."""
    mat = scipy.io.loadmat(mat_path, simplify_cells=True)
    ump = mat["UMparam"]
    uid = mat["UniqueIDConversion"]

    ks_dirs = ump["KSDir"]
    if isinstance(ks_dirs, str):
        ks_dirs = [ks_dirs]
    elif isinstance(ks_dirs, np.ndarray):
        ks_dirs = [str(s) for s in ks_dirs.flatten()]

    orig_clus_id = np.array(uid["OriginalClusID"]).flatten()
    recsesAll = np.array(uid["recsesAll"]).flatten()
    good_id = np.array(uid["GoodID"]).flatten().astype(bool)

    return ks_dirs, orig_clus_id, recsesAll, good_id


def _load_mat_hdf5(mat_path):
    """Load via h5py (MATLAB v7.3 HDF5). Returns (ks_dirs, orig_clus_id, recsesAll, good_id)."""
    with h5py.File(mat_path, "r") as f:
        ks_dirs = _hdf5_cellstr(f, f["UMparam"]["KSDir"])

        uid = f["UniqueIDConversion"]
        orig_clus_id = uid["OriginalClusID"][()].flatten()
        recsesAll = uid["recsesAll"][()].flatten()
        good_id = uid["GoodID"][()].flatten().astype(bool)

    return ks_dirs, orig_clus_id, recsesAll, good_id


def load_unitmatchemat(mat_path):
    """
    Load UnitMatch.mat and return:
        ks_dirs       – list of KS directory paths (one per session)
        orig_clus_id  – cluster IDs for every neuron
        recsesAll     – session index (1-based, MATLAB convention) for every neuron
        good_id       – boolean mask of good neurons
    """
    try:
        return _load_mat_scipy(mat_path)
    except Exception:
        return _load_mat_hdf5(mat_path)


# ── good-unit helpers ─────────────────────────────────────────────────────────


def build_good_units_per_session(ks_dirs, orig_clus_id, recsesAll, good_id):
    """
    Returns a list (one entry per session) of (N, 1) float arrays of cluster IDs,
    matching the shape expected by UnitMatchPy internals.
    recsesAll is 1-indexed (MATLAB convention).
    """
    good_units = []
    for i in range(len(ks_dirs)):
        mask = (recsesAll == (i + 1)) & good_id
        ids = orig_clus_id[mask].astype(float)
        good_units.append(ids.reshape(-1, 1))
    return good_units


def load_waveforms_for_good_units(wave_paths, good_units_per_session, param):
    """
    Directly loads RawSpikes waveforms for the specified good units, bypassing
    the TSV-based unit-label files used by util.load_good_waveforms.

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
            first_id = int(g_units[0])
            p_first = os.path.join(wave_path, f"Unit{first_id}_RawSpikes.npy")
            ref = np.load(p_first)  # (T, C, spikes) or similar
            buf = np.zeros((len(g_units), *ref.shape), dtype=ref.dtype)

            kept_ids = []
            kept_idx = []
            for j, uid in enumerate(g_units):
                p = os.path.join(wave_path, f"Unit{int(uid)}_RawSpikes.npy")
                if os.path.exists(p):
                    buf[j] = np.load(p)
                    kept_ids.append(uid)
                    kept_idx.append(j)
                else:
                    print(f"  Warning: missing {p}")

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


def get_save_dir(mat_path):
    """Return the DeepUnitMatch output directory for a given UnitMatch.mat path."""
    subfolder = os.path.relpath(os.path.dirname(mat_path), BASE_INPUT)
    subfolder = os.path.join(os.path.dirname(subfolder), "DeepUnitMatch")
    return os.path.join(BASE_OUTPUT, subfolder)


def get_umpy_save_dir(mat_path):
    """Return the UMPy output directory for a given UnitMatch.mat path."""
    subfolder = os.path.relpath(os.path.dirname(mat_path), BASE_INPUT)
    subfolder = os.path.join(os.path.dirname(subfolder), "UMPy")
    return os.path.join(BASE_OUTPUT, subfolder)


def results_exist(mat_path):
    """Return True when the DeepUnitMatch sentinel output file is present."""
    sentinel = os.path.join(get_save_dir(mat_path), "MatchingOverview.png")
    return os.path.isfile(sentinel)


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


def umpy_results_exist(mat_path):
    """Return True when the UMPy sentinel output file is present."""
    sentinel = os.path.join(get_umpy_save_dir(mat_path), "MatchingOverview.png")
    return os.path.isfile(sentinel)


# ── shared session loader ─────────────────────────────────────────────────────


def _prepare_session(mat_path):
    """
    Load and validate everything shared by both pipelines:
      mat → ks_dirs → params → probe check → waveforms.

    Returns a dict with all session data, or None on failure.
    Both run functions receive this dict and work on an independent copy of param
    so that mutations in one pipeline do not affect the other.
    """
    print("Loading UnitMatch.mat …")
    try:
        ks_dirs, orig_clus_id, recsesAll, good_id = load_unitmatchemat(mat_path)
    except Exception as e:
        print(f"  ERROR loading mat: {e}")
        traceback.print_exc()
        return None

    print(f"  {len(ks_dirs)} session(s):")
    for i, d in enumerate(ks_dirs):
        n_good = int(((recsesAll == (i + 1)) & good_id).sum())
        print(f"    [{i}] {d}  ({n_good} good units)")

    try:
        wave_paths, _, channel_pos = util.paths_from_KS(ks_dirs)
    except Exception as e:
        print(f"  ERROR in paths_from_KS: {e}")
        traceback.print_exc()
        return None

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

    good_units_per_session = build_good_units_per_session(
        ks_dirs, orig_clus_id, recsesAll, good_id
    )

    print("Loading waveforms …")
    try:
        waveform, session_id, session_switch, within_session, good_units, param = (
            load_waveforms_for_good_units(wave_paths, good_units_per_session, param)
        )
    except Exception as e:
        print(f"  ERROR loading waveforms: {e}")
        traceback.print_exc()
        return None

    param["good_units"] = good_units
    print(f"  {waveform.shape[0]} units across {param['n_sessions']} session(s)")

    return {
        "mat_path": mat_path,
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
    """Run the full DeepUnitMatch pipeline for one pre-loaded session."""
    mat_path = sess["mat_path"]
    print(f"\n--- DeepUnitMatch: {mat_path}")

    save_dir = get_save_dir(mat_path)
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

    # ── load model ───────────────────────────────────────────────────────────
    print("Loading DeepUnitMatch model …")
    model = test.load_trained_model(device=DEVICE)

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
        waveform, session_id, session_switch, _, good_units, param = (
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
                df, sim_mat, session_id[indices], data_dir, dist_thresh=50
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
            priors = np.array([1 - 2 / n_units, 2 / n_units])
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
    final_matches = test.directional_filter(probs, session_id, THRESH)
    n_matches = int(np.sum(final_matches)) // 2
    print(f"  {n_matches} matches found (threshold={THRESH})")

    # ── assign unique IDs ────────────────────────────────────────────────────
    UIDs = aid.assign_unique_id(probs, param, clus_info)

    # ── performance metrics (AUC against functional scores) ──────────────────
    functional_scores = {}
    try:
        isicorr = test.ISI_correlations(param)
        auc_isi = test.AUC(final_matches, isicorr, session_id)
        print(f"AUC (ISI correlations):            {auc_isi:.3f}")
        functional_scores["ISI_correlations"] = isicorr

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
            natimcorr = test.natim_correlations(param)
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
    mat_path = sess["mat_path"]
    print(f"\n--- UMPy: {mat_path}")

    save_dir = get_umpy_save_dir(mat_path)
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
        description="Run DeepUnitMatch and UMPy on UnitMatch.mat inputs."
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

    print(f"Scanning for UnitMatch.mat files under:\n  {BASE_INPUT}\n")

    mat_files = []
    for root, _, files in os.walk(BASE_INPUT):
        if "UnitMatch.mat" in files and os.path.basename(root) == "UnitMatch":
            mat_files.append(os.path.join(root, "UnitMatch.mat"))

    if not mat_files:
        print("No UnitMatch.mat files found.")
        return

    print(f"Found {len(mat_files)} file(s).\n")

    for i, mat_path in enumerate(mat_files):
        print(f"\n[{i + 1}/{len(mat_files)}] {mat_path}")

        run_deep = not results_exist(mat_path) or REDO
        run_ump = not umpy_results_exist(mat_path) or REDO

        if not run_deep and not run_ump:
            print("  Skipping both pipelines (results exist, REDO=False).")
            continue

        sess = _prepare_session(mat_path)
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
                f"  Skipping DeepUnitMatch (results exist, REDO=False): {get_save_dir(mat_path)}"
            )

        if run_ump:
            try:
                run_umpy(sess)
            except Exception as e:
                print(f"  UMPy FAILED: {e}")
                traceback.print_exc()
        else:
            print(
                f"  Skipping UMPy (results exist, REDO=False): {get_umpy_save_dir(mat_path)}"
            )

    print("\nAll done.")


if __name__ == "__main__":
    main()
