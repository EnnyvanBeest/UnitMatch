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
import datetime
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

import batch_lock
import pipeline_config as cfg
import pipeline_log as plog
import run_deepunitmatch_batch_onMerged as onm  # shared matching pipeline
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
# Paths come from pipeline_config.py (step 1: raw data -> non-merged output).
BASE_INPUT = cfg.RAW_KS_BASE
BASE_OUTPUT = cfg.UNMERGED_OUTPUT

DEVICE = "cuda" if test.torch.cuda.is_available() else "cpu"
print(f"Device: {DEVICE}")
THRESH = 0.5
# See batch_lock.sentinel_is_fresh() / run_deepunitmatch_batch_onMerged.py's
# REDO_FROM_DATE for what this does: a group is skipped once its
# MatchingOverview.png exists and is at least this new. None falls back to
# plain "skip if present"; a far-future date reproduces old REDO=True.
REDO_FROM_DATE = cfg.REDO_FROM_DATE
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
    """Return True when the DeepUnitMatch sentinel output file is present and fresh (see REDO_FROM_DATE)."""
    sentinel = os.path.join(get_save_dir(mat_path), "MatchingOverview.png")
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


def umpy_results_exist(mat_path):
    """Return True when the UMPy sentinel output file is present and fresh (see REDO_FROM_DATE)."""
    sentinel = os.path.join(get_umpy_save_dir(mat_path), "MatchingOverview.png")
    return batch_lock.sentinel_is_fresh(sentinel, REDO_FROM_DATE)


def get_group_lock_path(mat_path):
    """
    Lock file marking 'a run is currently processing this group' (DeepUnitMatch
    + UMPy together), so multiple machines pointed at the same BASE_INPUT/
    BASE_OUTPUT can split work across mat files without double-processing one.
    See batch_lock.py. Mirrors run_deepunitmatch_batch_onMerged.py's lock --
    this script previously had no cross-machine lock at all (only the
    REDO_FROM_DATE freshness check above), so running it from multiple
    machines at once would have them all discover the same mat_files and
    duplicate/race on the same output.
    """
    subfolder = os.path.relpath(os.path.dirname(mat_path), BASE_INPUT)
    return os.path.join(BASE_OUTPUT, os.path.dirname(subfolder), ".processing.lock")


# ── shared session loader ─────────────────────────────────────────────────────


LOG_STAGE = "step1"
SENTINEL = "MatchingOverview.png"


def group_key(mat_path):
    """'mouse/probe/location' of a raw UnitMatch.mat (<loc>/UnitMatch/UnitMatch.mat)."""
    return os.path.relpath(os.path.dirname(os.path.dirname(mat_path)), BASE_INPUT).replace(os.sep, "/")


def _prep_failed(mat_path, message, status="failed", tb=None):
    """Print + log why a group could not be prepared; returns None for the caller to return."""
    print(f"  {'SKIPPING' if status == 'skipped' else 'ERROR'}: {message}")
    plog.log_event(LOG_STAGE, group_key(mat_path), "prepare_session", status, message, tb)
    return None


def _prepare_session(mat_path):
    """
    Load one raw location (see _prepare_session_impl); every failure is
    logged and gives None. Units the DNN preprocessing rejects are removed for
    both methods (onm.restrict_to_snippet_units), as on the merged data.
    """
    try:
        sess = _prepare_session_impl(mat_path)
        if sess is None:
            return None
        sess.update(
            group=group_key(mat_path),  # for the shared run functions' logging
            log_stage=LOG_STAGE,
            merged_architecture=False,  # natural-image trial files: non-merged layout
        )
        return onm.restrict_to_snippet_units(sess)
    except Exception as e:
        return _prep_failed(mat_path, f"{type(e).__name__}: {e}", tb=traceback.format_exc())


def _prepare_session_impl(mat_path):
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
        traceback.print_exc()
        return _prep_failed(mat_path, f"loading UnitMatch.mat: {type(e).__name__}: {e}", tb=traceback.format_exc())

    print(f"  {len(ks_dirs)} session(s):")
    for i, d in enumerate(ks_dirs):
        n_good = int(((recsesAll == (i + 1)) & good_id).sum())
        print(f"    [{i}] {d}  ({n_good} good units)")

    try:
        wave_paths, _, channel_pos = util.paths_from_KS(ks_dirs)
    except Exception as e:
        traceback.print_exc()
        return _prep_failed(mat_path, f"paths_from_KS: {type(e).__name__}: {e}", tb=traceback.format_exc())

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
        return _prep_failed(
            mat_path,
            f"probe has {actual_n_xchannelpos} x-column position(s) across "
            f"{n_shanks} shank(s), expected a multiple of {param['n_xchannelpos']} "
            f"({param['n_xchannelpos'] * n_shanks}). "
            f"Set param['n_xchannelpos'] = {actual_n_xchannelpos // n_shanks} to process this probe type.",
            status="skipped",
        )

    good_units_per_session = build_good_units_per_session(
        ks_dirs, orig_clus_id, recsesAll, good_id
    )

    print("Loading waveforms …")
    try:
        waveform, session_id, session_switch, within_session, good_units, param = (
            load_waveforms_for_good_units(wave_paths, good_units_per_session, param)
        )
    except Exception as e:
        traceback.print_exc()
        return _prep_failed(mat_path, f"loading waveforms: {type(e).__name__}: {e}", tb=traceback.format_exc())

    param["good_units"] = good_units
    print(f"  {waveform.shape[0]} units across {param['n_sessions']} session(s)")
    if param["n_sessions"] < 2:
        plog.log_event(LOG_STAGE, group_key(mat_path), "prepare_session", "done",
                       f"single session ({waveform.shape[0]} units): no across-session pairs")

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


# ── matching: the pipeline shared with the merged-data runs ──────────────────
# DeepUnitMatch and UMPy run exactly as on the merged data (run_dum_core /
# run_umpy_core in run_deepunitmatch_batch_onMerged.py: one matching
# pipeline for both methods, unique-ID match sets, logging). The session
# dict says where the natural-image trial files are (non-merged layout) and
# how to log (stage "step1", group "mouse/probe/location").


def run_deep_unit_match(sess, model):
    """DeepUnitMatch on one non-merged location (its unique IDs define the merges in step 2)."""
    onm.run_dum_core(sess, get_save_dir(sess["mat_path"]), model, label="DeepUnitMatch")


def run_umpy(sess):
    """UMPy on one non-merged location."""
    onm.run_umpy_core(sess, get_umpy_save_dir(sess["mat_path"]), label="UMPy")


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
    onm.WRITE_MATLAB_COMPAT = args.write_matlab_compat
    print("Loading DeepUnitMatch model …")
    model = test.load_trained_model(device=DEVICE)

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

        run_deep = not results_exist(mat_path)
        run_ump = not umpy_results_exist(mat_path)

        if not run_deep and not run_ump:
            print("  Skipping both pipelines (results exist and are fresh).")
            continue

        lock_path = get_group_lock_path(mat_path)
        with batch_lock.try_lock(lock_path) as acquired:
            if not acquired:
                print(f"  Skipping (already being processed by another run): {lock_path}")
                continue

            # re-check now that we hold the lock: another machine may have
            # finished this group while we were scanning/waiting for the lock
            run_deep = not results_exist(mat_path)
            run_ump = not umpy_results_exist(mat_path)
            if not run_deep and not run_ump:
                print("  Skipping both pipelines (completed by another run).")
                continue

            sess = _prepare_session(mat_path)
            if sess is None:
                continue

            if run_deep:
                try:
                    run_deep_unit_match(sess, model)
                except Exception as e:
                    print(f"  DeepUnitMatch FAILED: {e}")
                    traceback.print_exc()
            else:
                print(
                    f"  Skipping DeepUnitMatch (results exist and are fresh): {get_save_dir(mat_path)}"
                )

            if run_ump:
                try:
                    run_umpy(sess)
                except Exception as e:
                    print(f"  UMPy FAILED: {e}")
                    traceback.print_exc()
            else:
                print(
                    f"  Skipping UMPy (results exist and are fresh): {get_umpy_save_dir(mat_path)}"
                )

    print("\nAll done.")


if __name__ == "__main__":
    main()
