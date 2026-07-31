# Batch wrapper: runs pyDANT (Density-based Across-day Neuron Tracking,
# https://github.com/jiumao2/pyDANT) on the same merged-dataset groups/
# sessions that run_deepunitmatch_batch_onMerged.py already runs UMPy/
# DeepUnitMatch on, so DANT can be compared to both (and to EMD) using the
# exact same functional-score AUC framework (test.AUC over a boolean
# final_matches matrix).
#
# Unlike EMD (MATLAB, pairwise day-by-day), pyDANT is a pure Python package
# that clusters units from an arbitrary number of sessions in one call, so
# this is a single-phase script -- no stage/external-tool/aggregate split
# needed the way EMD required.
#
# ── with vs without functional (activity-derived) data ──────────────────────
# pyDANT's own default config uses "AutoCorr" (autocorrelogram) alongside
# "Waveform" in both its motion-estimation and final-clustering feature
# lists -- i.e. it uses activity-derived data as *matching* input by default.
# UnitMatch's own philosophy deliberately avoids that, specifically so
# functional scores (ISI correlations, refpop correlations, etc.) stay
# independent evidence for post-hoc AUC validation rather than partially
# circular with what the matcher itself already used. So this script runs
# BOTH variants per group, sharing the same input arrays (only the feature
# selection differs):
#   DANT/               -- pyDANT's own recommended default (Waveform + AutoCorr)
#   DANT_no_functional/  -- Waveform only, matching UnitMatch's own philosophy
# Both are equally valid pyDANT configurations (this is not a hack -- pyDANT's
# settings.json schema already supports listing just ["Waveform"]); the two
# output folders let you see directly how much of DANT's tracking performance
# comes from activity-derived features vs waveform/anatomy alone.
#
# ── pyDANT's data contract (see UnitMatchPy/../pyDANT research notes) ───────
# runDANT(user_settings) expects a *flat* folder (user_settings["path_to_data"])
# containing, for ALL sessions pre-concatenated together:
#   waveform_all.npy       (n_unit, n_channel, n_sample) -- one MEAN waveform
#                           per unit (no cross-validation halves, unlike
#                           UnitMatch's own convention -- we average over the
#                           cv axis below)
#   session_index.npy      (n_unit,) int, 1-INDEXED, must start at 1 with no
#                           gaps (we use session_id + 1)
#   channel_locations.npy  (n_channel, 2) = (x, y) probe geometry, shared
#                           across all sessions in the group
#   spike_times/Unit{k}.npy  one file per unit, k = 0-indexed row in
#                           waveform_all, spike times in MILLISECONDS (NOT
#                           seconds, NOT samples -- verified directly against
#                           Preprocess.py/computeAutoCorr, which compare
#                           spike_times to the autoCorr/ISI window+binwidth
#                           settings, both in ms, with no internal unit
#                           conversion)
# There is no good-unit filtering inside pyDANT at all -- the caller is
# expected to have already restricted to good units, which _prepare_session()
# already does for us. This input is identical for both variants above, so
# it's built once per group and shared.
#
# Output written to user_settings["output_folder"] includes ClusterMatrix.npy
# -- an (n_unit, n_unit) BOOLEAN adjacency matrix (same row order as
# waveform_all.npy/session_index.npy, i.e. the same order we control) that is
# directly usable as our final_matches, no translation needed (unlike EMD,
# which required mapping MATLAB positions back to cluster IDs).
#
# channel_shanks.npy (optional; only required for pyDANT's separate
# runDANTMultiShank entry point) is intentionally omitted here -- we run the
# single-shank runDANT(), treating the whole probe as one shank. This is a
# scope simplification, not a pyDANT requirement; runDANTMultiShank could be
# wired in later for genuinely multi-shank probes if needed.
#
# peth.npy (optional; only used if "PETH" is listed in motionEstimation/
# clustering features) is also omitted -- we don't have trial-aligned PETH
# data readily prepared for every group, so "PETH" is dropped from both
# variants' feature lists rather than left in pointing at data we don't
# provide (pyDANT's own default motionEstimation feature list includes PETH;
# we drop it from the "DANT" variant here for that reason, not because it's
# itself the "no functional data" variant -- AutoCorr is still functional and
# still included in "DANT").
#
# IMPORTANT Windows/joblib gotcha (confirmed directly in pyDANT's source: no
# module there guards its own top-level code with `if __name__ == "__main__":`,
# and Preprocess.py/ComputeWaveformFeatures.py/MotionEstimation.py all use
# joblib.Parallel with the process-based "loky" backend by default): the
# entire batch loop below MUST run only under this script's own
# `if __name__ == "__main__":` guard, or worker processes spawned on Windows
# will re-import and re-execute this file's top-level code, recursively
# spawning more workers. main() is only ever invoked that way below.

import os
import sys
import datetime
import traceback

import numpy as np
import matplotlib

matplotlib.use("Agg")  # non-interactive backend for batch runs
import matplotlib.pyplot as plt

# Upstream compatibility gap: pyDANT's IterativeClustering.py calls np.atanh,
# which numpy only added in 2.0.0 (as an Array-API-standard-compatible alias
# for the older arctanh) -- so this crashes with AttributeError on any numpy
# < 2.0 (e.g. our own 1.26.4), even though it works fine for anyone on a
# recent numpy. Monkeypatched here rather than hand-edited in site-packages,
# or pinning numpy>=2.0, so the fix travels with this script to every machine
# regardless of which numpy version is installed there.
if not hasattr(np, "atanh"):
    np.atanh = np.arctanh

# Upstream bug: pyDANT's ComputeWaveformFeatures.py prints a micro sign (μm)
# unconditionally after every successful motion-estimation iteration. Windows'
# legacy console codepage (cp1252) can't encode it and crashes with
# UnicodeEncodeError -- confirmed directly (reproduced end-to-end against real
# data: the run completes fine with this reconfigured, crashes without it).
# Reconfiguring stdout/stderr to UTF-8 is the standard fix for this class of
# Windows-console issue and has no effect on Linux/Mac, where stdout is
# already UTF-8.
if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8")
    sys.stderr.reconfigure(encoding="utf-8")

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _HERE)
sys.path.insert(0, os.path.dirname(_HERE))
sys.path.insert(0, os.path.join(_HERE, "DeepUnitMatch"))

import batch_lock
import UnitMatchPy.overlord as ov
import UnitMatchPy.save_utils as su
from DeepUnitMatch.testing import test

from run_deepunitmatch_batch_onMerged import (
    BASE_INPUT,
    BASE_OUTPUT,
    find_merged_groups,
    _prepare_session,
)

from pyDANT import runDANT

FS = 30000.0  # Neuropixels sample rate -- same convention used throughout test.py

# n_jobs for pyDANT's own internal joblib parallelism (locations/autocorr/ISI
# preprocessing, motion-corrected-waveform computation, motion bootstrap CI).
# -1 uses all cores, matching pyDANT's own settings.json default -- fine here
# since groups (and variants within a group) are processed one at a time,
# same "one full-resource job at a time" pattern the rest of this pipeline uses.
N_JOBS = -1

# See batch_lock.sentinel_is_fresh() / run_deepunitmatch_batch_onMerged.py's
# REDO_FROM_DATE for what this does.
REDO_FROM_DATE = datetime.datetime(2026, 7, 30, 16, 0, 0)

# Two variants swept per group -- see module docstring for the "with vs
# without functional data" rationale. Both are otherwise identical (same
# motion-estimation loop structure, same clustering settings) except for
# which features are allowed in.
DANT_VARIANTS = {
    "DANT": {
        # A single combined feature set from the very first iteration, not a
        # staged "AutoCorr-only warmup then add Waveform" (pyDANT's own
        # default motionEstimation.features is staged that way, e.g.
        # [["AutoCorr","PETH"], ["Waveform","AutoCorr","PETH"]] -- we
        # originally mirrored that pattern). Confirmed via direct repro
        # against real data that the staged version can produce a
        # degenerate first-iteration clustering (AutoCorr alone finding
        # "matches" that turn out to be almost entirely same-session, i.e.
        # useless for drift estimation) which crashes pyDANT's own
        # computeMotion with zero cross-session pairs to work with, for
        # datasets where AutoCorr alone isn't very discriminative. Using
        # Waveform+AutoCorr together from iteration 0 avoids that
        # degenerate state entirely (confirmed: 0 -> ~12,000 cross-session
        # candidate pairs on the same real dataset) without changing what
        # the final clustering step (clustering_features, below) uses.
        "motionEstimation_features": [["Waveform", "AutoCorr"]],
        "clustering_features": ["Waveform", "AutoCorr"],
    },
    "DANT_no_functional": {
        "motionEstimation_features": [["Waveform"]],
        "clustering_features": ["Waveform"],
    },
}


# ── path helpers ─────────────────────────────────────────────────────────────


def get_dant_shared_dir(merged_dir):
    """Shared input arrays (identical across variants) live here, sibling to each variant's own output folder."""
    subfolder = os.path.relpath(os.path.dirname(merged_dir), BASE_INPUT)
    return os.path.join(BASE_OUTPUT, subfolder, "DANT_shared")


def get_dant_shared_input_dir(merged_dir):
    return os.path.join(get_dant_shared_dir(merged_dir), "_dant_input")


def get_dant_dir(merged_dir, variant):
    """Return the output directory for one DANT variant -- a sibling of DeepUnitMatch/UMPy/EMD."""
    subfolder = os.path.relpath(os.path.dirname(merged_dir), BASE_INPUT)
    return os.path.join(BASE_OUTPUT, subfolder, variant)


def get_dant_raw_output_dir(merged_dir, variant):
    return os.path.join(get_dant_dir(merged_dir, variant), "_dant_raw_output")


def dant_results_exist(merged_dir, variant):
    """Return True when this variant's sentinel output file is present and fresh (see REDO_FROM_DATE)."""
    sentinel = os.path.join(get_dant_dir(merged_dir, variant), "MatchingOverview.png")
    return batch_lock.sentinel_is_fresh(sentinel, REDO_FROM_DATE)


def get_dant_lock_path(merged_dir):
    """
    One lock per group (covering the shared input build + both variants), so
    multiple machines pointed at the same BASE_INPUT/BASE_OUTPUT can split
    work across groups without double-processing one. See batch_lock.py.
    """
    return os.path.join(get_dant_shared_dir(merged_dir), ".processing.lock")


# ── build pyDANT's input arrays from _prepare_session's output ──────────────


def build_dant_input(sess, input_dir):
    """
    Write waveform_all.npy / session_index.npy / channel_locations.npy /
    spike_times/Unit{k}.npy into input_dir, in exactly the layout
    pyDANT.Preprocess.preprocess() expects. Row/unit ordering matches
    sess["session_id"]/sess["good_units"] exactly (global index k = position
    in np.concatenate(good_units)), so pyDANT's output arrays -- which
    preserve that same row order -- can be read back directly against our
    own session_id/good_units without any remapping. Shared across variants
    (identical regardless of which features a variant uses).
    """
    os.makedirs(input_dir, exist_ok=True)

    param = sess["param"]
    waveform = sess["waveform"]  # (n_units, spike_width, n_channels, cv=2)
    session_id = sess["session_id"]
    channel_pos = sess["channel_pos"]  # list of (n_channels, 2 or 3), one per session
    good_units = sess["good_units"]
    ks_dirs = param["KS_dirs"]
    n_units = waveform.shape[0]

    # waveform_all.npy: average over cv -> (n_units, spike_width, n_channels)
    # -> pyDANT wants (n_units, n_channels, n_samples)
    waveform_mean = waveform.mean(axis=3)
    waveform_all = np.transpose(waveform_mean, (0, 2, 1)).astype(np.float64)
    np.save(os.path.join(input_dir, "waveform_all.npy"), waveform_all)

    # session_index.npy: pyDANT is 1-indexed with no gaps; session_id from
    # _prepare_session is already 0-indexed with no gaps (one entry per
    # successfully-loaded session), so a plain +1 is enough.
    session_index = (session_id + 1).astype(np.int64)
    np.save(os.path.join(input_dir, "session_index.npy"), session_index)

    # channel_locations.npy: (n_channel, 2) = (x, y/depth). Same probe assumed
    # across sessions in a group (same assumption run_EMD_batch_onMerged.m
    # makes) -- use the first session's positions. channel_pos entries have a
    # synthetic leading "1" column inserted by util.paths_from_KS when the
    # raw file is 2-column (x, z); 3-column raw files already have a real
    # leading (shank/depth-offset) column. Either way, columns 1:3 are (x, z).
    cp = channel_pos[0]
    if cp.shape[1] >= 3:
        channel_locations = cp[:, 1:3]
    else:
        channel_locations = cp[:, 0:2]
    np.save(
        os.path.join(input_dir, "channel_locations.npy"),
        channel_locations.astype(np.float64),
    )

    # spike_times/Unit{k}.npy: k = 0-indexed row in waveform_all/session_index
    # (= global position in np.concatenate(good_units)), in MILLISECONDS.
    original_ids = np.concatenate(good_units).flatten().astype(np.int64)
    spike_times_dir = os.path.join(input_dir, "spike_times")
    os.makedirs(spike_times_dir, exist_ok=True)

    session_spike_cache = {}
    for k in range(n_units):
        s = int(session_id[k])
        if s not in session_spike_cache:
            times = np.load(os.path.join(ks_dirs[s], "spike_times.npy")).flatten()
            clusters = np.load(os.path.join(ks_dirs[s], "spike_clusters.npy")).flatten()
            session_spike_cache[s] = (times, clusters)
        times, clusters = session_spike_cache[s]
        cid = original_ids[k]
        unit_times_samples = times[clusters == cid].astype(np.float64)
        unit_times_ms = (unit_times_samples / FS) * 1000.0
        np.save(os.path.join(spike_times_dir, f"Unit{k}.npy"), unit_times_ms)

    return n_units


def make_dant_settings(input_dir, output_dir, variant):
    """
    Mirrors pyDANT's own settings.json defaults, with "PETH" dropped from the
    feature lists (we don't have peth.npy prepared for every group -- see
    module docstring) and the feature lists themselves swapped in from
    DANT_VARIANTS[variant].
    """
    cfg = DANT_VARIANTS[variant]
    return {
        "path_to_data": input_dir,
        "output_folder": output_dir,
        "save_intermediate_results": True,
        "n_jobs": N_JOBS,
        "centering_waveforms": False,
        "spikeLocation": {
            "location_algorithm": "monopolar_triangulation",
            "n_nearest_channels": 20,
        },
        "waveformCorrection": {
            "n_nearest_channels": 38,
            "linear_correction": False,
            "n_templates": 2,
            "path_to_motion": "",
        },
        "autoCorr": {"window": 300, "binwidth": 1, "gaussian_sigma": 5},
        "ISI": {"window": 100, "binwidth": 1, "gaussian_sigma": 1},
        "motionEstimation": {
            "features": cfg["motionEstimation_features"],
            "max_iter": 15,
            "repeat_last_feature_set": True,
            "stop_early": True,
        },
        "clustering": {
            "max_distance": 100,
            "features": cfg["clustering_features"],
            "n_iter": 10,
            "weight_tol": 1e-8,
        },
        "autoCuration": {"auto_split": True},
    }


# ── per-variant scoring/saving (shared by every variant) ────────────────────


def score_and_save_variant(merged_dir, variant, sess, n_units):
    dant_dir = get_dant_dir(merged_dir, variant)
    output_dir = get_dant_raw_output_dir(merged_dir, variant)
    os.makedirs(dant_dir, exist_ok=True)

    input_dir = get_dant_shared_input_dir(merged_dir)
    dant_settings = make_dant_settings(input_dir, output_dir, variant)
    print(f"  [{variant}] Running pyDANT (features: clustering={dant_settings['clustering']['features']}) ...")

    # We keep pyDANT's own defaults untouched (e.g. clustering.max_distance),
    # so it can legitimately fail on some datasets -- e.g. zero cross-session
    # candidate pairs within max_distance on the raw, uncorrected first
    # motion-estimation pass, which crashes its LinearDiscriminantAnalysis fit
    # with a "0 samples" error. Rather than leaving this dataset/variant with
    # no output at all (which would (a) bias the cross-model comparison in
    # DANT's favour, since only the datasets it can solve would ever be
    # counted, and (b) leave no sentinel file, so it gets retried and fails
    # the same way on every future run), record it as zero matches -- the
    # same fair-comparison treatment any model finding zero matches gets.
    failure_reason = None
    try:
        runDANT(dant_settings)

        cluster_matrix_path = os.path.join(output_dir, "ClusterMatrix.npy")
        if not os.path.isfile(cluster_matrix_path):
            raise RuntimeError(f"pyDANT did not produce {cluster_matrix_path}")
        final_matches = np.load(cluster_matrix_path).astype(bool)
        if final_matches.shape != (n_units, n_units):
            raise RuntimeError(
                f"ClusterMatrix.npy shape {final_matches.shape} != expected ({n_units}, {n_units})"
            )
    except Exception as e:
        print(f"  [{variant}] pyDANT failed ({e}); recording as zero matches for a fair comparison.")
        traceback.print_exc()
        failure_reason = str(e)
        final_matches = np.zeros((n_units, n_units), dtype=bool)
        with open(os.path.join(dant_dir, "DANT_FAILURE.txt"), "w") as f:
            f.write(failure_reason + "\n")

    n_matches = int(np.sum(final_matches)) // 2
    print(f"  [{variant}] {n_matches} matches found")

    param = sess["param"]
    session_id = sess["session_id"]
    session_switch = sess["session_switch"]
    good_units = sess["good_units"]
    channel_pos = sess["channel_pos"]
    waveform = sess["waveform"]

    clus_info = {
        "good_units": good_units,
        "session_switch": session_switch,
        "session_id": session_id,
        "original_ids": np.concatenate(good_units),
    }
    extracted_wave_properties = ov.extract_parameters(
        waveform, channel_pos, clus_info, param
    )

    # ── performance metrics (AUC against functional scores) ──────────────────
    # Same seven functional scores computed identically for DUM/UMPy/EMD, so
    # DANT (both variants) lands in the same comparison framework.
    functional_scores = {}
    try:
        isicorr = test.ISI_correlations(param)
        auc_isi = test.AUC(final_matches, isicorr, session_id)
        print(f"  [{variant}] AUC (ISI correlations):            {auc_isi:.3f}")
        functional_scores["ISI_correlations"] = isicorr

        isikl = test.ISI_KL_divergence(param)
        auc_isikl = test.AUC(final_matches, -isikl, session_id)
        print(f"  [{variant}] AUC (ISI KL divergence):           {auc_isikl:.3f}")
        functional_scores["ISI_KL_divergence"] = isikl

        isiwass = test.ISI_wasserstein_distance(param)
        auc_isiwass = test.AUC(final_matches, -isiwass, session_id)
        print(f"  [{variant}] AUC (ISI Wasserstein distance):    {auc_isiwass:.3f}")
        functional_scores["ISI_wasserstein_distance"] = isiwass

        refpopcorr = test.refpop_correlations(param, matches=final_matches)
        auc_refpop = test.AUC(final_matches, refpopcorr, session_id)
        print(f"  [{variant}] AUC (ref. pop. correlation):       {auc_refpop:.3f}")
        functional_scores["refpop_correlations"] = refpopcorr

        frdiff = test.FR_diff(param)
        auc_fr = test.AUC(final_matches, -frdiff, session_id)
        print(f"  [{variant}] AUC (firing rate difference):      {auc_fr:.3f}")
        functional_scores["FR_diff"] = frdiff

        cvdiff = test.ISI_CV_diff(param)
        auc_cv = test.AUC(final_matches, -cvdiff, session_id)
        print(f"  [{variant}] AUC (ISI CV difference):           {auc_cv:.3f}")
        functional_scores["ISI_CV_diff"] = cvdiff

        try:
            natimcorr = test.natim_correlations(param, merged_architecture=True)
            auc_natim = test.AUC(final_matches, natimcorr, session_id)
            print(f"  [{variant}] AUC (nat. image correlations):     {auc_natim:.3f}")
            functional_scores["natim_correlations"] = natimcorr
        except Exception:
            pass  # natim data may not be available for every session
    except Exception as e:
        print(f"  [{variant}] WARNING: functional score computation failed: {e}")
        functional_scores = {}

    su.save_to_output(
        dant_dir,
        {},
        np.argwhere(final_matches),
        final_matches.astype(float),
        extracted_wave_properties["avg_centroid"],
        extracted_wave_properties["avg_waveform"],
        extracted_wave_properties["avg_waveform_per_tp"],
        extracted_wave_properties["max_site"],
        final_matches.astype(float),
        final_matches,
        clus_info,
        param,
        UIDs=None,
        matches_curated=None,
        save_match_table=True,
        functional_scores=functional_scores if functional_scores else None,
    )
    su.save_auc_summary(
        dant_dir, test.auc_summary_from_functional_scores(functional_scores, final_matches, session_id)
    )

    fig, ax = plt.subplots(figsize=(5, 5))
    im = ax.imshow(final_matches, cmap="viridis", aspect="auto")
    title = f"{variant} matches (n={n_matches})"
    if failure_reason:
        title += "\n(pyDANT run failed -- recorded as zero matches, see DANT_FAILURE.txt)"
    ax.set_title(title)
    ax.set_xlabel("Unit")
    ax.set_ylabel("Unit")
    fig.colorbar(im, ax=ax)
    fig.tight_layout()
    fig.savefig(os.path.join(dant_dir, "MatchingOverview.png"), dpi=150)
    plt.close(fig)

    if functional_scores:
        score_meta = {
            "ISI_correlations": ("ISI correlations", "viridis"),
            "ISI_KL_divergence": ("ISI KL divergence", "magma"),
            "ISI_wasserstein_distance": ("ISI Wasserstein distance", "magma"),
            "refpop_correlations": ("Ref. pop. correlations", "viridis"),
            "FR_diff": ("Firing rate difference", "magma"),
            "ISI_CV_diff": ("ISI CV difference", "magma"),
            "natim_correlations": ("Nat. image correlations", "viridis"),
        }
        keys = [k for k in score_meta if k in functional_scores]
        fig, axes = plt.subplots(1, len(keys), figsize=(5 * len(keys), 5))
        if len(keys) == 1:
            axes = [axes]
        for ax, key in zip(axes, keys):
            title, cmap = score_meta[key]
            im = ax.imshow(functional_scores[key], cmap=cmap, aspect="auto")
            ax.set_title(title)
            fig.colorbar(im, ax=ax)
        fig.tight_layout()
        fig.savefig(os.path.join(dant_dir, "FunctionalScores.png"), dpi=150)
        plt.close(fig)

    print(f"  [{variant}] Results saved to: {dant_dir}")


# ── main per-group pipeline ──────────────────────────────────────────────────


def run_group(merged_dir, pending_variants):
    print(f"\n--- DANT: {merged_dir}  (pending: {pending_variants})")

    sess = _prepare_session(merged_dir)
    if sess is None:
        print("  ERROR: _prepare_session failed, skipping.")
        return

    input_dir = get_dant_shared_input_dir(merged_dir)
    # Shared input is identical for every variant -- build once per group,
    # reused whichever variant(s) are actually pending.
    if not os.path.isfile(os.path.join(input_dir, "waveform_all.npy")):
        print("  Building shared pyDANT input arrays ...")
        n_units = build_dant_input(sess, input_dir)
    else:
        n_units = sess["waveform"].shape[0]
    print(f"  {n_units} units across {sess['param']['n_sessions']} session(s)")

    for variant in pending_variants:
        try:
            score_and_save_variant(merged_dir, variant, sess, n_units)
        except Exception as e:
            print(f"  [{variant}] FAILED: {e}")
            traceback.print_exc()


def main():
    print(f"Scanning for merged-data groups under:\n  {BASE_INPUT}\n")
    groups = find_merged_groups()
    if not groups:
        print("No merged-data groups found.")
        return
    print(f"Found {len(groups)} group(s).\n")

    for i, merged_dir in enumerate(groups):
        print(f"\n[{i + 1}/{len(groups)}] {merged_dir}")

        pending = [v for v in DANT_VARIANTS if not dant_results_exist(merged_dir, v)]
        if not pending:
            print("  Skipping (all variants exist and are fresh).")
            continue

        lock_path = get_dant_lock_path(merged_dir)
        with batch_lock.try_lock(lock_path, redo_from_date=REDO_FROM_DATE) as acquired:
            if not acquired:
                print(f"  Skipping (already being processed by another run): {lock_path}")
                continue

            # re-check now that we hold the lock: another machine may have
            # finished some/all pending variants while we were waiting
            pending = [v for v in DANT_VARIANTS if not dant_results_exist(merged_dir, v)]
            if not pending:
                print("  Skipping (completed by another run).")
                continue

            try:
                run_group(merged_dir, pending)
            except Exception as e:
                print(f"  Group FAILED: {e}")
                traceback.print_exc()

    print("\nAll done.")


if __name__ == "__main__":
    main()
