# Batch wrapper: runs a "spatial model only" ablation of DeepUnitMatch and
# UMPy on the *merged* dataset (see run_deepunitmatch_batch_onMerged.py for how
# that tree is built and how sessions/good units are derived) -- i.e. what
# happens to matching performance if every non-spatial predictor is removed
# and only the centroid-distance metric is left to decide matches.
#
# DUM (DUM_spatialonly/): run_deep_unit_match_core(..., spatial_only=True,
# total_score_for_threshold=True). The DNN is never loaded or run; sim_matrix
# is instead a constant all-ones array, so every downstream step (the
# drift-correction pre-pass, best-match de-duplication inside
# mf.drift_n_sessions, the final per-session-pair Bayes step) runs completely
# unmodified with a uniformly uninformative 'similarity' input, leaving
# centroid distance as the only real signal driving matches. Two places would
# otherwise choke on that constant input, so both are patched:
#   - test.get_threshold (inside the pre-pass' test.get_matches, setting its
#     Prob > threshold cutoff) has no well-defined on/off-diagonal KDE
#     intersection for constant input, so it special-cases that to return
#     threshold=0 -- every pair then passes the check, leaving get_matches'
#     spatial (position-based) filter as the only thing selecting
#     drift-correction candidates, so drift correction still runs off spatial
#     proximity alone.
#   - mf.get_threshold, used per session-pair to set the adaptive Bayes prior
#     (n_expected_matches = count of pairs above threshold), is fundamentally
#     unable to split a *constant* sim_mat into a meaningful subset -- any
#     threshold value makes every pair (or none) exceed it, driving the prior
#     to an extreme that forces match probability to 0 or 1 for everyone,
#     regardless of centroid distance. total_score_for_threshold=True fixes
#     this by feeding that calculation mf.get_total_score(scores_to_incl, ...)
#     (the normalised sum of 'similarity' + 'distance', exactly how UMPy forms
#     its own total_score) instead of raw sim_mat -- the constant 'similarity'
#     term drops out under min-max normalisation, leaving the real
#     centroid_dist signal to drive the threshold.
# See run_deep_unit_match_core's docstring in run_deepunitmatch_batch_onMerged.py
# for the exact mechanism of both fixes.
#
# UMPy (UMPy_spatialonly/): run_umpy_core(..., to_use=["centroid_dist"]),
# using ov.extract_metric_scores' existing to_use argument to drop every
# metric except centroid_dist from total_score/candidate_pairs/drift-
# correction/the Bayes predictors -- so, like DUM here, drift correction,
# threshold/prior estimation, and the final match decision are all driven by
# centroid distance alone.
#
# Both variants therefore still get real drift correction, but bootstrapped
# from different spatial signals: UMPy_spatialonly uses its own iterative
# niter=2 loop over waveform-derived centroid_dist, while DUM_spatialonly
# falls back to get_matches' raw electrode/unit-position filter (see above).
# They're not expected to land on identical numbers, but comparing the two is
# still a useful sanity check that the two Bayes implementations behave
# similarly once stripped down to spatial information alone.
#
# Each such folder sits alongside the DeepUnitMatch/ and UMPy/ subfolders that
# run_deepunitmatch_batch_onMerged.py writes for the same dataset.

import os
import sys
import datetime
import traceback
import argparse

import matplotlib

matplotlib.use("Agg")  # non-interactive backend for batch runs

# ── project paths ───────────────────────────────────────────────────────────
_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _HERE)
sys.path.insert(0, os.path.dirname(_HERE))
sys.path.insert(0, os.path.join(_HERE, "DeepUnitMatch"))

import batch_lock
import run_deepunitmatch_batch_onMerged as base_batch

# See batch_lock.sentinel_is_fresh() / run_deepunitmatch_batch_onMerged.py's
# REDO_FROM_DATE for what this does. A dataset/method combo is skipped once
# its MatchingOverview.png exists and is at least this new.
REDO_FROM_DATE = datetime.datetime(2026, 8, 21, 0, 0, 0)


# ── path helpers ─────────────────────────────────────────────────────────────


def get_spatialonly_save_dir(merged_dir, method):
    """Output dir for a given merged-data group + method ('DUM' or 'UMPy')."""
    subfolder = os.path.relpath(os.path.dirname(merged_dir), base_batch.BASE_INPUT)
    return os.path.join(base_batch.BASE_OUTPUT, subfolder, f"{method}_spatialonly")


def spatialonly_results_exist(merged_dir, method):
    """Return True when the sentinel output file is present and fresh for this dataset/method combo."""
    sentinel = os.path.join(
        get_spatialonly_save_dir(merged_dir, method), "MatchingOverview.png"
    )
    return batch_lock.sentinel_is_fresh(sentinel, REDO_FROM_DATE)


def get_group_lock_path(merged_dir):
    """
    Lock file marking 'a run is currently processing the spatial-only ablation
    for this group', so multiple machines pointed at the same BASE_INPUT/
    BASE_OUTPUT can split work across groups without double-processing one.
    See batch_lock.py. Named distinctly from the other batch scripts' locks so
    they can all run on the same group concurrently.
    """
    subfolder = os.path.relpath(os.path.dirname(merged_dir), base_batch.BASE_INPUT)
    return os.path.join(base_batch.BASE_OUTPUT, subfolder, ".processing_spatialonly.lock")


# ── entry point ───────────────────────────────────────────────────────────────


def parse_args():
    parser = argparse.ArgumentParser(
        description="Run a spatial-only (centroid-distance-only) ablation of DeepUnitMatch and UMPy on the merged dataset."
    )
    parser.add_argument(
        "--write-matlab-compat",
        action="store_true",
        help="Also write a MATLAB-compatible UnitMatch.mat from the Python outputs.",
    )
    return parser.parse_args()


def main():
    args = parse_args()
    base_batch.WRITE_MATLAB_COMPAT = args.write_matlab_compat

    print(f"Scanning for merged-data groups under:\n  {base_batch.BASE_INPUT}\n")
    groups = base_batch.find_merged_groups()
    if not groups:
        print("No merged-data groups found.")
        return
    print(f"Found {len(groups)} group(s).\n")

    for i, merged_dir in enumerate(groups):
        print(f"\n[{i + 1}/{len(groups)}] {merged_dir}")

        pending = [
            method
            for method in ("DUM", "UMPy")
            if not spatialonly_results_exist(merged_dir, method)
        ]
        if not pending:
            print("  Skipping (results exist and are fresh).")
            continue

        lock_path = get_group_lock_path(merged_dir)
        with batch_lock.try_lock(lock_path) as acquired:
            if not acquired:
                print(f"  Skipping (already being processed by another run): {lock_path}")
                continue

            # re-check now that we hold the lock: another machine may have
            # finished this group while we were scanning/waiting for the lock
            pending = [
                method
                for method in ("DUM", "UMPy")
                if not spatialonly_results_exist(merged_dir, method)
            ]
            if not pending:
                print("  Skipping (completed by another run).")
                continue

            sess = base_batch._prepare_session(merged_dir)
            if sess is None:
                continue

            for method in pending:
                save_dir = get_spatialonly_save_dir(merged_dir, method)
                label = f"{method}_spatialonly"
                try:
                    if method == "DUM":
                        base_batch.run_deep_unit_match_core(
                            sess,
                            save_dir,
                            model=None,
                            label=label,
                            spatial_only=True,
                            total_score_for_threshold=True,
                        )
                    else:
                        base_batch.run_umpy_core(
                            sess, save_dir, label=label, to_use=["centroid_dist"]
                        )
                except Exception as e:
                    print(f"  {label} FAILED: {e}")
                    traceback.print_exc()

    print("\nAll done.")


if __name__ == "__main__":
    main()
