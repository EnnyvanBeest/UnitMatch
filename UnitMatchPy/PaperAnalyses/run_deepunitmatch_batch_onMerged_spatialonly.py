# Batch wrapper: runs a "spatial model only" ablation on the *merged* dataset
# (see run_deepunitmatch_batch_onMerged.py for how that tree is built and how
# sessions/good units are derived) -- i.e. what happens to matching
# performance if every non-spatial predictor is removed and only the
# centroid-distance metric is left to decide matches.
#
# DeepUnitMatch and UMPy share one matching pipeline (run_umpy_core), so
# without their non-spatial scores they are the same model: a single
# condition, saved as UMPy_spatialonly/ -- run_umpy_core(...,
# to_use=["centroid_dist"]), which drops every score except centroid_dist
# from total_score/candidate_pairs/drift correction/the Bayes predictors.
# (The previous DUM_spatialonly condition used the DeepUnitMatch-specific
# legacy pipeline with a constant similarity; it no longer applies.)
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
import pipeline_config as cfg
import run_deepunitmatch_batch_onMerged as base_batch

# See batch_lock.sentinel_is_fresh() / run_deepunitmatch_batch_onMerged.py's
# REDO_FROM_DATE for what this does. A dataset/method combo is skipped once
# its MatchingOverview.png exists and is at least this new.
REDO_FROM_DATE = cfg.REDO_FROM_DATE  # see pipeline_config.py

# One condition: DUM and UMPy are identical once only centroid_dist is left.
METHODS = ("UMPy",)


# ── path helpers ─────────────────────────────────────────────────────────────


def get_spatialonly_save_dir(merged_dir, method):
    """Output dir for a given merged-data group + method (see METHODS)."""
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
        description="Run a spatial-only (centroid-distance-only) ablation of the shared DeepUnitMatch/UMPy matching pipeline on the merged dataset."
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
            for method in METHODS
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
                for method in METHODS
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
                    base_batch.run_umpy_core(
                        sess, save_dir, label=label, to_use=["centroid_dist"]
                    )
                except Exception as e:
                    print(f"  {label} FAILED: {e}")
                    traceback.print_exc()

    print("\nAll done.")


if __name__ == "__main__":
    main()
