# Batch wrapper: runs DeepUnitMatch and UMPy on the *merged* dataset (see
# run_deepunitmatch_batch_onMerged.py for how that tree is built and how
# sessions/good units are derived) with drift correction disabled entirely,
# to see how much each pipeline's normal performance actually depends on it.
# Everything else about each pipeline's default scoring stays exactly as it
# normally is (DUM still uses the DNN similarity score; UMPy still uses its
# own six metrics) -- this ablation isolates drift correction specifically,
# independent of the score-source questions the spatialonly/scoreswap batch
# scripts explore.
#
# DUM_nodrift/: run_deep_unit_match_core(..., apply_drift_correction=False).
# mf.drift_n_sessions is never called, so avg_centroid/avg_waveform_per_tp
# stay exactly as extract_parameters produced them -- centroid_dist and every
# downstream step use uncorrected waveforms. The drift-correction pre-pass
# itself still runs (it also builds the supervised match/non-match labels the
# final Bayes step needs, regardless of drift correction), only the actual
# correction step is skipped.
#
# UMPy_nodrift/: run_umpy_core(..., niter=1). ov.extract_metric_scores' own
# niter loop is what applies drift correction (niter=2 does one pass); niter=1
# runs the loop once and skips it entirely, computing total_score/
# candidate_pairs/scores_to_include from the raw, uncorrected waveforms.
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
from DeepUnitMatch.testing import test

# See batch_lock.sentinel_is_fresh() / run_deepunitmatch_batch_onMerged.py's
# REDO_FROM_DATE for what this does. A dataset/condition combo is skipped once
# its MatchingOverview.png exists and is at least this new.
REDO_FROM_DATE = datetime.datetime(2026, 8, 22, 0, 0, 0)

CONDITIONS = ("DUM_nodrift", "UMPy_nodrift")


# ── path helpers ─────────────────────────────────────────────────────────────


def get_nodrift_save_dir(merged_dir, condition):
    """Output dir for a given merged-data group + condition (one of CONDITIONS)."""
    subfolder = os.path.relpath(os.path.dirname(merged_dir), base_batch.BASE_INPUT)
    return os.path.join(base_batch.BASE_OUTPUT, subfolder, condition)


def nodrift_results_exist(merged_dir, condition):
    """Return True when the sentinel output file is present and fresh for this dataset/condition combo."""
    sentinel = os.path.join(
        get_nodrift_save_dir(merged_dir, condition), "MatchingOverview.png"
    )
    return batch_lock.sentinel_is_fresh(sentinel, REDO_FROM_DATE)


def get_group_lock_path(merged_dir):
    """
    Lock file marking 'a run is currently processing the no-drift-correction
    ablation for this group', so multiple machines pointed at the same
    BASE_INPUT/BASE_OUTPUT can split work across groups without
    double-processing one. See batch_lock.py. Named distinctly from the other
    batch scripts' locks so they can all run on the same group concurrently.
    """
    subfolder = os.path.relpath(os.path.dirname(merged_dir), base_batch.BASE_INPUT)
    return os.path.join(base_batch.BASE_OUTPUT, subfolder, ".processing_nodrift.lock")


# ── entry point ───────────────────────────────────────────────────────────────


def parse_args():
    parser = argparse.ArgumentParser(
        description="Run DeepUnitMatch and UMPy on the merged dataset with drift correction disabled."
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

    model = None  # DUM's default trained model, loaded once and reused across all groups

    for i, merged_dir in enumerate(groups):
        print(f"\n[{i + 1}/{len(groups)}] {merged_dir}")

        pending = [
            condition
            for condition in CONDITIONS
            if not nodrift_results_exist(merged_dir, condition)
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
                condition
                for condition in CONDITIONS
                if not nodrift_results_exist(merged_dir, condition)
            ]
            if not pending:
                print("  Skipping (completed by another run).")
                continue

            sess = base_batch._prepare_session(merged_dir)
            if sess is None:
                continue

            if model is None and "DUM_nodrift" in pending:
                print("Loading DeepUnitMatch model …")
                model = test.load_trained_model(device=base_batch.DEVICE)

            for condition in pending:
                save_dir = get_nodrift_save_dir(merged_dir, condition)
                try:
                    if condition == "DUM_nodrift":
                        base_batch.run_deep_unit_match_core(
                            sess,
                            save_dir,
                            model=model,
                            label=condition,
                            apply_drift_correction=False,
                        )
                    else:
                        base_batch.run_umpy_core(
                            sess, save_dir, label=condition, niter=1
                        )
                except Exception as e:
                    print(f"  {condition} FAILED: {e}")
                    traceback.print_exc()

    print("\nAll done.")


if __name__ == "__main__":
    main()
