# Batch wrapper: runs a "score swap" test of DeepUnitMatch and UMPy on the
# *merged* dataset (see run_deepunitmatch_batch_onMerged.py for how that tree
# is built and how sessions/good units are derived) -- i.e. feed each
# pipeline the OTHER pipeline's score, to separate "the score itself" from
# "everything else that differs between the two pipelines' architectures"
# (per-session-pair vs global vectorised processing, how drift correction is
# bootstrapped, how the adaptive prior/threshold is estimated, etc.).
#
# DUM_totalscore/: run_deep_unit_match_core(..., use_umpy_totalscore=True).
# The DNN is never loaded or run; sim_matrix is instead UMPy's own TotalScore
# -- computed by running ov.extract_metric_scores on this same session's
# extracted waveform properties (see that parameter's docstring in
# run_deepunitmatch_batch_onMerged.py). Every other step of DUM's own
# architecture (the drift-correction pre-pass, best-match de-duplication,
# adaptive-prior estimation, the final per-session-pair Bayes step, AUC/
# figures) runs completely unmodified, just fed this real, non-constant score
# in place of the DNN's.
#
# UMPy_simscore/: run_umpy_core(..., to_use=["similarity"], model=<DUM's
# default trained model>). This runs DUM's own preprocessing (get_snippets,
# kept_idx re-sync, test.inference) to get a DNN similarity score, which is
# folded into UMPy's own ov.extract_metric_scores as an extra 'similarity'
# entry and, via to_use=["similarity"], made the sole driver of UMPy's own
# total_score/candidate_pairs/drift-correction/Bayes predictors -- in place of
# UMPy's own six metrics. Every other step of UMPy's own architecture runs
# unmodified.
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
# REDO_FROM_DATE for what this does. A dataset/method combo is skipped once
# its MatchingOverview.png exists and is at least this new.
REDO_FROM_DATE = datetime.datetime(2026, 8, 22, 0, 0, 0)

CONDITIONS = ("DUM_totalscore", "UMPy_simscore")


# ── path helpers ─────────────────────────────────────────────────────────────


def get_scoreswap_save_dir(merged_dir, condition):
    """Output dir for a given merged-data group + condition (one of CONDITIONS)."""
    subfolder = os.path.relpath(os.path.dirname(merged_dir), base_batch.BASE_INPUT)
    return os.path.join(base_batch.BASE_OUTPUT, subfolder, condition)


def scoreswap_results_exist(merged_dir, condition):
    """Return True when the sentinel output file is present and fresh for this dataset/condition combo."""
    sentinel = os.path.join(
        get_scoreswap_save_dir(merged_dir, condition), "MatchingOverview.png"
    )
    return batch_lock.sentinel_is_fresh(sentinel, REDO_FROM_DATE)


def get_group_lock_path(merged_dir):
    """
    Lock file marking 'a run is currently processing the score-swap test for
    this group', so multiple machines pointed at the same BASE_INPUT/
    BASE_OUTPUT can split work across groups without double-processing one.
    See batch_lock.py. Named distinctly from the other batch scripts' locks so
    they can all run on the same group concurrently.
    """
    subfolder = os.path.relpath(os.path.dirname(merged_dir), base_batch.BASE_INPUT)
    return os.path.join(base_batch.BASE_OUTPUT, subfolder, ".processing_scoreswap.lock")


# ── entry point ───────────────────────────────────────────────────────────────


def parse_args():
    parser = argparse.ArgumentParser(
        description="Run DUM with UMPy's TotalScore, and UMPy with DUM's similarity score, on the merged dataset."
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

    model = None  # DUM's default trained model, loaded once and reused (only UMPy_simscore needs it)

    for i, merged_dir in enumerate(groups):
        print(f"\n[{i + 1}/{len(groups)}] {merged_dir}")

        pending = [
            condition
            for condition in CONDITIONS
            if not scoreswap_results_exist(merged_dir, condition)
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
                if not scoreswap_results_exist(merged_dir, condition)
            ]
            if not pending:
                print("  Skipping (completed by another run).")
                continue

            sess = base_batch._prepare_session(merged_dir)
            if sess is None:
                continue

            if model is None and "UMPy_simscore" in pending:
                print("Loading DeepUnitMatch model …")
                model = test.load_trained_model(device=base_batch.DEVICE)

            for condition in pending:
                save_dir = get_scoreswap_save_dir(merged_dir, condition)
                try:
                    if condition == "DUM_totalscore":
                        base_batch.run_deep_unit_match_core(
                            sess,
                            save_dir,
                            model=None,
                            label=condition,
                            use_umpy_totalscore=True,
                        )
                    else:
                        base_batch.run_umpy_core(
                            sess,
                            save_dir,
                            label=condition,
                            to_use=["similarity"],
                            model=model,
                        )
                except Exception as e:
                    print(f"  {condition} FAILED: {e}")
                    traceback.print_exc()

    print("\nAll done.")


if __name__ == "__main__":
    main()
