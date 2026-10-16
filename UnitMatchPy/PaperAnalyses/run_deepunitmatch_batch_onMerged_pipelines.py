# Batch wrapper: pipeline variants on the *merged* dataset, to find out which
# parts of the unified matching pipeline (run_umpy_core, used for both
# DeepUnitMatch/ and UMPy/) matter, for DeepUnitMatch and for UMPy, compared
# with the previous DeepUnitMatch-specific pipeline (DUM_legacy/).
#
# The two main differences between the pipelines, varied as a 2x2 with
# everything else unified (DeepUnitMatch with the paper model,
# pipeline_config.PAPER_MODEL_PATH):
#   A. candidate pairs (drift correction, Naive Bayes labels and prior):
#      per location = one total-score crossing (similarity 5:1 distance),
#      after per-session-pair alignment; per pair (as the DUM_legacy
#      pipeline) = per session pair, the method's similarity above that
#      pair's own same/different-unit crossing in both directions, within
#      max_dist, mutual best.
#   B. Naive Bayes + prior: one for the location, or one per session pair.
#
#                              Bayes per location     Bayes per pair
#   candidates per location    DeepUnitMatch/         DUM_pairBayes/
#   candidates per pair        DUM_pairCandidates/    DUM_pairCandidates_pairBayes/
#
# and C. the similarity:distance weight in the total score, 1:1 instead of
# 5:1 (DUM_weights1to1/).
#
# The same changes are applied to UMPy (UMPy_pairCandidates/,
# UMPy_pairBayes/, UMPy_pairCandidates_pairBayes/, UMPy_weights1to1/), so
# each change can be judged on whether it helps DeepUnitMatch specifically.
# UMPy's counterparts: the per-pair selection runs on UMPy's total score over its
# five waveform metrics (without the centroid distance, as the DNN
# similarity), with the same unit positions for the distance filter; 1:1 =
# each waveform metric weighted 0.2, so together they count as much as the
# distance (DUM's 5:1 mirrors UMPy's five metrics against one distance).
# Every run also writes candidate_pairs.json (number of candidate pairs: how
# stringent each selection is).
#
# Each folder sits alongside DeepUnitMatch/ and UMPy/ for the same location;
# sql.py / fast_testing.py --group pipelines.

import argparse
import os
import sys
import traceback

import matplotlib

matplotlib.use("Agg")  # non-interactive backend for batch runs

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _HERE)
sys.path.insert(0, os.path.dirname(_HERE))
sys.path.insert(0, os.path.join(_HERE, "DeepUnitMatch"))

import batch_lock
import pipeline_config as cfg
import run_deepunitmatch_batch_onMerged as base_batch

REDO_FROM_DATE = cfg.REDO_FROM_DATE  # see pipeline_config.py

# condition -> (method, run_umpy_core options); order = run order per location
UMPY_WEIGHTS_1TO1 = {score: 0.2 for score in base_batch.UMPY_SCORES if score != "centroid_dist"}
UMPY_WEIGHTS_1TO1["centroid_dist"] = 1.0

CONDITIONS = {
    "DUM_pairCandidates": ("DUM", dict(candidates="pair")),
    "DUM_pairBayes": ("DUM", dict(bayes="pair")),
    "DUM_pairCandidates_pairBayes": ("DUM", dict(candidates="pair", bayes="pair")),
    "DUM_weights1to1": ("DUM", dict(score_weights={"similarity": 1, "centroid_dist": 1})),
    "UMPy_pairCandidates": ("UMPy", dict(candidates="pair")),
    "UMPy_pairBayes": ("UMPy", dict(bayes="pair")),
    "UMPy_pairCandidates_pairBayes": ("UMPy", dict(candidates="pair", bayes="pair")),
    "UMPy_weights1to1": ("UMPy", dict(score_weights=UMPY_WEIGHTS_1TO1)),
}


def get_save_dir(merged_dir, condition):
    subfolder = os.path.relpath(os.path.dirname(merged_dir), base_batch.BASE_INPUT)
    return os.path.join(base_batch.BASE_OUTPUT, subfolder, condition)


def results_exist(merged_dir, condition):
    sentinel = os.path.join(get_save_dir(merged_dir, condition), base_batch.SENTINEL)
    return batch_lock.sentinel_is_fresh(sentinel, REDO_FROM_DATE)


def get_group_lock_path(merged_dir):
    subfolder = os.path.relpath(os.path.dirname(merged_dir), base_batch.BASE_INPUT)
    return os.path.join(base_batch.BASE_OUTPUT, subfolder, ".processing_pipelines.lock")


def parse_args():
    parser = argparse.ArgumentParser(description="Pipeline variants of DeepUnitMatch/UMPy on the merged dataset.")
    parser.add_argument("--only", nargs="+", choices=sorted(CONDITIONS), help="restrict to these conditions")
    parser.add_argument("--write-matlab-compat", action="store_true",
                        help="Also write a MATLAB-compatible UnitMatch.mat from the Python outputs.")
    return parser.parse_args()


def main():
    args = parse_args()
    base_batch.WRITE_MATLAB_COMPAT = args.write_matlab_compat
    conditions = [c for c in CONDITIONS if not args.only or c in args.only]

    groups = base_batch.find_merged_groups()
    print(f"Found {len(groups)} group(s); conditions {conditions}\n")
    model = None  # the paper model, loaded once when first needed

    for i, merged_dir in enumerate(groups):
        print(f"\n[{i + 1}/{len(groups)}] {merged_dir}")
        if all(results_exist(merged_dir, c) for c in conditions):
            print("  Skipping (results exist and are fresh).")
            continue

        with batch_lock.try_lock(get_group_lock_path(merged_dir)) as acquired:
            if not acquired:
                print("  Skipping (already being processed by another run).")
                continue
            pending = [c for c in conditions if not results_exist(merged_dir, c)]
            if not pending:
                continue

            sess = base_batch._prepare_session(merged_dir)
            if sess is None:
                continue  # logged by _prepare_session

            for condition in pending:
                method, options = CONDITIONS[condition]
                save_dir = get_save_dir(merged_dir, condition)
                try:
                    if method == "DUM":
                        if model is None:
                            model = base_batch.load_dum_model()
                        base_batch.run_dum_core(sess, save_dir, model, label=condition, **options)
                    else:
                        base_batch.run_umpy_core(sess, save_dir, label=condition, **options)
                except Exception as e:
                    print(f"  {condition} FAILED: {e}")  # logged by the run functions
                    traceback.print_exc()

    print("\nAll done.")


if __name__ == "__main__":
    main()
