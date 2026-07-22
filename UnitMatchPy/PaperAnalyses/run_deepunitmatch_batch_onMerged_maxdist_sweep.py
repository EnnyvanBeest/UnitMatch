# Batch wrapper: runs DeepUnitMatch and UMPy on the *merged* dataset (see
# run_deepunitmatch_batch_onMerged.py for how that tree is built and how
# sessions/good units are derived) once per candidate spatial-filter radius
# (param["max_dist"]), instead of the single default value, to test how
# sensitive matching performance is to that radius.
#
# max_dist plays two roles in the shared UnitMatchPy code (metric_functions.py
# / overlord.py), both relevant here:
#   1. Candidate-pair gate: pairs farther apart than max_dist are never
#      considered at all (overlord.extract_metric_scores' include_these_pairs,
#      and -- since the DUM adaptive-prior fix -- DUM's own per-pair
#      include_these_pairs_idx in run_deep_unit_match_core).
#   2. Normalisation scale for the centroid-distance score itself
#      (metric_functions.centroid_metrics / get_euclidean_metrics_chunked:
#      centroid_dist = 1 - (dist - min) / (max_dist - min)). Both of those
#      functions zero out centroid_dist entirely when max_dist is non-finite
#      (the normalising denominator becomes non-finite) -- so the MAX_DIST_INF
#      condition below isn't just "no candidate-pair gate", it also flattens
#      the distance-based Bayes predictor to a constant. That's a deliberate,
#      known tradeoff for this sweep (matching in that condition relies
#      entirely on non-spatial predictors), not a bug -- see the conversation
#      this script came out of for the full reasoning.
#
# neighbour_dist (used only inside metric_functions.get_threshold, by both
# UMPy's extract_metric_scores and DUM's adaptive-prior code) is a *separate*
# radius: the tight neighbourhood used to calibrate the match/no-match score
# threshold from "hard", spatially-plausible non-matches, deliberately
# narrower than max_dist's broad candidate net so the calibration isn't
# diluted by trivially-obvious non-matches. We keep that design intent across
# the sweep by *capping* it rather than coupling it 1:1 to max_dist: at each
# sweep point neighbour_dist = min(default neighbour_dist, max_dist) -- so it
# only ever shrinks (max_dist=20) and otherwise stays at its sensible default
# (max_dist=50/100/inf), rather than growing past the point where it would be
# diluted by mostly-obvious non-matches.
#
# Each such folder sits alongside the DeepUnitMatch/ and UMPy/ subfolders that
# run_deepunitmatch_batch_onMerged.py writes for the same dataset, using the
# exact same extraction/matching/saving pipeline (run_deep_unit_match_core /
# run_umpy_core), just against a modified param dict.

import os
import sys
import datetime
import traceback
import argparse

import numpy as np
import matplotlib

matplotlib.use("Agg")  # non-interactive backend for batch runs

# ── project paths ───────────────────────────────────────────────────────────
_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _HERE)
sys.path.insert(0, os.path.dirname(_HERE))
sys.path.insert(0, os.path.join(_HERE, "DeepUnitMatch"))

import batch_lock
import UnitMatchPy.default_params as default_params
import run_deepunitmatch_batch_onMerged as base_batch
from DeepUnitMatch.testing import test

MAX_DIST_VALUES = [20, 50, 100, float("inf")]
DEFAULT_NEIGHBOUR_DIST = default_params.get_default_param()["neighbour_dist"]

# See batch_lock.sentinel_is_fresh() / run_deepunitmatch_batch_onMerged.py's
# REDO_FROM_DATE for what this does. A dataset/max_dist/method combo is
# skipped once its MatchingOverview.png exists and is at least this new.
REDO_FROM_DATE = datetime.datetime(2026, 7, 22, 19, 0, 0)


# ── sweep-point helpers ──────────────────────────────────────────────────────


def neighbour_dist_for(max_dist):
    """
    Cap neighbour_dist at max_dist rather than coupling them 1:1 -- see the
    module docstring for why: neighbour_dist should stay at its conservative
    default except where max_dist itself is narrower (in which case a wider
    neighbour_dist would calibrate from pairs that aren't even valid
    candidates, which is nonsensical).
    """
    return min(DEFAULT_NEIGHBOUR_DIST, max_dist)


def format_dist(value):
    return "inf" if np.isinf(value) else str(int(value))


def sess_with_max_dist(sess, max_dist):
    """
    Shallow-copy sess with param["max_dist"]/param["neighbour_dist"] overridden.
    run_deep_unit_match_core/run_umpy_core each do their own
    copy.deepcopy(sess["param"]) internally, so this shallow override is
    isolated per sweep point without needing to touch the original sess.
    """
    sess_mod = dict(sess)
    sess_mod["param"] = dict(
        sess["param"],
        max_dist=max_dist,
        neighbour_dist=neighbour_dist_for(max_dist),
    )
    return sess_mod


# ── path helpers ─────────────────────────────────────────────────────────────


def get_maxdist_save_dir(merged_dir, max_dist, method):
    """Output dir for a given merged-data group + max_dist + method ('DUM' or 'UMPy')."""
    subfolder = os.path.relpath(os.path.dirname(merged_dir), base_batch.BASE_INPUT)
    return os.path.join(
        base_batch.BASE_OUTPUT, subfolder, f"{method}_maxdist={format_dist(max_dist)}"
    )


def maxdist_results_exist(merged_dir, max_dist, method):
    """Return True when the sentinel output file is present and fresh for this dataset/max_dist/method combo."""
    sentinel = os.path.join(
        get_maxdist_save_dir(merged_dir, max_dist, method), "MatchingOverview.png"
    )
    return batch_lock.sentinel_is_fresh(sentinel, REDO_FROM_DATE)


def get_group_lock_path(merged_dir):
    """
    Lock file marking 'a run is currently processing all pending max_dist
    sweep points for this group', so multiple machines pointed at the same
    BASE_INPUT/BASE_OUTPUT can split work across groups without
    double-processing one. See batch_lock.py. Named distinctly from the other
    batch scripts' locks so they can all run on the same group concurrently.
    """
    subfolder = os.path.relpath(os.path.dirname(merged_dir), base_batch.BASE_INPUT)
    return os.path.join(
        base_batch.BASE_OUTPUT, subfolder, ".processing_maxdist_sweep.lock"
    )


# ── entry point ───────────────────────────────────────────────────────────────


def parse_args():
    parser = argparse.ArgumentParser(
        description="Run DeepUnitMatch and UMPy on the merged dataset across a max_dist sweep."
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

    print(f"max_dist sweep values: {[format_dist(v) for v in MAX_DIST_VALUES]}")
    print(
        f"neighbour_dist per sweep point: "
        f"{[neighbour_dist_for(v) for v in MAX_DIST_VALUES]} "
        f"(default {DEFAULT_NEIGHBOUR_DIST})\n"
    )

    print(f"Scanning for merged-data groups under:\n  {base_batch.BASE_INPUT}\n")
    groups = base_batch.find_merged_groups()
    if not groups:
        print("No merged-data groups found.")
        return
    print(f"Found {len(groups)} group(s).\n")

    model = None  # DUM's default trained model, loaded once and reused across all sweep points

    for i, merged_dir in enumerate(groups):
        print(f"\n[{i + 1}/{len(groups)}] {merged_dir}")

        pending = [
            (max_dist, method)
            for max_dist in MAX_DIST_VALUES
            for method in ("DUM", "UMPy")
            if not maxdist_results_exist(merged_dir, max_dist, method)
        ]
        if not pending:
            print("  Skipping all sweep points (results exist and are fresh).")
            continue

        lock_path = get_group_lock_path(merged_dir)
        with batch_lock.try_lock(lock_path) as acquired:
            if not acquired:
                print(f"  Skipping (already being processed by another run): {lock_path}")
                continue

            # re-check now that we hold the lock: another machine may have
            # finished this group while we were scanning/waiting for the lock
            pending = [
                (max_dist, method)
                for max_dist in MAX_DIST_VALUES
                for method in ("DUM", "UMPy")
                if not maxdist_results_exist(merged_dir, max_dist, method)
            ]
            if not pending:
                print("  Skipping all sweep points (completed by another run).")
                continue

            sess = base_batch._prepare_session(merged_dir)
            if sess is None:
                continue

            if model is None and any(method == "DUM" for _, method in pending):
                print("Loading DeepUnitMatch model …")
                model = test.load_trained_model(device=base_batch.DEVICE)

            for max_dist, method in pending:
                save_dir = get_maxdist_save_dir(merged_dir, max_dist, method)
                label = f"{method}_maxdist={format_dist(max_dist)}"
                sess_mod = sess_with_max_dist(sess, max_dist)
                try:
                    if method == "DUM":
                        base_batch.run_deep_unit_match_core(
                            sess_mod, save_dir, model, label=label
                        )
                    else:
                        base_batch.run_umpy_core(sess_mod, save_dir, label=label)
                except Exception as e:
                    print(f"  {label} FAILED: {e}")
                    traceback.print_exc()

    print("\nAll done.")


if __name__ == "__main__":
    main()
