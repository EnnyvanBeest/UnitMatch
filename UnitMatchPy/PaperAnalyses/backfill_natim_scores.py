# Backfill script: walks every completed DeepUnitMatch/UMPy output under
# BASE_OUTPUT and adds the natim_correlations functional score + AUC wherever
# it's missing, without re-running the rest of the pipeline (DUM/UMPy
# matching, other functional scores, etc. are left untouched).
#
# Two independent bugs -- both now fixed in test.py / run_deepunitmatch_batch_onMerged.py
# -- meant natim_correlations silently failed for every merged-pipeline run so
# far, not only ones genuinely lacking natural-image data:
#   - DUM passed merged_architecture=True to natim_correlations(), which didn't
#     accept that keyword at all (TypeError, caught by the pipeline's own
#     surrounding try/except and silently dropped).
#   - UMPy never passed merged_architecture=True, so get_natim_responses()
#     looked for trial.imageIDs.npy etc. two directory levels above each
#     session (the convention for the original, non-merged KS directories)
#     instead of directly inside it, which is where generate_merged_dataset.py
#     actually copies those files for the merged tree. That produced a
#     FileNotFoundError per session (silently NaN'd out inside
#     get_natim_responses), and the resulting all-NaN functional score then
#     made test.AUC() raise "No matches found", again swallowed by the same
#     try/except.
#
# This script recomputes natim_correlations now that both call sites are
# fixed. For groups where natural-image data is genuinely absent from the
# merged tree (e.g. generate_merged_dataset.py itself never had access to the
# raw source server hosting it when it ran), this will still legitimately
# fail -- those are reported separately from the ones actually repaired, since
# no amount of re-running this script fixes that; it needs generate_merged_dataset.py
# re-run with proper access to backfill the missing trial files first.

import os
import sys
import json
import pickle
import traceback

import numpy as np
import pandas as pd

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _HERE)
sys.path.insert(0, os.path.dirname(_HERE))
sys.path.insert(0, os.path.join(_HERE, "DeepUnitMatch"))

import batch_lock
import UnitMatchPy.save_utils as su
from DeepUnitMatch.testing import test

from run_deepunitmatch_batch_onMerged import (
    BASE_INPUT,
    BASE_OUTPUT,
    find_merged_groups,
    get_save_dir,
    get_umpy_save_dir,
)


def get_backfill_lock_path(save_dir):
    """See batch_lock.py -- lets multiple machines split this backfill pass
    across groups without duplicating the (non-trivial: PSTH binning over
    every session) natim_correlations recomputation for the same folder."""
    return os.path.join(save_dir, ".processing_natim_backfill.lock")


def backfill_one(save_dir, label):
    sentinel = os.path.join(save_dir, "MatchingOverview.png")
    if not os.path.isfile(sentinel):
        print(f"  [{label}] pipeline never completed here, skipping: {save_dir}")
        return

    auc_path = os.path.join(save_dir, "AUC_summary.json")
    auc_summary = {}
    if os.path.isfile(auc_path):
        with open(auc_path) as f:
            auc_summary = json.load(f)
    if "natim_correlations" in auc_summary:
        print(f"  [{label}] already has natim_correlations, skipping.")
        return

    match_table_path = os.path.join(save_dir, "MatchTable.csv")
    param_path = os.path.join(save_dir, "UMparam.pickle")
    clus_info_path = os.path.join(save_dir, "ClusInfo.pickle")
    if not (os.path.isfile(match_table_path) and os.path.isfile(param_path) and os.path.isfile(clus_info_path)):
        print(f"  [{label}] missing MatchTable.csv/UMparam.pickle/ClusInfo.pickle, skipping: {save_dir}")
        return

    lock_path = get_backfill_lock_path(save_dir)
    with batch_lock.try_lock(lock_path) as acquired:
        if not acquired:
            print(f"  [{label}] already being backfilled by another run, skipping: {lock_path}")
            return

        # re-check now that we hold the lock: another machine may have
        # finished this one while we were waiting for the lock
        if os.path.isfile(auc_path):
            with open(auc_path) as f:
                auc_summary = json.load(f)
        if "natim_correlations" in auc_summary:
            print(f"  [{label}] completed by another run, skipping.")
            return

        with open(param_path, "rb") as f:
            param = pickle.load(f)
        with open(clus_info_path, "rb") as f:
            clus_info = pickle.load(f)
        session_id = clus_info["session_id"]
        n_units = param["n_units"]

        df = pd.read_csv(match_table_path)
        if len(df) != n_units * n_units:
            print(f"  [{label}] MatchTable.csv row count doesn't match n_units**2, skipping: {save_dir}")
            return
        # "Matches" is the same flattened (n_units, n_units) boolean matrix
        # make_match_table() built via np.reshape(..., n_units*n_units) --
        # reading it back this way guarantees the row order lines up with
        # whatever we reshape natimcorr into below, without needing a
        # separate Matches.npy/session_id reconstruction.
        final_matches = df["Matches"].to_numpy().reshape(n_units, n_units).astype(bool)

        try:
            natimcorr = test.natim_correlations(param, merged_architecture=True)
        except Exception as e:
            print(f"  [{label}] natim_correlations still unavailable (likely no natural-image "
                  f"data copied into the merged tree for this group): {e}")
            return

        try:
            auc_natim = test.AUC(final_matches, natimcorr, session_id)
        except Exception as e:
            print(f"  [{label}] natim_correlations computed but AUC failed (likely all-NaN, "
                  f"no usable data for any session in this group): {e}")
            return

        # success -- patch both files, leaving everything else untouched
        df["natim_correlations"] = np.reshape(natimcorr, (n_units * n_units))
        df.to_csv(match_table_path, index=False)

        auc_summary["natim_correlations"] = float(auc_natim)
        su.save_auc_summary(save_dir, auc_summary)

        print(f"  [{label}] backfilled natim_correlations (AUC={auc_natim:.3f}): {save_dir}")


def main():
    print(f"Scanning for merged-data groups under:\n  {BASE_INPUT}\n")
    groups = find_merged_groups()
    if not groups:
        print("No merged-data groups found.")
        return
    print(f"Found {len(groups)} group(s).\n")

    for i, merged_dir in enumerate(groups):
        print(f"\n[{i + 1}/{len(groups)}] {merged_dir}")
        for save_dir, label in [
            (get_save_dir(merged_dir), "DeepUnitMatch"),
            (get_umpy_save_dir(merged_dir), "UMPy"),
        ]:
            try:
                backfill_one(save_dir, label)
            except Exception as e:
                print(f"  [{label}] FAILED: {e}")
                traceback.print_exc()

    print("\nAll done.")


if __name__ == "__main__":
    main()
