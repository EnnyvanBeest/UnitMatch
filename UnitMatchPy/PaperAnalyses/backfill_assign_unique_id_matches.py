# Backfill script: for every completed DeepUnitMatch/UMPy output under
# BASE_OUTPUT, writes a sibling "<Model>_AssignUniqueID" folder whose
# final_matches come from AssignUniqueID's "intermediate" grouping
# (unique_id, the same UID already saved as the "UID 1"/"UID 2" columns in
# MatchTable.csv) instead of the raw threshold/directional-filter matches
# DUM/UMPy's own AUC_summary.json is currently based on.
#
# Why this matters for a fair comparison: EMD's and DANT's final_matches are
# already "clique-level" -- every pair within one final assigned group counts
# as a match (DANT's ClusterMatrix.npy marks ALL pairs within a post-curation
# cluster, not just the edges that were individually thresholded; EMD's own
# algorithm similarly outputs a complete partitioning). DUM/UMPy's current
# final_matches, by contrast, are "edge-level" -- test.directional_filter_matrix
# (DUM) / a raw probability threshold (UMPy) only mark pairs that were
# individually thresholded, with no requirement that the resulting matches be
# globally self-consistent, and assign_unique_id() is computed during the
# batch run (for the UID columns saved to MatchTable.csv) but its output was
# never used to define final_matches for AUC purposes. That's comparing
# different processing stages, not the same deliverable.
#
# assign_unique_id() adds real correctness checks beyond directional
# filtering (both-CV agreement, excluding same-session merges that would
# create excessive ISI refractory violations) and produces three grouping
# strictness levels -- liberal (any transitive chain), intermediate (chain
# merging but requiring full agreement among adjacent-session group
# members), conservative (a pair only merges if it agrees with every unit
# already in the other's group, across all sessions). "Intermediate" (the
# "UID 1"/"UID 2" columns) is UnitMatch's own default/recommended level, and
# is what this script uses. ("Liberal" would be the wrong choice here even
# though it sounds permissive-friendly: its transitive chaining can mark
# pairs as matched that were never individually verified, which is *more*
# permissive than the current raw threshold, not a cleanup of it.)
#
# Since assign_unique_id() is already computed and saved into every
# historical MatchTable.csv (as the "UID 1"/"UID 2" columns, alongside the
# already-computed functional-score columns), this needs no rerun of DUM/UMPy
# inference at all -- it's pure postprocessing of already-saved output,
# following the same pattern as backfill_natim_scores.py.

import os
import sys
import json
import pickle
import datetime
import traceback

import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")  # non-interactive backend for batch runs
import matplotlib.pyplot as plt

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

SOURCE_MODELS = {
    "DeepUnitMatch": get_save_dir,
    "UMPy": get_umpy_save_dir,
}
SUFFIX = "_AssignUniqueID"

# A group/source is skipped once its output exists and is at least this new
# (see batch_lock.sentinel_is_fresh). Set to None for plain "skip if present".
REDO_FROM_DATE = None


def get_uid_dir(merged_dir, source_model):
    """Output dir for one source model's AssignUniqueID-backfilled matches --
    a sibling of DeepUnitMatch/UMPy/EMD/DANT/..., named so plot_auc_summary.py
    picks it up as its own model automatically."""
    src_dir = SOURCE_MODELS[source_model](merged_dir)
    return src_dir + SUFFIX


def uid_results_exist(merged_dir, source_model):
    sentinel = os.path.join(get_uid_dir(merged_dir, source_model), "AUC_summary.json")
    return batch_lock.sentinel_is_fresh(sentinel, REDO_FROM_DATE)


def get_backfill_lock_path(merged_dir, source_model):
    """One lock per (group, source model) so multiple machines pointed at the
    same BASE_OUTPUT can split this backfill across groups without duplicating
    work on the same folder. Named distinctly from the natim backfill's lock
    so the two backfills can run concurrently."""
    src_dir = SOURCE_MODELS[source_model](merged_dir)
    return os.path.join(src_dir, ".processing_assign_unique_id_backfill.lock")


def process_source(merged_dir, source_model):
    src_dir = SOURCE_MODELS[source_model](merged_dir)
    sentinel = os.path.join(src_dir, "MatchingOverview.png")
    if not os.path.isfile(sentinel):
        print(f"  [{source_model}] pipeline never completed here, skipping: {src_dir}")
        return

    match_table_path = os.path.join(src_dir, "MatchTable.csv")
    param_path = os.path.join(src_dir, "UMparam.pickle")
    clus_info_path = os.path.join(src_dir, "ClusInfo.pickle")
    if not (os.path.isfile(match_table_path) and os.path.isfile(param_path) and os.path.isfile(clus_info_path)):
        print(f"  [{source_model}] missing MatchTable.csv/UMparam.pickle/ClusInfo.pickle, skipping: {src_dir}")
        return

    with open(param_path, "rb") as f:
        param = pickle.load(f)
    with open(clus_info_path, "rb") as f:
        clus_info = pickle.load(f)
    session_id = clus_info["session_id"]
    n_units = param["n_units"]

    df = pd.read_csv(match_table_path)
    if len(df) != n_units * n_units:
        print(f"  [{source_model}] MatchTable.csv row count doesn't match n_units**2, skipping: {src_dir}")
        return
    if "UID 1" not in df.columns or "UID 2" not in df.columns:
        print(
            f"  [{source_model}] MatchTable.csv has no UID columns (assign_unique_id "
            f"wasn't recorded for this run), skipping: {src_dir}"
        )
        return

    # Same flattened (n_units, n_units) row order make_match_table() built
    # every column with (np.reshape(..., n_units*n_units)) -- reading any
    # column back this way keeps everything aligned without needing a
    # separate reconstruction of session pairing.
    uid1 = df["UID 1"].to_numpy().reshape(n_units, n_units)
    uid2 = df["UID 2"].to_numpy().reshape(n_units, n_units)
    rec_ses1 = df["RecSes 1"].to_numpy().reshape(n_units, n_units)
    rec_ses2 = df["RecSes 2"].to_numpy().reshape(n_units, n_units)

    across_session = rec_ses1 != rec_ses2
    final_matches = (uid1 == uid2) & across_session
    n_matches = int(np.sum(final_matches)) // 2
    print(f"  [{source_model}] AssignUniqueID (intermediate) matches: {n_matches}")

    score_columns = [
        "ISI_correlations",
        "ISI_KL_divergence",
        "ISI_wasserstein_distance",
        "refpop_correlations",
        "FR_diff",
        "ISI_CV_diff",
        "natim_correlations",
    ]
    functional_scores = {
        col: df[col].to_numpy().reshape(n_units, n_units) for col in score_columns if col in df.columns
    }

    out_dir = get_uid_dir(merged_dir, source_model)
    os.makedirs(out_dir, exist_ok=True)

    auc_summary = test.auc_summary_from_functional_scores(functional_scores, final_matches, session_id)
    su.save_auc_summary(out_dir, auc_summary)

    fig, ax = plt.subplots(figsize=(5, 5))
    im = ax.imshow(final_matches, cmap="viridis", aspect="auto")
    ax.set_title(f"{source_model} + AssignUniqueID (intermediate) matches (n={n_matches})")
    ax.set_xlabel("Unit")
    ax.set_ylabel("Unit")
    fig.colorbar(im, ax=ax)
    fig.tight_layout()
    fig.savefig(os.path.join(out_dir, "MatchingOverview.png"), dpi=150)
    plt.close(fig)

    print(f"  [{source_model}] Wrote {out_dir}")


def process_group(merged_dir):
    pending = [m for m in SOURCE_MODELS if not uid_results_exist(merged_dir, m)]
    if not pending:
        print("  Skipping (all sources already backfilled and fresh).")
        return

    for source_model in pending:
        lock_path = get_backfill_lock_path(merged_dir, source_model)
        with batch_lock.try_lock(lock_path) as acquired:
            if not acquired:
                print(f"  [{source_model}] already being backfilled by another run, skipping: {lock_path}")
                continue

            # re-check now that we hold the lock: another machine may have
            # finished this one while we were waiting for the lock
            if uid_results_exist(merged_dir, source_model):
                print(f"  [{source_model}] completed by another run, skipping.")
                continue

            try:
                process_source(merged_dir, source_model)
            except Exception as e:
                print(f"  [{source_model}] FAILED: {e}")
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
        process_group(merged_dir)

    print("\nAll done.")


if __name__ == "__main__":
    main()
