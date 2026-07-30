# Backfill script: for every completed DeepUnitMatch/UMPy output under
# BASE_OUTPUT, writes sibling "<Model>_AssignUniqueID" (intermediate) and
# "<Model>_AssignUniqueID_Conservative" (conservative) folders whose
# final_matches come from AssignUniqueID's grouping (the "UID 1"/"UID 2" and
# "UID Cons 1"/"UID Cons 2" columns already saved in MatchTable.csv) instead
# of the raw threshold/directional-filter matches DUM/UMPy's own
# AUC_summary.json is currently based on.
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
# already in the other's group, across all sessions). We backfill both
# "intermediate" (UnitMatch's own default/recommended level) and
# "conservative" (the strictest, full-clique-agreement level), so both ends
# of the tradeoff are available to compare. ("Liberal" is deliberately not
# included: its transitive chaining can mark pairs as matched that were never
# individually verified, which is *more* permissive than the current raw
# threshold, not a cleanup of it.)
#
# Empirically (checked directly against a real dataset before writing this
# comment): both intermediate AND conservative grouping *increase* the total
# matched-pair count relative to the current raw threshold, not decrease it
# -- transitive chaining adds pairs that were never individually verified
# faster than the extra correctness checks remove them. Average functional-
# score AUC drops slightly as a result. So "cleaned up" here means "more
# internally consistent", not "stricter" -- worth keeping in mind when
# interpreting the comparison, this was not the a priori assumption.
#
# Since assign_unique_id() is already computed and saved into every
# historical MatchTable.csv (as the UID columns, alongside the already-
# computed functional-score columns), this needs no rerun of DUM/UMPy
# inference at all -- it's pure postprocessing of already-saved output,
# following the same pattern as backfill_natim_scores.py.

import os
import sys
import pickle
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

# suffix -> (UID column name for unit A, UID column name for unit B)
UID_VARIANTS = {
    "AssignUniqueID": ("UID 1", "UID 2"),
    "AssignUniqueID_Conservative": ("UID Cons 1", "UID Cons 2"),
}

# A group/source/variant is skipped once its output exists and is at least
# this new (see batch_lock.sentinel_is_fresh). None = plain "skip if present".
REDO_FROM_DATE = None


def get_uid_dir(merged_dir, source_model, variant_suffix):
    """Output dir for one source model's AssignUniqueID-backfilled matches --
    a sibling of DeepUnitMatch/UMPy/EMD/DANT/..., named so plot_auc_summary.py
    picks it up as its own model automatically."""
    src_dir = SOURCE_MODELS[source_model](merged_dir)
    return f"{src_dir}_{variant_suffix}"


def uid_results_exist(merged_dir, source_model, variant_suffix):
    sentinel = os.path.join(get_uid_dir(merged_dir, source_model, variant_suffix), "AUC_summary.json")
    return batch_lock.sentinel_is_fresh(sentinel, REDO_FROM_DATE)


def get_backfill_lock_path(merged_dir, source_model):
    """One lock per (group, source model), covering both UID variants (they
    share the same expensive MatchTable.csv load), so multiple machines
    pointed at the same BASE_OUTPUT can split this backfill across groups
    without duplicating work. Named distinctly from the natim backfill's
    lock so the two backfills can run concurrently."""
    src_dir = SOURCE_MODELS[source_model](merged_dir)
    return os.path.join(src_dir, ".processing_assign_unique_id_backfill.lock")


def compute_and_save_variant(merged_dir, source_model, variant_suffix, df, n_units, session_id):
    col_a, col_b = UID_VARIANTS[variant_suffix]
    if col_a not in df.columns or col_b not in df.columns:
        print(f"  [{source_model}/{variant_suffix}] MatchTable.csv missing '{col_a}'/'{col_b}', skipping.")
        return

    # Same flattened (n_units, n_units) row order make_match_table() built
    # every column with (np.reshape(..., n_units*n_units)) -- reading any
    # column back this way keeps everything aligned without needing a
    # separate reconstruction of session pairing.
    uid_a = df[col_a].to_numpy().reshape(n_units, n_units)
    uid_b = df[col_b].to_numpy().reshape(n_units, n_units)
    rec_ses1 = df["RecSes 1"].to_numpy().reshape(n_units, n_units)
    rec_ses2 = df["RecSes 2"].to_numpy().reshape(n_units, n_units)

    across_session = rec_ses1 != rec_ses2
    final_matches = (uid_a == uid_b) & across_session
    n_matches = int(np.sum(final_matches)) // 2
    print(f"  [{source_model}/{variant_suffix}] matches: {n_matches}")

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

    out_dir = get_uid_dir(merged_dir, source_model, variant_suffix)
    os.makedirs(out_dir, exist_ok=True)

    auc_summary = test.auc_summary_from_functional_scores(functional_scores, final_matches, session_id)
    su.save_auc_summary(out_dir, auc_summary)

    fig, ax = plt.subplots(figsize=(5, 5))
    im = ax.imshow(final_matches, cmap="viridis", aspect="auto")
    ax.set_title(f"{source_model} + {variant_suffix} matches (n={n_matches})")
    ax.set_xlabel("Unit")
    ax.set_ylabel("Unit")
    fig.colorbar(im, ax=ax)
    fig.tight_layout()
    fig.savefig(os.path.join(out_dir, "MatchingOverview.png"), dpi=150)
    plt.close(fig)

    print(f"  [{source_model}/{variant_suffix}] Wrote {out_dir}")


def process_source(merged_dir, source_model, pending_variants):
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

    for variant_suffix in pending_variants:
        try:
            compute_and_save_variant(merged_dir, source_model, variant_suffix, df, n_units, session_id)
        except Exception as e:
            print(f"  [{source_model}/{variant_suffix}] FAILED: {e}")
            traceback.print_exc()


def process_group(merged_dir):
    pending_by_source = {
        source_model: [
            v for v in UID_VARIANTS if not uid_results_exist(merged_dir, source_model, v)
        ]
        for source_model in SOURCE_MODELS
    }
    pending_by_source = {k: v for k, v in pending_by_source.items() if v}
    if not pending_by_source:
        print("  Skipping (all sources/variants already backfilled and fresh).")
        return

    for source_model, pending_variants in pending_by_source.items():
        lock_path = get_backfill_lock_path(merged_dir, source_model)
        with batch_lock.try_lock(lock_path) as acquired:
            if not acquired:
                print(f"  [{source_model}] already being backfilled by another run, skipping: {lock_path}")
                continue

            # re-check now that we hold the lock: another machine may have
            # finished some/all pending variants while we were waiting
            pending_variants = [v for v in pending_variants if not uid_results_exist(merged_dir, source_model, v)]
            if not pending_variants:
                print(f"  [{source_model}] completed by another run, skipping.")
                continue

            process_source(merged_dir, source_model, pending_variants)


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
