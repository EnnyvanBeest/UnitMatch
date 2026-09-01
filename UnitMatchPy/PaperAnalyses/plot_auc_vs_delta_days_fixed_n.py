# Companion to plot_auc_vs_delta_days.py: that script pools every across-
# session pair a model itself decided was a "Match" into its ΔDay bins, so a
# model that simply calls more pairs matches (a more lenient threshold) gets
# a different, not-directly-comparable AUC from a model that calls fewer --
# AUC is sensitive to how many/which positives are in the pool, independent
# of whether the *ranking* of functional-score similarity is actually
# better. This script controls for that the same way figure2.ipynb/
# figure-model-comparison.ipynb already do for two models at a time
# (fast_testing.py's get_matches_1model(fixed_n=...)/"N_set_by_UM"),
# generalised here to every rank-eligible model (see RANK_ELIGIBLE_MODELS
# below) and applied per ΔDay bin like plot_auc_vs_delta_days.py: for every
# session pair, every model's positive set is *re-derived* as the top N
# candidate pairs by that model's own "UM Probabilities" score, where N is
# REFERENCE_MODEL's own conflict-resolved match count for that same pair (an
# unordered-unit-pair count, not identity -- see _build_reference_n_lookup()).
# This mirrors fast_testing.py's get_matches_1model(fixed_n=...) closely: the
# ranking operates over the *full* across-session candidate pool for that
# pair, not just the pipeline's own already-decided matches, so a model can
# both lose matches it had decided on (if they rank low) and gain ones it
# hadn't (if they rank in the top N despite being below that model's own
# decision threshold) -- unlike an earlier, more conservative version of this
# script that only ever discarded a model's own decided matches and never
# added new ones. Per-pair scores are averaged across both directions before
# ranking (mirrors the notebooks' avg_across_directions()), one-to-many
# conflicts are resolved on that averaged score before ranking (mirrors
# remove_conflicts2(), and runs in the same order fast_testing uses: resolve
# conflicts, then take the top N), and ties are broken by a random draw
# (mirrors get_ordered_matches()'s random_tiebreaker column) -- see
# _rank_based_fixed_n_matches()'s docstring for what's still not reproduced
# (directional_filter's require-both-directions-above-threshold check).
#
# "UM Probabilities" is only a real, continuously-rankable score for
# DeepUnitMatch/UMPy-family models -- EMD/DANT/DANT_no_functional save
# final_matches.astype(float) into that same column (confirmed in
# run_emd_batch_onMerged.py), i.e. a bare 0/1 decision with nothing to rank
# among the rejected candidates. Padding those models up to N would mean
# picking additional "matches" uniformly at random from an undifferentiated
# pile of non-matches -- not a meaningful re-derivation, and not comparable
# to what's done for the ranked models -- so RANK_ELIGIBLE_MODELS excludes
# them entirely rather than represent them with a different, weaker method.
#
# Even restricted to ranked models, REFERENCE_MODEL itself is structurally
# advantaged: N for a pair is defined as *its* own decided count, so its own
# top-N-by-probability selection for that pair will typically reproduce (or
# closely approximate) its own real decision, while every other model's
# selection is actively reshaped to hit that same N. This is an inherent
# property of pinning N to one fixed model rather than, say, the per-pair
# minimum across included models -- known and not corrected for here.
#
# Reuses plot_auc_vs_delta_days.py for everything except the re-derivation
# step itself (session-date lookup, binning, AUC/summary/plotting functions
# are all identical once Matches has been replaced) via that module's
# matches_transform hook on collect_binned_pairs().
#
# Output goes to BASE_OUTPUT/auc_vs_delta_days_fixed_n_report, a sibling of
# plot_auc_vs_delta_days.py's own auc_vs_delta_days_report.

import os
import sys

import numpy as np
import pandas as pd

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _HERE)
sys.path.insert(0, os.path.dirname(_HERE))

import plot_auc_vs_delta_days as dvd
import plot_auc_summary as auc_summary_mod

BASE_OUTPUT = dvd.BASE_OUTPUT
OUTPUT_DIR = os.path.join(BASE_OUTPUT, "auc_vs_delta_days_fixed_n_report")

MIN_MATCHES_TO_INCLUDE = dvd.MIN_MATCHES_TO_INCLUDE
MIN_MICE_PER_BIN = dvd.MIN_MICE_PER_BIN

# Only models with a real per-pair "UM Probabilities" score can be re-ranked
# to hit a fixed N (see module docstring) -- the DeepUnitMatch/UMPy-family
# entries of plot_auc_vs_delta_days.MODELS_TO_INCLUDE. EMD, DANT and
# DANT_no_functional are deliberately left out.
RANK_ELIGIBLE_MODELS = {
    "DeepUnitMatch_AssignUniqueID_Conservative",
    "DeepUnitMatch_AssignUniqueID",
    "UMPy_AssignUniqueID_Conservative",
    "UMPy_AssignUniqueID",
}

# Every model's per-session-pair positive count is re-derived relative to
# this model's own match count for that pair (see module docstring). Chosen
# over bare "UMPy" (the notebooks' hardcoded default) because bare UMPy isn't
# itself one of MODELS_TO_INCLUDE here -- its UID-clustering variants are --
# and the conservative variant is the closest analogue to a single, stable
# "how many matches did UMPy actually commit to" count.
REFERENCE_MODEL = "UMPy_AssignUniqueID_Conservative"

# Which unit pairs survive the top-N selection is tie-broken randomly per
# (dataset, session pair); fixed so re-running this script reproduces the
# same figures.
RANDOM_SEED = 0

# The two RANK_ELIGIBLE_MODELS compared by the compact key-scores AUC-diff
# summary (dvd.plot_model_diff_summary(), restricted to
# dvd.KEY_SCORES_FOR_DIFF_SUMMARY) -- the Conservative UID variants are left
# out of this particular comparison since it's meant as a quick DUM-vs-UMPy
# headline, not a full model sweep (see main()'s own len(models) == 2 guard
# for the fuller all-scores version this module doesn't otherwise produce).
DIFF_SUMMARY_MODEL_A = "DeepUnitMatch_AssignUniqueID"
DIFF_SUMMARY_MODEL_B = "UMPy_AssignUniqueID"


def _build_reference_n_lookup(session_date_lookup, base_output=BASE_OUTPUT, reference_model=REFERENCE_MODEL):
    """
    {dataset: {(RecSes1, RecSes2) with RecSes1 < RecSes2: N}} -- N = how many
    unordered unit pairs reference_model committed to as a match for that
    session pair, in that dataset, AFTER resolving one-to-many conflicts
    (dvd._resolve_match_conflicts(), ranked by "UM Probabilities" same as
    everywhere else) -- not the raw pipeline "Matches" column's count, which
    has no uniqueness constraint at write time (see
    dvd._resolve_match_conflicts()'s docstring) and so can overcount relative
    to what reference_model actually settled on for a unit. This mirrors
    fast_testing.py's own reference N: it comes from a prior
    get_matches_1model(fixed_n=False) run, whose remove_conflicts2() always
    resolves conflicts before N = len(match_indices) is recorded.

    Datasets/pairs missing here get treated as N=0 by
    _rank_based_fixed_n_matches() (every other model's matches for that pair
    are then dropped entirely) -- reference_model not having found anything
    for a pair is itself information: no other model's matches for that pair
    are better-supported than reference_model's null result, so this
    comparison shouldn't count them either.

    Same discovery walk as plot_auc_vs_delta_days.collect_binned_pairs(), but
    restricted to reference_model's own folders, and only reading the columns
    needed to count matches (no functional scores).

    Also excludes pairs touching an unresolved within-session oversplit
    residual (dvd.exclude_within_session_duplicates's own default in
    collect_binned_pairs()) -- otherwise N here could count a match that
    collect_binned_pairs() itself will go on to drop for every model,
    including reference_model's own re-derived selection.
    """
    variant = dvd._uid_variant_source_model(reference_model)

    lookup = {}
    n_found = 0
    for root, dirs, files in os.walk(base_output, topdown=True):
        dirs[:] = [d for d in dirs if not auc_summary_mod._should_skip_dir(d)]
        if "AUC_summary.json" not in files:
            continue
        dirs[:] = []

        if os.path.basename(root) != reference_model:
            continue

        dataset = os.path.relpath(os.path.dirname(root), base_output).replace("\\", "/")
        if not session_date_lookup.get(dataset):
            continue

        match_table_dir = os.path.join(os.path.dirname(root), variant[0]) if variant is not None else root
        match_table_path = os.path.join(match_table_dir, "MatchTable.csv")
        if not os.path.isfile(match_table_path):
            continue

        header = pd.read_csv(match_table_path, nrows=0).columns
        usecols = ["ID1", "ID2", "RecSes 1", "RecSes 2", "Matches"]
        dtype_map = {"ID1": "int32", "ID2": "int32", "RecSes 1": "int32", "RecSes 2": "int32", "Matches": "int8"}
        if dvd.MATCH_CONFLICT_RANK_COL in header:
            usecols.append(dvd.MATCH_CONFLICT_RANK_COL)
            dtype_map[dvd.MATCH_CONFLICT_RANK_COL] = "float32"
        if variant is not None:
            _, (uid_col_a, uid_col_b) = variant
            if uid_col_a not in header or uid_col_b not in header:
                print(f"  {match_table_path} missing '{uid_col_a}'/'{uid_col_b}', skipping reference lookup for {dataset}.")
                continue
            usecols += [uid_col_a, uid_col_b]
            dtype_map[uid_col_a] = dtype_map[uid_col_b] = "int32"

        df = pd.read_csv(match_table_path, usecols=usecols, dtype=dtype_map)
        if variant is not None:
            _, (uid_col_a, uid_col_b) = variant
            df["Matches"] = (df[uid_col_a] == df[uid_col_b]).astype(np.int8)

        df = df.assign(Matches=dvd._resolve_match_conflicts(df))

        pos = df[(df["RecSes 1"] < df["RecSes 2"]) & (df["Matches"] == 1)]
        dup_flagged = dvd._build_within_session_duplicate_lookup(base_output, {dataset}).get(dataset, set())
        pos, _ = dvd._drop_within_session_duplicate_pairs(pos, dup_flagged)
        counts = pos.groupby(["RecSes 1", "RecSes 2"]).size()
        lookup[dataset] = {(int(r1), int(r2)): int(n) for (r1, r2), n in counts.items()}
        n_found += 1

    print(f"Built reference-N lookup ({reference_model}) for {n_found} dataset(s).")
    return lookup


def _rank_based_fixed_n_matches(df, ref_counts_for_dataset, rng, rank_col="UM Probabilities"):
    """
    Returns a new int8 array, same length as df, replacing df["Matches"]
    entirely: for each session pair, the top ref_counts_for_dataset[(r1, r2)]
    candidate unit pairs by rank_col (0 if that pair isn't in
    ref_counts_for_dataset at all) become the new positives -- everything
    else, including pairs the pipeline itself called a match, becomes a
    negative. Mirrors fast_testing.get_matches_1model(fixed_n=...): the
    ranking runs over every across-session candidate pair for that session
    pair (df must already be restricted to RecSes 1 != RecSes 2), not just
    the ones already flagged Matches==1, so this can both add matches a
    model's own threshold rejected and drop ones it accepted.

    df must have "ID2" and rank_col columns loaded (see extra_usecols on
    collect_binned_pairs()). rank_col is averaged across both directions of
    each unordered unit pair before ranking (mirrors the notebooks'
    avg_across_directions()) since a per-direction "UM Probabilities" value
    need not itself be symmetric; ties are broken by a per-call random draw
    (mirrors get_ordered_matches()'s random_tiebreaker column).

    Not reproduced here: the notebooks additionally run directional_filter
    (require both directions individually above threshold) before ranking.
    Averaging across directions already softly penalises one-way-only
    candidates rather than excluding them outright -- a simplification
    relative to the notebooks' own fixed_n comparison. One-to-many conflicts
    (remove_conflicts2's one-partner-per-unit constraint) ARE resolved here,
    on the direction-averaged score, before top-N ranking -- same order
    fast_testing.get_matches_1model(fixed_n=...) uses (remove_conflicts2()
    runs, THEN .head(fixed_n) picks winners from what's left), rather than
    collect_binned_pairs()'s own post-hoc resolve_match_conflicts pass, which
    would otherwise run after N is already fixed (on raw per-direction
    "UM Probabilities" values, a different score source) and could shrink the
    selected set below N. main() passes resolve_match_conflicts=False to
    collect_binned_pairs() for this reason -- conflict resolution happens
    exactly once, here.
    """
    r1 = df["RecSes 1"].to_numpy()
    r2 = df["RecSes 2"].to_numpy()
    id1 = df["ID1"].to_numpy()
    id2 = df["ID2"].to_numpy()
    score = df[rank_col].to_numpy()

    lo = r1 < r2
    sess_lo, sess_hi = np.where(lo, r1, r2), np.where(lo, r2, r1)
    unit_lo, unit_hi = np.where(lo, id1, id2), np.where(lo, id2, id1)

    canon = pd.DataFrame({
        "sess_lo": sess_lo, "sess_hi": sess_hi,
        "unit_lo": unit_lo, "unit_hi": unit_hi,
        "score": score,
        "row": np.arange(len(df)),
    })
    # One id per distinct unordered unit pair -- 2 rows per id (this row + its
    # mirror direction) in the normal case, 1 if only one direction is
    # present in this bin/dataset's data.
    canon["canon_id"] = canon.groupby(["sess_lo", "unit_lo", "sess_hi", "unit_hi"], sort=False).ngroup()

    pairs = canon.groupby("canon_id", sort=False).agg(
        sess_lo=("sess_lo", "first"), sess_hi=("sess_hi", "first"),
        unit_lo=("unit_lo", "first"), unit_hi=("unit_hi", "first"),
        avg_score=("score", "mean"),
    )
    # NaN scores (shouldn't normally occur -- UM Probabilities is computed for
    # every candidate pair -- but handled defensively) never win a ranking or
    # a conflict.
    pairs["avg_score"] = pairs["avg_score"].fillna(-np.inf)

    # Resolve one-to-many conflicts per session pair before ranking (see
    # docstring above): a candidate survives only if its avg_score is
    # simultaneously maximal among every other candidate sharing its "lo"
    # endpoint unit and every other candidate sharing its "hi" endpoint unit,
    # within the same session pair -- same both-endpoints-maximal rule as
    # dvd._resolve_match_conflicts()/fast_testing.remove_conflicts2(), just
    # applied to canonical (direction-averaged) candidates instead of raw
    # rows. Losing candidates are dropped from ranking entirely (rather than
    # zeroed-but-still-eligible, as remove_conflicts2 does) so a conflict can
    # never consume one of the N slots -- a minor, deliberate simplification
    # that only differs from fast_testing in the rare case where a session
    # pair has fewer than N conflict-free candidates.
    max_for_lo = pairs.groupby(["sess_lo", "sess_hi", "unit_lo"])["avg_score"].transform("max")
    max_for_hi = pairs.groupby(["sess_lo", "sess_hi", "unit_hi"])["avg_score"].transform("max")
    pairs = pairs[(pairs["avg_score"] == max_for_lo) & (pairs["avg_score"] == max_for_hi)]

    avg_score = pairs["avg_score"].to_numpy()
    tiebreak = rng.random(len(pairs))
    rows_by_canon = canon.groupby("canon_id", sort=False)["row"].apply(np.asarray)

    matches_fixed = np.zeros(len(df), dtype=np.int8)
    selected_canon_ids = []
    for (s_lo, s_hi), g in pairs.groupby(["sess_lo", "sess_hi"]):
        n_ref = ref_counts_for_dataset.get((int(s_lo), int(s_hi)), 0)
        if n_ref <= 0:
            continue
        local = g.index.to_numpy()
        pos = pairs.index.get_indexer(local)
        # Descending by avg_score, ties broken by the random draw.
        order = np.lexsort((tiebreak[pos], -avg_score[pos]))
        selected_canon_ids.append(local[order[:n_ref]])

    if selected_canon_ids:
        sel_ids = np.concatenate(selected_canon_ids)
        sel_rows = np.concatenate(rows_by_canon.loc[sel_ids].to_numpy())
        matches_fixed[sel_rows] = 1
    return matches_fixed


def main():
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    session_date_lookup = dvd.get_session_date_lookup()
    n_datasets_with_dates = sum(1 for m in session_date_lookup.values() if m)
    print(f"Have session dates for {n_datasets_with_dates}/{len(session_date_lookup)} dataset(s).\n")

    ref_n_lookup = _build_reference_n_lookup(session_date_lookup)
    rng = np.random.default_rng(RANDOM_SEED)

    def matches_transform(dataset, df):
        return _rank_based_fixed_n_matches(df, ref_n_lookup.get(dataset, {}), rng)

    with open(os.path.join(OUTPUT_DIR, "fixed_n_run_info.txt"), "w") as f:
        f.write(
            f"reference_model = {REFERENCE_MODEL}\n"
            f"rank_eligible_models = {sorted(RANK_ELIGIBLE_MODELS)}\n"
            f"excluded_models = {sorted(dvd.MODELS_TO_INCLUDE - RANK_ELIGIBLE_MODELS)} "
            f"(no continuous per-pair score to rank by -- see module docstring)\n"
            f"random_seed = {RANDOM_SEED}\n"
            f"datasets_with_reference_matches = {len(ref_n_lookup)}\n"
        )

    (
        score_bins, rate_bins, count_bins, dataset_level_rows, dataset_level_auc_rows,
        auc_long_df, dataset_bin_auc_rows, dataset_bin_rate_rows,
    ) = dvd.collect_binned_pairs(
        session_date_lookup,
        models_to_include=RANK_ELIGIBLE_MODELS,
        matches_transform=matches_transform,
        extra_usecols={"ID2": "int32", "UM Probabilities": "float32"},
        # _rank_based_fixed_n_matches() already resolves one-to-many conflicts
        # itself, before picking the top N (see its docstring) -- running
        # collect_binned_pairs()'s own post-hoc pass here too would re-decide
        # conflicts from raw per-direction "UM Probabilities" values (a
        # different score source than the direction-averaged one used for
        # ranking) after N is already fixed, and could shrink a model's
        # selected set below N.
        resolve_match_conflicts=False,
    )

    auc_df = dvd.summarise_auc(score_bins)
    auc_csv = os.path.join(OUTPUT_DIR, "auc_vs_delta_days_fixed_n.csv")
    auc_df.to_csv(auc_csv, index=False)
    print(f"Wrote {auc_csv}")

    rate_df = dvd.summarise_match_rate(rate_bins)
    rate_csv = os.path.join(OUTPUT_DIR, "match_rate_vs_delta_days_fixed_n.csv")
    rate_df.to_csv(rate_csv, index=False)
    print(f"Wrote {rate_csv}")

    count_df = dvd.summarise_dataset_counts(count_bins)
    count_csv = os.path.join(OUTPUT_DIR, "dataset_counts_vs_delta_days_fixed_n.csv")
    count_df.to_csv(count_csv, index=False)
    print(f"Wrote {count_csv}")

    models = sorted(set(auc_df["model"]) | set(rate_df["model"]) | set(count_df["model"]))
    colour_for = dvd.build_family_colours(models)

    # Model-pair statistics (overall mixed-effects test + per-bin paired
    # t-test), restricted to DIFF_SUMMARY_MODEL_A/B -- see
    # plot_auc_vs_delta_days.py main()'s own version of this same block.
    # Computed once from dataset_bin_auc_rows/dataset_bin_rate_rows and
    # reused by every plot below, pooled or mouse-averaged alike.
    stats_available = DIFF_SUMMARY_MODEL_A in models and DIFF_SUMMARY_MODEL_B in models
    if stats_available:
        print(f"Computing model-comparison stats: {DIFF_SUMMARY_MODEL_A} vs {DIFF_SUMMARY_MODEL_B}")
        rate_overall = dvd.test_overall_model_effect(dataset_bin_rate_rows, "rate", DIFF_SUMMARY_MODEL_A, DIFF_SUMMARY_MODEL_B)
        rate_bin_pvals = dvd.test_per_bin_model_effect(dataset_bin_rate_rows, "rate", DIFF_SUMMARY_MODEL_A, DIFF_SUMMARY_MODEL_B)
        all_scores_for_stats = sorted(auc_df["score"].unique())
        score_overall = {
            s: dvd.test_overall_model_effect(dataset_bin_auc_rows, "auc", DIFF_SUMMARY_MODEL_A, DIFF_SUMMARY_MODEL_B, score=s)
            for s in all_scores_for_stats
        }
        score_bin_pvals = {
            s: dvd.test_per_bin_model_effect(dataset_bin_auc_rows, "auc", DIFF_SUMMARY_MODEL_A, DIFF_SUMMARY_MODEL_B, score=s)
            for s in all_scores_for_stats
        }
    else:
        rate_overall, rate_bin_pvals = None, {}
        score_overall, score_bin_pvals = {}, {}

    for score in sorted(auc_df["score"].unique()):
        out_path = os.path.join(OUTPUT_DIR, f"summary_{score}_fixed_n.png")
        result = dvd.plot_score_summary(
            score, auc_df, rate_df, count_df, colour_for, out_path, title_suffix=" (fixed N)",
            stats_model_a=DIFF_SUMMARY_MODEL_A if stats_available else None,
            stats_model_b=DIFF_SUMMARY_MODEL_B if stats_available else None,
            rate_overall=rate_overall, rate_bin_pvals=rate_bin_pvals,
            auc_overall=score_overall.get(score), auc_bin_pvals=score_bin_pvals.get(score, {}),
        )
        if result:
            print(f"  Plotted {score} -> {result}")

    # Compact key-scores AUC-diff summary (DIFF_SUMMARY_MODEL_A vs
    # DIFF_SUMMARY_MODEL_B, dvd.KEY_SCORES_FOR_DIFF_SUMMARY only) -- the
    # fixed-N counterpart to plot_auc_vs_delta_days.py main()'s own version
    # of this same compact plot.
    if stats_available:
        key_scores_present = [s for s in dvd.KEY_SCORES_FOR_DIFF_SUMMARY if s in auc_df["score"].unique()]
        if key_scores_present:
            colour_for_score = dvd.build_qualitative_colours(auc_df["score"].unique())
            key_diff_df = dvd.compute_auc_diff(DIFF_SUMMARY_MODEL_A, DIFF_SUMMARY_MODEL_B, auc_df, scores=key_scores_present)
            key_diff_csv = os.path.join(
                OUTPUT_DIR, f"auc_diff_vs_delta_days_fixed_n_key_scores_{DIFF_SUMMARY_MODEL_A}_vs_{DIFF_SUMMARY_MODEL_B}.csv"
            )
            key_diff_df.to_csv(key_diff_csv, index=False)
            print(f"Wrote {key_diff_csv}")

            key_diff_out_path = os.path.join(
                OUTPUT_DIR, f"summary_diff_key_scores_fixed_n_{DIFF_SUMMARY_MODEL_A}_vs_{DIFF_SUMMARY_MODEL_B}.png"
            )
            result = dvd.plot_model_diff_summary(
                DIFF_SUMMARY_MODEL_A, DIFF_SUMMARY_MODEL_B, auc_df, rate_df, count_df, colour_for, colour_for_score,
                key_diff_out_path, scores=key_scores_present, title_suffix=" (key scores, fixed N)",
                rate_overall=rate_overall, rate_bin_pvals=rate_bin_pvals,
                score_overall={s: score_overall[s] for s in key_scores_present},
                score_bin_pvals={s: score_bin_pvals[s] for s in key_scores_present},
            )
            if result:
                print(f"  Plotted key-scores AUC-diff summary (fixed N) -> {result}")
    else:
        print(
            f"Skipping key-scores AUC-diff summary (fixed N): need both {DIFF_SUMMARY_MODEL_A} and "
            f"{DIFF_SUMMARY_MODEL_B} present in {sorted(models)}."
        )

    auc_mouse_df = dvd.summarise_auc_mouse_averaged(dataset_bin_auc_rows)
    auc_mouse_csv = os.path.join(OUTPUT_DIR, "auc_vs_delta_days_fixed_n_mouse_averaged.csv")
    auc_mouse_df.to_csv(auc_mouse_csv, index=False)
    print(f"Wrote {auc_mouse_csv}")

    rate_mouse_df = dvd.summarise_match_rate_mouse_averaged(dataset_bin_rate_rows)
    rate_mouse_csv = os.path.join(OUTPUT_DIR, "match_rate_vs_delta_days_fixed_n_mouse_averaged.csv")
    rate_mouse_df.to_csv(rate_mouse_csv, index=False)
    print(f"Wrote {rate_mouse_csv}")

    for score in sorted(auc_mouse_df["score"].unique()):
        out_path = os.path.join(OUTPUT_DIR, f"summary_{score}_fixed_n_mouse_averaged.png")
        result = dvd.plot_score_summary(
            score, auc_mouse_df, rate_mouse_df, count_df, colour_for, out_path,
            min_count=MIN_MICE_PER_BIN, count_col="n_mice", title_suffix=" (fixed N, mouse-averaged)",
            stats_model_a=DIFF_SUMMARY_MODEL_A if stats_available else None,
            stats_model_b=DIFF_SUMMARY_MODEL_B if stats_available else None,
            rate_overall=rate_overall, rate_bin_pvals=rate_bin_pvals,
            auc_overall=score_overall.get(score), auc_bin_pvals=score_bin_pvals.get(score, {}),
        )
        if result:
            print(f"  Plotted {score} (mouse-averaged) -> {result}")

    if any(model in models for model in dvd.CIRCULAR_SCORES):
        print(
            f"Note: excluding {sorted(dvd.SCORES_EXCLUDED_FROM_QUALITY_AVERAGE)} from every model's pooled/mouse-"
            f"balanced quality average (circular for at least one included model -- see dvd.CIRCULAR_SCORES)."
        )

    # ── quality vs quantity, same as plot_auc_vs_delta_days.py, on the fixed-N bins ──

    summary_df = dvd.compute_quality_quantity_summary(auc_df, rate_df)
    summary_csv = os.path.join(OUTPUT_DIR, "quality_vs_quantity_fixed_n.csv")
    summary_df.to_csv(summary_csv, index=False)
    print(f"Wrote {summary_csv}")
    print(summary_df.sort_values("quality", ascending=False).to_string(index=False))

    pooled_png = os.path.join(OUTPUT_DIR, "quality_vs_quantity_fixed_n_pooled.png")
    result = dvd.plot_quality_vs_quantity(
        summary_df, colour_for, pooled_png,
        title=f"Quality vs quantity, fixed N (pooled across datasets/ΔDay, N capped to {REFERENCE_MODEL})",
    )
    if result:
        print(f"  Plotted pooled quality-vs-quantity -> {result}")

    trajectories_png = os.path.join(OUTPUT_DIR, "quality_vs_quantity_fixed_n_by_bin.png")
    result = dvd.plot_quality_vs_quantity_trajectories(auc_df, rate_df, colour_for, trajectories_png)
    if result:
        print(f"  Plotted quality-vs-quantity by ΔDay bin -> {result}")

    per_score_png = os.path.join(OUTPUT_DIR, "quality_vs_quantity_fixed_n_per_score.png")
    result = dvd.plot_quality_vs_quantity_per_score(auc_df, rate_df, colour_for, per_score_png)
    if result:
        print(f"  Plotted quality-vs-quantity per score -> {result}")

    # auc_long_df (straight from each model's own AUC_summary.json) is only
    # used here to decide dataset inclusion -- that threshold is about
    # whether a dataset had enough raw trackable matches to trust at all, a
    # property of the pipeline's own original match counts regardless of any
    # N-fixing choice made afterwards, and is the same threshold
    # collect_binned_pairs() already applied internally to build score_bins/
    # rate_bins/count_bins above. It must NOT be used as the AUC ("quality")
    # source below -- those values were never touched by matches_transform,
    # so mixing them with the (transformed) P_track values in dataset_level_rows
    # would silently compare quantity computed one way against quality
    # computed another. dataset_level_auc_rows is the fixed-N-consistent
    # per-score AUC source instead (see its docstring in collect_binned_pairs()).
    passing_datasets = auc_summary_mod.get_passing_datasets(auc_long_df, min_matches=MIN_MATCHES_TO_INCLUDE)
    ptrack_df = pd.DataFrame(dataset_level_rows)
    functional_auc_df = pd.DataFrame(dataset_level_auc_rows)
    combined_df_long = pd.concat([functional_auc_df, ptrack_df], ignore_index=True)
    combined_df_long = combined_df_long[
        combined_df_long["model"].isin(models) & combined_df_long["dataset"].isin(passing_datasets)
    ]
    combined_csv = os.path.join(OUTPUT_DIR, "quality_quantity_fixed_n_long_by_dataset.csv")
    combined_df_long.to_csv(combined_csv, index=False)
    print(f"Wrote {combined_csv}")

    mouse_summary_df = dvd.mouse_balanced_quality_quantity(combined_df_long)
    mouse_summary_csv = os.path.join(OUTPUT_DIR, "quality_vs_quantity_fixed_n_mouse_balanced.csv")
    mouse_summary_df.to_csv(mouse_summary_csv, index=False)
    print(f"Wrote {mouse_summary_csv}")
    print(mouse_summary_df.sort_values("quality", ascending=False).to_string(index=False))

    mouse_png = os.path.join(OUTPUT_DIR, "quality_vs_quantity_fixed_n_mouse_balanced.png")
    result = dvd.plot_quality_vs_quantity(
        mouse_summary_df, colour_for, mouse_png,
        title="Quality vs quantity, fixed N (mouse-balanced: average of per-mouse averages)",
    )
    if result:
        print(f"  Plotted mouse-balanced quality-vs-quantity -> {result}")

    png_path, stats_path = dvd.test_ptrack_mouse_balanced(combined_df_long, OUTPUT_DIR)
    if png_path:
        print(f"  Plotted P_track (mouse-aware) -> {png_path}")
    print(f"Wrote P_track mixed-model summary to {stats_path}")

    # Proper dataset-level mixed-model + Holm-Bonferroni stats for every
    # functional score at fixed N -- there's no equivalent to
    # auc_summary_report/mixed_model_summaries.txt for this data otherwise,
    # since those were computed on the pipeline's own untransformed matches.
    score_png_paths, score_stats_path = dvd.test_functional_scores_mouse_balanced(
        combined_df_long, OUTPUT_DIR, filename="functional_score_mixed_model_summaries_fixed_n.txt",
    )
    for p in score_png_paths:
        print(f"  Plotted {os.path.basename(p)} (fixed N, dataset-level) -> {p}")
    print(f"Wrote fixed-N functional-score mixed-model summary to {score_stats_path}")


if __name__ == "__main__":
    main()
