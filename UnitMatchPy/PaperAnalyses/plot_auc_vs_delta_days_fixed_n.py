# Companion to plot_auc_vs_delta_days.py: that script compares every match each
# model's tracking found, so a model that tracks more units gets a different,
# not directly comparable AUC from one that tracks fewer -- AUC is sensitive to
# how many positives are in the pool, whether or not the matches are better.
# This script controls for that and stays on the tracking results.
#
# N per session pair is the number of DeepUnitMatch UID matches, and
# DeepUnitMatch keeps exactly those. UMPy gets the same N by lowering the
# threshold of its tracking: assign_unique_id() is rerun on UMPy's saved
# probabilities at RETRACK_THRESHOLDS, and each session pair uses the highest
# threshold at which UMPy reaches N, keeping its N most likely matches there.
#
# This is fast_testing.py's fixed N with tracking added: it takes the N most
# likely candidates of each session pair, which amounts to lowering the
# threshold per session pair until there are N. UMPy's matches thus come from
# different thresholds for different session pairs, not from one clustering.
# Session pairs that stay below N at the lowest threshold keep the matches they
# have (see retracking_thresholds_fixed_n.csv and fixed_n_run_info.txt).
#
# "Most likely" is the model's "UM Probabilities", averaged across both
# directions of the unit pair (as fast_testing's avg_across_directions), ties
# at random. One-to-many conflicts are resolved first, as in
# plot_auc_vs_delta_days.py, so each unit keeps at most one partner per session.
#
# Reuses plot_auc_vs_delta_days.py for everything except this selection, via
# its matches_transform hook on collect_binned_pairs().
#
# Output goes to BASE_OUTPUT/auc_vs_delta_days_fixed_n_report, a sibling of
# plot_auc_vs_delta_days.py's own auc_vs_delta_days_report.

import contextlib
import io
import os
import pickle
import sys

import numpy as np
import pandas as pd

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _HERE)
sys.path.insert(0, os.path.dirname(_HERE))

import plot_auc_vs_delta_days as dvd
import plot_auc_summary as auc_summary_mod
import UnitMatchPy.assign_unique_id as aid

BASE_OUTPUT = dvd.BASE_OUTPUT
OUTPUT_DIR = os.path.join(BASE_OUTPUT, "auc_vs_delta_days_fixed_n_report")

MIN_MATCHES_TO_INCLUDE = dvd.MIN_MATCHES_TO_INCLUDE
MIN_MICE_PER_BIN = dvd.MIN_MICE_PER_BIN

DIFF_SUMMARY_MODEL_A = "DeepUnitMatch_AssignUniqueID"
DIFF_SUMMARY_MODEL_B = "UMPy_AssignUniqueID"
COMPARED_MODELS = {DIFF_SUMMARY_MODEL_A, DIFF_SUMMARY_MODEL_B}

N_MODEL = DIFF_SUMMARY_MODEL_A  # its UID matches set N
RETRACKED_MODEL = DIFF_SUMMARY_MODEL_B  # its tracking threshold is lowered to reach N

# Thresholds tried below the model's own match_threshold, highest first.
RETRACK_THRESHOLDS = (0.4, 0.3, 0.2, 0.15, 0.1, 0.07, 0.05, 0.03, 0.02, 0.01)

# The rerun of the tracking is by far the slowest step. Delete this file to redo
# it, e.g. after the matching itself was rerun; a different RETRACK_THRESHOLDS
# or RANDOM_SEED is noticed and redoes it automatically.
RETRACK_CACHE_PATH = os.path.join(OUTPUT_DIR, "retracking_cache_fixed_n.pickle")

# Ties in the selection are broken at random; fixed so re-running this script
# reproduces the same figures.
RANDOM_SEED = 0

_flagged_units = {}  # dataset -> within-session duplicates, slow to read and the same for every model


def _read_match_table(dataset, model, base_output=BASE_OUTPUT):
    """
    The source MatchTable.csv of a "<Source>_AssignUniqueID" model with "Matches"
    its UID matches, with one row per ordered unit pair and the default index,
    as collect_binned_pairs() reads it. None if the model has no output here.
    """
    source_model, (uid_col_a, uid_col_b) = dvd._uid_variant_source_model(model)
    dataset_dir = os.path.join(base_output, *dataset.split("/"))
    path = os.path.join(dataset_dir, source_model, "MatchTable.csv")
    if not os.path.isfile(path) or not os.path.isfile(os.path.join(dataset_dir, model, "AUC_summary.json")):
        return None
    columns = {"ID1": "int32", "ID2": "int32", "RecSes 1": "int32", "RecSes 2": "int32",
               uid_col_a: "int32", uid_col_b: "int32", dvd.MATCH_CONFLICT_RANK_COL: "float32"}
    df = pd.read_csv(path, usecols=list(columns), dtype=columns)
    df["Matches"] = (df[uid_col_a] == df[uid_col_b]).astype(np.int8)
    return df


def _across_sessions(df, dataset, base_output=BASE_OUTPUT):
    """The across-session rows of df without within-session duplicates, as collect_binned_pairs() uses them."""
    if dataset not in _flagged_units:
        _flagged_units[dataset] = dvd._build_within_session_duplicate_lookup(base_output, {dataset}).get(dataset, set())
    df, _ = dvd._drop_within_session_duplicate_pairs(df[df["RecSes 1"] != df["RecSes 2"]], _flagged_units[dataset])
    return df


def _build_n_lookup(session_date_lookup, base_output=BASE_OUTPUT):
    """
    {dataset: {(RecSes1, RecSes2) with RecSes1 < RecSes2: N}}: number of
    conflict-resolved UID matches of N_MODEL. Session pairs without a match are
    absent, i.e. N = 0.
    """
    lookup = {}
    for dataset, dates in session_date_lookup.items():
        if not dates:
            continue
        df = _read_match_table(dataset, N_MODEL, base_output)
        if df is None:
            continue
        df = _across_sessions(df, dataset, base_output)
        df = df.assign(Matches=dvd._resolve_match_conflicts(df))
        n_per_pair = df[(df["RecSes 1"] < df["RecSes 2"]) & (df["Matches"] == 1)].groupby(["RecSes 1", "RecSes 2"]).size()
        lookup[dataset] = {(int(r1), int(r2)): int(n) for (r1, r2), n in n_per_pair.items()}
    print(f"Built N per session pair for {len(lookup)} dataset(s).")
    return lookup


def _best_n_uid_matches(df, n_for_dataset, rng, rank_col=dvd.MATCH_CONFLICT_RANK_COL):
    """
    Returns a new int8 array replacing df["Matches"] (df restricted to RecSes 1 !=
    RecSes 2, "Matches" already UID matches): for each session pair, the
    n_for_dataset[(r1, r2)] most likely conflict-resolved UID matches, in both
    directions. Everything else becomes a non-match, so a selected pair counts
    the same for (i, j) and (j, i).

    A unit pair is a candidate if it survived conflict resolution in the
    forward direction (RecSes 1 < RecSes 2), the direction N is counted in.
    Ranked by rank_col averaged over both directions, ties at random.
    """
    r1, r2 = df["RecSes 1"].to_numpy(), df["RecSes 2"].to_numpy()
    id1, id2 = df["ID1"].to_numpy(), df["ID2"].to_numpy()
    resolved = dvd._resolve_match_conflicts(df).astype(bool)
    selected = np.zeros(len(df), dtype=np.int8)

    rows = np.flatnonzero(df["Matches"].to_numpy())  # both directions of every UID match
    if len(rows) == 0:
        return selected

    forward = r1[rows] < r2[rows]
    keys = ["sess_lo", "sess_hi", "unit_lo", "unit_hi"]
    candidates = pd.DataFrame({
        "sess_lo": np.where(forward, r1[rows], r2[rows]), "sess_hi": np.where(forward, r2[rows], r1[rows]),
        "unit_lo": np.where(forward, id1[rows], id2[rows]), "unit_hi": np.where(forward, id2[rows], id1[rows]),
        "score": df[rank_col].to_numpy()[rows],
        "kept": forward & resolved[rows],
        "row": rows,
    })
    pairs = candidates.groupby(keys, sort=False).agg(score=("score", "mean"), kept=("kept", "max")).reset_index()
    pairs = pairs[pairs["kept"]].copy()
    pairs["tiebreak"] = rng.random(len(pairs))
    pairs = pairs.sort_values(["score", "tiebreak"], ascending=[False, True])

    rank = pairs.groupby(["sess_lo", "sess_hi"], sort=False).cumcount().to_numpy()
    n_allowed = np.array([n_for_dataset.get((int(a), int(b)), 0) for a, b in zip(pairs["sess_lo"], pairs["sess_hi"])])
    chosen = pairs.loc[rank < n_allowed, keys]

    selected[candidates.merge(chosen, on=keys)["row"].to_numpy()] = 1
    return selected


def _retrack_to_n(dataset, n_for_dataset, rng, base_output=BASE_OUTPUT):
    """
    Lowers RETRACKED_MODEL's tracking threshold per session pair until it has
    n_for_dataset[(r1, r2)] matches. Returns (rows of its MatchTable.csv that
    are matches, in both directions; one info dict per session pair).
    """
    source_model, _ = dvd._uid_variant_source_model(RETRACKED_MODEL)
    table = _read_match_table(dataset, RETRACKED_MODEL, base_output)
    if table is None:
        return np.array([], dtype=np.int64), []

    # A MatchTable is an (n_units, n_units) matrix with ID1 cycling fastest
    # (save_utils.make_match_table), so a row number gives both units.
    n_units = int(np.sqrt(len(table)))
    probs = table[dvd.MATCH_CONFLICT_RANK_COL].to_numpy().reshape(n_units, n_units)
    clus_info = {
        "original_ids": table["ID1"].to_numpy()[:n_units].reshape(-1, 1),
        "session_id": table["RecSes 1"].to_numpy()[:n_units] - 1,
    }
    with open(os.path.join(base_output, *dataset.split("/"), source_model, "UMparam.pickle"), "rb") as f:
        param = pickle.load(f)

    df = _across_sessions(table, dataset, base_output)
    rows = df.index.to_numpy()
    unit_1, unit_2 = rows % n_units, rows // n_units
    r1, r2 = df["RecSes 1"].to_numpy(), df["RecSes 2"].to_numpy()
    n_sessions = int(max(r1.max(), r2.max())) + 1
    pair_code = np.minimum(r1, r2).astype(np.int64) * n_sessions + np.maximum(r1, r2)

    thresholds = [param["match_threshold"]] + [t for t in RETRACK_THRESHOLDS if t < param["match_threshold"]]
    remaining = dict(n_for_dataset)
    selected_rows, info = [], []
    for threshold in thresholds:
        if not remaining:
            break
        with contextlib.redirect_stdout(io.StringIO()):
            uids = aid.assign_unique_id(probs, {**param, "match_threshold": threshold, "use_data_driven_prob_thrs": False}, clus_info)[1]
        matches = (uids[unit_1] == uids[unit_2]).astype(np.int8)
        chosen = _best_n_uid_matches(df.assign(Matches=matches), remaining, rng).astype(bool)
        forward = chosen & (r1 < r2)
        found = dict(zip(*np.unique(pair_code[forward], return_counts=True)))

        done = []
        for (s1, s2), n in remaining.items():
            n_found = int(found.get(s1 * n_sessions + s2, 0))
            if n_found >= n or threshold == thresholds[-1]:
                done.append((s1, s2))
                info.append({"dataset": dataset, "RecSes 1": s1, "RecSes 2": s2, "N": n,
                             "n_matches": n_found, "threshold": threshold})
        codes = [s1 * n_sessions + s2 for s1, s2 in done]
        selected_rows.append(rows[chosen & np.isin(pair_code, codes)])
        for pair in done:
            del remaining[pair]

    return np.concatenate(selected_rows) if selected_rows else np.array([], dtype=np.int64), info


def _get_retracked_matches(n_lookup, rng, base_output=BASE_OUTPUT):
    """
    ({dataset: matched rows of RETRACKED_MODEL}, one row per session pair with
    N > 0 saying which threshold it needed), cached in RETRACK_CACHE_PATH.
    """
    settings = {"thresholds": RETRACK_THRESHOLDS, "seed": RANDOM_SEED}
    if os.path.isfile(RETRACK_CACHE_PATH):
        with open(RETRACK_CACHE_PATH, "rb") as f:
            cache = pickle.load(f)
        if cache["settings"] == settings:
            print(f"Loading cached retracking from {RETRACK_CACHE_PATH}")
            return cache["rows"], cache["info"]

    rows, info = {}, []
    for i, (dataset, n_for_dataset) in enumerate(n_lookup.items()):
        if n_for_dataset:
            rows[dataset], dataset_info = _retrack_to_n(dataset, n_for_dataset, rng, base_output)
            info += dataset_info
        print(f"  retracked {i + 1}/{len(n_lookup)}: {dataset}", flush=True)

    info = pd.DataFrame(info, columns=["dataset", "RecSes 1", "RecSes 2", "N", "n_matches", "threshold"])
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    with open(RETRACK_CACHE_PATH, "wb") as f:
        pickle.dump({"settings": settings, "rows": rows, "info": info}, f)
    return rows, info


def main():
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    session_date_lookup = dvd.get_session_date_lookup()
    n_datasets_with_dates = sum(1 for m in session_date_lookup.values() if m)
    print(f"Have session dates for {n_datasets_with_dates}/{len(session_date_lookup)} dataset(s).\n")

    n_lookup = _build_n_lookup(session_date_lookup)
    rng = np.random.default_rng(RANDOM_SEED)
    retracked_rows, retrack_info = _get_retracked_matches(n_lookup, rng)

    retrack_csv = os.path.join(OUTPUT_DIR, "retracking_thresholds_fixed_n.csv")
    retrack_info.to_csv(retrack_csv, index=False)
    print(f"Wrote {retrack_csv}")

    def matches_transform(dataset, df, model):
        if model == N_MODEL:
            return _best_n_uid_matches(df, n_lookup.get(dataset, {}), rng)
        return np.isin(df.index.to_numpy(), retracked_rows.get(dataset, [])).astype(np.int8)

    n_short = int((retrack_info["n_matches"] < retrack_info["N"]).sum())
    with open(os.path.join(OUTPUT_DIR, "fixed_n_run_info.txt"), "w") as f:
        f.write(
            f"N per session pair = number of UID matches of {N_MODEL}\n"
            f"{RETRACKED_MODEL} retracked at thresholds {list(RETRACK_THRESHOLDS)}\n"
            f"session pairs with N > 0 = {len(retrack_info)}\n"
            f"of which {RETRACKED_MODEL} did not reach N at the lowest threshold = {n_short}\n"
            f"random_seed = {RANDOM_SEED}\n"
        )
    print(f"{n_short}/{len(retrack_info)} session pair(s) did not reach N at the lowest threshold.")

    (
        score_bins, rate_bins, count_bins, dataset_level_rows, dataset_level_auc_rows,
        auc_long_df, dataset_bin_auc_rows, dataset_bin_rate_rows,
    ) = dvd.collect_binned_pairs(
        session_date_lookup,
        models_to_include=COMPARED_MODELS,
        matches_transform=matches_transform,
        extra_usecols={"ID2": "int32", "UM Probabilities": "float32"},
        # _best_n_uid_matches() already resolves one-to-many conflicts itself,
        # before picking the best N -- a second pass here could shrink a
        # model's selected set below N.
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

    # Compact key-scores paper figure (DIFF_SUMMARY_MODEL_A vs
    # DIFF_SUMMARY_MODEL_B, dvd.KEY_SCORES_FOR_DIFF_SUMMARY only) -- the
    # fixed-N counterpart to plot_auc_vs_delta_days.py main()'s own version
    # of this same figure (dvd.plot_key_scores_diff_figure()).
    if stats_available:
        # Mouse-level and dataset-level versions.
        dvd.make_key_scores_figures(
            DIFF_SUMMARY_MODEL_A, DIFF_SUMMARY_MODEL_B, dataset_bin_auc_rows, dataset_bin_rate_rows, colour_for,
            OUTPUT_DIR, file_tag="_fixed_n",
        )
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
        title="Quality vs quantity, fixed N (pooled across datasets/ΔDay, N = smaller UID match count)",
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
