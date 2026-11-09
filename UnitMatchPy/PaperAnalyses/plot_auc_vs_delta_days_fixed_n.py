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
import UnitMatchPy.assign_unique_id as aid

BASE_OUTPUT = dvd.BASE_OUTPUT
OUTPUT_DIR = os.path.join(BASE_OUTPUT, "auc_vs_delta_days_fixed_n_report")

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

_flagged_units = {}  # dataset -> within-session duplicates, the same for every model


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
    """
    The fixed-N paper figures only (dvd.make_key_scores_figures()), see the
    top of this file. Shares everything else with plot_auc_vs_delta_days.main().
    """
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
        _, _, _, dataset_level_rows, dataset_level_auc_rows, _, dataset_bin_auc_rows, dataset_bin_rate_rows,
    ) = dvd.collect_binned_pairs(
        session_date_lookup,
        models_to_include=COMPARED_MODELS,
        matches_transform=matches_transform,
        extra_usecols={"ID2": "int32", "UM Probabilities": "float32"},
        # _best_n_uid_matches() already resolves one-to-many conflicts itself,
        # before picking the best N -- a second pass here could shrink a
        # model's selected set below N.
        resolve_match_conflicts=False,
        scores_to_use=dvd.KEY_SCORES_FOR_DIFF_SUMMARY,
        pooled=False,
    )
    dvd.write_long_by_dataset(
        dataset_level_rows, dataset_level_auc_rows, os.path.join(OUTPUT_DIR, "quality_quantity_fixed_n_long_by_dataset.csv"),
    )
    dvd.make_key_scores_figures(
        DIFF_SUMMARY_MODEL_A, DIFF_SUMMARY_MODEL_B, dataset_bin_auc_rows, dataset_bin_rate_rows,
        dvd.build_family_colours([DIFF_SUMMARY_MODEL_A, DIFF_SUMMARY_MODEL_B]), OUTPUT_DIR, file_tag="_fixed_n",
    )


if __name__ == "__main__":
    main()
