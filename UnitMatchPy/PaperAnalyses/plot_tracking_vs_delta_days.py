# Paper figures: P(track) and the AUC of the key functional scores vs signed
# ΔDay, DeepUnitMatch vs UMPy tracking (UID) matches, at free and fixed N.
#
# Free N:  each model's own UID matches, one-to-many conflicts resolved.
# Fixed N: per session pair, N = number of DeepUnitMatch UID matches; UMPy's
#          tracking threshold is lowered until it has N (see
#          plot_auc_vs_delta_days_fixed_n.py, whose selection is reused here).
#
# Per session pair: P(track) (units of RecSes 1 tracked into RecSes 2), number
# of matched pairs, and one AUC per score. P(track) uses every session pair;
# AUCs only session pairs where DeepUnitMatch or UMPy has >= 20 matches. Mouse
# value = median over its session pairs per ΔDay bin; curves are mean +/- SEM
# over mice. Statistics: see plot_auc_vs_delta_days.test_model_by_days()
# (titles) and test_per_bin_model_effect() (stars).
#
# Short entry point for the two figures. Reads every MatchTable.csv once for
# both, and reuses plot_auc_vs_delta_days.py for the helpers, statistics and
# plotting, so the figures are computed exactly as there.
#
# Output: BASE_OUTPUT/tracking_vs_delta_days_report.

import os
import sys

import numpy as np
import pandas as pd

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _HERE)

import plot_auc_vs_delta_days as dvd
import plot_auc_vs_delta_days_fixed_n as fixed_n
import plot_auc_summary as auc_summary_mod

MODEL_A = fixed_n.N_MODEL  # DeepUnitMatch_AssignUniqueID, sets N
MODEL_B = fixed_n.RETRACKED_MODEL  # UMPy_AssignUniqueID, retracked to N
SCORES = dvd.KEY_SCORES_FOR_DIFF_SUMMARY
OUTPUT_DIR = os.path.join(dvd.BASE_OUTPUT, "tracking_vs_delta_days_report")


def read_match_table(dataset, model, dates, flagged_units):
    """
    Across-session rows of the model's source MatchTable.csv with "Matches" its
    UID matches and a "delta_days" column, without pairs touching a
    within-session duplicate or a session without a date. Keeps the original
    row numbers as index (the fixed-N selection refers to them). None if absent.
    """
    source_model, (uid_a, uid_b) = dvd._uid_variant_source_model(model)
    dataset_dir = os.path.join(dvd.BASE_OUTPUT, *dataset.split("/"))
    path = os.path.join(dataset_dir, source_model, "MatchTable.csv")
    if not os.path.isfile(path) or not os.path.isfile(os.path.join(dataset_dir, model, "AUC_summary.json")):
        return None

    header = pd.read_csv(path, nrows=0).columns
    scores = [s for s in SCORES + [auc_summary_mod.FR_DIFF_SCORE] if s in header]
    dtypes = {"ID1": "int32", "ID2": "int32", "RecSes 1": "int32", "RecSes 2": "int32",
              uid_a: "int32", uid_b: "int32", dvd.MATCH_CONFLICT_RANK_COL: "float32"}
    dtypes.update({s: "float32" for s in scores})
    df = pd.read_csv(path, usecols=list(dtypes), dtype=dtypes)
    if auc_summary_mod.FR_DIFF_SCORE in df.columns:
        auc_summary_mod.add_fr_diff_norm_column(df, os.path.join(dataset_dir, source_model))

    df["Matches"] = (df[uid_a] == df[uid_b]).astype(np.int8)
    df = df[df["RecSes 1"] != df["RecSes 2"]]
    df, _ = dvd._drop_within_session_duplicate_pairs(df, flagged_units)
    date1, date2 = df["RecSes 1"].map(dates), df["RecSes 2"].map(dates)
    have_dates = date1.notna() & date2.notna()
    df = df[have_dates].copy()
    df["delta_days"] = (pd.to_datetime(date2[have_dates]) - pd.to_datetime(date1[have_dates])).dt.days.astype(int)
    return df


def session_pair_rows(df, matches, mouse, dataset, model, rng):
    """(P(track) rows, AUC rows), one per ordered session pair (and score), in the format of make_key_scores_figures()."""
    matches = np.asarray(matches).astype(bool)
    id1 = df["ID1"].to_numpy()
    days = df["delta_days"].to_numpy()
    scores = [s for s in SCORES if s in df.columns and s not in dvd.SCORES_EXCLUDED_PER_MOUSE.get(mouse, ())]
    metric = {s: (-1 if s in dvd.NEGATE_FOR_AUC else 1) * df[s].to_numpy() for s in scores}

    rate_rows, auc_rows = [], []
    for (r1, r2), idx in df.groupby(["RecSes 1", "RecSes 2"], sort=False).indices.items():
        delta = int(days[idx[0]])
        base = {"mouse": mouse, "dataset": dataset, "model": model, "RecSes 1": int(r1), "RecSes 2": int(r2),
                "delta_days": delta, "bin": dvd._signed_bin_label(dvd._signed_bin_ids(np.array([delta]))[0])}
        units = id1[idx]
        rate = len(np.unique(units[matches[idx]])) / len(np.unique(units))
        rate_rows.append({**base, "rate": rate, dvd.N_MATCHES_COL: float(matches[idx].sum())})
        for s in scores:
            valid = idx[np.isfinite(metric[s][idx])]
            if len(valid) == 0:
                continue
            valid = dvd._cap_negatives(valid, matches, dvd.BOOTSTRAP_MAX_N, rng)
            auc = dvd._auc_from_flat(matches[valid], metric[s][valid], punish_missing=dvd.PUNISH_MISSING_MATCHES)
            if np.isfinite(auc):
                auc_rows.append({**base, "score": s, "auc": auc})
    return rate_rows, auc_rows


def keep_auc_rows(rate_rows, auc_rows, min_matches=dvd.MIN_MATCHES_PER_SESSION_PAIR):
    """AUC rows of session pairs where at least one model has >= min_matches matches."""
    rates = pd.DataFrame(rate_rows)
    best = rates.groupby(["dataset", "RecSes 1", "RecSes 2"])[dvd.N_MATCHES_COL].max()
    keep = set(best.index[best >= min_matches])
    kept = [r for r in auc_rows if (r["dataset"], r["RecSes 1"], r["RecSes 2"]) in keep]
    print(f"AUCs: kept {len(keep)}/{len(best)} session pairs with >= {min_matches} matches in a model.")
    return kept


def main():
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    date_lookup = dvd.get_session_date_lookup()
    datasets = sorted(d for d, dates in date_lookup.items() if dates)
    flagged = dvd._build_within_session_duplicate_lookup(dvd.BASE_OUTPUT, datasets)

    # Fixed N: N per session pair from DeepUnitMatch, UMPy's retracked matches (cached).
    n_lookup = fixed_n._build_n_lookup(date_lookup)
    selection_rng = np.random.default_rng(fixed_n.RANDOM_SEED)
    retracked_rows, retrack_info = fixed_n._get_retracked_matches(n_lookup, selection_rng)
    retrack_info.to_csv(os.path.join(OUTPUT_DIR, "retracking_thresholds_fixed_n.csv"), index=False)
    n_short = int((retrack_info["n_matches"] < retrack_info["N"]).sum())
    print(f"Fixed N: {MODEL_B} short of N in {n_short}/{len(retrack_info)} session pairs.")

    auc_rng = np.random.default_rng(0)  # subsampling non-matches of very large session pairs
    rows = {"free_n": ([], []), "fixed_n": ([], [])}
    for i, dataset in enumerate(datasets):
        mouse = auc_summary_mod.mouse_from_dataset(dataset)
        for model in (MODEL_A, MODEL_B):
            df = read_match_table(dataset, model, date_lookup[dataset], flagged.get(dataset, set()))
            if df is None or df.empty:
                continue
            if model == MODEL_A:
                fixed = fixed_n._best_n_uid_matches(df, n_lookup.get(dataset, {}), selection_rng)
            else:
                fixed = np.isin(df.index.to_numpy(), retracked_rows.get(dataset, []))
            for name, matches in (("free_n", dvd._resolve_match_conflicts(df)), ("fixed_n", fixed)):
                rate_rows, auc_rows = session_pair_rows(df, matches, mouse, dataset, model, auc_rng)
                rows[name][0].extend(rate_rows)
                rows[name][1].extend(auc_rows)
        print(f"  {i + 1}/{len(datasets)} {dataset}", flush=True)

    colours = dvd.build_family_colours([MODEL_A, MODEL_B])
    for name, (rate_rows, auc_rows) in rows.items():
        print(f"\n{name}:")
        auc_rows = keep_auc_rows(rate_rows, auc_rows)
        pd.DataFrame(rate_rows).to_csv(os.path.join(OUTPUT_DIR, f"session_pairs_ptrack_{name}.csv"), index=False)
        pd.DataFrame(auc_rows).to_csv(os.path.join(OUTPUT_DIR, f"session_pairs_auc_{name}.csv"), index=False)
        dvd.make_key_scores_figures(MODEL_A, MODEL_B, auc_rows, rate_rows, colours, OUTPUT_DIR,
                                    file_tag=f"_{name}", levels=("mouse",))


if __name__ == "__main__":
    main()
