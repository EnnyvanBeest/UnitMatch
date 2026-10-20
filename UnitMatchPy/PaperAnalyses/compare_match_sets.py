# Compare the match sets of two matching runs (default: DeepUnitMatch vs
# DUM_legacy) on every location where both exist, and test whether the
# matches only one of them finds are real.
#
# Per location (across-session unit pairs only; a pair is a match of a run if
# that run called it in either direction):
#   - overlap: pairs both runs match, only run A, only run B; fraction of A's
#     matches that B also has (and vice versa).
#   - AUC of each group (both / A only / B only) against the pairs neither run
#     matches, per functional score (direction-averaged). AUC > 0.5 = the
#     group's pairs look more like the same neuron than non-matches do, i.e.
#     the extra matches carry real same-neuron evidence.
# Across mice: per mouse the median over its locations, then a Wilcoxon
# signed-rank test of AUC - 0.5 (is the group better than non-matches?).
#
# refpop_correlations is left out: each run builds it from its own matches,
# so it is not comparable between the two runs' groups.
#
# Read-only on the matching outputs; writes a per-location CSV and a text
# summary to --out-dir (default cfg.REPORTS_DIR).
#
#   python compare_match_sets.py                       # DeepUnitMatch vs DUM_legacy
#   python compare_match_sets.py --a DeepUnitMatch --b UMPy

import argparse
import datetime
import glob
import os

import numpy as np
import pandas as pd
from scipy.stats import wilcoxon
from sklearn.metrics import roc_auc_score

import pipeline_config as cfg

# functional score -> +1 if higher means more similar, -1 if lower does
SCORES = {
    "ISI_correlations": 1, "ISI_KL_divergence": -1, "ISI_wasserstein_distance": -1,
    "FR_diff": -1, "ISI_CV_diff": -1, "natim_correlations": 1,
}
MIN_PAIRS = 5  # a group's AUC at a location needs at least this many pairs


def unordered_pairs(path, columns):
    """Across-session rows of a MatchTable as one row per unordered unit pair
    (Matches: either direction; scores: mean of both directions)."""
    header = pd.read_csv(path, nrows=0).columns
    use = ["RecSes 1", "RecSes 2", "ID1", "ID2", "Matches"] + [c for c in columns if c in header]
    df = pd.read_csv(path, usecols=use)
    df = df[df["RecSes 1"] != df["RecSes 2"]].reset_index(drop=True)
    lo = (df["RecSes 1"] < df["RecSes 2"]).to_numpy()
    pairs = pd.DataFrame({
        "s1": np.where(lo, df["RecSes 1"], df["RecSes 2"]), "u1": np.where(lo, df["ID1"], df["ID2"]),
        "s2": np.where(lo, df["RecSes 2"], df["RecSes 1"]), "u2": np.where(lo, df["ID2"], df["ID1"]),
    })
    values = ["Matches"] + use[5:]
    for c in values:
        pairs[c] = df[c].to_numpy()
    agg = {"Matches": "max", **{c: "mean" for c in values[1:]}}
    return pairs.groupby(["s1", "u1", "s2", "u2"]).agg(agg)


def compare_location(path_a, path_b):
    a = unordered_pairs(path_a, list(SCORES))
    b = unordered_pairs(path_b, [])
    pairs = a.join(b[["Matches"]].rename(columns={"Matches": "Matches_b"}), how="inner")
    in_a, in_b = pairs["Matches"] > 0, pairs["Matches_b"] > 0
    groups = {"both": in_a & in_b, "a_only": in_a & ~in_b, "b_only": in_b & ~in_a}
    neither = ~in_a & ~in_b
    row = {"n_pairs": len(pairs), "n_a": int(in_a.sum()), "n_b": int(in_b.sum()),
           **{f"n_{g}": int(m.sum()) for g, m in groups.items()}}
    row["frac_a_in_b"] = row["n_both"] / row["n_a"] if row["n_a"] else np.nan
    row["frac_b_in_a"] = row["n_both"] / row["n_b"] if row["n_b"] else np.nan
    for score, sign in SCORES.items():
        if score not in pairs:
            continue
        for g, m in groups.items():
            pos = pairs.loc[m, score].dropna()
            neg = pairs.loc[neither, score].dropna()
            if len(pos) >= MIN_PAIRS and len(neg) >= MIN_PAIRS:
                y = np.r_[np.ones(len(pos)), np.zeros(len(neg))]
                row[f"AUC_{score}_{g}"] = roc_auc_score(y, sign * np.r_[pos.to_numpy(), neg.to_numpy()])
    return row


def main():
    parser = argparse.ArgumentParser(description="Overlap of two runs' match sets and AUCs of the shared/unique matches.")
    parser.add_argument("--a", default="DeepUnitMatch")
    parser.add_argument("--b", default="DUM_legacy")
    parser.add_argument("--out-dir", default=cfg.REPORTS_DIR)
    args = parser.parse_args()

    rows = []
    locations = sorted(glob.glob(os.path.join(cfg.ANALYSIS_OUTPUT, "*", "*", "*", args.a, "MatchTable.csv")))
    for path_a in locations:
        loc_dir = os.path.dirname(os.path.dirname(path_a))
        path_b = os.path.join(loc_dir, args.b, "MatchTable.csv")
        if not os.path.isfile(path_b):
            continue
        loc = os.path.relpath(loc_dir, cfg.ANALYSIS_OUTPUT).replace("\\", "/")
        row = compare_location(path_a, path_b)
        if row["n_pairs"] == 0:
            continue  # single session: no across-session pairs
        rows.append({"location": loc, "mouse": loc.split("/")[0], **row})
        print(f"{loc:32s} {args.a} {row['n_a']:5d}  {args.b} {row['n_b']:5d}  both {row['n_both']:5d}")
    df = pd.DataFrame(rows)

    lines = [f"Match sets: A = {args.a}, B = {args.b} ({len(df)} locations, {df['mouse'].nunique()} mice)", ""]
    tot = df[["n_a", "n_b", "n_both", "n_a_only", "n_b_only"]].sum()
    lines += [f"Total across-session matches: A {tot.n_a}, B {tot.n_b}; both {tot.n_both}, "
              f"A only {tot.n_a_only}, B only {tot.n_b_only}",
              f"Fraction of A's matches also in B: {tot.n_both / tot.n_a:.3f}; of B's also in A: {tot.n_both / tot.n_b:.3f}",
              f"Per location (median): A in B {df['frac_a_in_b'].median():.3f}, B in A {df['frac_b_in_a'].median():.3f}", "",
              "AUC of each group vs pairs neither run matches (per mouse: median over locations;",
              "Wilcoxon signed-rank of AUC - 0.5 across mice; 'A only' = matches only A finds):",
              f"{'score':26s} {'group':8s} {'n_mice':>6s} {'median AUC':>10s} {'mice > 0.5':>10s} {'p':>9s}"]
    labels = {"both": "both", "a_only": "A only", "b_only": "B only"}
    for score in SCORES:
        for g, label in labels.items():
            col = f"AUC_{score}_{g}"
            if col not in df:
                continue
            per_mouse = df.groupby("mouse")[col].median().dropna()
            if len(per_mouse) < 2:
                continue
            diff = per_mouse - 0.5
            p = wilcoxon(diff).pvalue if np.any(diff != 0) else np.nan
            lines.append(f"{score:26s} {label:8s} {len(per_mouse):6d} {per_mouse.median():10.3f} "
                         f"{int((diff > 0).sum()):>5d}/{len(per_mouse):<4d} {p:9.3g}")
    text = "\n".join(lines)
    print("\n" + text)

    os.makedirs(args.out_dir, exist_ok=True)
    stamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
    base = os.path.join(args.out_dir, f"match_sets_{args.a}_vs_{args.b}_{stamp}")
    df.to_csv(base + ".csv", index=False)
    with open(base + ".txt", "w", encoding="utf-8") as f:
        f.write(text + "\n")
    print(f"\nWritten to {base}.csv / .txt")


if __name__ == "__main__":
    main()
