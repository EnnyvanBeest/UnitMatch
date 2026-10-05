# Summary of the training comparison (train_paper_models.py --jobs comparison):
# which way of training -- A_original (pre-fix code), B_fixed (fixed v1) or
# C_v2 (fixed v1 + v2 options) -- matches neurons best on held-out mice.
#
# Criterion (fixed before looking at the results):
#   Per held-out mouse, from every evaluated location's AUC_summary.json
#   (whole-location AUC with each model's own matches, P > 0.5):
#     ISI correlations, reference-population correlations, natural-image
#     correlations, firing-rate difference (AUC, higher = better), and the
#     number of across-session matches.
#   Each value: median over the mouse's locations, then mean over the
#   replicates (3-mouse training subsets) in which that mouse was held out.
#   Variants are compared pairwise across mice (paired Wilcoxon signed-rank).
#   A variant is preferred over B_fixed only if it is significantly better on
#   at least one AUC and not significantly worse on any AUC or on the number
#   of matches; otherwise B_fixed (the simplest) is kept.
#
# Read-only on the results; writes a CSV + text report to cfg.REPORTS_DIR.

import datetime
import json
import os
import sys

import numpy as np
import pandas as pd
from scipy.stats import wilcoxon

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _HERE)
import pipeline_config as cfg

VARIANTS = ["A_original", "B_fixed", "C_v2"]
AUC_SCORES = ["ISI_correlations", "refpop_correlations", "natim_correlations", "FR_diff"]
COUNT_SCORE = "n_matches_across_sessions"
REPLICATES = [1, 2, 3]
ALPHA = 0.05


def collect():
    with open(os.path.join(cfg.TRAINING_STATE_ROOT, "manifest.json")) as f:
        manifest = json.load(f)
    rows = []
    for mouse in sorted(os.listdir(cfg.ANALYSIS_OUTPUT)):
        mouse_dir = os.path.join(cfg.ANALYSIS_OUTPUT, mouse)
        if not os.path.isdir(mouse_dir):
            continue
        for root, dirs, files in os.walk(mouse_dir):
            name = os.path.basename(root)
            if not name.startswith("cmp_m3_") or "AUC_summary.json" not in files:
                continue
            dirs[:] = []
            rep, variant = int(name.split("_")[2]), "_".join(name.split("_")[3:])
            if variant not in VARIANTS or mouse in manifest[f"m3_{rep}"]:
                continue  # unknown condition, or a training mouse (should not have been evaluated)
            with open(os.path.join(root, "AUC_summary.json")) as f:
                summary = json.load(f)
            location = os.path.relpath(os.path.dirname(root), cfg.ANALYSIS_OUTPUT).replace(os.sep, "/")
            for score in AUC_SCORES + [COUNT_SCORE]:
                value = summary.get(score)
                if isinstance(value, (int, float)) and np.isfinite(value):
                    rows.append({"mouse": mouse, "location": location, "rep": rep, "variant": variant,
                                 "score": score, "value": float(value)})
    return pd.DataFrame(rows)


def per_mouse(df):
    loc_median = df.groupby(["score", "variant", "rep", "mouse"])["value"].median()
    return loc_median.groupby(["score", "variant", "mouse"]).mean().unstack("variant")


def compare(table, a, b):
    """Rows of (score, median difference a-b, n mice a > b, n, p) across mice with both variants."""
    out = []
    for score in AUC_SCORES + [COUNT_SCORE]:
        if score not in table.index.get_level_values(0):
            continue
        pair = table.loc[score][[a, b]].dropna()
        diff = pair[a] - pair[b]
        p = wilcoxon(diff).pvalue if len(diff) > 1 and np.any(diff != 0) else np.nan
        out.append({"comparison": f"{a} vs {b}", "score": score, "median_diff": float(diff.median()) if len(diff) else np.nan,
                    "n_better": int((diff > 0).sum()), "n_mice": int(len(diff)), "p": p})
    return out


def verdict(results, variant):
    rows = [r for r in results if r["comparison"] == f"{variant} vs B_fixed"]
    sig = lambda r: np.isfinite(r["p"]) and r["p"] < ALPHA
    worse = [r["score"] for r in rows if sig(r) and r["median_diff"] < 0]
    better = [r["score"] for r in rows if sig(r) and r["median_diff"] > 0 and r["score"] in AUC_SCORES]
    if worse:
        return False, f"{variant}: significantly worse than B_fixed on {worse}"
    if better:
        return True, f"{variant}: significantly better than B_fixed on {better}, not worse on anything"
    return False, f"{variant}: no significant improvement over B_fixed"


def main():
    df = collect()
    if df.empty:
        print("No evaluated comparison results found yet.")
        return
    lines = ["Training comparison (held-out mice)", ""]
    coverage = df.groupby(["rep", "variant"])["location"].nunique().unstack("variant")
    lines += ["Evaluated locations per replicate and variant:", coverage.to_string(), ""]

    table = per_mouse(df)
    lines += ["Per-variant medians across mice:",
              table.groupby(level=0).median().round(4).to_string(), ""]

    results = compare(table, "B_fixed", "A_original") + compare(table, "C_v2", "B_fixed") + compare(table, "C_v2", "A_original")
    res = pd.DataFrame(results)
    lines += ["Paired comparisons across mice (Wilcoxon signed-rank):",
              res.assign(median_diff=res["median_diff"].round(4), p=res["p"].map(lambda p: f"{p:.3g}")).to_string(index=False), ""]

    keep_c, why_c = verdict(results, "C_v2")
    lines += ["Decision (criterion fixed in advance):", "  " + why_c,
              f"  -> use {'C_v2' if keep_c else 'B_fixed'} for the paper models", "",
              "B_fixed vs A_original shows what the training fixes alone changed."]
    text = "\n".join(lines)
    print(text)

    os.makedirs(cfg.REPORTS_DIR, exist_ok=True)
    stamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
    df.to_csv(os.path.join(cfg.REPORTS_DIR, f"training_comparison_{stamp}.csv"), index=False)
    with open(os.path.join(cfg.REPORTS_DIR, f"training_comparison_{stamp}.txt"), "w") as f:
        f.write(text + "\n")
    print(f"\nWritten to {cfg.REPORTS_DIR}")


if __name__ == "__main__":
    main()
