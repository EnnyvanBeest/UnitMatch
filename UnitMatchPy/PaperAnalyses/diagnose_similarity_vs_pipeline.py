# Diagnostic: is a model's poor result in the fixed-N comparison due to the
# model (its DNN similarity) or to the matching pipeline on top of it?
#
# For every model and location in matchtables_<group>.db it computes
#   1. within-day separation (model only, no pipeline): ROC AUC of the raw DNN
#      similarity for the same unit (ID1 == ID2, the two halves of the
#      recording) vs. different units, within each session;
#   2. the spread of the similarity (mean/sd of across-session pairs, mean of
#      the same-unit pairs), to see if a model's similarities are compressed;
#   3. the fixed-N functional AUCs (as fast_testing.py, N per session pair set
#      by DeepUnitMatch) twice: ranking pairs by the pipeline's match
#      probability ("UM Probabilities", as in the figures) and by the raw DNN
#      similarity ("similarity", no pipeline);
#   4. how much of the probability ranking is ties: per session pair, the
#      number of across-session pairs with probability >= 0.999 relative to N
#      (> 1 means the fixed-N selection is partly random tie-breaking).
#
# Reading: if a model separates units well within a day and its
# similarity-ranked AUCs are fine but its probability-ranked AUCs drop, the
# drop comes from the pipeline. If the similarity-ranked AUCs drop as well,
# it is the model.
#
# Read-only on the database; writes CSVs and a figure to --out.
#
#   python diagnose_similarity_vs_pipeline.py
#   python diagnose_similarity_vs_pipeline.py --mice AL031 AV008 --workers 2
#   python diagnose_similarity_vs_pipeline.py --models DeepUnitMatch DUM_AE_only

import os
import sys
import argparse
import sqlite3
import traceback
from concurrent.futures import ProcessPoolExecutor, as_completed

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _HERE)

import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from sklearn.metrics import roc_auc_score

import pipeline_config as cfg
import fast_testing as ft
from utils.helpers import create_um_lookup, get_locations_from_sqlite

DEFAULT_MODELS = ["DeepUnitMatch", "DUM_untrained", "DUM_AE_only", "DUM_finetuned_only"]
SATURATED = 0.999
RANKINGS = {"prob": "UM Probabilities", "sim": "similarity"}
AUC_COLS = ["AUC_isi", "AUC_fr_norm", "AUC_refpop_model", "AUC_natim"]


def read_location(db_path, mouse, probe, loc, models):
    con = sqlite3.connect(f"file:{db_path}?mode=ro", uri=True)
    try:
        table = f"{mouse}_{probe}_{loc}"
        present = {c[1] for c in con.execute(f"PRAGMA table_info({ft.quote_ident(table)})")}
        base = ["ID1", "ID2", "RecSes1", "RecSes2", "ISI_correlations", "ISI_KL_divergence",
                "ISI_wasserstein_distance", "FR_diff", "ISI_CV_diff",
                "natim_correlations_DeepUnitMatch", "refpop_correlations_DeepUnitMatch",
                "refpop_correlations_UMPy"]
        per_model = [f"{p}_{m}" for m in models
                     for p in ("UM Probabilities", "similarity", "refpop_correlations")]
        wanted = list(dict.fromkeys(base + per_model))
        missing = [c for c in per_model if c not in present and not c.startswith("refpop")]
        cols = [c for c in wanted if c in present]
        mt = pd.read_sql_query(
            f"SELECT {', '.join(ft.quote_ident(c) for c in cols)} FROM {ft.quote_ident(table)};", con)
    finally:
        con.close()
    for c in mt.columns:
        mt[c] = pd.to_numeric(mt[c], errors="coerce")
    for c in wanted:  # e.g. no natim/refpop scores at this location: NaN, as in fast_testing
        if c not in mt:
            mt[c] = np.nan
    return mt, missing


def within_day_stats(mt, model):
    """Same unit vs. different units within each session, on the raw similarity."""
    sim = f"similarity_{model}"
    within = mt.loc[mt["RecSes1"] == mt["RecSes2"], ["RecSes1", "ID1", "ID2", sim]].dropna()
    across = mt.loc[mt["RecSes1"] != mt["RecSes2"], sim].dropna()
    aucs = []
    for _, s in within.groupby("RecSes1"):
        same = (s["ID1"] == s["ID2"]).to_numpy()
        if 0 < same.sum() < len(same):
            aucs.append(roc_auc_score(same, s[sim].to_numpy()))
    same_vals = within.loc[within["ID1"] == within["ID2"], sim]
    return {
        "withinday_AUC": np.mean(aucs) if aucs else np.nan,
        "sim_same_mean": same_vals.mean(),
        "sim_across_mean": across.mean(),
        "sim_across_sd": across.std(),
        # distance of the same-unit pairs from the bulk, in units of its spread
        "sim_same_z": (same_vals.mean() - across.mean()) / across.std() if across.std() > 0 else np.nan,
    }


def saturation(mt, model, n_lookup):
    """Per session pair: pairs with probability >= SATURATED, relative to the fixed N."""
    prob = f"UM Probabilities_{model}"
    rows = []
    for _, r in n_lookup.iterrows():
        d = mt.loc[(mt["RecSes1"] == r["r1"]) & (mt["RecSes2"] == r["r2"]), prob]
        n_sat = int((d >= SATURATED).sum())
        rows.append({"r1": r["r1"], "r2": r["r2"], "N": r["N"], "n_saturated": n_sat,
                     "saturated_over_N": n_sat / r["N"] if r["N"] > 0 else np.nan})
    return pd.DataFrame(rows)


def process_location(db_path, location, models, n_lookup):
    mouse, probe, loc = location
    try:
        mt, missing = read_location(db_path, mouse, probe, loc, models)
        if missing:
            print(f"{mouse}/{probe}/{loc}: missing columns {missing}")
        out_loc, out_pairs = [], []
        for model in models:
            if f"similarity_{model}" not in mt or f"UM Probabilities_{model}" not in mt:
                continue
            row = {"mouse": mouse, "probe": probe, "loc": loc, "model": model}
            row.update(within_day_stats(mt, model))
            if n_lookup is not None:
                sat = saturation(mt, model, n_lookup)
                row["saturated_over_N"] = sat["saturated_over_N"].mean()
            out_loc.append(row)

            if n_lookup is None:
                continue
            for rank_name, prefix in RANKINGS.items():
                # alias the ranking column (and its refpop column) under a
                # "UM Probabilities_*" name, which all_results_1model expects
                alias = f"{model}__{rank_name}rank"
                d = mt.copy()
                d[f"UM Probabilities_{alias}"] = d[f"{prefix}_{model}"]
                d[f"refpop_correlations_{alias}"] = d[f"refpop_correlations_{model}"]
                res = ft.all_results_1model(d, mouse, probe, loc,
                                            model_name=f"UM Probabilities_{alias}", fixed_n=n_lookup)
                if res is None:
                    continue
                res = res[["mouse", "probe", "loc", "r1", "r2", "N"] + AUC_COLS].copy()
                res["model"], res["ranking"] = model, rank_name
                out_pairs.append(res)
        print(f"done {mouse}/{probe}/{loc}")
        return out_loc, out_pairs
    except Exception:
        print(f"FAILED {mouse}/{probe}/{loc}\n{traceback.format_exc()}")
        return [], []


def plot(per_mouse, models, out_dir):
    panels = [("withinday_AUC", "within-day AUC\n(raw similarity, no pipeline)"),
              ("sim_same_z", "same-unit similarity\n(z vs across-session pairs)"),
              ("saturated_over_N", "pairs with P >= 0.999 / N")]
    panels += [(f"{c}|{r}", f"{c}\nranked by {'probability' if r == 'prob' else 'raw similarity'}")
               for c in ("AUC_isi", "AUC_refpop_model", "AUC_natim") for r in ("prob", "sim")]
    fig, axes = plt.subplots(1, len(panels), figsize=(3.2 * len(panels), 4), squeeze=False)
    x = np.arange(len(models))
    for ax, (col, title) in zip(axes[0], panels):
        if col not in per_mouse:
            ax.set_visible(False)
            continue
        wide = per_mouse.pivot(index="mouse", columns="model", values=col).reindex(columns=models)
        for _, r in wide.iterrows():
            ax.plot(x, r.to_numpy(), color="0.8", lw=0.8, zorder=1)
        ax.scatter(np.repeat(x[None], len(wide), 0).ravel(), wide.to_numpy().ravel(), s=12,
                   color="tab:red", alpha=0.6, zorder=2)
        ax.errorbar(x, wide.mean(), wide.sem(), fmt="o", color="k", mfc="w", zorder=3)
        ax.set_xticks(x, models, rotation=45, ha="right")
        ax.set_title(title, fontsize=9)
        ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    path = os.path.join(out_dir, "similarity_vs_pipeline.png")
    fig.savefig(path, dpi=150)
    plt.close(fig)
    return path


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--group", default="models", choices=sorted(cfg.MODEL_GROUPS))
    p.add_argument("--models", nargs="+", default=DEFAULT_MODELS)
    p.add_argument("--mice", nargs="+", default=None, help="restrict to these mice")
    p.add_argument("--ref-model", default="DeepUnitMatch", help="model whose free-N count sets N")
    p.add_argument("--out", default=os.path.join(cfg.REPORTS_DIR, "similarity_vs_pipeline"))
    p.add_argument("--workers", type=int, default=min((os.cpu_count() or 2) - 1, 6))
    args = p.parse_args()

    db_path = cfg.database_path(args.group)
    os.makedirs(args.out, exist_ok=True)
    ref_path = os.path.join(cfg.RESULTS_DIR, args.group, f"UM Probabilities_{args.ref_model}_results.csv")
    um_lookup = create_um_lookup(pd.read_csv(ref_path)) if os.path.isfile(ref_path) else {}
    if not um_lookup:
        print(f"{ref_path} not found: fixed-N AUCs skipped (run fast_testing.py --group {args.group} first)")

    locations = get_locations_from_sqlite(db_path)
    if args.mice:
        locations = [l for l in locations if l[0] in set(args.mice)]
    print(f"{len(locations)} locations, models {args.models}, db {db_path}")

    loc_rows, pair_rows = [], []
    with ProcessPoolExecutor(max_workers=max(args.workers, 1)) as ex:
        futures = [ex.submit(process_location, db_path, l, args.models, um_lookup.get(l))
                   for l in locations]
        for f in as_completed(futures):
            a, b = f.result()
            loc_rows += a
            pair_rows += b

    per_loc = pd.DataFrame(loc_rows)
    per_loc.to_csv(os.path.join(args.out, "per_location.csv"), index=False)
    per_mouse = per_loc.groupby(["mouse", "model"]).mean(numeric_only=True).reset_index()

    if pair_rows:
        pairs = pd.concat(pair_rows, ignore_index=True)
        pairs.to_csv(os.path.join(args.out, "per_session_pair_fixedN.csv"), index=False)
        auc_mouse = pairs.groupby(["mouse", "model", "ranking"])[AUC_COLS].mean().unstack("ranking")
        auc_mouse.columns = [f"{c}|{r}" for c, r in auc_mouse.columns]
        per_mouse = per_mouse.merge(auc_mouse.reset_index(), on=["mouse", "model"], how="left")
    per_mouse.to_csv(os.path.join(args.out, "per_mouse.csv"), index=False)

    summary = per_mouse.groupby("model").mean(numeric_only=True).reindex(args.models).T
    summary.to_csv(os.path.join(args.out, "summary.csv"))
    with pd.option_context("display.width", 200, "display.float_format", "{:.3f}".format):
        print(summary)
    print("figure:", plot(per_mouse, args.models, args.out))


if __name__ == "__main__":
    main()
