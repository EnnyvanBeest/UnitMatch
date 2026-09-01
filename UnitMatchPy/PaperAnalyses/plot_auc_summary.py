# Quick analysis script: collect every AUC_summary.json written by the batch
# scripts (save_utils.save_auc_summary is called once per model/method per
# dataset group by run_deepunitmatch_batch_onMerged.py, run_emd_batch_onMerged.py,
# run_dant_batch_onMerged.py, and the maxdist/extramodels sweep scripts -- see
# their get_*_dir/get_*save_dir helpers), plot each functional score (and the
# match count) as one point per dataset grouped by model and coloured by
# mouse, and fit a linear mixed-effects model per score (fixed effect: model,
# random intercept: mouse) so model differences aren't confounded with
# mouse-to-mouse variability (some mice contribute more datasets than others).
#
# Every AUC_summary.json lives at BASE_OUTPUT/<dataset>/<model>/AUC_summary.json,
# where <dataset> is the group's path relative to BASE_OUTPUT (e.g.
# "AL032/19011111882/2") and <model> is the immediate parent folder name
# ("DeepUnitMatch", "UMPy", "EMD", "DANT", "DANT_no_functional",
# "DUM_maxdist=20", "n_output=64_after_ae_and_finetune", ...). The mouse is
# taken as the first path component of <dataset>.
#
# Two parallel analyses are produced for every score:
#   - dataset-level: every dataset is its own point, mouse-to-mouse variation
#     is handled by a random intercept in a linear mixed-effects model.
#   - mouse-averaged: each mouse's datasets are averaged first into a single
#     point per mouse per model, so the simpler paired-t-test/ANOVA machinery
#     (no random effects needed) can be used instead -- at the cost of
#     discarding within-mouse dataset-to-dataset variability.
# Both get pairwise post-hoc comparisons across all models (Holm-Bonferroni
# corrected across pairs within a score) drawn as significance brackets on
# the plots for pairs whose corrected p < 0.05.
#
# Output (plots, CSVs, stats text files) is written to BASE_OUTPUT/auc_summary_report
# alongside the rest of the pipeline's output, not into this git checkout.

import os
import sys
import json
import pickle
import warnings
from itertools import combinations

import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")  # non-interactive backend
import matplotlib.pyplot as plt
import statsmodels.formula.api as smf
from statsmodels.stats.multitest import multipletests
from scipy.stats import ttest_rel

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _HERE)
sys.path.insert(0, os.path.join(_HERE, "DeepUnitMatch"))

import run_deepunitmatch_batch_onMerged as base_batch
from DeepUnitMatch.testing import test as dumtest

# ── settings ─────────────────────────────────────────────────────────────────
# Point this at whichever BASE_OUTPUT tree you want to summarise (defaults to
# the merged-data comparison tree, where DUM/UMPy/EMD/DANT/sweeps all live).
BASE_OUTPUT = base_batch.BASE_OUTPUT
OUTPUT_DIR = os.path.join(BASE_OUTPUT, "auc_summary_report")

# Restrict to specific model/method subfolder names (e.g.
# {"DeepUnitMatch", "UMPy", "EMD", "DANT", "DANT_no_functional"}), or None to
# include every AUC_summary.json found under BASE_OUTPUT (main models,
# maxdist sweep, extra checkpoints, everything).
# MODELS_TO_INCLUDE = None
#MODELS_TO_INCLUDE = {"DeepUnitMatch_AssignUniqueID_Conservative", "DeepUnitMatch_AssignUniqueID", "UMPy_AssignUniqueID_Conservative", "UMPy_AssignUniqueID", "EMD", "DANT", "DANT_no_functional"}
MODELS_TO_INCLUDE = {"DeepUnitMatch_AssignUniqueID", "UMPy_AssignUniqueID"}

#MODELS_TO_INCLUDE = {"DUM_NewModelAug2026","DeepUnitMatch","DeepUnitMatch_AssignUniqueID_Conservative", "DeepUnitMatch_AssignUniqueID", "UMPy","UMPy_AssignUniqueID_Conservative", "UMPy_AssignUniqueID", "EMD", "DANT", "DANT_no_functional"}
#MODELS_TO_INCLUDE = {"DeepUnitMatch","n_output=256_after_ae_and_finetune","n_output=128_after_ae_and_finetune","n_output=32_after_ae_and_finetune","n_output=8_after_ae_and_finetune"}
# Reference level for the "vs reference" summaries (every other model's
# coefficient/mean-difference is "difference from this model"). Falls back to
# whichever model sorts first for a given score if this one has no data there.
REFERENCE_MODEL = "UMPy_AssignUniqueID"

# Counts are heavily right-skewed; fit stats on log1p(count) instead of the
# raw value (plots still show raw counts).
LOG_TRANSFORM_FOR_STATS = {"n_matches_across_sessions"}

ALPHA = 0.05

# A dataset is kept only if at least one included model found this many
# across-session matches (n_matches_across_sessions) -- datasets where every
# model stayed below this never demonstrated enough tracked matches to trust
# a score/AUC comparison on, regardless of which model is being looked at.
# Also used by plot_auc_vs_delta_days.py (via get_passing_datasets()) to
# apply the same dataset-level cutoff there.
MIN_MATCHES_TO_INCLUDE = 20

# ── FR_diff_norm: on-the-fly firing-rate normalization ──────────────────────
#
# fast_testing.py (the non-merged/SQL pipeline) now normalizes FR_diff by
# each source unit's own mean firing rate (fast_testing.normalize_FR_diff)
# instead of using the raw |firing-rate difference|. This section reproduces
# that same normalization (divide FR_diff by FR, from
# DeepUnitMatch.testing.test.get_FR) for the merged-tree/BASE_OUTPUT pipeline
# these plotting scripts read from -- computed on the fly, from each model's
# own UMparam.pickle, every time a script here runs, rather than persisted
# into MatchTable.csv/AUC_summary.json by the batch pipeline itself (which
# would require re-running matching for every dataset -- not done yet).
#
# Only DeepUnitMatch/UMPy model folders ever save a UMparam.pickle (written
# unconditionally by UnitMatchPy.save_utils.save_to_output, called from both
# branches of run_deepunitmatch_batch_onMerged.py) -- EMD/DANT/DANT_no_functional
# never do, so FR_diff_norm is silently unavailable for those models until
# the pipeline itself is changed to compute and persist it directly (at which
# point this on-the-fly path becomes unnecessary and can be removed).
FR_DIFF_SCORE = "FR_diff"
FR_DIFF_NORM_SCORE = "FR_diff_norm"


def load_unit_fr_lookup(model_dir):
    """
    Load model_dir/UMparam.pickle and return {(RecSes (1-based), unit ID):
    mean firing rate}, built from DeepUnitMatch.testing.test.get_FR(param) --
    the same cross-validated per-unit firing rate fast_testing.normalize_FR_diff()
    divides FR_diff by.

    Keyed off param["good_units"]'s own per-session cluster-ID lists (with a
    running offset) rather than assuming MatchTable.csv row order matches
    concatenation order: UnitMatchPy.save_utils.make_match_table() builds
    ID1/RecSes 1 from exactly this same
    "np.concatenate(param['good_units'])" + per-unit session index
    construction (see clus_info["original_ids"]/["session_id"]), which is
    also exactly get_FR's own internal per-unit indexing order -- so this
    mapping is derived the same way the table itself was, not inferred from
    table contents.

    Returns None if UMparam.pickle is missing (e.g. EMD/DANT model folders --
    see the module comment above) or fails to load/compute -- callers should
    skip FR_diff_norm for that model/dataset in that case, not treat it as an
    error.
    """
    pickle_path = os.path.join(model_dir, "UMparam.pickle")
    if not os.path.isfile(pickle_path):
        return None
    try:
        with open(pickle_path, "rb") as f:
            param = pickle.load(f)
        FR = dumtest.get_FR(param)  # (2, n_units), fold x unit
    except Exception as e:
        print(f"  WARNING: could not load/compute FR from {pickle_path}: {e}")
        return None

    fr_mean = FR.mean(axis=0)  # fold-averaged, same as fast_testing.normalize_FR_diff
    lookup = {}
    offset = 0
    for session_idx, gu in enumerate(param["good_units"]):
        ids = np.asarray(gu).flatten()
        for pos, uid in enumerate(ids):
            lookup[(session_idx + 1, int(uid))] = fr_mean[offset + pos]
        offset += len(ids)
    return lookup


def add_fr_diff_norm_column(df, model_dir, diff_col=FR_DIFF_SCORE, recses_col="RecSes 1", id_col="ID1"):
    """
    Adds an "FR_diff_norm" column to df (a MatchTable.csv-shaped DataFrame
    already containing diff_col), normalizing diff_col by the source unit's
    own mean firing rate -- see load_unit_fr_lookup()/the module comment
    above for why this is computed on the fly here rather than read from a
    persisted column.

    Returns True if the column was added, False if UMparam.pickle wasn't
    available for model_dir -- callers should treat False as "skip
    FR_diff_norm for this model", not an error.
    """
    if diff_col not in df.columns:
        return False
    fr_lookup = load_unit_fr_lookup(model_dir)
    if fr_lookup is None:
        return False

    fr_df = pd.DataFrame(
        [(recses, uid, fr) for (recses, uid), fr in fr_lookup.items()],
        columns=[recses_col, id_col, "_fr_tmp"],
    )
    merged_fr = df[[recses_col, id_col]].merge(fr_df, on=[recses_col, id_col], how="left")["_fr_tmp"]
    df[FR_DIFF_NORM_SCORE] = (df[diff_col].to_numpy() / merged_fr.to_numpy()).astype("float32")
    return True


def compute_fr_diff_norm_auc(model_dir):
    """
    On-the-fly counterpart to the scalar FR_diff AUC that AUC_summary.json
    already stores: recomputes the same across-session-pairs AUC (same
    recall/FPR/trapezoid logic as DeepUnitMatch.testing.test.AUC(), just over
    flat MatchTable.csv rows instead of the square (n_units, n_units) matrix
    -- duplicated here rather than imported from plot_auc_vs_delta_days.py's
    own _auc_from_flat() to avoid a circular import, since that module
    already imports this one), but on FR_diff normalized by each source
    unit's own mean firing rate (see add_fr_diff_norm_column()) instead of
    the raw difference.

    Returns None if MatchTable.csv/UMparam.pickle aren't available, or if
    there's nothing to rank (no across-session match or non-match) --
    mirrors auc_summary_from_functional_scores()'s own try/except AUC()
    handling: the caller should just omit the score, not treat this as an
    error.
    """
    match_table_path = os.path.join(model_dir, "MatchTable.csv")
    if not os.path.isfile(match_table_path):
        return None

    header = pd.read_csv(match_table_path, nrows=0).columns
    if FR_DIFF_SCORE not in header:
        return None

    usecols = ["ID1", "RecSes 1", "RecSes 2", "Matches", FR_DIFF_SCORE]
    dtype_map = {
        "ID1": "int32", "RecSes 1": "int32", "RecSes 2": "int32",
        "Matches": "int8", FR_DIFF_SCORE: "float32",
    }
    df = pd.read_csv(match_table_path, usecols=usecols, dtype=dtype_map)
    df = df[df["RecSes 1"] != df["RecSes 2"]]

    if not add_fr_diff_norm_column(df, model_dir, diff_col=FR_DIFF_SCORE):
        return None

    matches_bool = df["Matches"].to_numpy().astype(bool)
    # FR_diff_norm is still "lower = more similar" like FR_diff itself (a
    # ratio built from the same difference) -- negate for the higher-is-
    # better AUC convention, same sign NEGATE_FOR_AUC applies to FR_diff.
    metric = -df[FR_DIFF_NORM_SCORE].to_numpy()
    valid = np.isfinite(metric)
    matches_bool, metric = matches_bool[valid], metric[valid]

    P = int(matches_bool.sum())
    N = len(matches_bool) - P
    if P < 1 or N < 1:
        return None

    order = np.argsort(metric)[::-1]
    m_sorted = matches_bool[order]
    tp = np.cumsum(m_sorted)
    fp = np.cumsum(~m_sorted)
    recall = tp / P
    fpr = fp / N
    if hasattr(np, "trapezoid"):
        return float(np.trapezoid(recall, fpr))
    return float(np.trapz(recall, fpr))


def savefig_with_svg(fig, out_path, **kwargs):
    """Save out_path (a .png path) plus an editable .svg copy alongside it (e.g. for Inkscape)."""
    fig.savefig(out_path, **kwargs)
    fig.savefig(os.path.splitext(out_path)[0] + ".svg", **{k: v for k, v in kwargs.items() if k != "dpi"})


# ── collect ──────────────────────────────────────────────────────────────────

# Subtrees that never contain an AUC_summary.json themselves but sit alongside
# folders that do, and can be large (many small files/folders) -- an unpruned
# os.walk descends into these for no reason, which is expensive on a network
# share. DANT_shared/_dant_input/spike_times holds one .npy per unit; EMD's
# own output folder holds _stage/<session>/ and one result_<i>_<j>/ per
# session pair (grows with n_sessions choose 2) right next to its own
# AUC_summary.json.
_SKIP_DIR_NAMES = {"DANT_shared", "_stage"}
_SKIP_DIR_PREFIXES = ("result_",)


def _should_skip_dir(name):
    return (
        name in _SKIP_DIR_NAMES
        or name.startswith(".")
        or any(name.startswith(p) for p in _SKIP_DIR_PREFIXES)
    )


def mouse_from_dataset(dataset):
    """First path component of the BASE_OUTPUT-relative dataset path, e.g.
    'AL032/19011111882/2' -> 'AL032'."""
    return dataset.replace("\\", "/").split("/")[0]


def collect_auc_summaries(base_output=BASE_OUTPUT, models_to_include=MODELS_TO_INCLUDE):
    """
    Walk base_output for every AUC_summary.json and return a long-form
    DataFrame with columns: mouse, dataset, model, score, value.

    Prunes the walk aggressively (see _should_skip_dir / the "found it,
    stop descending" below) since base_output is typically a network share
    and a naive full recursive walk re-visits many large subtrees (EMD's
    per-session-pair result folders, DANT's per-unit spike-time cache, ...)
    that can never contain an AUC_summary.json of their own.
    """
    rows = []
    n_found = 0
    for root, dirs, files in os.walk(base_output, topdown=True):
        dirs[:] = [d for d in dirs if not _should_skip_dir(d)]

        if "AUC_summary.json" not in files:
            continue
        dirs[:] = []  # this is a leaf model dir -- don't descend further (e.g. EMD's own result_i_j/_stage)

        model = os.path.basename(root)
        if models_to_include is not None and model not in models_to_include:
            continue
        dataset = os.path.relpath(os.path.dirname(root), base_output)
        mouse = mouse_from_dataset(dataset)

        with open(os.path.join(root, "AUC_summary.json")) as f:
            try:
                summary = json.load(f)
            except json.JSONDecodeError:
                print(f"  WARNING: could not parse {root}\\AUC_summary.json, skipping.")
                continue

        n_found += 1
        if n_found % 25 == 0:
            print(f"  ...found {n_found} AUC_summary.json files so far", flush=True)

        if FR_DIFF_SCORE in summary and FR_DIFF_NORM_SCORE not in summary:
            fr_norm_auc = compute_fr_diff_norm_auc(root)
            if fr_norm_auc is not None:
                summary = {**summary, FR_DIFF_NORM_SCORE: fr_norm_auc}

        for score, value in summary.items():
            rows.append(
                {"mouse": mouse, "dataset": dataset, "model": model, "score": score, "value": value}
            )

    if not rows:
        raise RuntimeError(f"No AUC_summary.json files found under {base_output}")

    return pd.DataFrame(rows)


def get_passing_datasets(df_long, min_matches=MIN_MATCHES_TO_INCLUDE):
    """
    Datasets where at least one included model's n_matches_across_sessions
    reaches min_matches -- see MIN_MATCHES_TO_INCLUDE.
    """
    n_matches = df_long[df_long["score"] == "n_matches_across_sessions"]
    best = n_matches.groupby("dataset")["value"].max()
    return set(best[best >= min_matches].index)


def filter_low_match_datasets(df_long, min_matches=MIN_MATCHES_TO_INCLUDE):
    """
    Drop every row belonging to a dataset that fails get_passing_datasets()
    -- not just its n_matches_across_sessions rows -- so a dataset no model
    could find enough matches in doesn't contribute a spuriously extreme
    point to any score's comparison plot either.
    """
    passing = get_passing_datasets(df_long, min_matches)
    all_datasets = set(df_long["dataset"].unique())
    dropped = all_datasets - passing
    if dropped:
        print(
            f"Dropping {len(dropped)}/{len(all_datasets)} dataset(s) with < {min_matches} "
            f"matches in every included model: {sorted(dropped)}"
        )
    return df_long[df_long["dataset"].isin(passing)]


def compute_matching_failures(df_long):
    """
    Per model, count datasets where matching effectively failed: either the
    model produced no AUC_summary.json at all for that dataset (crashed or
    was skipped before save_utils.save_auc_summary ran), or it ran but found
    zero across-session matches (n_matches_across_sessions == 0 -- always
    present whenever AUC_summary.json exists, since
    test.auc_summary_from_functional_scores computes it unconditionally,
    unlike the per-score AUCs which are silently *omitted* from the JSON on
    failure rather than recorded as 0/NaN -- see collect_auc_summaries()).

    A model that fails outright on a dataset never shows up as a bad value in
    the per-score plots above (it's just absent from that score's points), so
    it's reported here separately instead of being folded into the AUC
    comparison.

    Returns a DataFrame with columns: model, n_no_output, n_zero_matches,
    n_failed, n_datasets, frac_failed.
    """
    all_datasets = set(df_long["dataset"].unique())
    n_total = len(all_datasets)
    rows = []
    for model in sorted(df_long["model"].unique()):
        attempted = df_long[(df_long["model"] == model) & (df_long["score"] == "n_matches_across_sessions")]
        attempted_datasets = set(attempted["dataset"])
        no_output = all_datasets - attempted_datasets
        zero_matches = set(attempted.loc[attempted["value"] == 0, "dataset"])
        rows.append(
            {
                "model": model,
                "n_no_output": len(no_output),
                "n_zero_matches": len(zero_matches),
                "n_failed": len(no_output) + len(zero_matches),
                "n_datasets": n_total,
                "frac_failed": (len(no_output) + len(zero_matches)) / n_total if n_total else np.nan,
            }
        )
    return pd.DataFrame(rows)


def plot_matching_failures(fail_df, output_dir):
    """
    Stacked bar per model: datasets where the model produced no output at all
    (crashed/skipped, bottom colour flipped to be on top for visibility)
    stacked on datasets where it ran but found zero across-session matches --
    see compute_matching_failures().
    """
    if fail_df.empty or fail_df["n_failed"].sum() == 0:
        return None

    models = fail_df["model"].tolist()
    x = np.arange(len(models))
    n_total = fail_df["n_datasets"].iloc[0]

    fig, ax = plt.subplots(figsize=(max(6, 1.2 * len(models)), 5))
    ax.bar(x, fail_df["n_zero_matches"], label="ran, zero matches found", color="tab:orange")
    ax.bar(x, fail_df["n_no_output"], bottom=fail_df["n_zero_matches"], label="no output produced", color="tab:red")

    for xi, n_failed in zip(x, fail_df["n_failed"]):
        if n_failed > 0:
            ax.text(xi, n_failed + 0.02 * n_total, str(int(n_failed)), ha="center", va="bottom", fontsize=9)

    ax.set_xticks(x)
    ax.set_xticklabels(models, rotation=30, ha="right")
    ax.set_ylabel(f"# datasets with failed matching (of {n_total})")
    ax.set_title("Matching failures per model")
    ax.grid(axis="y", alpha=0.3)
    ax.legend(fontsize=8)
    fig.tight_layout()

    out_path = os.path.join(output_dir, "matching_failures.png")
    savefig_with_svg(fig, out_path, dpi=150)
    plt.close(fig)
    return out_path


def compute_n_matches_vs_auc_summary(df_long):
    """
    Per (score, model): mean +/- s.d. of AUC, plus a *geometric* mean +/- s.d.
    of n_matches_across_sessions, across that model's datasets, and the
    dataset count each is over. Companion table for plot_n_matches_vs_auc()'s
    one-point-per-model summary (see there for why this is averaged rather
    than one point per dataset).

    n_matches_across_sessions is heavily right-skewed (a handful of datasets
    can have orders of magnitude more matches than the rest -- the same
    reason LOG_TRANSFORM_FOR_STATS fits mixed models on log1p(value) for it
    elsewhere in this script): a linear mean/s.d. is dominated by those
    outlier datasets and produces error bars wide enough to span almost the
    whole log axis, which defeats the point of averaging. Averaging
    log1p(n_matches) and mapping back with expm1 instead gives a "typical
    dataset" centre and a spread that reflects multiplicative variation.
    """
    scores = sorted(s for s in df_long["score"].unique() if s != "n_matches_across_sessions")
    n_matches = df_long[df_long["score"] == "n_matches_across_sessions"][["mouse", "dataset", "model", "value"]]
    n_matches = n_matches.rename(columns={"value": "n_matches"})
    n_matches["log_n_matches"] = np.log1p(n_matches["n_matches"])

    rows = []
    for score in scores:
        sub = df_long[df_long["score"] == score][["mouse", "dataset", "model", "value"]]
        sub = sub.merge(n_matches, on=["mouse", "dataset", "model"], how="inner")
        sub = sub.dropna(subset=["value", "n_matches"])

        auc = sub.groupby("model")["value"].agg(auc_mean="mean", auc_std="std")
        log_n = sub.groupby("model")["log_n_matches"].agg(log_mean="mean", log_std="std")
        n_datasets = sub.groupby("model").size().rename("n_datasets")
        summary = pd.concat([auc, log_n, n_datasets], axis=1)

        # A model with only one dataset for this score has an undefined
        # (NaN) s.d. -- treat that as zero spread rather than propagating
        # NaN into the plotted error bar.
        log_std = summary["log_std"].fillna(0.0)
        summary["n_matches_geomean"] = np.expm1(summary["log_mean"])
        summary["n_matches_geo_low"] = np.expm1(summary["log_mean"] - log_std).clip(lower=0.0)
        summary["n_matches_geo_high"] = np.expm1(summary["log_mean"] + log_std)
        summary = summary.drop(columns=["log_mean", "log_std"])

        summary.insert(0, "score", score)
        rows.append(summary.reset_index())

    if not rows:
        return pd.DataFrame(
            columns=[
                "score", "model", "auc_mean", "auc_std",
                "n_matches_geomean", "n_matches_geo_low", "n_matches_geo_high", "n_datasets",
            ]
        )
    return pd.concat(rows, ignore_index=True)


def plot_n_matches_vs_auc(summary_df, output_dir):
    """
    One figure, one subplot per functional score: x = geometric mean
    n_matches_across_sessions (log scale), y = mean AUC for that score, one
    point per model with error bars showing +/- 1 s.d. across that model's
    datasets (log-space for x, linear for y -- see
    compute_n_matches_vs_auc_summary() for why x uses a geometric mean/s.d.).
    Averaged rather than one point per dataset (which is unreadable once
    there are hundreds of datasets) so match yield and match quality stay
    readable together per model -- e.g. a model that only reaches a high AUC
    by finding very few (easy) matches sits at low x, whereas a model doing
    well on both axes sits at high x and high y.

    A model's mean here is implicitly conditioned on "datasets where it
    found at least one match" (see compute_matching_failures() for the
    zero-match/no-output failure modes, which contribute no data point at
    all here).
    """
    if summary_df.empty:
        return None

    scores = sorted(summary_df["score"].unique())
    models = sorted(summary_df["model"].unique())
    cmap = plt.get_cmap("tab10" if len(models) <= 10 else "tab20")
    colour_for = {m: cmap(i % cmap.N) for i, m in enumerate(models)}

    ncols = min(3, len(scores))
    nrows = int(np.ceil(len(scores) / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(5 * ncols, 4.5 * nrows), squeeze=False)

    for i, score in enumerate(scores):
        ax = axes[i // ncols][i % ncols]
        sub = summary_df[summary_df["score"] == score]
        for _, row in sub.iterrows():
            x, y = row["n_matches_geomean"], row["auc_mean"]
            y_std = row["auc_std"] if np.isfinite(row["auc_std"]) else 0.0
            xerr = [[x - row["n_matches_geo_low"]], [row["n_matches_geo_high"] - x]]
            ax.errorbar(
                x, y, xerr=xerr, yerr=y_std,
                fmt="o", color=colour_for[row["model"]], ecolor=colour_for[row["model"]],
                elinewidth=1, capsize=3, markersize=7, markeredgecolor="k", markeredgewidth=0.5,
                label=row["model"],
            )
        ax.set_xscale("log")
        ax.set_xlabel("n_matches_across_sessions (geometric mean ± s.d.)")
        ax.set_ylabel(f"AUC ({score})")
        ax.set_title(score)
        ax.grid(alpha=0.3)

    for j in range(len(scores), nrows * ncols):
        axes[j // ncols][j % ncols].set_visible(False)

    handles = [
        plt.Line2D(
            [0], [0], marker="o", linestyle="", color=colour_for[m],
            label=m, markeredgecolor="k", markeredgewidth=0.3,
        )
        for m in models
    ]
    fig.legend(handles=handles, title="model", bbox_to_anchor=(1.0, 0.5), loc="center left", fontsize=8)
    fig.tight_layout()

    out_path = os.path.join(output_dir, "n_matches_vs_auc.png")
    savefig_with_svg(fig, out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return out_path


def average_over_mice(df_long):
    """
    Collapse dataset-level rows to one value per (mouse, model, score) by
    averaging across that mouse's datasets. This removes the pseudoreplication
    from unequal numbers of datasets per mouse, at the cost of losing
    within-mouse dataset-to-dataset variability -- the "simple" companion to
    the dataset-level mixed-effects analysis, where mouse no longer needs to
    be modelled as a random effect because each mouse now contributes exactly
    one point per model.
    """
    return df_long.groupby(["mouse", "model", "score"], as_index=False)["value"].mean()


# ── significance annotation ─────────────────────────────────────────────────


def _p_to_stars(p):
    if p < 0.001:
        return "***"
    if p < 0.01:
        return "**"
    if p < ALPHA:
        return "*"
    return None


def add_significance_brackets(ax, model_order, pvals_adj, y_values):
    """
    Draw a horizontal bracket + star label between each pair of models (by
    x-position) whose corrected p-value is below ALPHA, stacked bottom-to-top
    by increasing span so overlapping brackets don't collide.

    pvals_adj: dict {frozenset({model_a, model_b}): p_adj}.
    """
    if not pvals_adj:
        return

    sig_pairs = []
    for pair, p in pvals_adj.items():
        stars = _p_to_stars(p)
        if stars is None:
            continue
        m1, m2 = tuple(pair)
        if m1 not in model_order or m2 not in model_order:
            continue
        i1, i2 = sorted((model_order.index(m1), model_order.index(m2)))
        sig_pairs.append((i1, i2, stars))

    if not sig_pairs:
        return

    sig_pairs.sort(key=lambda t: t[1] - t[0])  # narrower spans stacked lower

    y_values = np.asarray(y_values, dtype=float)
    y_values = y_values[~np.isnan(y_values)]
    y_top = y_values.max() if y_values.size else 1.0
    y_bottom = y_values.min() if y_values.size else 0.0
    y_range = (y_top - y_bottom) or abs(y_top) or 1.0
    step = 0.09 * y_range
    base = y_top + 0.08 * y_range

    for level, (i1, i2, stars) in enumerate(sig_pairs):
        y = base + level * step
        tick = step * 0.15
        ax.plot([i1, i1, i2, i2], [y, y + tick, y + tick, y], color="black", linewidth=1, clip_on=False)
        ax.text((i1 + i2) / 2, y + tick, stars, ha="center", va="bottom", fontsize=11, clip_on=False)

    ax.set_ylim(top=base + len(sig_pairs) * step + step)


# ── plotting ─────────────────────────────────────────────────────────────────


def plot_score(df_long, score, output_dir, pvals_adj=None):
    """
    One figure per score: x = model, y = value, one point per dataset
    (jittered horizontally so overlapping datasets stay visible), coloured
    by mouse, with a black bar marking each model's mean, and significance
    brackets for any pairwise comparison in pvals_adj below ALPHA.
    """
    sub = df_long[df_long["score"] == score].dropna(subset=["value"])
    if sub.empty:
        return None

    models = sorted(sub["model"].unique())
    mice = sorted(sub["mouse"].unique())
    cmap = plt.get_cmap("tab10" if len(mice) <= 10 else "tab20")
    colour_for = {m: cmap(i % cmap.N) for i, m in enumerate(mice)}

    rng = np.random.default_rng(0)
    fig, ax = plt.subplots(figsize=(max(6, 1.2 * len(models)), 5))
    for xi, model in enumerate(models):
        model_rows = sub[sub["model"] == model]
        jitter = rng.uniform(-0.15, 0.15, size=len(model_rows))
        ax.scatter(
            xi + jitter,
            model_rows["value"],
            c=[colour_for[m] for m in model_rows["mouse"]],
            s=35,
            alpha=0.85,
            edgecolors="k",
            linewidths=0.3,
            zorder=3,
        )
        ax.hlines(model_rows["value"].mean(), xi - 0.25, xi + 0.25, colors="black", linewidth=2, zorder=4)

    ax.set_xticks(range(len(models)))
    ax.set_xticklabels(models, rotation=30, ha="right")
    ax.set_ylabel(score)
    ax.set_title(score)
    ax.grid(axis="y", alpha=0.3)

    handles = [
        plt.Line2D(
            [0], [0], marker="o", linestyle="", color=colour_for[m],
            label=m, markeredgecolor="k", markeredgewidth=0.3,
        )
        for m in mice
    ]
    ax.legend(handles=handles, title="mouse", bbox_to_anchor=(1.02, 1), loc="upper left", fontsize=8)

    add_significance_brackets(ax, models, pvals_adj, sub["value"].values)
    fig.tight_layout()

    out_path = os.path.join(output_dir, f"{score}.png")
    savefig_with_svg(fig, out_path, dpi=150)
    plt.close(fig)
    return out_path


def plot_score_mouse_avg(df_mouse, score, output_dir, pvals_adj=None):
    """
    Mouse-averaged companion to plot_score: one point per mouse (averaged
    over that mouse's datasets), thin lines connecting each mouse's points
    across models to make the paired comparisons visible, and significance
    brackets for any pairwise comparison in pvals_adj below ALPHA.
    """
    sub = df_mouse[df_mouse["score"] == score].dropna(subset=["value"])
    if sub.empty:
        return None

    models = sorted(sub["model"].unique())
    mice = sorted(sub["mouse"].unique())
    cmap = plt.get_cmap("tab10" if len(mice) <= 10 else "tab20")
    colour_for = {m: cmap(i % cmap.N) for i, m in enumerate(mice)}

    fig, ax = plt.subplots(figsize=(max(6, 1.2 * len(models)), 5))

    wide = sub.pivot(index="mouse", columns="model", values="value")
    for mouse in wide.index:
        ys = [wide.loc[mouse, m] if m in wide.columns else np.nan for m in models]
        xs = list(range(len(models)))
        ax.plot(xs, ys, color=colour_for[mouse], alpha=0.5, linewidth=1, zorder=2)
        ax.scatter(xs, ys, color=[colour_for[mouse]] * len(xs), s=45, edgecolors="k", linewidths=0.4, zorder=3)

    for xi, model in enumerate(models):
        vals = sub[sub["model"] == model]["value"]
        ax.hlines(vals.mean(), xi - 0.25, xi + 0.25, colors="black", linewidth=2, zorder=4)

    ax.set_xticks(range(len(models)))
    ax.set_xticklabels(models, rotation=30, ha="right")
    ax.set_ylabel(score)
    ax.set_title(f"{score} (averaged per mouse)")
    ax.grid(axis="y", alpha=0.3)

    handles = [
        plt.Line2D(
            [0], [0], marker="o", linestyle="", color=colour_for[m],
            label=m, markeredgecolor="k", markeredgewidth=0.3,
        )
        for m in mice
    ]
    ax.legend(handles=handles, title="mouse", bbox_to_anchor=(1.02, 1), loc="upper left", fontsize=8)

    add_significance_brackets(ax, models, pvals_adj, sub["value"].values)
    fig.tight_layout()

    out_path = os.path.join(output_dir, f"{score}_mouse_averaged.png")
    savefig_with_svg(fig, out_path, dpi=150)
    plt.close(fig)
    return out_path


# ── statistics: dataset-level (mixed model) ─────────────────────────────────


def fit_mixed_model(df_long, score, reference_model=REFERENCE_MODEL):
    """
    Linear mixed-effects model: value ~ model (fixed) + (1 | mouse) (random
    intercept per mouse), fit with statsmodels' MixedLM. A random intercept
    per mouse absorbs consistent per-animal offsets so the fixed-effect
    comparison reflects differences between models rather than which mice
    happened to contribute more datasets to which model.

    Returns the fitted result, or None if there are fewer than 2 models or
    fewer than 2 mice with data for this score.
    """
    sub = df_long[df_long["score"] == score].dropna(subset=["value"]).copy()
    if sub["model"].nunique() < 2 or sub["mouse"].nunique() < 2:
        return None

    if score in LOG_TRANSFORM_FOR_STATS:
        sub["value"] = np.log1p(sub["value"])

    ref = reference_model if reference_model in sub["model"].unique() else sorted(sub["model"].unique())[0]
    other_levels = sorted(m for m in sub["model"].unique() if m != ref)
    sub["model"] = pd.Categorical(sub["model"], categories=[ref] + other_levels)

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        try:
            result = smf.mixedlm("value ~ C(model)", sub, groups=sub["mouse"]).fit()
        except Exception as e:
            print(f"  WARNING: mixed model failed for '{score}': {e}")
            return None
    return result


def pairwise_mixed_pvalues(df_long, score):
    """
    Post-hoc: refit the mixed model above for every pair of models (each
    pair using only mice/datasets that have that pair) to get a p-value for
    that specific pairwise difference, then Holm-Bonferroni correct across
    all pairs computed for this score. Returns {frozenset({m1, m2}): p_adj}.
    """
    sub = df_long[df_long["score"] == score].dropna(subset=["value"]).copy()
    models = sorted(sub["model"].unique())
    if len(models) < 2:
        return {}

    if score in LOG_TRANSFORM_FOR_STATS:
        sub["value"] = np.log1p(sub["value"])

    raw_pvals = {}
    for m1, m2 in combinations(models, 2):
        pair_sub = sub[sub["model"].isin([m1, m2])].copy()
        if pair_sub["mouse"].nunique() < 2:
            continue
        pair_sub["model"] = pd.Categorical(pair_sub["model"], categories=[m1, m2])
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            try:
                result = smf.mixedlm("value ~ C(model)", pair_sub, groups=pair_sub["mouse"]).fit()
            except Exception:
                continue
        coef_names = [c for c in result.pvalues.index if c not in ("Intercept", "Group Var")]
        if not coef_names:
            continue
        raw_pvals[frozenset((m1, m2))] = result.pvalues[coef_names[0]]

    return _holm_correct(raw_pvals)


def _holm_correct(raw_pvals):
    """dict {pair: p_raw} -> dict {pair: p_adj} via Holm-Bonferroni across all pairs."""
    if not raw_pvals:
        return {}
    pairs = list(raw_pvals.keys())
    pvals = [raw_pvals[p] for p in pairs]
    _, p_adj, _, _ = multipletests(pvals, method="holm")
    return dict(zip(pairs, p_adj))


# ── statistics: mouse-averaged (paired t-test) ──────────────────────────────


def fit_paired_ttest_vs_reference(df_mouse, score, reference_model=REFERENCE_MODEL):
    """
    Simple companion to fit_mixed_model for the mouse-averaged data: a paired
    t-test of each other model against reference_model, paired on mouse (only
    mice with both models present are used). No random-effects machinery
    needed here since each mouse already contributes exactly one point per
    model.

    Returns a report string, or None if fewer than 2 models have data.
    """
    sub = df_mouse[df_mouse["score"] == score].dropna(subset=["value"]).copy()
    if sub["model"].nunique() < 2:
        return None

    if score in LOG_TRANSFORM_FOR_STATS:
        sub["value"] = np.log1p(sub["value"])

    ref = reference_model if reference_model in sub["model"].unique() else sorted(sub["model"].unique())[0]
    wide = sub.pivot(index="mouse", columns="model", values="value")
    if ref not in wide.columns:
        return None

    lines = [f"Paired t-test vs reference model '{ref}' (paired on mouse):"]
    for model in sorted(c for c in wide.columns if c != ref):
        paired = wide[[ref, model]].dropna()
        if len(paired) < 2:
            lines.append(f"  {model}: skipped (fewer than 2 mice with both models present)")
            continue
        stat, p = ttest_rel(paired[model], paired[ref])
        mean_diff = (paired[model] - paired[ref]).mean()
        lines.append(f"  {model}: mean diff = {mean_diff:+.4f}, t = {stat:.3f}, p = {p:.4g}, n_mice = {len(paired)}")
    return "\n".join(lines)


def pairwise_paired_ttest_pvalues(df_mouse, score):
    """
    Post-hoc for the mouse-averaged data: paired t-test for every pair of
    models (paired on mouse, using only mice with both models present),
    Holm-Bonferroni corrected across all pairs computed for this score.
    Returns {frozenset({m1, m2}): p_adj}.
    """
    sub = df_mouse[df_mouse["score"] == score].dropna(subset=["value"]).copy()
    models = sorted(sub["model"].unique())
    if len(models) < 2:
        return {}

    if score in LOG_TRANSFORM_FOR_STATS:
        sub["value"] = np.log1p(sub["value"])

    wide = sub.pivot(index="mouse", columns="model", values="value")

    raw_pvals = {}
    for m1, m2 in combinations(models, 2):
        if m1 not in wide.columns or m2 not in wide.columns:
            continue
        paired = wide[[m1, m2]].dropna()
        if len(paired) < 2:
            continue
        _, p = ttest_rel(paired[m1], paired[m2])
        raw_pvals[frozenset((m1, m2))] = p

    return _holm_correct(raw_pvals)


def _format_posthoc_lines(pvals_adj):
    lines = ["  Post-hoc pairwise comparisons (Holm-Bonferroni corrected across pairs):"]
    if not pvals_adj:
        lines.append("    (skipped: fewer than 2 comparable models/mice)")
        return lines
    for pair, p in sorted(pvals_adj.items(), key=lambda kv: kv[1]):
        m1, m2 = tuple(pair)
        stars = _p_to_stars(p) or "ns"
        lines.append(f"    {m1} vs {m2}: p_adj = {p:.4g} ({stars})")
    return lines


# ── main ─────────────────────────────────────────────────────────────────────


def main():
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    print(f"Scanning for AUC_summary.json under:\n  {BASE_OUTPUT}\n")
    df_long = collect_auc_summaries()

    n_datasets = df_long[["mouse", "dataset"]].drop_duplicates().shape[0]
    n_mice = df_long["mouse"].nunique()
    models_found = sorted(df_long["model"].unique())
    print(
        f"Found {len(df_long)} score values across {n_datasets} datasets, "
        f"{n_mice} mice, {len(models_found)} models: {models_found}\n"
    )

    csv_path = os.path.join(OUTPUT_DIR, "auc_summary_long.csv")
    df_long.to_csv(csv_path, index=False)
    print(f"Wrote long-form table to {csv_path}")

    fail_df = compute_matching_failures(df_long)
    fail_csv_path = os.path.join(OUTPUT_DIR, "matching_failures.csv")
    fail_df.to_csv(fail_csv_path, index=False)
    print(f"Wrote matching-failure counts to {fail_csv_path}")
    fail_png = plot_matching_failures(fail_df, OUTPUT_DIR)
    if fail_png:
        print(f"Plotted matching failures -> {fail_png}")

    # Everything below (n_matches-vs-AUC, score comparisons, mixed models,
    # mouse-averaged stats) uses only datasets where at least one included
    # model found enough matches to trust -- see filter_low_match_datasets().
    # compute_matching_failures() above ran on the full set so it still
    # reports on datasets this drops.
    df_long = filter_low_match_datasets(df_long)

    n_matches_auc_summary = compute_n_matches_vs_auc_summary(df_long)
    n_matches_auc_csv = os.path.join(OUTPUT_DIR, "n_matches_vs_auc_summary.csv")
    n_matches_auc_summary.to_csv(n_matches_auc_csv, index=False)
    print(f"Wrote n_matches-vs-AUC summary to {n_matches_auc_csv}")
    n_matches_auc_png = plot_n_matches_vs_auc(n_matches_auc_summary, OUTPUT_DIR)
    if n_matches_auc_png:
        print(f"Plotted n_matches vs AUC -> {n_matches_auc_png}")

    df_mouse = average_over_mice(df_long)
    mouse_csv_path = os.path.join(OUTPUT_DIR, "auc_summary_mouse_averaged.csv")
    df_mouse.to_csv(mouse_csv_path, index=False)
    print(f"Wrote mouse-averaged table to {mouse_csv_path}")

    # ── dataset-level: mixed model + post-hoc pairwise ──
    stats_lines = []
    for score in sorted(df_long["score"].unique()):
        pvals_adj = pairwise_mixed_pvalues(df_long, score)
        png_path = plot_score(df_long, score, OUTPUT_DIR, pvals_adj=pvals_adj)
        if png_path:
            print(f"  Plotted {score} -> {png_path}")

        result = fit_mixed_model(df_long, score)
        header = f"\n{'=' * 70}\n{score}"
        if score in LOG_TRANSFORM_FOR_STATS:
            header += " (fit on log1p(value))"
        stats_lines.append(header + f"\n{'=' * 70}")
        if result is None:
            stats_lines.append("  (skipped: fewer than 2 models or fewer than 2 mice with data)")
            continue
        stats_lines.append(str(result.summary()))
        stats_lines.append("")
        stats_lines.extend(_format_posthoc_lines(pvals_adj))

    stats_path = os.path.join(OUTPUT_DIR, "mixed_model_summaries.txt")
    with open(stats_path, "w") as f:
        f.write("\n".join(stats_lines))
    print(f"\nWrote dataset-level (mixed-model) stats to {stats_path}")

    # ── mouse-averaged: paired t-test + post-hoc pairwise ──
    mouse_stats_lines = []
    for score in sorted(df_mouse["score"].unique()):
        pvals_adj_mouse = pairwise_paired_ttest_pvalues(df_mouse, score)
        png_path = plot_score_mouse_avg(df_mouse, score, OUTPUT_DIR, pvals_adj=pvals_adj_mouse)
        if png_path:
            print(f"  Plotted {score} (mouse-averaged) -> {png_path}")

        header = f"\n{'=' * 70}\n{score}"
        if score in LOG_TRANSFORM_FOR_STATS:
            header += " (fit on log1p(value))"
        mouse_stats_lines.append(header + f"\n{'=' * 70}")
        simple_summary = fit_paired_ttest_vs_reference(df_mouse, score)
        if simple_summary is None:
            mouse_stats_lines.append("  (skipped: fewer than 2 models with data)")
            continue
        mouse_stats_lines.append(simple_summary)
        mouse_stats_lines.append("")
        mouse_stats_lines.extend(_format_posthoc_lines(pvals_adj_mouse))

    mouse_stats_path = os.path.join(OUTPUT_DIR, "mouse_averaged_stats.txt")
    with open(mouse_stats_path, "w") as f:
        f.write("\n".join(mouse_stats_lines))
    print(f"Wrote mouse-averaged stats to {mouse_stats_path}")


if __name__ == "__main__":
    main()
