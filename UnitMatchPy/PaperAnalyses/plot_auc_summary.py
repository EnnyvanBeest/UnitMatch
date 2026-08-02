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

import run_deepunitmatch_batch_onMerged as base_batch

# ── settings ─────────────────────────────────────────────────────────────────
# Point this at whichever BASE_OUTPUT tree you want to summarise (defaults to
# the merged-data comparison tree, where DUM/UMPy/EMD/DANT/sweeps all live).
BASE_OUTPUT = base_batch.BASE_OUTPUT
OUTPUT_DIR = os.path.join(BASE_OUTPUT, "auc_summary_report")

# Restrict to specific model/method subfolder names (e.g.
# {"DeepUnitMatch", "UMPy", "EMD", "DANT", "DANT_no_functional"}), or None to
# include every AUC_summary.json found under BASE_OUTPUT (main models,
# maxdist sweep, extra checkpoints, everything).
MODELS_TO_INCLUDE = None
#MODELS_TO_INCLUDE = {"DeepUnitMatch","DeepUnitMatch_AssignUniqueID_Conservative", "DeepUnitMatch_AssignUniqueID", "UMPy","UMPy_AssignUniqueID_Conservative", "UMPy_AssignUniqueID", "EMD", "DANT", "DANT_no_functional"}
#MODELS_TO_INCLUDE = {"DeepUnitMatch","n_output=256_after_ae_and_finetune","n_output=128_after_ae_and_finetune","n_output=32_after_ae_and_finetune","n_output=8_after_ae_and_finetune"}
# Reference level for the "vs reference" summaries (every other model's
# coefficient/mean-difference is "difference from this model"). Falls back to
# whichever model sorts first for a given score if this one has no data there.
REFERENCE_MODEL = "EMD"

# Counts are heavily right-skewed; fit stats on log1p(count) instead of the
# raw value (plots still show raw counts).
LOG_TRANSFORM_FOR_STATS = {"n_matches_across_sessions"}

ALPHA = 0.05


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

        for score, value in summary.items():
            rows.append(
                {"mouse": mouse, "dataset": dataset, "model": model, "score": score, "value": value}
            )

    if not rows:
        raise RuntimeError(f"No AUC_summary.json files found under {base_output}")

    return pd.DataFrame(rows)


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
    fig.savefig(out_path, dpi=150)
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
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
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
    fig.savefig(out_path, dpi=150)
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
    fig.savefig(out_path, dpi=150)
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
