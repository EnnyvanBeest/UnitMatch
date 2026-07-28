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

import os
import sys
import json
import warnings

import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")  # non-interactive backend
import matplotlib.pyplot as plt
import statsmodels.formula.api as smf

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _HERE)

import run_deepunitmatch_batch_onMerged as base_batch

# ── settings ─────────────────────────────────────────────────────────────────
# Point this at whichever BASE_OUTPUT tree you want to summarise (defaults to
# the merged-data comparison tree, where DUM/UMPy/EMD/DANT/sweeps all live).
BASE_OUTPUT = base_batch.BASE_OUTPUT
OUTPUT_DIR = os.path.join(_HERE, "auc_summary_report")

# Restrict to specific model/method subfolder names (e.g.
# {"DeepUnitMatch", "UMPy", "EMD", "DANT", "DANT_no_functional"}), or None to
# include every AUC_summary.json found under BASE_OUTPUT (main models,
# maxdist sweep, extra checkpoints, everything).
MODELS_TO_INCLUDE = {"DeepUnitMatch", "UMPy", "EMD", "DANT", "DANT_no_functional"}

# Reference level for the mixed model (every other model's coefficient is
# then "difference from this model, averaged over mice"). Falls back to
# whichever model sorts first for a given score if this one has no data there.
REFERENCE_MODEL = "UMPy"

# Counts are heavily right-skewed; fit the mixed model on log1p(count)
# instead of the raw value (plots still show raw counts).
LOG_TRANSFORM_FOR_STATS = {"n_matches_across_sessions"}


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


# ── plotting ─────────────────────────────────────────────────────────────────


def plot_score(df_long, score, output_dir):
    """
    One figure per score: x = model, y = value, one point per dataset
    (jittered horizontally so overlapping datasets stay visible), coloured
    by mouse, with a black bar marking each model's mean.
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
    fig.tight_layout()

    out_path = os.path.join(output_dir, f"{score}.png")
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    return out_path


# ── statistics ───────────────────────────────────────────────────────────────


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

    stats_lines = []
    for score in sorted(df_long["score"].unique()):
        png_path = plot_score(df_long, score, OUTPUT_DIR)
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

    stats_path = os.path.join(OUTPUT_DIR, "mixed_model_summaries.txt")
    with open(stats_path, "w") as f:
        f.write("\n".join(stats_lines))
    print(f"\nWrote mixed-effects model summaries to {stats_path}")


if __name__ == "__main__":
    main()
