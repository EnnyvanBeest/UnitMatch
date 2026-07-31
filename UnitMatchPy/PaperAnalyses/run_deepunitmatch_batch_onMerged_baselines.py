# Batch wrapper: runs DeepUnitMatch on the *merged* dataset (see
# run_deepunitmatch_batch_onMerged.py for how that tree is built and how
# sessions/good units are derived) once per baseline checkpoint under
# DeepUnitMatch/BaselineModels, instead of the single default model.
#
# These three checkpoints (see BaselineModels/make_baseline_models.py for
# provenance) isolate what each of the two training stages (AE pretraining,
# CLIP finetuning) contributes on its own, relative to the production
# default model (which goes through both):
#   untrained_baseline/ckpt_epoch_*        -> BASE_OUTPUT/<dataset>/DUM_untrained
#       Random weights, evaluated as-is -- neither stage.
#   unfinetuned_baseline/ckpt_epoch_*      -> BASE_OUTPUT/<dataset>/DUM_unfinetuned
#       AE-pretrained only, never CLIP-finetuned.
#   finetuned_only_baseline/ckpt_epoch_*   -> BASE_OUTPUT/<dataset>/DUM_finetuned_only
#       CLIP-finetuned on top of a frozen random backbone -- never AE-pretrained.
#
# Each such folder sits alongside the DeepUnitMatch/ and UMPy/ subfolders that
# run_deepunitmatch_batch_onMerged.py writes for the same dataset, using the
# exact same inference/matching/saving pipeline (run_deep_unit_match_core),
# just against a different checkpoint.

import os
import sys
import datetime
import traceback
import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")  # non-interactive backend for batch runs

# ── project paths ───────────────────────────────────────────────────────────
_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _HERE)
sys.path.insert(0, os.path.dirname(_HERE))
sys.path.insert(0, os.path.join(_HERE, "DeepUnitMatch"))

import batch_lock
import run_deepunitmatch_batch_onMerged as base_batch
from DeepUnitMatch.testing import test

BASELINE_MODELS_ROOT = os.path.join(
    os.path.dirname(_HERE), "DeepUnitMatch", "BaselineModels"
)

# (folder name under BaselineModels, output subfolder name under BASE_OUTPUT)
BASELINE_DIRS = (
    ("untrained_baseline", "DUM_untrained"),
    ("unfinetuned_baseline", "DUM_unfinetuned"),
    ("finetuned_only_baseline", "DUM_finetuned_only"),
)

# See batch_lock.sentinel_is_fresh() / run_deepunitmatch_batch_onMerged.py's
# REDO_FROM_DATE for what this does: a dataset/model combo is skipped once its
# MatchingOverview.png exists and is at least this new. None falls back to
# plain "skip if present"; a far-future date reproduces old REDO=True.
REDO_FROM_DATE = datetime.datetime(2026, 7, 22, 19, 0, 0)


# ── baseline-model discovery ─────────────────────────────────────────────────


def discover_baseline_models():
    """
    Scan DeepUnitMatch/BaselineModels for the untrained/unfinetuned/
    finetuned_only checkpoints.

    Returns a list of dicts: {"checkpoint": path, "n_output": int, "subfolder_name": str}
    """
    models = []

    for baseline_dir, subfolder_name in BASELINE_DIRS:
        d = Path(BASELINE_MODELS_ROOT) / baseline_dir
        if not d.is_dir():
            print(f"  WARNING: {d} not found, skipping.")
            continue
        ckpts = sorted(d.glob("ckpt_epoch_*"))
        if not ckpts:
            print(f"  WARNING: no checkpoint found in {d}, skipping.")
            continue
        models.append(
            {
                "checkpoint": str(ckpts[0]),
                "n_output": 256,
                "subfolder_name": subfolder_name,
            }
        )

    return models


# ── path helpers ─────────────────────────────────────────────────────────────


def get_baseline_save_dir(merged_dir, model_info):
    """Output dir for a given merged-data group + baseline model, mirroring get_save_dir/get_umpy_save_dir."""
    subfolder = os.path.relpath(os.path.dirname(merged_dir), base_batch.BASE_INPUT)
    return os.path.join(base_batch.BASE_OUTPUT, subfolder, model_info["subfolder_name"])


def baseline_results_exist(merged_dir, model_info):
    """Return True when the sentinel output file is present and fresh for this dataset/model combo (see REDO_FROM_DATE)."""
    sentinel = os.path.join(get_baseline_save_dir(merged_dir, model_info), "MatchingOverview.png")
    return batch_lock.sentinel_is_fresh(sentinel, REDO_FROM_DATE)


def get_group_lock_path(merged_dir):
    """
    Lock file marking 'a run is currently processing all pending baseline
    models for this group', so multiple machines pointed at the same
    BASE_INPUT/BASE_OUTPUT can split work across groups without
    double-processing one. See batch_lock.py. Named distinctly from the other
    batch scripts' locks so they can all run on the same group concurrently.
    """
    subfolder = os.path.relpath(os.path.dirname(merged_dir), base_batch.BASE_INPUT)
    return os.path.join(base_batch.BASE_OUTPUT, subfolder, ".processing_baselines.lock")


# ── model cache ───────────────────────────────────────────────────────────────


def get_model_for_checkpoint(cache, model_info):
    """Load (and cache) the model for a given checkpoint so it's not reloaded per dataset."""
    key = model_info["checkpoint"]
    if key not in cache:
        print(
            f"Loading model: {model_info['subfolder_name']} "
            f"(n_output={model_info['n_output']}, checkpoint={key}) …"
        )
        cache[key] = test.load_trained_model(
            device=base_batch.DEVICE,
            read_path=key,
            n_output=model_info["n_output"],
        )
    return cache[key]


# ── entry point ───────────────────────────────────────────────────────────────


def parse_args():
    parser = argparse.ArgumentParser(
        description="Run DeepUnitMatch on the merged dataset with each baseline (untrained/unfinetuned/finetuned_only) model."
    )
    parser.add_argument(
        "--write-matlab-compat",
        action="store_true",
        help="Also write a MATLAB-compatible UnitMatch.mat from the Python outputs.",
    )
    return parser.parse_args()


def main():
    args = parse_args()
    base_batch.WRITE_MATLAB_COMPAT = args.write_matlab_compat

    baseline_models = discover_baseline_models()
    if not baseline_models:
        print("No baseline models found under:\n  " + BASELINE_MODELS_ROOT)
        return
    print(f"Found {len(baseline_models)} baseline model(s):")
    for m in baseline_models:
        print(f"  {m['subfolder_name']}  (n_output={m['n_output']})  <- {m['checkpoint']}")

    print(f"\nScanning for merged-data groups under:\n  {base_batch.BASE_INPUT}\n")
    groups = base_batch.find_merged_groups()
    if not groups:
        print("No merged-data groups found.")
        return
    print(f"Found {len(groups)} group(s).\n")

    model_cache = {}

    for i, merged_dir in enumerate(groups):
        print(f"\n[{i + 1}/{len(groups)}] {merged_dir}")

        pending = [
            m for m in baseline_models if not baseline_results_exist(merged_dir, m)
        ]
        if not pending:
            print("  Skipping all baseline models (results exist and are fresh).")
            continue

        lock_path = get_group_lock_path(merged_dir)
        with batch_lock.try_lock(lock_path) as acquired:
            if not acquired:
                print(f"  Skipping (already being processed by another run): {lock_path}")
                continue

            # re-check now that we hold the lock: another machine may have
            # finished this group while we were scanning/waiting for the lock
            pending = [
                m for m in baseline_models if not baseline_results_exist(merged_dir, m)
            ]
            if not pending:
                print("  Skipping all baseline models (completed by another run).")
                continue

            sess = base_batch._prepare_session(merged_dir)
            if sess is None:
                continue

            for model_info in pending:
                save_dir = get_baseline_save_dir(merged_dir, model_info)
                try:
                    model = get_model_for_checkpoint(model_cache, model_info)
                    base_batch.run_deep_unit_match_core(
                        sess, save_dir, model, label=model_info["subfolder_name"]
                    )
                except Exception as e:
                    print(f"  {model_info['subfolder_name']} FAILED: {e}")
                    traceback.print_exc()

    print("\nAll done.")


if __name__ == "__main__":
    main()
