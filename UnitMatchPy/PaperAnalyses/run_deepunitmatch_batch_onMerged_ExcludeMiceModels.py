# Batch wrapper: runs DeepUnitMatch on the *merged* dataset (see
# run_deepunitmatch_batch_onMerged.py for how that tree is built and how
# sessions/good units are derived) once per trained checkpoint under
# DeepUnitMatch/ExtraModels/exclude_mice_exp, instead of the single default
# model.
#
# That folder holds a "leave-some-mice-out" training experiment: for each
# training-set size N in {1, 6, 12} mice, several independently trained models
# (model_m<N>_<idx>) live under
#   exclude_mice_exp/<N>_mice/model_m<N>_<idx>/{after_ae,after_ae_and_finetune}/ckpt_epoch_*
# mirroring the two-stage (AE pretraining, then CLIP finetuning) checkpoints
# used elsewhere in DeepUnitMatch/ExtraModels and DeepUnitMatch/BaselineModels.
# Results are saved to
#   BASE_OUTPUT/<dataset>/exclude_mice_m<N>_<idx>_after_ae
#   BASE_OUTPUT/<dataset>/exclude_mice_m<N>_<idx>_after_ae_and_finetune
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

EXCLUDE_MICE_ROOT = os.path.join(
    os.path.dirname(_HERE), "DeepUnitMatch", "ExtraModels", "exclude_mice_exp"
)


# See batch_lock.sentinel_is_fresh() / run_deepunitmatch_batch_onMerged.py's
# REDO_FROM_DATE for what this does: a dataset/model combo is skipped once its
# MatchingOverview.png exists and is at least this new. None falls back to
# plain "skip if present" (appropriate here since these models are new and no
# prior run needs to be superseded); a far-future date would reproduce old
# REDO=True behaviour.
REDO_FROM_DATE = None
RUN_UNFINETUNED_MODELS = False  # if False, skip the after_ae (not fine-tuned) checkpoints


# ── model discovery ──────────────────────────────────────────────────────────


def discover_exclude_mice_models():
    """
    Scan DeepUnitMatch/ExtraModels/exclude_mice_exp for trained checkpoints.

    Returns a list of dicts: {"checkpoint": path, "n_output": int, "subfolder_name": str}
    """
    models = []

    root = Path(EXCLUDE_MICE_ROOT)
    for mice_dir in sorted(root.glob("*_mice")):
        if not mice_dir.is_dir():
            continue

        for model_dir in sorted(mice_dir.glob("model_*")):
            if not model_dir.is_dir():
                continue
            # "model_m12_1" -> "m12_1"
            model_name = model_dir.name.split("model_", 1)[-1]

            for stage in ("after_ae", "after_ae_and_finetune"):
                if stage == "after_ae" and not RUN_UNFINETUNED_MODELS:
                    continue
                stage_dir = model_dir / stage
                if not stage_dir.is_dir():
                    continue
                ckpts = sorted(stage_dir.glob("ckpt_epoch_*"))
                if not ckpts:
                    print(f"  WARNING: no checkpoint found in {stage_dir}, skipping.")
                    continue
                models.append(
                    {
                        "checkpoint": str(ckpts[0]),
                        "n_output": 256,
                        "subfolder_name": f"exclude_mice_{model_name}_{stage}",
                    }
                )

    return models


# ── path helpers ─────────────────────────────────────────────────────────────


def get_exclude_mice_save_dir(merged_dir, model_info):
    """Output dir for a given merged-data group + exclude-mice model, mirroring get_save_dir/get_umpy_save_dir."""
    subfolder = os.path.relpath(os.path.dirname(merged_dir), base_batch.BASE_INPUT)
    return os.path.join(base_batch.BASE_OUTPUT, subfolder, model_info["subfolder_name"])


def exclude_mice_results_exist(merged_dir, model_info):
    """Return True when the sentinel output file is present and fresh for this dataset/model combo (see REDO_FROM_DATE)."""
    sentinel = os.path.join(get_exclude_mice_save_dir(merged_dir, model_info), "MatchingOverview.png")
    return batch_lock.sentinel_is_fresh(sentinel, REDO_FROM_DATE)


def get_group_lock_path(merged_dir):
    """
    Lock file marking 'a run is currently processing all pending exclude-mice
    models for this group', so multiple machines pointed at the same
    BASE_INPUT/BASE_OUTPUT can split work across groups without
    double-processing one. See batch_lock.py. Named distinctly from the other
    batch scripts' locks so they can all run on the same group concurrently.
    """
    subfolder = os.path.relpath(os.path.dirname(merged_dir), base_batch.BASE_INPUT)
    return os.path.join(base_batch.BASE_OUTPUT, subfolder, ".processing_exclude_mice_models.lock")


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
        description="Run DeepUnitMatch on the merged dataset with each exclude-mice trained model."
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

    exclude_mice_models = discover_exclude_mice_models()
    if not exclude_mice_models:
        print("No exclude-mice models found under:\n  " + EXCLUDE_MICE_ROOT)
        return
    print(f"Found {len(exclude_mice_models)} exclude-mice model(s):")
    for m in exclude_mice_models:
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
            m for m in exclude_mice_models if not exclude_mice_results_exist(merged_dir, m)
        ]
        if not pending:
            print("  Skipping all exclude-mice models (results exist and are fresh).")
            continue

        lock_path = get_group_lock_path(merged_dir)
        with batch_lock.try_lock(lock_path) as acquired:
            if not acquired:
                print(f"  Skipping (already being processed by another run): {lock_path}")
                continue

            # re-check now that we hold the lock: another machine may have
            # finished this group while we were scanning/waiting for the lock
            pending = [
                m for m in exclude_mice_models if not exclude_mice_results_exist(merged_dir, m)
            ]
            if not pending:
                print("  Skipping all exclude-mice models (completed by another run).")
                continue

            sess = base_batch._prepare_session(merged_dir)
            if sess is None:
                continue

            for model_info in pending:
                save_dir = get_exclude_mice_save_dir(merged_dir, model_info)
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
