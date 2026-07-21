# Batch wrapper: runs DeepUnitMatch on the *merged* dataset (see
# run_deepunitmatch_batch_onMerged.py for how that tree is built and how
# sessions/good units are derived) once per "extra" trained DNN checkpoint
# under DeepUnitMatch/ExtraModels, instead of the single default model.
#
# Two experiment families are discovered under DeepUnitMatch/ExtraModels:
#   loss_sensitivity_experiment/W_ij=<x>/ckpt_epoch_*
#       -> saved to BASE_OUTPUT/<dataset>/DUM_W_ij=<x>
#   n_output_experiment/n_output=<x>/{after_ae,after_ae_and_finetune}/ckpt_epoch_*
#       -> saved to BASE_OUTPUT/<dataset>/n_output=<x>_after_ae
#          and      BASE_OUTPUT/<dataset>/n_output=<x>_after_ae_and_finetune
# ('n_output=256-chinesecharacters' is a different training dataset, not a
# point on the n_output sweep, and is skipped.)
#
# Each such folder sits alongside the DeepUnitMatch/ and UMPy/ subfolders that
# run_deepunitmatch_batch_onMerged.py writes for the same dataset, using the
# exact same inference/matching/saving pipeline (run_deep_unit_match_core),
# just against a different checkpoint.

import os
import sys
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

import run_deepunitmatch_batch_onMerged as base_batch
from DeepUnitMatch.testing import test

EXTRA_MODELS_ROOT = os.path.join(
    os.path.dirname(_HERE), "DeepUnitMatch", "ExtraModels"
)
EXCLUDED_N_OUTPUT_DIRS = {"256-chinesecharacters"}

REDO = False  # if True, rerun even when output already exists


# ── extra-model discovery ────────────────────────────────────────────────────


def discover_extra_models():
    """
    Scan DeepUnitMatch/ExtraModels for trained checkpoints.

    Returns a list of dicts: {"checkpoint": path, "n_output": int, "subfolder_name": str}
    """
    models = []

    loss_root = Path(EXTRA_MODELS_ROOT) / "loss_sensitivity_experiment"
    for d in sorted(loss_root.glob("W_ij=*")):
        if not d.is_dir():
            continue
        x = d.name.split("=", 1)[1]
        ckpts = sorted(d.glob("ckpt_epoch_*"))
        if not ckpts:
            print(f"  WARNING: no checkpoint found in {d}, skipping.")
            continue
        models.append(
            {
                "checkpoint": str(ckpts[0]),
                "n_output": 256,
                "subfolder_name": f"DUM_W_ij={x}",
            }
        )

    nout_root = Path(EXTRA_MODELS_ROOT) / "n_output_experiment"
    for d in sorted(nout_root.glob("n_output=*")):
        if not d.is_dir():
            continue
        x = d.name.split("=", 1)[1]
        if x in EXCLUDED_N_OUTPUT_DIRS:
            print(f"  Skipping {d.name} (excluded).")
            continue
        try:
            n_output = int(x)
        except ValueError:
            print(f"  WARNING: could not parse n_output from {d.name}, skipping.")
            continue

        for stage in ("after_ae", "after_ae_and_finetune"):
            stage_dir = d / stage
            if not stage_dir.is_dir():
                continue
            ckpts = sorted(stage_dir.glob("ckpt_epoch_*"))
            if not ckpts:
                print(f"  WARNING: no checkpoint found in {stage_dir}, skipping.")
                continue
            models.append(
                {
                    "checkpoint": str(ckpts[0]),
                    "n_output": n_output,
                    "subfolder_name": f"n_output={x}_{stage}",
                }
            )

    return models


# ── path helpers ─────────────────────────────────────────────────────────────


def get_extra_save_dir(merged_dir, model_info):
    """Output dir for a given merged-data group + extra model, mirroring get_save_dir/get_umpy_save_dir."""
    subfolder = os.path.relpath(os.path.dirname(merged_dir), base_batch.BASE_INPUT)
    return os.path.join(base_batch.BASE_OUTPUT, subfolder, model_info["subfolder_name"])


def extra_results_exist(merged_dir, model_info):
    """Return True when the sentinel output file is present for this dataset/model combo."""
    sentinel = os.path.join(get_extra_save_dir(merged_dir, model_info), "MatchingOverview.png")
    return os.path.isfile(sentinel)


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
        description="Run DeepUnitMatch on the merged dataset with each extra trained model."
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

    extra_models = discover_extra_models()
    if not extra_models:
        print("No extra models found under:\n  " + EXTRA_MODELS_ROOT)
        return
    print(f"Found {len(extra_models)} extra model(s):")
    for m in extra_models:
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
            m for m in extra_models if REDO or not extra_results_exist(merged_dir, m)
        ]
        if not pending:
            print("  Skipping all extra models (results exist, REDO=False).")
            continue

        sess = base_batch._prepare_session(merged_dir)
        if sess is None:
            continue

        for model_info in pending:
            save_dir = get_extra_save_dir(merged_dir, model_info)
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
