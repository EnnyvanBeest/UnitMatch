# Batch wrapper: runs DeepUnitMatch on the *merged* dataset (see
# run_deepunitmatch_batch_onMerged.py for how that tree is built and how
# sessions/good units are derived) using the new checkpoint trained by
# train_deepunitmatch_from_merged.py (architecture changes: symmetrized/
# correctly-weighted CLIP loss, trainable backbone during finetuning, real
# per-channel probe geometry, amplitude/noise augmentation -- see that
# script and utils/mymodel.py's ChannelPositionalBias for details), instead
# of the original production checkpoint at DeepUnitMatch/utils/model.
#
# Mirrors run_deepunitmatch_batch_onMerged_extramodels.py's pattern (reuse
# find_merged_groups/_prepare_session/run_deep_unit_match_core from the base
# script, save into a distinctly-named subfolder alongside the DeepUnitMatch/
# and UMPy/ output the base script already writes for the same dataset) but
# simplified to a single model instead of a discovery loop over many, since
# there is exactly one new checkpoint here. UMPy output is untouched -- it
# doesn't depend on which DNN checkpoint is used, and the base script has
# presumably already produced it for these groups.

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
from DeepUnitMatch.testing import test_v2 as test

# The exp_name train_deepunitmatch_from_merged.py was run with -- finetune
# checkpoints live under DeepUnitMatch/ModelExp/experiments/<EXP_NAME>/ckpt/.
EXP_NAME = "NewModelAug2026"
CKPT_DIR = os.path.join(
    os.path.dirname(_HERE), "DeepUnitMatch", "ModelExp", "experiments", EXP_NAME, "ckpt"
)
N_OUTPUT = 256
SUBFOLDER_NAME = f"DUM_{EXP_NAME}"

# See batch_lock.sentinel_is_fresh() / run_deepunitmatch_batch_onMerged.py's
# REDO_FROM_DATE for what this does. None (the default here) means "skip a
# group once its sentinel exists" -- this model's output is new, so there's
# no stale pre-fix run to redo yet. Set to a specific date if you retrain
# this same EXP_NAME again later and want old output redone.
REDO_FROM_DATE = None


def latest_checkpoint(ckpt_dir):
    """Most recent ckpt_epoch_N file in ckpt_dir, same convention train_AE.py/
    train_finetune.py use for --cont."""
    ckpts = [f for f in os.listdir(ckpt_dir) if f.startswith("ckpt_epoch_")]
    if not ckpts:
        raise FileNotFoundError(f"No ckpt_epoch_* files found in {ckpt_dir}")
    ckpts.sort(key=lambda x: int(x.split("_")[-1]))
    return os.path.join(ckpt_dir, ckpts[-1])


# ── path helpers ─────────────────────────────────────────────────────────────


def get_save_dir(merged_dir):
    """Output dir for a given merged-data group with this new model, mirroring get_save_dir/get_umpy_save_dir."""
    subfolder = os.path.relpath(os.path.dirname(merged_dir), base_batch.BASE_INPUT)
    return os.path.join(base_batch.BASE_OUTPUT, subfolder, SUBFOLDER_NAME)


def results_exist(merged_dir):
    """Return True when the sentinel output file is present (and fresh, see REDO_FROM_DATE)."""
    sentinel = os.path.join(get_save_dir(merged_dir), "MatchingOverview.png")
    return batch_lock.sentinel_is_fresh(sentinel, REDO_FROM_DATE)


def get_group_lock_path(merged_dir):
    """
    Lock file marking 'a run is currently processing this group with the new
    model', so multiple machines pointed at the same BASE_INPUT/BASE_OUTPUT
    can split work across groups without double-processing one. Named
    distinctly from the base/extramodels scripts' locks so all three can run
    on the same group concurrently without blocking each other.
    """
    subfolder = os.path.relpath(os.path.dirname(merged_dir), base_batch.BASE_INPUT)
    return os.path.join(base_batch.BASE_OUTPUT, subfolder, f".processing_{SUBFOLDER_NAME}.lock")


# ── entry point ───────────────────────────────────────────────────────────────


def parse_args():
    parser = argparse.ArgumentParser(
        description="Run DeepUnitMatch on the merged dataset with the new trained model."
    )
    parser.add_argument(
        "--checkpoint",
        default=None,
        help=f"Explicit checkpoint path (default: latest ckpt_epoch_* under {CKPT_DIR})",
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

    checkpoint = args.checkpoint or latest_checkpoint(CKPT_DIR)
    print(f"Using checkpoint: {checkpoint}")
    print("Loading model …")
    model = test.load_trained_model(
        device=base_batch.DEVICE, read_path=checkpoint, n_output=N_OUTPUT
    )

    print(f"\nScanning for merged-data groups under:\n  {base_batch.BASE_INPUT}\n")
    groups = base_batch.find_merged_groups()
    if not groups:
        print("No merged-data groups found.")
        return
    print(f"Found {len(groups)} group(s).\n")

    for i, merged_dir in enumerate(groups):
        print(f"\n[{i + 1}/{len(groups)}] {merged_dir}")

        if results_exist(merged_dir):
            print("  Skipping (results exist and are fresh).")
            continue

        lock_path = get_group_lock_path(merged_dir)
        with batch_lock.try_lock(lock_path) as acquired:
            if not acquired:
                print(f"  Skipping (already being processed by another run): {lock_path}")
                continue

            # re-check now that we hold the lock: another machine may have
            # finished this group while we were scanning/waiting for the lock
            if results_exist(merged_dir):
                print("  Skipping (completed by another run).")
                continue

            sess = base_batch._prepare_session(merged_dir)
            if sess is None:
                continue

            save_dir = get_save_dir(merged_dir)
            try:
                base_batch.run_deep_unit_match_core(
                    sess, save_dir, model, label=SUBFOLDER_NAME
                )
            except Exception as e:
                print(f"  {SUBFOLDER_NAME} FAILED: {e}")
                traceback.print_exc()

    print("\nAll done.")


if __name__ == "__main__":
    main()
