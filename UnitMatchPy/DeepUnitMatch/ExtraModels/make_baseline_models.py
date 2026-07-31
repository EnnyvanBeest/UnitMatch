# One-off copier for the three baseline checkpoints used to sanity-check the
# trained DeepUnitMatch models in run_deepunitmatch_batch_onMerged_extramodels.py.
#
# Source: the canonical DeepUnitMatch training archive
# (ModelExp/{AE_experiments,experiments}/<exp_name>/ckpt/...), e.g.
# .../Dropbox/DeepUnitMatch_Data/Data_June2026/ModelExp -- the same archive
# that produced the production default model (verified byte-identical to
# ModelExp/experiments/original/ckpt/model == DeepUnitMatch/utils/model).
#
#   untrained_baseline/ckpt_epoch_0
#       = ModelExp/AE_experiments/untrained/ckpt/ckpt_epoch_0
#       Randomly initialized encoder, evaluated directly for matching --
#       never trained on anything ("fully untrained").
#
#   unfinetuned_baseline/ckpt_epoch_299
#       = ModelExp/AE_experiments/unfinetuned/ckpt/ckpt_epoch_299
#       Autoencoder-pretrained for 300 epochs, evaluated directly for
#       matching -- never went through the CLIP-loss finetune stage.
#
#   finetuned_only_baseline/ckpt_epoch_49
#       = ModelExp/experiments/untrainedAE/ckpt/ckpt_epoch_49
#       Started from a randomly initialized encoder (no AE pretraining --
#       ModelExp/AE_experiments/untrainedAE/ckpt/ckpt_epoch_0, a *different*
#       random init than untrained_baseline above) and CLIP-loss finetuned
#       for 50 epochs via train/train_finetune.py's --finetune-omitted /
#       from_scratch=True path: same freeze policy as the production model
#       (only FcBlock/projector/clip_loss are trained; the conv backbone
#       stays at its random init throughout).
#
# All three were already produced (not by this script) as part of the
# original DeepUnitMatch training run family. Re-run this only if the
# archive moves or the checkpoints need refreshing:
#   python make_baseline_models.py --archive-root "<path to the ModelExp parent dir>"

import argparse
import shutil
import sys
from pathlib import Path

_HERE = Path(__file__).resolve().parent

DEFAULT_ARCHIVE_ROOT = (
    r"C:\Users\EnnyB\Dropbox\DeepUnitMatch_Data\Data_June2026"
)

# (relative path under <archive_root>/ModelExp, destination dir name, destination filename)
BASELINES = [
    ("AE_experiments/untrained/ckpt/ckpt_epoch_0", "untrained_baseline", "ckpt_epoch_0"),
    ("AE_experiments/unfinetuned/ckpt/ckpt_epoch_299", "unfinetuned_baseline", "ckpt_epoch_299"),
    ("experiments/untrainedAE/ckpt/ckpt_epoch_49", "finetuned_only_baseline", "ckpt_epoch_49"),
]


def copy_baselines(archive_root):
    model_exp = Path(archive_root) / "ModelExp"
    for rel_src, dest_dir, dest_name in BASELINES:
        src = model_exp / rel_src
        if not src.exists():
            raise FileNotFoundError(f"Expected source checkpoint not found: {src}")

        out_dir = _HERE / dest_dir
        out_dir.mkdir(exist_ok=True)
        out_path = out_dir / dest_name
        shutil.copyfile(src, out_path)
        print(f"  Copied {src} -> {out_path}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Copy the untrained/unfinetuned/finetuned_only baseline "
        "checkpoints from the DeepUnitMatch training archive into ExtraModels."
    )
    parser.add_argument(
        "--archive-root",
        default=DEFAULT_ARCHIVE_ROOT,
        help="Directory containing the ModelExp/ tree (default: %(default)s).",
    )
    args = parser.parse_args()

    print(f"Copying baseline checkpoints from: {args.archive_root}")
    copy_baselines(args.archive_root)
    print("Done.")
