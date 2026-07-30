# One-off generator for the two baseline checkpoints used to sanity-check the
# trained DeepUnitMatch models in run_deepunitmatch_batch_onMerged_extramodels.py:
#
#   untrained_baseline/ckpt_epoch_0
#       SpatioTemporalCNN_V2(n_output=256) with only its default
#       (Kaiming/LayerNorm) initialization -- never seen any data. Same
#       checkpoint shape as an autoencoder-only checkpoint (see
#       testing/test.py:load_trained_model), so it loads via the "encoder"
#       branch (no "clip_loss" key).
#
#   unfinetuned_baseline/ckpt_epoch_299
#       Copy of ExtraModels/n_output_experiment/n_output=256/after_ae/ckpt_epoch_299
#       -- the n_output=256 encoder after autoencoder pretraining but *before*
#       the CLIP-loss finetune step that produces the production default
#       model (DeepUnitMatch/utils/model). Kept as its own top-level folder
#       (rather than relying on the n_output sweep) so it's always
#       discovered regardless of RUN_UNFINETUNED_N_OUTPUT_MODELS.
#
# Run once: python make_baseline_models.py

import shutil
import sys
from pathlib import Path

import torch

_HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(_HERE.parent))  # DeepUnitMatch/ root, for `utils`

from utils.mymodel import SpatioTemporalCNN_V2

N_OUTPUT = 256  # matches DeepUnitMatch/utils/model, the production default
SEED = 0


def make_untrained_checkpoint():
    out_dir = _HERE / "untrained_baseline"
    out_dir.mkdir(exist_ok=True)
    out_path = out_dir / "ckpt_epoch_0"
    if out_path.exists():
        print(f"  Already exists, skipping: {out_path}")
        return

    torch.manual_seed(SEED)
    model = SpatioTemporalCNN_V2(n_channel=30, n_time=60, n_output=N_OUTPUT).double()
    torch.save({"encoder": model.state_dict(), "epoch": 0}, out_path)
    print(f"  Wrote {out_path}")


def make_unfinetuned_checkpoint():
    src = (
        _HERE
        / "n_output_experiment"
        / f"n_output={N_OUTPUT}"
        / "after_ae"
        / "ckpt_epoch_299"
    )
    out_dir = _HERE / "unfinetuned_baseline"
    out_dir.mkdir(exist_ok=True)
    out_path = out_dir / "ckpt_epoch_299"
    if out_path.exists():
        print(f"  Already exists, skipping: {out_path}")
        return
    if not src.exists():
        raise FileNotFoundError(f"Expected source checkpoint not found: {src}")

    shutil.copyfile(src, out_path)
    print(f"  Copied {src} -> {out_path}")


if __name__ == "__main__":
    print("Creating untrained baseline checkpoint...")
    make_untrained_checkpoint()
    print("Creating unfinetuned baseline checkpoint...")
    make_unfinetuned_checkpoint()
    print("Done.")
