# Wrapper: preprocess the full merged-data corpus and train a new
# DeepUnitMatch AE + finetuned model from ALL locations found in it.
#
# Read-only against BASE_INPUT (the existing merged_data_v2 raw KS-style
# tree already used for inference by run_deepunitmatch_batch_onMerged.py).
# Preprocessed snippets (param_fun.get_snippets output) are written to
# PREPROCESSED_OUTPUT, a NEW sibling folder -- this script never writes into
# or deletes anything under merged_data_v2 or any other existing
# preprocessed-output folder. Model checkpoints go to the local
# DeepUnitMatch/ModelExp/ folder, matching train_AE.py / train_finetune.py's
# own conventions, and utils/model (the production checkpoint) is never
# touched -- promoting a new checkpoint to production is a deliberate,
# separate step (see utils/model_PreJuly2026 for how the previous production
# checkpoint was preserved before this training run).
#
# Stages (run independently via --stage, default is to run all of them in
# order):
#   preprocess     For every merged-data location, run the same
#                  param_fun.get_snippets call the inference pipeline
#                  already uses per-comparison, except the output is
#                  persisted under PREPROCESSED_OUTPUT instead of being
#                  discarded in a temp dir afterward. Idempotent: a location
#                  with a ".done" sentinel is skipped unless --redo is
#                  passed, so an interrupted run can simply be re-launched.
#   train_ae       Reconstruction-pretrain the encoder over every
#                  preprocessed location pooled together.
#   train_finetune Contrastive-finetune the pretrained encoder over every
#                  preprocessed location pooled together.
#
# Example:
#   python train_deepunitmatch_from_merged.py --stage preprocess --limit 1   # smoke test
#   python train_deepunitmatch_from_merged.py --stage preprocess             # full corpus
#   python train_deepunitmatch_from_merged.py --stage train_ae
#   python train_deepunitmatch_from_merged.py --stage train_finetune

import os
import sys
import argparse
import numpy as np
from torch.utils.data import ConcatDataset

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _HERE)
sys.path.insert(0, os.path.dirname(_HERE))
sys.path.insert(0, os.path.join(_HERE, "DeepUnitMatch"))

import run_deepunitmatch_batch_onMerged as batch_pipeline  # reuses find_merged_groups / _prepare_session
from DeepUnitMatch.utils import param_fun_v2 as param_fun
from DeepUnitMatch.utils.AE_npdataset import AE_NeuropixelsDataset
from DeepUnitMatch.utils.npdataset_v2 import NeuropixelsDataset_cortexlab
from DeepUnitMatch.train import train_AE as train_ae_mod
from DeepUnitMatch.train import train_finetunev2 as train_finetune_mod

BASE_INPUT = batch_pipeline.BASE_INPUT  # raw KS-style merged tree -- read-only
# New sibling folder for this training run's preprocessed snippets. Never
# reuses merged_data_v2 (the raw input) or any other existing output folder.
PREPROCESSED_OUTPUT = os.path.join(
    os.path.dirname(BASE_INPUT), "merged_data_vNewModelAug2026"
)

DEFAULT_EXP_NAME = "NewModelAug2026"


def get_location_out_dir(merged_dir):
    """Mirror BASE_INPUT's subfolder layout under PREPROCESSED_OUTPUT."""
    subfolder = os.path.relpath(os.path.dirname(merged_dir), BASE_INPUT)
    return os.path.join(PREPROCESSED_OUTPUT, subfolder)


def preprocess_all(limit=None, redo=False):
    groups = batch_pipeline.find_merged_groups()
    if limit:
        groups = groups[:limit]
    print(f"Found {len(groups)} merged-data location(s) under {BASE_INPUT}")

    n_done, n_skipped, n_failed = 0, 0, 0
    for i, merged_dir in enumerate(groups):
        out_dir = get_location_out_dir(merged_dir)
        sentinel = os.path.join(out_dir, ".done")
        print(f"\n[{i + 1}/{len(groups)}] {merged_dir}")

        if os.path.exists(sentinel) and not redo:
            print(f"  Already preprocessed -> {out_dir} (skipping; pass --redo to reprocess)")
            n_skipped += 1
            continue

        sess = batch_pipeline._prepare_session(merged_dir)
        if sess is None:
            print("  SKIPPING (failed to load raw session data)")
            n_failed += 1
            continue

        unit_ids = np.concatenate(sess["param"]["good_units"]).squeeze()
        os.makedirs(out_dir, exist_ok=True)
        try:
            # Same call the inference pipeline already makes per-comparison
            # (run_deepunitmatch_batch_onMerged.run_deep_unit_match_core);
            # the only difference is out_dir is a persistent, per-location
            # folder under PREPROCESSED_OUTPUT instead of a temp dir that
            # gets discarded after one matching run.
            param_fun.get_snippets(
                sess["waveform"],
                sess["channel_pos"],
                sess["session_id"],
                save_path=out_dir,
                unit_ids=unit_ids,
                param=sess["param"],
            )
        except Exception as e:
            print(f"  ERROR in get_snippets: {e}")
            n_failed += 1
            continue

        with open(sentinel, "w") as f:
            f.write("ok")
        n_done += 1

    print(
        f"\nPreprocessing complete: {n_done} newly done, {n_skipped} already done, "
        f"{n_failed} failed."
    )


def list_preprocessed_locations():
    """Every location folder under PREPROCESSED_OUTPUT with a .done sentinel."""
    locations = []
    if not os.path.isdir(PREPROCESSED_OUTPUT):
        return locations
    for root, _, files in os.walk(PREPROCESSED_OUTPUT):
        if ".done" in files:
            locations.append(root)
    return sorted(locations)


class MultiLocationFinetuneDataset(NeuropixelsDataset_cortexlab):
    """
    Pools sessions from multiple param_fun.get_snippets output locations
    into one dataset, so TrainExperimentBatchSampler builds one batch per
    (location, session) across the WHOLE corpus instead of just one
    location -- NeuropixelsDataset_cortexlab itself only ever looks at a
    single data_dir whose immediate children are session folders, which is
    the right shape for one location but not for many pooled together.

    Deliberately does not call NeuropixelsDataset_cortexlab.__init__: that
    constructor is built around a single data_dir, so it's simpler to fill
    in the same attributes (mode, unit_order, experiment_unit_map,
    all_files) it expects directly, reusing its inherited
    select_good_units_files/__getitem__/_augment for everything else.
    """

    def __init__(self, location_dirs, mode="train"):
        self.mode = mode
        self.unit_order = "filesystem"
        self.experiment_unit_map = {}
        global_id = 0
        for location_dir in location_dirs:
            data_dir = os.path.join(location_dir, "processed_waveforms")
            if not os.path.isdir(data_dir):
                continue
            sessions = sorted(
                os.listdir(data_dir), key=lambda s: int(s) if s.isdigit() else s
            )
            for session in sessions:
                session_dir = os.path.join(data_dir, session)
                files = self.select_good_units_files(session_dir, load_pre_merge=False)
                if files:
                    self.experiment_unit_map[global_id] = files
                    global_id += 1

        self.all_files = [
            (exp, file)
            for exp, files in self.experiment_unit_map.items()
            for file in files
        ]
        print(
            f"MultiLocationFinetuneDataset: {len(self.all_files)} unit-files across "
            f"{len(self.experiment_unit_map)} session(s) from {len(location_dirs)} location(s)."
        )


def run_train_ae(args):
    locations = list_preprocessed_locations()
    if not locations:
        raise RuntimeError(
            f"No preprocessed locations found under {PREPROCESSED_OUTPUT}. "
            f"Run --stage preprocess first."
        )
    print(f"Pretraining AE on {len(locations)} preprocessed location(s).")
    dataset = ConcatDataset(
        [AE_NeuropixelsDataset(loc, batch_size=args.batchsize) for loc in locations]
    )
    train_ae_mod.run_training(
        exp_name=args.exp_name,
        dataset=dataset,
        lr=args.lr_ae,
        total_epoch=args.epochs_ae,
        cont=args.cont,
        batchsize=args.batchsize,
    )


def run_train_finetune(args):
    locations = list_preprocessed_locations()
    if not locations:
        raise RuntimeError(
            f"No preprocessed locations found under {PREPROCESSED_OUTPUT}. "
            f"Run --stage preprocess first."
        )
    # NB: train_finetune.run_finetune looks for the AE checkpoint under
    # ModelExp/AE_experiments/<exp_name> -- i.e. it reuses exp_name for both
    # stages by convention, not a separate "which AE checkpoint" argument.
    # This wrapper passes the same args.exp_name to both run_train_ae and
    # run_train_finetune so that lines up automatically.
    dataset = MultiLocationFinetuneDataset(locations, mode="train")
    train_finetune_mod.run_finetune(
        exp_name=args.exp_name,
        dataset=dataset,
        lr_backbone=args.lr_backbone,
        lr_enc=args.lr_enc,
        lr_proj=args.lr_proj,
        total_epoch=args.epochs_finetune,
        cont=args.cont,
        batchsize=args.batchsize,
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--stage",
        choices=["preprocess", "train_ae", "train_finetune", "all"],
        default="all",
        help="Which stage to run (default: all, in order)",
    )
    parser.add_argument("--exp_name", default=DEFAULT_EXP_NAME)
    parser.add_argument(
        "--limit", type=int, default=None,
        help="Only preprocess the first N locations (for a quick smoke test)",
    )
    parser.add_argument(
        "--redo", action="store_true",
        help="Reprocess locations even if their .done sentinel is present",
    )
    parser.add_argument(
        "--cont", action="store_true",
        help="Resume training from the latest checkpoint under this exp_name",
    )
    parser.add_argument("--batchsize", type=int, default=32)
    parser.add_argument("--lr_ae", type=float, default=1e-5)
    parser.add_argument("--epochs_ae", type=int, default=300)
    parser.add_argument("--lr_backbone", type=float, default=2e-6)
    parser.add_argument("--lr_enc", type=float, default=2e-5)
    parser.add_argument("--lr_proj", type=float, default=1.1e-4)
    parser.add_argument("--epochs_finetune", type=int, default=50)
    args = parser.parse_args()

    print(f"Raw input (read-only) : {BASE_INPUT}")
    print(f"Preprocessed output   : {PREPROCESSED_OUTPUT}")
    print(f"Experiment name       : {args.exp_name}")

    if args.stage in ("preprocess", "all"):
        preprocess_all(limit=args.limit, redo=args.redo)
    if args.stage in ("train_ae", "all"):
        run_train_ae(args)
    if args.stage in ("train_finetune", "all"):
        run_train_finetune(args)


if __name__ == "__main__":
    main()
