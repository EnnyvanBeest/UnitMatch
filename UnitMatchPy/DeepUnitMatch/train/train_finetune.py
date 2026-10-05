import logging
import os
import argparse
import subprocess
import numpy as np
import tqdm
import torch
import torch.optim as optim
from torch.utils.data import DataLoader
from torch.utils.tensorboard import SummaryWriter
from pathlib import Path

from utils import metric
from utils.losses import *
from utils.npdataset import (
    NeuropixelsDataset_cortexlab,
    TrainExperimentBatchSampler,
    ValidationExperimentBatchSampler,
    split_sessions,
)
import json

# Fraction of recording sessions held out (never trained on) to validate the
# contrastive fine-tuning; fixed seed so every run uses the same split.
VAL_FRACTION = 0.05
SPLIT_SEED = 0


def make_session_loaders(dataset, batchsize, ckpt_folder):
    """
    Train/validation loaders over disjoint sets of whole sessions (see
    npdataset.split_sessions); the split is saved as session_split.json in
    ckpt_folder. Returns (train_loader, val_loader); val_loader is None if
    there are too few sessions to hold any out.
    """
    train_sessions, val_sessions = split_sessions(dataset, VAL_FRACTION, SPLIT_SEED)
    with open(os.path.join(ckpt_folder, "session_split.json"), "w") as f:
        json.dump(
            {
                "val_fraction": VAL_FRACTION,
                "seed": SPLIT_SEED,
                "train_sessions": [str(s) for s in train_sessions],
                "val_sessions": [str(s) for s in val_sessions],
            },
            f,
            indent=1,
        )
    print(f"Sessions: {len(train_sessions)} for training, {len(val_sessions)} held out for validation")
    train_loader = DataLoader(
        dataset,
        batch_sampler=TrainExperimentBatchSampler(dataset, batchsize, shuffle=True, experiments=train_sessions),
    )
    val_loader = None
    if val_sessions:
        val_loader = DataLoader(
            dataset,
            batch_sampler=ValidationExperimentBatchSampler(dataset, shuffle=True, experiments=val_sessions),
        )
    return train_loader, val_loader
from utils.mymodel import *
from testing.test import load_encoder_state, latest_checkpoint

logger = logging.getLogger(__name__)
device = torch.device("cuda" if torch.cuda.is_available() else "cpu")


def _encode_pair(model, batch):
    """
    Encode both halves of a batch. Batches are 5-tuples (waveforms only) or,
    with dataset.with_geometry, 9-tuples that also carry each half's channel
    positions/validity, which are then passed to the encoder.
    """
    if len(batch) == 9:
        estimates, candidates, _, pos_e, valid_e, pos_c, valid_c, _, _ = batch
        geometry = [pos_e, valid_e, pos_c, valid_c]
    else:
        estimates, candidates = batch[0], batch[1]
        geometry = None
    if torch.cuda.is_available():
        estimates, candidates = estimates.cuda(), candidates.cuda()
        if geometry is not None:
            geometry = [g.cuda() for g in geometry]
    if geometry is None:
        return model(estimates), model(candidates)
    pos_e, valid_e, pos_c, valid_c = geometry
    return model(estimates, pos_e, valid_e), model(candidates, pos_c, valid_c)


def validation(epoch, model, projector, val_loader, clip_loss, writer):
    if val_loader is None:
        return
    # held-out sessions are evaluated without augmentation
    dataset = val_loader.dataset
    previous_mode = getattr(dataset, "mode", None)
    dataset.mode = "val"
    try:
        _validation(epoch, model, projector, val_loader, clip_loss, writer)
    finally:
        dataset.mode = previous_mode


def _validation(epoch, model, projector, val_loader, clip_loss, writer):
    model.eval()
    projector.eval()
    clip_loss.eval()
    losses = metric.AverageMeter()
    experiment_accuracies = []
    if torch.cuda.is_available():
        model = model.cuda()
        projector = projector.cuda()
        clip_loss = clip_loss.cuda()

    with torch.no_grad():
        progress_bar = tqdm.tqdm(
            total=len(val_loader), desc="Epoch {:3d}".format(epoch)
        )
        for batch in val_loader:
            bsz = batch[0].shape[0]
            # Forward pass
            enc_estimates, enc_candidates = _encode_pair(model, batch)  # [bsz, n_output]
            proj_estimates = projector(enc_estimates)
            proj_candidates = projector(enc_candidates)
            loss_clip = clip_loss(proj_estimates, proj_candidates)
            loss = loss_clip
            losses.update(loss.item(), bsz)

            probs = clip_prob(enc_estimates, enc_candidates)
            predicted_indices = torch.argmax(
                probs, dim=1
            )  # Get the index of the max probability for each batch element
            ground_truth_indices = torch.arange(
                bsz, device=device
            )  # Diagonal indices as ground truth
            correct_predictions = (
                (predicted_indices == ground_truth_indices).sum().item()
            )  # Count correct predictions
            accuracy = correct_predictions / bsz
            experiment_accuracies.append(accuracy)
            progress_bar.update(1)
        progress_bar.close()

    print(
        "Epoch: %d" % (epoch),
        "Validation Loss: %.9f" % (losses.avg),
        "Validation Accuracy: %.9f" % (np.mean(experiment_accuracies)),
    )
    # print('clip temp tau', clip_loss.temp_tau)
    writer.add_scalar("Validation/Loss", losses.avg, epoch)
    writer.add_scalar("Validation/Accuracy", np.mean(experiment_accuracies), epoch)
    return


def train(epoch, model, projector, optimizer, train_loader, clip_loss, writer):
    model.train()
    projector.train()  # dropout on (validation switches it off)
    clip_loss.train()
    losses = metric.AverageMeter()
    iteration = len(train_loader) * epoch
    if torch.cuda.is_available():
        model = model.cuda()
        projector = projector.cuda()
        clip_loss = clip_loss.cuda()

    progress_bar = tqdm.tqdm(total=len(train_loader), desc="Epoch {:3d}".format(epoch))
    for batch in train_loader:
        bsz = batch[0].shape[0]
        optimizer.zero_grad()
        enc_estimates, enc_candidates = _encode_pair(model, batch)  # [bsz, n_output]
        proj_estimates = projector(enc_estimates)
        proj_candidates = projector(enc_candidates)
        loss_clip = clip_loss(proj_estimates, proj_candidates)
        loss = loss_clip
        losses.update(loss.item(), bsz)
        loss.backward()
        optimizer.step()
        progress_bar.update(1)
        iteration += 1
        if iteration % 50 == 0:
            writer.add_scalar("Train/Loss", losses.avg, iteration)

    progress_bar.close()
    print(" Epoch: %d" % (epoch), "Loss: %.9f" % (losses.avg))
    return


def run(args):
    save_folder = os.path.join("ModelExp", "experiments", args.exp_name)
    ckpt_folder = os.path.join(save_folder, "ckpt")
    log_folder = os.path.join(save_folder, "log")
    os.makedirs(ckpt_folder, exist_ok=True)
    os.makedirs(log_folder, exist_ok=True)
    writer = SummaryWriter(log_dir=log_folder)

    train_data_root = args.train_root
    np_dataset = NeuropixelsDataset_cortexlab(
        data_dir=train_data_root, batch_size=args.batchsize, mode="train"
    )
    train_loader, val_loader = make_session_loaders(np_dataset, args.batchsize, ckpt_folder)

    print("train dataset length: %d" % (len(np_dataset)))

    print(
        f"To open Tensorboard, run this: tensorboard --logdir {os.path.join(os.getcwd(), log_folder)}"
    )

    model = SpatioTemporalCNN_V2(n_channel=30, n_time=60, n_output=256).to(device)

    model = model.double()
    finetune_folder = os.path.join("ModelExp", "AE_experiments", args.finetune)
    ckpt_finetune_folder = os.path.join(finetune_folder, "ckpt")
    ckpt_lst = os.listdir(ckpt_finetune_folder)
    ckpt_lst.sort(key=lambda x: int(x.split("_")[-1]))
    read_path = os.path.join(ckpt_finetune_folder, ckpt_lst[-1])
    print("load checkpoint from %s" % (read_path))
    checkpoint = torch.load(read_path)
    model.load_state_dict(checkpoint["encoder"])
    for name, param in model.named_parameters():
        if "FcBlock" not in name:
            param.requires_grad = False

    projector = Projector(
        input_dim=256, output_dim=128, hidden_dim=128, n_hidden_layers=1, dropout=0.1
    ).to(device)
    projector = projector.double()

    # clip_loss = ClipLoss1D().to(device)
    clip_loss = CustomClipLoss().to(device)

    encoder_fc_params = [
        param
        for name, param in model.named_parameters()
        if "FcBlock" in name and param.requires_grad
    ]
    projector_params = list(
        projector.parameters()
    )  # Assuming projector is defined elsewhere
    clip_loss_params = list(
        clip_loss.parameters()
    )  # Assuming clip_loss is defined elsewhere

    # Combine parameters from different parts with their respective learning rates
    optimizer_params = [
        {
            "params": encoder_fc_params,
            "lr": args.lr_enc,
        },  # Smaller learning rate for FcBlock
        {
            "params": projector_params + clip_loss_params,
            "lr": args.lr_proj,
        },  # Larger learning rate for projector and clip_loss
    ]

    optimizer = optim.Adam(optimizer_params)

    start_epoch = 0
    if args.cont:
        # load latest checkpoint, if one exists yet -- --cont on a folder with
        # no checkpoints saved yet (e.g. a freshly created exp) just starts
        # from epoch 0 instead of crashing on ckpt_lst[-1].
        ckpt_lst = [f for f in os.listdir(ckpt_folder) if f.startswith("ckpt_epoch_")]
        if ckpt_lst:
            ckpt_lst.sort(key=lambda x: int(x.split("_")[-1]))
            read_path = os.path.join(ckpt_folder, ckpt_lst[-1])
            print("load checkpoint from %s" % (read_path))
            checkpoint = torch.load(read_path)
            model.load_state_dict(checkpoint["model"])
            optimizer.load_state_dict(checkpoint["optimizer"])
            clip_loss.load_state_dict(checkpoint["clip_loss"])
            if "projector" in checkpoint:
                projector.load_state_dict(checkpoint["projector"])
            else:
                print("WARNING: checkpoint has no projector state (saved before it was stored); projector restarts")
            start_epoch = checkpoint["epoch"] + 1
        else:
            print(f"--cont given but no checkpoint found in {ckpt_folder}; starting from epoch 0")

    if args.total_epoch == 0:
        # don't train, just save the untrained checkpoint
        state = {
            "model": model.state_dict(),
            "optimizer": optimizer.state_dict(),
            "clip_loss": clip_loss.state_dict(),
            "projector": projector.state_dict(),
            "epoch": 0,
        }
        save_file = os.path.join(ckpt_folder, "ckpt_epoch_0")
        torch.save(state, save_file)

    for epoch in range(start_epoch, args.total_epoch):
        train(epoch, model, projector, optimizer, train_loader, clip_loss, writer)
        if epoch % args.save_freq == 0:
            state = {
                "model": model.state_dict(),
                "optimizer": optimizer.state_dict(),
                "clip_loss": clip_loss.state_dict(),
                "projector": projector.state_dict(),
                "epoch": epoch,
            }
            save_file = os.path.join(ckpt_folder, "ckpt_epoch_%s" % (str(epoch)))
            torch.save(state, save_file)

        # validate and test
        validation(epoch, model, projector, val_loader, clip_loss, writer)
    # test(epoch, model, test_loader, writer)

    return


def run_finetune(
    exp_name,
    dataset,
    lr_enc=2 * 1e-5,
    lr_proj=1.1 * 1e-4,
    save_freq=1,
    total_epoch=50,
    cont=False,
    batchsize=40,
    launch_tensorboard=True,
    n_output=256,
    negative_weight=10.0,
    ae_exp_name=None,
    random_backbone=False,
    lr_backbone=None,
    with_geometry=False,
    amp_noise_jitter=False,
    symmetric_loss=False,
):
    """
    Contrastive fine-tuning of an autoencoder-pretrained encoder.

    launch_tensorboard: Whether to kill any running tensorboard.exe and launch
        a new one for this run (default: True, matching prior behaviour). Set
        False for unattended/batch/parallel callers.
    n_output: encoder output size (must match the AE checkpoint).
    negative_weight: W_ij for different neurons in the loss (W_ii = 1).
    ae_exp_name: AE experiment to start from (ModelExp/AE_experiments/<name>);
        default exp_name. Lets several fine-tunings share one AE.
    random_backbone: don't load AE weights -- fine-tune on a randomly
        initialised, frozen backbone (the "fine-tuned only" baseline).

    Optional v2 features (defaults = v1: frozen backbone, waveform only,
    roll augmentation only, one-directional loss):
    lr_backbone: if given, the conv/spatial backbone is trained too, with this
        (small) learning rate; None keeps it frozen.
    with_geometry: give the encoder each channel's position relative to the
        peak (ChannelPositionalBias); snippets must contain ChannelPos/ChannelValid.
    amp_noise_jitter: add random gain/noise to training snippets.
    symmetric_loss: average the loss over both directions.

    The options are stored in every checkpoint ("config"); inference reads
    "with_geometry" from there to feed channel positions to such models.
    """

    from argparse import Namespace

    args = Namespace(
        exp_name=exp_name,
        lr_enc=lr_enc,
        lr_proj=lr_proj,
        save_freq=save_freq,
        total_epoch=total_epoch,
        cont=cont,
        batchsize=batchsize,
        dataset=dataset,
    )
    config = {
        "n_output": n_output,
        "negative_weight": negative_weight,
        "ae_exp_name": ae_exp_name or exp_name,
        "random_backbone": random_backbone,
        "lr_enc": lr_enc,
        "lr_proj": lr_proj,
        "lr_backbone": lr_backbone,
        "with_geometry": with_geometry,
        "amp_noise_jitter": amp_noise_jitter,
        "symmetric_loss": symmetric_loss,
        "batchsize": batchsize,
        "total_epoch": total_epoch,
    }

    save_path = Path(__file__).parent.parent
    save_folder = os.path.join(save_path, "ModelExp", "experiments", args.exp_name)
    if os.path.exists(save_folder) and not args.cont:
        # Check if the folder already exists and if we are not continuing training
        raise ValueError(
            "Running this experiment will overwrite a previous one. Please choose a different name."
        )

    ckpt_folder = os.path.join(save_folder, "ckpt")
    log_folder = os.path.join(save_folder, "log")
    os.makedirs(ckpt_folder, exist_ok=True)
    os.makedirs(log_folder, exist_ok=True)
    writer = SummaryWriter(log_dir=log_folder)
    with open(os.path.join(ckpt_folder, "config.json"), "w") as f:
        json.dump(config, f, indent=1)

    np_dataset = dataset
    np_dataset.with_geometry = with_geometry
    np_dataset.amp_noise_jitter = amp_noise_jitter
    train_loader, val_loader = make_session_loaders(np_dataset, args.batchsize, ckpt_folder)

    print("train dataset length: %d" % (len(np_dataset)))

    if launch_tensorboard:
        import psutil

        # Kill any existing process on port 6006
        for proc in psutil.process_iter(["pid", "name"]):
            try:
                if proc.name() == "tensorboard.exe":
                    proc.kill()
            except:
                pass

        tensorboard_process = subprocess.Popen(
            ["tensorboard", "--logdir", log_folder, "--port", "6006"],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
        print("To view the tensorboard, click this link: http://localhost:6006")

    model = SpatioTemporalCNN_V2(n_channel=30, n_time=60, n_output=n_output).to(device)
    model = model.double()
    if random_backbone:
        print("random_backbone: fine-tuning on a randomly initialised (frozen) backbone")
    else:
        ae_ckpt_folder = os.path.join(save_path, "ModelExp", "AE_experiments", config["ae_exp_name"], "ckpt")
        read_path = latest_checkpoint(ae_ckpt_folder)
        print("load checkpoint from %s" % (read_path))
        checkpoint = torch.load(read_path, weights_only=False)
        load_encoder_state(model, checkpoint["encoder"], source=read_path)

    backbone_params = []
    for name, param in model.named_parameters():
        if "FcBlock" in name:
            continue
        if lr_backbone is None:
            param.requires_grad = False
        elif with_geometry or not name.startswith("pos_bias."):
            backbone_params.append(param)
        else:
            param.requires_grad = False  # positional bias unused without geometry

    projector = Projector(
        input_dim=n_output, output_dim=128, hidden_dim=128, n_hidden_layers=1, dropout=0.1
    ).to(device)
    projector = projector.double()

    clip_loss = CustomClipLoss(negative_weight=negative_weight, symmetric=symmetric_loss).to(device)

    encoder_fc_params = [
        param
        for name, param in model.named_parameters()
        if "FcBlock" in name and param.requires_grad
    ]
    projector_params = list(projector.parameters())
    clip_loss_params = list(clip_loss.parameters())

    # Combine parameters from different parts with their respective learning rates
    optimizer_params = [
        {
            "params": encoder_fc_params,
            "lr": args.lr_enc,
        },  # Smaller learning rate for FcBlock
        {
            "params": projector_params + clip_loss_params,
            "lr": args.lr_proj,
        },  # Larger learning rate for projector and clip_loss
    ]
    if backbone_params:
        optimizer_params.append({"params": backbone_params, "lr": lr_backbone})

    optimizer = optim.Adam(optimizer_params)

    def state(epoch):
        return {
            "model": model.state_dict(),
            "optimizer": optimizer.state_dict(),
            "clip_loss": clip_loss.state_dict(),
            "projector": projector.state_dict(),
            "config": config,
            "epoch": epoch,
        }

    start_epoch = 0
    if args.cont:
        # load latest checkpoint, if one exists yet -- --cont on a folder with
        # no checkpoints saved yet (e.g. a freshly created exp) just starts
        # from epoch 0 instead of crashing on ckpt_lst[-1].
        ckpt_lst = [f for f in os.listdir(ckpt_folder) if f.startswith("ckpt_epoch_")]
        if ckpt_lst:
            read_path = latest_checkpoint(ckpt_folder)
            print("load checkpoint from %s" % (read_path))
            checkpoint = torch.load(read_path, weights_only=False)
            model.load_state_dict(checkpoint["model"])
            optimizer.load_state_dict(checkpoint["optimizer"])
            clip_loss.load_state_dict(checkpoint["clip_loss"])
            if "projector" in checkpoint:
                projector.load_state_dict(checkpoint["projector"])
            else:
                print("WARNING: checkpoint has no projector state (saved before it was stored); projector restarts")
            start_epoch = checkpoint["epoch"] + 1
        else:
            print(f"--cont given but no checkpoint found in {ckpt_folder}; starting from epoch 0")

    if args.total_epoch == 0:
        # don't train, just save the untrained checkpoint
        torch.save(state(0), os.path.join(ckpt_folder, "ckpt_epoch_0"))

    for epoch in range(start_epoch, args.total_epoch):
        train(epoch, model, projector, optimizer, train_loader, clip_loss, writer)
        if epoch % args.save_freq == 0 or epoch == args.total_epoch - 1:  # always keep the final epoch
            torch.save(state(epoch), os.path.join(ckpt_folder, "ckpt_epoch_%s" % (str(epoch))))

        # validate
        validation(epoch, model, projector, val_loader, clip_loss, writer)

    return


if __name__ == "__main__":
    arg_parser = argparse.ArgumentParser()
    arg_parser.add_argument(
        "--exp_name",
        "-e",
        type=str,
        required=True,
        help="The checkpoints and logs will be save in ./checkpoint/$EXP_NAME",
    )
    arg_parser.add_argument(
        "--finetune",
        "-f",
        type=str,
        required=True,
        help="Load the AE encoder from the path ./checkpoint/$finetune",
    )
    arg_parser.add_argument(
        "--lr_enc",
        "-le",
        type=float,
        default=2 * 1e-5,
        help="Learning rate for encoder",
    )
    arg_parser.add_argument(
        "--lr_proj",
        "-lp",
        type=float,
        default=1.1 * 1e-4,
        help="Learning rate for projector",
    )
    arg_parser.add_argument(
        "--save_freq", "-s", type=int, default=1, help="frequency of saving model"
    )
    arg_parser.add_argument(
        "--total_epoch",
        "-t",
        type=int,
        default=50,
        help="total epoch number for training",
    )
    arg_parser.add_argument(
        "--cont",
        "-c",
        action="store_true",
        help="whether to load saved checkpoints from $EXP_NAME and continue training",
    )
    arg_parser.add_argument(
        "--batchsize", "-b", type=int, default=40, help="batch size"
    )
    arg_parser.add_argument(
        "--train_root",
        type=str,
        default=r"\path\to\your\data",
        help="root directory of training data",
    )
    args = arg_parser.parse_args()

    run(args)
