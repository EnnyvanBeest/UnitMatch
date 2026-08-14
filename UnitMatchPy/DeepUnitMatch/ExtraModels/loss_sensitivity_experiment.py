"""Contrastive-loss off-diagonal weight sensitivity experiments.

This script is intentionally self-contained in Rebuttal/ so the package training
code remains unchanged. It mirrors DeepUnitMatch/train/train_finetune.py while
exposing CustomClipLoss.negative_weight as an experiment sweep parameter.
"""

from __future__ import annotations

import argparse
import csv
import h5py
import json
import random
import re
import sys
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence

import numpy as np
import torch
import torch.optim as optim
from torch.utils.data import DataLoader, Dataset
import tqdm


REBUTTAL_DIR = Path(__file__).resolve().parent
REPO_ROOT = REBUTTAL_DIR.parent
DEEPUNITMATCH_DIR = REPO_ROOT / "DeepUnitMatch"
DEFAULT_DATA_ROOT = Path(r"F:\UCL\kenneth_harris\code\Data_DUM")

for import_path in (DEEPUNITMATCH_DIR, REPO_ROOT):
    import_path_str = str(import_path)
    if import_path_str not in sys.path:
        sys.path.insert(0, import_path_str)

try:
    from torch.utils.tensorboard import SummaryWriter
except Exception:  # pragma: no cover - only used when tensorboard is absent.
    class SummaryWriter:  # type: ignore[no-redef]
        def __init__(self, *args, **kwargs):
            pass

        def add_scalar(self, *args, **kwargs):
            pass

        def close(self):
            pass


from DeepUnitMatch.utils.losses import CustomClipLoss, Projector, clip_prob
from DeepUnitMatch.utils.metric import AverageMeter
from DeepUnitMatch.utils.mymodel import SpatioTemporalCNN_V2
from DeepUnitMatch.utils.npdataset import (
    NeuropixelsDataset,
    NeuropixelsDataset_cortexlab,
    TrainExperimentBatchSampler,
    ValidationExperimentBatchSampler,
)
from DeepUnitMatch.utils.train_utils import read_good_ids


def set_seed(seed: int, deterministic: bool = False) -> None:
    """Seed Python, NumPy, and PyTorch for comparable sweeps."""

    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    if torch.cuda.is_available():
        torch.cuda.manual_seed_all(seed)

    if deterministic:
        torch.backends.cudnn.deterministic = True
        torch.backends.cudnn.benchmark = False
        try:
            torch.use_deterministic_algorithms(True, warn_only=True)
        except TypeError:
            torch.use_deterministic_algorithms(True)


def resolve_device(device_name: str) -> torch.device:
    if device_name == "auto":
        return torch.device("cuda" if torch.cuda.is_available() else "cpu")
    return torch.device(device_name)


def _checkpoint_epoch(path: Path) -> int:
    match = re.search(r"(\d+)(?:\.[^.]+)?$", path.name)
    if match is None:
        return -1
    return int(match.group(1))


def latest_checkpoint(ckpt_dir: Path) -> Path:
    ckpt_dir = Path(ckpt_dir)
    if not ckpt_dir.exists():
        raise FileNotFoundError(f"Checkpoint directory does not exist: {ckpt_dir}")

    ckpts = [path for path in ckpt_dir.iterdir() if path.is_file()]
    if not ckpts:
        raise FileNotFoundError(f"No checkpoint files found in: {ckpt_dir}")

    return max(ckpts, key=lambda path: (_checkpoint_epoch(path), path.stat().st_mtime))


def resolve_ae_checkpoint(ae_checkpoint: Optional[str], ae_exp_name: Optional[str]) -> Path:
    """Resolve either an explicit AE checkpoint path or an AE experiment name."""

    if ae_checkpoint:
        path = Path(ae_checkpoint).expanduser()
        if not path.is_absolute():
            path = (Path.cwd() / path).resolve()
        if not path.is_file():
            raise FileNotFoundError(f"AE checkpoint not found: {path}")
        return path

    if not ae_exp_name:
        raise ValueError("Provide either --ae_checkpoint or --ae_exp_name.")

    as_path = Path(ae_exp_name).expanduser()
    candidate_dirs: List[Path] = []
    if as_path.exists():
        resolved = as_path.resolve()
        candidate_dirs.append(resolved / "ckpt" if (resolved / "ckpt").exists() else resolved)

    candidate_dirs.extend(
        [
            DEEPUNITMATCH_DIR / "ModelExp" / "AE_experiments" / ae_exp_name / "ckpt",
            REPO_ROOT / "ModelExp" / "AE_experiments" / ae_exp_name / "ckpt",
            Path.cwd() / "ModelExp" / "AE_experiments" / ae_exp_name / "ckpt",
        ]
    )

    for ckpt_dir in candidate_dirs:
        if ckpt_dir.exists():
            return latest_checkpoint(ckpt_dir)

    searched = "\n".join(str(path) for path in candidate_dirs)
    raise FileNotFoundError(
        f"Could not resolve AE experiment '{ae_exp_name}'. Searched:\n{searched}"
    )


def format_setting_value(value: float) -> str:
    return ("%g" % value).replace("-", "m").replace(".", "p")


def _looks_like_cortexlab_root(root: Path) -> bool:
    """Return True for root/mouse/probe/location/experiment folders."""

    if not root.is_dir():
        return False

    for mouse_dir in root.iterdir():
        if not mouse_dir.is_dir():
            continue
        for probe_dir in mouse_dir.iterdir():
            if not probe_dir.is_dir():
                continue
            for location_dir in probe_dir.iterdir():
                if not location_dir.is_dir():
                    continue
                for experiment_dir in location_dir.iterdir():
                    if (
                        experiment_dir.is_dir()
                        and (experiment_dir / "metadata.json").is_file()
                        and (experiment_dir / "processed_waveforms").is_dir()
                    ):
                        return True
    return False


def resolve_dataset_root(train_root: str, dataset_kind: str) -> tuple[str, str]:
    """Resolve Data_DUM and public processed layouts to the expected reader root.

    Data_DUM has an outer folder containing ALL_DATA_nomatchtables plus
    matchtables.db. The training readers need the ALL_DATA_nomatchtables child.
    """

    root = Path(train_root).expanduser()
    if not root.is_absolute():
        root = (Path.cwd() / root).resolve()

    if dataset_kind not in {"auto", "cortexlab", "processed", "flat_cortexlab"}:
        raise ValueError(f"Unknown dataset kind: {dataset_kind}")

    all_data_root = root / "ALL_DATA_nomatchtables"
    if dataset_kind in {"auto", "cortexlab"} and all_data_root.is_dir():
        return str(all_data_root), "cortexlab"

    if dataset_kind == "auto":
        if (root / "processed_waveforms").is_dir():
            return str(root), "processed"
        if _looks_like_cortexlab_root(root):
            return str(root), "cortexlab"
        raise ValueError(
            "Could not infer dataset kind. Expected either a Data_DUM folder "
            "containing ALL_DATA_nomatchtables, a cortexlab root containing "
            "mouse/probe/location/experiment folders, or a public processed "
            "root containing processed_waveforms/."
        )

    if dataset_kind == "cortexlab" and not _looks_like_cortexlab_root(root):
        raise ValueError(
            f"{root} does not look like a cortexlab root. For Data_DUM, pass "
            "the outer Data_DUM folder or its ALL_DATA_nomatchtables child."
        )

    if dataset_kind == "flat_cortexlab":
        return str(root), dataset_kind

    if dataset_kind == "processed" and not (root / "processed_waveforms").is_dir():
        raise ValueError(
            f"{root} does not look like a processed root because it has no "
            "processed_waveforms/ child."
        )

    return str(root), dataset_kind


def write_json(path: Path, payload: Dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")


def write_rows_csv(path: Path, rows: Sequence[Dict], fieldnames: Optional[Sequence[str]] = None) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        return

    if fieldnames is None:
        fields: List[str] = []
        for row in rows:
            for key in row.keys():
                if key not in fields:
                    fields.append(key)
        fieldnames = fields

    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(fieldnames), extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


class NestedCortexlabPairDataset(Dataset):
    """Pair dataset for Data_DUM's root/mouse/probe/location/experiment layout."""

    def __init__(self, root: str, batch_size: int = 40, mode: str = "val"):
        self.root = Path(root).resolve()
        self.batch_size = batch_size
        self.mode = mode
        self.mouse_names = sorted(
            path.name for path in self.root.iterdir() if path.is_dir()
        )
        experiment_unit_map = read_good_ids(
            str(self.root),
            self.batch_size,
            self.mouse_names,
            finetune=True,
        )

        self.experiment_paths = sorted(experiment_unit_map.keys())
        self.experiment_unit_map = {
            index: experiment_unit_map[experiment_path]
            for index, experiment_path in enumerate(self.experiment_paths)
        }
        self.experiment_path_by_index = {
            index: experiment_path
            for index, experiment_path in enumerate(self.experiment_paths)
        }
        self.all_files = [
            (experiment_index, file_path)
            for experiment_index, files in self.experiment_unit_map.items()
            for file_path in files
        ]

        if len(self.all_files) < 1:
            print("No data in contrastive dataset. Try a smaller batch size?")
        else:
            print(
                f"Initialised nested contrastive dataset with "
                f"{len(self.all_files)} files across "
                f"{len(self.experiment_unit_map)} experiments."
            )

    def __len__(self) -> int:
        return len(self.all_files)

    def __getitem__(self, index: int):
        experiment_index, neuron_file = self.all_files[index]
        with h5py.File(neuron_file, "r") as handle:
            waveform = handle["waveform"][()]
            max_site_pos = handle["MaxSitepos"][()]

        if waveform.shape != (60, 30, 2):
            waveform = np.zeros((60, 30, 2))

        if self.mode == "train":
            waveform_fh = self._augment_original(waveform[..., 0].copy())
            waveform_sh = self._augment_original(waveform[..., 1].copy())
        else:
            waveform_fh = waveform[..., 0]
            waveform_sh = waveform[..., 1]

        return (
            waveform_fh,
            waveform_sh,
            max_site_pos,
            self.experiment_path_by_index[experiment_index],
            neuron_file,
        )

    def _augment_original(self, data: np.ndarray) -> np.ndarray:
        roll_choice = random.choice(["roll_up", "roll_down", "none"])

        if roll_choice == "roll_up":
            channels = data.shape[1]
            odd_indices = np.arange(0, channels - 1, 2)
            even_indices = np.arange(1, channels - 1, 2)
            if len(odd_indices) > 1:
                data[:, odd_indices[:-1]] = data[:, odd_indices[1:]]
            if len(even_indices) > 1:
                data[:, even_indices[:-1]] = data[:, even_indices[1:]]
        elif roll_choice == "roll_down":
            channels = data.shape[1]
            odd_indices = np.arange(2, channels, 2)
            even_indices = np.arange(3, channels, 2)
            if len(odd_indices) > 0:
                data[:, odd_indices] = data[:, odd_indices - 2]
            if len(even_indices) > 0:
                data[:, even_indices] = data[:, even_indices - 2]

        return data


def make_paired_dataset(
    train_root: str,
    dataset_kind: str,
    batchsize: int,
    mode: str,
):
    if dataset_kind == "processed":
        return NeuropixelsDataset(train_root, batch_size=batchsize, mode=mode)
    if dataset_kind == "cortexlab":
        return NestedCortexlabPairDataset(train_root, batch_size=batchsize, mode=mode)
    if dataset_kind == "flat_cortexlab":
        return NeuropixelsDataset_cortexlab(train_root, batch_size=batchsize, mode=mode)
    raise ValueError(f"Unknown dataset kind: {dataset_kind}")


def make_paired_loaders(
    train_root: str,
    dataset_kind: str,
    batchsize: int,
    num_workers: int,
    validation_mode: str,
):
    train_root, dataset_kind = resolve_dataset_root(train_root, dataset_kind)
    train_dataset = make_paired_dataset(
        train_root=train_root,
        dataset_kind=dataset_kind,
        batchsize=batchsize,
        mode="train",
    )
    val_dataset = make_paired_dataset(
        train_root=train_root,
        dataset_kind=dataset_kind,
        batchsize=batchsize,
        mode=validation_mode,
    )

    train_sampler = TrainExperimentBatchSampler(train_dataset, batchsize, shuffle=True)
    val_sampler = ValidationExperimentBatchSampler(val_dataset, shuffle=False)

    train_loader = DataLoader(
        train_dataset,
        batch_sampler=train_sampler,
        num_workers=num_workers,
        pin_memory=torch.cuda.is_available(),
    )
    val_loader = DataLoader(
        val_dataset,
        batch_sampler=val_sampler,
        num_workers=num_workers,
        pin_memory=torch.cuda.is_available(),
    )
    return train_loader, val_loader


def _move_waveforms(tensor: torch.Tensor, device: torch.device) -> torch.Tensor:
    return tensor.to(device=device, dtype=torch.float64, non_blocking=True)


def _accuracy_from_probabilities(probs: torch.Tensor) -> float:
    if probs.numel() == 0:
        return float("nan")
    target = torch.arange(probs.size(0), device=probs.device)
    predicted = torch.argmax(probs, dim=1)
    return (predicted == target).double().mean().item()


def _load_encoder_state(
    model: SpatioTemporalCNN_V2,
    ae_checkpoint: Path,
    device: torch.device,
) -> None:
    checkpoint = torch.load(ae_checkpoint, map_location=device)
    if "encoder" in checkpoint:
        encoder_state = checkpoint["encoder"]
    elif "model" in checkpoint:
        encoder_state = checkpoint["model"]
    else:
        encoder_state = checkpoint

    load_error: Optional[RuntimeError] = None
    try:
        model.load_state_dict(encoder_state)
        return
    except RuntimeError as exc:
        load_error = exc

    stripped = {
        key.replace("encoder.", "", 1): value
        for key, value in encoder_state.items()
        if key.startswith("encoder.")
    }
    if not stripped:
        raise load_error if load_error is not None else RuntimeError(
            f"Could not load encoder state from {ae_checkpoint}"
        )
    model.load_state_dict(stripped)


def _freeze_encoder(model: SpatioTemporalCNN_V2, freeze_non_fc: bool) -> None:
    if not freeze_non_fc:
        return
    for name, parameter in model.named_parameters():
        parameter.requires_grad = "FcBlock" in name


def _make_optimizer(
    model: SpatioTemporalCNN_V2,
    projector: Projector,
    clip_loss: CustomClipLoss,
    lr_enc: float,
    lr_proj: float,
) -> optim.Optimizer:
    encoder_params = [
        parameter
        for parameter in model.parameters()
        if parameter.requires_grad
    ]
    head_params = list(projector.parameters()) + list(clip_loss.parameters())

    optimizer_params = []
    if encoder_params:
        optimizer_params.append({"params": encoder_params, "lr": lr_enc})
    optimizer_params.append({"params": head_params, "lr": lr_proj})
    return optim.Adam(optimizer_params)


def train_contrastive_epoch(
    epoch: int,
    model: SpatioTemporalCNN_V2,
    projector: Projector,
    optimizer: optim.Optimizer,
    train_loader: DataLoader,
    clip_loss: CustomClipLoss,
    writer: SummaryWriter,
    device: torch.device,
    log_every: int,
    progress: bool,
) -> float:
    model.train()
    projector.train()
    clip_loss.train()
    losses = AverageMeter()
    iterator = tqdm.tqdm(
        train_loader,
        desc=f"Epoch {epoch:3d} train",
        disable=not progress,
    )

    for step, batch in enumerate(iterator, start=1):
        estimates, candidates = batch[0], batch[1]
        estimates = _move_waveforms(estimates, device)
        candidates = _move_waveforms(candidates, device)
        bsz = estimates.shape[0]

        optimizer.zero_grad(set_to_none=True)
        enc_estimates = model(estimates)
        enc_candidates = model(candidates)
        proj_estimates = projector(enc_estimates)
        proj_candidates = projector(enc_candidates)
        loss = clip_loss(proj_estimates, proj_candidates)
        losses.update(loss.item(), bsz)
        loss.backward()
        optimizer.step()

        iteration = epoch * len(train_loader) + step
        if log_every > 0 and iteration % log_every == 0:
            writer.add_scalar("Train/Loss", losses.avg, iteration)
        if progress:
            iterator.set_postfix(loss=f"{losses.avg:.6f}")

    return float(losses.avg)


def validate_contrastive(
    epoch: int,
    model: SpatioTemporalCNN_V2,
    projector: Projector,
    val_loader: DataLoader,
    clip_loss: CustomClipLoss,
    writer: SummaryWriter,
    device: torch.device,
    progress: bool,
) -> Dict[str, float]:
    model.eval()
    projector.eval()
    clip_loss.eval()

    losses = AverageMeter()
    encoder_acc = AverageMeter()
    projector_acc = AverageMeter()
    iterator = tqdm.tqdm(
        val_loader,
        desc=f"Epoch {epoch:3d} val",
        disable=not progress,
    )

    with torch.no_grad():
        for batch in iterator:
            estimates, candidates = batch[0], batch[1]
            estimates = _move_waveforms(estimates, device)
            candidates = _move_waveforms(candidates, device)
            bsz = estimates.shape[0]

            enc_estimates = model(estimates)
            enc_candidates = model(candidates)
            proj_estimates = projector(enc_estimates)
            proj_candidates = projector(enc_candidates)
            loss = clip_loss(proj_estimates, proj_candidates)
            losses.update(loss.item(), bsz)

            enc_probs = clip_prob(enc_estimates, enc_candidates)
            proj_probs = clip_loss.get_probabilities(proj_estimates, proj_candidates)
            encoder_acc.update(_accuracy_from_probabilities(enc_probs), bsz)
            projector_acc.update(_accuracy_from_probabilities(proj_probs), bsz)

            if progress:
                iterator.set_postfix(
                    loss=f"{losses.avg:.6f}",
                    enc_acc=f"{encoder_acc.avg:.4f}",
                    proj_acc=f"{projector_acc.avg:.4f}",
                )

    writer.add_scalar("Validation/Loss", losses.avg, epoch)
    writer.add_scalar("Validation/EncoderAccuracy", encoder_acc.avg, epoch)
    writer.add_scalar("Validation/ProjectorAccuracy", projector_acc.avg, epoch)

    return {
        "val_loss": float(losses.avg),
        "val_encoder_accuracy": float(encoder_acc.avg),
        "val_projector_accuracy": float(projector_acc.avg),
    }


def save_contrastive_checkpoint(
    ckpt_dir: Path,
    epoch: int,
    model: SpatioTemporalCNN_V2,
    projector: Projector,
    optimizer: optim.Optimizer,
    clip_loss: CustomClipLoss,
    n_output: int,
    negative_weight: float,
) -> Path:
    ckpt_dir.mkdir(parents=True, exist_ok=True)
    checkpoint_path = ckpt_dir / f"ckpt_epoch_{epoch}"
    torch.save(
        {
            "model": model.state_dict(),
            "projector": projector.state_dict(),
            "optimizer": optimizer.state_dict(),
            "clip_loss": clip_loss.state_dict(),
            "epoch": epoch,
            "n_output": n_output,
            "negative_weight": negative_weight,
        },
        checkpoint_path,
    )
    return checkpoint_path


def _load_contrastive_checkpoint(
    checkpoint_path: Path,
    model: SpatioTemporalCNN_V2,
    projector: Projector,
    optimizer: optim.Optimizer,
    clip_loss: CustomClipLoss,
    device: torch.device,
) -> int:
    checkpoint = torch.load(checkpoint_path, map_location=device)
    model.load_state_dict(checkpoint["model"])
    if "projector" in checkpoint:
        projector.load_state_dict(checkpoint["projector"])
    optimizer.load_state_dict(checkpoint["optimizer"])
    clip_loss.load_state_dict(checkpoint["clip_loss"])
    return int(checkpoint["epoch"]) + 1


def _summarize_rows(rows: Sequence[Dict], setting: Dict) -> Dict:
    summary = dict(setting)
    summary["epochs_recorded"] = len(rows)
    if not rows:
        return summary

    best_encoder = max(rows, key=lambda row: float(row["val_encoder_accuracy"]))
    best_projector = max(rows, key=lambda row: float(row["val_projector_accuracy"]))
    final = rows[-1]
    summary.update(
        {
            "best_encoder_epoch": best_encoder["epoch"],
            "best_val_encoder_accuracy": best_encoder["val_encoder_accuracy"],
            "best_projector_epoch": best_projector["epoch"],
            "best_val_projector_accuracy": best_projector["val_projector_accuracy"],
            "final_epoch": final["epoch"],
            "final_train_loss": final.get("train_loss", ""),
            "final_val_loss": final["val_loss"],
            "final_val_encoder_accuracy": final["val_encoder_accuracy"],
            "final_val_projector_accuracy": final["val_projector_accuracy"],
            "last_checkpoint": final.get("checkpoint", ""),
        }
    )
    return summary


def train_contrastive_setting(
    *,
    exp_dir: Path,
    train_root: str,
    dataset_kind: str,
    ae_checkpoint: Path,
    n_output: int,
    negative_weight: float,
    total_epoch: int,
    batchsize: int,
    lr_enc: float,
    lr_proj: float,
    save_freq: int,
    seed: int,
    device: torch.device,
    num_workers: int,
    freeze_non_fc: bool,
    projector_output_dim: int,
    projector_hidden_dim: int,
    projector_dropout: float,
    deterministic: bool = False,
    resume: bool = False,
    validation_mode: str = "val",
    log_every: int = 50,
    progress: bool = True,
) -> Dict:
    """Train one contrastive setting and return a summary row."""

    set_seed(seed, deterministic=deterministic)
    train_root, dataset_kind = resolve_dataset_root(train_root, dataset_kind)
    exp_dir.mkdir(parents=True, exist_ok=True)
    ckpt_dir = exp_dir / "ckpt"
    log_dir = exp_dir / "log"
    metrics_path = exp_dir / "metrics.csv"

    write_json(
        exp_dir / "config.json",
        {
            "train_root": train_root,
            "dataset_kind": dataset_kind,
            "ae_checkpoint": str(ae_checkpoint),
            "n_output": n_output,
            "negative_weight": negative_weight,
            "total_epoch": total_epoch,
            "batchsize": batchsize,
            "lr_enc": lr_enc,
            "lr_proj": lr_proj,
            "save_freq": save_freq,
            "seed": seed,
            "device": str(device),
            "num_workers": num_workers,
            "freeze_non_fc": freeze_non_fc,
            "projector_output_dim": projector_output_dim,
            "projector_hidden_dim": projector_hidden_dim,
            "projector_dropout": projector_dropout,
            "validation_mode": validation_mode,
        },
    )

    train_loader, val_loader = make_paired_loaders(
        train_root=train_root,
        dataset_kind=dataset_kind,
        batchsize=batchsize,
        num_workers=num_workers,
        validation_mode=validation_mode,
    )

    model = SpatioTemporalCNN_V2(n_channel=30, n_time=60, n_output=n_output).to(device)
    model = model.double()
    _load_encoder_state(model, ae_checkpoint, device=device)
    _freeze_encoder(model, freeze_non_fc=freeze_non_fc)

    projector = Projector(
        input_dim=n_output,
        output_dim=projector_output_dim,
        hidden_dim=projector_hidden_dim,
        n_hidden_layers=1,
        dropout=projector_dropout,
    ).to(device)
    projector = projector.double()

    clip_loss = CustomClipLoss(negative_weight=negative_weight).to(device)
    optimizer = _make_optimizer(
        model=model,
        projector=projector,
        clip_loss=clip_loss,
        lr_enc=lr_enc,
        lr_proj=lr_proj,
    )

    start_epoch = 0
    if resume and ckpt_dir.exists():
        existing_ckpts = [path for path in ckpt_dir.iterdir() if path.is_file()]
        if existing_ckpts:
            start_epoch = _load_contrastive_checkpoint(
                latest_checkpoint(ckpt_dir),
                model=model,
                projector=projector,
                optimizer=optimizer,
                clip_loss=clip_loss,
                device=device,
            )

    writer = SummaryWriter(log_dir=str(log_dir))
    rows: List[Dict] = []

    if total_epoch == 0 and start_epoch == 0:
        val_metrics = validate_contrastive(
            0,
            model=model,
            projector=projector,
            val_loader=val_loader,
            clip_loss=clip_loss,
            writer=writer,
            device=device,
            progress=progress,
        )
        checkpoint_path = save_contrastive_checkpoint(
            ckpt_dir=ckpt_dir,
            epoch=0,
            model=model,
            projector=projector,
            optimizer=optimizer,
            clip_loss=clip_loss,
            n_output=n_output,
            negative_weight=negative_weight,
        )
        rows.append(
            {
                "epoch": 0,
                "negative_weight": negative_weight,
                "n_output": n_output,
                "train_loss": "",
                **val_metrics,
                "checkpoint": str(checkpoint_path),
            }
        )

    for epoch in range(start_epoch, total_epoch):
        train_loss = train_contrastive_epoch(
            epoch,
            model=model,
            projector=projector,
            optimizer=optimizer,
            train_loader=train_loader,
            clip_loss=clip_loss,
            writer=writer,
            device=device,
            log_every=log_every,
            progress=progress,
        )
        val_metrics = validate_contrastive(
            epoch,
            model=model,
            projector=projector,
            val_loader=val_loader,
            clip_loss=clip_loss,
            writer=writer,
            device=device,
            progress=progress,
        )

        checkpoint_path = ""
        if save_freq > 0 and (epoch % save_freq == 0 or epoch == total_epoch - 1):
            checkpoint_path = str(
                save_contrastive_checkpoint(
                    ckpt_dir=ckpt_dir,
                    epoch=epoch,
                    model=model,
                    projector=projector,
                    optimizer=optimizer,
                    clip_loss=clip_loss,
                    n_output=n_output,
                    negative_weight=negative_weight,
                )
            )

        rows.append(
            {
                "epoch": epoch,
                "negative_weight": negative_weight,
                "n_output": n_output,
                "train_loss": train_loss,
                **val_metrics,
                "checkpoint": checkpoint_path,
            }
        )
        write_rows_csv(metrics_path, rows)

    writer.close()
    write_rows_csv(metrics_path, rows)
    return _summarize_rows(
        rows,
        {
            "negative_weight": negative_weight,
            "n_output": n_output,
            "experiment_dir": str(exp_dir),
        },
    )


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Sweep CustomClipLoss off-diagonal negative weights."
    )
    parser.add_argument("--exp_name", type=str, default="loss_sensitivity")
    parser.add_argument(
        "--train_root",
        type=str,
        default=str(DEFAULT_DATA_ROOT),
        help=(
            "Dataset root. Defaults to the local Data_DUM folder and resolves "
            "Data_DUM/ALL_DATA_nomatchtables automatically."
        ),
    )
    parser.add_argument(
        "--dataset_kind",
        type=str,
        default="cortexlab",
        choices=("auto", "cortexlab", "processed", "flat_cortexlab"),
    )
    parser.add_argument("--ae_checkpoint", type=str, default=None)
    parser.add_argument("--ae_exp_name", type=str, default=None)
    parser.add_argument("--n_output", type=int, default=256)
    parser.add_argument(
        "--negative_weights",
        type=float,
        nargs="+",
        default=[1.0, 5.0, 10.0, 15.0, 20.0],
        # default=[10.0],
    )
    
    
    parser.add_argument("--total_epoch", type=int, default=50)
    parser.add_argument("--batchsize", type=int, default=40)
    parser.add_argument("--lr_enc", type=float, default=2e-5)
    parser.add_argument("--lr_proj", type=float, default=1.1e-4)
    parser.add_argument("--save_freq", type=int, default=1)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--device", type=str, default="auto")
    parser.add_argument("--num_workers", type=int, default=0)
    parser.add_argument("--projector_output_dim", type=int, default=128)
    parser.add_argument("--projector_hidden_dim", type=int, default=128)
    parser.add_argument("--projector_dropout", type=float, default=0.1)
    parser.add_argument("--validation_mode", type=str, default="val", choices=("train", "val"))
    parser.add_argument("--log_every", type=int, default=50)
    parser.add_argument("--deterministic", action="store_true")
    parser.add_argument("--train_all_encoder", action="store_true")
    parser.add_argument("--cont", action="store_true")
    parser.add_argument("--overwrite", action="store_true")
    parser.add_argument("--no_progress", action="store_true")
    parser.add_argument(
        "--output_dir",
        type=str,
        default=str(REBUTTAL_DIR / "results" / "loss_sensitivity"),
    )
    return parser


def main(argv: Optional[Sequence[str]] = None) -> None:
    parser = build_arg_parser()
    args = parser.parse_args(argv)

    device = resolve_device(args.device)
    ae_checkpoint = resolve_ae_checkpoint(args.ae_checkpoint, args.ae_exp_name)
    train_root, dataset_kind = resolve_dataset_root(args.train_root, args.dataset_kind)
    exp_root = Path(args.output_dir).expanduser()
    if not exp_root.is_absolute():
        exp_root = (Path.cwd() / exp_root).resolve()
    exp_root = exp_root / args.exp_name

    if exp_root.exists() and any(exp_root.iterdir()) and not args.cont and not args.overwrite:
        raise ValueError(
            f"Output directory already exists: {exp_root}. "
            "Use --cont to resume or --overwrite to write new metrics into it."
        )

    exp_root.mkdir(parents=True, exist_ok=True)
    write_json(
        exp_root / "sweep_config.json",
        {
            **vars(args),
            "device_resolved": str(device),
            "ae_checkpoint_resolved": str(ae_checkpoint),
            "train_root_resolved": train_root,
            "dataset_kind_resolved": dataset_kind,
        },
    )

    summaries: List[Dict] = []
    for negative_weight in args.negative_weights:
        label = format_setting_value(negative_weight)
        setting_dir = exp_root / f"negative_weight_{label}"
        print(f"\nRunning negative_weight={negative_weight:g} in {setting_dir}")
        summary = train_contrastive_setting(
            exp_dir=setting_dir,
            train_root=train_root,
            dataset_kind=dataset_kind,
            ae_checkpoint=ae_checkpoint,
            n_output=args.n_output,
            negative_weight=negative_weight,
            total_epoch=args.total_epoch,
            batchsize=args.batchsize,
            lr_enc=args.lr_enc,
            lr_proj=args.lr_proj,
            save_freq=args.save_freq,
            seed=args.seed,
            device=device,
            num_workers=args.num_workers,
            freeze_non_fc=not args.train_all_encoder,
            projector_output_dim=args.projector_output_dim,
            projector_hidden_dim=args.projector_hidden_dim,
            projector_dropout=args.projector_dropout,
            deterministic=args.deterministic,
            resume=args.cont,
            validation_mode=args.validation_mode,
            log_every=args.log_every,
            progress=not args.no_progress,
        )
        summaries.append(summary)
        write_rows_csv(exp_root / "summary.csv", summaries)

    print(f"\nSaved summary to {exp_root / 'summary.csv'}")


if __name__ == "__main__":
    main()
