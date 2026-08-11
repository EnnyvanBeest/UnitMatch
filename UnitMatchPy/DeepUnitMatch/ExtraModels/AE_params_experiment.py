"""Autoencoder latent feature dimension experiments.

This script sweeps SpatioTemporalAutoEncoder_V2.n_output while keeping changes
inside Rebuttal/. Each setting can run AE pretraining and optional contrastive
finetuning so the rebuttal can report both reconstruction and matching effects.
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Dict, List, Optional, Sequence

import numpy as np
import torch
import torch.optim as optim
from torch.utils.data import DataLoader, random_split
import tqdm

from loss_sensitivity_experiment import (
    DEFAULT_DATA_ROOT,
    REBUTTAL_DIR,
    SummaryWriter,
    latest_checkpoint,
    resolve_dataset_root,
    resolve_device,
    set_seed,
    train_contrastive_setting,
    write_json,
    write_rows_csv,
)

from DeepUnitMatch.utils.AE_npdataset import AE_NeuropixelsDataset, AE_NeuropixelsDataset_cortexlab
from DeepUnitMatch.utils.losses import AELoss
from DeepUnitMatch.utils.metric import AverageMeter
from DeepUnitMatch.utils.mymodel import SpatioTemporalAutoEncoder_V2, count_parameters


def make_ae_dataset(train_root: str, dataset_kind: str, batchsize: int):
    if dataset_kind == "processed":
        return AE_NeuropixelsDataset(train_root, batch_size=batchsize)
    if dataset_kind == "cortexlab":
        return AE_NeuropixelsDataset_cortexlab(root=train_root, batch_size=batchsize)
    raise ValueError(f"Unknown dataset kind: {dataset_kind}")


def split_lengths(n_items: int) -> List[int]:
    """Use the original 90/5/5 split, with guards for small smoke-test data."""

    if n_items <= 0:
        raise ValueError("Dataset is empty.")
    if n_items == 1:
        return [1, 0, 0]
    if n_items == 2:
        return [1, 1, 0]

    val_size = max(1, int(0.05 * n_items))
    test_size = max(1, int(0.05 * n_items))
    train_size = n_items - val_size - test_size

    if train_size < 1:
        train_size = 1
        if test_size > 1:
            test_size -= 1
        else:
            val_size -= 1
    return [train_size, val_size, test_size]


def make_ae_loaders(
    train_root: str,
    dataset_kind: str,
    batchsize: int,
    num_workers: int,
    seed: int,
):
    train_root, dataset_kind = resolve_dataset_root(train_root, dataset_kind)
    dataset = make_ae_dataset(train_root, dataset_kind, batchsize)
    lengths = split_lengths(len(dataset))
    generator = torch.Generator().manual_seed(seed)
    train_dataset, val_dataset, test_dataset = random_split(
        dataset,
        lengths,
        generator=generator,
    )

    def make_loader(subset, shuffle: bool):
        if len(subset) == 0:
            return None
        return DataLoader(
            subset,
            batch_size=batchsize,
            shuffle=shuffle,
            num_workers=num_workers,
            pin_memory=torch.cuda.is_available(),
        )

    return (
        make_loader(train_dataset, shuffle=True),
        make_loader(val_dataset, shuffle=False),
        make_loader(test_dataset, shuffle=False),
        lengths,
    )


def _move_waveforms(tensor: torch.Tensor, device: torch.device) -> torch.Tensor:
    return tensor.to(device=device, dtype=torch.float64, non_blocking=True)


def train_ae_epoch(
    epoch: int,
    model: SpatioTemporalAutoEncoder_V2,
    optimizer: optim.Optimizer,
    train_loader: DataLoader,
    ae_loss: AELoss,
    writer: SummaryWriter,
    device: torch.device,
    log_every: int,
    progress: bool,
) -> float:
    model.train()
    losses = AverageMeter()
    iterator = tqdm.tqdm(
        train_loader,
        desc=f"Epoch {epoch:3d} AE train",
        disable=not progress,
    )

    for step, data in enumerate(iterator, start=1):
        data = _move_waveforms(data, device)
        bsz = data.shape[0]

        optimizer.zero_grad(set_to_none=True)
        reconstruction = model(data)
        loss = ae_loss(reconstruction, data)
        losses.update(loss.item(), bsz)
        loss.backward()
        optimizer.step()

        iteration = epoch * len(train_loader) + step
        if log_every > 0 and iteration % log_every == 0:
            writer.add_scalar("Train/Loss", losses.avg, iteration)
        if progress:
            iterator.set_postfix(loss=f"{losses.avg:.6f}")

    return float(losses.avg)


def evaluate_ae(
    epoch: int,
    model: SpatioTemporalAutoEncoder_V2,
    loader: Optional[DataLoader],
    ae_loss: AELoss,
    writer: Optional[SummaryWriter],
    device: torch.device,
    split_name: str,
    progress: bool,
) -> float:
    if loader is None:
        return float("nan")

    model.eval()
    losses = AverageMeter()
    iterator = tqdm.tqdm(
        loader,
        desc=f"Epoch {epoch:3d} AE {split_name}",
        disable=not progress,
    )

    with torch.no_grad():
        for data in iterator:
            data = _move_waveforms(data, device)
            bsz = data.shape[0]
            reconstruction = model(data)
            loss = ae_loss(reconstruction, data)
            losses.update(loss.item(), bsz)
            if progress:
                iterator.set_postfix(loss=f"{losses.avg:.6f}")

    if writer is not None:
        writer.add_scalar(f"{split_name.capitalize()}/Loss", losses.avg, epoch)
    return float(losses.avg)


def save_ae_checkpoint(
    ckpt_dir: Path,
    epoch: int,
    model: SpatioTemporalAutoEncoder_V2,
    optimizer: optim.Optimizer,
    n_output: int,
    parameter_counts: Dict[str, int],
) -> Path:
    ckpt_dir.mkdir(parents=True, exist_ok=True)
    checkpoint_path = ckpt_dir / f"ckpt_epoch_{epoch}"
    torch.save(
        {
            "model": model.state_dict(),
            "encoder": model.encoder.state_dict(),
            "optimizer": optimizer.state_dict(),
            "epoch": epoch,
            "n_output": n_output,
            "parameter_counts": parameter_counts,
        },
        checkpoint_path,
    )
    return checkpoint_path


def load_ae_checkpoint(
    checkpoint_path: Path,
    model: SpatioTemporalAutoEncoder_V2,
    optimizer: optim.Optimizer,
    device: torch.device,
) -> int:
    checkpoint = torch.load(checkpoint_path, map_location=device)
    model.load_state_dict(checkpoint["model"])
    optimizer.load_state_dict(checkpoint["optimizer"])
    return int(checkpoint["epoch"]) + 1


def _best_val_loss(rows: Sequence[Dict]) -> Dict:
    valid_rows = [row for row in rows if not np.isnan(float(row["val_loss"]))]
    if not valid_rows:
        return rows[-1] if rows else {}
    return min(valid_rows, key=lambda row: float(row["val_loss"]))


def train_ae_setting(
    *,
    exp_dir: Path,
    train_root: str,
    dataset_kind: str,
    n_output: int,
    total_epoch: int,
    batchsize: int,
    lr: float,
    save_freq: int,
    seed: int,
    device: torch.device,
    num_workers: int,
    deterministic: bool = False,
    resume: bool = False,
    log_every: int = 50,
    progress: bool = True,
) -> Dict:
    set_seed(seed, deterministic=deterministic)
    train_root, dataset_kind = resolve_dataset_root(train_root, dataset_kind)
    exp_dir.mkdir(parents=True, exist_ok=True)
    ckpt_dir = exp_dir / "ckpt"
    log_dir = exp_dir / "log"
    metrics_path = exp_dir / "metrics.csv"

    train_loader, val_loader, test_loader, lengths = make_ae_loaders(
        train_root=train_root,
        dataset_kind=dataset_kind,
        batchsize=batchsize,
        num_workers=num_workers,
        seed=seed,
    )
    if train_loader is None:
        raise ValueError("AE training split is empty.")

    model = SpatioTemporalAutoEncoder_V2(
        n_channel=30,
        n_time=60,
        n_output=n_output,
    ).to(device)
    model = model.double()

    parameter_counts = {
        "encoder_parameters": count_parameters(model.encoder),
        "decoder_parameters": count_parameters(model.decoder),
        "autoencoder_parameters": count_parameters(model),
    }
    ae_loss = AELoss(lambda1=0.0, lambda2=1.0).to(device)
    optimizer = optim.Adam(model.parameters(), lr=lr)

    write_json(
        exp_dir / "config.json",
        {
            "train_root": train_root,
            "dataset_kind": dataset_kind,
            "n_output": n_output,
            "total_epoch": total_epoch,
            "batchsize": batchsize,
            "lr": lr,
            "save_freq": save_freq,
            "seed": seed,
            "device": str(device),
            "num_workers": num_workers,
            "split_lengths": {
                "train": lengths[0],
                "val": lengths[1],
                "test": lengths[2],
            },
            **parameter_counts,
        },
    )

    start_epoch = 0
    if resume and ckpt_dir.exists():
        existing_ckpts = [path for path in ckpt_dir.iterdir() if path.is_file()]
        if existing_ckpts:
            start_epoch = load_ae_checkpoint(
                latest_checkpoint(ckpt_dir),
                model=model,
                optimizer=optimizer,
                device=device,
            )

    writer = SummaryWriter(log_dir=str(log_dir))
    rows: List[Dict] = []

    if total_epoch == 0 and start_epoch == 0:
        val_loss = evaluate_ae(
            0,
            model=model,
            loader=val_loader,
            ae_loss=ae_loss,
            writer=writer,
            device=device,
            split_name="validation",
            progress=progress,
        )
        checkpoint_path = save_ae_checkpoint(
            ckpt_dir=ckpt_dir,
            epoch=0,
            model=model,
            optimizer=optimizer,
            n_output=n_output,
            parameter_counts=parameter_counts,
        )
        rows.append(
            {
                "epoch": 0,
                "n_output": n_output,
                "train_loss": "",
                "val_loss": val_loss,
                "checkpoint": str(checkpoint_path),
            }
        )

    for epoch in range(start_epoch, total_epoch):
        train_loss = train_ae_epoch(
            epoch,
            model=model,
            optimizer=optimizer,
            train_loader=train_loader,
            ae_loss=ae_loss,
            writer=writer,
            device=device,
            log_every=log_every,
            progress=progress,
        )
        val_loss = evaluate_ae(
            epoch,
            model=model,
            loader=val_loader,
            ae_loss=ae_loss,
            writer=writer,
            device=device,
            split_name="validation",
            progress=progress,
        )

        checkpoint_path = ""
        if save_freq > 0 and (epoch % save_freq == 0 or epoch == total_epoch - 1):
            checkpoint_path = str(
                save_ae_checkpoint(
                    ckpt_dir=ckpt_dir,
                    epoch=epoch,
                    model=model,
                    optimizer=optimizer,
                    n_output=n_output,
                    parameter_counts=parameter_counts,
                )
            )

        rows.append(
            {
                "epoch": epoch,
                "n_output": n_output,
                "train_loss": train_loss,
                "val_loss": val_loss,
                "checkpoint": checkpoint_path,
            }
        )
        write_rows_csv(metrics_path, rows)

    test_epoch = int(rows[-1]["epoch"]) if rows else max(start_epoch - 1, 0)
    test_loss = evaluate_ae(
        test_epoch,
        model=model,
        loader=test_loader,
        ae_loss=ae_loss,
        writer=writer,
        device=device,
        split_name="test",
        progress=progress,
    )
    writer.close()
    write_rows_csv(metrics_path, rows)

    if rows:
        best = _best_val_loss(rows)
        final = rows[-1]
        last_checkpoint = final.get("checkpoint") or str(latest_checkpoint(ckpt_dir))
    else:
        best = {}
        final = {}
        last_checkpoint = str(latest_checkpoint(ckpt_dir))

    return {
        "n_output": n_output,
        **parameter_counts,
        "ae_best_epoch": best.get("epoch", ""),
        "ae_best_val_loss": best.get("val_loss", ""),
        "ae_final_epoch": final.get("epoch", ""),
        "ae_final_train_loss": final.get("train_loss", ""),
        "ae_final_val_loss": final.get("val_loss", ""),
        "ae_test_loss": test_loss,
        "ae_checkpoint": last_checkpoint,
        "ae_experiment_dir": str(exp_dir),
    }


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Sweep SpatioTemporalAutoEncoder_V2 latent feature dimensions."
    )
    parser.add_argument("--exp_name", type=str, default="AE_params")
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
        choices=("auto", "cortexlab", "processed"),
    )
    parser.add_argument(
        "--n_outputs",
        type=int,
        nargs="+",
        default=[8, 32, 128, 256],
        # default=[256,],
    )
    parser.add_argument("--ae_total_epoch", type=int, default=300) # 300
    parser.add_argument("--finetune_total_epoch", type=int, default=50) # 50
    parser.add_argument("--ae_batchsize", type=int, default=32)
    parser.add_argument("--finetune_batchsize", type=int, default=40)
    parser.add_argument("--ae_lr", type=float, default=1e-5)
    parser.add_argument("--lr_enc", type=float, default=2e-5)
    parser.add_argument("--lr_proj", type=float, default=1.1e-4)
    parser.add_argument("--negative_weight", type=float, default=10.0)
    parser.add_argument("--ae_save_freq", type=int, default=1)
    parser.add_argument("--finetune_save_freq", type=int, default=1)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--device", type=str, default="auto")
    parser.add_argument("--num_workers", type=int, default=0)
    parser.add_argument(
        "--projector_output_dim",
        type=int,
        default=None,
        help=(
            "Projection-head output dimension. If omitted, uses "
            "min(n_output, 128) separately for each swept n_output."
        ),
    )
    # parser.add_argument("--projector_hidden_dim", type=int, default=128)
    parser.add_argument(
        "--projector_hidden_dim",
        type=int,
        default=None,
        help=(
            "Projection-head hidden dimension. If omitted, uses "
            "min(n_output, 128) separately for each swept n_output."
        ),
    )
    
    parser.add_argument("--projector_dropout", type=float, default=0.1)
    parser.add_argument("--validation_mode", type=str, default="val", choices=("train", "val"))
    parser.add_argument("--log_every", type=int, default=50)
    parser.add_argument("--skip_finetune", action="store_true")
    parser.add_argument("--train_all_encoder", action="store_true")
    parser.add_argument("--deterministic", action="store_true")
    parser.add_argument("--cont", action="store_true")
    parser.add_argument("--overwrite", action="store_true")
    parser.add_argument("--no_progress", action="store_true")
    parser.add_argument(
        "--output_dir",
        type=str,
        default=str(REBUTTAL_DIR / "results" / "AE_params"),
    )
    return parser


def main(argv: Optional[Sequence[str]] = None) -> None:
    parser = build_arg_parser()
    args = parser.parse_args(argv)

    device = resolve_device(args.device)
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
            "train_root_resolved": train_root,
            "dataset_kind_resolved": dataset_kind,
        },
    )

    summaries: List[Dict] = []
    for n_output in args.n_outputs:
        setting_dir = exp_root / f"n_output_{n_output}"
        ae_dir = setting_dir / "ae"
        finetune_dir = setting_dir / "finetune"
        print(f"\nRunning n_output={n_output} in {setting_dir}")

        ae_summary = train_ae_setting(
            exp_dir=ae_dir,
            train_root=train_root,
            dataset_kind=dataset_kind,
            n_output=n_output,
            total_epoch=args.ae_total_epoch,
            batchsize=args.ae_batchsize,
            lr=args.ae_lr,
            save_freq=args.ae_save_freq,
            seed=args.seed,
            device=device,
            num_workers=args.num_workers,
            deterministic=args.deterministic,
            resume=args.cont,
            log_every=args.log_every,
            progress=not args.no_progress,
        )

        combined = dict(ae_summary)
        if not args.skip_finetune:
            projector_output_dim = (
                min(n_output, 128)
                if args.projector_output_dim is None
                else args.projector_output_dim
            )
            if projector_output_dim > n_output:
                raise ValueError(
                    "projector_output_dim must not be larger than n_output. "
                    f"Got projector_output_dim={projector_output_dim}, "
                    f"n_output={n_output}."
                )
            projector_hidden_dim = (
                min(n_output, 128)
                if args.projector_hidden_dim is None
                else args.projector_hidden_dim
            )
            if projector_hidden_dim > n_output:
                raise ValueError(
                    "projector_hidden_dim must not be larger than n_output. "
                    f"Got projector_hidden_dim={projector_hidden_dim}, "
                    f"n_output={n_output}."
                )
                
            finetune_summary = train_contrastive_setting(
                exp_dir=finetune_dir, 
                train_root=train_root,
                dataset_kind=dataset_kind,
                ae_checkpoint=Path(ae_summary["ae_checkpoint"]),
                n_output=n_output,
                negative_weight=args.negative_weight,
                total_epoch=args.finetune_total_epoch,
                batchsize=args.finetune_batchsize,
                lr_enc=args.lr_enc,
                lr_proj=args.lr_proj,
                save_freq=args.finetune_save_freq,
                seed=args.seed,
                device=device,
                num_workers=args.num_workers,
                freeze_non_fc=not args.train_all_encoder,
                projector_output_dim=projector_output_dim,
                projector_hidden_dim=projector_hidden_dim,
                projector_dropout=args.projector_dropout,
                deterministic=args.deterministic,
                resume=args.cont,
                validation_mode=args.validation_mode,
                log_every=args.log_every,
                progress=not args.no_progress,
            )
            combined.update(
                {
                    "finetune_negative_weight": args.negative_weight,
                    "finetune_projector_output_dim": projector_output_dim,
                    "finetune_projector_hidden_dim": projector_hidden_dim,
                    "finetune_best_encoder_epoch": finetune_summary.get("best_encoder_epoch", ""),
                    "finetune_best_val_encoder_accuracy": finetune_summary.get(
                        "best_val_encoder_accuracy", ""
                    ),
                    "finetune_best_projector_epoch": finetune_summary.get(
                        "best_projector_epoch", ""
                    ),
                    "finetune_best_val_projector_accuracy": finetune_summary.get(
                        "best_val_projector_accuracy", ""
                    ),
                    "finetune_final_epoch": finetune_summary.get("final_epoch", ""),
                    "finetune_final_train_loss": finetune_summary.get("final_train_loss", ""),
                    "finetune_final_val_loss": finetune_summary.get("final_val_loss", ""),
                    "finetune_final_val_encoder_accuracy": finetune_summary.get(
                        "final_val_encoder_accuracy", ""
                    ),
                    "finetune_final_val_projector_accuracy": finetune_summary.get(
                        "final_val_projector_accuracy", ""
                    ),
                    "finetune_checkpoint": finetune_summary.get("last_checkpoint", ""),
                    "finetune_experiment_dir": finetune_summary.get("experiment_dir", ""),
                }
            )

        summaries.append(combined)
        write_rows_csv(exp_root / "summary.csv", summaries)

    print(f"\nSaved summary to {exp_root / 'summary.csv'}")


if __name__ == "__main__":
    main()
