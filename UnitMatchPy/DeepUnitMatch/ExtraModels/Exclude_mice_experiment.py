"""Train DeepUnitMatch with 1/18, 6/18, or 12/18 of the mice.

The experiment follows ``AE_params_experiment.py``: each mouse subset is used
for autoencoder pretraining and, unless disabled, contrastive finetuning.  The
latent feature dimension is fixed at the original value of 256 so that dataset
amount is the only swept experimental factor.

Only the nested CortexLab layout is supported because filtering is performed at
the ``root/mouse/probe/location/experiment`` level.
"""

from __future__ import annotations

import argparse
from contextlib import contextmanager
from pathlib import Path
import sys
from typing import Dict, Iterator, List, Optional, Sequence

import torch
from torch.utils.data import DataLoader, random_split


REBUTTAL_DIR = Path(__file__).resolve().parent
REPO_ROOT = REBUTTAL_DIR.parent
for import_path in (REBUTTAL_DIR, REPO_ROOT):
    import_path_str = str(import_path)
    if import_path_str not in sys.path:
        sys.path.insert(0, import_path_str)

import AE_params_experiment as ae_experiment
import loss_sensitivity_experiment as contrastive_experiment

from DeepUnitMatch.utils.AE_npdataset import AE_NeuropixelsDataset_cortexlab
from DeepUnitMatch.utils.npdataset import (
    TrainExperimentBatchSampler,
    ValidationExperimentBatchSampler,
)
from DeepUnitMatch.utils.train_utils import read_good_ids


N_OUTPUT = 256

mice_names = [
    "AL031",
    "AL032",
    "AL036",
    "AV008",
    "AV009",
    "AV015",
    "AV021",
    "AV049",
    "CB015",
    "CB016",
    "CB017",
    "CB018",
    "CB020",
    "EB014",
    "EB019",
    "FT033",
    "FT039",
    "JF084",
]

dict_1mouse = {
    "model_m1_1": ["AL032"],
    "model_m1_2": ["AV015"],
    "model_m1_3": ["CB016"],
}

dict_6mice = {
    "Model_m6_1": ["AL036", "AV008", "AV021", "CB015", "FT033", "FT039"],
    "Model_m6_2": ["AL031", "AV008", "CB018", "CB020", "EB014", "JF084"],
    "Model_m6_3": ["AL032", "AV015", "AV049", "CB016", "CB017", "EB019"],
}

dict_12mice = {
    "Model_m12_1": [
        "AL036",
        "AV009",
        "AV021",
        "CB015",
        "FT033",
        "FT039",
        "AL031",
        "AV008",
        "CB018",
        "CB020",
        "EB014",
        "JF084",
    ],
    "Model_m12_2": [
        "AL031",
        "AV009",
        "CB018",
        "CB020",
        "EB014",
        "JF084",
        "AL032",
        "AV015",
        "AV049",
        "CB016",
        "CB017",
        "EB019",
    ],
    "Model_m12_3": [
        "AL036",
        "AV009",
        "AV021",
        "CB015",
        "FT033",
        "FT039",
        "AL032",
        "AV015",
        "AV049",
        "CB016",
        "CB017",
        "EB019",
    ],
}

MOUSE_SPLITS = {
    1: dict_1mouse,
    6: dict_6mice,
    12: dict_12mice,
}


def _all_model_names() -> List[str]:
    return [
        model_name
        for split in MOUSE_SPLITS.values()
        for model_name in split
    ]


def _selected_models(
    requested_model_names: Optional[Sequence[str]],
) -> Iterator[tuple[str, List[str]]]:
    requested = set(requested_model_names or _all_model_names())
    unknown = requested.difference(_all_model_names())
    if unknown:
        raise ValueError(
            "Unknown model name(s): "
            f"{', '.join(sorted(unknown))}. Valid names are "
            f"{', '.join(_all_model_names())}."
        )

    for split in MOUSE_SPLITS.values():
        for model_name, selected_mice in split.items():
            if model_name in requested:
                yield model_name, list(selected_mice)


def validate_mouse_splits() -> None:
    """Validate the hard-coded rebuttal splits before starting a long run."""

    roster = set(mice_names)
    if len(roster) != 18 or len(roster) != len(mice_names):
        raise ValueError("mice_names must contain 18 unique mouse names.")

    seen_model_names = set()
    for expected_count, split in MOUSE_SPLITS.items():
        if len(split) != 3:
            raise ValueError(
                f"The {expected_count}-mouse case must contain three models."
            )
        for model_name, selected_mice in split.items():
            if model_name in seen_model_names:
                raise ValueError(f"Duplicate model name: {model_name}")
            seen_model_names.add(model_name)

            if len(selected_mice) != expected_count:
                raise ValueError(
                    f"{model_name} must contain {expected_count} mice, "
                    f"but contains {len(selected_mice)}."
                )
            if len(set(selected_mice)) != len(selected_mice):
                raise ValueError(f"{model_name} contains duplicate mice.")

            unknown_mice = set(selected_mice).difference(roster)
            if unknown_mice:
                raise ValueError(
                    f"{model_name} contains unknown mice: "
                    f"{', '.join(sorted(unknown_mice))}."
                )


def validate_dataset_mice(train_root: str, selected_mice: Sequence[str]) -> None:
    root = Path(train_root)
    missing = [
        mouse_name
        for mouse_name in selected_mice
        if not (root / mouse_name).is_dir()
    ]
    if missing:
        raise FileNotFoundError(
            f"Mouse directories missing from {root}: {', '.join(missing)}"
        )


class MouseFilteredAEDataset(AE_NeuropixelsDataset_cortexlab):
    """AE dataset restricted to explicit mouse directories."""

    def __init__(
        self,
        root: str,
        mouse_names: Sequence[str],
        batch_size: int = 32,
    ):
        self.root = str(Path(root).resolve())
        self.mouse_names = list(mouse_names)
        self.batch_size = batch_size
        self.np_file_names = read_good_ids(
            self.root,
            self.batch_size,
            self.mouse_names,
            finetune=False,
        )
        self.n_neurons = len(self.np_file_names)
        if self.n_neurons < 1:
            raise ValueError(
                "No AE waveforms were found for mice "
                f"{', '.join(self.mouse_names)} with batch size {batch_size}."
            )


class MouseFilteredPairDataset(
    contrastive_experiment.NestedCortexlabPairDataset
):
    """Contrastive dataset restricted to explicit mouse directories."""

    def __init__(
        self,
        root: str,
        mouse_names: Sequence[str],
        batch_size: int = 40,
        mode: str = "val",
    ):
        self.root = Path(root).resolve()
        self.batch_size = batch_size
        self.mode = mode
        self.mouse_names = list(mouse_names)

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

        if not self.all_files:
            raise ValueError(
                "No contrastive waveforms were found for mice "
                f"{', '.join(self.mouse_names)} with batch size {batch_size}."
            )
        print(
            f"Initialised mouse-filtered contrastive dataset with "
            f"{len(self.all_files)} files across "
            f"{len(self.experiment_unit_map)} experiments from "
            f"{len(self.mouse_names)} mice."
        )


def make_mouse_filtered_ae_loaders(
    train_root: str,
    dataset_kind: str,
    batchsize: int,
    num_workers: int,
    seed: int,
    selected_mice: Sequence[str],
):
    train_root, dataset_kind = contrastive_experiment.resolve_dataset_root(
        train_root,
        dataset_kind,
    )
    if dataset_kind != "cortexlab":
        raise ValueError(
            "Mouse exclusion requires the nested cortexlab dataset layout."
        )

    validate_dataset_mice(train_root, selected_mice)
    dataset = MouseFilteredAEDataset(
        train_root,
        mouse_names=selected_mice,
        batch_size=batchsize,
    )
    lengths = ae_experiment.split_lengths(len(dataset))
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


def make_mouse_filtered_pair_loaders(
    train_root: str,
    dataset_kind: str,
    batchsize: int,
    num_workers: int,
    validation_mode: str,
    selected_mice: Sequence[str],
):
    train_root, dataset_kind = contrastive_experiment.resolve_dataset_root(
        train_root,
        dataset_kind,
    )
    if dataset_kind != "cortexlab":
        raise ValueError(
            "Mouse exclusion requires the nested cortexlab dataset layout."
        )

    validate_dataset_mice(train_root, selected_mice)
    train_dataset = MouseFilteredPairDataset(
        train_root,
        mouse_names=selected_mice,
        batch_size=batchsize,
        mode="train",
    )
    val_dataset = MouseFilteredPairDataset(
        train_root,
        mouse_names=selected_mice,
        batch_size=batchsize,
        mode=validation_mode,
    )

    train_sampler = TrainExperimentBatchSampler(
        train_dataset,
        batchsize,
        shuffle=True,
    )
    val_sampler = ValidationExperimentBatchSampler(
        val_dataset,
        shuffle=False,
    )

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


@contextmanager
def use_mouse_filtered_loaders(
    selected_mice: Sequence[str],
) -> Iterator[None]:
    """Inject filtered loaders into the two established training routines."""

    original_ae_loader = ae_experiment.make_ae_loaders
    original_pair_loader = contrastive_experiment.make_paired_loaders

    def make_ae_loaders(
        train_root: str,
        dataset_kind: str,
        batchsize: int,
        num_workers: int,
        seed: int,
    ):
        return make_mouse_filtered_ae_loaders(
            train_root=train_root,
            dataset_kind=dataset_kind,
            batchsize=batchsize,
            num_workers=num_workers,
            seed=seed,
            selected_mice=selected_mice,
        )

    def make_paired_loaders(
        train_root: str,
        dataset_kind: str,
        batchsize: int,
        num_workers: int,
        validation_mode: str,
    ):
        return make_mouse_filtered_pair_loaders(
            train_root=train_root,
            dataset_kind=dataset_kind,
            batchsize=batchsize,
            num_workers=num_workers,
            validation_mode=validation_mode,
            selected_mice=selected_mice,
        )

    ae_experiment.make_ae_loaders = make_ae_loaders
    contrastive_experiment.make_paired_loaders = make_paired_loaders
    try:
        yield
    finally:
        ae_experiment.make_ae_loaders = original_ae_loader
        contrastive_experiment.make_paired_loaders = original_pair_loader


def train_mouse_setting(
    *,
    model_name: str,
    selected_mice: Sequence[str],
    setting_dir: Path,
    train_root: str,
    dataset_kind: str,
    args: argparse.Namespace,
    device: torch.device,
) -> Dict:
    """Train the AE and optional contrastive model for one mouse subset."""

    mouse_count = len(selected_mice)
    selection = {
        "model_name": model_name,
        "mouse_count": mouse_count,
        "dataset_fraction": f"{mouse_count}/18",
        "mice": list(selected_mice),
        "excluded_mice": [
            mouse_name
            for mouse_name in mice_names
            if mouse_name not in selected_mice
        ],
        "n_output": N_OUTPUT,
    }
    contrastive_experiment.write_json(
        setting_dir / "mouse_selection.json",
        selection,
    )

    ae_dir = setting_dir / "ae"
    finetune_dir = setting_dir / "finetune"
    with use_mouse_filtered_loaders(selected_mice):
        ae_summary = ae_experiment.train_ae_setting(
            exp_dir=ae_dir,
            train_root=train_root,
            dataset_kind=dataset_kind,
            n_output=N_OUTPUT,
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

        combined = {
            "model_name": model_name,
            "mouse_count": mouse_count,
            "dataset_fraction": f"{mouse_count}/18",
            "mice": ",".join(selected_mice),
            **ae_summary,
        }
        if args.skip_finetune:
            return combined

        finetune_summary = contrastive_experiment.train_contrastive_setting(
            exp_dir=finetune_dir,
            train_root=train_root,
            dataset_kind=dataset_kind,
            ae_checkpoint=Path(ae_summary["ae_checkpoint"]),
            n_output=N_OUTPUT,
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
            projector_output_dim=args.projector_output_dim,
            projector_hidden_dim=args.projector_hidden_dim,
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
            "finetune_projector_output_dim": args.projector_output_dim,
            "finetune_projector_hidden_dim": args.projector_hidden_dim,
            "finetune_best_encoder_epoch": finetune_summary.get(
                "best_encoder_epoch",
                "",
            ),
            "finetune_best_val_encoder_accuracy": finetune_summary.get(
                "best_val_encoder_accuracy",
                "",
            ),
            "finetune_best_projector_epoch": finetune_summary.get(
                "best_projector_epoch",
                "",
            ),
            "finetune_best_val_projector_accuracy": finetune_summary.get(
                "best_val_projector_accuracy",
                "",
            ),
            "finetune_final_epoch": finetune_summary.get("final_epoch", ""),
            "finetune_final_train_loss": finetune_summary.get(
                "final_train_loss",
                "",
            ),
            "finetune_final_val_loss": finetune_summary.get(
                "final_val_loss",
                "",
            ),
            "finetune_final_val_encoder_accuracy": finetune_summary.get(
                "final_val_encoder_accuracy",
                "",
            ),
            "finetune_final_val_projector_accuracy": finetune_summary.get(
                "final_val_projector_accuracy",
                "",
            ),
            "finetune_checkpoint": finetune_summary.get(
                "last_checkpoint",
                "",
            ),
            "finetune_experiment_dir": finetune_summary.get(
                "experiment_dir",
                "",
            ),
        }
    )
    return combined


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Train n_output=256 models with the predefined 1-, 6-, and "
            "12-mouse dataset subsets."
        )
    )
    parser.add_argument("--exp_name", type=str, default="dataset_amount")
    parser.add_argument(
        "--train_root",
        type=str,
        default=str(contrastive_experiment.DEFAULT_DATA_ROOT),
        help=(
            "Nested dataset root. The outer Data_DUM folder and its "
            "ALL_DATA_nomatchtables child are both accepted."
        ),
    )
    parser.add_argument(
        "--dataset_kind",
        type=str,
        default="cortexlab",
        choices=("auto", "cortexlab"),
    )
    parser.add_argument(
        "--model_names",
        type=str,
        nargs="+",
        default=None,
        help=(
            "Optional subset of the nine predefined model names. By default, "
            "all three replicates for all three dataset amounts are trained."
        ),
    )
    parser.add_argument("--ae_total_epoch", type=int, default=300)
    parser.add_argument("--finetune_total_epoch", type=int, default=50)
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
    parser.add_argument("--projector_output_dim", type=int, default=128)
    parser.add_argument("--projector_hidden_dim", type=int, default=128)
    parser.add_argument("--projector_dropout", type=float, default=0.1)
    parser.add_argument(
        "--validation_mode",
        type=str,
        default="val",
        choices=("train", "val"),
    )
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
        default=str(REBUTTAL_DIR / "results" / "Exclude_mice"),
    )
    return parser


def main(argv: Optional[Sequence[str]] = None) -> None:
    parser = build_arg_parser()
    args = parser.parse_args(argv)
    validate_mouse_splits()
    selected_models = list(_selected_models(args.model_names))
    if args.projector_output_dim > N_OUTPUT:
        raise ValueError(
            "projector_output_dim must not be larger than n_output. "
            f"Got projector_output_dim={args.projector_output_dim}, "
            f"n_output={N_OUTPUT}."
        )
    if args.projector_hidden_dim > N_OUTPUT:
        raise ValueError(
            "projector_hidden_dim must not be larger than n_output. "
            f"Got projector_hidden_dim={args.projector_hidden_dim}, "
            f"n_output={N_OUTPUT}."
        )

    device = contrastive_experiment.resolve_device(args.device)
    train_root, dataset_kind = contrastive_experiment.resolve_dataset_root(
        args.train_root,
        args.dataset_kind,
    )
    if dataset_kind != "cortexlab":
        raise ValueError(
            "Mouse exclusion requires the nested cortexlab dataset layout."
        )
    validate_dataset_mice(
        train_root,
        [
            mouse_name
            for _, selected_mice in selected_models
            for mouse_name in selected_mice
        ],
    )

    exp_root = Path(args.output_dir).expanduser()
    if not exp_root.is_absolute():
        exp_root = (Path.cwd() / exp_root).resolve()
    exp_root = exp_root / args.exp_name

    if (
        exp_root.exists()
        and any(exp_root.iterdir())
        and not args.cont
        and not args.overwrite
    ):
        raise ValueError(
            f"Output directory already exists: {exp_root}. "
            "Use --cont to resume or --overwrite to write new metrics into it."
        )

    exp_root.mkdir(parents=True, exist_ok=True)
    contrastive_experiment.write_json(
        exp_root / "sweep_config.json",
        {
            **vars(args),
            "n_output": N_OUTPUT,
            "device_resolved": str(device),
            "train_root_resolved": train_root,
            "dataset_kind_resolved": dataset_kind,
            "all_mice": mice_names,
            "mouse_splits": MOUSE_SPLITS,
        },
    )

    summaries: List[Dict] = []
    for model_name, selected_mice in selected_models:
        mouse_count = len(selected_mice)
        setting_dir = exp_root / f"{mouse_count}_mice" / model_name
        print(
            f"\nRunning {model_name}: {mouse_count}/18 mice "
            f"({', '.join(selected_mice)}) in {setting_dir}"
        )
        summary = train_mouse_setting(
            model_name=model_name,
            selected_mice=selected_mice,
            setting_dir=setting_dir,
            train_root=train_root,
            dataset_kind=dataset_kind,
            args=args,
            device=device,
        )
        summaries.append(summary)
        contrastive_experiment.write_rows_csv(
            exp_root / "summary.csv",
            summaries,
        )

    print(f"\nSaved summary to {exp_root / 'summary.csv'}")


if __name__ == "__main__":
    main()
