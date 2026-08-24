"""End-to-end mouse-count cross-validation: train + immediately evaluate.

Redo of Exclude_mice_experiment.py (which trained on {1, 6, 12}-mouse subsets
and left inference to a separate discovery-based script,
PaperAnalyses/run_deepunitmatch_batch_onMerged_ExcludeMiceModels.py), with
three differences:
  - 3 replicates per group size (including 1-mouse, matching
    Exclude_mice_experiment.py's original single-mouse subsets), plus a new
    18-mouse (all mice, single replicate) control:
    3x1 + 3x3 + 3x6 + 3x12 + 1x18 = 13 models total.
  - Training and inference are bundled into a single resumable pipeline per
    model ("xval") instead of two separate scripts.
  - The whole thing is safe to launch as several concurrent python processes
    (same machine or several machines) pointed at the same shared drive: each
    of the 13 xvals is claimed via batch_lock.try_lock() before a process
    starts it, so nobody duplicates another instance's run.

Data
----
Training snippets are derived from the same *merged* tree used for inference
(DeepUM_NatMeth2026V2_merged/merged_data_v2), not DeepUM_NatMeth2026V2 itself:
that root is the pre-merge UnitMatch.mat/UMparam.pickle tree
generate_merged_dataset.py reads to *build* the merged tree, and has no
metadata.json/processed_waveforms per experiment for the AE/finetune dataset
loaders to read. The merged tree is already known-good for this purpose --
PaperAnalyses/train_deepunitmatch_from_merged.py preprocesses and trains from
it today. The preprocessing stage below mirrors that script, but uses the v1
(non-_v2) param_fun/npdataset/losses/mymodel modules -- the same ones
Exclude_mice_experiment.py and run_deepunitmatch_batch_onMerged.py's inference
use -- so the resulting checkpoints stay compatible with
DeepUnitMatch.testing.test (not test_v2).

Preprocessed snippets are cached once per merged-data location under a shared
folder (PREPROCESSED_CACHE, a new sibling of merged_data_v2) and reused across
every xval that needs that location's mouse -- there's heavy overlap between
the 10 mouse subsets, so this avoids redoing the same extraction 10 times.

Parallelism / "main log"
-------------------------
Coordination state (per-xval status.json + the outer per-xval lock) lives
under SHARED_STATE_ROOT on the same shared drive, not the local checkout --
multiple *machines* need to see the same state for the "don't double-start an
xval" guarantee to hold. Run with --status at any time to print every xval's
state (queued/preprocessing/training_ae/training_finetune/inference/done/
failed) -- that table, built by scanning each xval's own status.json, *is*
the main log the task asked for; there's no single shared mutable file that
would itself need locking on the network share.

Known limitation: training checkpoints themselves are written by
DeepUnitMatch/train/train_AE.py and train_finetune.py to a *local* path
(DeepUnitMatch/ModelExp/... on whichever machine runs the training) -- that's
inherent to those scripts and out of scope to change here. So while the outer
lock reliably prevents two machines from starting the same xval concurrently,
if a machine is interrupted mid-training and its lock later goes stale and is
reclaimed by a *different* machine, that machine restarts the xval from
scratch rather than resuming the first machine's partial progress. Resuming
mid-training only works if the same machine picks its own xval back up.

Example
-------
  python xval_end_to_end.py                  # process all 13 xvals, in order
  python xval_end_to_end.py --only m3_1 m6_2  # just these two
  python xval_end_to_end.py --status          # print progress, don't run anything
"""

from __future__ import annotations

import argparse
import datetime
import json
import os
import random
import socket
import sys
import time
import traceback
from pathlib import Path
from typing import Dict, List, Sequence

import numpy as np

import matplotlib

matplotlib.use("Agg")  # non-interactive backend, must be set before pyplot is imported

# ── project paths ───────────────────────────────────────────────────────────
_HERE = Path(__file__).resolve().parent               # .../DeepUnitMatch/ExtraModels
DEEPUNITMATCH_DIR = _HERE.parent                       # .../DeepUnitMatch
REPO_ROOT = DEEPUNITMATCH_DIR.parent                   # .../UnitMatchPy
PAPER_ANALYSES_DIR = REPO_ROOT / "PaperAnalyses"
for _p in (_HERE, DEEPUNITMATCH_DIR, REPO_ROOT, PAPER_ANALYSES_DIR):
    _ps = str(_p)
    if _ps not in sys.path:
        sys.path.insert(0, _ps)

import batch_lock
import run_deepunitmatch_batch_onMerged as base_batch
from torch.utils.data import ConcatDataset

from DeepUnitMatch.testing import test
from DeepUnitMatch.train import train_AE as train_ae_mod
from DeepUnitMatch.train import train_finetune as train_finetune_mod
from DeepUnitMatch.utils import npdataset
from DeepUnitMatch.utils import param_fun
from DeepUnitMatch.utils.AE_npdataset import AE_NeuropixelsDataset

# ── roster / sweep spec ──────────────────────────────────────────────────────

# Same 18-mouse roster Exclude_mice_experiment.py used (there, split unevenly
# across dict_1mouse/dict_6mice/dict_12mice with inconsistent key casing --
# collapsed here into one canonical, lowercase-safe list).
CANONICAL_MICE = [
    "AL031", "AL032", "AL036",
    "AV008", "AV009", "AV015", "AV021", "AV049",
    "CB015", "CB016", "CB017", "CB018", "CB020",
    "EB014", "EB019",
    "FT033", "FT039",
    "JF084",
]

GROUP_SPEC = [(1, 3), (3, 3), (6, 3), (12, 3), (18, 1)]  # (mice per model, replicates)
SPLIT_SEED = 0  # fixed so every process/machine derives the same manifest

# Hand-picked mice for the 1-mouse xvals (every other group size is drawn at
# random -- see generate_manifest) -- order gives m1_1/m1_2/m1_3 respectively.
SINGLE_MOUSE_XVAL_MICE = ["AL032", "AV009", "FT039"]
assert all(m in CANONICAL_MICE for m in SINGLE_MOUSE_XVAL_MICE)

N_OUTPUT = 256

# ── shared paths (network drive, so every machine agrees) ───────────────────
_MERGED_PARENT = Path(os.path.dirname(base_batch.BASE_INPUT))  # .../DeepUM_NatMeth2026V2_merged
PREPROCESSED_CACHE = _MERGED_PARENT / "xval_end_to_end_preprocessed_v1"
LOCKS_DIR = PREPROCESSED_CACHE / ".locks"
SHARED_STATE_ROOT = _MERGED_PARENT / "xval_end_to_end_state"

# One xval = AE pretrain (up to --ae_total_epoch epochs) + finetune (up to
# --finetune_total_epoch epochs) + inference over every merged-data group --
# realistically many hours to a few days. Much longer than batch_lock's
# default 24h stale window (sized for a single inference group), so a
# legitimately still-running xval isn't mistaken for an abandoned one.
TRAINING_STALE_AFTER_SECONDS = 5 * 24 * 3600


# ── mouse-split manifest ─────────────────────────────────────────────────────


def generate_manifest() -> "dict[str, List[str]]":
    """Deterministically assign mice to each of the 13 xval replicates.

    For group sizes 3/6/12, each replicate independently draws `group_size`
    mice at random (without replacement within a replicate; replicates may
    overlap each other, same as Exclude_mice_experiment.py's hand-picked
    splits did), seeded so every process/machine derives the identical
    manifest without needing to share it. The 18-mouse group has only one
    possible draw: every mouse. The 1-mouse group is the one exception --
    it uses the hand-picked SINGLE_MOUSE_XVAL_MICE instead of a random draw,
    and doesn't consume rng state, so it doesn't perturb the 3/6/12-mouse
    draws.
    """
    rng = random.Random(SPLIT_SEED)
    manifest: "dict[str, List[str]]" = {}
    for group_size, n_replicates in GROUP_SPEC:
        for replicate in range(1, n_replicates + 1):
            name = f"m{group_size}_{replicate}"
            if group_size == 1:
                mice = [SINGLE_MOUSE_XVAL_MICE[replicate - 1]]
            elif group_size >= len(CANONICAL_MICE):
                mice = list(CANONICAL_MICE)
            else:
                mice = rng.sample(CANONICAL_MICE, group_size)
            manifest[name] = sorted(mice)
    return manifest


# ── status / "main log" ──────────────────────────────────────────────────────


def xval_dir(name: str) -> Path:
    return SHARED_STATE_ROOT / name


def status_path(name: str) -> Path:
    return xval_dir(name) / "status.json"


def write_status(name: str, **fields) -> None:
    path = status_path(name)
    path.parent.mkdir(parents=True, exist_ok=True)
    payload = {}
    if path.is_file():
        try:
            payload = json.loads(path.read_text())
        except Exception:
            payload = {}
    payload.update(fields)
    payload["updated"] = datetime.datetime.now().isoformat(timespec="seconds")
    payload["host"] = socket.gethostname()
    payload["pid"] = os.getpid()
    path.write_text(json.dumps(payload, indent=2, sort_keys=True))


def read_status(name: str) -> dict:
    path = status_path(name)
    if not path.is_file():
        return {"state": "queued"}
    try:
        return json.loads(path.read_text())
    except Exception:
        return {"state": "unknown"}


def xval_is_done(name: str) -> bool:
    return read_status(name).get("state") == "done"


def print_status_table(manifest: "dict[str, List[str]]") -> None:
    print(f"{'xval':<8} {'mice':<6} {'state':<18} {'updated':<20} host")
    for name in manifest:
        st = read_status(name)
        print(
            f"{name:<8} {len(manifest[name]):<6} {st.get('state', 'queued'):<18} "
            f"{st.get('updated', ''):<20} {st.get('host', '')}"
        )


def xval_lock_path(name: str) -> str:
    return str(xval_dir(name) / ".processing.lock")


# ── merged-group discovery / mouse filtering ─────────────────────────────────


def group_by_mouse(groups: Sequence[str]) -> "dict[str, List[str]]":
    """merged_dir -> mouse name is the first path component under BASE_INPUT."""
    by_mouse: "dict[str, List[str]]" = {}
    for merged_dir in groups:
        rel = os.path.relpath(os.path.dirname(merged_dir), base_batch.BASE_INPUT)
        mouse = rel.split(os.sep)[0]
        by_mouse.setdefault(mouse, []).append(merged_dir)
    return by_mouse


# ── shared preprocessing (v1 snippets, cached per merged-data location) ─────


def get_location_out_dir(merged_dir: str) -> str:
    subfolder = os.path.relpath(os.path.dirname(merged_dir), base_batch.BASE_INPUT)
    return str(PREPROCESSED_CACHE / subfolder)


def location_done(merged_dir: str) -> bool:
    return os.path.isfile(os.path.join(get_location_out_dir(merged_dir), ".done"))


def get_location_lock_path(merged_dir: str) -> str:
    subfolder = os.path.relpath(os.path.dirname(merged_dir), base_batch.BASE_INPUT)
    safe = subfolder.replace(os.sep, "_").replace("/", "_")
    return str(LOCKS_DIR / f"{safe}.lock")


def preprocess_one_location(merged_dir: str) -> None:
    out_dir = get_location_out_dir(merged_dir)
    sess = base_batch._prepare_session(merged_dir)
    if sess is None:
        raise RuntimeError(f"could not load session data for {merged_dir}")
    unit_ids = np.concatenate(sess["param"]["good_units"]).squeeze()
    os.makedirs(out_dir, exist_ok=True)
    param_fun.get_snippets(
        sess["waveform"],
        sess["channel_pos"],
        sess["session_id"],
        save_path=out_dir,
        unit_ids=unit_ids,
        param=sess["param"],
    )
    (Path(out_dir) / ".done").write_text("ok")


def ensure_locations_preprocessed(
    merged_dirs: Sequence[str], max_wait_s: int, poll_s: int = 30
) -> bool:
    """
    Make sure every location in merged_dirs has cached snippets, doing the
    work ourselves for whichever ones aren't locked by another process, and
    briefly waiting for the rest (some other process may already be on them --
    preprocessing a single location is comparatively quick). Returns False
    (without raising) if some locations are still missing once max_wait_s has
    elapsed, so the caller can defer this xval to a later pass instead of
    hanging indefinitely.
    """
    deadline = time.time() + max_wait_s
    failed: set = set()  # locations that errored this call -- don't tight-loop retrying them
    while True:
        missing = [d for d in merged_dirs if not location_done(d)]
        if not missing:
            return True

        did_work = False
        for merged_dir in missing:
            if merged_dir in failed or location_done(merged_dir):
                continue
            with batch_lock.try_lock(get_location_lock_path(merged_dir)) as acquired:
                if not acquired:
                    continue
                if location_done(merged_dir):
                    continue
                did_work = True
                print(f"    Preprocessing location: {merged_dir}")
                try:
                    preprocess_one_location(merged_dir)
                except Exception as e:
                    print(f"    ERROR preprocessing {merged_dir}: {e}")
                    traceback.print_exc()
                    failed.add(merged_dir)

        missing = [d for d in merged_dirs if not location_done(d)]
        if not missing:
            return True
        if all(d in failed for d in missing):
            print(f"    Giving up: {len(missing)} location(s) failed to preprocess: {missing}")
            return False
        if time.time() > deadline:
            print(
                f"    Timed out after {max_wait_s}s waiting for {len(missing)} "
                f"location(s) still being preprocessed elsewhere: {missing}"
            )
            return False
        if not did_work:
            time.sleep(poll_s)


# ── training dataset pooling (v1) ────────────────────────────────────────────


class MultiLocationFinetuneDatasetV1(npdataset.NeuropixelsDataset_cortexlab):
    """
    Pools sessions from multiple preprocessed locations into one v1 finetune
    dataset, so TrainExperimentBatchSampler builds one batch per session
    across every selected mouse's locations. Mirrors
    train_deepunitmatch_from_merged.MultiLocationFinetuneDataset, but against
    the v1 (non-_v2) npdataset.NeuropixelsDataset_cortexlab so the trained
    checkpoint stays compatible with DeepUnitMatch.testing.test.

    Deliberately doesn't call NeuropixelsDataset_cortexlab.__init__: that
    constructor is built around a single data_dir, so it's simpler to fill in
    the attributes it expects (mode, unit_order, experiment_unit_map,
    all_files) directly, reusing its inherited
    select_good_units_files/__getitem__/_augment_original for everything else.
    """

    def __init__(self, location_dirs: Sequence[str], mode: str = "train"):
        self.mode = mode
        self.unit_order = "filesystem"
        self.experiment_unit_map: Dict[int, List[str]] = {}
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
            f"  MultiLocationFinetuneDatasetV1: {len(self.all_files)} unit-files "
            f"across {len(self.experiment_unit_map)} session(s) from "
            f"{len(location_dirs)} location(s)."
        )


def finetune_ckpt_dir(exp_name: str) -> str:
    return str(DEEPUNITMATCH_DIR / "ModelExp" / "experiments" / exp_name / "ckpt")


def latest_finetune_checkpoint(exp_name: str) -> str:
    ckpt_dir = finetune_ckpt_dir(exp_name)
    ckpts = [f for f in os.listdir(ckpt_dir) if f.startswith("ckpt_epoch_")]
    if not ckpts:
        raise FileNotFoundError(f"No checkpoint found in {ckpt_dir}")
    ckpts.sort(key=lambda x: int(x.split("_")[-1]))
    return os.path.join(ckpt_dir, ckpts[-1])


def run_ae_stage(exp_name: str, locations: Sequence[str], args: argparse.Namespace) -> None:
    dataset = ConcatDataset(
        [AE_NeuropixelsDataset(loc, batch_size=args.ae_batchsize) for loc in locations]
    )
    train_ae_mod.run_training(
        exp_name=exp_name,
        dataset=dataset,
        lr=args.ae_lr,
        save_freq=args.ae_save_freq,
        total_epoch=args.ae_total_epoch,
        cont=True,  # safe no-op if already fully trained; resumes a partial local run
        batchsize=args.ae_batchsize,
        launch_tensorboard=False,
    )


def run_finetune_stage(exp_name: str, locations: Sequence[str], args: argparse.Namespace) -> None:
    dataset = MultiLocationFinetuneDatasetV1(locations, mode="train")
    # exp_name must match run_ae_stage's: run_finetune looks up the AE
    # checkpoint under ModelExp/AE_experiments/<exp_name>.
    train_finetune_mod.run_finetune(
        exp_name=exp_name,
        dataset=dataset,
        lr_enc=args.lr_enc,
        lr_proj=args.lr_proj,
        save_freq=args.finetune_save_freq,
        total_epoch=args.finetune_total_epoch,
        cont=True,
        batchsize=args.finetune_batchsize,
        launch_tensorboard=False,
    )


# ── inference (reuses run_deepunitmatch_batch_onMerged's pipeline) ──────────


def get_xval_save_dir(merged_dir: str, xval_name: str) -> str:
    subfolder = os.path.relpath(os.path.dirname(merged_dir), base_batch.BASE_INPUT)
    return os.path.join(base_batch.BASE_OUTPUT, subfolder, f"xval_{xval_name}")


def xval_group_done(merged_dir: str, xval_name: str) -> bool:
    sentinel = os.path.join(get_xval_save_dir(merged_dir, xval_name), "MatchingOverview.png")
    return batch_lock.sentinel_is_fresh(sentinel)


def run_inference_stage(xval_name: str, checkpoint_path: str, eval_groups: Sequence[str]) -> None:
    print(f"  Loading model for inference: {checkpoint_path}")
    model = test.load_trained_model(
        device=base_batch.DEVICE, read_path=checkpoint_path, n_output=N_OUTPUT
    )
    for merged_dir in eval_groups:
        if xval_group_done(merged_dir, xval_name):
            continue
        save_dir = get_xval_save_dir(merged_dir, xval_name)
        try:
            sess = base_batch._prepare_session(merged_dir)
            if sess is None:
                continue
            base_batch.run_deep_unit_match_core(
                sess, save_dir, model, label=f"xval_{xval_name}"
            )
        except Exception as e:
            print(f"    xval_{xval_name} inference FAILED on {merged_dir}: {e}")
            traceback.print_exc()


# ── one full xval pipeline ───────────────────────────────────────────────────


def run_one_xval(
    name: str,
    mice: List[str],
    args: argparse.Namespace,
    groups_by_mouse: "dict[str, List[str]]",
) -> None:
    exp_name = f"xval_{name}"
    lock_path = xval_lock_path(name)

    with batch_lock.try_lock(lock_path, stale_after=TRAINING_STALE_AFTER_SECONDS) as acquired:
        if not acquired:
            print(f"  Skipping {name} (already being processed by another run).")
            return

        # re-check now that we hold the lock: another process may have
        # finished this xval while we were scanning/waiting for the lock
        if xval_is_done(name):
            print(f"  Skipping {name} (completed by another run).")
            return

        write_status(name, state="preprocessing", mice=mice)
        train_locations = [
            d for mouse in mice for d in groups_by_mouse.get(mouse, [])
        ]
        if not train_locations:
            write_status(name, state="failed", error="no training locations found for selected mice")
            print(f"  {name}: no merged-data locations found for mice {mice}; skipping.")
            return

        xval_dir(name).mkdir(parents=True, exist_ok=True)
        (xval_dir(name) / "mouse_selection.json").write_text(
            json.dumps(
                {
                    "xval_name": name,
                    "mice": mice,
                    "n_mice": len(mice),
                    "excluded_mice": [m for m in CANONICAL_MICE if m not in mice],
                    "n_train_locations": len(train_locations),
                },
                indent=2,
                sort_keys=True,
            )
        )

        ok = ensure_locations_preprocessed(train_locations, max_wait_s=args.preprocess_wait_s)
        if not ok:
            write_status(name, state="queued")  # not done -- retry on a later pass
            print(f"  {name}: preprocessing incomplete for now; will retry on a later run.")
            return

        # AE_NeuropixelsDataset/MultiLocationFinetuneDatasetV1 read from each
        # location's *cached, preprocessed* folder (processed_waveforms/ under
        # PREPROCESSED_CACHE) -- not from the raw merged_dir itself, which has
        # no processed_waveforms/ of its own until preprocess_one_location()
        # (called above, keyed by merged_dir) has written one there.
        cache_locations = [get_location_out_dir(d) for d in train_locations]

        write_status(name, state="training_ae")
        print(f"  {name}: training AE on {len(cache_locations)} location(s) from {len(mice)} mice.")
        try:
            run_ae_stage(exp_name, cache_locations, args)
        except Exception as e:
            write_status(name, state="failed", error=f"AE training: {e}")
            print(f"  {name}: AE TRAINING FAILED: {e}")
            traceback.print_exc()
            return

        write_status(name, state="training_finetune")
        print(f"  {name}: finetuning on {len(cache_locations)} location(s).")
        try:
            run_finetune_stage(exp_name, cache_locations, args)
        except Exception as e:
            write_status(name, state="failed", error=f"finetune training: {e}")
            print(f"  {name}: FINETUNE TRAINING FAILED: {e}")
            traceback.print_exc()
            return

        try:
            checkpoint_path = latest_finetune_checkpoint(exp_name)
        except FileNotFoundError as e:
            write_status(name, state="failed", error=str(e))
            print(f"  {name}: {e}")
            return

        eval_groups = [d for group in groups_by_mouse.values() for d in group]
        write_status(name, state="inference", checkpoint=checkpoint_path)
        print(f"  {name}: running inference on {len(eval_groups)} merged-data group(s).")
        run_inference_stage(name, checkpoint_path, eval_groups)

        write_status(name, state="done", checkpoint=checkpoint_path)
        print(f"  {name}: DONE.")


# ── entry point ───────────────────────────────────────────────────────────────


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Train + evaluate DeepUnitMatch on 3/6/12/18-mouse subsets, end to end."
    )
    parser.add_argument(
        "--only", nargs="+", default=None,
        help="Only run these xval names (default: all 13, in a fixed order). "
        "Names are m1_1..m1_3, m3_1..m3_3, m6_1..m6_3, m12_1..m12_3, m18_1.",
    )
    parser.add_argument(
        "--status", action="store_true",
        help="Print every xval's current state (from its status.json) and exit.",
    )
    parser.add_argument("--ae_total_epoch", type=int, default=300)
    parser.add_argument("--finetune_total_epoch", type=int, default=50)
    parser.add_argument("--ae_batchsize", type=int, default=32)
    parser.add_argument("--finetune_batchsize", type=int, default=40)
    parser.add_argument("--ae_lr", type=float, default=1e-5)
    parser.add_argument("--lr_enc", type=float, default=2e-5)
    parser.add_argument("--lr_proj", type=float, default=1.1e-4)
    parser.add_argument("--ae_save_freq", type=int, default=1)
    parser.add_argument("--finetune_save_freq", type=int, default=1)
    parser.add_argument(
        "--preprocess_wait_s", type=int, default=1800,
        help="How long to wait for locations another process is currently "
        "preprocessing before deferring this xval to a later pass.",
    )
    parser.add_argument(
        "--write-matlab-compat", action="store_true",
        help="Also write a MATLAB-compatible UnitMatch.mat from the Python inference outputs.",
    )
    return parser


def main() -> None:
    args = build_arg_parser().parse_args()
    base_batch.WRITE_MATLAB_COMPAT = args.write_matlab_compat

    manifest = generate_manifest()

    if args.status:
        print_status_table(manifest)
        return

    print(f"Scanning for merged-data groups under:\n  {base_batch.BASE_INPUT}\n")
    groups = base_batch.find_merged_groups()
    if not groups:
        print("No merged-data groups found.")
        return

    groups_by_mouse_all = group_by_mouse(groups)
    groups_by_mouse = {
        m: groups_by_mouse_all[m] for m in CANONICAL_MICE if m in groups_by_mouse_all
    }
    missing_mice = [m for m in CANONICAL_MICE if m not in groups_by_mouse_all]
    if missing_mice:
        print(f"WARNING: no merged-data groups found for: {missing_mice}")
    n_groups = sum(len(v) for v in groups_by_mouse.values())
    print(
        f"Found {n_groups} merged-data group(s) across "
        f"{len(groups_by_mouse)}/{len(CANONICAL_MICE)} canonical mice.\n"
    )

    names = args.only if args.only else list(manifest.keys())
    unknown = [n for n in names if n not in manifest]
    if unknown:
        raise ValueError(f"Unknown xval name(s): {unknown}. Valid names: {list(manifest.keys())}")

    for name in names:
        if xval_is_done(name):
            print(f"[{name}] already done, skipping.")
            continue
        print(f"\n[{name}] mice={manifest[name]}")
        run_one_xval(name, manifest[name], args, groups_by_mouse)

    print("\nAll done (for this pass).")


if __name__ == "__main__":
    main()
