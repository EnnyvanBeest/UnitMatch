# Train and evaluate DeepUnitMatch models on the merged data (pipeline step 3).
#
# One driver for every model the paper needs, so all of them are trained on
# the same data with the same code. Job lists:
#
#   comparison  Which way of training to use for the paper models. The
#               three 3-mouse subsets of the cross-validation (m3_1..m3_3,
#               xval_end_to_end.generate_manifest), each trained as
#                 A_original  training code before the fixes (frozen snapshot,
#                             DeepUnitMatch/ExtraModels/legacy_v1)
#                 B_fixed     current code (fixed v1: frozen backbone)
#                 C_v2        current code + v2 options (trainable backbone,
#                             channel positions, amplitude/noise jitter,
#                             symmetric loss)
#               and evaluated on the held-out mice. The three variants of a
#               subset share one autoencoder. Summarise with
#               compare_training_variants.py.
#   paper       every model in the paper (main model, W_ij and n_output
#               sweeps, baselines, mouse-count cross-validation), all with
#               the winning variant -- set PAPER_VARIANT first (see paper_jobs).
#
# Stages (run in this order; each can run on several machines at once --
# every unit of work is claimed with a lock on the share):
#   preprocess  DNN snippets of the training locations -> cfg.TRAINING_CACHE
#   train       autoencoders + fine-tuning; every finished model is copied to
#               cfg.MODELS_ROOT/<name>/ (checkpoints are first written to the
#               training machine's local DeepUnitMatch/ModelExp)
#   evaluate    each trained model on its evaluation locations, through the
#               matching pipeline shared with UMPy (run_dum_core) ->
#               cfg.ANALYSIS_OUTPUT/<mouse>/<probe>/<location>/<model name>/
#   status      table of every job's progress
#
# Examples:
#   python train_paper_models.py --jobs comparison --stage preprocess
#   python train_paper_models.py --jobs comparison --stage train      # GPU machine(s)
#   python train_paper_models.py --jobs comparison --stage evaluate   # any machine(s)
#   python train_paper_models.py --jobs comparison --stage status
#   python compare_training_variants.py

import argparse
import datetime
import json
import os
import shutil
import socket
import sys
import traceback

_HERE = os.path.dirname(os.path.abspath(__file__))
_REPO = os.path.dirname(_HERE)  # .../UnitMatchPy
_DUM = os.path.join(_REPO, "DeepUnitMatch")
for _p in (os.path.join(_DUM, "ExtraModels"), _DUM, _REPO, _HERE):
    if _p not in sys.path:
        sys.path.insert(0, _p)

import numpy as np
import torch
from torch.utils.data import ConcatDataset

import batch_lock
import pipeline_config as cfg
import pipeline_log as plog
import run_deepunitmatch_batch_onMerged as base_batch
from DeepUnitMatch.testing import test
from DeepUnitMatch.utils import param_fun
from DeepUnitMatch.utils.AE_npdataset import AE_NeuropixelsDataset
from train import train_AE as train_ae_mod
from train import train_finetune as finetune_fixed
from utils import npdataset as npdataset_fixed
from legacy_v1 import train_finetune as finetune_legacy
from legacy_v1 import npdataset as npdataset_legacy
from xval_end_to_end import generate_manifest

LOG_STAGE = "training"
MODELEXP = os.path.join(_DUM, "ModelExp")  # where train_AE / train_finetune write
TRAINING_STALE_AFTER_SECONDS = 5 * 24 * 3600  # a training job can run for days

# Training settings shared by every job (as in the Methods)
AE_EPOCHS = 300
AE_LR = 1e-5
AE_BATCHSIZE = 32
FT_EPOCHS = 50
FT_BATCHSIZE = 40

VARIANTS = {
    "A_original": {"code": "legacy"},
    "B_fixed": {"code": "fixed"},
    "C_v2": {
        "code": "fixed",
        "lr_backbone": 2e-6,
        "with_geometry": True,
        "amp_noise_jitter": True,
        "symmetric_loss": True,
    },
}


# ── jobs ─────────────────────────────────────────────────────────────────────


def make_job(name, ae, train_mice, variant, evaluate="heldout", ft_options=None,
             kind="finetune", ae_n_output=256, same_as=None):
    """
    name        model name (= output folder name for its evaluation)
    ae          autoencoder experiment name (shared between jobs that use it);
                None for a fine-tuning on a random backbone
    train_mice  mice whose locations are used for training (AE and fine-tuning)
    variant     key of VARIANTS
    evaluate    "heldout" (all mice not trained on), "all", or a list of mice
    ft_options  extra run_finetune arguments (n_output, negative_weight, random_backbone)
    kind        "finetune" (AE + fine-tuning), "ae_only" (the AE itself is the
                model) or "untrained" (randomly initialised encoder)
    ae_n_output encoder size of the AE (jobs sharing an AE must agree)
    same_as     name of a job trained identically: its model is reused
                (copied) instead of training again
    """
    return {
        "name": name,
        "ae": ae,
        "train_mice": sorted(train_mice),
        "variant": variant,
        "evaluate": evaluate,
        "ft_options": ft_options or {},
        "kind": kind,
        "ae_n_output": ae_n_output,
        "same_as": same_as,
    }


def comparison_jobs(manifest):
    jobs = []
    for rep in (1, 2, 3):
        subset = manifest[f"m3_{rep}"]
        ae = f"cmp_m3_{rep}"
        for variant in VARIANTS:
            jobs.append(make_job(f"cmp_m3_{rep}_{variant}", ae, subset, variant))
    return jobs


# Training variant of every paper model: set to the winner of the comparison
# (compare_training_variants.py) before training any paper job.
PAPER_VARIANT = None


def paper_jobs(manifest):
    """
    Every model in the paper, all trained with PAPER_VARIANT:
      DUM_paper            main model: all mice (W_ij = 10, n_output = 256),
                           evaluated on all mice
      DUM_W_ij=<w>         W_ij 1/5/15/20; same AE as DUM_paper
      DUM_n_output=<n>     encoder size 8/32/128, each with its own AE
      DUM_untrained        randomly initialised encoder (no training)
      DUM_AE_only          the autoencoder of DUM_paper, not fine-tuned
      DUM_finetuned_only   fine-tuning without AE pretraining (random backbone)
      xval_m<k>_<r>        trained on k = 1/3/6/12 mice (3 subsets each, from
                           xval_end_to_end.generate_manifest). The 3-mouse
                           models are the comparison models of PAPER_VARIANT
                           (same mice, same training), reused rather than
                           retrained; the 18-mouse point is DUM_paper.
    Every paper model is evaluated (inference) on all mice; which mice a
    model was trained on is in MODELS_ROOT/<name>/model.json and
    TRAINING_STATE_ROOT/manifest.json, so figures can mark training vs.
    held-out points.
    """
    if PAPER_VARIANT is None:
        raise SystemExit("Set PAPER_VARIANT (top of train_paper_models.py) to the training variant "
                         "chosen with compare_training_variants.py first.")
    v = PAPER_VARIANT
    everyone = manifest["m18_1"]
    jobs = [make_job("DUM_paper", "paper_ae", everyone, v, evaluate="all")]
    for w in (1, 5, 15, 20):
        jobs.append(make_job(f"DUM_W_ij={w}", "paper_ae", everyone, v, evaluate="all",
                             ft_options={"negative_weight": float(w)}))
    for n in (8, 32, 128):
        jobs.append(make_job(f"DUM_n_output={n}", f"paper_ae_n{n}", everyone, v, evaluate="all",
                             ft_options={"n_output": n}, ae_n_output=n))
    jobs.append(make_job("DUM_untrained", "paper_ae_untrained", everyone, v, evaluate="all", kind="untrained"))
    jobs.append(make_job("DUM_AE_only", "paper_ae", everyone, v, evaluate="all", kind="ae_only"))
    jobs.append(make_job("DUM_finetuned_only", None, everyone, v, evaluate="all",
                         ft_options={"random_backbone": True}))
    for k in (1, 3, 6, 12):
        for rep in (1, 2, 3):
            key = f"m{k}_{rep}"
            same_as = f"cmp_m3_{rep}_{v}" if k == 3 else None
            jobs.append(make_job(f"xval_{key}", f"xval_{key}_ae", manifest[key], v, evaluate="all",
                                 same_as=same_as))
    return jobs


JOB_LISTS = {"comparison": comparison_jobs, "paper": paper_jobs}


# ── locations / paths ────────────────────────────────────────────────────────


def merged_groups():
    """{'mouse/probe/location': merged_dir} of every merged-data location."""
    return {base_batch.group_key(d): d for d in base_batch.find_merged_groups()}


def mouse_of(group):
    return group.split("/")[0]


def cache_dir(group):
    return os.path.join(cfg.TRAINING_CACHE, *group.split("/"))


def state_dir(name):
    return os.path.join(cfg.TRAINING_STATE_ROOT, name)


def published_dir(name):
    return os.path.join(cfg.MODELS_ROOT, name)


def read_status(name):
    try:
        with open(os.path.join(state_dir(name), "status.json")) as f:
            return json.load(f)
    except (OSError, json.JSONDecodeError):
        return {}


def write_status(name, **fields):
    status = read_status(name)
    status.update(fields)
    os.makedirs(state_dir(name), exist_ok=True)
    with open(os.path.join(state_dir(name), "status.json"), "w") as f:
        json.dump(status, f, indent=1)


def lock_path(kind, name):
    return os.path.join(cfg.TRAINING_STATE_ROOT, ".locks", f"{kind}__{name.replace('/', '__')}.lock")


def publish(name, ckpt_path, extra=None):
    """Copy a finished checkpoint (+ metadata) to MODELS_ROOT/<name>/."""
    out = published_dir(name)
    os.makedirs(out, exist_ok=True)
    shutil.copy2(ckpt_path, os.path.join(out, "model.pt"))
    meta = {"source": ckpt_path, "machine": socket.gethostname(),
            "time": datetime.datetime.now().isoformat(timespec="seconds")}
    meta.update(extra or {})
    with open(os.path.join(out, "model.json"), "w") as f:
        json.dump(meta, f, indent=1)
    return os.path.join(out, "model.pt")


# ── stage 1: preprocessing ───────────────────────────────────────────────────


def preprocess(groups):
    """DNN snippets of every location (param_fun.get_snippets) into cfg.TRAINING_CACHE."""
    for i, (group, merged_dir) in enumerate(sorted(groups.items())):
        sentinel = os.path.join(cache_dir(group), ".done")
        if os.path.isfile(sentinel):
            continue
        with batch_lock.try_lock(lock_path("preprocess", group)) as acquired:
            if not acquired or os.path.isfile(sentinel):
                continue
            print(f"[{i + 1}/{len(groups)}] preprocessing {group}")
            try:
                sess = base_batch._prepare_session(merged_dir)
                if sess is None:
                    raise RuntimeError("session preparation failed (see pipeline log)")
                out = cache_dir(group)
                if os.path.isdir(out):
                    shutil.rmtree(out)  # unfinished earlier attempt
                os.makedirs(out)
                param_fun.get_snippets(
                    sess["waveform"], sess["channel_pos"], sess["session_id"], save_path=out,
                    unit_ids=np.concatenate(sess["param"]["good_units"]).squeeze(),
                    param=sess["param"],
                )
                with open(sentinel, "w") as f:
                    f.write("ok")
                plog.log_event(LOG_STAGE, group, "preprocess", "done")
            except Exception as e:
                traceback.print_exc()
                plog.log_event(LOG_STAGE, group, "preprocess", "failed", f"{type(e).__name__}: {e}",
                               traceback.format_exc())


def training_locations(mice, groups):
    """Preprocessed locations of these mice (missing ones are reported, not silently dropped)."""
    locations, missing = [], []
    for group in sorted(groups):
        if mouse_of(group) in mice:
            (locations if os.path.isfile(os.path.join(cache_dir(group), ".done")) else missing).append(group)
    if missing:
        raise RuntimeError(f"not preprocessed yet: {missing} -- run --stage preprocess first")
    return [cache_dir(g) for g in locations]


def finetune_dataset(base_cls, location_dirs):
    """All sessions of several locations in one training dataset (one batch never mixes sessions)."""

    class MultiLocationDataset(base_cls):
        def __init__(self):  # deliberately not base_cls.__init__ (single data_dir)
            self.mode = "train"
            self.unit_order = "filesystem"
            self.experiment_unit_map = {}
            next_id = 0
            for location_dir in location_dirs:
                data_dir = os.path.join(location_dir, "processed_waveforms")
                for session in sorted(os.listdir(data_dir), key=lambda s: int(s) if s.isdigit() else s):
                    files = self.select_good_units_files(os.path.join(data_dir, session), load_pre_merge=False)
                    if files:
                        self.experiment_unit_map[next_id] = files
                        next_id += 1
            self.all_files = [(e, f) for e, files in self.experiment_unit_map.items() for f in files]
            print(f"training data: {len(self.all_files)} units in {len(self.experiment_unit_map)} sessions "
                  f"from {len(location_dirs)} locations")

    return MultiLocationDataset()


# ── stage 2: training ────────────────────────────────────────────────────────


def ensure_local_ae(ae):
    """Fine-tuning reads the AE from the local ModelExp; fetch it from the share if trained elsewhere."""
    local_ckpt_dir = os.path.join(MODELEXP, "AE_experiments", ae, "ckpt")
    if os.path.isdir(local_ckpt_dir) and any(f.startswith("ckpt_epoch_") for f in os.listdir(local_ckpt_dir)):
        return
    shared = os.path.join(published_dir(ae), "model.pt")
    if not os.path.isfile(shared):
        raise RuntimeError(f"autoencoder {ae} not trained yet")
    os.makedirs(local_ckpt_dir, exist_ok=True)
    epoch = read_status(ae).get("epochs", AE_EPOCHS) - 1
    shutil.copy2(shared, os.path.join(local_ckpt_dir, f"ckpt_epoch_{epoch}"))


def train_autoencoder(ae, mice, groups, n_output=256, epochs=AE_EPOCHS):
    """Train (or, with epochs=0, just initialise) an autoencoder; True once it is available."""
    if read_status(ae).get("trained"):
        return True
    with batch_lock.try_lock(lock_path("train", ae), stale_after=TRAINING_STALE_AFTER_SECONDS) as acquired:
        if not acquired:
            print(f"  autoencoder {ae} is being trained elsewhere")
            return False
        if read_status(ae).get("trained"):
            return True
        print(f"=== training autoencoder {ae} on {mice} (n_output {n_output}, {epochs} epochs)")
        if epochs == 0:
            torch.manual_seed(0)  # reproducible random initialisation (untrained baseline)
        dataset = ConcatDataset([AE_NeuropixelsDataset(d, batch_size=AE_BATCHSIZE)
                                 for d in training_locations(mice, groups)])
        train_ae_mod.run_training(exp_name=ae, dataset=dataset, lr=AE_LR, total_epoch=epochs,
                                  cont=True, batchsize=AE_BATCHSIZE, launch_tensorboard=False,
                                  n_output=n_output)
        ckpt = test.latest_checkpoint(os.path.join(MODELEXP, "AE_experiments", ae, "ckpt"))
        publish(ae, ckpt, {"kind": "autoencoder", "train_mice": sorted(mice), "n_output": n_output,
                           "epochs": epochs})
        write_status(ae, trained=datetime.datetime.now().isoformat(timespec="seconds"),
                     machine=socket.gethostname(), epochs=max(epochs, 1), train_mice=sorted(mice),
                     n_output=n_output)
        plog.log_event(LOG_STAGE, ae, "autoencoder", "done")
        return True


def _mark_trained(name, job, ckpt_path, kind):
    publish(name, ckpt_path, {"kind": kind, "job": job})
    write_status(name, trained=datetime.datetime.now().isoformat(timespec="seconds"),
                 machine=socket.gethostname(), job=job)
    plog.log_event(LOG_STAGE, name, "train", "done", kind if kind != "finetuned" else "")


def train_job(job, groups):
    name = job["name"]
    if read_status(name).get("trained"):
        return
    kind = job.get("kind", "finetune")

    if job.get("same_as"):
        # trained identically before (e.g. a comparison model): reuse it
        source = job["same_as"]
        if not read_status(source).get("trained"):
            print(f"  {name}: waiting for {source} to be trained")
            return
        _mark_trained(name, job, os.path.join(published_dir(source), "model.pt"), f"same as {source}")
        return

    if kind in ("ae_only", "untrained"):
        epochs = 0 if kind == "untrained" else AE_EPOCHS
        if train_autoencoder(job["ae"], job["train_mice"], groups, job["ae_n_output"], epochs):
            _mark_trained(name, job, os.path.join(published_dir(job["ae"]), "model.pt"), kind)
        return

    if job["ae"] is not None and not train_autoencoder(job["ae"], job["train_mice"], groups, job["ae_n_output"]):
        return  # AE still training elsewhere: try this job later
    with batch_lock.try_lock(lock_path("train", name), stale_after=TRAINING_STALE_AFTER_SECONDS) as acquired:
        if not acquired or read_status(name).get("trained"):
            return
        print(f"=== fine-tuning {name} ({job['variant']})")
        try:
            if job["ae"] is not None:
                ensure_local_ae(job["ae"])
            spec = dict(VARIANTS[job["variant"]])
            code = spec.pop("code")
            locations = training_locations(job["train_mice"], groups)
            if code == "legacy" and (job["ft_options"] or job["ae"] is None):
                raise RuntimeError("the original training code has no options (n_output, W_ij, random "
                                   "backbone): sweeps and baselines need a variant with the fixed code")
            if code == "legacy":
                # the original code reads the AE from an experiment of the same
                # name, so its fine-tuning experiment is named after the AE
                exp = job["ae"]
                finetune_legacy.run_finetune(
                    exp, finetune_dataset(npdataset_legacy.NeuropixelsDataset_cortexlab, locations),
                    total_epoch=FT_EPOCHS, cont=True, batchsize=FT_BATCHSIZE, launch_tensorboard=False,
                )
            else:
                exp = name
                spec.update(job["ft_options"])
                finetune_fixed.run_finetune(
                    exp, finetune_dataset(npdataset_fixed.NeuropixelsDataset_cortexlab, locations),
                    ae_exp_name=job["ae"], total_epoch=FT_EPOCHS, cont=True, batchsize=FT_BATCHSIZE,
                    launch_tensorboard=False, **spec,
                )
            ckpt = test.latest_checkpoint(os.path.join(MODELEXP, "experiments", exp, "ckpt"))
            _mark_trained(name, job, ckpt, "finetuned")
        except Exception as e:
            traceback.print_exc()
            plog.log_event(LOG_STAGE, name, "train", "failed", f"{type(e).__name__}: {e}", traceback.format_exc())


# ── stage 3: evaluation ──────────────────────────────────────────────────────


def evaluation_groups(job, groups):
    all_mice = set(mouse_of(g) for g in groups)
    if job["evaluate"] == "heldout":
        mice = all_mice - set(job["train_mice"])
    elif job["evaluate"] == "all":
        mice = all_mice
    else:
        mice = set(job["evaluate"])
    return sorted(g for g in groups if mouse_of(g) in mice)


def eval_dir(group, name):
    return os.path.join(cfg.ANALYSIS_OUTPUT, *group.split("/"), name)


def eval_done(group, name):
    return os.path.isfile(os.path.join(eval_dir(group, name), base_batch.SENTINEL))


def evaluate(jobs, groups):
    trained = [j for j in jobs if read_status(j["name"]).get("trained")]
    if len(trained) < len(jobs):
        print(f"{len(jobs) - len(trained)} job(s) not trained yet; evaluating the {len(trained)} trained ones")
    todo = {}
    for job in trained:
        for group in evaluation_groups(job, groups):
            if not eval_done(group, job["name"]):
                todo.setdefault(group, []).append(job)
    models = {}
    for i, (group, group_jobs) in enumerate(sorted(todo.items())):
        with batch_lock.try_lock(lock_path("evaluate", group)) as acquired:
            if not acquired:
                continue
            group_jobs = [j for j in group_jobs if not eval_done(group, j["name"])]
            if not group_jobs:
                continue
            print(f"[{i + 1}/{len(todo)}] {group}: {[j['name'] for j in group_jobs]}")
            sess = base_batch._prepare_session(groups[group])
            if sess is None:
                continue  # logged by _prepare_session
            for job in group_jobs:
                name = job["name"]
                if name not in models:
                    models[name] = test.load_trained_model(
                        device=base_batch.DEVICE, read_path=os.path.join(published_dir(name), "model.pt"))
                try:
                    base_batch.run_dum_core(sess, eval_dir(group, name), models[name], label=name)
                except Exception as e:
                    print(f"  {name} FAILED: {e}")  # logged by run_dum_core
                    traceback.print_exc()


# ── status ───────────────────────────────────────────────────────────────────


def print_status(jobs, groups):
    print(f"{'job':28s} {'variant':11s} {'AE':4s} {'trained':20s} evaluated")
    for job in jobs:
        st = read_status(job["name"])
        eval_groups = evaluation_groups(job, groups)
        n_done = sum(eval_done(g, job["name"]) for g in eval_groups)
        ae = "-" if job["ae"] is None else ("yes" if read_status(job["ae"]).get("trained") else "no")
        print(f"{job['name']:28s} {job['variant']:11s} {ae:4s} "
              f"{st.get('trained', '-'):20s} {n_done}/{len(eval_groups)}")


# ── entry point ──────────────────────────────────────────────────────────────


def main():
    parser = argparse.ArgumentParser(description="Train and evaluate DeepUnitMatch models on the merged data.")
    parser.add_argument("--jobs", choices=sorted(JOB_LISTS), required=True)
    parser.add_argument("--stage", choices=["preprocess", "train", "evaluate", "status"], required=True)
    parser.add_argument("--only", nargs="+", help="restrict to these job names")
    args = parser.parse_args()

    groups = merged_groups()
    manifest = generate_manifest()
    missing_mice = sorted({m for mice in manifest.values() for m in mice} - {mouse_of(g) for g in groups})
    if missing_mice:
        print(f"WARNING: mice of the subset manifest without merged data: {missing_mice}")
    os.makedirs(cfg.TRAINING_STATE_ROOT, exist_ok=True)
    with open(os.path.join(cfg.TRAINING_STATE_ROOT, "manifest.json"), "w") as f:
        json.dump(manifest, f, indent=1)

    jobs = JOB_LISTS[args.jobs](manifest)
    if args.only:
        jobs = [j for j in jobs if j["name"] in set(args.only)]

    if args.stage == "preprocess":
        # only training locations: evaluation makes its own snippets (run_dum_core)
        training_mice = {m for j in jobs for m in j["train_mice"]}
        preprocess({g: d for g, d in groups.items() if mouse_of(g) in training_mice})
    elif args.stage == "train":
        for job in jobs:
            train_job(job, groups)
    elif args.stage == "evaluate":
        evaluate(jobs, groups)
    print_status(jobs, groups)


if __name__ == "__main__":
    main()
