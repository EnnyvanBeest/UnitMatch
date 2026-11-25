# Tables describing the dataset of the paper: what was recorded, where, when,
# for how long, and during which experiments.
#
# Reads only (network shares); writes four CSV files to --out
# (default RESULTS_DIR/dataset_table):
#   dataset_sessions.csv   one row per recording session (the unit of matching)
#   dataset_locations.csv  one row per recording location (mouse/probe/loc)
#   dataset_mice.csv       one row per mouse; sex/strain/date of birth/implant
#                          date are left empty to fill in by hand (e.g. from Alyx)
#   expdef_names.csv       every expDef found, with its count and category, to
#                          check/extend EXPDEF_CATEGORIES below
#
# Sources:
#   - sessions of each location and the units per session: the DeepUnitMatch
#     output on the merged data (ANALYSIS_OUTPUT/<loc>/DeepUnitMatch/UMparam.pickle);
#     only the sessions the pipeline kept (KS_dirs); those in KS_dirs_all but not
#     in KS_dirs were dropped (no units left) and are left out
#   - used_in_analysis: False for locations with a single session (no session
#     pair, so in none of the across-session analyses)
#   - used_in_model: the paper model (DUM_paper) was trained on all sessions of
#     its training mice (model.json), except for the 5% of sessions held out to
#     validate the fine-tuning (session_split.json, --split; these are False).
#     The autoencoder step used a random 95% of all units, also of those sessions.
#     The split numbers sessions in the order fine-tuning loaded them from the
#     training cache, which is rebuilt here; left empty if that fails.
#   - original recording of each session: metadata_index.json (merged session
#     folder i = exp_key i + 1)
#   - duration, start time, probe type/serial, sampling rate: the SpikeGLX
#     .ap.meta of the recording (found via dat_path in the Kilosort params.py)
#   - recorded shanks and depth range: channel_positions.npy of the merged session
#   - expDefs: the Rigbox experiments (Subjects/<mouse>/<date>/<n>/) that
#     overlap the recording by >= 30 s (Block start/end, else Timeline start +
#     lastTimestamp). For the recordings of a day that overlap no experiment, a
#     whole-hour offset of the ephys clock is tried (some ephys computers were
#     set to another time zone) and reported in expdef_match. Recordings that
#     still overlap nothing get expdef_match = "none during recording" and an
#     empty expdefs; that day's experiments are in expdefs_same_day. The
#     SpikeGLX run name (recording_name) often says what was shown as well.
#
# Results of the slow network reads are cached in --out/_cache.json, so a rerun
# (e.g. after editing EXPDEF_CATEGORIES) is fast; delete it to read everything again.
#
# Usage (UMPy env):  python build_dataset_table.py [--out DIR] [--workers 16]

import argparse
import concurrent.futures as cf
import datetime as dt
import json
import os
import pickle
import re
import sys
import traceback

import numpy as np
import pandas as pd
import scipy.io as sio

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _HERE)
import pipeline_config as cfg

# expDef file name (without .m) -> category in the paper's table. Names not
# listed here come out as "UNMAPPED"; expdef_names.csv lists all of them.
EXPDEF_CATEGORIES = {
    "spontaneousActivity": "spontaneous activity",
    "JF_choiceworldStimuli_wheel_left_center_all": "visual learning task",
    "imageWorld2": "visual stimulation: natural images",
    "imageWorld_AllInOne": "visual stimulation: natural images",
    "imageWorld_lefthemifield": "visual stimulation: natural images",
    "sparseNoise": "visual stimulation: sparse noise",
    "AP_sparseNoise": "visual stimulation: sparse noise",
    "sparseNoiseAsync_NS2": "visual stimulation: sparse noise",
    "fullFieldGratingWorld": "visual stimulation: gratings",
    "AVPassive_ckeckerboard_postactive": "audiovisual stimulation (passive)",
    "AVPassive_ckeckerboard_postactive_forCelian": "audiovisual stimulation (passive)",
    "gratingAndSoundWorld": "audiovisual stimulation (passive)",
    "filmWorldRepeatDelayZero": "visual stimulation: natural movies",
    "playTwoVid": "visual stimulation: natural movies",
    "multiSpaceWorld_checker_training": "audiovisual task",
    "noGo_stage1": "visual task (go/no-go)",
    "noGo_stage4_q": "visual task (go/no-go)",
    # mpep experiments (no expDef): the stimulus x-file from Protocol.mat
    "stimSparseNoiseAsync_Luminance": "visual stimulation: sparse noise",
    "stimMergePureTonesAndFlickeringGratingsPlaids": "audiovisual stimulation (passive)",
    "stimCosAmplModulatedSineWavesOutput": "auditory stimulation",
    # Timeline experiments with a <date>_<n>_<mouse>_VRBehavior_trial.mat log
    "(VR behaviour log)": "virtual reality",
    "(no expDef, Timeline only)": "unknown (Timeline only)",
}

# Keyword in the SpikeGLX recording name (case-insensitive) -> category. Used
# only where no experiment during the recording says what was shown (only
# Timeline-only / unmapped experiments, or none); category_source is then
# "recording name". Any informative expDef wins over the name.
RECORDING_NAME_CATEGORIES = {
    "natim": "visual stimulation: natural images",
    "naturalimages": "visual stimulation: natural images",
    "audiovis": "audiovisual stimulation (passive)",
    "avpassive": "audiovisual stimulation (passive)",
    "sparsenoise": "visual stimulation: sparse noise",
    "gratings": "visual stimulation: gratings",
    "spont": "spontaneous activity",  # also "spontaneous"
}
UNINFORMATIVE_CATEGORIES = {"unknown (Timeline only)", "UNMAPPED"}

# SpikeGLX imDatPrb_type -> probe name
PROBE_TYPES = {
    "0": "Neuropixels 1.0",
    "21": "Neuropixels 2.0 single-shank",
    "2003": "Neuropixels 2.0 single-shank",
    "2004": "Neuropixels 2.0 single-shank",
    "24": "Neuropixels 2.0 four-shank",
    "2013": "Neuropixels 2.0 four-shank",
    "2014": "Neuropixels 2.0 four-shank",
}

SHANK_PITCH_UM = 250  # Neuropixels 2.0 four-shank
MIN_OVERLAP = dt.timedelta(seconds=30)  # experiment counted if it overlaps the recording this much
WINDOW_SLACK = dt.timedelta(minutes=2)  # experiments without end time: may start just before the recording
# ephys clock offsets tried for unmatched recordings: only the ones found in
# this dataset, +7/+8 h (CB020, FT033: an ephys computer on US Pacific time,
# +7 h after the US switch to winter time). Trying every whole hour gave chance
# matches onto other experiments of the day (-1 h AL032, +4 h EB014).
CLOCK_OFFSETS_H = [7, 8]
EXP_CACHE = "experiments_v5"  # cache key of read_experiments (v2: end times; v3: mpep x-files; v4: VR logs; v5: expDef from Block)


# ── cache ────────────────────────────────────────────────────────────────────


class Cache:
    def __init__(self, path):
        self.path = path
        self.data = {}
        if os.path.isfile(path):
            with open(path, "r", encoding="utf-8") as f:
                self.data = json.load(f)

    def get(self, kind, key, compute):
        store = self.data.setdefault(kind, {})
        if key not in store:
            store[key] = compute(key)
        return store[key]

    def save(self):
        tmp = self.path + ".tmp"
        with open(tmp, "w", encoding="utf-8") as f:
            json.dump(self.data, f)
        os.replace(tmp, self.path)


# ── SpikeGLX meta ────────────────────────────────────────────────────────────


def find_meta(ks_dir):
    """The .ap.meta of the recording sorted in ks_dir (dat_path in params.py), or None."""
    params = os.path.join(ks_dir, "params.py")
    if not os.path.isfile(params):
        return None
    with open(params, "r", errors="replace") as f:
        found = re.findall(r"dat_path\s*=\s*r?['\"](.*?)['\"]", f.read())
    if not found:
        return None
    # dat_path may be relative ("../\\server..."), on a drive letter or a list; only the name is used
    name = os.path.basename(found[0].replace("\\\\", "\\").replace("\\", "/"))
    name = re.sub(r"\.c?bin$", ".meta", name)
    cur = ks_dir
    for _ in range(6):  # output -> pyKS -> imec folder -> run folder -> ...
        cur = os.path.dirname(cur)
        candidate = os.path.join(cur, name)
        if os.path.isfile(candidate):
            return candidate
    return None


def read_meta(ks_dir):
    path = find_meta(ks_dir)
    if path is None:
        return {"meta_path": None}
    meta = {}
    with open(path, "r", errors="replace") as f:
        for line in f:
            if "=" in line:
                k, v = line.rstrip("\n").split("=", 1)
                meta[k.lstrip("~")] = v
    return {
        "meta_path": path,
        "start": meta.get("fileCreateTime"),
        "duration_s": float(meta["fileTimeSecs"]) if "fileTimeSecs" in meta else None,
        "probe_type_code": meta.get("imDatPrb_type"),
        "probe_sn": meta.get("imDatPrb_sn"),
        "sampling_rate": float(meta["imSampRate"]) if "imSampRate" in meta else None,
    }


# ── Rigbox experiments ───────────────────────────────────────────────────────


def matlab_datenum(value):
    value = float(np.asarray(value).ravel()[0])
    return (dt.datetime.fromordinal(int(value)) + dt.timedelta(days=value % 1)
            - dt.timedelta(days=366))


def mat_times(path, var):
    """(start, end) of var (block / Timeline) in a .mat file, v5 or v7.3; end may be None.

    Block: startDateTime / endDateTime. Timeline: startDateTime + lastTimestamp (s).
    """
    fields = ("startDateTime", "endDateTime", "lastTimestamp")
    try:
        s = sio.loadmat(path, squeeze_me=True, struct_as_record=False, variable_names=[var])[var]
        values = {k: getattr(s, k, None) for k in fields}
    except NotImplementedError:  # v7.3 (HDF5)
        import h5py
        with h5py.File(path, "r") as f:
            values = {k: (f[var][k][()] if k in f[var] else None) for k in fields}
    ok = {k: v for k, v in values.items() if v is not None and np.asarray(v).size > 0}
    if "startDateTime" not in ok:
        return None, None
    start = matlab_datenum(ok["startDateTime"])
    end = None
    if "endDateTime" in ok:
        end = matlab_datenum(ok["endDateTime"])
    elif "lastTimestamp" in ok:
        end = start + dt.timedelta(seconds=float(np.asarray(ok["lastTimestamp"]).ravel()[0]))
    return start, end


def expdef_name(exp_dir, files):
    """Base name of the experiment's defFunction (expDef), or None."""
    for suffix in ("_parameters.json", "_parameters.mat"):
        f = [x for x in files if x.endswith(suffix)]
        if not f:
            continue
        try:
            if suffix.endswith(".json"):
                with open(os.path.join(exp_dir, f[0]), "r", errors="replace") as fh:
                    definition = json.load(fh).get("defFunction")
            else:
                par = sio.loadmat(os.path.join(exp_dir, f[0]), squeeze_me=True,
                                  struct_as_record=False)["parameters"]
                definition = getattr(par, "defFunction", None)
        except Exception:
            continue
        if definition:
            return os.path.splitext(os.path.basename(str(definition).replace("\\", "/")))[0]
    # no (readable) parameters file: the Block stores the expDef as well
    for f in (x for x in files if x.endswith("_Block.mat")):
        try:
            definition = getattr(sio.loadmat(os.path.join(exp_dir, f), squeeze_me=True,
                                             struct_as_record=False, variable_names=["block"])["block"],
                                 "expDef", None)
        except Exception:  # incl. v7.3 Blocks; rare without a parameters file
            continue
        if isinstance(definition, str) and definition:
            return os.path.splitext(os.path.basename(definition.replace("\\", "/")))[0]
    return None


def mpep_xfile(exp_dir, files):
    """Stimulus x-file of an mpep experiment (e.g. stimSparseNoiseAsync_Luminance), or None.

    Protocol.mat (Protocol.xfile, v5 or v7.3), else the first line of the
    <mouse>_<date>_<n>.p stimulus file.
    """
    xfile = None
    if "Protocol.mat" in files:
        path = os.path.join(exp_dir, "Protocol.mat")
        try:
            xfile = sio.loadmat(path, squeeze_me=True, struct_as_record=False)["Protocol"].xfile
        except NotImplementedError:  # v7.3 (HDF5): char array
            import h5py
            try:
                with h5py.File(path, "r") as f:
                    xfile = "".join(chr(c) for c in np.asarray(f["Protocol"]["xfile"][()]).ravel())
            except Exception:
                xfile = None
        except Exception:
            xfile = None
    if not xfile:
        p = [x for x in files if x.endswith(".p")]
        if p:
            try:
                with open(os.path.join(exp_dir, p[0]), "rb") as fh:
                    line = fh.readline().decode(errors="replace").strip()
                xfile = line if line.endswith(".x") else None
            except Exception:
                xfile = None
    return re.sub(r"\.x$", "", str(xfile).strip()) if xfile else None


def read_experiments(date_dir):
    """Rigbox experiments of one day: [{n, expdef, start, end}] (ISO strings or None)."""
    experiments = []
    if not os.path.isdir(date_dir):
        return experiments
    for sub in sorted((s for s in os.listdir(date_dir) if s.isdigit()), key=int):
        exp_dir = os.path.join(date_dir, sub)
        if not os.path.isdir(exp_dir):
            continue
        files = os.listdir(exp_dir)
        block = [x for x in files if x.endswith("_Block.mat")]
        timeline = [x for x in files if x.endswith("_Timeline.mat")]
        if not block and not timeline:
            continue  # empty / aborted experiment folder
        name = expdef_name(exp_dir, files) or mpep_xfile(exp_dir, files)
        if name is None:  # other Timeline-based experiments (EB mice), recognised by their log file
            lower = [x.lower() for x in files]
            if any("vrbehavior" in x or "vr_behavior" in x for x in lower):
                name = "(VR behaviour log)"
            elif "runtraining.mat" in lower:
                name = "RunTraining"
        if name is None and timeline and not block:
            name = "(no expDef, Timeline only)"
        start = end = None
        for f, var in ([(b, "block") for b in block] + [(t, "Timeline") for t in timeline]):
            try:
                start, end = mat_times(os.path.join(exp_dir, f), var)
            except Exception:
                start = end = None
            if start is not None:
                break
        experiments.append({"n": int(sub), "expdef": name or "(unknown)",
                            "start": start.isoformat() if start else None,
                            "end": end.isoformat() if end else None})
    return experiments


def date_dir_of(ks_dir):
    r"""\\server\Subjects\<mouse>\<date> of a recording path."""
    parts = ks_dir.split("\\")
    return "\\".join(parts[:6]) if len(parts) >= 6 else None


def recording_name(meta_path, ks_dir):
    """SpikeGLX run name (e.g. AV009_2022-04-20_ActivePassive_day2), else the folder above the sorter."""
    if meta_path:
        return re.sub(r"_g\d+_t\d+\.imec\d*\.ap\.meta$", "", os.path.basename(meta_path))
    parts = ks_dir.split("\\")
    return parts[-3] if len(parts) >= 3 else ""


def experiments_during(meta, experiments, offset_h=0):
    """Experiments that overlap the recording (by >= MIN_OVERLAP; or, without
    an end time, start inside it), with the ephys clock shifted by offset_h hours."""
    if not (meta.get("start") and meta.get("duration_s")):
        return []
    t0 = dt.datetime.fromisoformat(meta["start"]) + dt.timedelta(hours=offset_h)
    t1 = t0 + dt.timedelta(seconds=meta["duration_s"])
    inside = []
    for e in experiments:
        if not e["start"]:
            continue
        s = dt.datetime.fromisoformat(e["start"])
        if e.get("end"):
            overlap = min(t1, dt.datetime.fromisoformat(e["end"])) - max(t0, s)
            if overlap >= MIN_OVERLAP:
                inside.append(e)
        elif t0 - WINDOW_SLACK <= s <= t1:
            inside.append(e)
    return inside


def clock_offsets(metas, experiments):
    """Whole-hour offset of the ephys clock per recording of one day.

    Recordings that overlap an experiment get 0. For the others (e.g. an ephys
    computer set to another time zone) one offset is chosen together: the one
    that matches most of them; 0 if none matches.
    """
    offsets = {k: 0 for k in metas}
    unmatched = {k: m for k, m in metas.items() if not experiments_during(m, experiments, 0)}
    if not experiments or not unmatched:
        return offsets

    def n_matched(h):
        return sum(bool(experiments_during(m, experiments, h)) for m in unmatched.values())

    best = max(CLOCK_OFFSETS_H, key=lambda h: (n_matched(h), -abs(h)))
    if n_matched(best) > 0:
        offsets.update({k: best for k in unmatched})
    return offsets


# ── merged-data geometry ─────────────────────────────────────────────────────


def channel_geometry(merged_session_dir):
    path = os.path.join(merged_session_dir, "channel_positions.npy")
    if not os.path.isfile(path):
        return {}
    pos = np.load(path)
    shanks = sorted({int(round(x / SHANK_PITCH_UM)) for x in pos[:, 0]})
    return {"shanks": " ".join(map(str, shanks)),
            "depth_min_um": float(pos[:, 1].min()), "depth_max_um": float(pos[:, 1].max())}


# ── main ─────────────────────────────────────────────────────────────────────


def category_of(names):
    cats = []
    for n in names:
        c = EXPDEF_CATEGORIES.get(n, "UNMAPPED")
        if c not in cats:
            cats.append(c)
    return "; ".join(cats)


def categories_with_name_fallback(expdefs, name):
    """(categories, source). The expDefs win: if any of them is informative,
    only those categories are used (a Timeline-only experiment alongside them
    is left out). Only if none is (Timeline only / unmapped / no experiment)
    does the recording name decide, if it says anything."""
    cats = [c for c in category_of(expdefs).split("; ") if c]
    informative = [c for c in cats if c not in UNINFORMATIVE_CATEGORIES]
    if informative:
        return "; ".join(informative), "expDef"
    lowered = re.sub(r"[^a-z]", "", name.lower())
    from_name = list(dict.fromkeys(c for k, c in RECORDING_NAME_CATEGORIES.items() if k in lowered))
    if from_name:
        return "; ".join(from_name), "recording name"
    return "; ".join(cats), "expDef" if cats else ""


def unique_in_order(items):
    return "; ".join(dict.fromkeys(items))


def load_location(key):
    """Sessions of one location from the DeepUnitMatch output on the merged data."""
    mouse, probe, loc = key.split("/")
    path = os.path.join(cfg.ANALYSIS_OUTPUT, mouse, probe, loc, "DeepUnitMatch", "UMparam.pickle")
    with open(path, "rb") as f:
        um = pickle.load(f)
    all_dirs = [str(p) for p in um.get("KS_dirs_all", um["KS_dirs"])]
    used = [str(p) for p in um["KS_dirs"]]
    n_units = dict(zip(used, um.get("n_units_per_session", [None] * len(used))))
    return [{"merged_dir": d, "folder": int(os.path.basename(d.rstrip("\\/"))),
             "used": d in used, "n_units": n_units.get(d)} for d in all_dirs]


def training_sessions(loc_key):
    """[[session folder, n snippet files]] of one location in the training cache, in
    the order fine-tuning loaded them (numeric; sessions without files are skipped)."""
    loc_dir = os.path.join(cfg.TRAINING_CACHE, *loc_key.split("/"))
    if not os.path.isfile(os.path.join(loc_dir, ".done")):
        return None
    data_dir = os.path.join(loc_dir, "processed_waveforms")
    out = []
    for s in sorted(os.listdir(data_dir), key=lambda s: int(s) if s.isdigit() else s):
        n = len(os.listdir(os.path.join(data_dir, s)))
        if n:
            out.append([s, n])
    return out


def finetune_split(loc_keys, train_mice, split_path, listing):
    """{(loc_key, session index): True (fine-tuned on) / False (held out for validation)}.

    session_split.json numbers the sessions in the order train_paper_models'
    finetune_dataset loaded them: locations of the training mice sorted by
    'mouse/probe/loc', sessions numerically within each. Returns ({}, reason)
    if that order can't be rebuilt exactly.
    """
    if not os.path.isfile(split_path):
        return {}, f"no session split at {split_path}"
    with open(split_path, "r", encoding="utf-8") as f:
        split = json.load(f)
    train_ids = {int(s) for s in split["train_sessions"]}
    val_ids = {int(s) for s in split["val_sessions"]}
    order = []
    for key in sorted(k for k in loc_keys if k.split("/")[0] in train_mice):
        sessions = listing.get(key)
        if not isinstance(sessions, list):
            return {}, f"training cache of {key} not readable"
        order += [(key, int(s)) for s, _ in sessions]
    if len(order) != len(train_ids) + len(val_ids):
        return {}, (f"training cache has {len(order)} sessions, the split {len(train_ids) + len(val_ids)}: "
                    f"cannot rebuild the order")
    return {k: (i in train_ids) for i, k in enumerate(order)}, None


def parse_args():
    p = argparse.ArgumentParser(description="Build the dataset description tables of the paper.")
    p.add_argument("--out", default=os.path.join(cfg.RESULTS_DIR, "dataset_table"))
    p.add_argument("--workers", type=int, default=16, help="parallel network reads")
    p.add_argument("--mice", nargs="+", help="only these mice (for a quick test)")
    p.add_argument("--split", default=os.path.join(
        os.path.dirname(_HERE), "DeepUnitMatch", "ModelExp", "experiments", "DUM_paper", "ckpt",
        "session_split.json"), help="session_split.json of the paper model's fine-tuning "
                                    "(on the training machine)")
    return p.parse_args()


def main():
    args = parse_args()
    os.makedirs(args.out, exist_ok=True)
    cache = Cache(os.path.join(args.out, "_cache.json"))

    with open(cfg.METADATA_INDEX_PATH, "r", encoding="utf-8") as f:
        index = json.load(f)
    all_loc_keys = list(index)  # the fine-tuning order needs every location, also with --mice
    if args.mice:
        index = {k: v for k, v in index.items() if v["mouse"] in args.mice}
    print(f"{len(index)} locations in {cfg.METADATA_INDEX_PATH}")

    # sessions of every location
    sessions, problems = [], []
    for key, entry in index.items():
        try:
            used = sorted((s for s in load_location(key) if s["used"]), key=lambda s: s["folder"])
            for i, s in enumerate(used):
                exp_key = str(s["folder"] + 1)
                sessions.append(dict(
                    mouse=entry["mouse"], probe=entry["probe"], loc=entry["loc"], **s,
                    loc_key=key, index_in_loc=i,  # = session number in the matching / training cache
                    date=entry["dates"].get(exp_key, ""),
                    ks_dir=entry["exp_paths"].get(exp_key, ""),
                ))
        except Exception as e:
            problems.append(f"{key}: {e!r}")
    # sessions the pipeline dropped (no units left) are left out. A location with
    # one session has no session pair, so it is in none of the across-session
    # analyses (fast_testing results), but it was used to train the model.
    per_loc = {}
    for s in sessions:
        per_loc[s["loc_key"]] = per_loc.get(s["loc_key"], 0) + 1
    for s in sessions:
        s["used_in_analysis"] = per_loc[s["loc_key"]] >= 2
    print(f"{len(sessions)} sessions ({sum(s['used_in_analysis'] for s in sessions)} in the "
          f"across-session analyses; {sum(n < 2 for n in per_loc.values())} single-session location(s))")

    def safe(fn, arg):
        try:
            return fn(arg)
        except Exception as e:
            return {"error": repr(e)}

    # slow network reads, in parallel and cached
    def run_cached(kind, keys, fn, label):
        todo = [k for k in dict.fromkeys(keys) if k and k not in cache.data.get(kind, {})]
        print(f"{label}: {len(todo)} to read ({len(set(keys)) - len(todo)} cached)")
        with cf.ThreadPoolExecutor(args.workers) as pool:
            futures = {pool.submit(safe, fn, k): k for k in todo}
            for i, fut in enumerate(cf.as_completed(futures), 1):
                cache.data.setdefault(kind, {})[futures[fut]] = fut.result()
                if i % 200 == 0:
                    print(f"  {i}/{len(todo)}")
                    cache.save()
        cache.save()

    run_cached("meta", [s["ks_dir"] for s in sessions], read_meta, "SpikeGLX meta files")
    run_cached(EXP_CACHE, [date_dir_of(s["ks_dir"]) for s in sessions if s["ks_dir"]],
               read_experiments, "Rigbox day folders")
    run_cached("geometry", [s["merged_dir"] for s in sessions], channel_geometry, "channel positions")

    # training of the paper model: mice from model.json, held-out fine-tuning sessions from the split
    train_mice = set()
    model_json = os.path.join(cfg.MODELS_ROOT, "DUM_paper", "model.json")
    if os.path.isfile(model_json):
        with open(model_json, "r", encoding="utf-8") as f:
            train_mice = set((json.load(f).get("job") or {}).get("train_mice", []))
    run_cached("training_sessions", [k for k in all_loc_keys if k.split("/")[0] in train_mice],
               training_sessions, "training cache")
    finetuned, split_problem = finetune_split(all_loc_keys, train_mice, args.split,
                                              cache.data.get("training_sessions", {}))
    if split_problem:
        problems.append(f"used_in_model left empty: {split_problem}")

    def experiments_of(date_dir):
        experiments = cache.data[EXP_CACHE].get(date_dir or "", [])
        if isinstance(experiments, dict):  # read error
            problems.append(f"{date_dir}: experiments {experiments.get('error')}")
            return []
        return experiments

    # ephys clock offset per day, from all recordings of that day (every location)
    metas_by_day = {}
    for s in sessions:
        if s["ks_dir"]:
            metas_by_day.setdefault(date_dir_of(s["ks_dir"]), {})[s["ks_dir"]] = cache.data["meta"].get(s["ks_dir"], {})
    offsets = {}
    for d, m in metas_by_day.items():
        offsets.update(clock_offsets(m, experiments_of(d)))

    rows = []
    for s in sessions:
        meta = cache.data["meta"].get(s["ks_dir"], {}) if s["ks_dir"] else {}
        day = date_dir_of(s["ks_dir"]) if s["ks_dir"] else None
        experiments = experiments_of(day)
        offset_h = offsets.get(s["ks_dir"], 0)
        geom = cache.data["geometry"].get(s["merged_dir"], {})
        matched = experiments_during(meta, experiments, offset_h)
        expdefs = [e["expdef"] for e in matched]
        if matched:
            how = "recording window" + (f" (ephys clock {offset_h:+d} h)" if offset_h else "")
        elif experiments:
            how = "none during recording"
        else:
            how = "no experiments that day"
        code = meta.get("probe_type_code")
        rec_name = recording_name(meta.get("meta_path"), s["ks_dir"]) if s["ks_dir"] else ""
        categories, category_source = categories_with_name_fallback(expdefs, rec_name)
        rows.append(dict(
            mouse=s["mouse"], probe=s["probe"], loc=s["loc"], session=s["folder"],
            used_in_model=(finetuned.get((s["loc_key"], s["index_in_loc"])) if finetuned
                           else (None if s["mouse"] in train_mice else False)),
            used_in_analysis=s["used_in_analysis"],
            date=s["date"],
            recording_name=rec_name,
            start_time=(meta.get("start") or "")[11:19],
            duration_min=round(meta["duration_s"] / 60, 1) if meta.get("duration_s") else np.nan,
            probe_type=PROBE_TYPES.get(code, f"imDatPrb_type {code}" if code else ""),
            probe_sn=meta.get("probe_sn", ""),
            sampling_rate_hz=meta.get("sampling_rate", np.nan),
            shanks=geom.get("shanks", ""), depth_min_um=geom.get("depth_min_um", np.nan),
            depth_max_um=geom.get("depth_max_um", np.nan),
            n_good_units=s["n_units"],
            expdefs=unique_in_order(expdefs), expdef_categories=categories,
            category_source=category_source, expdef_match=how,
            expdefs_same_day=unique_in_order(e["expdef"] for e in experiments),
            notes="; ".join(x for x in (meta.get("error"), geom.get("error"),
                                        None if meta.get("meta_path") else "no .ap.meta found") if x),
            ks_dir=s["ks_dir"], meta_path=meta.get("meta_path") or "",
        ))
    df = pd.DataFrame(rows)
    df["date_dt"] = pd.to_datetime(df["date"], errors="coerce")
    first_loc = df.groupby(["mouse", "probe", "loc"])["date_dt"].transform("min")
    first_mouse = df.groupby("mouse")["date_dt"].transform("min")
    df["day_at_location"] = (df["date_dt"] - first_loc).dt.days
    df["day_since_first_recording"] = (df["date_dt"] - first_mouse).dt.days
    df = df.sort_values(["mouse", "probe", "loc", "date_dt", "session"])
    cols = ["mouse", "probe", "loc", "session", "used_in_model", "used_in_analysis", "date", "day_at_location",
            "day_since_first_recording", "recording_name", "start_time", "duration_min",
            "probe_type", "probe_sn", "sampling_rate_hz", "shanks", "depth_min_um", "depth_max_um",
            "n_good_units", "expdefs", "expdef_categories", "category_source", "expdef_match",
            "expdefs_same_day", "notes",
            "ks_dir", "meta_path"]
    df[cols].to_csv(os.path.join(args.out, "dataset_sessions.csv"), index=False)

    # per location
    used = df

    def join_unique(series, sep="; "):
        items = [x for v in series.dropna() for x in str(v).split(sep) if x]
        return sep.join(dict.fromkeys(items))

    locations = used.groupby(["mouse", "probe", "loc"]).agg(
        probe_type=("probe_type", join_unique),
        probe_sn=("probe_sn", join_unique),
        shanks=("shanks", lambda v: " ".join(map(str, sorted({int(x) for s in v.dropna() for x in str(s).split()})))),
        depth_min_um=("depth_min_um", "min"), depth_max_um=("depth_max_um", "max"),
        used_in_analysis=("used_in_analysis", "all"),
        n_sessions=("session", "size"),
        n_sessions_used_in_model=("used_in_model", lambda v: int((v == True).sum())),  # noqa: E712
        n_days=("date", "nunique"),
        first_date=("date", "min"), last_date=("date", "max"),
        span_days=("day_at_location", "max"),
        median_good_units=("n_good_units", "median"),
        total_recording_h=("duration_min", lambda v: round(v.sum() / 60, 1)),
        expdef_categories=("expdef_categories", join_unique),
    ).reset_index()
    locations["brain_region"] = ""  # fill in by hand
    locations.to_csv(os.path.join(args.out, "dataset_locations.csv"), index=False)

    # per mouse
    mice = used.groupby("mouse").agg(
        probe_types=("probe_type", join_unique),
        n_locations=("loc", lambda v: used.loc[v.index, ["probe", "loc"]].drop_duplicates().shape[0]),
        n_sessions=("session", "size"),
        first_date=("date", "min"), last_date=("date", "max"),
        span_days=("day_since_first_recording", "max"),
        expdef_categories=("expdef_categories", join_unique),
    ).reset_index()
    for col in ("sex", "strain", "date_of_birth", "implant_date", "data_source"):
        mice[col] = ""  # fill in by hand (e.g. Alyx)
    mice["in_paper_model_training"] = mice["mouse"].isin(train_mice)
    mice.to_csv(os.path.join(args.out, "dataset_mice.csv"), index=False)

    # expDef names, to check the category mapping
    names = [n for v in used["expdefs"] for n in str(v).split("; ") if n]
    expdef_table = (pd.Series(names).value_counts().rename_axis("expdef").reset_index(name="n_sessions"))
    expdef_table["category"] = expdef_table["expdef"].map(EXPDEF_CATEGORIES).fillna("UNMAPPED")
    expdef_table.to_csv(os.path.join(args.out, "expdef_names.csv"), index=False)

    # summary
    print(f"\nWritten to {args.out}")
    print(f"  {len(used)} sessions, {len(locations)} locations, {len(mice)} mice")
    print(f"  used_in_analysis: {int(used['used_in_analysis'].sum())} sessions; used_in_model: "
          f"{used['used_in_model'].value_counts(dropna=False).to_dict()}")
    print("  expdef_match:", used["expdef_match"].value_counts().to_dict())
    print(f"  sessions without .ap.meta: {(used['meta_path'] == '').sum()}")
    unmapped = expdef_table[expdef_table["category"] == "UNMAPPED"]
    if len(unmapped):
        print(f"  UNMAPPED expDefs ({len(unmapped)}), add them to EXPDEF_CATEGORIES:")
        print(unmapped.to_string(index=False))
    if problems:
        print(f"\n{len(problems)} problem(s):")
        for p in problems[:50]:
            print("  " + p)


if __name__ == "__main__":
    try:
        main()
    except Exception:
        traceback.print_exc()
        raise SystemExit(1)
