import numpy as np
import pandas as pd
import os
import sys
import pickle
import shutil
import datetime
import h5py
import scipy.io
from pathlib import Path
import matplotlib.pyplot as plt

sys.path.insert(0, os.getcwd())
sys.path.insert(0, os.path.join(os.getcwd(), "testing"))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import traceback

import batch_lock
import pipeline_config as cfg
import pipeline_log as plog

LOG_STAGE = "merge"


def group_key(source_dir):
    """'mouse/probe/location' of a step-1 output folder (<X>/DeepUnitMatch)."""
    return os.path.relpath(os.path.dirname(source_dir), DUM_NONMERGED_DATAPATH).replace(os.sep, "/")


# Input: step-1 output on the non-merged data (one DeepUnitMatch/ folder per
# recording location, holding UMparam.pickle + MatchTable.csv).
DUM_NONMERGED_DATAPATH = cfg.UNMERGED_OUTPUT
# Output: the merged KS-style tree, mirroring DUM_NONMERGED_DATAPATH's layout.
MERGED_DATAPATH = cfg.MERGED_DATA
# Raw KS root the original UnitMatch.mat files were generated from (same as
# run_deepunitmatch_batch.py's BASE_INPUT). Needed to correctly map a session's
# position in the full KS_dirs list to its RecSes number in MatchTable.csv --
# see original_index_to_recses() below.
RAW_KS_BASE = cfg.RAW_KS_BASE
# A candidate pair is merged when merging does not increase contamination:
# C(merged spike train) / C(larger unit alone) <= MAX_C_RATIO. "<=" rather than
# "<" because estimate_C clamps zero contamination to 0.01, so two units with
# no refractory violations at all (merged or not) give a ratio of exactly 1 and
# should still be merged.
MAX_C_RATIO = 1.0
# Merged-away units keep a backup of their original waveform here (a sibling
# of RawWaveforms, so downstream loaders that read RawWaveforms never see it).
BACKUP_WAVEFORM_DIRNAME = "RawWaveforms_premerge_backup"
# Label written into the synthetic cluster_bc_unitType.tsv for units that were
# absorbed by a merge, so downstream code excludes them via the label file
# instead of relying on their waveform file being missing.
MERGED_UNIT_LABEL = "MERGED"

# Written into a session's target folder once it has been fully copied/merged.
# Lets re-runs skip sessions that are already done instead of re-copying
# waveforms and redoing the merge over the network every time.
MERGE_COMPLETE_MARKER = "merge_complete.flag"
# See batch_lock.sentinel_is_fresh() / run_deepunitmatch_batch_onMerged.py's
# REDO_FROM_DATE for what this does: a session is skipped once its
# merge_complete.flag exists and is at least this new. None falls back to
# plain "skip if present"; a far-future date reproduces old REDO=True.
#
# This matters more here than in the other batch scripts: this script's
# within-session merge decisions are read from run_deepunitmatch_batch.py's
# DeepUnitMatch/MatchTable.csv (UID 1 == UID 2 pairs, see run_merging_process
# below), which depends on DUM's Naive Bayes match probabilities. A fix to
# DUM's matching (e.g. its adaptive prior/threshold) changes that MatchTable,
# so previously-merged sessions must be treated as stale even though their own
# merge_complete.flag file didn't change.
REDO_FROM_DATE = cfg.REDO_FROM_DATE


def _decode_hdf5_str(f, ref_or_ds):
    """Decode a MATLAB char array stored as uint16 in an HDF5 file."""
    obj = f[ref_or_ds] if isinstance(ref_or_ds, h5py.Reference) else ref_or_ds
    chars = obj[()].flatten()
    return "".join(chr(int(c)) for c in chars)


def _load_uid_conversion_scipy(mat_path):
    """Load via scipy (MATLAB < v7.3). Returns (orig_clus_id, recsesAll, good_id)."""
    mat = scipy.io.loadmat(mat_path, simplify_cells=True)
    uid = mat["UniqueIDConversion"]
    orig_clus_id = np.array(uid["OriginalClusID"]).flatten()
    recsesAll = np.array(uid["recsesAll"]).flatten()
    good_id = np.array(uid["GoodID"]).flatten().astype(bool)
    return orig_clus_id, recsesAll, good_id


def _load_uid_conversion_hdf5(mat_path):
    """Load via h5py (MATLAB v7.3 HDF5). Returns (orig_clus_id, recsesAll, good_id)."""
    with h5py.File(mat_path, "r") as f:
        uid = f["UniqueIDConversion"]
        orig_clus_id = uid["OriginalClusID"][()].flatten()
        recsesAll = uid["recsesAll"][()].flatten()
        good_id = uid["GoodID"][()].flatten().astype(bool)
    return orig_clus_id, recsesAll, good_id


def load_uid_conversion(mat_path):
    """Load UnitMatch.mat's UniqueIDConversion. Returns (orig_clus_id, recsesAll, good_id)."""
    try:
        return _load_uid_conversion_scipy(mat_path)
    except Exception:
        return _load_uid_conversion_hdf5(mat_path)


def original_index_to_recses(recsesAll, good_id, n_ks_dirs):
    """
    Map each 0-based position in the full KS_dirs list to the 1-based *compacted*
    RecSes number used in MatchTable.csv.

    MatchTable.csv's RecSes numbering only counts sessions that had at least one
    good unit in the original DeepUnitMatch comparison -- sessions with zero good
    units are dropped, not just left empty. So a session's position in the full
    KS_dirs list does not equal its RecSes number whenever an earlier session had
    zero good units. This uses the original UnitMatch.mat's recsesAll/GoodID
    (the authoritative, uncompacted source) to build the correct mapping.
    Sessions with zero good units are simply absent from the returned dict.
    """
    good_session_idx = [
        i for i in range(n_ks_dirs) if ((recsesAll == (i + 1)) & good_id).sum() > 0
    ]
    return {
        orig_idx: recses for recses, orig_idx in enumerate(good_session_idx, start=1)
    }


def step1_sessions(data, recsesAll, good_id):
    """
    (full KS_dirs list, {0-based position in it -> RecSes}) of one step-1 run.

    Step-1 runs save the mapping themselves: param["KS_dirs_all"] is the full
    list and param["session_index"][k] the position of RecSes k+1 in it. This
    also covers sessions the run dropped for other reasons than having no good
    units (e.g. all units with non-finite waveforms). Older outputs without
    it fall back to the GoodID reconstruction (original_index_to_recses).
    """
    if "session_index" in data:
        return data["KS_dirs_all"], {
            int(orig_idx): recses for recses, orig_idx in enumerate(data["session_index"], start=1)
        }
    return data["KS_dirs"], original_index_to_recses(recsesAll, good_id, len(data["KS_dirs"]))


def write_synthetic_bc_unit_type_tsv(
    orig_clus_id, recsesAll, good_id, original_session_idx, out_path
):
    """
    Write a cluster_bc_unitType.tsv-compatible file for one session, derived from
    the original UnitMatch.mat's GoodID -- always, not just when the session's own
    bombcell output is missing. This matches run_deepunitmatch_batch.py, which
    never reads cluster_bc_unitType.tsv either: it defines good units purely from
    UniqueIDConversion. A raw KS session's own bombcell labels are session-wide
    and can't be used directly here, since the same KS session may be a candidate
    in more than one probe/depth-group's UnitMatch.mat (e.g. two different depths
    recorded on different days sharing an earlier day's session); GoodID is
    correctly scoped to *this* group's comparison, the session-wide TSV is not.

    The mat only records a good/not-good call, not the full bombcell MUA/NOISE
    distinction, so every non-good unit is labelled 'MUA' -- downstream code only
    ever selects GOOD / NON-SOMA GOOD rows, so this coarser labelling of the rest
    has no effect on which units get used.

    original_session_idx is 0-based, matching position in the full KS_dirs list
    (i.e. recsesAll == original_session_idx + 1).
    """
    mask = recsesAll == (original_session_idx + 1)
    df = pd.DataFrame(
        {
            "cluster_id": orig_clus_id[mask].astype(int),
            "bc_unitType": np.where(good_id[mask], "GOOD", "MUA"),
        }
    )
    df.to_csv(out_path, sep="\t", index=False)


def estimate_C(spike_train, t_r=0.002, t_c=0.0001):
    n_v = 0
    for i in range(spike_train.shape[0]):
        spike1 = spike_train[i]
        for spike2 in spike_train[i + 1 :]:
            diff = spike2 - spike1
            if diff < (t_r - t_c):
                if diff > t_c:
                    n_v += 1
            else:
                break
    N = len(spike_train)
    T = max(spike_train) - min(spike_train) - 2 * N * t_c
    if (1 - (2 * n_v * T) / (N**2 * t_r)) > 0:
        C = 0.5 * (1 - np.sqrt(1 - (2 * n_v * T) / (N**2 * t_r)))
        if C == 0:
            C = 0.01
    else:
        C = 1
    return C


def merge_waveforms(
    id1, id2, weight1, weight2, target_waveform_dir, plot_waveforms=False
):

    wave1 = os.path.join(target_waveform_dir, f"Unit{str(int(id1))}_RawSpikes.npy")
    wave2 = os.path.join(target_waveform_dir, f"Unit{str(int(id2))}_RawSpikes.npy")

    # Back up the waveform files (outside RawWaveforms). A unit that absorbs
    # several others keeps the backup of its *original* waveform.
    backup_dir = os.path.join(os.path.dirname(target_waveform_dir), BACKUP_WAVEFORM_DIRNAME)
    os.makedirs(backup_dir, exist_ok=True)
    for wave in (wave1, wave2):
        backup = os.path.join(backup_dir, os.path.basename(wave))
        if not os.path.exists(backup):
            shutil.copy(wave, backup)

    waveform1 = np.load(wave1)
    waveform2 = np.load(wave2)
    new_waveform = (weight1 * waveform1 + weight2 * waveform2) / (weight1 + weight2)

    # Overwrite the waveform of the first unit with the new merged waveform
    save_path = os.path.join(
        target_waveform_dir,
        f"Unit{str(int(id1))}_RawSpikes.npy",
    )
    np.save(save_path, new_waveform)

    # Back the waveform file of the second unit
    os.remove(wave2)

    print(
        f"Merged waveforms of units {id1} and {id2} into unit {id1}, and deleted unit {id2}."
    )

    if plot_waveforms:
        _, ax = plt.subplots(1, 3, figsize=(10, 4))

        ax[0].imshow(np.nanmean(waveform1, axis=2), aspect="auto", cmap="viridis")
        ax[1].imshow(np.nanmean(waveform2, axis=2), aspect="auto", cmap="viridis")
        ax[2].imshow(np.nanmean(new_waveform, axis=2), aspect="auto", cmap="viridis")
        ax[2].set_title(f"Merged Waveform of Unit {id1}")
        ax[2].set_xlabel("Time (samples)")
        ax[2].set_ylabel("Channels")


def revert_waveform_merges(target_waveform_dir):
    # Find all backup waveform files (kept next to RawWaveforms, see merge_waveforms)
    backup_dir = os.path.join(os.path.dirname(target_waveform_dir), BACKUP_WAVEFORM_DIRNAME)
    if not os.path.isdir(backup_dir):
        print("No backups found, nothing to revert.")
        return
    backup_files = [f for f in os.listdir(backup_dir) if f.endswith("_RawSpikes.npy")]

    for backup_file in backup_files:
        backup_path = os.path.join(backup_dir, backup_file)
        original_path = os.path.join(target_waveform_dir, backup_file)
        original_file = backup_file

        # Restore the original waveform file from the backup
        shutil.copy(backup_path, original_path)
        print(f"Restored {original_file} from backup.")

        os.remove(backup_path)  # Optionally remove the backup file after restoration

    print("All merges have been reverted.")


def session_merge_candidates(mt, recses):
    """
    Merge candidates of one session: groups (arrays of unit IDs) of different
    units in session `recses` (RecSes number in MatchTable.csv) that were
    assigned the same unique ID ("UID 1" == "UID 2", intermediate IDs).
    """
    if recses is None:  # session step 1 didn't match (no usable good units): nothing to merge
        return []
    rows = mt[
        (mt["RecSes 1"] == recses)
        & (mt["RecSes 2"] == recses)
        & (mt["ID1"] != mt["ID2"])
        & (mt["UID 1"] == mt["UID 2"])
    ]
    return [
        np.unique(np.concatenate([g["ID1"].values, g["ID2"].values]))
        for _, g in rows.groupby("UID 1")
    ]


def decide_merges(unit_groups, spk_times, spk_clusters, max_c_ratio=None, verbose=True):
    """
    Which candidate units to merge, by refractory-period contamination.

    Within each group, repeatedly take the pair whose merged spike train is
    least contaminated relative to the larger unit alone
    (C(merged) / C(larger unit), estimate_C) and merge it if that ratio is
    <= max_c_ratio (MAX_C_RATIO); stop when no pair qualifies. A merged unit
    keeps the lower cluster ID and can absorb further units of its group.

    Returns (merges, merged_spk_clusters): merges is a list of dicts
    (kept, absorbed, n_kept, n_absorbed = spike counts just before that
    merge, c_ratio) in the order applied; spk_clusters itself is not modified.
    """
    if max_c_ratio is None:
        max_c_ratio = MAX_C_RATIO
    spk_clusters = spk_clusters.copy()
    merges = []
    for unit_indices in unit_groups:
        unit_indices = np.asarray(unit_indices)
        while len(unit_indices) > 1:
            ratios = []
            for unit_id1 in unit_indices:
                times1 = spk_times[np.where(spk_clusters == unit_id1)]
                for unit_id2 in unit_indices:
                    if unit_id1 < unit_id2:  # Avoid duplicate pairs and self-comparison
                        times2 = spk_times[np.where(spk_clusters == unit_id2)]
                        non_merged = times1 if len(times1) > len(times2) else times2
                        merged = np.sort(np.concatenate([times1, times2]))
                        ratios.append(
                            {
                                "unit_id1": unit_id1,
                                "unit_id2": unit_id2,
                                "C_ratio": estimate_C(merged) / estimate_C(non_merged),
                            }
                        )
            # Pair whose merge increases contamination the least. (Positional
            # .iloc: after sort_values the index labels keep their pre-sort
            # order, so label-based [0] would pick the first pair computed.)
            best = pd.DataFrame(ratios).sort_values("C_ratio").iloc[0]
            if best["C_ratio"] > max_c_ratio:
                if verbose:
                    print(f"No more units to merge based on C ratio threshold (lowest C ratio = {best['C_ratio']:.4f}).")
                break
            # (the row mixes int and float columns, so pandas returns the ids as floats)
            id_a, id_b = int(best["unit_id1"]), int(best["unit_id2"])
            kept, absorbed = min(id_a, id_b), max(id_a, id_b)
            if verbose:
                print(f"Merging units {id_a} and {id_b} with C ratio {best['C_ratio']:.4f}")
            # spike counts *before* relabelling (waveform weights)
            n_kept = int(np.sum(spk_clusters == kept))
            n_absorbed = int(np.sum(spk_clusters == absorbed))
            spk_clusters[spk_clusters == absorbed] = kept
            merges.append({"kept": kept, "absorbed": absorbed, "n_kept": n_kept,
                           "n_absorbed": n_absorbed, "c_ratio": float(best["C_ratio"])})
            unit_indices = unit_indices[unit_indices != absorbed]
    return merges, spk_clusters


def get_group_lock_path(target_dir):
    """
    Lock file marking 'a run is currently merging this group', so multiple
    machines pointed at the same merged-data output can split work across
    groups without double-processing one. See batch_lock.py.
    """
    return os.path.join(target_dir, ".processing.lock")


def run_merging_process(UMparam_files, source_dirs, MAX_C_RATIO=1.0):

    for UMparam_file, source_dir in zip(UMparam_files, source_dirs):
        try:
            with open(UMparam_file, "rb") as file:
                data = pickle.load(file)

            # Same <mouse>/<probe>/<location>/DeepUnitMatch layout under the merged root.
            target_dir = os.path.join(
                MERGED_DATAPATH, os.path.relpath(source_dir, DUM_NONMERGED_DATAPATH)
            )

            # every session folder of the location, also those step 1 dropped
            ks_dirs_all = data.get("KS_dirs_all", data["KS_dirs"])
            target_KSDirs = [
                os.path.join(target_dir, str(idx)) for idx in range(len(ks_dirs_all))
            ]
            already_done = [
                batch_lock.sentinel_is_fresh(
                    os.path.join(d, MERGE_COMPLETE_MARKER), REDO_FROM_DATE
                )
                for d in target_KSDirs
            ]
            if all(already_done):
                print(
                    f"Skipping UMparam_file: {UMparam_file} (all {len(target_KSDirs)} session(s) already processed)."
                )
                continue

            lock_path = get_group_lock_path(target_dir)
            with batch_lock.try_lock(lock_path) as acquired:
                if not acquired:
                    print(f"  Skipping (already being processed by another run): {lock_path}")
                    continue

                # re-check now that we hold the lock: another machine may have
                # finished this group while we were scanning/waiting for the lock
                already_done = [
                    batch_lock.sentinel_is_fresh(
                        os.path.join(d, MERGE_COMPLETE_MARKER), REDO_FROM_DATE
                    )
                    for d in target_KSDirs
                ]
                if all(already_done):
                    print(f"  Skipping UMparam_file: {UMparam_file} (completed by another run).")
                    continue

                print(f"Processing UMparam_file: {UMparam_file}...")

                # could load both UM and DUM matchtables here to get the merged units, but for now just use the UM matchtable to find which units to merge
                mt = pd.read_csv(os.path.join(os.path.split(UMparam_file)[0], "MatchTable.csv"))

                # source_dir is DUM_NONMERGED_DATAPATH/<X>/DeepUnitMatch; the matching
                # original UnitMatch.mat lives at RAW_KS_BASE/<X>/UnitMatch/UnitMatch.mat.
                x = os.path.dirname(os.path.relpath(source_dir, DUM_NONMERGED_DATAPATH))
                mat_path = os.path.join(RAW_KS_BASE, x, "UnitMatch", "UnitMatch.mat")
                try:
                    orig_clus_id, recsesAll, good_id = load_uid_conversion(mat_path)
                    ks_dirs_all, idx_to_recses = step1_sessions(data, recsesAll, good_id)
                except Exception as e:
                    print(f"  WARNING: could not read {mat_path} ({e}); skipping {source_dir}.")
                    plog.log_event(LOG_STAGE, group_key(source_dir), "merge", "failed",
                                   f"could not read {mat_path}: {e}")
                    continue

                group_absorbed = 0
                for idx, KSDir in enumerate(ks_dirs_all):
                    target_KSDir = target_KSDirs[idx]
                    marker_path = os.path.join(target_KSDir, MERGE_COMPLETE_MARKER)

                    if already_done[idx]:
                        print(
                            f"Processing KSDir: {KSDir} (idx {idx})... already done, skipping."
                        )
                        continue

                    print(f"Processing KSDir: {KSDir} (idx {idx})...")

                    if not os.path.exists(target_KSDir):
                        os.makedirs(target_KSDir)

                    # Move KSDir files that are necessary but won't be modified to the new target directory
                    for file_name in [
                        "spike_times.npy",
                        "channel_positions.npy",
                        "cluster_info.tsv",
                    ]:
                        source_file = os.path.join(KSDir, file_name)
                        target_file = os.path.join(target_KSDir, file_name)
                        if os.path.exists(source_file):
                            shutil.copy2(source_file, target_file)  # preserves metadata

                    # Move natural images info if exists
                    natural_images_folder = os.path.dirname(os.path.dirname(KSDir))
                    for file_name in ['trial.imageIDs.npy', 'trial.offsetTimes.npy', 'trial.onsetTimes.npy']:
                            source_file = os.path.join(natural_images_folder, file_name)
                            target_file = os.path.join(target_KSDir, file_name)
                            if os.path.exists(source_file):
                                shutil.copy2(source_file, target_file)  # preserves metadata

                    # cluster_bc_unitType.tsv is always *derived* from the original
                    # UnitMatch.mat's GoodID rather than copied from the KS session's own
                    # bombcell output. This matches run_deepunitmatch_batch.py, which never
                    # reads cluster_bc_unitType.tsv either -- it defines good units purely
                    # from UniqueIDConversion. It also avoids a real scoping mismatch: a raw
                    # KS session's bombcell labels are session-wide, but the same session
                    # can be a candidate for multiple probe/depth-group comparisons (see
                    # write_synthetic_bc_unit_type_tsv docstring) -- GoodID is the one
                    # source that's correctly scoped to *this* group's comparison.
                    write_synthetic_bc_unit_type_tsv(
                        orig_clus_id,
                        recsesAll,
                        good_id,
                        idx,
                        os.path.join(target_KSDir, "cluster_bc_unitType.tsv"),
                    )

                    source_waveform_dir = next(
                        (
                            str(path)
                            for path in Path(KSDir).rglob("RawWaveforms")
                            if path.is_dir()
                        ),
                        None,
                    )
                    target_waveform_dir = os.path.join(target_KSDir, "qMetrics", "RawWaveforms")

                    spk_clusters = np.load(os.path.join(KSDir, "spike_clusters.npy"))
                    spk_times = np.load(os.path.join(KSDir, "spike_times.npy"))
                    spk_times = spk_times / 30000  # convert to seconds

                    # Copy the RawWaveforms directory to the new target directory.
                    # Start from a clean copy: a session being (re)processed must not
                    # keep merged waveforms or backups from an earlier, interrupted or
                    # outdated run (these only ever live in the *target* tree).
                    backup_waveform_dir = os.path.join(
                        os.path.dirname(target_waveform_dir), BACKUP_WAVEFORM_DIRNAME
                    )
                    for stale_dir in (target_waveform_dir, backup_waveform_dir):
                        if os.path.isdir(stale_dir):
                            shutil.rmtree(stale_dir)
                    os.makedirs(target_waveform_dir)

                    for file in Path(source_waveform_dir).iterdir():
                        if file.is_file():
                            shutil.copy2(
                                file, os.path.join(target_waveform_dir, file.name)
                            )  # preserves metadata

                    # Find which units need to be merged based on the matching results.
                    # RecSes in MatchTable.csv is the *compacted* session count (sessions
                    # with zero good units dropped), so idx+1 is only correct when no
                    # earlier session had zero good units -- use the recovered mapping
                    # instead of assuming idx+1 == RecSes.
                    # (None: zero good units in the original comparison -> no candidates)
                    recses = idx_to_recses.get(idx)
                    unit_groups = session_merge_candidates(mt, recses)
                    merges, spk_clusters = decide_merges(unit_groups, spk_times, spk_clusters)
                    absorbed_units = []
                    for m in merges:
                        total = m["n_kept"] + m["n_absorbed"]
                        merge_waveforms(
                            m["kept"],
                            m["absorbed"],
                            m["n_kept"] / total,
                            m["n_absorbed"] / total,
                            target_waveform_dir,
                            plot_waveforms=False,
                        )
                        absorbed_units.append(m["absorbed"])

                    # Absorbed units no longer exist: relabel them in the synthetic
                    # label file so every downstream loader excludes them.
                    if absorbed_units:
                        tsv_path = os.path.join(target_KSDir, "cluster_bc_unitType.tsv")
                        labels = pd.read_csv(tsv_path, sep="\t")
                        labels.loc[
                            labels["cluster_id"].isin(absorbed_units), "bc_unitType"
                        ] = MERGED_UNIT_LABEL
                        labels.to_csv(tsv_path, sep="\t", index=False)
                        print(f"  {len(absorbed_units)} unit(s) absorbed by merges in this session.")
                        group_absorbed += len(absorbed_units)

                    # Save the updated spike_clusters array to the new target directory
                    target_spike_clusters_file = os.path.join(
                        target_KSDir, "spike_clusters.npy"
                    )
                    np.save(target_spike_clusters_file, spk_clusters)

                    # Mark this session complete so a re-run can skip it. If anything
                    # above raised, execution never reaches here, so the session stays
                    # unmarked and gets fully reprocessed next time.
                    with open(marker_path, "w") as f:
                        f.write("ok")

                plog.log_event(
                    LOG_STAGE, group_key(source_dir), "merge", "done",
                    f"{len(ks_dirs_all)} session(s), {group_absorbed} unit(s) absorbed by merges",
                )
            print(f"Processing UMparam_file: {UMparam_file} done.")
        except Exception as e:
            # one failing location must not stop the others: log it and move on
            traceback.print_exc()
            plog.log_event(
                LOG_STAGE, group_key(source_dir), "merge", "failed",
                f"{type(e).__name__}: {e}", traceback.format_exc(),
            )


def _same_unit_pairs(unit_groups):
    """Unordered pairs of units that a list of unit groups puts together."""
    pairs = set()
    for group in unit_groups:
        group = [int(u) for u in group]
        pairs.update(frozenset((a, b)) for i, a in enumerate(group) for b in group[i + 1 :])
    return pairs


def _merged_pairs(merges):
    """Unordered pairs of units that end up in the same merged unit."""
    final = {}
    for m in merges:  # in order: a unit absorbed later may itself have absorbed others
        final[m["absorbed"]] = m["kept"]
    def root(u):
        while u in final:
            u = final[u]
        return u
    by_root = {}
    for u in set(final) | set(final.values()):
        by_root.setdefault(root(u), []).append(u)
    return _same_unit_pairs(by_root.values())


def compare_candidates(out_csv):
    """
    REPORT ONLY (nothing is merged or copied): for every location, apply the
    merge rule to DeepUnitMatch's merge candidates (as in the real run) and
    to UMPy's (same rule, UMPy's unique IDs), using the step-1 outputs and
    the raw spike times, and record candidates, merges and their overlap.
    Rows are appended to out_csv as locations finish.
    """
    rows = []
    for umparam_file, source_dir in _find_step1_outputs():
        group = group_key(source_dir)
        umpy_dir = os.path.join(os.path.dirname(source_dir), "UMPy")
        try:
            with open(umparam_file, "rb") as f:
                data = pickle.load(f)
            mt = {
                "DUM": pd.read_csv(os.path.join(source_dir, "MatchTable.csv"), usecols=["RecSes 1", "RecSes 2", "ID1", "ID2", "UID 1", "UID 2"]),
                "UM": pd.read_csv(os.path.join(umpy_dir, "MatchTable.csv"), usecols=["RecSes 1", "RecSes 2", "ID1", "ID2", "UID 1", "UID 2"]),
            }
            x = os.path.dirname(os.path.relpath(source_dir, DUM_NONMERGED_DATAPATH))
            orig_clus_id, recsesAll, good_id = load_uid_conversion(os.path.join(RAW_KS_BASE, x, "UnitMatch", "UnitMatch.mat"))
            ks_dirs_all, idx_to_recses = step1_sessions(data, recsesAll, good_id)
            counts = {k: 0 for k in ["cand_DUM", "cand_UM", "cand_both", "merged_DUM", "merged_UM", "merged_both",
                                     "absorbed_DUM", "absorbed_UM"]}
            for idx, ks_dir in enumerate(ks_dirs_all):
                recses = idx_to_recses.get(idx)
                groups = {k: session_merge_candidates(mt[k], recses) for k in mt}
                if not any(groups.values()):
                    continue
                spk_clusters = np.load(os.path.join(ks_dir, "spike_clusters.npy"))
                spk_times = np.load(os.path.join(ks_dir, "spike_times.npy")) / 30000
                cand = {k: _same_unit_pairs(g) for k, g in groups.items()}
                merges = {k: decide_merges(g, spk_times, spk_clusters, verbose=False)[0] for k, g in groups.items()}
                merged = {k: _merged_pairs(m) for k, m in merges.items()}
                for k in mt:
                    counts[f"cand_{k}"] += len(cand[k])
                    counts[f"merged_{k}"] += len(merged[k])
                    counts[f"absorbed_{k}"] += len(merges[k])
                counts["cand_both"] += len(cand["DUM"] & cand["UM"])
                counts["merged_both"] += len(merged["DUM"] & merged["UM"])
            row = {"group": group, **counts}
            print(f"{group}: {counts}")
        except Exception as e:
            traceback.print_exc()
            row = {"group": group, "error": f"{type(e).__name__}: {e}"}
        rows.append(row)
        pd.DataFrame(rows).to_csv(out_csv, index=False)

    df = pd.DataFrame(rows)
    tot = df.drop(columns=["group"] + (["error"] if "error" in df else [])).sum()
    pct = lambda a, b: 100 * a / b if b else float("nan")
    print("\nTotals over all locations:")
    print(f"  candidate pairs: DUM {tot['cand_DUM']:.0f}, UM {tot['cand_UM']:.0f}, both {tot['cand_both']:.0f}")
    print(f"  merged pairs:    DUM {tot['merged_DUM']:.0f}, UM {tot['merged_UM']:.0f}, both {tot['merged_both']:.0f} "
          f"({pct(tot['merged_both'], tot['merged_DUM']):.1f}% of DUM's merges also made from UM's candidates)")
    print(f"  units absorbed:  DUM {tot['absorbed_DUM']:.0f}, UM {tot['absorbed_UM']:.0f}")
    print(f"Per-location table: {out_csv}")


def _find_step1_outputs():
    """(UMparam.pickle, its DeepUnitMatch folder) of every step-1 location."""
    found = []
    for root, _, files in os.walk(DUM_NONMERGED_DATAPATH):
        if "UMparam.pickle" in files and os.path.basename(root) == "DeepUnitMatch":
            found.append((os.path.join(root, "UMparam.pickle"), root))
    return sorted(found)


def main():
    import argparse

    parser = argparse.ArgumentParser(description="Build the merged dataset from the step-1 output.")
    parser.add_argument(
        "--compare-candidates", action="store_true",
        help="Report only: overlap of the merges from DeepUnitMatch's vs UMPy's unique IDs "
             "(nothing is merged; report in pipeline_config.REPORTS_DIR).",
    )
    args = parser.parse_args()
    if args.compare_candidates:
        os.makedirs(cfg.REPORTS_DIR, exist_ok=True)
        stamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
        compare_candidates(os.path.join(cfg.REPORTS_DIR, f"merge_candidates_DUM_vs_UM_{stamp}.csv"))
        return

    # find all UMparam.pickle from DeepUnitMatch folders, to copy their original files to new directory
    UMparam_files = []
    source_dirs = []
    for root, _, files in os.walk(DUM_NONMERGED_DATAPATH):
        for item in files:
            if item.endswith("UMparam.pickle") & ("DeepUnitMatch" in root):
                UMparam_files.append(os.path.join(root, item))
                source_dirs.append(root)

                # need to loop over DUM outputs
    run_merging_process(UMparam_files, source_dirs, MAX_C_RATIO=MAX_C_RATIO)


if __name__ == "__main__":
    main()
