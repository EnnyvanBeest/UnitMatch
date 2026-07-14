import numpy as np
import pandas as pd
import os
import sys
import pickle
import shutil
import h5py
import scipy.io
from pathlib import Path
import matplotlib.pyplot as plt

sys.path.insert(0, os.getcwd())
sys.path.insert(0, os.path.join(os.getcwd(), "testing"))

DUM_NONMERGED_DATAPATH = r'\\znas\Lab\Share\UNITMATCHTABLES_ENNY_CELIAN_JULIE\DeepUM_NatMeth2026V2'
# Raw KS root the original UnitMatch.mat files were generated from (same as
# run_deepunitmatch_batch.py's BASE_INPUT). Needed to correctly map a session's
# position in the full KS_dirs list to its RecSes number in MatchTable.csv --
# see original_index_to_recses() below.
RAW_KS_BASE = r'\\znas.cortexlab.net\Lab\Share\UNITMATCHTABLES_ENNY_CELIAN_JULIE\FullAnimal_KSChanMap'
MAX_C_RATIO = 1.0

# Written into a session's target folder once it has been fully copied/merged.
# Lets re-runs skip sessions that are already done instead of re-copying
# waveforms and redoing the merge over the network every time.
MERGE_COMPLETE_MARKER = 'merge_complete.flag'
REDO = False  # if True, reprocess sessions even if already marked complete


def _decode_hdf5_str(f, ref_or_ds):
    """Decode a MATLAB char array stored as uint16 in an HDF5 file."""
    obj = f[ref_or_ds] if isinstance(ref_or_ds, h5py.Reference) else ref_or_ds
    chars = obj[()].flatten()
    return ''.join(chr(int(c)) for c in chars)


def _load_uid_conversion_scipy(mat_path):
    """Load via scipy (MATLAB < v7.3). Returns (orig_clus_id, recsesAll, good_id)."""
    mat = scipy.io.loadmat(mat_path, simplify_cells=True)
    uid = mat['UniqueIDConversion']
    orig_clus_id = np.array(uid['OriginalClusID']).flatten()
    recsesAll    = np.array(uid['recsesAll']).flatten()
    good_id      = np.array(uid['GoodID']).flatten().astype(bool)
    return orig_clus_id, recsesAll, good_id


def _load_uid_conversion_hdf5(mat_path):
    """Load via h5py (MATLAB v7.3 HDF5). Returns (orig_clus_id, recsesAll, good_id)."""
    with h5py.File(mat_path, 'r') as f:
        uid = f['UniqueIDConversion']
        orig_clus_id = uid['OriginalClusID'][()].flatten()
        recsesAll    = uid['recsesAll'][()].flatten()
        good_id      = uid['GoodID'][()].flatten().astype(bool)
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
    good_session_idx = [i for i in range(n_ks_dirs) if ((recsesAll == (i + 1)) & good_id).sum() > 0]
    return {orig_idx: recses for recses, orig_idx in enumerate(good_session_idx, start=1)}


def write_synthetic_bc_unit_type_tsv(orig_clus_id, recsesAll, good_id, original_session_idx, out_path):
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
    df = pd.DataFrame({
        'cluster_id': orig_clus_id[mask].astype(int),
        'bc_unitType': np.where(good_id[mask], 'GOOD', 'MUA'),
    })
    df.to_csv(out_path, sep='\t', index=False)

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

def merge_waveforms(id1, id2, weight1, weight2, target_waveform_dir, plot_waveforms=False):

    wave1 = os.path.join(target_waveform_dir, f"Unit{str(int(id1))}_RawSpikes.npy")
    wave2 = os.path.join(target_waveform_dir, f"Unit{str(int(id2))}_RawSpikes.npy")

    # Back up the waveform files
    shutil.copy(wave1, wave1.replace("_RawSpikes.npy", "_RawSpikes_backup.npy"))
    shutil.copy(wave2, wave2.replace("_RawSpikes.npy", "_RawSpikes_backup.npy"))

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

    print(f"Merged waveforms of units {id1} and {id2} into unit {id1}, and deleted unit {id2}.")

    if plot_waveforms:
        _, ax = plt.subplots(1, 3, figsize=(10, 4))
        
        ax[0].imshow(np.nanmean(waveform1, axis=2), aspect="auto", cmap="viridis")
        ax[1].imshow(np.nanmean(waveform2, axis=2), aspect="auto", cmap="viridis")
        ax[2].imshow(np.nanmean(new_waveform, axis=2), aspect="auto", cmap="viridis")
        ax[2].set_title(f"Merged Waveform of Unit {id1}")
        ax[2].set_xlabel("Time (samples)")
        ax[2].set_ylabel("Channels")

def revert_waveform_merges(target_waveform_dir):
    # Find all backup waveform files in the target directory
    backup_files = [f for f in os.listdir(target_waveform_dir) if f.endswith("_RawSpikes_backup.npy")]

    for backup_file in backup_files:
        original_file = backup_file.replace("_RawSpikes_backup.npy", "_RawSpikes.npy")
        backup_path = os.path.join(target_waveform_dir, backup_file)
        original_path = os.path.join(target_waveform_dir, original_file)

        # Restore the original waveform file from the backup
        shutil.copy(backup_path, original_path)
        print(f"Restored {original_file} from backup.")
        
        os.remove(backup_path)  # Optionally remove the backup file after restoration

    print("All merges have been reverted.")


def run_merging_process(UMparam_files, source_dirs, MAX_C_RATIO=1.0):

    for UMparam_file, source_dir in zip(UMparam_files, source_dirs):
        with open(UMparam_file, 'rb') as file:
            data = pickle.load(file)

        target_dir = source_dir.replace(
            r"DeepUM_NatMeth2026V2",
            r"DeepUM_NatMeth2026V2_merged\merged_data_v2"
        )

        target_KSDirs = [os.path.join(target_dir, str(idx)) for idx in range(len(data['KS_dirs']))]
        already_done = [
            (not REDO) and os.path.isfile(os.path.join(d, MERGE_COMPLETE_MARKER))
            for d in target_KSDirs
        ]
        if all(already_done):
            print(f"Skipping UMparam_file: {UMparam_file} (all {len(target_KSDirs)} session(s) already processed).")
            continue

        print(f"Processing UMparam_file: {UMparam_file}...")

        # could load both UM and DUM matchtables here to get the merged units, but for now just use the UM matchtable to find which units to merge
        mt = pd.read_csv(os.path.join(os.path.split(UMparam_file)[0], 'MatchTable.csv'))
        merged_units_idx = (mt['RecSes 1'] == mt['RecSes 2']) & (mt['ID1'] != mt['ID2']) & (mt['UID 1'] == mt['UID 2'])
        # merged_units_idx = (mt['RecSes 1'] == mt['RecSes 2']) & (mt['ID1'] != mt['ID2']) & (mt['UM Probabilities'] > 0.5)

        # source_dir is DUM_NONMERGED_DATAPATH/<X>/DeepUnitMatch; the matching
        # original UnitMatch.mat lives at RAW_KS_BASE/<X>/UnitMatch/UnitMatch.mat.
        x = os.path.dirname(os.path.relpath(source_dir, DUM_NONMERGED_DATAPATH))
        mat_path = os.path.join(RAW_KS_BASE, x, 'UnitMatch', 'UnitMatch.mat')
        try:
            orig_clus_id, recsesAll, good_id = load_uid_conversion(mat_path)
            idx_to_recses = original_index_to_recses(recsesAll, good_id, len(data['KS_dirs']))
        except Exception as e:
            print(f"  WARNING: could not read {mat_path} ({e}); skipping {source_dir}.")
            continue

        for idx, KSDir in enumerate(data['KS_dirs']):
            target_KSDir = target_KSDirs[idx]
            marker_path = os.path.join(target_KSDir, MERGE_COMPLETE_MARKER)

            if already_done[idx]:
                print(f"Processing KSDir: {KSDir} (idx {idx})... already done, skipping.")
                continue

            print(f"Processing KSDir: {KSDir} (idx {idx})...")

            if not os.path.exists(target_KSDir):
                os.makedirs(target_KSDir)

            # Move KSDir files that are necessary but won't be modified to the new target directory
            for file_name in ['spike_times.npy', 'channel_positions.npy', 'cluster_info.tsv']:
                source_file = os.path.join(KSDir, file_name)
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
                orig_clus_id, recsesAll, good_id, idx,
                os.path.join(target_KSDir, 'cluster_bc_unitType.tsv'))

            source_waveform_dir = next(
                (str(path) for path in Path(KSDir).rglob("RawWaveforms") if path.is_dir()),
                None,
            )   
            target_waveform_dir = os.path.join(target_KSDir, 'qMetrics','RawWaveforms')

            spk_clusters = np.load(os.path.join(KSDir, 'spike_clusters.npy'))
            spk_times = np.load(os.path.join(KSDir, 'spike_times.npy'))
            spk_times = spk_times / 30000  # convert to seconds

            # Copy the RawWaveforms directory to the new target directory
            if not os.path.exists(target_waveform_dir):
                os.makedirs(target_waveform_dir)

            for file in Path(source_waveform_dir).iterdir():
                if file.is_file():
                    shutil.copy2(file, os.path.join(target_waveform_dir, file.name))  # preserves metadata

            # Find which units need to be merged based on the matching results.
            # RecSes in MatchTable.csv is the *compacted* session count (sessions
            # with zero good units dropped), so idx+1 is only correct when no
            # earlier session had zero good units -- use the recovered mapping
            # instead of assuming idx+1 == RecSes.
            recses = idx_to_recses.get(idx)
            if recses is None:
                # This session had zero good units in the original comparison,
                # so it can't have any self-merges to apply.
                sess_idx = pd.Series(False, index=mt.index)
            else:
                sess_idx = (mt['RecSes 1'] == recses) & (mt['RecSes 2'] == recses)
            uids = mt['UID 1'][merged_units_idx & sess_idx].unique()
            for uid in uids:
                # Get the indices of the units to be merged
                units_to_merge = mt[(mt['UID 1'] == uid) & (merged_units_idx & sess_idx)]
                unit_indices = np.unique(np.concatenate([units_to_merge['ID1'].values, units_to_merge['ID2'].values]))

                FLAG = True
                while (len(unit_indices) > 1) & FLAG:

                    # Compute C ratio for all pairs of units to be merged
                    ratios = []
                    for unit_id1 in unit_indices:
                        times1 = spk_times[np.where(spk_clusters == unit_id1)]
                        for unit_id2 in unit_indices:
                            if unit_id1 < unit_id2:  # Avoid duplicate pairs and self-comparison
                                times2 = spk_times[np.where(spk_clusters == unit_id2)]

                                if len(times1) > len(times2):
                                    non_merged = times1
                                else:
                                    non_merged = times2
                                merged = np.concatenate([times1, times2])
                                merged = np.sort(merged)

                                C = estimate_C(merged)
                                C_pm = estimate_C(non_merged)
                                ratios.append({
                                    'unit_id1': unit_id1,
                                    'unit_id2': unit_id2,
                                    'C_ratio': C / C_pm
                                })
                    ratios = pd.DataFrame(ratios)
                    ratios = ratios.sort_values('C_ratio')

                    if ratios.C_ratio[0] > MAX_C_RATIO:
                        print(f"Merging units {ratios.unit_id1[0]} and {ratios.unit_id2[0]} with C ratio {ratios.C_ratio[0]:.4f}")
                        new_unit_id = min(ratios.unit_id1[0], ratios.unit_id2[0])
                        unit_id_to_merge = max(ratios.unit_id1[0], ratios.unit_id2[0])

                        # Merge the spike_clusters
                        spk_clusters[spk_clusters == unit_id_to_merge] = new_unit_id

                        # Merge the waveforms for the merged units
                        fr_1 = len(spk_times[spk_clusters == new_unit_id])
                        fr_2 = len(spk_times[spk_clusters == unit_id_to_merge])
                        weight_1 = fr_1 / (fr_1 + fr_2)
                        weight_2 = fr_2 / (fr_1 + fr_2)
                        merge_waveforms(new_unit_id, unit_id_to_merge, weight_1, weight_2, target_waveform_dir, plot_waveforms=False)

                        unit_indices = unit_indices[unit_indices != unit_id_to_merge]
                        
                    else:
                        
                        print(f"No more units to merge based on C ratio threshold (C ratio = {ratios.C_ratio[0]:.4f}).")
                        FLAG = False

            # Save the updated spike_clusters array to the new target directory
            target_spike_clusters_file = os.path.join(target_KSDir, 'spike_clusters.npy')
            np.save(target_spike_clusters_file, spk_clusters)

            # Mark this session complete so a re-run can skip it. If anything
            # above raised, execution never reaches here, so the session stays
            # unmarked and gets fully reprocessed next time.
            with open(marker_path, 'w') as f:
                f.write('ok')

        print(f"Processing UMparam_file: {UMparam_file} done.")

def main():
    # find all UMparam.pickle from DeepUnitMatch folders, to copy their original files to new directory
    UMparam_files = []
    source_dirs = []
    for root, _, files in os.walk(DUM_NONMERGED_DATAPATH):
        for item in files:
            if item.endswith("UMparam.pickle") & ('DeepUnitMatch' in root):
                UMparam_files.append(os.path.join(root, item))
                source_dirs.append(root)

                # need to loop over DUM outputs
    run_merging_process(UMparam_files, source_dirs, MAX_C_RATIO=MAX_C_RATIO)       
    
if __name__ == '__main__':
    main()