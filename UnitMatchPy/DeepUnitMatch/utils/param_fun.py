import numpy as np
from scipy.signal import detrend
import h5py, os
import shutil
from pathlib import Path

# Track which processed_waveforms session directories we've already cleaned this run
_CLEANED_WAVEFORM_DIRS = set()

def _resolve_processed_root(save_path=None) -> Path:
    """
    Return the path to the temporary output folder "processed_waveforms".

    Note: save_path is treated as a *base* directory; we only create/remove the
    "processed_waveforms" subfolder within it (never save_path itself).
    """
    # Historical default: store outputs under the DeepUnitMatch package directory.
    base_path = Path(__file__).parent.parent if save_path is None else Path(save_path)

    # If the user mistakenly passes a path that already ends with "processed_waveforms",
    # normalize to its parent so we still only operate on ".../<base>/processed_waveforms".
    if base_path.name == "processed_waveforms":
        base_path = base_path.parent

    return base_path / "processed_waveforms"

def get_default_param(param = None):
    """
    Create param, a dictionary with the default parameters.
    If a dictionary is given, it will add values to it without overwriting existing values.
    Do not need to give a dictionary.
    """
    tmp = {'nTime' : 82, 'nChannels' : 384, 'ChannelRadius' : 110,
           'RnChannels' : 30, 'RnTime' : 60,
           'n_xchannelpos' : 2,  # number of unique x-axis (column) positions on the probe shank (2 for NP2.0, 4 for NP1.0)
        }
    # if no dictionary is given just returns the default parameters
    if param is None:
        out = tmp
    else:
        # Add default parameters to param dictionary, does not overwrite pre existing param values
        out = tmp | param
        # spike_width is an older alias for nTime; ensure nTime follows it when available
        if 'spike_width' in out:
            out['nTime'] = out['spike_width']
    if out['RnChannels'] % 2 != 0:
        print('RnChannels is not even, please check')
    return out

def detrend_waveform(waveform):
    """
    This function accepts the raw waveforms (nTime, nChannels, CV)
    The output is the same shape, and linearly detrended across time.
    """ 
    return detrend(waveform, axis = 1)

def get_spatialfp(waveform):
    """
    Input: waveform np array (nTime, nChannels, first/second half)
    By taking the maximum value along the time axis of the absolute values
    Output: np array (nChannels, first/second half)
    spatialfp = spatial footprint
    """
    SpatialFP = np.max(np.abs(waveform), axis = 0)
    return SpatialFP

def get_max_site(SpatialFP):
    """
    Input: SpatialFP (nChannels, first/second half)
    By taking the index of maximum argument along the channel axis,
    Output: (first/second half), which gives the maximum site for each unit in first/second half
    """
    MaxSite = np.argmax(SpatialFP, axis = 0)
    return MaxSite

def sort_good_channels(goodChannelMap, goodpos, n_xchannelpos=2):
    '''
    Sorts good channels by depth (z-axis), interleaved across x-axis column positions.

    Works for any number of probe columns (n_xchannelpos). Channels are grouped by
    their x-axis position, each group is sorted by depth, then all groups are
    interleaved round-robin (largest group first to fill the leading positions).
    '''
    unique_x_values = np.unique(goodpos[:, 0])
    unique_x_values.sort()

    n_unique_x = len(unique_x_values)
    # Multi-shank probes have large x-gaps (>50 µm) between shanks; each shank
    # contributes n_xchannelpos columns, so the total must be a positive multiple.
    x_gaps = np.diff(unique_x_values)
    n_shanks = int(np.sum(x_gaps > 50)) + 1
    if n_unique_x != n_xchannelpos * n_shanks:
        print(f"Found {n_unique_x} unique x-axis positions across {n_shanks} shank(s); "
              f"expected a multiple of {n_xchannelpos} ({n_xchannelpos * n_shanks}).")
        return [-1], [-1]

    # Group channels by x-column and sort each group by depth (z)
    groups_ch = []
    groups_pos = []
    for x_val in unique_x_values:
        idx = np.where(goodpos[:, 0] == x_val)[0]
        z_order = np.argsort(goodpos[idx, 1])
        groups_ch.append(goodChannelMap[idx][z_order])
        groups_pos.append(goodpos[idx][z_order])

    # Put the largest group first so it occupies the lowest-index slots
    # (preserves the original NP2.0 behaviour where the bigger column gets even positions)
    col_order = np.argsort([len(g) for g in groups_ch])[::-1]
    groups_ch  = [groups_ch[i]  for i in col_order]
    groups_pos = [groups_pos[i] for i in col_order]

    # Round-robin interleave across columns, stepping through depth together
    sorted_goodChannelMap = np.empty_like(goodChannelMap)
    sorted_goodpos = np.empty_like(goodpos)
    pos = 0
    for depth_idx in range(max(len(g) for g in groups_ch)):
        for g_ch, g_pos in zip(groups_ch, groups_pos):
            if depth_idx < len(g_ch):
                sorted_goodChannelMap[pos] = g_ch[depth_idx]
                sorted_goodpos[pos] = g_pos[depth_idx]
                pos += 1

    return sorted_goodChannelMap, sorted_goodpos

def extract_Rwaveforms(waveform, ChannelPos,ChannelMap, param):
    """
    Using waveforms, ChannelPos and param, to find the max channel for each unit and cv, this function also
    returns good idx's / positions, by selecting channels within ChannelRadius (default 150 um) 
    Input: waveform (nTime, nChannels, CV), ChannelPos (nChannels, 2), param (dictionary)
    """

    nChannels = param['nChannels']
    nTime = param['nTime']
    RnChannels = param['RnChannels']
    RnTime = param['RnTime']
    ChannelRadius = param['ChannelRadius']
    # original time 0-82, new time 11-71
    start_time,end_time = (nTime - RnTime) // 2, (nTime + RnTime) // 2
    if waveform.ndim==2:
        waveform = np.stack([waveform, waveform], axis=2)
    waveform = waveform[start_time:end_time,:,:] # selecting the middle 60 time points (11-71)
    waveform = detrend_waveform(waveform) # detrend the waveform
    MeanCV = np.mean(waveform, axis = 2) # average of each cv
    SpatialFootprint = get_spatialfp(MeanCV) # choose max time 
    MaxSiteMean = get_max_site(SpatialFootprint) # argument of MaxSite
    MaxSitepos = ChannelPos[MaxSiteMean,:] #gives the 2-d positions of the max sites

    # Finds the indices where the distance from the max site mean is small
    goodidx = np.empty(nChannels, dtype=bool)
    for i in range(ChannelPos.shape[0]): #looping over each site
        dist = np.linalg.norm(ChannelPos[MaxSiteMean,:] - ChannelPos[i,:])
        good = dist < ChannelRadius
        goodidx[i] = good

    goodChannelMap = ChannelMap[goodidx] #selecting the good channels
    goodpos = ChannelPos * np.tile(goodidx, (2,1)).T
    goodpos = goodpos[goodidx,:]
    sorted_goodChannelMap, sorted_goodpos = sort_good_channels(
        goodChannelMap, goodpos, n_xchannelpos=param.get('n_xchannelpos', 2)
    )
    if sorted_goodChannelMap[0]==-1 and sorted_goodpos[0]==-1:
        # we have a spike sorting error so need to 0 this recording
        return np.array([-1,-1]), np.array([-1,-1]), [0], [0], np.zeros((1,1,1))
    Rwaveform = waveform[:, sorted_goodChannelMap, :] #selecting the good channels
    
    ## this part is tricks to make the data proper for DNN training
    GlobalMean = np.mean(Rwaveform) # mean of all channels and time points
    Rwaveform = Rwaveform - GlobalMean # subtracting the global mean, zero mean is good for DNN
    # padding the data to make it proper for DNN
    NewGlobalMean = np.mean(Rwaveform)
    z_sorted_goodpos = np.unique(sorted_goodpos[:,1])
    mean_z_sorted_goodpos = np.mean(z_sorted_goodpos)
    z_MaxSitepos = MaxSitepos[1]
    num_good_channels = np.sum(goodidx)
    # print('num_good_channels',num_good_channels)
    padding_needed = RnChannels - num_good_channels
    pad_before = 0
    pad_after = 0
    if z_MaxSitepos < mean_z_sorted_goodpos:
        pad_before = padding_needed  # Pad at the beginning if MaxSitepos is below the mean z position
    else:
        pad_after = padding_needed  # Pad at the end if MaxSitepos is above the mean z position
    Rwaveform = np.pad(Rwaveform, ((0, 0), (pad_before, pad_after), (0, 0)), 'constant', constant_values=(NewGlobalMean, NewGlobalMean))
    
    return MaxSiteMean, MaxSitepos, sorted_goodChannelMap, sorted_goodpos, Rwaveform

def save_waveforms_hdf5(file_name, Rwaveform, MaxSitepos, session, save_path=None):
    """
    Saves the preprocessed, reduced waveform and the max site position as a HDF5 file.
    """
    processed_root = _resolve_processed_root(save_path)
    dest_directory = processed_root / str(session)
    dest_path = dest_directory / file_name

    # If the destination directory already exists, remove existing waveform files once per session
    try:
        if dest_directory.is_dir() and dest_directory not in _CLEANED_WAVEFORM_DIRS:
            for fpath in dest_directory.iterdir():
                if fpath.is_file():
                    fpath.unlink()
        _CLEANED_WAVEFORM_DIRS.add(dest_directory)
    except Exception as e:
        print(f"Warning: could not clean destination directory {dest_directory}: {e}")

    dest_directory.mkdir(parents=True, exist_ok=True)

    new_data = {
        "waveform": Rwaveform,          # (60,30,2) 
        "MaxSitepos": MaxSitepos
    }
    with h5py.File(dest_path, 'w') as f:
        for key, value in new_data.items():
            f.create_dataset(key, data=value)

def get_snippets(waveforms, ChannelPos, session_id, save_path=None, unit_ids=None, param=None):
    """
    Convert raw waveforms from a pair of sessions (waveforms) to snippets with shape (60,30,2).

    Arguments:
    - waveforms: n_units x n_timepoints x n_channels x n_repeats (CV)
    - ChannelPos: list of channel positions for each session. Each session's channel positions should be a 384 x 3 or 384 x 2 array.
    - session_id: list of session identifiers, 1 integer per unit.
    - unit_ids: Optional original unit IDs (length n_units). If provided, snippet files are saved as
      Unit{unit_id}*_RawSpikes.npy so downstream inference can match UnitMatch ordering/labels.
    - save_path: Optional path to save the processed waveform files. If None, defaults to 'processed_waveforms' folder in current directory.
    """
    # If this save path already contains outputs from a previous run, remove the entire
    # processed_waveforms folder so we don't mix sessions. This folder is generated by
    # this code and treated as temporary.
    processed_root = _resolve_processed_root(save_path)
    if processed_root.exists():
        try:
            if processed_root.is_dir():
                shutil.rmtree(processed_root)
            else:
                processed_root.unlink()
        except Exception as e:
            print(f"Warning: could not remove temporary directory {processed_root}: {e}")

    processed_waveforms = []
    positions = []
    kept_indices = []
    params = get_default_param(param)
    ChannelMap = np.arange(384).astype(np.int32)

    if ChannelPos[0].shape[1] == 3:
        # Just need 2D probe positions
        ChannelPos = ChannelPos[0][:, 1:]
    else:
        # Just need 1 array as these are assumed to be the same for all repeats (CV)
        ChannelPos = ChannelPos[0]

    if unit_ids is not None:
        unit_ids = np.asarray(unit_ids).squeeze()
        if unit_ids.shape[0] != len(waveforms):
            raise ValueError("unit_ids must have the same length as waveforms.")

    expected_n_time = params['nTime']
    expected_n_channels = params['nChannels']
    for i, unit in enumerate(waveforms):
        n_time, n_channels, n_repeats = unit.shape
        if n_time != expected_n_time or n_channels != expected_n_channels or n_repeats != 2:
            print(f"Expected waveform shape to be ({expected_n_time},{expected_n_channels},2) - instead got {unit.shape}")
            continue
        if np.any(np.isnan(unit)) or np.any(np.isinf(unit)):
            if not np.any(np.isnan(unit[:,:,0])) and not np.any(np.isinf(unit[:,:,0])):
                unit = unit[:,:,0]
                MaxSiteMean, MaxSitepos, sorted_goodChannelMap, sorted_goodpos, Rwaveform = extract_Rwaveforms(unit, ChannelPos, ChannelMap, params)
            elif not np.any(np.isnan(unit[:,:,1])) and not np.any(np.isinf(unit[:,:,1])):
                unit = unit[:,:,1]
                MaxSiteMean, MaxSitepos, sorted_goodChannelMap, sorted_goodpos, Rwaveform = extract_Rwaveforms(unit, ChannelPos, ChannelMap, params)
            else:
                print(f"Corrupted data where unexpected in {i}th unit")
                continue
        else:
            # Clean data - go ahead as normal
            MaxSiteMean, MaxSitepos, sorted_goodChannelMap, sorted_goodpos, Rwaveform = extract_Rwaveforms(unit, ChannelPos, ChannelMap, params)

        if Rwaveform.shape != (params['RnTime'], params['RnChannels'], 2):
            print(f"Incorrect shape for {i}th unit: {Rwaveform.shape}")
            continue

        unit_id = int(unit_ids[i]) if unit_ids is not None else i
        save_waveforms_hdf5(f"Unit{unit_id}_RawSpikes.npy", Rwaveform, MaxSitepos, session_id[i], save_path=save_path)

        processed_waveforms.append(Rwaveform)
        positions.append(MaxSitepos)
        kept_indices.append(i)

    return np.array(processed_waveforms), np.array(positions), np.array(kept_indices, dtype=int)
