import numpy as np
from scipy.signal import detrend
import h5py, os
from pathlib import Path
from utils.myutil import get_default_param

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

def sort_good_channels(goodChannelMap, goodpos):
    '''
    Sorts the good channels by their y-axis value and then by their z-axis value.
    '''
    # Step 1: Identify the unique y-axis values and sort them
    unique_y_values = np.unique(goodpos[:, 0])
    unique_y_values.sort()

    # Safety check: ensure there are exactly two unique y-axis values
    if len(unique_y_values) != 2:
        print(unique_y_values)
        print(f"Found {len(unique_y_values)} instead of 2 unique x values.")
        # TODO: adapt this code to be robust to Neuropixels 1.0 recordings as well as 2.0.
        # For now we are only using Neuropixels 2.0 data -> if we enter this block it means there was a mistake in spike sorting
        # Therefore the fix for now is to use return 0s so that we can handle this in extract_Rwaveforms

        # print(f"Channel Map: {goodChannelMap}")
        # print(f"Pos: {goodpos}")
        # raise ValueError(f"There should be exactly two unique y-axis values for Neuropixels 2.0 shank - instead got {len(unique_y_values)}: [{unique_y_values}]")
        return [-1],[-1]

    channels_y_min_indices = np.where(goodpos[:, 0] == unique_y_values[0])[0]
    channels_y_max_indices = np.where(goodpos[:, 0] == unique_y_values[1])[0]

    channels_y_min = goodChannelMap[channels_y_min_indices]
    channels_y_max = goodChannelMap[channels_y_max_indices]

    pos_y_min = goodpos[channels_y_min_indices]
    pos_y_max = goodpos[channels_y_max_indices]
    
    # Step 3: Sort each group by the z-axis value
    z_min_sorted_indices = np.argsort(goodpos[goodpos[:, 0] == unique_y_values[0], 1])
    z_max_sorted_indices = np.argsort(goodpos[goodpos[:, 0] == unique_y_values[1], 1])

    channels_y_min_sorted = channels_y_min[z_min_sorted_indices]
    channels_y_max_sorted = channels_y_max[z_max_sorted_indices]

    pos_y_min_sorted = pos_y_min[z_min_sorted_indices]
    pos_y_max_sorted = pos_y_max[z_max_sorted_indices]

    # Step 4: Interleave the channels from the two groups
    sorted_goodChannelMap = np.empty_like(goodChannelMap)
    sorted_goodChannelMap[::2] = channels_y_min_sorted[:len(sorted_goodChannelMap)//2]  # Even indices
    sorted_goodChannelMap[1::2] = channels_y_max_sorted[:len(sorted_goodChannelMap)//2]  # Odd indices
    sorted_goodpos = np.empty_like(goodpos)
    sorted_goodpos[::2, :] = pos_y_min_sorted
    sorted_goodpos[1::2, :] = pos_y_max_sorted
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
    sorted_goodChannelMap,sorted_goodpos = sort_good_channels(goodChannelMap, goodpos)
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

def save_waveforms_hdf5(file_name, Rwaveform, MaxSitepos, session):
    """
    Saves the preprocessed, reduced waveform and the max site position as a HDF5 file.
    """
    current_dir = Path(__file__).parent.parent
    dest_path = os.path.join(current_dir, "processed_waveforms", str(session), file_name)
    dest_directory = os.path.dirname(dest_path)
    os.makedirs(dest_directory, exist_ok=True)

    new_data = {
        "waveform": Rwaveform,          # (60,30,2) 
        "MaxSitepos": MaxSitepos
    }
    with h5py.File(dest_path, 'w') as f:
        for key, value in new_data.items():
            f.create_dataset(key, data=value)

def get_snippets(waveforms, ChannelPos, session_id):
    """
    Convert raw waveforms from a pair of sessions (waveforms) to snippets with shape (60,30,2).

    Arguments:
    - waveforms: n_units x n_timepoints x n_channels x n_repeats (CV)
    - ChannelPos: list of channel positions for each session. Each session's channel positions should be a 384 x 3 or 384 x 2 array.
    """
    processed_waveforms = []
    positions = []
    params = get_default_param()
    ChannelMap = np.arange(384).astype(np.int32)

    if ChannelPos[0].shape[1] == 3:
        # Just need 2D probe positions
        ChannelPos = ChannelPos[0][:, 1:]
    else:
        # Just need 1 array as these are assumed to be the same for all repeats (CV)
        ChannelPos = ChannelPos[0]

    for i, unit in enumerate(waveforms):
        n_time, n_channels, n_repeats = unit.shape
        if n_time != 82 or n_channels != 384 or n_repeats != 2:
            print(f"Expected waveform shape to be (82,384,2) - instead got {unit.shape}")
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
        
        save_waveforms_hdf5(f'Unit{i}.npy', Rwaveform, MaxSitepos, session_id[i])
    
        processed_waveforms.append(Rwaveform)
        positions.append(MaxSitepos)

        if Rwaveform.shape!=(60,30,2):
            print(f"Incorrect shape for {i}th unit: {Rwaveform.shape}")

    return np.array(processed_waveforms), np.array(positions)
