# utility function for loading files etc
import numpy as np
import pandas as pd
import os
import matplotlib.pyplot as plt
import mat73
import sqlite3
import time
import random

def load_tsv(path):
    """
    Loadsa .tsv file as a numpy array, with the headers removed

    Parameters
    ----------
    path : str
        The path to the tsv to load

    Returns
    -------
    ndarray
        The tsv as a ndarray
    """
    df  = pd.read_csv(path, sep='\t', skiprows = 0)
    return df.values

def get_session_number(unit_id, session_switch):
    """
    Finds the session number of a unit given its id and the session_switch array

    Parameters
    ----------
    unit_id : int
        The UnitMatch unit id
    session_switch : ndarray
        A array which marks at which units the a new session starts

    Returns
    -------
    int
        The session number of the unit
    """
    for i in range(len(session_switch) - 1):
        if (session_switch[i] <= unit_id < session_switch[i+1]):
            return i

def get_session_data(n_units_per_session):
    """
    Calculates information on the sessions using the number of units per session

    Parameters
    ----------
    n_units_per_session : ndarray
        An array where each value is how many units appeared in the session

    Returns
    -------
    ndarrays
        The calculated session information
    """
    n_sessions = len(n_units_per_session)
    #Total number of units                  
    n_units = n_units_per_session.sum()

    sessionid = np.zeros(n_units, dtype = int)
    #What units the a new session starts
    session_switch = np.cumsum(n_units_per_session)
    session_switch = np.insert(session_switch, 0, 0)
    for i in range(n_sessions):
        #The session id for each unit
        sessionid[session_switch[i]:session_switch[i+1]] = int(i)

    return n_units, sessionid, session_switch, n_sessions

def get_within_session(session_id, param):
    """
    Creates an array with 1 if the units are in the same session and a 0 otherwise

    Parameters
    ----------
    session_id : ndarray
        The session id for each unit
    param : dict
        the param dictionary

    Returns
    -------
    ndarray
        A n_unit * n_unit array which marks units in the same session
    """
    n_units = param['n_units']

    tmp1 = np.expand_dims(session_id , axis=1)
    tmp2 = np.expand_dims(session_id, axis=0)

    within_session = np.ones((n_units, n_units))
    within_session[tmp1 == tmp2] = 0

    return within_session

def load_good_waveforms(wave_paths, unit_label_paths, param, good_units_only=True):
    if len(wave_paths) == len(unit_label_paths):
        n_sessions = len(wave_paths)
    else:
        print('Warning: gave different number of paths for waveforms and labels!')
        return

    good_units = []
    n_units_per_session_all = []
    all_units = []

    if good_units_only:
        for i in range(len(unit_label_paths)):
            if os.path.split(unit_label_paths[0])[1] == 'cluster_bc_unitType.tsv':
                unit_label = load_tsv(unit_label_paths[i])
                tmp_idx = np.argwhere(np.isin(unit_label[:, 1], ['GOOD', 'NON-SOMA GOOD']))
            else:
                unit_label = load_tsv(unit_label_paths[i])
                tmp_idx = np.argwhere(unit_label[:, 1]=='good')

            n_units_per_session_all.append(unit_label.shape[0])
            good_unit_idx = unit_label[tmp_idx, 0]
            good_units.append(good_unit_idx)
            all_units.append(unit_label[:, 0])

    waveforms = []
    if good_units_only:
        for ls in range(len(wave_paths)):
            try:
                p_file = os.path.join(wave_paths[ls], f'Unit{int(good_units[ls][0].squeeze())}_RawSpikes.npy')
                tmp = np.load(p_file)
                tmp_waveform = np.zeros((len(good_units[ls]), tmp.shape[0], tmp.shape[1], tmp.shape[2]))

                for i in range(len(good_units[ls])):
                    p_file_good = os.path.join(wave_paths[ls], f'Unit{int(good_units[ls][i].squeeze())}_RawSpikes.npy')
                    tmp_waveform[i] = np.load(p_file_good)
                waveforms.append(tmp_waveform)
            except Exception as e:
                print(f'Error loading waveform for session {ls}: {e}')
            finally:
                del tmp_waveform
                del tmp
    else:
        for ls in range(len(wave_paths)):
            try:
                p_file = os.path.join(wave_paths[ls], f'Unit0_RawSpikes.npy')
                tmp = np.load(p_file)
                tmp_waveform = np.zeros((len(os.listdir(wave_paths[ls])), tmp.shape[0], tmp.shape[1], tmp.shape[2]))

                for i in range(len(os.listdir(wave_paths[ls]))):
                    p_file_good = os.path.join(wave_paths[ls], f'Unit{i}_RawSpikes.npy')
                    tmp_waveform[i] = np.load(p_file_good)
                waveforms.append(tmp_waveform)
                print(f'UnitMatch is treating all the units as good and including all units from {wave_paths[ls]}, we recommended using curated data!')
            except Exception as e:
                print(f'Error loading waveform for session {ls}: {e}')
            finally:
                del tmp_waveform
                del tmp

    n_units_per_session = np.zeros(n_sessions, dtype='int')
    waveform = np.array([])

    for i in range(n_sessions):
        if i == 0:
            waveform = waveforms[i]
        else:
            waveform = np.concatenate((waveform, waveforms[i]), axis=0)

        n_units_per_session[i] = waveforms[i].shape[0]

    param['n_units'], session_id, session_switch, param['n_sessions'] = get_session_data(n_units_per_session)
    within_session = get_within_session(session_id, param)
    param['n_channels'] = waveform.shape[2]
    param['n_units_per_session'] = n_units_per_session_all

    if param['spike_width'] != waveform.shape[1]:
        param['spike_width'] = waveform.shape[1]
        param['peak_loc'] = np.floor(waveform.shape[1] / 2).astype(int)
        param['waveidx'] = np.arange(param['peak_loc'] - 8, param['peak_loc'] + 15, dtype=int)

    return waveform, session_id, session_switch, within_session, good_units, param

def get_good_units(unit_label_paths, good = True):
    """
    This function is used if you want to find good units then load them in
    (first half of load_good_waveforms)

    Parameters
    ----------
    unit_label_paths : list
        A list were each entry is a path to either BombCell good units (cluster_bc_unitType.tsv)
        or the KiloSort good units (cluster_group.tsv') for each session
    good : bool, optional
        If True will only load in units marked good
        If False will load all units labeled in the given .tsv, by default True

    Returns
    -------
    ndarray
        A list of all the good unit ids
    """
    good_units = []
    for i in range(len(unit_label_paths)):
    #see if bombcell unit labels
        if os.path.split(unit_label_paths[0])[1] == 'cluster_bc_unitType.tsv':
            unit_label = load_tsv(unit_label_paths[i])
            if good == True:
                tmp_idx = np.argwhere(np.isin(unit_label[:,1],['GOOD','NON-SOMA GOOD']))
            else:
                tmp_idx = unit_label[:,0].astype(np.int32)
            tmp_idx = tmp_idx[:, np.newaxis] # keep the array shape consistent between different methods
        else:
            unit_label = load_tsv(unit_label_paths[i])
            if good == True:
                tmp_idx = np.argwhere(unit_label[:,1] == 'good')
            else:
                tmp_idx = unit_label[:,0].astype(np.int32) # every unit index in the first column
                
        good_unit_idx = unit_label[tmp_idx, 0]
        good_units.append(good_unit_idx)
    return good_units

def load_good_units(good_units, wave_paths, param):
    if len(wave_paths) == len(good_units):
        n_sessions = len(wave_paths)
    else:
        print('Warning: gave different number of paths for waveforms and labels!')
        return

    waveforms = []
    for ls in range(len(wave_paths)):
        try:
            p_file = os.path.join(wave_paths[ls], f'Unit{int(good_units[ls][0].squeeze())}_RawSpikes.npy')
            tmp = np.load(p_file)
            tmp_waveform = np.zeros((len(good_units[ls]), tmp.shape[0], tmp.shape[1], tmp.shape[2]))

            for i in range(len(good_units[ls])):
                tmp_path_good = os.path.join(wave_paths[ls], f'Unit{int(good_units[ls][i].squeeze())}_RawSpikes.npy')
                tmp_waveform[i] = np.load(tmp_path_good)
            waveforms.append(tmp_waveform)
        except Exception as e:
            print(f'Error loading waveform for session {ls}: {e}')
        finally:
            del tmp_waveform
            del tmp

    n_units_per_session = np.zeros(n_sessions, dtype='int')
    waveform = np.array([])

    for i in range(n_sessions):
        if i == 0:
            waveform = waveforms[i]
        else:
            waveform = np.concatenate((waveform, waveforms[i]), axis=0)

        n_units_per_session[i] = waveforms[i].shape[0]

    param['n_units'], session_id, session_switch, param['n_sessions'] = get_session_data(n_units_per_session)
    within_session = get_within_session(session_id, param)
    param['n_channels'] = waveform.shape[2]
    param['n_units_per_session'] = n_units_per_session

    return waveform, session_id, session_switch, within_session, param

def evaluate_output(output_prob, param, within_session, session_switch, match_threshold = 0.5):
    """
    This function evaluates summary values for the UnitMatch results by finding:
    The number of units matched to themselves across cv
    The false negative %, how many did not match to themselves across cv
    the false positive % in two ways, how many miss-matches are there in the off-diagonal per session
    and how many  false match out of how many matches we should get

    Parameters
    ----------
    output_prob : ndarray (n_units, n_units)
        The output match probability array
    param : dict
        The param dictionary
    within_session : ndarray
        The array which marks units pairs in the same session
    session_switch : ndarray
        The array which marks when a new session starts
    match_threshold : float, optional
        The threshold value which decides matches, by default 0.5
    """

    output_threshold = np.zeros_like(output_prob)
    output_threshold[output_prob > match_threshold] = 1

    # get the number of diagonal matches
    n_diag = np.sum(output_threshold[np.eye(param['n_units']).astype(bool)])
    self_match = n_diag / param['n_units'] *100
    print(f'The percentage of units matched to themselves is: {self_match:.2f}%')
    print(f'The percentage of false -ve\'s then is: {100 - self_match:.2f}% \n')

    #off-diagonal miss-matches
    n_off_diag = np.zeros_like(output_prob)
    n_off_diag = output_threshold
    n_off_diag[within_session == 1] = 0 
    n_off_diag[np.eye(param['n_units']) == 1] = 0 
    false_positive_est =  n_off_diag.sum() / (param['n_units']) 
    print(f'The rate of miss-match(es) per expected match {false_positive_est:.2f}')


    #compute matlab FP per session per session
    false_positive_est_per_session = np.zeros(param['n_sessions'])
    for did in range(param['n_sessions']):
        tmp_diag = output_threshold[session_switch[did]:session_switch[did + 1], session_switch[did]:session_switch[did + 1]]
        n_units = tmp_diag.shape[0]
        tmp_diag[np.eye(n_units) == 1] = 0 
        false_positive_est_per_session[did] = tmp_diag.sum() / (n_units ** 2 - n_units) * 100
        print(f'The percentage of false +ve\'s is {false_positive_est_per_session[did]:.2f}% for session {did +1}')

    print('\nThis assumes that the spike sorter has made no mistakes')

def curate_matches(matches_GUI, is_match, not_match, mode='and'):
    """
    There are two options, 'and' 'or'. 
    'And' gives a match if both CV give it as a match
    'Or gives a match if either CV gives it as a match

    Parameters
    ----------
    matches_GUI : ndarray or None
        The array of matches calculated for the GUI or None if not available
    is_match : list
        A list of pairs manually curated as a match in the GUI
    not_match : list
        A list of pairs manually curated as NOT a match in the GUI
    mode : str, optional
        either 'and' or 'or' depending on preferred rules of CV concatenation, by default 'and'

    Returns
    -------
    ndarray
        The curated list of matches
    """
    if matches_GUI is not None:
        matches_a = matches_GUI[0]
        matches_b = matches_GUI[1]
    else:
        matches_a = np.zeros((0, 2))
        matches_b = np.zeros((0, 2))

    # if both arrays are empty leave function
    if np.logical_and(len(is_match) == 0, len(not_match) == 0):
        print('There are no curated matches/none matches')
        return None

    # if one array is empty make it have corrected shape
    if len(is_match) == 0:
        is_match = np.zeros((0, 2))
    else:
        is_match = np.array(is_match)

    if len(not_match) == 0:
        not_match = np.zeros((0, 2))
    else:
        not_match = np.array(not_match)

    if mode == 'and':
        matches_tmp = np.concatenate((matches_a, matches_b), axis=0)
        matches_tmp, counts = np.unique(matches_tmp, return_counts=True, axis=0)
        matches = matches_tmp[counts == 2]
    elif mode == 'or':
        matches = np.unique(np.concatenate((matches_a, matches_b), axis=0), axis=0)
    else:
        print('please make mode = \'and\' or \'or\' ')
        return None

    # add matches in IS Matches
    matches = np.unique(np.concatenate((matches, is_match), axis=0), axis=0)
    print(matches.shape)
    # remove Matches in NotMatch
    matches_tmp = np.concatenate((matches, not_match), axis=0)
    matches_tmp, counts = np.unique(matches_tmp, return_counts=True, axis=0)
    matches = matches_tmp[counts == 1]

    return matches

def isin_2d(test_arr, parent_arr):
    return  (test_arr[:, None] == parent_arr).all(-1).any(-1)

def fill_missing_pos(KS_dir, n_channels):
    """
    KiloSort (especially in 4.0) may not include channel positions for inactive channels, 
    as UnitMatch require the full channel_pos array this function will extrapolate it from the given channel positions.
    If there is only one missing positions will fill that in else it will go through each found missing position
    and attempt to fill them in iteratively, finally it will attempt to find a pattern in the x-coord and fill in positions that way

    Parameters
    ----------
    KS_dir : str
        The path to the KiloSort directory
    n_channels : int
        The number of channels 

    Returns
    -------
    ndarray
        The full channel_pos array
    """
    print('The channel_positions.npy file does not match with the raw waveforms \n \
           we have attempted to fill in the missing positions, please check the attempt worked and examine the channel_positions and RawWaveforms shape')
    path_tmp = os.path.join(KS_dir, 'channel_positions.npy')
    pos = np.load(path_tmp)

    path_tmp = os.path.join(KS_dir, 'channel_map.npy')
    channel_map = np.load(path_tmp).squeeze()

    channel_pos = np.full((n_channels,2), np.nan)
    channel_pos[channel_map,:] = pos

    channel_pos_new = []
    #get the unique x positions
    x_unique = np.unique(channel_pos[:,0])
    x_unique = x_unique[~np.isnan(x_unique)]

    #go through each column
    for x in x_unique:
        #get the known y-values for that column
        y_column = channel_pos[np.argwhere(channel_pos[:,0] == x), 1].squeeze()

        #test to see if any other columns have the same set of y positions
        same_x_pattern = np.unique(channel_pos[np.in1d(channel_pos[:,1], y_column), 0])
        same_y_pattern = channel_pos[np.in1d(channel_pos[:,0], same_x_pattern), 1]

        #find the mode difference, i.e the steps between y-positions 
        y_steps, y_step_counts = np.unique(np.diff(np.unique(same_y_pattern)), return_counts= True)
        y_steps = y_steps[np.argmax(y_step_counts)].astype(int)

        #find the min/max y-positions to fill in all positions for the column
        ymin = np.min(same_y_pattern).astype(int)
        ymax = np.max(same_y_pattern).astype(int)
        ypos = np.arange(ymin, ymax+y_steps, y_steps)

        channel_pos_column = np.stack((np.full_like(ypos,x), ypos)).T
        channel_pos_new.append(channel_pos_column)

    if pos.shape[0] == n_channels - 1:
        print('Likely to be correctly filled')
        missed_pos_idx = np.argwhere(isin_2d(np.vstack(channel_pos_new), pos) == False)
        channel_pos[np.isnan(channel_pos)[:,0],:] = np.vstack(channel_pos_new)[missed_pos_idx]
        return channel_pos

    channel_pos_fill = channel_pos.copy()
    #find all the positions not in the original channel_positions
    missed_pos_idxs = np.argwhere(isin_2d(np.vstack(channel_pos_new), pos) == False)
    #find all the empty channel positions indexs
    empty_idx = np.unique(np.argwhere(np.isnan(channel_pos))[:,0])

    #check to see if there are the same amount of empty positions as found positions
    if missed_pos_idx.shape[0] == empty_idx.shape[0]:
        #go through each estimated missing positions
        for idx in missed_pos_idxs:
            missed_pos = np.vstack(channel_pos_new)[idx].squeeze()
            #find the index of the nearest channel to the estimated missing positions
            distance_to_pos = np.sqrt(np.sum((channel_pos - missed_pos)**2, axis = 1) )
            nearest_idx = np.nanargmin(distance_to_pos)
            #assume that the nearest missing index to the nearest positions is the correct index!
            estimated_missing_idx = empty_idx[np.argmin(np.abs(empty_idx - nearest_idx))]
            #fill in this index
            channel_pos_fill[estimated_missing_idx,:] = missed_pos
            #delete this index from the empty idx_list
            empty_idx = np.delete(empty_idx, np.argmin(np.abs(empty_idx - nearest_idx)))

        #If all of the original positions are correct
        if np.sum(channel_pos == channel_pos_fill) //2 == pos.shape[0]:
            print('Likely to be correctly filled')
            return channel_pos_fill

    n_unique_x = x_unique.shape[0]
    channel_pos_fill = np.zeros_like(channel_pos)
    test_points = []
    starts = []
    x_points = []
    #Often the channel positions is filled in in a periodic fashion in x coord
    for i, x in enumerate(x_unique):
        #find which sequence of positions this x-column fills 
        x_point = np.argwhere(channel_pos[:,0] == x).squeeze()[0]
        start = x_point % n_unique_x
        points = np.arange(start, n_channels, n_unique_x)
        #fill in the positions for this column
        for j, point in enumerate(points):
            channel_pos_fill[point,:] = channel_pos_new[i][j,:]

        test_points.append(points)
        starts.append(start)
        x_points.append(x_point)

    #If all of the original positions are correct
    if np.sum(channel_pos == channel_pos_fill) //2 == pos.shape[0]:
        print('Likely to be correctly filled')
        return channel_pos_fill
    else:
        print('Error in filling channel positions, please fill in manually \n \
            Have returned the known channel positions and NaN')
        return channel_pos

def paths_from_KS(KS_dirs, custom_raw_waveform_paths=None, custom_bombcell_paths=None):
    """
    This function will find specific paths to required files from a KiloSort directory
    or use custom paths if provided

    Parameters
    ----------
    KS_dirs : list
        The list of paths to the KiloSort directory for each session
    custom_raw_waveform_paths : list, optional
        Custom paths to raw waveform directories for each session. If provided,
        these will be used instead of searching within KS_dirs
    custom_bombcell_paths : list, optional
        Custom paths to BombCell/unit label files for each session. If provided,
        these will be used instead of searching within KS_dirs

    Returns
    -------
    list
        The lists to the files for each session
    """
    n_sessions = len(KS_dirs)

    #load in the number of channels
    tmp = os.getcwd()

    # Use custom raw waveform paths if provided, otherwise search in KS directories
    if custom_raw_waveform_paths is not None:
        if len(custom_raw_waveform_paths) != n_sessions:
            raise ValueError(f"Number of custom raw waveform paths ({len(custom_raw_waveform_paths)}) must match number of KS directories ({n_sessions})")
        wave_paths = custom_raw_waveform_paths
    elif custom_bombcell_paths is not None:
        if len(custom_bombcell_paths) != n_sessions:
            raise ValueError(f"Number of custom raw waveform paths ({len(custom_bombcell_paths)}) must match number of KS directories ({n_sessions})")
        wave_paths = []
        for i in range(n_sessions):
             wave_paths.append( os.path.join(custom_bombcell_paths[i], 'RawWaveforms'))
    else:
        wave_paths = []
        for i in range(n_sessions):
            #check if it is in KS directory
            if os.path.exists(os.path.join(KS_dirs[i], 'RawWaveforms')):
                wave_paths.append( os.path.join(KS_dirs[i], 'RawWaveforms'))
            #Raw waveforms curated via bombcell
            elif os.path.exists(os.path.join(KS_dirs[i], 'qMetrics', 'RawWaveforms')):
                wave_paths.append( os.path.join(KS_dirs[i],'qMetrics', 'RawWaveforms'))
            elif os.path.exists(os.path.join(KS_dirs[i], 'bombcell', 'RawWaveforms')):
                wave_paths.append( os.path.join(KS_dirs[i],'bombcell', 'RawWaveforms'))
            else:
                raise Exception('Could not find RawWaveforms folder')
    #load in a waveform from each session to get the number of channels!
    n_channels = []
    for i in range(n_sessions):
        path_tmp = wave_paths[i]
        file = os.listdir(path_tmp)
        waveform_tmp = np.load(os.path.join(path_tmp,file[0]))
        n_channels.append(waveform_tmp.shape[1])

    os.chdir(tmp)

    #Load channel_pos
    channel_pos = []
    for i in range(n_sessions):
        path_tmp = os.path.join(KS_dirs[i], 'channel_positions.npy')
        pos_tmp = np.load(path_tmp)
        if pos_tmp.shape[0] != n_channels[i]:
            print('Attmepting to fill in missing channel positions')
            pos_tmp = fill_missing_pos(KS_dirs[i], n_channels[i])

        #  Want 3-D positions, however at the moment code only needs 2-D so add 1's to 0 axis position
        pos_tmp = np.insert(pos_tmp, 0, np.ones(pos_tmp.shape[0]), axis = 1)
        channel_pos.append(pos_tmp)

    # Use custom BombCell paths if provided, otherwise search in KS directories
    if custom_bombcell_paths is not None:
        if len(custom_bombcell_paths) != n_sessions:
            raise ValueError(f"Number of custom BombCell paths ({len(custom_bombcell_paths)}) must match number of KS directories ({n_sessions})")
        unit_label_paths = custom_bombcell_paths
    else:
        unit_label_paths = []
        # load Good unit Paths
        for i in range(n_sessions):
            if os.path.exists(os.path.join(KS_dirs[i], 'cluster_bc_unitType.tsv')):
               unit_label_paths.append( os.path.join(KS_dirs[i], 'cluster_bc_unitType.tsv')) 
               print('Using BombCell: cluster_bc_unitType')
            else:
                unit_label_paths.append( os.path.join(KS_dirs[i], 'cluster_group.tsv'))
                print('Using cluster_group.tsv')
    
    return wave_paths, unit_label_paths, channel_pos

def get_probe_geometry(channel_pos, param, verbose = False):
    """
    From the channel positions, estimate the number of shanks and there spacing used to 
    identify a unit to a probe.
    needmin_new_shank_distance from param.

    Parameters
    ----------
    channel_pos : ndarray
        The channel positions of a session 
    param : dict
        The param dict
    verbose : bool, optional
        If True will print the calculated results, by default False

    Returns
    -------
    param
        The param dictionary updated with the calculated params
    """
    min_new_shank_distance = param['min_new_shank_distance']
    x_val = np.unique(channel_pos[:,0])
    x_val = np.sort(x_val) #make sure they are in ascending order 
    x_spacing = np.diff(x_val)
    too_close = np.argwhere(x_spacing < min_new_shank_distance)

    #remove when shanks are two close
    x_val = np.delete(x_val, too_close + 1)

    n_shanks = x_val.size
    if n_shanks == 1:
        shank_spacing = 100 #make it cover the full possible range
    else:
        shank_spacing = np.abs(np.diff(x_val))[0]

    #NOTE this shank_dist is the distance within a centroid would be assigned to a shank!
    shank_dist = shank_spacing * 0.9 #the area from the start which is assigned to each probe
    #e.g 0 -> 0.9to shank 1, 0.9-1.8 to shank 2...
    #NOTE current shank identification doesn't work for 9+ shanks

    param['no_shanks'] = n_shanks
    param['shank_dist'] = shank_dist
    if verbose == True:
        print(f'We have found {n_shanks} with spacing ~ {shank_spacing}')
    return param

def read_datapaths(mice, base):
    """
    Input should be a list of mouse names as strings, e.g. ["AL031", ...]
    Output is a dictionary uniquely identifying each (mouse, probe, location) and the relevant recordings.
    Structure of output given below...

    raw_waveforms_dict["mouse"] : a list of mouse names (can be repeats)
    raw_waveforms_dict["probe"] : the corresponding probe for each entry in the mouse list
    raw_waveforms_dict["loc"] : the corresponding locations on the probe recordings were taken from
    raw_waveforms_dict["recordings"] : corresponding list where each entry is a numpy array listing all the paths to the recordings folders.
    """
    if type(mice) == str:
        # Sanitise inputs so that a single string can be passed in rather than a list.
        mice = [mice]

    # Initialise output dictionary
    raw_waveforms_dict = {}
    raw_waveforms_dict["mouse"] = []
    raw_waveforms_dict["probe"] = []
    raw_waveforms_dict["loc"] = []
    raw_waveforms_dict["recordings"] = []

    # Find Unitmatch.mat for each recording
    for mouse in mice:
        name_path = os.path.join(base, mouse)
        probes = os.listdir(name_path)
        for probe in probes:
            name_probe = os.path.join(name_path, probe)
            locations = os.listdir(name_probe)
            for location in locations:
                name_probe_location = os.path.join(name_probe, location)
                if not os.path.isdir(name_probe_location):
                    continue
                if not os.path.exists(os.path.join(name_probe_location, "UnitMatch")):
                    print(f"No UnitMatch folder where it was expected for mouse {mouse}")
                else:
                    datapath = os.path.join(name_probe_location, "UnitMatch")
                    try:
                        f = mat73.loadmat(os.path.join(datapath, "UnitMatch.mat"), verbose=False)
                    except:
                        pass
                    # find the directory to look for the raw waveforms
                    paths = f["UMparam"]["KSDir"]

                    # build the dictionary containing all relevant information
                    raw_waveforms_dict["mouse"].append(mouse)
                    raw_waveforms_dict["probe"].append(probe)
                    raw_waveforms_dict["loc"].append(location)
                    raw_waveforms_dict["recordings"].append(np.array(paths))

    return raw_waveforms_dict

def get_valid_locations(root):
    """Pre-filter valid locations to avoid repeated file system checks"""
    valid_locations = []
    
    for mouse in os.listdir(root):
        mouse_path = os.path.join(root, mouse)
        if not os.path.isdir(mouse_path):
            continue
            
        for probe in os.listdir(mouse_path):
            probe_path = os.path.join(mouse_path, probe)
            if not os.path.isdir(probe_path):
                continue
                
            for loc in os.listdir(probe_path):
                loc_path = os.path.join(probe_path, loc)
                if not os.path.isdir(loc_path):
                    continue
                    
                mt_path = os.path.join(loc_path, "merged_mt.csv")
                if os.path.exists(mt_path):
                    valid_locations.append((mouse, probe, loc, mt_path))
    
    return valid_locations

def get_exp_id(experiment_path:str, mouse:str):
    """
    experiment_path should be a full absolute path to the desired experiment folder on the server
    (NOT LOCAL PATH)
    """
    exp_name = os.path.basename(os.path.dirname(os.path.dirname(experiment_path)))
    experiment_id = exp_name[exp_name.find(mouse):]
    experiment_id = experiment_id.replace(mouse, '')
    experiment_id = experiment_id.replace("\\", "_")
    return experiment_id

def filter_good_units_and_merge(mt_path, mouse, KS_dirs, waveforms, session_switch, param):
    """
    Handles merges but dimensions not guaranteed to match SQL database.
    """

    filtered_waveforms = []
    filtered_ses_switch = [0]
    filtered_within_session = []
    filtered_session_id = []
    good_units = []

    exp_names = set()
    for i, experiment in enumerate(KS_dirs):
        experiment_id = get_exp_id(experiment, mouse)
        experiment_id = experiment_id[:60]
        if experiment_id in exp_names:
            experiment_id += f"_{i+1}"
        exp_names.add(experiment_id)

        good_ids = []
        merges = []
        if not os.path.exists(os.path.join(os.path.dirname(mt_path), experiment_id, "processed_waveforms")):
            good_units.append(good_ids)
            print(f"Warning: Could not find processed_waveforms for {experiment_id}...")
            continue

        wf_files = os.listdir(os.path.join(os.path.dirname(mt_path), experiment_id, "processed_waveforms"))
        for file in wf_files:
            id = get_unit_id(file)
            if id is not None:
                if type(id) == int:
                    good_ids.append(id)
                elif '+' in id:
                    id1 = int(id.split('+')[0])
                    id2 = int(id.split('+')[1])
                    merges.append((id1, id2))

        indices = np.arange(session_switch[i], session_switch[i + 1])
        session_waveforms = waveforms[indices]
        for merge in merges:
            if max(merge) in good_ids:
                good_ids.remove(max(merge))
            session_waveforms[min(merge)] = 0.5 * (session_waveforms[merge[0]] + session_waveforms[merge[1]])
            if min(merge) not in good_ids:
                good_ids.append(min(merge))
        
        good_ids.sort()
        filtered = session_waveforms[good_ids]
        if np.isnan(filtered).any():
            filtered = fill_missing_values(filtered)
        filtered_waveforms.append(filtered)
        filtered_ses_switch.append(filtered_ses_switch[-1] + filtered.shape[0])
        filtered_session_id.extend([i] * filtered.shape[0])
        good_units.append(good_ids)

    filtered_session_id = np.array(filtered_session_id)
    filtered_within_session = (filtered_session_id[:, None] != filtered_session_id).astype(int)
    filtered_waveforms = np.concatenate(filtered_waveforms, axis=0)
    filtered_session_id = np.array(filtered_session_id)
    filtered_ses_switch = np.array(filtered_ses_switch)
    param['n_units'], param['n_sessions'] = filtered_waveforms.shape[0], len(filtered_ses_switch) - 1

    return filtered_waveforms, filtered_ses_switch, filtered_within_session, filtered_session_id, good_units, param

def filter_good_units(mouse, probe, loc, conn, waveforms, session_switch, param):
    """
    Guaranteed to work and match shape of match tables in SQL database, but ignores merges. 
    """

    filtered_waveforms = []
    filtered_ses_switch = [0]
    filtered_within_session = []
    filtered_session_id = []
    good_units = []

    mt = pd.read_sql_query(f"SELECT RecSes1,ID1 FROM '{mouse}_{probe}_{loc}'", conn)
    ids = {}
    sessions = mt['RecSes1'].unique()
    sessions.sort()

    for i, session in enumerate(sessions):
        df = mt.loc[mt['RecSes1'] == session]
        good_ids = df['ID1'].unique()
        good_ids.sort()
        ids[session] = good_ids
        indices = np.arange(session_switch[i], session_switch[i + 1])
        session_waveforms = waveforms[indices]
        filtered = session_waveforms[good_ids]
        if np.isnan(filtered).any():
            filtered = fill_missing_values(filtered)
        filtered_waveforms.append(filtered)
        filtered_ses_switch.append(filtered_ses_switch[-1] + filtered.shape[0])
        filtered_session_id.extend([i] * filtered.shape[0])
        good_units.append(good_ids)
    
    filtered_session_id = np.array(filtered_session_id)
    filtered_within_session = (filtered_session_id[:, None] != filtered_session_id).astype(int)
    filtered_waveforms = np.concatenate(filtered_waveforms, axis=0)
    filtered_session_id = np.array(filtered_session_id)
    filtered_ses_switch = np.array(filtered_ses_switch)
    param['n_units'], param['n_sessions'] = filtered_waveforms.shape[0], len(filtered_ses_switch) - 1

    return filtered_waveforms, filtered_ses_switch, filtered_within_session, filtered_session_id, good_units, param

def get_unit_id(filepath:str):
    "This version of get_unit_id is mostly for internal use. The DeepUnitMatch.ipynb notebook uses the one in DeepUnitMatch/utils/helpers.py"
    fp = os.path.basename(filepath)
    if fp[:4] == "Unit" and fp[-14:] == "_RawSpikes.npy":
        fp = fp.replace("Unit", "")
        id = fp.replace("_RawSpikes.npy", "")
        if '+' in id:
            return id
        if '#' in id:
            return None
        try:
            return int(id)
        except:
            print(id)
            raise ValueError(f"Invalid filepath format for this waveform: {filepath}", 
                         "Filename for waveform XX should be UnitXX_RawSpikes.npy")
    else:
        raise ValueError(f"Invalid filepath format for this waveform: {filepath}", 
                         "Filename for waveform XX should be UnitXX_RawSpikes.npy")

def retry_db_operation(operation_func, max_retries=None, initial_delay=0.1, max_delay=30.0, 
                      backoff_multiplier=2.0, jitter_range=0.1):
    """
    Retry database operations with exponential backoff to handle SQLite locking issues.
    
    Parameters
    ----------
    operation_func : callable
        Function to execute that may fail due to database locking
    max_retries : int, optional
        Maximum number of retries. If None, retry indefinitely
    initial_delay : float, default 0.1
        Initial delay in seconds before first retry
    max_delay : float, default 30.0
        Maximum delay between retries in seconds
    backoff_multiplier : float, default 2.0
        Multiplier for exponential backoff
    jitter_range : float, default 0.1
        Random jitter as fraction of delay to add randomness
    
    Returns
    -------
    Any
        Result of the operation_func if successful
        
    Raises
    ------
    sqlite3.Error
        If max_retries is exceeded and operation still fails
    """
    attempt = 0
    delay = initial_delay
    
    while True:
        try:
            return operation_func()
        except (sqlite3.OperationalError, sqlite3.DatabaseError) as e:
            # Check if it's a database lock error
            error_str = str(e).lower()
            if ("database is locked" in error_str or 
                "cannot commit" in error_str or 
                "database disk image is malformed" in error_str or
                "attempt to write a readonly database" in error_str):
                
                attempt += 1
                
                # If max_retries is set and exceeded, raise the error
                if max_retries is not None and attempt > max_retries:
                    print(f"Database operation failed after {max_retries} retries: {e}")
                    raise e
                
                # Calculate delay with jitter
                jitter = random.uniform(-jitter_range, jitter_range) * delay
                actual_delay = min(delay + jitter, max_delay)
                
                if attempt <= 5 or attempt % 10 == 0:  # Reduce log spam
                    print(f"Database locked, retrying in {actual_delay:.2f} seconds (attempt {attempt})...")
                
                time.sleep(actual_delay)
                
                # Exponential backoff
                delay = min(delay * backoff_multiplier, max_delay)
            else:
                # For other types of database errors, re-raise immediately
                raise e

def create_robust_db_connection(db_path, timeout=30.0):
    """
    Create a robust SQLite database connection with retry logic.
    
    Parameters
    ----------
    db_path : str
        Path to the SQLite database file
    timeout : float, default 30.0
        Database timeout in seconds
        
    Returns
    -------
    tuple
        (connection, cursor) objects
    """
    def _connect():
        conn = sqlite3.connect(db_path, timeout=timeout)
        # Set WAL mode for better concurrency (allows readers while writing)
        conn.execute("PRAGMA journal_mode=WAL")
        # Set busy timeout
        conn.execute(f"PRAGMA busy_timeout={int(timeout * 1000)}")
        # Optimize for concurrent access
        conn.execute("PRAGMA synchronous=NORMAL")  # Faster than FULL, still safe with WAL
        conn.execute("PRAGMA cache_size=10000")    # Larger cache
        conn.execute("PRAGMA temp_store=memory")   # Use memory for temporary tables
        cursor = conn.cursor()
        return conn, cursor
    
    return retry_db_operation(_connect)
    
def add_col_to_sql(conn, cursor, mouse, probe, loc, col_name, data, indices=None, overwrite=False):
    table_name = f"{mouse}_{probe}_{loc}"
    
    # Check if column already exists
    cursor.execute(f"PRAGMA table_info({table_name})")
    columns = [column[1] for column in cursor.fetchall()]
    
    if col_name in columns:
        if not overwrite:
            raise ValueError(f"Column '{col_name}' already exists in table '{table_name}'. Use overwrite=True to replace it.")
    else:
        cursor.execute(f"ALTER TABLE {table_name} ADD COLUMN {col_name} REAL")

    cursor.execute("BEGIN TRANSACTION")
    cursor.execute("CREATE TEMP TABLE temp_update (rowid INTEGER PRIMARY KEY, value REAL)")

    cursor.execute(f"SELECT rowid FROM {table_name} ORDER BY rowid")
    rowids = [item[0] for item in cursor.fetchall()]

    if indices is None:
        assert len(rowids) == len(data), f"Length of data ({len(data)}) does not match number of rows ({len(rowids)})."
        insert_data = [(int(rowid), float(d)) for rowid, d in zip(rowids, data)]
    else:
        assert len(indices) == len(data), f"Length of data ({len(data)}) does not match length of indices ({len(indices)})."
        insert_data = [(int(rowids[i]), float(d)) for i, d in zip(indices, data)]
    
    cursor.executemany("INSERT INTO temp_update VALUES (?, ?)", insert_data)
    cursor.execute(f"""
        UPDATE {table_name} 
        SET {col_name} = (SELECT value FROM temp_update WHERE temp_update.rowid = {table_name}.rowid)
        WHERE rowid IN (SELECT rowid FROM temp_update)
    """)
    cursor.execute("DROP TABLE temp_update")
    cursor.execute("COMMIT")

def fill_missing_values(waveforms):
    """
    Fills missing values when 1 cv repeat is missing by just copying the existing values.

    Arguments
    ---------
    waveforms : ndarray (n_units, n_timepoints, n_channels, n_cv)
        The waveform array with possible NaN values
    """

    # find the units that have NaN values
    indices = np.argwhere(np.isnan(waveforms))
    units = np.unique(indices[:, 0])
    for i in units:
        waveforms[i] = np.nanmean(waveforms[i], axis=-1, keepdims=True)
    return waveforms
