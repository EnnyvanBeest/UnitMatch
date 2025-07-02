# utility function for loading files etc
import numpy as np
import pandas as pd
import os
import matplotlib.pyplot as plt

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

    for i in range(len(unit_label_paths)):
        if os.path.split(unit_label_paths[0])[1] == 'cluster_bc_unitType.tsv':
            unit_label = load_tsv(unit_label_paths[i])
            tmp_idx = np.argwhere(np.isin(unit_label[:, 1], ['GOOD', 'NON-SOMA GOOD']))
        else:
            unit_label = load_tsv(unit_label_paths[i])
            tmp_idx = np.argwhere(unit_label[:, 1] == 'good')

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
                p_file = os.path.join(wave_paths[ls], f'Unit{int(all_units[ls][0])}_RawSpikes.npy')
                tmp = np.load(p_file)
                tmp_waveform = np.zeros((len(os.listdir(wave_paths[ls])), tmp.shape[0], tmp.shape[1], tmp.shape[2]))

                for i in range(len(os.listdir(wave_paths[ls]))):
                    p_file_good = os.path.join(wave_paths[ls], f'Unit{int(all_units[ls][i])}_RawSpikes.npy')
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
    print(f'The percentage of units matched to themselves is: {self_match}%')
    print(f'The percentage of false -ve\'s then is: {100 - self_match}% \n')

    #off-diagonal miss-matches
    n_off_diag = np.zeros_like(output_prob)
    n_off_diag = output_threshold
    n_off_diag[within_session == 1] = 0 
    n_off_diag[np.eye(param['n_units']) == 1] = 0 
    false_positive_est =  n_off_diag.sum() / (param['n_units']) 
    print(f'The rate of miss-match(es) per expected match {false_positive_est}')


    #compute matlab FP per session per session
    false_positive_est_per_session = np.zeros(param['n_sessions'])
    for did in range(param['n_sessions']):
        tmp_diag = output_threshold[session_switch[did]:session_switch[did + 1], session_switch[did]:session_switch[did + 1]]
        n_units = tmp_diag.shape[0]
        tmp_diag[np.eye(n_units) == 1] = 0 
        false_positive_est_per_session[did] = tmp_diag.sum() / (n_units ** 2 - n_units) * 100
        print(f'The percentage of false +ve\'s is {false_positive_est_per_session[did]}% for session {did +1}')

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