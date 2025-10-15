import numpy as np
import UnitMatchPy.param_functions as pf

def re_scale(vector):
    """
    scale a score between 0 and 1 to use for a probability distribution.

    Parameters
    ----------
    vector : ndarray 
        A array of scores

    Returns
    -------
    ndarray
        re-scaled score array
    """
    score = ( np.nanquantile(vector,0.99) - vector) / (np.nanquantile(vector,0.99) - np.nanmin(vector))
    score[score<0] = 0
    score[np.isnan(score)] = 0 # set NaN values to 0
    return score

def get_simple_metric(waveform_parameter, outlier = False):
    """
    This function is suitable for, spatial decay, spatial decay fit, amplitude and other (n_units,2) parameters,

    Parameters
    ----------
    waveform_parameter : ndarray (n_units , 2)
        A simple parameters where each unit/cv has one score
    outlier : bool, optional
        If True will remove more outliers before creating a score, by default False

    Returns
    -------
    ndarray
        A (n_units, n_units) array of scores
    """
    n_units = waveform_parameter.shape[0]
    #select each cv then broadcast to n_units * n_nits as to compare each pair of units
    x1 = np.broadcast_to(waveform_parameter[:,0], (n_units, n_units)).T
    x2 = np.broadcast_to(waveform_parameter[:,1], (n_units, n_units))

    # takes the difference weighted by the mean of the value of the two cross-validation halves
    # to make the difference between 99-100 less significant than 5-6
    diff = np.abs(x1 - x2) / np.nanmean(np.abs( np.stack((x1,x2), axis = -1)), axis = 2)

    if outlier == True: 
        diff[diff<0] = np.nanquantile(diff, 0.9999) 
        diff[ diff > np.nanquantile(diff,0.9999)] = np.nanquantile(diff, 0.9999)

    sqrt_diff = np.sqrt(diff) 
    metric = re_scale(sqrt_diff)

    return metric

def get_wave_corr(avg_waveform, param):
    """
    Calculates the correlation between weighted average waveform, and rescales it into a score

    Parameters
    ----------
    avg_waveform : ndarray (n_units, spike_width, 2)
        The average waveform
    param : dict
        The param dictionary

    Returns
    -------
    ndarray
        The(n_units, n_units) score array for the waveform correlation
    """
    waveidx = param['waveidx']
    n_units = param['n_units']

    x1 = avg_waveform[waveidx,:,0].T 
    x2 = avg_waveform[waveidx,:,1].T

    #calculate the correlation
    wave_corr_tmp = np.corrcoef(x1,x2)
    wave_corr = wave_corr_tmp[:n_units,n_units:]


    wave_corr = np.arctanh(wave_corr) # apply Fisher z transformation
    #make into a score
    wave_corr = (wave_corr - np.nanquantile(wave_corr,0.005)) / (np.nanquantile(wave_corr,0.995) - np.nanquantile(wave_corr, 0.005))
    wave_corr[np.isnan(wave_corr)] = 0
    wave_corr[wave_corr<0] = 0
    wave_corr[wave_corr>1] = 1 

    return wave_corr

def get_waveforms_mse(avg_waveform, param):
    """
    Calculates the waveform mean square error, and rescales it into a score

    Parameters
    ----------
    avg_waveform : ndarray (n_units, spike_width, 2)
        The average waveform
    param : dict
        The param dictionary

    Returns
    -------
    ndarray
        The(n_units, n_units) score array for the waveform mean square error
    """
    waveidx = param['waveidx']
    projected_waveform_norm = avg_waveform[waveidx,:,:]
    projected_waveform_norm =  (projected_waveform_norm - np.nanmin(projected_waveform_norm,axis = 0)) / (np.nanmax(projected_waveform_norm, axis=0) - np.nanmin(projected_waveform_norm, axis = 0))

    x1 = np.tile(np.expand_dims(projected_waveform_norm[:,:,0], axis = 2), (1,1,projected_waveform_norm.shape[1]))
    x2 = np.swapaxes(np.tile(np.expand_dims(projected_waveform_norm[:,:,1], axis = 2), (1,1,projected_waveform_norm.shape[1])), 1,2 )
    raw_wave_mse = np.nanmean( (x1 - x2)**2, axis = 0 ).squeeze()
    raw_wave_mse_norm = np.sqrt(raw_wave_mse)

    waveform_mse = re_scale(raw_wave_mse_norm)
    return waveform_mse

def flip_dim(avg_waveform_per_tp, param):
    """
    Creates a version of the weighted average waveform per time point, where the x axis is flipped, due to the effect
    where the average position tends to wards the center of the x coords when the wave is decaying

    Parameters
    ----------
    avg_waveform_per_tp : ndarray
        The average waveform per time point
    param : _type_
        The param dictionary

    Returns
    -------
    ndarray
        The average waveform per time point with the x axis flipped and not flipped 
    """
    n_units = param['n_units']
    spike_width = param['spike_width']

    flip_dim = np.array((1,)) # BE CAREFUL HERE, Which dimension is the x-axis  
    avg_waveform_per_tp_flip = np.full((3, n_units,spike_width ,2 , len(flip_dim)+1), np.nan)

    for i in range(len(flip_dim)):
        tmpdat = avg_waveform_per_tp[flip_dim[i]  ,:,:,:]

        new_vals = np.nanmin(tmpdat, axis =1, keepdims=True) + np.nanmax(tmpdat, axis = 1, keepdims=True) - tmpdat
        avg_waveform_per_tp_flip[:,:,:,:,i] = avg_waveform_per_tp
        avg_waveform_per_tp_flip[flip_dim[i],:,:,:,i] = new_vals

    avg_waveform_per_tp_flip[:,:,:,:,-1] = avg_waveform_per_tp

    return avg_waveform_per_tp_flip

def get_Euclidean_dist(avg_waveform_per_tp_flip,param):
    """
    Calculated the Euclidean distance between the units at each time point and for the flipped axis case

    Parameters
    ----------
    avg_waveform_per_tp_flip : ndarray
        The average waveform per time point with the x axis flipped and not flipped 
    param : _type_
        The param dictionary

    Returns
    -------
    ndarray
        The euclidean distance between each unit for each time point
    """
    # Memory-efficient implementation using loops instead of massive tiling
    waveidx = param['waveidx']
    n_units = param['n_units']
    
    # Extract the relevant slices once
    data_cv0 = avg_waveform_per_tp_flip[:,:,waveidx,0,:]  # (3, n_units, len(waveidx), n_flips)
    data_cv1 = avg_waveform_per_tp_flip[:,:,waveidx,1,:]  # (3, n_units, len(waveidx), n_flips)
    
    # Initialise output array to match original
    euclid_dist = np.full((n_units, len(waveidx), data_cv0.shape[-1], n_units), np.nan)
    
    # Process unit by unit
    for i in range(n_units):
        for j in range(n_units):
            # Get data for units i and j
            unit_i_data = data_cv0[:, i, :, :]  # (3, len(waveidx), n_flips)
            unit_j_data = data_cv1[:, j, :, :]  # (3, len(waveidx), n_flips)
            
            # Check for NaN values
            w = np.isnan(np.abs(unit_i_data[0,:,:] - unit_j_data[0,:,:]))
            
            # Compute Euclidean distance between units i and j
            diff = unit_i_data - unit_j_data  # (3, len(waveidx), n_flips)
            tmp_euclid = np.linalg.norm(diff, axis=0)  # (len(waveidx), n_flips)
            tmp_euclid[w] = np.nan
            
            # Store result
            euclid_dist[i, :, :, j] = tmp_euclid
    
    return euclid_dist

def centroid_metrics(euclid_dist, param):
    """
    This function calculates the score for the centroid distance and centroid variance.

    Parameters
    ----------
    euclid_dist : ndarray
        The euclidean distance between each unit for each time point
    param : _type_
        The param dictionary

    Returns
    -------
    ndarrays
        An array (n_units, n_units) of scores for centroid distance and variance
    """
    max_dist = param['max_dist']
    waveidx = param['waveidx']
    new_peak_loc = param['peak_loc']    

    centroid_dist = np.nanmin( euclid_dist[:,new_peak_loc - waveidx ==0,:,:].squeeze(), axis =1 ).squeeze()


    centroid_dist = 1 - ((centroid_dist - np.nanmin(centroid_dist)) / (max_dist - np.nanmin(centroid_dist)))
    centroid_dist[centroid_dist<0] = 0
    centroid_dist[np.isnan(centroid_dist)] = 0

    #Centroid Var
    # need ddof = 1 to match with ML
    centroid_var = np.nanmin( np.nanvar(euclid_dist, axis = 1, ddof = 1 ).squeeze(), axis =1 ).squeeze()
    centroid_var = np.sqrt(centroid_var)
    centroid_var = re_scale(centroid_var)
    centroid_var[np.isnan(centroid_var)] = 0

    return centroid_dist, centroid_var

def get_recentered_euclidean_dist(avg_waveform_per_tp_flip, avg_centroid, param):
    """
    Find a Euclidean distance where the location per time has been centered around the average position

    Parameters
    ----------
    avg_waveform_per_tp_flip : ndarray
        The average waveform per time point with the x axis flipped and not flipped
    avg_centroid : ndarray
        The average location for each unit
    param : dict
        The param dictionary

    Returns
    -------
    ndarray
        the re-centered euclidean distance
    """

    waveidx = param['waveidx']
    n_units = param['n_units']
    spike_width = param['spike_width']

    # recentered projected location , aka subtract the avg location, the we can unique info
    avg_centroid_broadcast = np.tile(np.expand_dims(avg_centroid, axis= (3,4)), (1,1,1,spike_width,2))
    avg_waveform_per_tp_fip_recenterd = np.swapaxes( np.swapaxes(avg_waveform_per_tp_flip, 2,3) - avg_centroid_broadcast,2,3)
    x1 = np.tile( np.expand_dims(avg_waveform_per_tp_fip_recenterd[:,:,waveidx,0,:], axis = -1), (1,1,1,1,n_units)).squeeze()
    x2 = np.swapaxes(np.tile( np.expand_dims(avg_waveform_per_tp_fip_recenterd[:,:,waveidx,1,:], axis = -1), (1,1,1,1,n_units)).squeeze(), 1, 4)

    w = np.isnan( np.abs(x1[0,:,:,:,:] - x2[0,:,:,:,:])).squeeze()

    tmp_euclid = np.linalg.norm(x1-x2, axis = 0)
    tmp_euclid[w] = np.nan
    euclid_dist_2 = tmp_euclid
    del x1
    del x2
    del w
    del tmp_euclid
    return euclid_dist_2

def recentered_metrics(euclid_dist_2):
    """
    Calculates the re-centered centroid distance score metric

    Parameters
    ----------
    euclid_dist_2 : ndarray
        The re-centered euclidean distance

    Returns
    -------
    ndarray
        The centroid distance re_centered score metric
    """
    centroid_dist_recentered = np.nanmin( np.nanmean(euclid_dist_2, axis =1), axis =1)
    centroid_dist_recentered = re_scale(centroid_dist_recentered)
    centroid_dist_recentered[np.isnan(centroid_dist_recentered)] = 0
    return centroid_dist_recentered

def dist_angle(avg_waveform_per_tp_flip, param):
    """
    This function uses the weighted average location per time point, to find metrics based of off:
    The distance traveled by the unit at each time point
    The angle between units at each time point 

    Parameters
    ----------
    avg_waveform_per_tp_flip : ndarray
        The average waveform per time point with the x axis flipped and not flipped
    param : dict
        The param dictionary

    Returns
    -------
    ndarrays
        The score arrays for the trajectory distance and angle
    """
    waveidx = param['waveidx']
    n_units = param['n_units']
    min_angle_dist = param['min_angle_dist']
    
    #Distance between time steps and angle
    x1 = avg_waveform_per_tp_flip[:,:,waveidx[1]:waveidx[-1] +1,:,:]
    x2 = avg_waveform_per_tp_flip[:,:,waveidx[0]:waveidx[-2] +1,:,:] # Difference between python and ML indexing



    traj_dist = np.linalg.norm(x1-x2, axis= 0)

    loc_angle = np.full(np.append(traj_dist.shape, 3), np.nan)
    #only select points which have enough movement to get a angle
    good_angle = np.zeros_like(traj_dist)
    good_angle[traj_dist>=min_angle_dist] = 1

    count_id = 0
    for dim_id1 in range(avg_waveform_per_tp_flip.shape[0]):
        for dim_id2 in np.arange(1,avg_waveform_per_tp_flip.shape[0]):
            if dim_id2 <= dim_id1:
                continue
            ang = np.abs( x1[dim_id1,:,:,:,:] - x2[dim_id1,:,:,:,:]) / (np.abs(x1[dim_id2,:,:,:,:] - x2[dim_id2,:,:,:,:]) + 1e-10)
            
            loc_angle[:,:,:,:,count_id] = np.arctan(ang) * good_angle # only selects angles for units where there is sufficient distance between time poitns
            count_id +=1


    loc_angle = np.nansum(loc_angle, axis=4)
    x1 = np.tile (np.expand_dims(loc_angle[:,:,0,:], axis = -1), (1,1,1,n_units))
    x2 = np.swapaxes(np.tile(np.expand_dims(loc_angle[:,:,1,:], axis = -1), (1,1,1,n_units)), 0,3)


    angle_subtraction = np.abs(x1-x2)
    angle_subtraction[np.isnan(np.abs(x1 - x2))] = 2 * np.pi # make nan values punished 

    traj_angle_sim = np.nanmin( np.nansum(angle_subtraction, axis=1), axis = 1)
    traj_angle_sim = re_scale(traj_angle_sim)
    traj_angle_sim[np.isnan(traj_angle_sim)] = 0

    x1 = np.tile (np.expand_dims(traj_dist[:,:,0,:], axis = -1), (1,1,1,n_units))
    x2 = np.swapaxes(np.tile(np.expand_dims(traj_dist[:,:,1,:], axis = -1), (1,1,1,n_units)), 0,3)

    traj_dist_compared = np.abs(x1-x2)

    traj_dist_sim = np.nanmin( np.nansum(traj_dist_compared , axis = 1), axis= 1)
    traj_dist_sim = np.sqrt(traj_dist_sim)
    traj_dist_sim = re_scale(traj_dist_sim)
    traj_dist_sim[np.isnan(traj_dist_sim)] = 0

    return traj_angle_sim, traj_dist_sim


def get_threshold(total_score, within_session, euclid_dist, param, is_first_pass = True):
    """
    Uses the total_score and Euclidean distance, to determine a threshold for putative matches.

    If it is the first pass through the data i.e no drift correction has been done, we would expect the
    Total score for the matches to be smaller than expected, therefore we calculate the difference in mean
    for within and and between session to lower the threshold

    Parameters
    ----------
    total_score : ndaray 
        The total scores for each metric
    within_session : ndarray
        An array marking what session each unit is in
    euclid_dist : ndarray
        The euclidean distance between each units centroid
    param : dict
        The param dictionary
    is_first_pass : bool, optional
        If it is the first pass i.e has between session drift correction been applied, by default True

    Returns
    -------
    float
        A threshold for matches / non-matches
    """
    # take between session out
    score_vector = param['score_vector']
    Bins = param['bins']

    tmp = total_score.copy()
    tmp[euclid_dist > param['neighbour_dist'] ] = np.nan

    tmp[within_session == 1] = np.nan

    hd, __ = np.histogram(np.diag(tmp), Bins)
    hd = hd /  param['n_units']
    hnd, __ = np.histogram( (tmp - tmp *np.eye(param['n_units'])), Bins)
    hnd = hnd / np.nansum( (tmp - tmp *np.eye(param['n_units'])) )
    hp, __ = np.histogram(tmp, Bins)
    hp = hp / np.nansum(tmp)

    thrs_opt = score_vector[np.argwhere( (pf.smooth(hd,3) > pf.smooth(hnd,3) ) * (score_vector > 0.6) == True)][0]
    # if ThrsOpt.size == 0:
    #     ThrsOpt = 0.6 # give default threshold if above doestn return value
    # fit tmp to a normal ditn
    fit = tmp[ ~np.isnan(tmp) * (tmp < thrs_opt)]
    muw = np.mean(fit)
    stdw = np.std(fit)

    if param['n_sessions'] > 1:
        if is_first_pass == True:
            # take within session out

            tmp = total_score.copy()
            tmp[euclid_dist > param['neighbour_dist'] ] = np.nan

            tmp[within_session == 0] = np.nan

            ha, __ = np.histogram(tmp, Bins)
            ha = ha / np.nansum(tmp)
            fit = tmp[ ~np.isnan(tmp) * (tmp < thrs_opt)]
            mua = np.mean(fit)
            stda = np.std(fit)


            # for first pass only (i.e before drift correction)
            #This is to decrease the threshold for the first pass so that the threshold is lowered
            # as without drift correction even matches should have a lower total score
            if (~np.isnan(mua) and mua<muw):
                thrs_opt = thrs_opt - np.abs(muw - mua)

    return thrs_opt

def get_good_matches(pairs, total_score):
    """
    This function takes in a list of potential matches (n, 2) and return a list of matches where each unit can only appear once.
    This mean one unit will not be matched to multiple units providing a better estimate at the cost of loosing some units.
    The "best" match is decided by the match which has the highest total score.

    Parameters
    ----------
    pairs : ndarray
        A list of matches
    total_score : ndarray
        The total score array, the sum of all of the individual metric arrays

    Returns
    -------
    ndarray
        A list of matches where each unit appears once
    """
    # Need to make sure the first and second unit in the matches only appears once
    for pair_id in range(2):

        idx, count = np.unique(pairs[:,pair_id], return_counts = True)
        ls = np.argwhere(count != 1)
        tmp_vals = idx[ls] # returns the unit idx for PairID, where there is more than one potential match

        #Go through each case where there is more than 1 match
        for i in range(len(tmp_vals)):
            # find the unit idx pair, e.g if unit 2 matches with unit 272 and 278 this will find (2,272) then (2,278)
            tmp_pair = np.argwhere(pairs[:,pair_id] == tmp_vals[i]) # idx of pair for the multiple matched unit
            tmp = pairs[tmp_pair,:].squeeze() # unit idx pair, for each double match e.g (2,272) and (2,278)

            scores = np.zeros(len(tmp))
            for z in range(len(tmp)):
                scores[z] = total_score[tmp[z,0], tmp[z,1]] #Get their score

            best_match = np.argmax(scores)
            #set the worse matches to -1
            for z in range(len(tmp)):
                if z != best_match:
                    # cannot remove yet, as it will change all of the found indices, so set the value to -1, the at the end can remove all apperaances of -1
                    pairs[tmp_pair[z], :] = np.full_like(pairs[tmp_pair[z], :], -1)
    
    good_pairs = np.delete(pairs, np.argwhere(pairs[:,0] == -1), axis = 0)

    return good_pairs
   

def drift_correction_basic(candidate_pairs, session_switch, avg_centroid, avg_waveform_per_tp):
    """
    Uses the median difference in position, between putative matches to gain a value of drift between sessions,
    this is then applied to the AvgCentroid and the AvgWaveformPerTP.
    This function works on a pair of sessions.

    Parameters
    ----------
    candidate_pairs : ndarray
        The current pairs which may be a match
    session_switch : ndarray
        A array which contains at what unit numbers does a new session start
    avg_centroid : ndarray
        The average centroid for each unit
    avg_waveform_per_tp : ndarray
        The average waveform per time point for each unit

    Returns
    -------
    ndarray, ndarray, ndarray
        The array of drift correction values, the avg_centroid and avg_waveform_per_tp updated with the drift correction
    """
    #Drift.. currently only doing drift correction between 2 days/sessions
    best_pairs = np.argwhere(candidate_pairs == 1)

    # #Just to make it like the matlab code, as in matlab the axes are swapped
    best_pairs[:, [0,1]] = best_pairs[:, [1,0]]
    idx = np.argwhere( ((best_pairs[:,0] < session_switch[1]) * (best_pairs[:,1] >= session_switch[1])) == True)

    drift = np.nanmedian( np.nanmean( avg_centroid[:, best_pairs[idx,0].squeeze(),:], axis = 2) - np.nanmean( avg_centroid[:,best_pairs[idx,1].squeeze(),:], axis = 2), axis = 1)

    ##need to add the drift to the location
    avg_waveform_per_tp[0,session_switch[1]:,:,:] += drift[0]
    avg_waveform_per_tp[1,session_switch[1]:,:,:] += drift[1]
    avg_waveform_per_tp[2,session_switch[1]:,:,:] += drift[2]

    avg_centroid[0,session_switch[1]:,:] += drift[0]
    avg_centroid[1,session_switch[1]:,:] += drift[1]
    avg_centroid[2,session_switch[1]:,:] += drift[2]

    return drift, avg_centroid, avg_waveform_per_tp

def apply_drift_correction_basic(pairs, sid, session_switch, avg_centroid, avg_waveform_per_tp):
    """
    This function applies the basic style drift correction to a pair of sessions, as part of a n_session drift correction  

    Parameters
    ----------
    pairs : ndarray
        A list of potential matches
    sid : int
        The session id
    session_switch : ndarray
        What units does the session switch
    avg_centroid : ndarray
        The average centroid for each unit
    avg_waveform_per_tp : ndarray
        The average waveform per time point for each unit

    Returns
    -------
    ndarray, ndarray, ndarray
        The array of drift correction values, the avg_centroid and avg_waveform_per_tp updated with the drift correction, for all sessions
    """

    drift = np.nanmedian( np.nanmean( avg_centroid[:, pairs[:,0],:], axis = 2) - np.nanmean( avg_centroid[:,pairs[:,1],:], axis = 2), axis = 1)


    ##need to add the drift to the location
    avg_waveform_per_tp[0,session_switch[sid+1]:session_switch[sid+2],:,:] += drift[0]
    avg_waveform_per_tp[1,session_switch[sid+1]:session_switch[sid+2],:,:] += drift[1]
    avg_waveform_per_tp[2,session_switch[sid+1]:session_switch[sid+2],:,:] += drift[2]

    avg_centroid[0,session_switch[sid+1]:session_switch[sid+2],:] += drift[0]
    avg_centroid[1,session_switch[sid+1]:session_switch[sid+2],:] += drift[1]
    avg_centroid[2,session_switch[sid+1]:session_switch[sid+2],:] += drift[2]

    return drift, avg_waveform_per_tp, avg_centroid

def apply_drift_correction_per_shank(pairs, sid, session_switch, avg_centroid, avg_waveform_per_tp, param):
    """
    This is the same as "basic" drift correction, however treats each shank separately, if there is enough units per shank

    Parameters
    ----------
    pairs : ndarray
        A list of potential matches
    sid : int
        The session id
    session_switch : ndarray
        What units does the session switch
    avg_centroid : ndarray
        The average centroid for each unit
    avg_waveform_per_tp : ndarray
        The average waveform per time point for each unit

    Returns
    -------
    ndarray, ndarray, ndarray
        The array of drift correction values, the avg_centroid and avg_waveform_per_tp updated with the drift correction, for all sessions
    """
    shank_id = shank_ID_per_session(avg_centroid ,session_switch ,sid , param)
    n_shanks = param['no_shanks']
    shank_dist = param['shank_dist']

    
    centroid_a = np.nanmean( avg_centroid[:,pairs[:,0],:], axis = 2)
    centroid_b = np.nanmean( avg_centroid[:,pairs[:,1],:], axis = 2)

    max_dist = 0
    min_dist = 0
    drift_per_shank = np.zeros([n_shanks,3])

    for i in range(n_shanks):
        max_dist += shank_dist #beginning of loop shift maximum distance

        #get which shank each unit is on
        correct_shank_a = np.logical_and(centroid_a[1,:] < max_dist,  centroid_a[1,:] > min_dist)
        correct_shank_b = np.logical_and(centroid_b[1,:] < max_dist,  centroid_b[1,:] > min_dist)

        if np.all(correct_shank_a == correct_shank_b) != True:
            print(f'These pairs may be bad {np.argwhere(correct_shank_a != correct_shank_b)}')
            #delete pairs which are on different shanks ##sam's fix
            bad_idxs = np.argwhere(correct_shank_a != correct_shank_b)
            correct_shank_a = np.delete(correct_shank_a, bad_idxs)
            correct_shank_b = np.delete(correct_shank_b, bad_idxs)

            
        drifts = centroid_a[:,correct_shank_a] - centroid_b[:,correct_shank_b]
        drift =  np.nanmedian(drifts, axis = 1)
        drift_per_shank[i,:] = drift

        #need to get idx for each shank, to apply correct drift correction

        shank_session_idx = session_switch[sid] + np.argwhere( shank_id == i) 

        avg_waveform_per_tp[0,shank_session_idx,:,:] += drift[0]
        avg_waveform_per_tp[1,shank_session_idx,:,:] += drift[1]
        avg_waveform_per_tp[2,shank_session_idx,:,:] += drift[2]

        avg_centroid[0,shank_session_idx,:] += drift[0]
        avg_centroid[1,shank_session_idx,:] += drift[1]
        avg_centroid[2,shank_session_idx,:] += drift[2]

        min_dist += shank_dist

    return drift_per_shank, avg_waveform_per_tp, avg_centroid

def shank_ID_per_session(avg_centroid ,session_switch ,sid , param):
    """
    This function use the average centroid, to assign each unit in a session to a shank

    Parameters
    ----------
    avg_centroid : ndarray
        The average centroid for each unit
    session_switch : ndarray
        What units does the session switch
    sid : int
        The session id
    param : dict
        The param dictionary

    Returns
    -------
    ndarray
        An array assigning each unit to a shank
    """

    n_shanks = param['no_shanks']
    shank_dist = param['shank_dist']
    max_dist = 0
    min_dist = 0

    #load centroid position for 1 recording session
    centroid_pos = np.nanmean( avg_centroid[:, session_switch[sid]:session_switch[sid + 1],:], axis = 2)
    shank_id = np.zeros(centroid_pos.shape[1])

    #goes through each shank and assigns a unit a shank id depending on the area they are in
    for i in range(n_shanks):
        max_dist += shank_dist
        #Test to see if the centroid position is in the region of the i'th shank
        shank_idx = np.logical_and(centroid_pos[1,:] < max_dist,  centroid_pos[1,:] > min_dist)
        shank_id[shank_idx] = i
        min_dist += shank_dist
       
    return shank_id

def test_matches_per_shank(pairs, avg_centroid, sid, param):
    """
    Checks to see how many matches there are per shank, and returns false if there are less than match num threshold for one shank


    Parameters
    ----------
    pairs : ndarray
        The list of current matches
    avg_centroid : ndarray
        The average centroid for each unit
    sid : int
        The session id
    param : dict
        The param dictionary

    Returns
    -------
    bool
        Will return True if there is enough units per shank to do per shank drift correction
    """

    do_per_shank_correction = True

    centroid_a = np.nanmean( avg_centroid[:,pairs[:,0],:], axis = 2)
    centroid_b = np.nanmean( avg_centroid[:,pairs[:,1],:], axis = 2)
    shank_id_tmp = np.zeros(centroid_a.shape[1])

    units_per_shank_thrs = param['units_per_shank_thrs']

    max_dist = 0
    min_dist = 0
    n_shanks = param['no_shanks']
    shank_dist = param['shank_dist']

    if n_shanks != 1:
        for i in range(n_shanks):
            max_dist += shank_dist

            correct_shank_a = np.logical_and(centroid_a[1,:] < max_dist,  centroid_a[1,:] > min_dist)
            correct_shank_b = np.logical_and(centroid_b[1,:] < max_dist,  centroid_b[1,:] > min_dist)

            if np.all(correct_shank_a == correct_shank_b) != True:
                #delete pairs which are on different shanks
                bad_idx = np.argwhere(correct_shank_a != correct_shank_b)[0]
                print(f'These pairs may be bad {bad_idx}, so are excluded from drift correction')
                correct_shank_a[bad_idx]=False
                correct_shank_b[bad_idx]=False
                
            shank_id_tmp[correct_shank_a] = i

            min_dist += shank_dist

        __, counts = np.unique(shank_id_tmp, return_counts=True)

        if np.any(counts < units_per_shank_thrs):
            do_per_shank_correction = False
            print(f'Session pair {sid+1}/{sid+2} has {counts} matches per shank, which is below threshold to do per shank drift correction')
    else:
        do_per_shank_correction = False
    return do_per_shank_correction


def drift_n_sessions(candidate_pairs, session_switch, avg_centroid, avg_waveform_per_tp, total_score, param, best_match = True, best_drift = True):
    """
    This function applies drift correction between n_sessions, currently this is done by aligning session 2 to session 1,
    then session 3 to session 2 etc.   
    This function, calls another function to apply the drift correction, to be able to apply different type of drift correction 
    easily.

    Parameters
    ----------
    candidate_pairs : ndarray
        An array of likely matches
    session_switch : ndarray
        What units does the session switch
    avg_centroid : ndarray
        The average centroid for each unit
    avg_waveform_per_tp : ndarray
        The average waveform per time point for each unit
    total_score : ndarray
        The summed array for each individual metric
    param : dict
        The param dictionary
    best_match : bool, optional
        If True will only consider the best match for each unit, by default True
    best_drift : bool, optional
        If True will do per shank drift correction if there is sufficient units by default True

    Returns
    -------
    ndarray, ndarray, ndarray
        The drift values, the avg_centroid and the avg_waveform_per_tp arrays drift corrected
    """
    n_sessions = param['n_sessions']
    best_pairs = np.argwhere(candidate_pairs == 1)

    #make it like the matlab code, (small unit idx, larger unit idx)
    best_pairs[:, [0,1]] = best_pairs[:, [1,0]]

    drifts = np.zeros( (n_sessions - 1, 3))

    for did in range(n_sessions - 1):
        idx = np.argwhere( ( (best_pairs[:,0] >= session_switch[did]) * (best_pairs[:,0] < session_switch[did + 1]) *
                            (best_pairs[:,1] >= session_switch[did + 1]) * (best_pairs[:,1] < session_switch[did + 2]) ) == True)

        pairs = best_pairs[idx,:].squeeze()
        if pairs.ndim == 1:
            pairs = np.expand_dims(pairs, axis = 0)
        if best_match == True:
            pairs = get_good_matches(pairs, total_score)

        #Test to see if there are enough matches to do drift correction per shank
        if test_matches_per_shank(pairs, avg_centroid, did, param) == True and best_drift == True:
            drifts = np.zeros( (n_sessions - 1, param['no_shanks'], 3))
            drifts[did,:,:], avg_waveform_per_tp, avg_centroid = apply_drift_correction_per_shank(pairs, did, session_switch, avg_centroid, avg_waveform_per_tp, param)
            print(f'Done drift correction per shank for session pair {did+1} and {did+2}')
        elif len(pairs)>0: #if there exist pairs across sessions:
            drifts = np.zeros( (n_sessions - 1, 3))
            drifts[did,:], avg_waveform_per_tp, avg_centroid = apply_drift_correction_basic(pairs, did, session_switch, avg_centroid, avg_waveform_per_tp)
        elif len(pairs)==0:
            print('No pairs across sessions to perform drift correction.')
    return drifts, avg_centroid, avg_waveform_per_tp



def get_total_score(scores_to_include, param):
    """
    Calculates the total score array (n_units, n_units) normalised sum of all metrics,
    and the predictors array (n_units, n_units, n_scores) a stacked version fo score_to_include needed for probability calculations

    Parameters
    ----------
    scores_to_include : dict
        The dictionary of scores used 
    param : dict
        The param dictionary

    Returns
    -------
    ndarray
        The total scores for each unit, and the array version of scores_to_include
    """
    total_score = np.zeros((param['n_units'],param['n_units']))
    predictors =  np.zeros((param['n_units'],param['n_units'], 0))


    for sid in scores_to_include:
        tmp = scores_to_include[f'{sid}']
        predictors = np.concatenate((predictors, np.expand_dims(tmp, axis = 2)), axis = 2)
        total_score += tmp

    total_score = (total_score - np.min(total_score)) / (np.max(total_score) - np.min(total_score))

    return total_score, predictors
