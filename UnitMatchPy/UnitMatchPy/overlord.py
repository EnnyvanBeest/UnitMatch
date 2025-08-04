import UnitMatchPy.param_functions as pf
import UnitMatchPy.metric_functions as mf
import numpy as np

def extract_parameters(waveform, channel_pos, clus_info, param):
    """
    This function runs all of the extract parameters functions needed to run UnitMatch

    Parameters
    ----------
    waveform : ndarray (n_units, spike_width, n_channels, 2)
        The average waveforms needed for UnitMatch
    channel_pos : list
        The complete channel positions for each session
    clus_info : dict
        The clus_info dictionary
    param : dict
        The param dictionary

    Returns
    -------
    dict
        The extracted waveform properties as a dictionary of arrays
    """
    waveform = pf.detrend_waveform(waveform)

    max_site, good_idx, good_pos, max_site_mean = pf.get_max_sites(waveform, channel_pos, clus_info, param)

    spatial_decay_fit , spatial_decay,  d_10, avg_centroid, avg_waveform, peak_time = pf.decay_and_average_waveform(waveform, channel_pos, good_idx, max_site, max_site_mean, clus_info, param)

    amplitude, waveform, avg_waveform = pf.get_amplitude_shift_waveform(waveform, avg_waveform, peak_time, param)

    waveform_duration, avg_waveform_per_tp, good_wave_idxs = pf.get_avg_waveform_per_tp(waveform,channel_pos, d_10, max_site_mean, amplitude, avg_waveform, clus_info, param)

    extracted_wave_properties = {'spatial_decay_fit' : spatial_decay_fit, 'spatial_decay' : spatial_decay , 'avg_centroid' : avg_centroid,
                                'waveform_duration' : waveform_duration, 'avg_waveform_per_tp' : avg_waveform_per_tp, 'good_wave_idxs' : good_wave_idxs,
                                 'amplitude' : amplitude, 'avg_waveform' : avg_waveform, 'max_site' : max_site, 'max_site_mean' : max_site_mean}
    
    return extracted_wave_properties

def extract_metric_scores(extracted_wave_properties, session_switch, within_session, param, niter  = 2, to_use = None):
    """
    This function runs all of the metric calculations and drift correction to calculate the probability
    distribution needed for UnitMatch.

    Parameters
    ----------
    extracted_wave_properties : dict
        The extracted properties from extract_parameters()
    session_switch : ndarray
        An array which indicates when anew recording session starts
    within_session : ndarray
        The array which gives each unit a label depending on their session
    param : dict
        The param dictionary
    niter : int, optional
        The number of pass through the function, 1 mean no drift correction
            2 is one pass of drift correction, by default 2
    to_use : list, optional
        A list of scores to include in the total score calculation. If None (default), all scores are included.

    Returns
    -------
    ndarrays
        The total scores and candidate pairs needed for probability analysis
    """

    #unpack need arrays from the ExtractedWaveProperties dictionary
    amplitude = extracted_wave_properties['amplitude']
    spatial_decay = extracted_wave_properties['spatial_decay']
    spatial_decay_fit = extracted_wave_properties['spatial_decay_fit']
    avg_waveform = extracted_wave_properties['avg_waveform']
    avg_waveform_per_tp = extracted_wave_properties['avg_waveform_per_tp']
    avg_centroid = extracted_wave_properties['avg_centroid']

    #These scores are NOT effected by the drift correction
    amp_score = mf.get_simple_metric(amplitude)
    spatial_decay_score = mf.get_simple_metric(spatial_decay)
    spatial_decay_fit_score = mf.get_simple_metric(spatial_decay_fit, outlier = True)
    wave_corr_score = mf.get_wave_corr(avg_waveform, param)
    wave_mse_score = mf.get_waveforms_mse(avg_waveform, param)

    #effected by drift
    for i in range(niter):
        avg_waveform_per_tp_flip = mf.flip_dim(avg_waveform_per_tp, param)
        euclid_dist = mf.get_Euclidean_dist(avg_waveform_per_tp_flip, param)

        centroid_dist, centroid_var = mf.centroid_metrics(euclid_dist, param)

        euclid_dist_rc = mf.get_recentered_euclidean_dist(avg_waveform_per_tp_flip, avg_centroid, param)

        centroid_dist_recentered = mf.recentered_metrics(euclid_dist_rc)
        traj_angle_score, traj_dist_score = mf.dist_angle(avg_waveform_per_tp_flip, param)


        # Average Euc Dist
        euclid_dist = np.nanmin(euclid_dist[:,param['peak_loc'] - param['waveidx'] == 0, :,:].squeeze(), axis = 1 )

        # TotalScore
        include_these_pairs = np.argwhere( euclid_dist < param['max_dist']) #array indices of pairs to include
        include_these_pairs_idx = np.zeros_like(euclid_dist)
        include_these_pairs_idx[euclid_dist < param['max_dist']] = 1 

        # Make a dictionary of score to include
        centroid_overlord_score = (centroid_dist_recentered + centroid_var) / 2
        waveform_score = (wave_corr_score + wave_mse_score) / 2
        trajectory_score = (traj_angle_score + traj_dist_score) / 2

        scores_to_include = {'amp_score' : amp_score, 'spatial_decay_score' : spatial_decay_score, 'centroid_overlord_score' : centroid_overlord_score,
                        'centroid_dist' : centroid_dist, 'waveform_score' : waveform_score, 'trajectory_score': trajectory_score }
        if to_use is not None:
            to_exclude = [key for key in scores_to_include if key not in to_use]
            for key in to_exclude:
                scores_to_include.pop(key, None)
        
        total_score, predictors = mf.get_total_score(scores_to_include, param)

        #Initial thresholding
        if (i < niter - 1):
            #get the thershold for a match
            thrs_opt = mf.get_threshold(total_score, within_session, euclid_dist, param, is_first_pass = True)

            param['n_expected_matches'] = np.sum( (total_score > thrs_opt).astype(int))
            prior_match = 1 - ( param['n_expected_matches'] / len(include_these_pairs))
            candidate_pairs = total_score > thrs_opt

            drifts, avg_centroid, avg_waveform_per_tp = mf.drift_n_sessions(candidate_pairs, session_switch, avg_centroid, avg_waveform_per_tp, total_score, param)

    thrs_opt = mf.get_threshold(total_score, within_session, euclid_dist, param, is_first_pass = False)
    param['n_expected_matches'] = np.sum( (total_score > thrs_opt).astype(int))
    prior_match = 1 - ( param['n_expected_matches'] / len(include_these_pairs))
    thrs_opt = np.quantile(total_score[include_these_pairs_idx.astype(bool)], prior_match)
    candidate_pairs = total_score > thrs_opt

    return total_score, candidate_pairs, scores_to_include, predictors 
