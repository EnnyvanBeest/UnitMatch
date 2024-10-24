import numpy as np

def get_default_param(param = None):
    """
    Creates the default param dictionary for UnitMatch.

    Parameters
    ----------
    param : dict, optional
        A given dictionary which will take priority over default params, by default None

    Returns
    -------
    dict
        The param dictionary
    """
    tmp = {'spike_width' : 82, # will automatically update after loading in the waveform
            'waveidx' : np.arange(33,56), # what idxs cover the spike
            'channel_radius' : 150, # the max distance channels which a unit will consider
            'peak_loc' : 40, # were the spike peak will be aligned to
            'max_dist' : 100, # max distance for possible matches
            'neighbour_dist' : 50,
            'stepsz' : 0.01, # size of of a step in probability distributions
            'smooth_prob' : 9, # probability smoothing size
            'min_angle_dist' : 0.1, # smallest distance for and angle to be consider
            'min_new_shank_distance' : 100, #The smallest distance which separates 2 shanks
            'units_per_shank_thrs' : 15, # threshold for doing per shank drift correction
            'match_threshold' : 0.5 # probability threshold to consider as a match
        }
    
    tmp['score_vector'] = np.arange(tmp['stepsz']/2 ,1 ,tmp['stepsz'])
    tmp['bins'] = np.arange(0, 1 + tmp['stepsz'], tmp['stepsz'])

    #if no dictionary is given just returns the default parameters
    if param == None:
        out = tmp
    else:
        #Add default parameters to param dictionary, does not overwrite pre existing param values
        out = tmp | param
    return out
