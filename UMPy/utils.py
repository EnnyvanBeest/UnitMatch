# utility function for loading files etc
import numpy as np
import pandas as pd
import os
import matplotlib.pyplot as plt

def load_tsv(path):
     """
     Loads a tsv, as a numpy array with the headers removed.
     """
     df  = pd.read_csv(path, sep='\t', skiprows = 0)
     return df.values

def get_session_data(n_units_perday):
    """
    Input the number of units per day/session as a numpy array, will return:
    the total number of units, sessionid and array where each unit is given a number according to what session it is a member of
    the index's of when the session switches in form [0, end of session 1, end of session 2....end of final session]
    """  
    n_days = len(n_units_perday)
    n_units = n_units_perday.sum()

    sessionid = np.zeros(n_units)
    SessionSwitch = np.cumsum(n_units_perday)
    SessionSwitch = np.insert(SessionSwitch, 0, 0)
    for i in range(n_days):
        sessionid[SessionSwitch[i]:SessionSwitch[i+1]] = i

    return n_units, sessionid, SessionSwitch, n_days

def get_within_session(sessionid, param):
    """
    Uses the session id to great a n_units * n_units array, where it is 0 if the units are from the same session
    and it is one if the units are from a different session
    """
    n_units = param['n_units']

    tmp1 = np.expand_dims(sessionid , axis=1)
    tmp2 = np.expand_dims(sessionid, axis=0)

    WithinSession = np.ones((n_units, n_units))
    WithinSession[tmp1 == tmp2] = 0

    return WithinSession

def get_default_param(param = None):
    """
    Create param, a dictionary with the default parameters.
    If a dictionary is given, it will add values to it without overwriting existing values.
    Do not need to give a dictionary.
    """
    tmp = {'SpikeWidth' : 82, 'waveidx' : np.arange(33,56), 'ChannelRadius' : 150,
         'PeakLoc' : 40, 'MaxDist' : 100, 'NeighbourDist' : 50, 'stepsz' : 0.01, 
         'SmoothProb' : 9, 'MinAngleDist' : 0.1
        }
    tmp['ScoreVector'] = np.arange(tmp['stepsz']/2 ,1 ,tmp['stepsz'])
    tmp['Bins'] = np.arange(0, 1 + tmp['stepsz'], tmp['stepsz'])

    #if no dictionary is given just returns the default parameters
    if param == None:
        out = tmp
    else:    
        #Add default parameters to param dictionary, does not overwrite pre existing param values
        out = tmp | param
    return out

def load_good_waveforms(WaveF_paths, unit_label_paths, param):
    """"
    This is the reccomeneded way to read in data. It uses 
    """
    if len(WaveF_paths) == len(unit_label_paths):
        n_sessions = len(WaveF_paths)
    else:
        print('Warning: gave different number of paths for waveforms and labels!')
        return

    GoodUnits = []
    for i in range(len(unit_label_paths)):
        unit_label = load_tsv(unit_label_paths[i])
        tmp_idx = np.argwhere(unit_label[:,1] == 'good')
        goodunit_idx = unit_label[tmp_idx, 0]
        GoodUnits.append(goodunit_idx)

    waveforms = []
    #go through each session and load in units to waveforms list
    for ls in range(len(WaveF_paths)):
        #load in the first good unit, to get the shape of each waveform
        tmp = np.load(WaveF_paths[ls] + rf'\Unit{int(GoodUnits[ls][0].squeeze())}_RawSpikes.npy')
        tmp_waveform = np.zeros( (len(GoodUnits[ls]), tmp.shape[0], tmp.shape[1], tmp.shape[2]))

        for i in range(len(GoodUnits[ls])):
            #loads in all GoodUnits for that session
            path = WaveF_paths[ls] + rf'\Unit{int(GoodUnits[ls][i].squeeze())}_RawSpikes.npy'
            tmp_waveform[i] = np.load(path)
        #adds that session to the list
        waveforms.append(tmp_waveform)

    del tmp_waveform
    del tmp

    n_units_perday = np.zeros(n_sessions, dtype = 'int')
    waveform = np.array([])

    #add all of the individual waveforms to one waveform array
    for i in range(n_sessions):
        if i == 0:
            waveform = waveforms[i] 
        else:
            waveform = np.concatenate((waveform, waveforms[i]), axis = 0)

        n_units_perday[i] = waveforms[i].shape[0]

    param['n_units'], sessionid, SessionSwitch, param['n_days'] = get_session_data(n_units_perday)
    WithinSession = get_within_session(sessionid, param)
    param['n_channels'] = waveform.shape[2]
    return waveform, sessionid, SessionSwitch, WithinSession, param

def compare_units(WeightedAvgWaveF, AvgCentroid, unit1, unit2):
    """
    Basic helper function, plots the average wavefucntion (of cv 0) and the average centroid to quickly compare 2 units 
    """
    plt.plot(WeightedAvgWaveF[:,unit1,0])
    plt.plot(WeightedAvgWaveF[:,unit2,0])
    print(f'Average centroid of unit {unit1} is :{AvgCentroid[:,unit1,0]}')
    print(f'Average centroid of unit {unit2} is :{AvgCentroid[:,unit2,0]}')


def evaluate_output(output, param, WithinSession, SessionSwitch, match_threshold = 0.5):
    """"
    Input: output - the n_units * n_units probability matrix (each value is prob those units match)
    the param dictionary and optionally the threshold used to calculate if a unit is a match

    This function then print:
    The number of units matched to themselves across cv
    The false negative %, how many did not match to themselves across cv
    the false positive % in two ways, how many miss-matches are there in the off-diagonal per session
    and how many  false match out of how many matches we should get
    """

    thrs_output = np.zeros_like(output)
    thrs_output[output > match_threshold] = 1

    # get the number of diagonal matches
    UM_diag = np.sum(thrs_output[np.eye(param['n_units']).astype(bool)])
    self_match = UM_diag / param['n_units'] *100
    print(f'The percentage of units matched to themselves is: {self_match}%')
    print(f'The percentage of false -ve\'s then is: {100 - self_match}% \n')

    #off-diagonal miss-matches
    UM_offdiag = np.zeros_like(output)
    UM_offdiag = thrs_output
    UM_offdiag[WithinSession == 1] = 0 
    UM_offdiag[np.eye(param['n_units']) == 1] = 0 
    FP_est =  UM_offdiag.sum() / (param['n_units']) 
    print(f'The rate of miss-match(es) per expected match {FP_est}')


    #compute matlab FP per session per session
    FP_est_pd = np.zeros(param['n_days'])
    for did in range(param['n_days']):
        tmp_diag = thrs_output[SessionSwitch[did]:SessionSwitch[did + 1], SessionSwitch[did]:SessionSwitch[did + 1]]
        n_units = tmp_diag.shape[0]
        tmp_diag[np.eye(n_units) == 1] = 0 
        FP_est_pd[did] = tmp_diag.sum() / (n_units ** 2 - n_units) * 100
        print(f'The percentage of false +ve\'s is {FP_est_pd[did]}% for session {did +1}')

    print('\nThis assumes that the spike sorter has made no mistakes')


##########################################################################################################################
#The folowing functions are the old way of reading in units, is slower and will not work if unit are missing e.g 1,2,4

def load_waveforms(WaveF_paths, unit_label_paths, param):
    """
    This function uses a list of paths to the average waveforms and good units to load in all
    the waveforms and session related information.
    """

    #assuming the number of sessions is the same as length of WaveF_paths
    n_sessions = len(WaveF_paths)

    # load in individual session waveforms as list of np arrays
    waveforms = []
    for i in range(len(WaveF_paths)):
        #waveforms.append(util.get_waveform(WaveF_paths[i]))
        tmp = get_waveform(WaveF_paths[i])
        goodunit_idxtmp = get_good_unit_idx(unit_label_paths[i])
        waveforms.append(good_units(tmp, goodunit_idxtmp))
        del tmp
        del goodunit_idxtmp


    n_units_perday = np.zeros(n_sessions, dtype = 'int')
    waveform = np.array([])

    #add all of the individual waveforms to one waveform array
    for i in range(n_sessions):
        if i == 0:
            waveform = waveforms[i] 
        else:
            waveform = np.concatenate((waveform, waveforms[i]), axis = 0)

        n_units_perday[i] = waveforms[i].shape[0]

    param['n_units'], sessionid, SessionSwitch, param['n_days'] = get_session_data(n_units_perday)
    WithinSession = get_within_session(sessionid, param)
    param['n_channels'] = waveform.shape[2]

    return waveform, sessionid, SessionSwitch, WithinSession, param

def get_waveform(folder_path):
    '''
    Assuming the raw spike are saved as Unitxxx_RawSpikes.npy where xxx is the number id of the spike, 
    Requires the path to the folder where all the spike are saved.
    requires all spike to have same dimensions

    returns the waveform matrix (No. Units, spike dims), assumed (No.units, time, channel_no, first half/second half) 

    Could:
    - parallelize
    - open in blocks of n units
    - adapt to open other types of files
    '''    
    num_files = len(os.listdir(folder_path))

    tmp = np.load(folder_path + r'\Unit0_RawSpikes.npy')
    waveform = np.zeros((num_files, tmp.shape[0], tmp.shape[1], tmp.shape[2]))

    for i in range(num_files):
        path = folder_path + rf'\Unit{i}_RawSpikes.npy'
        waveform[i] = np.load(path)
    return waveform


def get_files_outdated(folder_path):
    ''' Similar to get_wave for (maybe more optimized), however:
    1 doesn't return numerically order spikes
    2 doesn't assume anything about how spikes are saved ( except numpy array of constant dims)
    '''
    for path, dirs, files in os.walk(folder_path, topdown=True):
        i = 0
        tmp = np.load(os.path.join(path, files[0]))
        waveform = np.zeros((len(files),tmp.shape[0], tmp.shape[1], tmp.shape[2] ))
        for f in files:
           waveform[i] =  np.load(os.path.join(path, f))
           i +=1 
    return waveform 

def get_good_unit_idx(unit_label_path):
    """ 
    Assuming until label path, is the path to a tsv file where the second row onwards is 2 columns, where the second one is the unit label
    """
    unit_label = load_tsv(unit_label_path)
    tmp_idx = np.argwhere(unit_label[:,1] == 'good')
    goodunit_idx = unit_label[tmp_idx, 0]
    return goodunit_idx

def good_units(waveform, goodunit_idx):
    """
    Using goodunit_idx, this function returns the good units of a waveform
    ** may want to edit, so it can select good units if the unit axes isn't the first axis and is adaptable to any shape of input 
    """
    waveform = waveform[goodunit_idx,:,:,:].squeeze()
    return waveform