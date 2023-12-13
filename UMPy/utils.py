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


def get_waveform(folder_path):
    '''
    Asuming the raw spike are saved as Unitxxx_RawSpikes.npy where xxx is the number id of the spike, 
    Requires the path to the folder were all the spike are saved.
    requires all spike to have same dimensions

    returns the waveform matrix (No. Units, spike dims), assumed (No.units, time, channel_no, firsthalf/secondhalf) 

    Could:
    - parallelise
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


def get_files_outdated(folderpath):
    ''' Similar to get_wavefor (maybe more optimsied), however:
    1 doesn return numerically order spikes
    2 doesnt assume anything about how spikes are saved ( excpet numpy array of constant dims)
    '''
    for path, dirs, files in os.walk(folderpath, topdown=True):
        i = 0
        tmp = np.load(os.path.join(path, files[0]))
        waveform = np.zeros((len(files),tmp.shape[0], tmp.shape[1], tmp.shape[2] ))
        for f in files:
           waveform[i] =  np.load(os.path.join(path, f))
           i +=1 
    return waveform 

def get_good_unit_idx(unit_label_path):
    """ Assuming until label path, is the path to a tsv file where the second rown onwardsits 2 coloumsn, where the second one is the unit label
    """
    unit_label = load_tsv(unit_label_path)
    goodunit_idx = np.argwhere(unit_label[:,1] == 'good')
    return goodunit_idx

def good_units(waveform, goodunit_idx):
    """
    Usinggoodunit_idx, this function returns the good units of a waveform
    ** may want to edit, so it can slect good units if the unit axes isnt the first axis and is adaptable to any shape of input 
    """
    waveform = waveform[goodunit_idx,:,:,:].squeeze()
    return waveform


def get_session_data(n_units_perday):
    """
    Input the number of units per day/session as a numpy array, will reuturn:
    the total number of units, sessionid and array where each unit is given a number according to what session it is a member of
    the indexs of when the session switches in form [0, end of session 1, end of session 2....end of final session]
    """  
    n_days = len(n_units_perday)
    n_units = n_units_perday.sum()

    sessionid = np.zeros(n_units)
    sessionswitch = np.cumsum(n_units_perday)
    sessionswitch = np.insert(sessionswitch, 0, 0)
    for i in range(n_days):
        sessionid[sessionswitch[i]:sessionswitch[i+1]] = i

    return n_units, sessionid, sessionswitch, n_days

def get_within_session(sessionid, param):
    """
    Uses the session id to great a n_units * n_units array, where it is 0 if the units are from the same session
    and it is one if the units are from a different session
    """
    n_units = param['n_units']

    tmp1 = np.expand_dims(sessionid , axis=1)
    tmp2 = np.expand_dims(sessionid, axis=0)

    within_session = np.ones((n_units, n_units))
    within_session[tmp1 == tmp2] = 0

    return within_session

def get_default_param():
    """
    Create param, a dictionary with the default parameters.
    """
    param = {'spike_width' : 82, 'waveidx' : np.arange(33,56), 'channel_radius' : 150,
         'peak_loc' : 40, 'maxdist' : 100, 'NeighbourDist' : 50, 'stepsz' : 0.01, 
         'SmoothProb' : 9, 
        }
    param['ScoreVector'] = np.arange(param['stepsz']/2 ,1 ,param['stepsz'])
    param['Bins'] = np.arange(0, 1 + param['stepsz'], param['stepsz'])


    return param

def load_waveforms(WaveF_paths, unit_label_paths, param):
    """
    This function uses a list of paths to the average waveforms and good units to load in all
    the waveforms and seesion related infomation.
    """

    #assuming unit length paths same length as wavef path
    n_sessions = len(WaveF_paths)

    # load in as list of np arrays
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

    for i in range(n_sessions):
        if i == 0:
            waveform = waveforms[i] 
        else:
            waveform = np.concatenate((waveform, waveforms[i]), axis = 0)

        n_units_perday[i] = waveforms[i].shape[0]

    param['n_units'], sessionid, sessionswitch, param['n_days'] = get_session_data(n_units_perday)
    within_session = get_within_session(sessionid, param)

    return waveform, sessionid, sessionswitch, within_session, param

def compare_units(WeightedAvgWaveF, avg_centroid, unit1, unit2):
    """
    Basic helper function, plots the average wavefucntion (of cv 0) and the average centroid to quickly comapre 2 units 
    """
    plt.plot(WeightedAvgWaveF[:,unit1,0])
    plt.plot(WeightedAvgWaveF[:,unit2,0])
    print(f'Average centroid of unit {unit1} is :{avg_centroid[:,unit1,0]}')
    print(f'Average centroid of unit {unit2} is :{avg_centroid[:,unit2,0]}')


def evaluate_output(output, param, within_session, sessionswitch, match_threshold = 0.5):
    """"
    Input: output - the n_units * n_units probabilty matrix (each value is prob those units match)
    the param dictionary and optioanly the threshold used to calculate if a unit is a match

    This function then print:
    The number of untis matched to themselves across cv
    The false negative %, how many did not mathc to themsleves across cv
    the false postive % in two ways, how many miss-matches are there in the off-diagonal per session
    and how many  false match out of how many mantches we should get
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
    UM_offdiag[within_session == 1] = 0 
    UM_offdiag[np.eye(param['n_units']) == 1] = 0 
    FP_est =  UM_offdiag.sum() / (param['n_units']) 
    print(f'The rate of miss-match(es) per expected match {FP_est}')


    #compute matlab FP per session per session
    FP_est_pd = np.zeros(param['n_days'])
    for did in range(param['n_days']):
        tmp_diag = thrs_output[sessionswitch[did]:sessionswitch[did + 1], sessionswitch[did]:sessionswitch[did + 1]]
        n_units = tmp_diag.shape[0]
        tmp_diag[np.eye(n_units) == 1] = 0 
        FP_est_pd[did] = tmp_diag.sum() / (n_units ** 2 - n_units) * 100
        print(f'The percentage of false +ve\'s is {FP_est_pd[did]}% for session {did +1}')

    print('\nThis assumes that the spike sorter has made no mistakes')

