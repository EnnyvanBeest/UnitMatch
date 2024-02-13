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

def get_session_number(unitid, SessionSwitch):

    for i in range(len(SessionSwitch) - 1):
        if (SessionSwitch[i] <= unitid < SessionSwitch[i+1]):
            return i

def get_session_data(nUnitsPerSession):
    """
    Input the number of units per day/session as a numpy array, will return:
    the total number of units, sessionid and array where each unit is given a number according to what session it is a member of
    the index's of when the session switches in form [0, end of session 1, end of session 2....end of final session]
    """  
    nSessions = len(nUnitsPerSession)                  
    nUnits = nUnitsPerSession.sum()

    sessionid = np.zeros(nUnits, dtype = int)
    SessionSwitch = np.cumsum(nUnitsPerSession)
    SessionSwitch = np.insert(SessionSwitch, 0, 0)
    for i in range(nSessions):
        sessionid[SessionSwitch[i]:SessionSwitch[i+1]] = int(i)

    return nUnits, sessionid, SessionSwitch, nSessions

def get_within_session(sessionid, param):
    """
    Uses the session id to great a nUnits * nUnits array, where it is 0 if the units are from the same session
    and it is one if the units are from a different session
    """
    nUnits = param['nUnits']

    tmp1 = np.expand_dims(sessionid , axis=1)
    tmp2 = np.expand_dims(sessionid, axis=0)

    WithinSession = np.ones((nUnits, nUnits))
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
         'SmoothProb' : 9, 'MinAngleDist' : 0.1, 'NoShanks' : 4, 'ShankDist' : 175,
         'MatchNumThreshold' : 15
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

def load_good_waveforms(WavePaths, UnitLabelPaths, param):
    """"
    This is the reccomeneded way to read in data. It uses 
    """
    if len(WavePaths) == len(UnitLabelPaths):
        nSessions = len(WavePaths)
    else:
        print('Warning: gave different number of paths for waveforms and labels!')
        return

    GoodUnits = []
    for i in range(len(UnitLabelPaths)):
        UnitLabel = load_tsv(UnitLabelPaths[i])
        TmpIdx = np.argwhere(UnitLabel[:,1] == 'good')
        goodunit_idx = UnitLabel[TmpIdx, 0]
        GoodUnits.append(goodunit_idx)

    waveforms = []
    #go through each session and load in units to waveforms list
    for ls in range(len(WavePaths)):
        #load in the first good unit, to get the shape of each waveform
        tmp = np.load(WavePaths[ls] + rf'\Unit{int(GoodUnits[ls][0].squeeze())}_RawSpikes.npy')
        tmpWaveform = np.zeros( (len(GoodUnits[ls]), tmp.shape[0], tmp.shape[1], tmp.shape[2]))

        for i in range(len(GoodUnits[ls])):
            #loads in all GoodUnits for that session
            path = WavePaths[ls] + rf'\Unit{int(GoodUnits[ls][i].squeeze())}_RawSpikes.npy'
            tmpWaveform[i] = np.load(path)
        #adds that session to the list
        waveforms.append(tmpWaveform)

    del tmpWaveform
    del tmp

    nUnitsPerSession = np.zeros(nSessions, dtype = 'int')
    waveform = np.array([])

    #add all of the individual waveforms to one waveform array
    for i in range(nSessions):
        if i == 0:
            waveform = waveforms[i] 
        else:
            waveform = np.concatenate((waveform, waveforms[i]), axis = 0)

        nUnitsPerSession[i] = waveforms[i].shape[0]

    param['nUnits'], sessionid, SessionSwitch, param['nSessions'] = get_session_data(nUnitsPerSession)
    WithinSession = get_within_session(sessionid, param)
    param['nChannels'] = waveform.shape[2]
    return waveform, sessionid, SessionSwitch, WithinSession, GoodUnits, param

def get_good_units(UnitLabelPaths, good = True):
    """
    Requires the paths to .tsv files, which contain the unit index's and if they area a good unit.
    Will return a list where each index of the list is a numpy array ofall the good index's.
    This function is set to only get index's for units labelled 'good', pass good = False to get ALL unit index's
    """
    GoodUnits = []
    for i in range(len(UnitLabelPaths)):
        UnitLabel = load_tsv(UnitLabelPaths[i])
        if good == True:
            TmpIdx = np.argwhere(UnitLabel[:,1] == 'good')
        else:
            TmpIdx = UnitLabel[:,0] # every unit index in the first column
        GoodUnitIdx = UnitLabel[TmpIdx, 0]
        GoodUnits.append(GoodUnitIdx)
    return GoodUnits

def load_good_units(GoodUnits, WavePaths, param):
    """
    Requires a list which contains a numpy array with the units to load per session, as well as a path to
    a file which contains all the the raw averaged units 
    """
    if len(WavePaths) == len(GoodUnits):
        nSessions = len(WavePaths)
    else:
        print('Warning: gave different number of paths for waveforms and labels!')
        return
    
    waveforms = []
    #go through each session and load in units to waveforms list
    for ls in range(len(WavePaths)):
        #load in the first good unit, to get the shape of each waveform
        tmp = np.load(WavePaths[ls] + rf'\Unit{int(GoodUnits[ls][0].squeeze())}_RawSpikes.npy')
        tmpWaveform = np.zeros( (len(GoodUnits[ls]), tmp.shape[0], tmp.shape[1], tmp.shape[2]))

        for i in range(len(GoodUnits[ls])):
            #loads in all GoodUnits for that session
            path = WavePaths[ls] + rf'\Unit{int(GoodUnits[ls][i].squeeze())}_RawSpikes.npy'
            tmpWaveform[i] = np.load(path)
        #adds that session to the list
        waveforms.append(tmpWaveform)

    del tmpWaveform
    del tmp

    nUnitsPerSession = np.zeros(nSessions, dtype = 'int')
    waveform = np.array([])

    #add all of the individual waveforms to one waveform array
    for i in range(nSessions):
        if i == 0:
            waveform = waveforms[i] 
        else:
            waveform = np.concatenate((waveform, waveforms[i]), axis = 0)

        nUnitsPerSession[i] = waveforms[i].shape[0]

    param['nUnits'], sessionid, SessionSwitch, param['nSessions'] = get_session_data(nUnitsPerSession)
    WithinSession = get_within_session(sessionid, param)
    param['nChannels'] = waveform.shape[2]
    return waveform, sessionid, SessionSwitch, WithinSession, param


def compare_units(AvgWaveform, AvgCentroid, unit1, unit2):
    """
    Basic helper function, plots the average wavefucntion (of cv 0) and the average centroid to quickly compare 2 units 
    """
    plt.plot(AvgWaveform[:,unit1,0])
    plt.plot(AvgWaveform[:,unit2,0])
    print(f'Average centroid of unit {unit1} is :{AvgCentroid[:,unit1,0]}')
    print(f'Average centroid of unit {unit2} is :{AvgCentroid[:,unit2,0]}')


def evaluate_output(output, param, WithinSession, SessionSwitch, MatchThreshold = 0.5):
    """"
    Input: output - the n_units * n_units probability matrix (each value is prob those units match)
    the param dictionary and optionally the threshold used to calculate if a unit is a match

    This function then print:
    The number of units matched to themselves across cv
    The false negative %, how many did not match to themselves across cv
    the false positive % in two ways, how many miss-matches are there in the off-diagonal per session
    and how many  false match out of how many matches we should get
    """

    OutputThreshold = np.zeros_like(output)
    OutputThreshold[output > MatchThreshold] = 1

    # get the number of diagonal matches
    nDiag = np.sum(OutputThreshold[np.eye(param['nUnits']).astype(bool)])
    SelfMatch = nDiag / param['nUnits'] *100
    print(f'The percentage of units matched to themselves is: {SelfMatch}%')
    print(f'The percentage of false -ve\'s then is: {100 - SelfMatch}% \n')

    #off-diagonal miss-matches
    nOffDiag = np.zeros_like(output)
    nOffDiag = OutputThreshold
    nOffDiag[WithinSession == 1] = 0 
    nOffDiag[np.eye(param['nUnits']) == 1] = 0 
    FPest =  nOffDiag.sum() / (param['nUnits']) 
    print(f'The rate of miss-match(es) per expected match {FPest}')


    #compute matlab FP per session per session
    FPestPerSession = np.zeros(param['nSessions'])
    for did in range(param['nSessions']):
        tmpDiag = OutputThreshold[SessionSwitch[did]:SessionSwitch[did + 1], SessionSwitch[did]:SessionSwitch[did + 1]]
        nUnits = tmpDiag.shape[0]
        tmpDiag[np.eye(nUnits) == 1] = 0 
        FPestPerSession[did] = tmpDiag.sum() / (nUnits ** 2 - nUnits) * 100
        print(f'The percentage of false +ve\'s is {FPestPerSession[did]}% for session {did +1}')

    print('\nThis assumes that the spike sorter has made no mistakes')


##########################################################################################################################
#The folowing functions are the old way of reading in units, is slower and will not work if unit are missing e.g 1,2,4

# def load_waveforms(WavePaths, UnitLabelPaths, param):
#     """
#     This function uses a list of paths to the average waveforms and good units to load in all
#     the waveforms and session related information.
#     """

#     #assuming the number of sessions is the same as length of WaveF_paths
#     nSessions = len(WavePaths)

#     # load in individual session waveforms as list of np arrays
#     waveforms = []
#     for i in range(len(WavePaths)):
#         #waveforms.append(util.get_waveform(WaveF_paths[i]))
#         tmp = get_waveform(WavePaths[i])
#         GoodUnitIdxTmp = get_good_unit_idx(UnitLabelPaths[i])
#         waveforms.append(good_units(tmp, GoodUnitIdxTmp))
#         del tmp
#         del GoodUnitIdxTmp


#     nUnitsPerSession = np.zeros(nSessions, dtype = 'int')
#     waveform = np.array([])

#     #add all of the individual waveforms to one waveform array
#     for i in range(nSessions):
#         if i == 0:
#             waveform = waveforms[i] 
#         else:
#             waveform = np.concatenate((waveform, waveforms[i]), axis = 0)

#         nUnitsPerSession[i] = waveforms[i].shape[0]

#     param['n_units'], sessionid, SessionSwitch, param['n_days'] = get_session_data(nUnitsPerSession)
#     WithinSession = get_within_session(sessionid, param)
#     param['n_channels'] = waveform.shape[2]

#     return waveform, sessionid, SessionSwitch, WithinSession, param

# def get_waveform(FolderPath):
#     '''
#     Assuming the raw spike are saved as Unitxxx_RawSpikes.npy where xxx is the number id of the spike, 
#     Requires the path to the folder where all the spike are saved.
#     requires all spike to have same dimensions

#     returns the waveform matrix (No. Units, spike dims), assumed (No.units, time, channel_no, first half/second half) 

#     Could:
#     - parallelize
#     - open in blocks of n units
#     - adapt to open other types of files
#     '''    
#     nFiles = len(os.listdir(FolderPath))

#     tmp = np.load(FolderPath + r'\Unit0_RawSpikes.npy')
#     waveform = np.zeros((nFiles, tmp.shape[0], tmp.shape[1], tmp.shape[2]))

#     for i in range(nFiles):
#         path = FolderPath + rf'\Unit{i}_RawSpikes.npy'
#         waveform[i] = np.load(path)
#     return waveform


# def get_files_outdated(FolderPath):
#     ''' Similar to get_wave for (maybe more optimized), however:
#     1 doesn't return numerically order spikes
#     2 doesn't assume anything about how spikes are saved ( except numpy array of constant dims)
#     '''
#     for path, dirs, files in os.walk(FolderPath, topdown=True):
#         i = 0
#         tmp = np.load(os.path.join(path, files[0]))
#         waveform = np.zeros((len(files),tmp.shape[0], tmp.shape[1], tmp.shape[2] ))
#         for f in files:
#            waveform[i] =  np.load(os.path.join(path, f))
#            i +=1 
#     return waveform 

# def get_good_unit_idx(UnitLabelPath):
#     """ 
#     Assuming until label path, is the path to a tsv file where the second row onwards is 2 columns, where the second one is the unit label
#     """
#     UnitLabel = load_tsv(UnitLabelPath)
#     TmpIdx = np.argwhere(UnitLabel[:,1] == 'good')
#     GoodUnitIdx = UnitLabel[TmpIdx, 0]
#     return GoodUnitIdx

# def good_units(waveform, GoodUnitIdx):
#     """
#     Using goodunit_idx, this function returns the good units of a waveform
#     ** may want to edit, so it can select good units if the unit axes isn't the first axis and is adaptable to any shape of input 
#     """
#     waveform = waveform[GoodUnitIdx,:,:,:].squeeze()
#     return waveform
    
# #################################
# # The following is how the above function would be used to read in data
# # read in data and select the good units and exact metadata

# #loads in waveforms, 
# waveform1 = util.get_waveform(WavePath1)
# waveform2 = util.get_waveform(WavePath2)

# #selects 'good' units 
# GoodUnitIdx1 = util.get_good_unit_idx(UnitLabelPath1)
# waveform1 = util.good_units(waveform1, GoodUnitIdx1)
# GoodUnitIdx2 = util.get_good_unit_idx(UnitLabelPath2)
# waveform2 = util.good_units(waveform2, GoodUnitIdx2)

# #joins the waveforms together, and keep track of length of each session
# waveform = np.concatenate((waveform1,waveform2), axis = 0 )
# nUnitsPerSession = np.asarray([waveform1.shape[0], waveform2.shape[0]])

# # assigns a session id to each units and notes when the sessions switch
# param['nUnits'], sessionid, SessionSwitch, param['nSessions'] = util.get_session_data(nUnitsPerSession)
# WithinSession = util.get_within_session(sessionid, param)
# param['nChannels'] = waveform.shape[2]
