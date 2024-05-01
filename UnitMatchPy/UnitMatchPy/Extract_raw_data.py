#Function for extracting and averaging raw data

import os
import numpy as np
from pathlib import Path
from scipy.ndimage import gaussian_filter
from mtscomp import decompress
from joblib import Parallel, delayed
import UnitMatchPy.utils as util

#Decompressed data functions
def Read_Meta(metaPath):
    "Readin Meta data as a dictionary"
    metaDict = {}
    with metaPath.open() as f:
        mdatList = f.read().splitlines()
        # convert the list entries into key value pairs
        for m in mdatList:
            csList = m.split(sep='=')
            if csList[0][0] == '~':
                currKey = csList[0][1:len(csList[0])]
            else:
                currKey = csList[0]
            metaDict.update({currKey: csList[1]})

    return(metaDict)

def get_sample_idx(SpikeTimes, UnitIDs, SampleAmount, units):
    """
    Needs spike times, unit ID's (from kilosort dir) and maximum number of samples per unit.
    Returns a (nUnits, SampleAmount) array with what spikes to sample for every unit, selected spikes evely spaced over time and
    fill with NaN if the unit has less spikes than SampleAmount
    """

    UniqueUnitIDs = np.unique(UnitIDs) 
    nUnitsALL = len(UniqueUnitIDs)

    SampleIdx = np.zeros((nUnitsALL, SampleAmount))
    #Process ALL unit
    for i, idx in enumerate(units):
        UnitTimes = SpikeTimes[UnitIDs == idx]
        if SampleAmount < len(UnitTimes):
            ChooseIdx = np.linspace(0,len(UnitTimes)-1, SampleAmount, dtype = int) # -1 so can't indx out of region
            SampleIdx[i,:] = UnitTimes[ChooseIdx]
        else:
            SampleIdx[i,:len(UnitTimes)] = UnitTimes
            SampleIdx[i,len(UnitTimes):] = np.nan
    
    return SampleIdx

def Extract_A_Unit(SampleIdx, Data, HalfWidth, SpikeWidth, nChannels, SampleAmount):
    """ 
    This function extracts and averages the raw data for A unit, and splits the unit into two half for cross verification.
    returns AvgWavforms shape (nChannels, SpikeWidth, 2)

    NOTE - Here SampleIdx is a array of shape (SampleAmount), i.e use SampleIdx[UnitIdx] to get the AvgWAveform for that unit
    """

    Channels = np.arange(0,nChannels)

    AllSampleWaveforms = np.zeros( (SampleAmount, SpikeWidth, nChannels))
    for i, idx in enumerate(SampleIdx[:]):
        if np.isnan(idx):
            continue 
        tmp = Data[ int(idx - HalfWidth - 1): int(idx + HalfWidth - 1), Channels] # -1, to better fit with ML
        tmp.astype(np.float32)
        #gaussina smooth, over time gaussina window = 5, sigma = window size / 5
        tmp = gaussian_filter(tmp, 1, radius = 2, axes = 0) #edges are handled differently to ML
        # window ~ radius *2 + 1
        tmp = tmp - np.mean(tmp[:20,:], axis = 0)
        AllSampleWaveforms[i] = tmp

    #median and split CV's
    nWavs = np.sum(~np.isnan(SampleIdx[:]))
    CVlim = np.floor(nWavs / 2).astype(int)

    #find median over samples
    AvgWaveforms = np.zeros((SpikeWidth, nChannels, 2))
    AvgWaveforms[:, :, 0] = np.median(AllSampleWaveforms[:CVlim, :, :], axis = 0) #median over samples
    AvgWaveforms[:, :, 1] = np.median(AllSampleWaveforms[CVlim:nWavs, :, :], axis = 0) #median over samples
    return AvgWaveforms

def Extract_A_UnitKS4(SampleIdx, Data, SamplesBefore, SamplesAfter, SpikeWidth, nChannels, SampleAmount):
    """ 
    This function extracts and averages the raw data for A unit, and splits the unit into two half for cross verification.
    returns AvgWavforms shape (nChannels, SpikeWidth, 2)

    NOTE - Here SampleIdx is a array of shape (SampleAmount), i.e use SampleIdx[UnitIdx] to get the AvgWAveform for that unit
    """

    Channels = np.arange(0,nChannels)

    AllSampleWaveforms = np.zeros( (SampleAmount, SpikeWidth, nChannels))
    for i, idx in enumerate(SampleIdx[:]):
        if np.isnan(idx):
            continue 
        tmp = Data[ int(idx - SamplesBefore - 1): int(idx + SamplesAfter - 1), Channels] # -1, to better fit with ML
        tmp.astype(np.float32)
        #gaussina smooth, over time gaussina window = 5, sigma = window size / 5
        tmp = gaussian_filter(tmp, 1, radius = 2, axes = 0) #edges are handled differently to ML
        # window ~ radius *2 + 1
        tmp = tmp - np.mean(tmp[:20,:], axis = 0)
        AllSampleWaveforms[i] = tmp

    #median and split CV's
    nWavs = np.sum(~np.isnan(SampleIdx[:]))
    CVlim = np.floor(nWavs / 2).astype(int)

    #find median over samples
    AvgWaveforms = np.zeros((SpikeWidth, nChannels, 2))
    AvgWaveforms[:, :, 0] = np.median(AllSampleWaveforms[:CVlim, :, :], axis = 0) #median over samples
    AvgWaveforms[:, :, 1] = np.median(AllSampleWaveforms[CVlim:nWavs, :, :], axis = 0) #median over samples
    return AvgWaveforms


def Save_AvgWaveforms(AvgWaveforms, SaveDir, GoodUnits, ExtractGoodUnitsOnly = False):
    """
    Will save the extracted average waveforms in a folder called 'RawWaveforms' in the given SaveDir
    Each waveform will be saved in a unique .npy file called 'UnitX_RawSpikes.npy.
    Supply GoodUnits, a array of which idx's are included, if they you are not extract all units
    from the recording session. 
    """
    CurrentDir = os.getcwd()
    os.chdir(SaveDir)
    DirList = os.listdir()
    if 'RawWaveforms' in DirList:
        TmpPath = os.path.join(SaveDir, 'RawWaveforms')
        
    else:
        os.mkdir('RawWaveforms')
        TmpPath = os.path.join(SaveDir, 'RawWaveforms')

    os.chdir(TmpPath)

    #first axis is each unit

    #ALL waveforms from 0->nUnits
    if ExtractGoodUnitsOnly == False:
        for i in range(AvgWaveforms.shape[0]):
            np.save(f'Unit{i}_RawSpikes.npy', AvgWaveforms[i,:,:,:])

    #If only extracting GoodUnits
    else:
        for i, idx in enumerate(GoodUnits):
            # ironically need idx[0], to selct value so saves with correct name
            np.save(f'Unit{idx[0]}_RawSpikes.npy', AvgWaveforms[i,:,:,:])
    
    os.chdir(CurrentDir)




# Load in necessary files from KS directory and raw data directory   
# extracting n Sessions
def get_raw_data_paths(RawDataDirPaths):
    """
    This function requires RawDatPaths, a list of pahts to the Raw data directories, e.g where .cbin, .ch .meta files are
    This function will return a list fo paths to the.cbin, .ch and .meta files
    """
    cbinPaths = []
    chPaths = []
    metaPaths = []

    for i in range(len(RawDataDirPaths)):
        for f in os.listdir(RawDataDirPaths[i]):
            name, ext = os.path.splitext(f)
            
            if ext == '.cbin':
                cbinPaths.append(os.path.join(RawDataDirPaths[i], name + ext))

            if ext == '.ch':
                chPaths.append(os.path.join(RawDataDirPaths[i], name + ext))

            if ext == '.meta':
                metaPaths.append(os.path.join(RawDataDirPaths[i], name + ext))
    
    return cbinPaths, chPaths, metaPaths


def extract_KSdata(KSdirs, ExtractGoodUnitsOnly = False):
    """
    This fucntion requires KSdirs, a lsit of KiloSort directories for each session.
    This function will then load in the spike_times, spike_ids and a Good_Units
    """
    nSessions = len(KSdirs)

    #Load Spike Times
    SpikeTimes = []
    for i in range(nSessions):
        PathTmp = os.path.join(KSdirs[i], 'spike_times.npy')
        SpikeTimestmp = np.load(PathTmp)
        SpikeTimes.append(SpikeTimestmp)

    
    #Load Spike ID's
    SpikeIDs = []
    for i in range(nSessions):
        PathTmp = os.path.join(KSdirs[i], 'spike_clusters.npy')
        SpikeIDstmp = np.load(PathTmp)
        SpikeIDs.append(SpikeIDstmp)


    if ExtractGoodUnitsOnly:
        #Good unit ID's
        UnitLabelPaths = []

        # load Good unit Paths
        for i in range(nSessions):
            UnitLabelPaths.append( os.path.join(KSdirs[i], 'cluster_group.tsv'))

        GoodUnits = util.get_good_units(UnitLabelPaths)

        return SpikeIDs, SpikeTimes, GoodUnits
    
    else:
        return SpikeIDs, SpikeTimes, [None for s in range(nSessions)]
