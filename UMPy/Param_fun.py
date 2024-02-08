# This file will containa ll the necessary function for extracting waveform parameters from average waveforms

import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import detrend
import scipy as sp


import pandas as pd


def detrend_waveform(waveform):
    """
    This function accepts the raw waveforms (nUnits, SpikeWidth, nChannels, CV)
    The output is the same shape, and linearly detrended across time.
    """ 
    return detrend(waveform, axis = 2)

def get_spatialfp(waveform):
    """
    Input: waveform np array (nUnits, time, nChannels, first/second half)
    By taking the maximum value along the time axis of the absolute values
    Output: np array (Units, nChannels, first/second half)
    """
    SpatialFP = np.max(np.abs(waveform), axis = 1)
    return SpatialFP

def get_max_site(SpatialFP):
    """
    Input: SpatialFP (nUnits, nChannels, first/second half)
    By taking the index of maximum argument along the channel axis,
    Output: (nUnits, first/second half), which gives the maximum site for each unit in first/second half
    """
    MaxSite = np.argmax(SpatialFP, axis = 1)
    return MaxSite

def get_max_sites(waveform, ChannelPos, ClusInfo, param):
    """
    Using waveforms, ChannelPos and param, to find the max channel for each unit and cv, this function also
    returns good idx's / positions, by selecting channels within ChannelRadius (default 150 um) 
    """

    nUnits = param['nUnits']
    nChannels = param['nChannels']
    ChannelRadius = param['ChannelRadius']
    waveidx = param['waveidx']
    SessionID = ClusInfo['SessionID']
    
    MeanCV = np.mean(waveform, axis = 3) # average of each cv
    SpatialFootprint = get_spatialfp(MeanCV) # choose max time 
    MaxSiteMean= get_max_site(SpatialFootprint) # argument of MaxSite

    # Finds the indices where the distance from the max site mean is small
    goodidx = np.empty((nUnits, nChannels))
    for i in range(ChannelPos[0].shape[0]): #looping over each site, assuming all ChannelPos, are the same no of channels
        for j in range(nUnits):
            dist = np.linalg.norm(ChannelPos[SessionID[j]][MaxSiteMean[j],:] - ChannelPos[SessionID[j]][i,:])
            good = dist < ChannelRadius
            goodidx[j,i] = good

    goodpos = np.zeros((nUnits, nChannels, 3))
    for i in range(nUnits):
        goodpos[i] = ChannelPos[SessionID[i]] * np.tile(goodidx[i], (3,1)).T #gives the 3-d positions of the channels if they are close to the max site  

    # Wave_filt, is the waveform, only at 'good' spatial points and a waveidx/good time points
    WaveformFilt = np.zeros_like(waveform)
    mask = np.zeros_like(waveform)
    for i in range(nUnits):
        mask[i, waveidx[0]:waveidx[-1], goodidx[i,:].astype(bool),:] = 1 

    WaveformFilt = mask * waveform   #applying the filter, so only allow good time points/spatial points

    SpatialFootprintFilt = get_spatialfp(WaveformFilt)

    MaxSite = get_max_site(SpatialFootprintFilt) #This is the MAX site of each individual cv

    return MaxSite, goodidx, goodpos, MaxSiteMean

def exponential_func(d, p_1, p_2):
    """ 
    function, which the waveform is fitted to, to get the spatial decay parameters 
    """
    return p_1 * np.exp(-p_2 * d)

def smooth(array, size = 2):
    """ 
    array is a 1-D np.array
    size is the size of the smoothing window
    function is ~equivalent to https://uk.mathworks.com/help/curvefit/smooth.html
    """
   
    tmp = np.ones_like(array)
    filt = np.ones(size)
    out = np.convolve(np.squeeze(array), filt,'same') / np.convolve(tmp, filt ,'same')      
    return out

def decay_and_average_Waveform(waveform,ChannelPos, goodidx, MaxSite, MaxSiteMean, ClusInfo, param):
    """
    This functions, extracts decay parameters of the units, and uses them to create weighted average waveforms 
    of each unit.
    """
    nUnits = param['nUnits']
    SpikeWidth = param['SpikeWidth']
    nChannels = param['nChannels']
    NewPeakLoc = param['PeakLoc']
    waveidx = param['waveidx']
    ChannelRadius = param['ChannelRadius']
    SessionID = ClusInfo['SessionID']


    SpatialDecayFit = np.zeros((nUnits,2))
    SpatialDecay = np.zeros((nUnits,2))
    d_10 = np.zeros((nUnits, 2))
    AvgCentroid = np.zeros((3,nUnits, 2))
    AvgWaveform = np.zeros((SpikeWidth, nUnits, 2))
    PeakTime = np.zeros((nUnits, 2))



    for i in range(nUnits):
        #use the good indices calculated, to get the nearby positions of each unit
        goodpos = ChannelPos[SessionID[i]][goodidx[i,:].astype(bool),:]
        for cv in range(2):
            Dist2MaxChan = np.linalg.norm( goodpos - ChannelPos[SessionID[i]][MaxSite[i,cv]], axis= 1 )
            TmpAmp = abs(waveform[i,NewPeakLoc,goodidx[i,:].astype(bool),cv])

            # need to remove 0 values, as divide by Dist2MaxChan, and need TmpAmp to be same size
            TmpAmp = TmpAmp[Dist2MaxChan != 0]
            Dist2MaxChan = Dist2MaxChan[Dist2MaxChan != 0]

            # there is variation in how different programming languages/options fit to a curve
            popt, pcurve = sp.optimize.curve_fit(exponential_func, Dist2MaxChan.T , TmpAmp, p0 = (np.max(TmpAmp), 0.05), method = 'trf' )
            #popt, pcurve = sp.optimize.least_squares(exponential_func, Dist2MaxChan, TmpAmp, p0 = (np.max(TmpAmp), 0.05) )

            SpatialDecayFit[i,cv] = popt[1]
            SpatialDecay[i,cv] = np.mean(TmpAmp / Dist2MaxChan)

            tmpmin = np.log(10) / popt[1] # min value which is above noise

            #this if statement shouldn't trigger for a good unit
            if tmpmin > ChannelRadius or tmpmin <0:
                tmpmin = ChannelRadius

            d_10[i,cv] = tmpmin # distance to where the amplitude decay to 10% of its peak

            # Find channel sites which are within d_10 of the max site for that unit and cv
            tmpidx = np.empty(nChannels)
            for s in range(nChannels): #looping over each site

########################################################################################################
                #can change to use max site of each cv, or max site of the mean of each cv
                #dist = np.linalg.norm(ChannelPos[SessionID[i]][MaxSite[i],:] - ChannelPos[SessionID[i]][s,:])
                dist = np.linalg.norm(ChannelPos[SessionID[i]][MaxSiteMean[i],:] - ChannelPos[SessionID[i]][s,:])
########################################################################################################

                good = dist < d_10[i,cv]
                tmpidx[s] = good

            loc = ChannelPos[SessionID[i]][tmpidx.astype(bool),:]

            #average centroid is the sum of spatial footprint * position / spatial foot print
            SpatialFootprint = np.max(np.abs(waveform[i,:,tmpidx.astype(bool),cv]), axis = 1)
            SpatialFootprint = np.expand_dims(SpatialFootprint, axis = -1)
            mu = np.sum( np.tile(SpatialFootprint[:], (1,3)) * loc, axis = 0) / np.sum(SpatialFootprint[:])

            AvgCentroid[:,i,cv] = mu

            #Weighted average waveform
            Dist2MaxProj = np.linalg.norm(loc - AvgCentroid[:,i,cv].T, axis = 1) #it is transposed here
            weight = (d_10[i,cv] - Dist2MaxProj) / d_10[i,cv]
            AvgWaveform[:,i,cv] = np.nansum( waveform[i,:,tmpidx.astype(bool),cv].T * np.tile(weight, (SpikeWidth, 1)), axis = 1 ) / np.sum(weight)
            
            #significant time points 
            wvdurtmp = np.argwhere( np.abs(AvgWaveform[:,i,cv]) - np.mean(AvgWaveform[0:20,i,cv]) > 2.5 * np.std(AvgWaveform[0:20,i,cv], axis = 0))
            if wvdurtmp.size == 0:
                wvdurtmp = waveidx
            wvdurtmp = wvdurtmp[np.isin(wvdurtmp, waveidx)] #okay over achiever, gonna cut you off there
            if wvdurtmp.size == 0:
                wvdurtmp = waveidx  

            # may want to ad smoothing here?
            #PeakTime[i,cv] = np.argmax( np.abs(smooth(WeightedAvgWaveF[wvdurtmp[0]:wvdurtmp[-1],i,cv], 2) ) )
            PeakTime[i,cv] = np.argmax( np.abs(AvgWaveform[wvdurtmp[0]:wvdurtmp[-1] + 1 ,i,cv])) #potential  +1 due to python/ML difference
            PeakTime[i,cv] = PeakTime[i,cv] + wvdurtmp[0] 

    return SpatialDecayFit , SpatialDecay,  d_10, AvgCentroid, AvgWaveform, PeakTime


def get_amplitude_shift_Waveform(waveform,AvgWaveform, PeakTime, param):
    """
    This function aligns the different cv as to maximize their corelation, then shift BOTH cv's by the SAME amount
    so cv 1 has its peak at param['PeakLoc']
    """
    nUnits = param['nUnits']
    SpikeWidth = param['SpikeWidth']
    NewPeakLoc = param['PeakLoc']

    Amplitude = np.zeros((nUnits,2))

    for i in range(nUnits):
        # Shift cv 2 so it has the there is maximum alignment between the 2 waveforms
        tmpcorr = np.correlate(AvgWaveform[:,i,0], AvgWaveform[:,i,1], 'full')
        MaxLag = np.argmax(tmpcorr) - (len(AvgWaveform[:,i,0]) - 1)

        AvgWaveform[:,i,:] = AvgWaveform[:,i,:]
        waveform[i,:,:,:] = waveform[i,:,:,:]

        AvgWaveform[:,i,1] = np.roll(AvgWaveform[:,i,1], MaxLag)
        waveform[i,:,:,1] = np.roll(waveform[i,:,:,1], MaxLag, axis = 0)

        if MaxLag>0:
            AvgWaveform[0:MaxLag ,i,1] = np.nan
            waveform[i,0:MaxLag ,:,1] = np.nan

        elif MaxLag<0:
            AvgWaveform[SpikeWidth - 1 + MaxLag: ,i,1] = np.nan
            waveform[i,SpikeWidth - 1 + MaxLag: ,:,1] = np.nan        


        for cv in [0,1]:
            #shift so both cv, have peak at the prescribed location
            if PeakTime[i,0] != NewPeakLoc: #yes the first CV 
                shift = int(-( PeakTime[i,0] - NewPeakLoc))
                if shift !=0:
                    AvgWaveform[:,i,cv] = np.roll(AvgWaveform[:,i,cv], shift)
                    waveform[i,:,:,cv] = np.roll(waveform[i,:,:,cv], shift, axis = 0)

                    if shift>0:
                        AvgWaveform[0:shift ,i,cv] = np.nan
                        waveform[i,0:shift ,:,cv] = np.nan
                    elif shift<0:
                        # altShift = int(-( PeakTime[i,cv] - NewPeakLoc)) # This altShift should not be used , include to match with (old?) ML code

                        AvgWaveform[SpikeWidth - 1 + shift: ,i,cv] = np.nan
                        waveform[i,SpikeWidth - 1 + shift: ,:,cv] = np.nan    

            tmpPeak = AvgWaveform[NewPeakLoc,i,cv]

            if np.isnan(tmpPeak) == True : # This is a catch for  bad units
                #tmpPeak = np.nanmax(WeightedAvgWaveF_copy[:,i,cv]) 
                tmpPeak = np.nan
                print(f'This unit {i}, CV {cv} is very likely a bad unit!')
            Amplitude[i,cv] = tmpPeak  

    return Amplitude, waveform, AvgWaveform


def avg_Waveform_PerTP(waveform,ChannelPos, d_10, MaxSiteMean, Amplitude, AvgWaveform, ClusInfo, param):
    """
    This function calculates the weighted average waveform per time point, as well as the good time points for each unit
    """
    nUnits = param['nUnits']
    SpikeWidth = param['SpikeWidth']
    nChannels = param['nChannels']
    waveidx = param['waveidx']
    SessionID = ClusInfo['SessionID']

    GoodSiteId = np.empty((nUnits,nChannels,2))
    WaveformDuration = np.full((nUnits,2), np.nan)
    AvgWaveformPerTP = np.full((3, nUnits,SpikeWidth,2), np.nan)
    WaveIdx = np.zeros((nUnits, SpikeWidth, 2))

    for i in range(nUnits):
        
        for cv in range(2):
            #dist = np.linalg.norm(ChannelPos[SessionID[i]][MaxSite[i,cv],:] - ChannelPos[SessionID[i]][:,:], axis = 1)
            dist = np.linalg.norm(ChannelPos[SessionID[i]][MaxSiteMean[i],:] - ChannelPos[SessionID[i]][:,:], axis = 1)

            test = dist < np.abs(d_10[i,cv])
            GoodSiteId[i,:,cv] = test
            Locs = ChannelPos[SessionID[i]][GoodSiteId[i,:,cv].astype(bool), :]

            # select time points where the values are above 25% of the amplitude
            wvdurtmp = np.argwhere( np.abs(np.sign(Amplitude[i,cv]) * AvgWaveform[waveidx,i,cv] )> (np.sign(Amplitude[i,cv]) * Amplitude[i,cv] * 0.25))

            if wvdurtmp.sum != 0:
                try:
                    wvdurtmp = np.arange(wvdurtmp[0],wvdurtmp[-1]+1) + waveidx[0]  
                    wvdurtmp = np.expand_dims(wvdurtmp, axis = 1)
                    WaveformDuration[i,cv] = len(wvdurtmp) 
                except:
                    print(f'unit{i} is very likely bad, no good time points in average wavefunction')  
                    WaveformDuration[i,cv] = np.nan    

        
            else:
                print(f'problem with wvdurtmp for unit {i}')
            
        #Projected location per time point 
            for idx in wvdurtmp:
                tmp = np.abs(waveform[i,idx,GoodSiteId[i,:,cv].astype(bool),cv])
                tmp = np.expand_dims(tmp, axis = 1)
                tmp = np.tile(tmp, (1,3))
                tmp = np.sum(tmp * Locs , axis = 0) / np.sum(tmp, axis = 0)
                AvgWaveformPerTP[:,i,idx,cv] = tmp.reshape((3,1))
                WaveIdx[i,idx,cv] = 1 
    
    return WaveformDuration, AvgWaveformPerTP, WaveIdx
        