# This file will containall the necessary function for extracting waveform parameters from average waveforms

import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import detrend
import scipy as sp


import pandas as pd


def detrend_waveform(waveform):
    """
    This function accepts the raw waveforms (no_units, SpikeWidth, no_channels, CV)
    The output is the same shape, and linearly detrended across time.
    """ 
    return detrend(waveform, axis = 2)

def get_spatialfp(waveform):
    """
    Input: waveform np array (No. units, time, channel No., first/second half)
    By taking the maximum value along the time axis of the absolute values
    Output: np array (No. units, channel No., first/second half)
    """
    spatial_fp = np.max(np.abs(waveform), axis = 1)
    return spatial_fp

def get_max_site(spatial_fp):
    """
    Input: spatial_fp (No. units, channel No., first/second half)
    By taking the index of maximum argument along the channel axis,
    Output: (No units, first/second half), which gives the maximum site for each unit in first/second half
    """
    MaxSite = np.argmax(spatial_fp, axis = 1)
    return MaxSite

def get_max_sites(waveform, ChannelPos, param):
    """
    Using waveforms, ChannelPos and param, to find the max channel for each unit and cv, this function also
    returns good idx's / positions, by selecting channels within ChannelRadius (default 150 um) 
    """

    n_units = param['n_units']
    n_channels = param['n_channels']
    ChannelRadius = param['ChannelRadius']
    waveidx = param['waveidx']
    
    mean_cv = np.mean(waveform, axis = 3) # average of each cv
    spatial_footprint = get_spatialfp(mean_cv) # choose max time 
    MaxSiteMean= get_max_site(spatial_footprint) # argument of MaxSite

    # Finds the indices where the distance between the avg max site is small
    good_idx = np.empty((n_units, n_channels))
    for i in range(ChannelPos.shape[0]): #looping over each site
        dist = np.linalg.norm(ChannelPos[MaxSiteMean,:] - ChannelPos[i,:], axis=1)
        good = dist < ChannelRadius
        good_idx[:,i] = good

    good_pos = np.zeros((n_units, n_channels, 3))
    for i in range(n_units):
        good_pos[i] = ChannelPos * np.tile(good_idx[i], (3,1)).T #gives the 3-d positions of the channels if they are close to the max site  

    # Wave_filt, is the waveform, only at 'good' spatial points
    wavef_filt = np.zeros_like(waveform)
    mask = np.zeros_like(waveform)
    for i in range(n_units):
        mask[i, waveidx[0]:waveidx[-1], good_idx[i,:].astype(bool),:] = 1 

    wavef_filt = mask * waveform

    filt_spfp = get_spatialfp(wavef_filt)

    filt_max_site = get_max_site(filt_spfp) #This is the MAX site of individual cv
    #waveforms, using a filter waveform to be robust to noise
    return filt_max_site, good_idx, good_pos, MaxSiteMean

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

def decay_and_average_WaveF(waveform,ChannelPos, good_idx, max_site, MaxSiteMean, param):
    """
    This functions, extracts decay parameters of the units, and uses them to create weighted average waveforms 
    of each unit.
    """
    n_units = param['n_units']
    SpikeWidth = param['SpikeWidth']
    n_channels = param['n_channels']
    NewPeakLoc = param['PeakLoc']
    waveidx = param['waveidx']
    ChannelRadius = param['ChannelRadius']


    SpatialDecayFit = np.zeros((n_units,2))
    SpatialDecay = np.zeros((n_units,2))
    d_10 = np.zeros((n_units, 2))
    AvgCentroid = np.zeros((3,n_units, 2))
    WeightedAvgWaveF = np.zeros((SpikeWidth, n_units, 2))
    PeakTime = np.zeros((n_units, 2))



    for i in range(n_units):
        #use the good indices calculated, to get the nearby positions of each unit
        near_pos = ChannelPos[good_idx[i,:].astype(bool),:]
        for cv in range(2):
            Dist2MaxChan = np.linalg.norm( near_pos - ChannelPos[max_site[i,cv]], axis= 1 )
            tmp_amp = abs(waveform[i,NewPeakLoc,good_idx[i,:].astype(bool),cv])

            # need to remove 0 values, as divide by Dist2MaxChan, and need tmp_amp to be same size
            tmp_amp = tmp_amp[Dist2MaxChan != 0]
            Dist2MaxChan = Dist2MaxChan[Dist2MaxChan != 0]

            # there is variation in how different languages/options fit to a curve
            popt, pcurve = sp.optimize.curve_fit(exponential_func, Dist2MaxChan.T , tmp_amp, p0 = (np.max(tmp_amp), 0.05), method = 'trf' )
            #popt, pcurve = sp.optimize.least_squares(exponential_func, Dist2MaxChan, tmp_amp, p0 = (np.max(tmp_amp), 0.05) )

            SpatialDecayFit[i,cv] = popt[1]
            SpatialDecay[i,cv] = np.mean(tmp_amp / Dist2MaxChan)

            tmpmin = np.log(10) / popt[1] # min value which is above noise

            #this if statement shouldn't trigger for a good unit
            if tmpmin > ChannelRadius or tmpmin <0:
                tmpmin = ChannelRadius
                # print(f'this index{i} and cv {cv} produces a highly questionable d_10 distance')

            d_10[i,cv] = tmpmin # distance to where the amplitude decay to 10% of its peak

            # Now need to get indices and location for specific cv
            tmp_idx = np.empty((n_channels))
            for s in range(n_channels): #looping over each site

########################################################################################################
                #can test distance to its max site, or the mean of it two halves max site, default is mean
                #dist = np.linalg.norm(ChannelPos[max_site[i],:] - ChannelPos[s,:])
                dist = np.linalg.norm(ChannelPos[MaxSiteMean[i],:] - ChannelPos[s,:])
########################################################################################################

                good = dist < d_10[i,cv]
                tmp_idx[s] = good

            loc = ChannelPos[tmp_idx.astype(bool),:]

            #average centroid, sum of spatial footprint * position / spatial foot print
            spatial_fp = np.max(np.abs(waveform[i,:,tmp_idx.astype(bool),cv]), axis = 1)
            spatial_fp = np.expand_dims(spatial_fp, axis = -1)
            mu = np.sum( np.tile(spatial_fp[:], (1,3)) * loc, axis = 0) / np.sum(spatial_fp[:])

            AvgCentroid[:,i,cv] = mu

            #Weighted average waveform
            Dist2MaxProj = np.linalg.norm(loc - AvgCentroid[:,i,cv].T, axis = 1) #it is transposed here
            weight = (d_10[i,cv] - Dist2MaxProj) / d_10[i,cv]
            WeightedAvgWaveF[:,i,cv] = np.nansum( waveform[i,:,tmp_idx.astype(bool),cv].T * np.tile(weight, (SpikeWidth, 1)), axis = 1 ) / np.sum(weight)
            
            #significant time points 
            wvdurtmp = np.argwhere( np.abs(WeightedAvgWaveF[:,i,cv]) - np.mean(WeightedAvgWaveF[0:20,i,cv]) > 2.5 * np.std(WeightedAvgWaveF[0:20,i,cv], axis = 0))
            if wvdurtmp.size == 0:
                wvdurtmp = waveidx
            wvdurtmp = wvdurtmp[np.isin(wvdurtmp, waveidx)] #okay over achiever, gonna cut you off there
            if wvdurtmp.size == 0:
                wvdurtmp = waveidx  

            # may want to ad smoothing here?
            #PeakTime[i,cv] = np.argmax( np.abs(smooth(WeightedAvgWaveF[wvdurtmp[0]:wvdurtmp[-1],i,cv], 2) ) )
            PeakTime[i,cv] = np.argmax( np.abs(WeightedAvgWaveF[wvdurtmp[0]:wvdurtmp[-1] + 1 ,i,cv])) # +1 due to python/ML difference
            PeakTime[i,cv] = PeakTime[i,cv] + wvdurtmp[0] 

    return SpatialDecayFit , SpatialDecay,  d_10, AvgCentroid, WeightedAvgWaveF, PeakTime


def get_amplitude_shift_WaveF(waveform,WeightedAvgWaveF, PeakTime, param):
    """
    This function aligns the different cv to maximize their corelation, the shifts BOTH cv by the SAME amount
    so cv 1 has its peak at param['PeakLoc']
    """
    n_units = param['n_units']
    SpikeWidth = param['SpikeWidth']
    NewPeakLoc = param['PeakLoc']

    Amplitude = np.zeros((n_units,2))

    for i in range(n_units):
        # Shift cv 2 so it has the there is maximum alignment between the 2 waveforms
        tmpcorr = np.correlate(WeightedAvgWaveF[:,i,0], WeightedAvgWaveF[:,i,1], 'full')
        max_lag = np.argmax(tmpcorr) - (len(WeightedAvgWaveF[:,i,0]) - 1)

        WeightedAvgWaveF[:,i,:] = WeightedAvgWaveF[:,i,:]
        waveform[i,:,:,:] = waveform[i,:,:,:]

        WeightedAvgWaveF[:,i,1] = np.roll(WeightedAvgWaveF[:,i,1], max_lag)
        waveform[i,:,:,1] = np.roll(waveform[i,:,:,1], max_lag, axis = 0)

        if max_lag>0:
            WeightedAvgWaveF[0:max_lag ,i,1] = np.nan
            waveform[i,0:max_lag ,:,1] = np.nan

        elif max_lag<0:
            WeightedAvgWaveF[SpikeWidth - 1 + max_lag: ,i,1] = np.nan
            waveform[i,SpikeWidth - 1 + max_lag: ,:,1] = np.nan        


        for cv in [0,1]:
            #shift so both cv, have peak at the prescribed location
            if PeakTime[i,0] != NewPeakLoc: #yes the first CV 
                shift = int(-( PeakTime[i,0] - NewPeakLoc))
                if shift !=0:
                    WeightedAvgWaveF[:,i,cv] = np.roll(WeightedAvgWaveF[:,i,cv], shift)
                    waveform[i,:,:,cv] = np.roll(waveform[i,:,:,cv], shift, axis = 0)

                    if shift>0:
                        WeightedAvgWaveF[0:shift ,i,cv] = np.nan
                        waveform[i,0:shift ,:,cv] = np.nan
                    elif shift<0:
                        # alt_shift = int(-( PeakTime[i,cv] - NewPeakLoc)) # This alt_shift should not be used , include to match with (old?) ML code

                        WeightedAvgWaveF[SpikeWidth - 1 + shift: ,i,cv] = np.nan
                        waveform[i,SpikeWidth - 1 + shift: ,:,cv] = np.nan    

            temp_peak = WeightedAvgWaveF[NewPeakLoc,i,cv]

            if np.isnan(temp_peak) == True : # This is a catch for  bad units
                #temp_peak = np.nanmax(WeightedAvgWaveF_copy[:,i,cv])
                temp_peak = np.nan
                print(f'This unit {i}, CV {cv} is very likely a bad unit!')
            Amplitude[i,cv] = temp_peak  


    return Amplitude, waveform, WeightedAvgWaveF
                                                                          
def avg_WaveF_PerTP(waveform,ChannelPos, d_10, MaxSiteMean, Amplitude, WeightedAvgWaveF, param):
    """
    This function calculates the weighted average wavefunction per time point, as well as the good time points for each unit
    """
    n_units = param['n_units']
    SpikeWidth = param['SpikeWidth']
    n_channels = param['n_channels']
    waveidx = param['waveidx']

    good_site_id = np.empty((n_units,n_channels,2))
    WaveformDuration = np.full((n_units,2), np.nan)
    WeightedAvgWaveF_PerTP = np.full((3, n_units,SpikeWidth,2), np.nan)
    WaveIdx = np.zeros((n_units, SpikeWidth, 2))

    for i in range(n_units):
        
        for cv in range(2):
            #dist = np.linalg.norm(ChannelPos[filt_max_site[i,cv],:] - ChannelPos[:,:], axis = 1)
            dist = np.linalg.norm(ChannelPos[MaxSiteMean[i],:] - ChannelPos[:,:], axis = 1)

            test = dist < np.abs(d_10[i,cv])
            good_site_id[i,:,cv] = test
            Locs = ChannelPos[good_site_id[i,:,cv].astype(bool), :]

            # select time points where the values are above 25% of the amplitude
            wvdurtmp = np.argwhere( np.abs(np.sign(Amplitude[i,cv]) * WeightedAvgWaveF[waveidx,i,cv] )> (np.sign(Amplitude[i,cv]) * Amplitude[i,cv] * 0.25))

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
                tmp_1 = np.abs(waveform[i,idx,good_site_id[i,:,cv].astype(bool),cv])
                tmp_1 = np.expand_dims(tmp_1, axis = 1)
                tmp_1 = np.tile(tmp_1, (1,3))
                tmp = np.sum(tmp_1 * Locs , axis = 0) / np.sum(tmp_1, axis = 0)
                WeightedAvgWaveF_PerTP[:,i,idx,cv] = tmp.reshape((3,1))
                WaveIdx[i,idx,cv] = 1 
    
    return WaveformDuration, WeightedAvgWaveF_PerTP, WaveIdx
        