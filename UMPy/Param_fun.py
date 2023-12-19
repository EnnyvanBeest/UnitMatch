# This file will containall the necessary function for extracting waveform parameters from average waveforms

import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import detrend
import scipy as sp


import pandas as pd


def detrendwaveform(waveform):
    """
    This function accepts the raw waveforms (no_units, spike_width, no_channels, CV)
    The output is the same shape, and linearly detrended acros time.
    """ 
    return detrend(waveform, axis = 2)

def get_spatialfp(waveform):
    """
    Input: waveform nparray (No. units, time, channel No., first/second half)
    By taking the maximum value along the time axis of the absolute values
    Output: nparray (No. units, channel No., first/second half)
    """
    spatial_fp = np.max(np.abs(waveform), axis = 1)
    return spatial_fp

def get_max_site(spatial_fp):
    """
    Input: spatial_fp (No. units, channel No., first/second half)
    By taking the index of maximum argument allong the channel axis,
    Output: (No units, first/second half), whcih gives the maximum site for each unit in first/second half
    """
    max_site = np.argmax(spatial_fp, axis = 1)
    return max_site

def get_max_sites(waveform, channel_pos, param):
    """
    Using waveforms, channel_pos and parm, to find the max channel for each unit and cv, this function also
    returns good idx's / positions, by selecting channels within channel_radius (default 150 um) 
    """

    n_units, n_channels = waveform.shape[0], waveform.shape[2]
    channel_radius = param['channel_radius']
    
    mean_cv = np.mean(waveform, axis = 3) # average of each half
    spatial_footprint = get_spatialfp(mean_cv) # choose max time 
    max_site_mean= get_max_site(spatial_footprint) # argument of max_site

    # Finds the indicies where the distance between tha avg max site is small
    good_idx = np.empty((n_units, n_channels))
    for i in range(channel_pos.shape[0]): #looping over each site
        dist = np.linalg.norm(channel_pos[max_site_mean,:] - channel_pos[i,:], axis=1)
        good = dist < channel_radius
        good_idx[:,i] = good

    # good index shape(n_units, n_channels), is true if distance is small enough

    good_pos = np.zeros((n_units, n_channels, 3))
    for i in range(n_units):
        good_pos[i] = channel_pos * np.tile(good_idx[i], (3,1)).T
    #gives the 3-d positions of the channels if they are close to the max site
    

    # Wave_filt, is the waveform, only at 'good' spatial points
    wavef_filt = np.zeros_like(waveform)
    mask = np.zeros_like(waveform)
    for i in range(n_units):
        mask[i,34:71,good_idx[i,:].astype(bool),:] = 1 

    wavef_filt = mask * waveform

    filt_spfp = get_spatialfp(wavef_filt)

    filt_max_site = get_max_site(filt_spfp) #This is the MAX site of individual cv
    #waveforms, using a filter waveform to be rhobust to noise
    return filt_max_site, good_idx, good_pos, max_site_mean

def expontetial_func(d, p_1, p_2):
    """ 
    function, which the waveform is fitted to, to get the spatial decay parameters 
    """
    return p_1 * np.exp(-p_2 * d)

def smooth(array, size = 2):
    """ 
    array is a 1-D np.array
    size is the size of the smoothing window
    function is ~equvilent to https://uk.mathworks.com/help/curvefit/smooth.html
    """
   
    tmp = np.ones_like(array)
    filt = np.ones(size)
    out = np.convolve(np.squeeze(array), filt,'same') / np.convolve(tmp, filt ,'same')      
    return out

def decay_and_average_WaveF(waveform,channel_pos, good_idx, max_site, max_site_mean, param):
    """
    This functions, extracts decay parameters of the units, and uses them to create weighted average waveforms 
    of each unit.
    """
    n_units, spike_width, n_channels,__ = waveform.shape
    new_peak_loc = param['peak_loc']
    waveidx = param['waveidx']
    channel_radius = param['channel_radius']


    spatialdecayfit = np.zeros((n_units,2))
    spatialdecay = np.zeros((n_units,2))
    d_10 = np.zeros((n_units, 2))
    avg_centroid = np.zeros((3,n_units, 2))
    WeightedAvgWaveF = np.zeros((spike_width, n_units, 2))
    PeakTime = np.zeros((n_units, 2))



    for i in range(n_units):
        #use the good indices calculated, to get the nearby postions of each unit
        near_pos = channel_pos[good_idx[i,:].astype(bool),:]
        for cv in range(2):
            Dist2MaxChan = np.linalg.norm( near_pos - channel_pos[max_site[i,cv]], axis= 1 )
            tmp_amp = abs(waveform[i,new_peak_loc,good_idx[i,:].astype(bool),cv])

            # need to remove 0 values, as divide by Dist2MaxChan, and need tmp_amp to be same size
            tmp_amp = tmp_amp[Dist2MaxChan != 0]
            Dist2MaxChan = Dist2MaxChan[Dist2MaxChan != 0]

            # there is variation in how different languages/options fit to a curve
            popt, pcurve = sp.optimize.curve_fit(expontetial_func, Dist2MaxChan.T , tmp_amp, p0 = (np.max(tmp_amp), 0.05), method = 'trf' )
            #popt, pcurve = sp.optimize.least_squares(expontetial_func, Dist2MaxChan, tmp_amp, p0 = (np.max(tmp_amp), 0.05) )

            spatialdecayfit[i,cv] = popt[1]
            spatialdecay[i,cv] = np.mean(tmp_amp / Dist2MaxChan)

            tmpmin = np.log(10) / popt[1]

            #this if statement shouldn't trigger for a good unit
            if tmpmin > channel_radius or tmpmin <0:
                tmpmin = channel_radius
                # print(f'this index{i} and cv {cv} produces a highly questionable d_10 distance')

            d_10[i,cv] = tmpmin # distance to where the amplitude decay to 10% of its peak

            # Now need to get indices and location for specfic cv
            tmp_idx = np.empty((n_channels))
            for s in range(n_channels): #looping over each site

########################################################################################################
                #can test distance to its max site, or the mean of it two halves max site, default is mean
                #dist = np.linalg.norm(channel_pos[max_site[i],:] - channel_pos[s,:])
                dist = np.linalg.norm(channel_pos[max_site_mean[i],:] - channel_pos[s,:])
########################################################################################################

                good = dist < d_10[i,cv]
                tmp_idx[s] = good

            loc = channel_pos[tmp_idx.astype(bool),:]

            #average centroid, sum of spatial footprint * postion / spatial foot print
            spatial_fp = np.max(np.abs(waveform[i,:,tmp_idx.astype(bool),cv]), axis = 1)
            spatial_fp = np.expand_dims(spatial_fp, axis = -1)
            mu = np.sum( np.tile(spatial_fp[:], (1,3)) * loc, axis = 0) / np.sum(spatial_fp[:])

            avg_centroid[:,i,cv] = mu

            #Weighted average waveform
            Dist2MaxProj = np.linalg.norm(loc - avg_centroid[:,i,cv].T, axis = 1) # is it transpose here
            weight = (d_10[i,cv] - Dist2MaxProj) / d_10[i,cv]
            WeightedAvgWaveF[:,i,cv] = np.nansum( waveform[i,:,tmp_idx.astype(bool),cv].T * np.tile(weight, (spike_width, 1)), axis = 1 ) / np.sum(weight)
            
            #significant timepoints 
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

    return spatialdecayfit , spatialdecay,  d_10, avg_centroid, WeightedAvgWaveF, PeakTime


def get_amplitude_shifWaveF(waveform,WeightedAvgWaveF, PeakTime, param):
    """
    This function alligns the different cv to maximise thier corelation, the shifts BOTH cv bythe SAME amount
    so cv 1 has its peak at param['peak_loc']
    """
    n_units, spike_width, __ , __= waveform.shape
    new_peak_loc = param['peak_loc']


    # stop unsing a copy of the waveform
    Amplitude = np.zeros((n_units,2))

    for i in range(n_units):
        # Shift cv 2 so it has the there is maximum allignement betwen the 2 waveforms
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
            WeightedAvgWaveF[spike_width - 1 + max_lag: ,i,1] = np.nan
            waveform[i,spike_width - 1 + max_lag: ,:,1] = np.nan        


        for cv in [0,1]:
            #shift so both cv, have peak at the prescribed location
            if PeakTime[i,0] != new_peak_loc: #yes the first CV 
                shift = int(-( PeakTime[i,0] - new_peak_loc))
                if shift !=0:
                    WeightedAvgWaveF[:,i,cv] = np.roll(WeightedAvgWaveF[:,i,cv], shift)
                    waveform[i,:,:,cv] = np.roll(waveform[i,:,:,cv], shift, axis = 0)

                    if shift>0:
                        WeightedAvgWaveF[0:shift ,i,cv] = np.nan
                        waveform[i,0:shift ,:,cv] = np.nan
                    elif shift<0:
                        # alt_shift = int(-( PeakTime[i,cv] - new_peak_loc)) # This alt_shift shuld not be used , include to match with (old?) ML code

                        WeightedAvgWaveF[spike_width - 1 + shift: ,i,cv] = np.nan
                        waveform[i,spike_width - 1 + shift: ,:,cv] = np.nan    

            temp_peak = WeightedAvgWaveF[new_peak_loc,i,cv]

            if np.isnan(temp_peak) == True : # This is a catch for  bad units
                #temp_peak = np.nanmax(WeightedAvgWaveF_copy[:,i,cv])
                temp_peak = np.nan
                print(f'This unit {i}, CV {cv} is very likely a bad unit!')
            Amplitude[i,cv] = temp_peak  


    return Amplitude, waveform, WeightedAvgWaveF
                                                                          
def avg_WaveF_PerTP(waveform,channel_pos, d_10, max_site_mean, Amplitude, WeightedAvgWaveF, param):
    """
    This function calculates the weighted average wavefunction per time point, as well as the good time points for each unit
    """
    n_units, spike_width, n_channels,__ = waveform.shape
    waveidx = param['waveidx']

    good_site_id = np.empty((n_units,n_channels,2))
    WaveformDuration = np.full((n_units,2), np.nan)
    WeightedAvgWaveF_PerTP = np.full((3, n_units,spike_width,2), np.nan)
    WaveIdx = np.zeros((n_units, spike_width, 2))

    for i in range(n_units):
        
        for cv in range(2):
            #dist = np.linalg.norm(channel_pos[filt_max_site[i,cv],:] - channel_pos[:,:], axis = 1)
            dist = np.linalg.norm(channel_pos[max_site_mean[i],:] - channel_pos[:,:], axis = 1)

            test = dist < np.abs(d_10[i,cv])
            good_site_id[i,:,cv] = test
            Locs = channel_pos[good_site_id[i,:,cv].astype(bool), :]


            wvdurtmp = np.argwhere( np.abs(np.sign(Amplitude[i,cv]) * WeightedAvgWaveF[waveidx,i,cv] )> (np.sign(Amplitude[i,cv]) * Amplitude[i,cv] * 0.25))

            if wvdurtmp.sum != 0:

                #wvdurtmp = wvdurtmp + waveidx[0] 
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
        