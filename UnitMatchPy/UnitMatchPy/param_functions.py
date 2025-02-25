# This file will contain all the necessary function for extracting waveform parameters from average waveforms
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import detrend
import scipy as sp
from scipy.ndimage import gaussian_filter
import pandas as pd

def detrend_waveform(waveform):
    """
    This function linearly de-trend the raw waveform over time.

    Parameters
    ----------
    waveform : ndarray (n_units, spike_width, n_channels, cv)
        The average waveform for each unit and each cv

    Returns
    -------
    ndarray 
        The linearly de-trended waveform  
    """
    return detrend(waveform, axis = 1)

def get_spatial_fp(waveform):
    """
    This function finds the maximum values over time to find the spatial foot print

    Parameters
    ----------
    waveform : ndarray (n_units, spike_width, n_channels, cv)
        The average waveform for each unit and each cv

    Returns
    -------
    ndarray 
        The spatial footprint (n_units, n_channels, cv)
    """
    spatial_fp = np.max(np.abs(waveform), axis = 1)
    return spatial_fp

def get_max_site(spatial_fp):
    """
    This function finds the maximum channels for each unit.

    Parameters
    ----------
    spatial_fp : ndarray (n_units, n_channels, cv)
        The spatial footprint for each unit and each cv

    Returns
    -------
    ndarray 
        The maximum channels for each unit ad cv (n_units, cv)
    """
    max_site = np.argmax(spatial_fp, axis = 1)
    return max_site

def get_max_sites(waveform, channel_pos, clus_info, param):
    """
    This functions finds the max sites for each unit as well as the 'good' indexes and positions to use.

    Parameters
    ----------
    waveform : ndarray (n_units, spike_width, n_channels, cv)
        The waveforms for each unit and cv
    channel_pos : (n_units, 3)
        The spatial position of each unit
    clus_info : dict
        The clus_info dictionary
    param : dict
        The param dictionary

    Returns
    -------
    ndarray
        The max site of each unit, the good idx/positions for each unit
    """

    n_units = param['n_units']
    n_channels = param['n_channels']
    channel_radius = param['channel_radius']
    waveidx = param['waveidx']
    session_id = clus_info['session_id']
    
    mean_cv = np.mean(waveform, axis = 3) # average of each cv
    spatial_fp = get_spatial_fp(mean_cv) # choose max time
    max_site_mean= get_max_site(spatial_fp) # argument of MaxSite

    # Finds the indices where the distance from the max site mean is small
    good_idx = np.empty((n_units, n_channels))
    for i in range(channel_pos[0].shape[0]): #looping over each site, assuming all channel_pos, are the same no of channels
        for j in range(n_units):
            dist = np.linalg.norm(channel_pos[session_id[j]][max_site_mean[j],:] - channel_pos[session_id[j]][i,:])
            good = dist < channel_radius
            good_idx[j,i] = good

    good_pos = np.zeros((n_units, n_channels, 3))
    for i in range(n_units):
        #gives the 3-d positions of the channels if they are close to the max site
        good_pos[i] = channel_pos[session_id[i]] * np.tile(good_idx[i], (3,1)).T

    # waveform_filt, is the waveform, only at 'good' spatial points and a waveidx/good time points
    waveform_filt = np.zeros_like(waveform)
    mask = np.zeros_like(waveform)
    for i in range(n_units):
        mask[i, waveidx[0]:waveidx[-1], good_idx[i,:].astype(bool),:] = 1 

    waveform_filt = mask * waveform #applying the filter, so only allow good time points/spatial points

    spatial_fp_filt = get_spatial_fp(waveform_filt)

    max_site = get_max_site(spatial_fp_filt) #This is the max site of each individual cv

    return max_site, good_idx, good_pos, max_site_mean

def exponential_func(d, p_1, p_2):
    """
    The exponential decay function the waveform is fitted to, to get the decay parameters
    """
    return p_1 * np.exp(-p_2 * d)

def smooth(array, size = 2):
    """
    This function smooths a 1-d array, it is approximately the same as
    https://uk.mathworks.com/help/curvefit/smooth.html

    Parameters
    ----------
    array : ndarray
        The 1D array to be smoothed
    size : int, optional
        The size of the smoothing window, by default 2

    Returns
    -------
    ndarray
        The smoothed array
    """
    tmp = np.ones_like(array)
    filt = np.ones(size)
    out = np.convolve(np.squeeze(array), filt,'same') / np.convolve(tmp, filt ,'same')      
    return out

def decay_and_average_waveform(waveform, channel_pos, good_idx, max_site, max_site_mean, clus_info, param):
    """
    This functions, extracts decay parameters of the units, and uses them to create weighted average waveforms 
    for each unit.

    Parameters
    ----------
    waveform : ndarray (n_units, spike_width. n_channels, cv)
        The waveforms for each unit and cv
    channel_pos : (n_units, 3)
        The spatial position of each unit
    good_idx : ndarray
        The good time points for each unit
    max_site : ndarray
        The maximum site for each unit
    max_site_mean : ndarray
        The maximum site for each unit, when taking the mean over each cv
     clus_info : dict
        The clus_info dictionary
    param : dict
        The param dictionary

    Returns
    -------
    ndarrays
        decay parameters, avg waveform and avg centroid for each unit
    """
    n_units = param['n_units']
    spike_width = param['spike_width']
    n_channels = param['n_channels']
    new_peak_loc = param['peak_loc']
    waveidx = param['waveidx']
    channel_radius = param['channel_radius']
    session_id = clus_info['session_id']

    spatial_decay_fit = np.zeros((n_units,2))
    spatial_decay = np.zeros((n_units,2))
    d_10 = np.zeros((n_units, 2))
    avg_centroid = np.zeros((3,n_units, 2))
    avg_waveform = np.zeros((spike_width, n_units, 2))
    peak_time = np.zeros((n_units, 2))

    for i in range(n_units):
        #use the good indices calculated, to get the nearby positions of each unit
        good_pos = channel_pos[session_id[i]][good_idx[i,:].astype(bool),:]
        for cv in range(2):
            dist_to_max_chan = np.linalg.norm( good_pos - channel_pos[session_id[i]][max_site[i,cv]], axis= 1 )
            tmp_amp = abs(waveform[i, new_peak_loc, good_idx[i,:].astype(bool), cv])

            # need to remove 0 values, as divide by Dist2MaxChan, and need TmpAmp to be same size
            tmp_amp = tmp_amp[dist_to_max_chan != 0]
            dist_to_max_chan = dist_to_max_chan[dist_to_max_chan != 0]

            # there is variation in how different programming languages/options fit to a curve
            popt, pcurve = sp.optimize.curve_fit(exponential_func, dist_to_max_chan.T , tmp_amp, p0 = (np.max(tmp_amp), 0.05), method = 'trf', maxfev=2000)
            #popt, pcurve = sp.optimize.least_squares(exponential_func, Dist2MaxChan, TmpAmp, p0 = (np.max(TmpAmp), 0.05) )

            spatial_decay_fit[i,cv] = popt[1]
            spatial_decay[i,cv] = np.mean(tmp_amp / dist_to_max_chan)

            tmp_min = np.log(10) / popt[1] # min value which is above noise

            #this if statement shouldn't trigger for a good unit
            if tmp_min > channel_radius or tmp_min <0:
                tmp_min = channel_radius

            d_10[i,cv] = tmp_min # distance to where the amplitude decay to 10% of its peak

            # Find channel sites which are within d_10 of the max site for that unit and cv
            tmp_idx = np.empty(n_channels)
            for s in range(n_channels): #looping over each site

########################################################################################################
                #can change to use max site of each cv, or max site of the mean of each cv
                #dist = np.linalg.norm(ChannelPos[SessionID[i]][MaxSite[i],:] - ChannelPos[SessionID[i]][s,:])
                dist = np.linalg.norm(channel_pos[session_id[i]][max_site_mean[i],:] - channel_pos[session_id[i]][s,:])
########################################################################################################

                good = dist < d_10[i,cv]
                tmp_idx[s] = good

            loc = channel_pos[session_id[i]][tmp_idx.astype(bool),:]

            #average centroid is the sum of spatial footprint * position / spatial foot print
            spatial_fp = np.max(np.abs(waveform[i,:,tmp_idx.astype(bool),cv]), axis = 1)
            spatial_fp = np.expand_dims(spatial_fp, axis = -1)
            mu = np.sum( np.tile(spatial_fp[:], (1,3)) * loc, axis = 0) / np.sum(spatial_fp[:])

            avg_centroid[:,i,cv] = mu

            #Weighted average waveform
            dist_to_max_proj = np.linalg.norm(loc - avg_centroid[:,i,cv].T, axis = 1) #it is transposed here
            weight = (d_10[i,cv] - dist_to_max_proj) / d_10[i,cv]
            avg_waveform[:,i,cv] = np.nansum( waveform[i,:,tmp_idx.astype(bool),cv].T * np.tile(weight, (spike_width, 1)), axis = 1 ) / np.sum(weight)
            
            #significant time points 
            wave_duration_tmp = np.argwhere( np.abs(avg_waveform[:,i,cv]) - np.mean(avg_waveform[0:20,i,cv]) > 2.5 * np.std(avg_waveform[0:20,i,cv], axis = 0))
            if wave_duration_tmp.size == 0:
                wave_duration_tmp = waveidx
            wave_duration_tmp = wave_duration_tmp[np.isin(wave_duration_tmp, waveidx)] #okay over achiever, gonna cut you off there
            if wave_duration_tmp.size == 0:
                wave_duration_tmp = waveidx  

            # may want to add smoothing here?
            #PeakTime[i,cv] = np.argmax( np.abs(smooth(WeightedAvgWaveF[wvdurtmp[0]:wvdurtmp[-1],i,cv], 2) ) )

            peak_time[i,cv] = np.argmax( np.abs(avg_waveform[wave_duration_tmp[0]:wave_duration_tmp[-1] + 1 ,i,cv]))
            peak_time[i,cv] = peak_time[i,cv] + wave_duration_tmp[0] 

    return spatial_decay_fit , spatial_decay,  d_10, avg_centroid, avg_waveform, peak_time


def get_amplitude_shift_waveform(waveform, avg_waveform, peak_time, param):
    """
    This function aligns the different cv as to maximize their corelation, then shift BOTH cv's by the SAME amount
    so cv 1 has its peak at param['peak_loc']

    Parameters
    ----------
    waveform : ndarray (n_units, spike_width, n_channels, cv)
        The waveforms for each unit and cv
    avg_waveform : ndarray (n_units, spike_width, cv)
        The weighted average waveform over channels
    peak_time : ndarray
        The peak time for each waveform
    param : dict
        The param dictionary

    Returns
    -------
    ndarrays
        The amplitude for each unit and waveforms re-aligned to each other
    """
    n_units = param['n_units']
    spike_width = param['spike_width']
    new_peak_loc = param['peak_loc']

    amplitude = np.zeros((n_units,2))

    for i in range(n_units):
        # Shift cv 2 so it has the there is maximum alignment between the 2 waveforms
        tmp_corr = np.correlate(avg_waveform[:,i,0], avg_waveform[:,i,1], 'full')
        max_lag = np.argmax(tmp_corr) - (len(avg_waveform[:,i,0]) - 1)

        avg_waveform[:,i,:] = avg_waveform[:,i,:]
        waveform[i,:,:,:] = waveform[i,:,:,:]

        avg_waveform[:,i,1] = np.roll(avg_waveform[:,i,1], max_lag)
        waveform[i,:,:,1] = np.roll(waveform[i,:,:,1], max_lag, axis = 0)

        if max_lag>0:
            avg_waveform[0:max_lag ,i,1] = np.nan
            waveform[i,0:max_lag ,:,1] = np.nan

        elif max_lag<0:
            avg_waveform[spike_width - 1 + max_lag: ,i,1] = np.nan
            waveform[i,spike_width - 1 + max_lag: ,:,1] = np.nan        


        for cv in [0,1]:
            #shift so both cv, have peak at the prescribed locations
            if peak_time[i,0] != new_peak_loc: #yes the first CV
                shift = int(-( peak_time[i,0] - new_peak_loc))
                if shift !=0:
                    avg_waveform[:,i,cv] = np.roll(avg_waveform[:,i,cv], shift)
                    waveform[i,:,:,cv] = np.roll(waveform[i,:,:,cv], shift, axis = 0)

                    if shift>0:
                        avg_waveform[0:shift ,i,cv] = np.nan
                        waveform[i,0:shift ,:,cv] = np.nan
                    elif shift<0:
                        avg_waveform[spike_width - 1 + shift: ,i,cv] = np.nan
                        waveform[i,spike_width - 1 + shift: ,:,cv] = np.nan    

            tmp_peak = avg_waveform[new_peak_loc,i,cv]

            if np.isnan(tmp_peak) == True : # This is a catch for bad units
                tmp_peak = np.nanmax(avg_waveform[:,i,cv]) 
                #tmp_peak = np.nan
                print(f'This unit {i}, CV {cv} is very likely a bad unit!')
            amplitude[i,cv] = tmp_peak

    return amplitude, waveform, avg_waveform


def get_avg_waveform_per_tp(waveform, channel_pos, d_10, max_site_mean, amplitude, avg_waveform, clus_info, param):
    """
    This function calculates the weighted average waveform per time point, as well as the good time points for each unit

    Parameters
    ----------
    waveform : ndarray (n_units, spike_width. n_channels, cv)
        The waveforms for each unit and cv
    channel_pos : (n_units, 3)
        The spatial position of each unit
    d_10 : ndarray
        The distance at which the signal decays to 10% of the peak
    max_site_mean : ndarray
        The maximum site for each unit, when taking the mean over each cv
    amplitude : ndarray
        The amplitude for each unit and cv
    avg_waveform : ndarray (n_units, spike_width, cv)
        The weighted average waveform over channels
    clus_info : dict
        The clus_info dictionary
    param : dict
        The param dictionary

    Returns
    -------
    ndarrays
        waveform duration, the average waveform per time point and good waveform indexes
    """
    n_units = param['n_units']
    spike_width = param['spike_width']
    n_channels = param['n_channels']
    waveidx = param['waveidx']
    session_id = clus_info['session_id']

    good_site_id = np.empty((n_units,n_channels,2))
    waveform_duration = np.full((n_units,2), np.nan)
    avg_waveform_per_tp = np.full((3, n_units,spike_width,2), np.nan)
    good_wave_idxs = np.zeros((n_units, spike_width, 2))

    for i in range(n_units):   
        for cv in range(2):
            #dist = np.linalg.norm(ChannelPos[SessionID[i]][MaxSite[i,cv],:] - ChannelPos[SessionID[i]][:,:], axis = 1)
            dist = np.linalg.norm(channel_pos[session_id[i]][max_site_mean[i],:] - channel_pos[session_id[i]][:,:], axis = 1)

            test = dist < np.abs(d_10[i,cv])
            good_site_id[i,:,cv] = test
            loc = channel_pos[session_id[i]][good_site_id[i,:,cv].astype(bool), :]

            # select time points where the values are above 25% of the amplitude
            wave_duration_tmp = np.argwhere( np.abs(np.sign(amplitude[i,cv]) * avg_waveform[waveidx,i,cv] )> (np.sign(amplitude[i,cv]) * amplitude[i,cv] * 0.25))

            if wave_duration_tmp.sum != 0:
                try:
                    wave_duration_tmp = np.arange(wave_duration_tmp[0],wave_duration_tmp[-1]+1) + waveidx[0]  
                    wave_duration_tmp = np.expand_dims(wave_duration_tmp, axis = 1)
                    waveform_duration[i,cv] = len(wave_duration_tmp) 
                except:
                    print(f'unit{i} is very likely bad, no good time points in average waveform')  
                    waveform_duration[i,cv] = np.nan

            else:
                print(f'problem with active time points (waveform_duration) for unit {i}')
            
        #Projected location per time point 
            for idx in wave_duration_tmp:
                tmp = np.abs(waveform[i,idx,good_site_id[i,:,cv].astype(bool),cv])
                tmp = np.expand_dims(tmp, axis = 1)
                tmp = np.tile(tmp, (1,3))
                tmp = np.sum(tmp * loc , axis = 0) / np.sum(tmp, axis = 0)
                avg_waveform_per_tp[:,i,idx,cv] = tmp.reshape((3,1))
                good_wave_idxs[i,idx,cv] = 1 
    
    #apply some smoothing the the trajectories
    for uid in range(n_units):
        for cv in range(2):
            for dim in range(avg_waveform_per_tp.shape[0]):
                avg_waveform_per_tp[dim, uid, good_wave_idxs[uid,:,cv].astype(bool), cv] = gaussian_filter(avg_waveform_per_tp[dim, uid, good_wave_idxs[uid,:,cv].astype(bool), cv], 1, radius = 2, axes = 0)
        
    return waveform_duration, avg_waveform_per_tp, good_wave_idxs
        