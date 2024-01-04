import Param_fun as pf
import Metrics_fun as mf
import Bayes_fun as bf
import utils as util
import numpy as np

def extract_parameters(waveform, ChannelPos, param):
    """
    This function runs all of the parameter extraction in one step!
    """
    
    waveform = pf.detrend_waveform(waveform)

    MaxSite, good_idx, good_pos, MaxSiteMean = pf.get_max_sites(waveform, ChannelPos, param)

    SpatialDecayFit , SpatialDecay,  d_10, AvgCentroid, WeightedAvgWaveF, PeakTime = pf.decay_and_average_WaveF(waveform,ChannelPos, good_idx, MaxSite, MaxSiteMean, param)

    Amplitude, waveform, WeightedAvgWaveF = pf.get_amplitude_shift_WaveF(waveform,WeightedAvgWaveF, PeakTime, param)

    WaveformDuration, WeightedAvgWaveF_PerTP, WaveIdx = pf.avg_WaveF_PerTP(waveform,ChannelPos, d_10, MaxSiteMean, Amplitude, WeightedAvgWaveF, param)

    ExtractedWaveProperties = {'SpatialDecayFit' : SpatialDecayFit, 'SpatialDecay' : SpatialDecay , 'AvgCentroid' : AvgCentroid,
                                'WaveformDuration' : WaveformDuration, 'WeightedAvgWaveF_PerTP' : WeightedAvgWaveF_PerTP, 'WaveIdx' : WaveIdx,
                                 'Amplitude' : Amplitude, 'WeightedAvgWaveF' : WeightedAvgWaveF }
    
    return ExtractedWaveProperties

def extract_metric_scores(ExtractedWaveProperties, SessionSwitch, WithinSession, param, n_iterations  = 2):
    """
    This function requires extract parameters and session information, to extract the standard set of scores.
    n_itterations is the number of passes through the extraction and then drift correction, if n_iterations = 1 no drift correection is done
    for n_iterations > 2 drift correction is done multiple times. n_iterations = 2 is reccomened.

    ** if needed can change so you input a list which choose what scores to include... 
    """

    #unpack need arrays from the ExtractedWaveProperties dictionary
    Amplitude = ExtractedWaveProperties['Amplitude']
    SpatialDecay = ExtractedWaveProperties['SpatialDecay']
    SpatialDecayFit = ExtractedWaveProperties['SpatialDecayFit']
    WeightedAvgWaveF = ExtractedWaveProperties['WeightedAvgWaveF']
    WeightedAvgWaveF_PerTP = ExtractedWaveProperties['WeightedAvgWaveF_PerTP']
    AvgCentroid = ExtractedWaveProperties['AvgCentroid']

    #These scores are NOT effected by the drift correction
    Amp_score = mf.get_simple_metric(Amplitude)
    SD_score = mf.get_simple_metric(SpatialDecay)
    SDF_score = mf.get_simple_metric(SpatialDecayFit, outlier = True)
    WVcorr_score = mf.get_WVcorr(WeightedAvgWaveF, param)
    WFMSE_score = mf.get_WaveFormMSE(WeightedAvgWaveF, param)

    #effectedt by drift
    for i in range(n_iterations):
        WAW_PerTP_flip = mf.flip_dim(WeightedAvgWaveF_PerTP, param)
        EuclDist = mf.get_Euclidean_dist(WAW_PerTP_flip,param)

        CentroidDist, CentroidVar = mf.Centroid_metrics(EuclDist, param)

        EuclDist_2 = mf.get_recentered_Euclidean_dist(WAW_PerTP_flip, AvgCentroid, param)

        CentroidDistRecentered = mf.recentered_metrics(EuclDist_2)
        TrajAngle_score, TrajDist_score = mf.dist_angle(WAW_PerTP_flip, param)


        # Average Euc Dist
        EuclDist = np.nanmin(EuclDist[:,param['PeakLoc'] - param['waveidx'] == 0, :,:].squeeze(), axis = 1 )

        # TotalScore
        IncludeThesePairs = np.argwhere( EuclDist < param['MaxDist']) #array indices of pairs to include
        IncludeThesePairs_idx = np.zeros_like(EuclDist)
        IncludeThesePairs_idx[EuclDist < param['MaxDist']] = 1 

        # Make a dictionary of score to include
        CentroidOverlord_score = (CentroidDistRecentered + CentroidVar) / 2
        Waveform_score = (WVcorr_score + WFMSE_score) / 2
        Trajectory_score = (TrajAngle_score + TrajDist_score) / 2

        Scores2Include = {'Amp_score' : Amp_score, 'SD_score' : SD_score, 'CentroidOverlord' : CentroidOverlord_score,
                        'CentroidDist' : CentroidDist, 'Waveform_score' : Waveform_score, 'Trajectory_score': Trajectory_score }

        TotalScore, Predictors = mf.get_total_score(Scores2Include, param)

        #Initial thresholding
        if (i < n_iterations - 1):
            #get the thershold for a match
            ThrsOpt = mf.get_threshold(TotalScore, WithinSession, EuclDist, param, is_first_pass = True)

            param['nExpectedMatches'] = np.sum( (TotalScore > ThrsOpt).astype(int))
            priorMatch = 1 - ( param['nExpectedMatches'] / len(IncludeThesePairs))
            CandidatePairs = TotalScore > ThrsOpt

            drifts, AvgCentroid, WeightedAvgWaveF_PerTP = mf.drift_n_days(CandidatePairs, SessionSwitch, AvgCentroid, WeightedAvgWaveF_PerTP, TotalScore, param)


    ThrsOpt = mf.get_threshold(TotalScore, WithinSession, EuclDist, param, is_first_pass = False)
    param['nExpectedMatches'] = np.sum( (TotalScore > ThrsOpt).astype(int))
    priorMatch = 1 - ( param['nExpectedMatches'] / len(IncludeThesePairs))

    ThrsOpt = np.quantile(TotalScore[IncludeThesePairs_idx.astype(bool)], priorMatch)
    CandidatePairs = TotalScore > ThrsOpt

    return TotalScore, CandidatePairs, Scores2Include, Predictors 
