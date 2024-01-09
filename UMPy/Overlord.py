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

    SpatialDecayFit , SpatialDecay,  d_10, AvgCentroid, WeightedAvgWaveF, PeakTime = pf.decay_and_average_Waveform(waveform,ChannelPos, good_idx, MaxSite, MaxSiteMean, param)

    Amplitude, waveform, WeightedAvgWaveF = pf.get_amplitude_shift_Waveform(waveform,WeightedAvgWaveF, PeakTime, param)

    WaveformDuration, WeightedAvgWaveF_PerTP, WaveIdx = pf.avg_Waveform_PerTP(waveform,ChannelPos, d_10, MaxSiteMean, Amplitude, WeightedAvgWaveF, param)

    ExtractedWaveProperties = {'SpatialDecayFit' : SpatialDecayFit, 'SpatialDecay' : SpatialDecay , 'AvgCentroid' : AvgCentroid,
                                'WaveformDuration' : WaveformDuration, 'WeightedAvgWaveF_PerTP' : WeightedAvgWaveF_PerTP, 'WaveIdx' : WaveIdx,
                                 'Amplitude' : Amplitude, 'WeightedAvgWaveF' : WeightedAvgWaveF }
    
    return ExtractedWaveProperties

def extract_metric_scores(ExtractedWaveProperties, SessionSwitch, WithinSession, param, niter  = 2):
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
    AmpScore = mf.get_simple_metric(Amplitude)
    SpatialDecayScore = mf.get_simple_metric(SpatialDecay)
    SpatialDecayFitScore = mf.get_simple_metric(SpatialDecayFit, outlier = True)
    WVcorrScore = mf.get_WVcorr(WeightedAvgWaveF, param)
    WVcorrScore = mf.get_WaveformMSE(WeightedAvgWaveF, param)

    #effected by drift
    for i in range(niter):
        WAW_PerTP_flip = mf.flip_dim(WeightedAvgWaveF_PerTP, param)
        EuclDist = mf.get_Euclidean_dist(WAW_PerTP_flip,param)

        CentroidDist, CentroidVar = mf.Centroid_metrics(EuclDist, param)

        EuclDistRC = mf.get_recentered_Euclidean_dist(WAW_PerTP_flip, AvgCentroid, param)

        CentroidDistRecentered = mf.recentered_metrics(EuclDistRC)
        TrajAngleScore, TrajDistScore = mf.dist_angle(WAW_PerTP_flip, param)


        # Average Euc Dist
        EuclDist = np.nanmin(EuclDist[:,param['PeakLoc'] - param['waveidx'] == 0, :,:].squeeze(), axis = 1 )

        # TotalScore
        IncludeThesePairs = np.argwhere( EuclDist < param['MaxDist']) #array indices of pairs to include
        IncludeThesePairs_idx = np.zeros_like(EuclDist)
        IncludeThesePairs_idx[EuclDist < param['MaxDist']] = 1 

        # Make a dictionary of score to include
        CentroidOverlordScore = (CentroidDistRecentered + CentroidVar) / 2
        WaveformScore = (WVcorrScore + WVcorrScore) / 2
        TrajectoryScore = (TrajAngleScore + TrajDistScore) / 2

        Scores2Include = {'AmpScore' : AmpScore, 'SpatialDecayScore' : SpatialDecayScore, 'CentroidOverlord' : CentroidOverlordScore,
                        'CentroidDist' : CentroidDist, 'WaveformScore' : WaveformScore, 'TrajectoryScore': TrajectoryScore }

        TotalScore, Predictors = mf.get_total_score(Scores2Include, param)

        #Initial thresholding
        if (i < niter - 1):
            #get the thershold for a match
            ThrsOpt = mf.get_threshold(TotalScore, WithinSession, EuclDist, param, IsFirstPass = True)

            param['nExpectedMatches'] = np.sum( (TotalScore > ThrsOpt).astype(int))
            priorMatch = 1 - ( param['nExpectedMatches'] / len(IncludeThesePairs))
            CandidatePairs = TotalScore > ThrsOpt

            drifts, AvgCentroid, WeightedAvgWaveF_PerTP = mf.drift_nSessions(CandidatePairs, SessionSwitch, AvgCentroid, WeightedAvgWaveF_PerTP, TotalScore, param)


    ThrsOpt = mf.get_threshold(TotalScore, WithinSession, EuclDist, param, IsFirstPass = False)
    param['nExpectedMatches'] = np.sum( (TotalScore > ThrsOpt).astype(int))
    priorMatch = 1 - ( param['nExpectedMatches'] / len(IncludeThesePairs))

    ThrsOpt = np.quantile(TotalScore[IncludeThesePairs_idx.astype(bool)], priorMatch)
    CandidatePairs = TotalScore > ThrsOpt

    return TotalScore, CandidatePairs, Scores2Include, Predictors 
