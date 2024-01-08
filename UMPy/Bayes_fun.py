import numpy as np
import Param_fun as pf

def get_ParameterKernels(Scores2Include, labels, Cond, param, addone = 1):
    """
    Requires Score2Include, a dictionary where the keys are the metric used and the values are
    nUnits * nUnits with the score for each unit.

    Smoothing and add one is done to try and compensate for the fact the histogram used as a prediction for the 
    probability distn has few values, therefore this smoothing hopes to make it more similar to the true distn
    by smoothing nearby peaks and trough to reduce shot noise
    """

    ScoreVector = param['ScoreVector']
    Bins = param['Bins']
    SmoothProb = param['SmoothProb']

    ParameterKernels = np.full((len(ScoreVector), len(Scores2Include), len(Cond)), np.nan)

    ScoreID = 0
    for sc in Scores2Include:
        Scorestmp = Scores2Include[sc]

        SmoothTmp = SmoothProb # Not doing the different ones for now (default the same)


        for Ck in range(len(Cond)):
            
            HistTmp , __ = np.histogram(Scorestmp[labels == Ck], Bins)
            ParameterKernels[:,ScoreID, Ck] = pf.smooth(HistTmp, SmoothTmp)
            ParameterKernels[:,ScoreID, Ck] /= np.sum(ParameterKernels[:,ScoreID,Ck])

            ParameterKernels[:,ScoreID, Ck] += addone* np.min(ParameterKernels[ParameterKernels[:,ScoreID, Ck] !=0, ScoreID, Ck], axis = 0)

        ScoreID +=1    

    return ParameterKernels

def apply_naive_bayes(ParameterKernels,Priors, Predictors, param, Cond):
    """
    Using the Paramater kernels, Priors and Predictors, calculate the probability each pair of units is a 
    match
    """
    ScoreVector = param['ScoreVector']
    
    nPairs = Predictors.shape[0] ** 2

    unravel = np.reshape(Predictors , (Predictors.shape[0] * Predictors.shape[1], Predictors.shape[2],1))
    x1 = np.tile(unravel, ( 1, 1, len(ScoreVector)))
    tmp = np.expand_dims(ScoreVector, axis = (0,1))
    x2 = np.tile(tmp, (x1.shape[0], x1.shape[1], 1))
    minidx = np.argmin( np.abs(x1 - x2), axis = 2)

    likelihood = np.full((nPairs, len(Cond)), np.nan)
    for Ck in range(len(Cond)):
        tmpp = np.zeros_like(minidx, np.float64)
        for yy in range(minidx.shape[1]):
            tmpp[:,yy] = ParameterKernels[minidx[:,yy],yy,Ck]
        likelihood[:,Ck] = np.prod(tmpp, axis=1)


    Probability = np.full((nPairs,2), np.nan )    
    for Ck in range(len(Cond)):
        Probability[:,Ck] = Priors[Ck] * likelihood[:,Ck] / np.nansum((Priors * likelihood), axis =1)
    
    return Probability