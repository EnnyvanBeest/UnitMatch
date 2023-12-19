import numpy as np
import Param_fun as pf

def get_ParameterKernels(Scores2Include, labels, Cond, param, addone = 1):
    """
    ******* I think it would probaly be nice to input alot of these parameters as a dictionary,
    ******* which can be automatical created with defaut values

    Mainly requires Score2Include, a dictionary where the keys are the metric used and the values are
    n_unts * n_units with the score for each unit.

    Smoothing and addone is done to try and compensate for the fact the histogram used as a prediction for the 
    probabilty distn has few values, therefore this smoothing hopesto make it more similar to the true distn
    by smoothing nearby peaks and trough to reduce shot nosie
    """

    ScoreVector = param['ScoreVector']
    Bins = param['Bins']
    SmoothProb = param['SmoothProb']

    ParamaterKernels = np.full((len(ScoreVector), len(Scores2Include), len(Cond)), np.nan)

    sc_no = 0
    for sc in Scores2Include:
        Scorestmp = Scores2Include[sc]

        SmoothTmp = SmoothProb # Notdoingthe different ones for now (default the same)


        for Ck in range(len(Cond)):
            
            hist_tmp , __ = np.histogram(Scorestmp[labels == Ck], Bins)
            ParamaterKernels[:,sc_no, Ck] = pf.smooth(hist_tmp, SmoothTmp)
            ParamaterKernels[:,sc_no, Ck] /= np.sum(ParamaterKernels[:,sc_no,Ck])

            ParamaterKernels[:,sc_no, Ck] += addone* np.min(ParamaterKernels[ParamaterKernels[:,sc_no, Ck] !=0, sc_no, Ck], axis = 0)

        sc_no +=1    

    return ParamaterKernels

def apply_naive_bayes(ParamaterKernels,Priors, Predictors, param, Cond):
    """
    Using the Paramter knernels, Priors and Predictors, calculatethe proabiblty each pair of units is a 
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
            tmpp[:,yy] = ParamaterKernels[minidx[:,yy],yy,Ck]
        likelihood[:,Ck] = np.prod(tmpp, axis=1)


    Probabilty = np.full((nPairs,2), np.nan )    
    for Ck in range(len(Cond)):
        Probabilty[:,Ck] = Priors[Ck] * likelihood[:,Ck] / np.nansum((Priors * likelihood), axis =1)
    
    return Probabilty