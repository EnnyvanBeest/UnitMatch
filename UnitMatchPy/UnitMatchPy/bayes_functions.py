import numpy as np
import UnitMatchPy.param_functions as pf

def get_parameter_kernels(scores_to_include, labels, cond, param, add_one = 1):
    """
    Calculate the probability distribution for the scores for this dataset.
    Smoothing and add one is done to try and compensate for the fact the histogram used as a prediction for the 
    probability distribution has few values, therefore this smoothing hopes to make it more similar to the true distribution
    by smoothing nearby peaks and trough to reduce shot noise.

    Parameters
    ----------
    scores_to_include : dict
        Keys are the metrics used, the values are (n_unit, n_unit) arrays containing the scores for that metric
    labels : ndarray (n_unit, n_unit)
        whether  that unit is a candidate pair or not
    cond : ndarray
        The unique value of labels array
    param : dict
        The param dictionary
    add_one : int, optional
        add a small value to smooth the probability distribution, by default 1

    Returns
    -------
    ndarray
        The probability distributions for each metric
    """
    print('Calculating the probability distributions of the metric scores')
    score_vector = param['score_vector']
    bins = param['bins']
    smooth_prob = param['smooth_prob']

    parameter_kernels = np.full((len(score_vector), len(scores_to_include), len(cond)), np.nan)

    score_id = 0
    for sc in scores_to_include:
        scores_tmp = scores_to_include[sc]

        smooth_tmp = smooth_prob # Not doing the different ones for now (default the same)


        for ck in range(len(cond)):           
            hist_tmp , __ = np.histogram(scores_tmp[labels == ck], bins)
            parameter_kernels[:,score_id, ck] = pf.smooth(hist_tmp, smooth_tmp)
            parameter_kernels[:,score_id, ck] /= np.sum(parameter_kernels[:,score_id,ck])
            parameter_kernels[:,score_id, ck] += add_one* np.min(parameter_kernels[parameter_kernels[:,score_id, ck] !=0, score_id, ck], axis = 0)

        score_id += 1

    return parameter_kernels

def apply_naive_bayes(parameter_kernels, priors, predictors, param, cond):
    """
    Using the parameter_kernels, priors and predictors to calculate the probability each pair of 
    units is a match

    Parameters
    ----------
    parameter_kernels : ndarray
        The probability distributions for each metric
    priors : ndarray 
        The prior probability a pair is a match or not
    predictors : ndarray (n_units, n_units, n_metrics)
        The combined array of all of the metrics for all of the units
    param : dict
        The param dictionary
    cond : ndarray
        The possibilities e.g is it a match or not

    Returns
    -------
    ndarray
        The probability the unit is or is not a match
    """
    print('Calculating the match probabilities')
    
    score_vector = param['score_vector']
    n_pairs = predictors.shape[0] ** 2

    unravel = np.reshape(predictors , (predictors.shape[0] * predictors.shape[1], predictors.shape[2],1))
    x1 = np.tile(unravel, ( 1, 1, len(score_vector)))
    tmp = np.expand_dims(score_vector, axis = (0,1))
    x2 = np.tile(tmp, (x1.shape[0], x1.shape[1], 1))
    min_idx = np.argmin( np.abs(x1 - x2), axis = 2)

    likelihood = np.full((n_pairs, len(cond)), np.nan)
    for ck in range(len(cond)):
        tmp_prob = np.zeros_like(min_idx, np.float64)
        for yy in range(min_idx.shape[1]):
            tmp_prob[:,yy] = parameter_kernels[min_idx[:,yy],yy,ck]
        likelihood[:,ck] = np.prod(tmp_prob, axis=1)


    prob = np.full((n_pairs,2), np.nan )
    for ck in range(len(cond)):
        prob[:,ck] = priors[ck] * likelihood[:,ck] / np.nansum((priors * likelihood), axis =1)
    
    return prob
