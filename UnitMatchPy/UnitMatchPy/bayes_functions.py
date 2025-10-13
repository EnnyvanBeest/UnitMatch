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
    # print('Calculating the probability distributions of the metric scores')
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
    # print('Calculating the match probabilities')
    
    score_vector = param['score_vector']
    n_pairs = predictors.shape[0] ** 2
    n_metrics = predictors.shape[2]

    predictors_flat = predictors.reshape(n_pairs, n_metrics)
    
    # Process in chunks
    chunk_size = min(1000, n_pairs)
    min_idx = np.zeros((n_pairs, n_metrics), dtype=int)
    
    for start_idx in range(0, n_pairs, chunk_size):
        end_idx = min(start_idx + chunk_size, n_pairs)
        chunk = predictors_flat[start_idx:end_idx]  # (chunk_size, n_metrics)
        
        # For each metric, find closest score indices for this chunk
        for metric_idx in range(n_metrics):
            metric_values = chunk[:, metric_idx]  # (chunk_size,)
            # Broadcast: (chunk_size, 1) - (1, n_scores) = (chunk_size, n_scores)
            distances = np.abs(metric_values[:, np.newaxis] - score_vector[np.newaxis, :])
            min_idx[start_idx:end_idx, metric_idx] = np.argmin(distances, axis=1)

    likelihood = np.full((n_pairs, len(cond)), np.nan)
    for ck in range(len(cond)):

        metric_indices = np.arange(n_metrics)[np.newaxis, :]  # (1, n_metrics)
        
        # Extract probabilities for all pairs and metrics at once
        tmp_prob = parameter_kernels[min_idx, metric_indices, ck]  # (n_pairs, n_metrics)
        
        # Compute product across metrics for each pair
        likelihood[:, ck] = np.prod(tmp_prob, axis=1)

    denominator = np.nansum(priors[np.newaxis, :] * likelihood, axis=1)  # (n_pairs,)
    
    # Compute final probabilities for all conditions at once
    prob = (priors[np.newaxis, :] * likelihood) / denominator[:, np.newaxis]  # (n_pairs, len(cond))
    
    return prob
