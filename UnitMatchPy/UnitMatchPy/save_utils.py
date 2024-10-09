import os
import pickle
import pandas as pd
import numpy as np

def make_match_table(scores_to_include, matches, output_prob, total_score, output_threshold, clus_info, param, UIDs = None, matches_curated = None):
    """
    Creates the match table showing information for every pair of units.

    Parameters
    ----------
    scores_to_include : dictionary
        The dictionary containing the used scores and the values for each pair of units
    matches : list
        The list of all matches
    output_prob : ndarray (n_units, n_units)
        The array containing the UnitMatch match probability for each pair of units for cv 1/2 and 2/1
    total_score : ndarray (n_units, n_units)
        The total score for each pair of units being a match for each pair of units for cv 1/2 and 2/1
    output_threshold : int
        The probability threshold value for deciding matches
    clus_info : dict
        The clus_info dictionary
    param : dict
        The param dictionary
    UIDs : list, optional
        The list of unique ids for different case from assign_unique_id, by default None
    matches_curated : list, optional
        A list of matches manually curated using the GUI, by default None

    Returns
    -------
    dataframe
        A pandas dataframe as a match table for each pair of units
    """
    # Making Match Table
    n_units = param['n_units']

    #Give UMID add it as well?
    # xx, yy =np.meshgrid(np.arange(nUnits), np.arange(nUnits))
    # UnitA = np.reshape(xx, (nUnits * nUnits)).astype(np.int16)
    # UnitB = np.reshape(yy, (nUnits * nUnits)).astype(np.int16)

    original_ids = clus_info['original_ids'].squeeze()
    xx, yy = np.meshgrid(original_ids, original_ids)
    unit_a_list = xx.reshape(n_units*n_units)
    unit_b_list = yy.reshape(n_units*n_units)

    session_id = clus_info['session_id']
    xx, yy = np.meshgrid(session_id, session_id)
    unit_a_session_list = xx.reshape(n_units*n_units) + 1 # Add one here so it counts from one not 0
    unit_b_session_list = yy.reshape(n_units*n_units) + 1

    all_matches = np.reshape(output_threshold, (n_units*n_units)).astype(np.int8) # Uses Matches currated .. as well
    total_score_list = np.reshape(total_score, (n_units*n_units))
    prob_list = np.reshape(output_prob, (n_units*n_units))


    #create the initial array with the important info
    #see if there is a curated list and if soo add it
    if matches_curated is not None:
        matches_curated_list = np.zeros((n_units,n_units))
        for match in matches:
            matches_curated_list[match[0], match[1]] = 1

        matches_curated_list = np.reshape(matches_curated, (n_units * n_units)).astype(np.int8)

        df = pd.DataFrame(np.array([unit_a_list, unit_b_list, unit_a_session_list, unit_b_session_list, all_matches, matches_curated_list, prob_list, total_score_list]).T, columns = ['ID1', 'ID2', 'RecSes 1', 'RecSes 2', 'Matches', 'Matches Currated', 'UM Probabilities', 'TotalScore'])

    else:
        df = pd.DataFrame(np.array([unit_a_list, unit_b_list, unit_a_session_list, unit_b_session_list, all_matches, prob_list, total_score_list]).T, columns = ['ID1', 'ID2', 'RecSes 1', 'RecSes 2', 'Matches', 'UM Probabilities', 'TotalScore'])

    #add a dictionary to the match table
    for key, value in scores_to_include.items():
        df[key] = np.reshape(value, (n_units *n_units)).T

    #if you have supplied UIDs create a data frame using them and merge it to the save table
    if UIDs is not None:
        unique_id_liberal = UIDs[0]
        unique_id = UIDs[1]
        unique_id_conservative = UIDs[2]
        original_unique_id = UIDs[3]

        xx, yy = np.meshgrid(unique_id_liberal, unique_id_liberal)
        unit_a_liberal_id = xx.reshape(n_units*n_units)
        unit_b_liberal_id = yy.reshape(n_units*n_units)

        xx, yy = np.meshgrid(original_unique_id, original_unique_id)
        unit_a_original_id = xx.reshape(n_units*n_units)
        unit_b_original_id = yy.reshape(n_units*n_units)

        xx, yy = np.meshgrid(unique_id_conservative, unique_id_conservative)
        unit_a_conservative_id = xx.reshape(n_units*n_units)
        unit_b_conservative_id = yy.reshape(n_units*n_units)

        xx, yy = np.meshgrid(unique_id, unique_id)
        unit_a_int_id = xx.reshape(n_units*n_units)
        unit_b_int_id = yy.reshape(n_units*n_units)

        unique_id_df = pd.DataFrame(np.array([unit_a_original_id, unit_b_original_id, unit_a_liberal_id, unit_b_liberal_id, unit_a_int_id, unit_b_int_id, unit_a_conservative_id, unit_b_conservative_id]).T, columns = ['UID1', 'UID2', 'UID Liberal 1', 'UID Liberal 2', 'UID int 1', 'UM UID int 2', 'UID Conservative 1', 'UID Conservative 2'])
        df = df.join(unique_id_df)

    return df

def save_to_output(save_dir, scores_to_include, matches, output_prob, avg_centroid, avg_waveform, avg_waveform_per_tp, max_site,
                   total_score, output_threshold, clus_info, param, UIDs = None, matches_curated = None, save_match_table = True):
    """
    Saves all useful information calculated by UnitMatch to a given save_dir

    Parameters
    ----------
    save_dir : str
        The absolute path to the save directory
    scores_to_include : dictionary
        The dictionary containing the used scores and the values for each pair of units
    matches : list
        The list of all matches
    output_prob : ndarray (n_units, n_units)
        The array containing the UnitMatch match probability for each pair of units for cv 1/2 and 2/1
    avg_centroid : _type_
        _description_
    avg_waveform : ndarray (n_units, spike_width, cv)
        The weighted average waveform over channels
    avg_waveform_per_tp : ndarray
        The average waveform per time point
    max_site : ndarray
        The maximum site for each unit
    total_score : ndarray (n_units, n_units)
        The total score for each pair of units being a match for each pair of units for cv 1/2 and 2/1
    output_threshold : int
        The probability threshold value for deciding matches
    clus_info : dict
        The clus_info dictionary
    param : dict
        The param dictionary
        _description_
    UIDs : list, optional
        The list of unique ids for different case from assign_unique_id, by default None
    matches_curated : list, optional
        A list of matches manually curated using the GUI, by default None
    save_match_table : bool, optional
        If True will save a match table containing the information for every pair of units in a table, by default True
    """

    #Choose a file where the save directory will be made
    #options for create and overwrite?
    if os.path.isdir(save_dir) == False:
        os.mkdir(save_dir)

    #save scores
    UM_scores_path = os.path.join(save_dir, 'UM Scores')
    np.savez(UM_scores_path, **scores_to_include)


    # #save ClusInfo
    clus_info_path = os.path.join(save_dir, 'ClusInfo.pickle')
    with open(clus_info_path, 'wb') as fp:
        pickle.dump(clus_info, fp)

    #Save param
    param_path = os.path.join(save_dir, 'UMparam.pickle')
    with open(param_path, 'wb') as fp:
        pickle.dump(param, fp)

    #Save output
    match_prob_path = os.path.join(save_dir, 'MatchProb')
    #MAY WANT TO CHANGE TOSAVE PROB FOR BOTH CV AND AVG?
    np.save(match_prob_path, output_prob)

    #Save Waveform info
    waveform_info = {"avg_centroid" : avg_centroid, "avg_waveform" : avg_waveform, "avg_waveform_per_tp" : avg_waveform_per_tp, 
                    "max_site" : max_site}
    waveform_info_path = os.path.join(save_dir, 'WaveformInfo')
    np.savez(waveform_info_path, **waveform_info)

    #save autimatuc matches
    matches_path = os.path.join(save_dir, 'Matches')
    np.save(matches_path, matches)

    if matches_curated != None:
        matches_curated_path = os.path.join(save_dir, 'Matches Currated')
        np.save(matches_curated_path, matches_curated)

    if save_match_table == True:
        df = make_match_table(scores_to_include, matches, output_prob, total_score, output_threshold, clus_info, param, UIDs = UIDs, matches_curated = None)
        match_table_path = os.path.join(save_dir, 'MatchTable.csv')
        df.to_csv(match_table_path, index = False)

def save_to_output_seperate_CV(save_dir, scores_to_include, matches, output_prob, avg_centroid, avg_waveform, avg_waveform_per_tp, max_site,
                   total_score, match_threshold, clus_info, param, UIDs = None, matches_curated = None, save_match_table = True):
    """
    Saves all useful information calculated by UnitMatch to a given save_dir.
    This function will split the n_unit * n_units arrays into separate arrays depending on which pair of CV is used

    Parameters
    ----------
    save_dir : str
        The absolute path to the save directory
    scores_to_include : dictionary
        The dictionary containing the used scores and the values for each pair of units
    matches : list
        The list of all matches
    output_prob : ndarray (n_units, n_units)
        The array containing the UnitMatch match probability for each pair of units for cv 1/2 and 2/1
    avg_centroid : _type_
        _description_
    avg_waveform : ndarray (n_units, spike_width, cv)
        The weighted average waveform over channels
    avg_waveform_per_tp : ndarray
        The average waveform per time point
    max_site : ndarray
        The maximum site for each unit
    total_score : ndarray (n_units, n_units)
        The total score for each pair of units being a match for each pair of units for cv 1/2 and 2/1
    output_threshold : int
        The probability threshold value for deciding matches
    clus_info : dict
        The clus_info dictionary
    param : dict
        The param dictionary
        _description_
    UIDs : list, optional
        The list of unique ids for different case from assign_unique_id, by default None
    matches_curated : list, optional
        A list of matches manually curated using the GUI, by default None
    save_match_table : bool, optional
        If True will save a match table containing the infomation for every pair of units in a table, by default True
    """


    #Start by separating the info into CV
    matches12_part1 = np.argwhere(np.tril(output_prob) > match_threshold) 
    matches12_part2 = np.argwhere(np.tril(output_prob).T > match_threshold)
    matches12 = np.unique(np.concatenate((matches12_part1,matches12_part2)), axis = 0)

    matches21_part1 = np.argwhere(np.triu(output_prob) > match_threshold) 
    matches21_part2 = np.argwhere(np.triu(output_prob).T > match_threshold)
    matches21 = np.unique(np.concatenate((matches21_part1,matches21_part2)), axis = 0)

    output12_tmp1 = np.tril(output_prob)
    output12_tmp2 = np.tril(output_prob).T
    np.fill_diagonal(output12_tmp2, 0)
    output12 = output12_tmp1 + output12_tmp2

    output21_tmp1= np.triu(output_prob)
    output21_tmp2 = np.triu(output_prob).T
    np.fill_diagonal(output21_tmp2, 0)
    output21 = output21_tmp1 + output21_tmp2

    scores_to_include12 = {}
    scores_to_include21 = {}
    for key, value in scores_to_include.items():
        tmp1 = np.tril(value)
        tmp2 = np.tril(value).T
        np.fill_diagonal(tmp2, 0)

        tmp3 = np.triu(value)
        tmp4 = np.triu(value).T
        np.fill_diagonal(tmp4, 0)

        scores_to_include12[key] = tmp1 + tmp2 
        scores_to_include21[key] = tmp3 + tmp4  

    #Now can save all of these like above
    #Choose a file where the save directory will be made
    #options for create and overwrite?
    if os.path.isdir(save_dir) == False:
        os.mkdir(save_dir)

    #save scores
    UM_scores_path_cv12 = os.path.join(save_dir, 'UM Scores CV12')
    np.savez(UM_scores_path_cv12, **scores_to_include12)
    UM_scores_path_cv21 = os.path.join(save_dir, 'UM Scores CV21')
    np.savez(UM_scores_path_cv21, **scores_to_include21)

    # #save ClusInfo
    clus_info_path = os.path.join(save_dir, 'ClusInfo.pickle')
    with open(clus_info_path, 'wb') as fp:
        pickle.dump(clus_info, fp)

    #Save param
    param_path = os.path.join(save_dir, 'UMparam.pickle')
    with open(param_path, 'wb') as fp:
        pickle.dump(param, fp)

    #Save output n_unit*n_units probability array
    match_prob_path_cv12 = os.path.join(save_dir, 'MatchProb CV12')
    np.save(match_prob_path_cv12, output12)
    match_prob_path_cv21 = os.path.join(save_dir, 'MatchProb CV21')
    np.save(match_prob_path_cv21, output21)

    #Save Waveform info
    waveform_info = {"avg_centroid" : avg_centroid, "avg_waveform" : avg_waveform, "avg_waveform_per_tp" : avg_waveform_per_tp,
                    "max_site" : max_site}
    wavefrom_info_path = os.path.join(save_dir, 'WaveformInfo')
    np.savez(wavefrom_info_path, **waveform_info)

    #save automatic matches
    matches_path_cv12 = os.path.join(save_dir, 'Matches CV12')
    np.save(matches_path_cv12, matches12)
    matches_path_cv21 = os.path.join(save_dir, 'Matches CV21')
    np.save(matches_path_cv21, matches21)

    if matches_curated != None:
        MatchesCuratedPath = os.path.join(save_dir, 'Matches Currated')
        np.save(MatchesCuratedPath, matches_curated)

    output_threshold = np.zeros_like(output_prob)
    output_threshold[output_prob > match_threshold] = 1

    if save_match_table == True:
        df = make_match_table(scores_to_include, matches, output_prob, total_score, output_threshold, clus_info, param, UIDs = UIDs, matches_curated = None)
        match_table_path = os.path.join(save_dir, 'MatchTable.csv')
        df.to_csv(match_table_path, index = False)


def load_output(save_dir, load_match_table = True):
    """
    Will load all the information in the save directory

    Parameters
    ----------
    save_dir : str
        The absolute path to the save directory
    load_match_table : bool, optional
        If True will load in the match table, by default True

    Returns
    -------
    All data saved in the UM save directory
    """

    #load scores
    UM_scores_path = os.path.join(save_dir, 'UM Scores.npz')
    UM_scores = dict(np.load(UM_scores_path))

    #load ClusInfo
    clus_info_path = os.path.join(save_dir, 'ClusInfo.pickle')
    with open(clus_info_path, 'rb') as fp:
        clus_info = pickle.load(fp)

    #load param
    param_path = os.path.join(save_dir, 'UMparam.pickle')
    with open(param_path, 'rb') as fp:
        param = pickle.load(fp)

    #load output
    match_prob_path = os.path.join(save_dir, 'MatchProb.npy')
    match_prob = np.load(match_prob_path)

    
    #Load Waveform info
    wavefrom_info_path = os.path.join(save_dir, 'WaveformInfo.npz')
    wavefrom_info =dict(np.load(wavefrom_info_path))

    if load_match_table == True:
        match_table_path = os.path.join(save_dir, 'MatchTable.csv') 
        match_table = pd.read_csv(match_table_path)
    
        return UM_scores, clus_info, param, match_prob, wavefrom_info, match_table
    return UM_scores, clus_info, param, wavefrom_info, match_prob

def load_output_separate_CV(save_dir, load_match_table = True):
    """
    Will load all the information in the save directory.
    This function will load in data if the CV are saved separately

    Parameters
    ----------
    save_dir : str
        The absolute path to the save directory
    load_match_table : bool, optional
        If True will load in the match table, by default True

    Returns
    -------
    All data saved in the UM save directory
    """

    UM_scores_path_cv12 = os.path.join(save_dir, 'UM Scores CV12.npz')
    UM_scores_cv12 = dict(np.load(UM_scores_path_cv12))
    UM_scores_path_cv21 = os.path.join(save_dir, 'UM Scores CV21.npz')
    UM_scores_cv21 = dict(np.load(UM_scores_path_cv21))

    # #Load ClusInfo
    clus_info_path = os.path.join(save_dir, 'ClusInfo.pickle')
    with open(clus_info_path, 'rb') as fp:
        clus_info = pickle.load(fp)

    #Load param
    param_path = os.path.join(save_dir, 'UMparam.pickle')
    with open(param_path, 'rb') as fp:
        param = pickle.load(fp)

    #Load output nUnit*nUnits probabilities array
    match_prob_path_cv12 = os.path.join(save_dir, 'MatchProb CV12.npy')
    match_prob_cv12 = np.load(match_prob_path_cv12)
    match_prob_path_cv21 = os.path.join(save_dir, 'MatchProb CV21.npy')
    match_prob_cv21 = np.load(match_prob_path_cv21)

    #Load Waveform info
    wavefrom_info_path = os.path.join(save_dir, 'WaveformInfo.npz')
    wavefrom_info =dict(np.load(wavefrom_info_path))

    #save automatic matches
    matches_path_cv12 = os.path.join(save_dir, 'Matches CV12.npy')
    matches_cv12 = np.load(matches_path_cv12)
    matches_path_cv21 = os.path.join(save_dir, 'Matches CV21.npy')
    matches_cv21 = np.load(matches_path_cv21)

    if load_match_table == True:
        match_table_path = os.path.join(save_dir, 'MatchTable.csv') 
        match_table = pd.read_csv(match_table_path)

        return UM_scores_cv12, UM_scores_cv21, clus_info, param, match_prob_cv12, match_prob_cv21,matches_cv12, matches_cv21, wavefrom_info, match_table
    return UM_scores_cv12, UM_scores_cv21, clus_info, param, match_prob_cv12, match_prob_cv21,matches_cv12, matches_cv21, wavefrom_info

def save_prob_for_phy(probability, param, clus_info):
    """
    Saves the within session UnitMatch probabilities for each session in their KiloSort directory,
    to be used with the UnitMatch Phy plugin 

    Parameters
    ----------
    Probability : ndarray (nUnits, nUnits)
        The calculates UnitMatch probability array
    param : dictionary
        The param dictionary
    ClusInfo : dictionary
        The ClusInfo dictionary
    """

    session_switch = clus_info['session_switch']
    n_units_per_session = param['n_units_per_session']

    for sid in range(session_switch.shape[0] - 1):
        #file to save the array in
        save_file_tmp = os.path.join(param['KSdirs'][sid], 'probability_templates.npy')

        matrix_prob = np.full((n_units_per_session[sid], n_units_per_session[sid]), np.nan) #Make the size of all the units

        session_output = probability[session_switch[sid]:session_switch[sid+1], session_switch[sid]:session_switch[sid+1]]
        #If Only good units where used add values know to a matrix of NaNs 
        if session_switch[sid+1] - session_switch[sid] != n_units_per_session[sid]:
            all_good_units = clus_info['good_units'][sid].squeeze().astype(int)
            for id, gid in enumerate(clus_info['good_units'][sid].astype(int)):
                matrix_prob[gid, all_good_units] = session_output[id,:]
        else:
            matrix_prob = session_output

        np.save(save_file_tmp, matrix_prob)
