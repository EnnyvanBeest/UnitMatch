import numpy as np

def check_is_in(test_array, parent_array):
    """
    Checks to see if the test_array is contained within the parent_array

    Parameters
    ----------
    test_array : ndarray
        sub array
    parent_array : ndarray
        parent array which may contain the test_array

    Returns
    -------
    bool
        True if the test_array is within the parent_array
    """
    is_in = (test_array[:, None] == parent_array).all(-1).any(-1)
    return is_in


def assign_unique_id(output_prob_array, param, clus_info):
    """
    Assign units to a common group depending on different criteria:
    Conservative - adds units which match with EVERY unit in the proposed group
    Intermediate - adds units which match with EVERY unit in same/adjacent sessions in the proposed group
    Liberal - adds all units which match with any unit in the proposed group
    Each unit will be given a unique group id for each case.

    Parameters
    ----------
    output : ndarray (n_units, n_uunits)
        The 2d probability matrix which gives the UnitMatch probability of a unit with every other unit
    param : dict
        The param dictionary
    clus_info : dict
        The clus_info dictionary

    Returns
    -------
    List
        A list of arrays which gives each unit its group ID for each case
    """
    all_cluster_ids = clus_info['original_ids'] # each units has unique ID

    #create arrays for the unique ids
    unique_id_liberal = np.arange(all_cluster_ids.shape[0])
    ori_unique_id = np.arange(all_cluster_ids.shape[0])
    unique_id_conservative = np.arange(all_cluster_ids.shape[0])
    unique_id = np.arange(all_cluster_ids.shape[0]) #Intermediate Case

    #use a data driven probability threshold
    if param.get('use_data_driven_prob_thrs', False):
        stepsz = 0.1
        bin_edges = np.arange(0,  1 + stepsz, stepsz)
        plot_vec = np.arange(stepsz / 2, 1, stepsz)

        hw, __ = np.histogram(np.diag(output_prob_array), bins = len(bin_edges), density = True)

        threshold = plot_vec[np.diff(hw) > 0.1]
    else:
        threshold = param['match_threshold']

    pairs = np.argwhere(output_prob_array > threshold)
    pairs = np.delete(pairs, np.argwhere(pairs[:,0] == pairs[:,1]), axis =0) #delete self-matches
    pairs = np.sort(pairs, axis = 1)# arange so smaller pairID is in column 1
    #Only keep one copy of pairs only if both CV agree its a match
    pairs_unique, count = np.unique(pairs, axis = 0, return_counts=True)
    pairs_unique_filt = np.delete(pairs_unique, count == 1, axis = 0) #if Count = 1 only 1 CV for that pair!

    #get the mean probability for each match
    prob_mean = np.nanmean(np.vstack((output_prob_array[pairs_unique_filt[:,0], pairs_unique_filt[:,1]], \
                                      output_prob_array[pairs_unique_filt[:,1], pairs_unique_filt[:,0]])), axis=0)
    #sort by the mean probability
    pairs_prob = np.hstack((pairs_unique_filt, prob_mean[:, np.newaxis])) 
    sorted_idxs = np.argsort(-pairs_prob[:,2], axis = 0) #start go in descending order
    pairs_prob_sorted = np.zeros_like(pairs_prob)
    pairs_prob_sorted = pairs_prob[sorted_idxs,:]

    #Create a list which has both copies of each match e.g (1,2) and (2,1) for easier comparison
    pairs_all = np.zeros((pairs_unique_filt.shape[0]*2,2))
    pairs_all[:pairs_unique_filt.shape[0],:] = pairs_unique_filt
    pairs_all[pairs_unique_filt.shape[0]:,:] = pairs_unique_filt[:, (1,0)]

    n_matches_conservative = 0
    n_matches_liberal = 0
    n_matches = 0
    #Go through each pair and assign to groups!!
    for pair in pairs_prob_sorted[:,:2]:
        pair = pair.astype(np.int16)

        #Get the conservative group ID for the current 2 units
        unit_a_conservative_id = unique_id_conservative[pair[0]]
        unit_b_conservative_id = unique_id_conservative[pair[1]]
        # get all units which have the same ID
        same_group_id_a = np.argwhere(unique_id_conservative == unit_a_conservative_id).squeeze()
        same_group_id_b = np.argwhere(unique_id_conservative == unit_b_conservative_id).squeeze()
        #reshape array to be a 1d array if needed
        if len(same_group_id_a.shape) == 0:
            same_group_id_a = same_group_id_a[np.newaxis]
        if len(same_group_id_b.shape) == 0:
            same_group_id_b = same_group_id_b[np.newaxis]

        #will need to check if pair[0] has match with SameGroupIdB and vice versa
        check_pairs_a = np.stack((same_group_id_b, np.broadcast_to(np.array(pair[0]), same_group_id_b.shape)), axis = -1)
        check_pairs_b = np.stack((same_group_id_a, np.broadcast_to(np.array(pair[1]), same_group_id_a.shape)), axis = -1)
        # delete the potential self-matches 
        check_pairs_a = np.delete(check_pairs_a, np.argwhere(check_pairs_a[:,0] == check_pairs_a[:,1]), axis =0)
        check_pairs_b = np.delete(check_pairs_b, np.argwhere(check_pairs_b[:,0] == check_pairs_b[:,1]), axis =0)

        if (np.logical_and(np.all(check_is_in(check_pairs_a, pairs_all)), np.all(check_is_in(check_pairs_b, pairs_all)))):
            #If each pairs matches with every unit in the other pairs group
            #can add as match to all classes
            all_pairs = np.vstack((check_pairs_a, check_pairs_b))
            all_group_idxs = np.unique(all_pairs)
            unique_id_conservative[all_group_idxs] = np.min(unique_id_conservative[all_group_idxs])
            n_matches_conservative +=1

        ##Intermediate matches
        #Now test to see if each pairs match with every unit in the other pair IF they are in the same/adjacent sessions 
        unit_a_id = unique_id[pair[0]]
        unit_b_id = unique_id[pair[1]]

        same_group_id_a = np.argwhere(unique_id == unit_a_id).squeeze()
        same_group_id_b = np.argwhere(unique_id == unit_b_id).squeeze()
        if len(same_group_id_a.shape) == 0:
            same_group_id_a = same_group_id_a[np.newaxis]
        if len(same_group_id_b.shape) == 0:
            same_group_id_b = same_group_id_b[np.newaxis]

        check_pairs_a = np.stack((same_group_id_b, np.broadcast_to(np.array(pair[0]), same_group_id_b.shape)), axis = -1)
        check_pairs_b = np.stack((same_group_id_a, np.broadcast_to(np.array(pair[1]), same_group_id_a.shape)), axis = -1)
        #delete potential self-matches
        check_pairs_a = np.delete(check_pairs_a, np.argwhere(check_pairs_a[:,0] == check_pairs_a[:,1]), axis =0)
        check_pairs_b = np.delete(check_pairs_b, np.argwhere(check_pairs_b[:,0] == check_pairs_b[:,1]), axis =0)

        #check to see if they are in the same or adjacent sessions
        near_session_a = np.abs(np.diff(clus_info['session_id'][check_pairs_a])) <= 1
        near_session_b = np.abs(np.diff(clus_info['session_id'][check_pairs_b])) <= 1

        check_pairs_near_a = check_pairs_a[near_session_a.squeeze()]
        check_pairs_near_b = check_pairs_b[near_session_b.squeeze()]

        #Catch the case where the units ARE NOT in adjacent session, so CheckPairsNear is empty
        if np.logical_and(check_pairs_near_a.size == 0, check_pairs_near_b.size == 0):
            all_pairs = np.vstack((check_pairs_a, check_pairs_b))
            all_group_idxs = np.unique(all_pairs)
            unique_id[all_group_idxs] = np.min(unique_id[all_group_idxs])
            n_matches +=1
        elif (np.logical_and(np.all(check_is_in(check_pairs_near_a, pairs_all)), np.all(check_is_in(check_pairs_near_b, pairs_all)))):
            all_pairs = np.vstack((check_pairs_a, check_pairs_b))
            all_group_idxs = np.unique(all_pairs)
            unique_id[all_group_idxs] = np.min(unique_id[all_group_idxs])
            n_matches +=1

        ## Liberal Matches
        same_group_id_a = np.argwhere(unique_id_liberal == unique_id_liberal[pair[0]]).squeeze()
        same_group_id_b = np.argwhere(unique_id_liberal == unique_id_liberal[pair[1]]).squeeze()

        all_pairs = np.hstack((same_group_id_a, same_group_id_b))
        all_group_idxs = np.unique(all_pairs)
        unique_id_liberal[all_group_idxs] = np.min(unique_id_liberal[all_group_idxs])
        n_matches_liberal +=1
     
    print(f'Number of Liberal Matches: {n_matches_liberal}')
    print(f'Number of Intermediate Matches: {n_matches}')
    print(f'Number of Conservative Matches: {n_matches_conservative}')

    return [unique_id_liberal, unique_id, unique_id_conservative, ori_unique_id]