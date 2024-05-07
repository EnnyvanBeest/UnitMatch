import numpy as np

def check_is_in(TestArray, ParentArray):
    """
    Test to se if a row in TestArray is in ParrentArray
    Arguments:
        TestArray -- ndarray (N, 2)
        ParentArray -- ndArray(M, 2)

    Returns:
        IsIn -- ndarray (N) dtype - bool
    """    
    IsIn = (TestArray[:, None] == ParentArray).all(-1).any(-1)
    return IsIn

def AssignUID(Output, param, ClusInfo):

    AllClusterIDs = ClusInfo['OriginalID'] # each units has unique ID

    #create arrays for the uniwue ids
    UniqueIDLiberal = np.arange(AllClusterIDs.shape[0])
    OriUniqueID = np.arange(AllClusterIDs.shape[0])
    UniqueIDConservative = np.arange(AllClusterIDs.shape[0])
    UniqueID = np.arange(AllClusterIDs.shape[0])

    GoodRecSesID = ClusInfo['SessionID']
    RecOpt = np.unique(ClusInfo['SessionID'])
    nRec = RecOpt.shape[0]

    #data Driven Threshold?
    if param.get('UseDataDrivenProbThrs', False):
        stepsz = 0.1
        binedges = np.arange(0,  1 + stepsz, stepsz)
        plotvec = np.arange(stepsz / 2, 1, stepsz)

        hw, __ = np.histogram(np.diag(Output), bins = len(binedges), density = True)

        Threshold = plotvec[np.diff(hw) > 0.1]
    else:
        Threshold = param['MatchThreshold']

    Pairs = np.argwhere(Output > Threshold)
    Pairs = np.delete(Pairs, np.argwhere(Pairs[:,0] == Pairs[:,1]), axis =0) #delete self-matches
    Pairs = np.sort(Pairs, axis = 1)# arange so smaller pairID is in column 1
    #Only keep one copy of pairs only if both CV agrree its a match
    PairsUnique, Count = np.unique(Pairs, axis = 0, return_counts=True)
    PairsUniqueFilt = np.delete(PairsUnique, Count == 1, axis = 0) #if Count = 1 only 1 CV for that pair!

    #get the mean probabilty for each match
    ProbMean = np.nanmean(np.vstack((Output[PairsUniqueFilt[:,0], PairsUniqueFilt[:,1]], Output[PairsUniqueFilt[:,1], PairsUniqueFilt[:,0]])), axis=0)
    #sort by the mean probabilty
    PairsProb = np.hstack((PairsUniqueFilt, ProbMean[:, np.newaxis])) 
    SortIdxs = np.argsort(-PairsProb[:,2], axis = 0) #start go in decending order
    PairsProbSorted = np.zeros_like(PairsProb)
    PairsProbSorted = PairsProb[SortIdxs,:]

    #Create a list which has both copies of each match e.g (1,2) and (2,1) for easier comparisson 
    PairsAll = np.zeros((PairsUniqueFilt.shape[0]*2,2))
    PairsAll[:PairsUniqueFilt.shape[0],:] = PairsUniqueFilt
    PairsAll[PairsUniqueFilt.shape[0]:,:] = PairsUniqueFilt[:, (1,0)]

    nMatchesConservative = 0
    nMatchesLiberal = 0
    nMatches = 0
    #Go through each pair and assign to groups!!
    for pair in PairsProbSorted[:,:2]:
        pair = pair.astype(np.int16)

        #Get the conservative group ID for thecurrent 2 units
        UnitAConservativeID = UniqueIDConservative[pair[0]]
        UnitBConservativeID = UniqueIDConservative[pair[1]]
        # get all units which have the same ID
        SameGroupIdA = np.argwhere(UniqueIDConservative == UnitAConservativeID).squeeze()
        SameGroupIdB = np.argwhere(UniqueIDConservative == UnitBConservativeID).squeeze()
        #reshape array to be a 1d array if needed
        if len(SameGroupIdA.shape) == 0:
            SameGroupIdA = SameGroupIdA[np.newaxis]
        if len(SameGroupIdB.shape) == 0:
            SameGroupIdB = SameGroupIdB[np.newaxis]

        #will need to check if pair[0] has match with SameGroupIdB and vice versa
        CheckPairsA = np.stack((SameGroupIdB, np.broadcast_to(np.array(pair[0]), SameGroupIdB.shape)), axis = -1)
        CheckPairsB = np.stack((SameGroupIdA, np.broadcast_to(np.array(pair[1]), SameGroupIdA.shape)), axis = -1)
        # delete the potential self-matches 
        CheckPairsA = np.delete(CheckPairsA, np.argwhere(CheckPairsA[:,0] == CheckPairsA[:,1]), axis =0)
        CheckPairsB = np.delete(CheckPairsB, np.argwhere(CheckPairsB[:,0] == CheckPairsB[:,1]), axis =0)

        if (np.logical_and(np.all(check_is_in(CheckPairsA, PairsAll)), np.all(check_is_in(CheckPairsB, PairsAll)))):
            #If each pairs matches with every unit in the other pairs group
            #can add as match to all classes
            UniqueIDConservative[pair[0]] = min(UniqueIDConservative[pair])
            UniqueIDConservative[pair[1]] = min(UniqueIDConservative[pair])
            nMatchesConservative +=1

            UniqueID[pair[0]] = min(UniqueID[pair])
            UniqueID[pair[1]] = min(UniqueID[pair])
            nMatches +=1

            UniqueIDLiberal[pair[0]] = min(UniqueIDLiberal[pair])
            UniqueIDLiberal[pair[0]] = min(UniqueIDLiberal[pair])
            nMatchesLiberal +=1
        else:
            #Now test to see if each pairs match with every unit in the other pair IF they are in the same/adjacent sessions 
            UnitAID = UniqueID[pair[0]]
            UnitBID = UniqueID[pair[1]]

            SameGroupIdA = np.argwhere(UniqueID == UnitAID).squeeze()
            SameGroupIdB = np.argwhere(UniqueID == UnitBID).squeeze()
            if len(SameGroupIdA.shape) == 0:
                SameGroupIdA = SameGroupIdA[np.newaxis]
            if len(SameGroupIdB.shape) == 0:
                SameGroupIdB = SameGroupIdB[np.newaxis]

            CheckPairsA = np.stack((SameGroupIdB, np.broadcast_to(np.array(pair[0]), SameGroupIdB.shape)), axis = -1)
            CheckPairsB = np.stack((SameGroupIdA, np.broadcast_to(np.array(pair[1]), SameGroupIdA.shape)), axis = -1)
            #delte potential self-matches
            CheckPairsA = np.delete(CheckPairsA, np.argwhere(CheckPairsA[:,0] == CheckPairsA[:,1]), axis =0)
            CheckPairsB = np.delete(CheckPairsB, np.argwhere(CheckPairsB[:,0] == CheckPairsB[:,1]), axis =0)

            #check to see if they are in the same or adjacent sessions
            InNearSessionA = np.abs(np.diff(ClusInfo['SessionID'][CheckPairsA])) <= 1
            InNearSessionB = np.abs(np.diff(ClusInfo['SessionID'][CheckPairsB])) <= 1

            CheckPairsNearA = CheckPairsA[InNearSessionA.squeeze()]
            CheckPairsNearB = CheckPairsB[InNearSessionB.squeeze()]

            if (np.logical_and(np.all(check_is_in(CheckPairsNearA, PairsAll)), np.all(check_is_in(CheckPairsNearB, PairsAll)))):
                UniqueID[pair[0]] = min(UniqueID[pair])
                UniqueID[pair[1]] = min(UniqueID[pair])
                nMatches +=1

                UniqueIDLiberal[pair[0]] = min(UniqueIDLiberal[pair])
                UniqueIDLiberal[pair[1]] = min(UniqueIDLiberal[pair])
                nMatchesLiberal +=1
            else:
                UniqueIDLiberal[pair[0]] = min(UniqueIDLiberal[pair])
                UniqueIDLiberal[pair[1]] = min(UniqueIDLiberal[pair])
                nMatchesLiberal +=1
                
    print(f'Number of Liberal Matches: {nMatchesLiberal}')
    print(f'Number of Intermediate Matches: {nMatches}')
    print(f'Number of Conservative Matches: {nMatchesConservative}')
    return [UniqueIDLiberal, UniqueID, UniqueIDConservative, OriUniqueID]