import os
import pickle
import pandas as pd
import numpy as np

def make_match_table(Scores2Include, Matches, Output, TotalScore, OutputThreshold, ClusInfo, param, MatchesCurated = None):
    # Making Match Table
    nUnits = param['nUnits']

    #Give UMID add it as well?
    # xx, yy =np.meshgrid(np.arange(nUnits), np.arange(nUnits))
    # UnitA = np.reshape(xx, (nUnits * nUnits)).astype(np.int16)
    # UnitB = np.reshape(yy, (nUnits * nUnits)).astype(np.int16)

    OriginalIDs = ClusInfo['OriginalID'].squeeze()
    xx, yy = np.meshgrid(OriginalIDs, OriginalIDs)
    UnitALS = xx.reshape(nUnits*nUnits)
    UnitBLS = yy.reshape(nUnits*nUnits)

    SessionIDs = ClusInfo['SessionID']
    xx, yy = np.meshgrid(SessionIDs, SessionIDs)
    UnitAsessionLS = xx.reshape(nUnits*nUnits) + 1 # Add onehere so it counts from one not 0
    UnitBsessionLS = yy.reshape(nUnits*nUnits) + 1

    ALLMatches = np.reshape(OutputThreshold, (nUnits*nUnits)).astype(np.int8) # Uses Matches currated .. as well
    TotalScoreLS = np.reshape(TotalScore, (nUnits*nUnits))
    ProbLS = np.reshape(Output, (nUnits*nUnits))

    #add euclidiean distance

    #create the intial array wiht theimportant info
    #see if there is a curated lsit and if soo add it
    if MatchesCurated != None:
        MatchesCuratedLS = np.zeros((nUnits,nUnits))
        for match in Matches:
            MatchesCuratedLS[match[0], match[1]] = 1

        MatchesCuratedLS = np.reshape(MatchesCurated, (nUnits * nUnits)).astype(np.int8)

        df = pd.DataFrame(np.array([UnitALS, UnitBLS, UnitAsessionLS, UnitBsessionLS, ALLMatches, MatchesCuratedLS, ProbLS, TotalScoreLS]).T, columns = ['ID1', 'ID2', 'RecSes 1', 'RecSes 2', 'Matches', 'Matches Currated', 'UM Probabilities', 'TotalScore'])

    else:
        df = pd.DataFrame(np.array([UnitALS, UnitBLS, UnitAsessionLS, UnitBsessionLS, ALLMatches, ProbLS, TotalScoreLS]).T, columns = ['ID1', 'ID2', 'RecSes 1', 'RecSes 2', 'Matches', 'UM Probabilities', 'TotalScore'])

    #add a dictionary to the match table
    for key, value in Scores2Include.items():
        df[key] = np.reshape(value, (nUnits *nUnits)).T

    return df

def save_to_output(SaveDir, Scores2Include, Matches, Output, AvgCentroid, AvgWaveform, AvgWaveformPerTP, MaxSite,
                   TotalScore, OutputThreshold, ClusInfo, param, MatchesCurated = None, SaveMatchTable = True):
    """ 
    Save all of the useful infomation into a SaveDIR.
    keyword arguments are MatchesCurated and SaveMatchTable, Supply MatchesCurrated of you have currated matches and
    pass SaveMatchTable = False, to not save the .csv match table as this takes the most time
    
    """

    #Choose a file where the save directory will be made
    #options for create and overwrite?
    if os.path.isdir(SaveDir) == False:
        os.mkdir(SaveDir)

    #save scores
    UMscoresPath = os.path.join(SaveDir, 'UM Scores')
    np.savez(UMscoresPath, **Scores2Include)


    # #save ClusInfo
    ClusInfoPath = os.path.join(SaveDir, 'ClusInfo.pickle')
    with open(ClusInfoPath, 'wb') as fp:
        pickle.dump(ClusInfo, fp)

    #Save param
    ParamPath = os.path.join(SaveDir, 'UMparam.pickle')
    with open(ParamPath, 'wb') as fp:
        pickle.dump(param, fp)

    #Save output
    MatchProbPath = os.path.join(SaveDir, 'MatchProb')
    #MAY WANT TO CHANGE TOSAVE PROB FOR BOTH CV AND AVG?
    np.save(MatchProbPath, Output)

    #Save Waveform info
    WaveformInfo = {"AvgCentroid" : AvgCentroid, "AvgWaveform" : AvgWaveform, "AvgWaveformPerTP" : AvgWaveformPerTP, 
                    "MaxSite" : MaxSite}
    WavefromInfoPath = os.path.join(SaveDir, 'WaveformInfo')
    np.savez(WavefromInfoPath, **WaveformInfo)

    #save autimatuc matches
    MatchesPath = os.path.join(SaveDir, 'Matches')
    np.save(MatchesPath, Matches)

    if MatchesCurated != None:
        MatchesCuratedPath = os.path.join(SaveDir, 'Matches Currated')
        np.save(MatchesCuratedPath, MatchesCurated)

    if SaveMatchTable == True:
        df = make_match_table(Scores2Include, Matches, Output, TotalScore, OutputThreshold, ClusInfo, param, MatchesCurated = None)
        MatchTablePath = os.path.join(SaveDir, 'MatchTable.csv')
        df.to_csv(MatchTablePath, index = False)

def save_to_output_seperate_CV(SaveDir, Scores2Include, Matches, Output, AvgCentroid, AvgWaveform, AvgWaveformPerTP, MaxSite,
                   TotalScore, MatchThreshold, ClusInfo, param, MatchesCurated = None, SaveMatchTable = True):

    #Start by sperating the info into CV
    Matches12part1 = np.argwhere(np.tril(Output) > MatchThreshold) 
    Matches12part2 = np.argwhere(np.tril(Output).T > MatchThreshold)
    Matches12 = np.unique(np.concatenate((Matches12part1,Matches12part2)), axis = 0)

    Matches21part1 = np.argwhere(np.triu(Output) > MatchThreshold) 
    Matches21part2 = np.argwhere(np.triu(Output).T > MatchThreshold)
    Matches21 = np.unique(np.concatenate((Matches21part1,Matches21part2)), axis = 0)

    Output12tmp1 = np.tril(Output)
    Output12tmp2 = np.tril(Output).T
    np.fill_diagonal(Output12tmp2, 0)
    Output12 = Output12tmp1 + Output12tmp2

    Output21tmp1= np.triu(Output)
    Output21tmp2 = np.triu(Output).T
    np.fill_diagonal(Output21tmp2, 0)
    Output21 = Output21tmp1 + Output21tmp2

    Scores2Include12 = {}
    Scores2Include21 = {}
    for key, value in Scores2Include.items():
        tmp1 = np.tril(value)
        tmp2 = np.tril(value).T
        np.fill_diagonal(tmp2, 0)

        tmp3 = np.triu(value)
        tmp4 = np.triu(value).T
        np.fill_diagonal(tmp4, 0)

        Scores2Include12[key] = tmp1 + tmp2 
        Scores2Include21[key] = tmp3 + tmp4  

    #Now can save all of these like above
    #Choose a file where the save directory will be made
    #options for create and overwrite?
    if os.path.isdir(SaveDir) == False:
        os.mkdir(SaveDir)

    #save scores
    UMscoresPathCV12 = os.path.join(SaveDir, 'UM Scores CV12')
    np.savez(UMscoresPathCV12, **Scores2Include12)
    UMscoresPathCV21 = os.path.join(SaveDir, 'UM Scores CV21')
    np.savez(UMscoresPathCV21, **Scores2Include21)

    # #save ClusInfo
    ClusInfoPath = os.path.join(SaveDir, 'ClusInfo.pickle')
    with open(ClusInfoPath, 'wb') as fp:
        pickle.dump(ClusInfo, fp)

    #Save param
    ParamPath = os.path.join(SaveDir, 'UMparam.pickle')
    with open(ParamPath, 'wb') as fp:
        pickle.dump(param, fp)

    #Save output nUnit*nUnits probabilite array
    MatchProbPathCV12 = os.path.join(SaveDir, 'MatchProb CV12')
    np.save(MatchProbPathCV12, Output12)
    MatchProbPathCV21 = os.path.join(SaveDir, 'MatchProb CV21')
    np.save(MatchProbPathCV21, Output21)

    #Save Waveform info
    WaveformInfo = {"AvgCentroid" : AvgCentroid, "AvgWaveform" : AvgWaveform, "AvgWaveformPerTP" : AvgWaveformPerTP, 
                    "MaxSite" : MaxSite}
    WavefromInfoPath = os.path.join(SaveDir, 'WaveformInfo')
    np.savez(WavefromInfoPath, **WaveformInfo)

    #save autimatuc matches
    MatchesPathCV12 = os.path.join(SaveDir, 'Matches CV12')
    np.save(MatchesPathCV12, Matches12)
    MatchesPathCV21 = os.path.join(SaveDir, 'Matches CV21')
    np.save(MatchesPathCV21, Matches21)

    if MatchesCurated != None:
        MatchesCuratedPath = os.path.join(SaveDir, 'Matches Currated')
        np.save(MatchesCuratedPath, MatchesCurated)

    OutputThreshold = np.zeros_like(Output)
    OutputThreshold[Output > MatchThreshold] = 1

    if SaveMatchTable == True:
        df = make_match_table(Scores2Include, Matches, Output, TotalScore, OutputThreshold, ClusInfo, param, MatchesCurated = None)
        MatchTablePath = os.path.join(SaveDir, 'MatchTable.csv')
        df.to_csv(MatchTablePath, index = False)


def load_output(SaveDir, LoadMatchTable = False):

    #load scores
    UMscoresPath = os.path.join(SaveDir, 'UM Scores.npz')
    UMScores = dict(np.load(UMscoresPath))

    #load ClusInfo
    ClusInfoPath = os.path.join(SaveDir, 'ClusInfo.pickle')
    with open(ClusInfoPath, 'rb') as fp:
        ClusInfo = pickle.load(fp)

    #load param
    ParamPath = os.path.join(SaveDir, 'UMparam.pickle')
    with open(ParamPath, 'rb') as fp:
        param = pickle.load(fp)

    #load output
    MatchProbPath = os.path.join(SaveDir, 'MatchProb.npy')
    MatchProb = np.load(MatchProbPath)

    
    #Load Waveform info
    WavefromInfoPath = os.path.join(SaveDir, 'WaveformInfo.npz')
    WavefromInfo =dict(np.load(WavefromInfoPath))

    if LoadMatchTable == True:
        MatchTablePath = os.path.join(SaveDir, 'MatchTable.csv') 
        MatchTable = pd.read_csv(MatchTablePath)
    
        return UMScores, ClusInfo, param, MatchProb, WavefromInfo, MatchTable
    return UMScores, ClusInfo, param, WavefromInfo, MatchProb

def load_output_seperate_CV(SaveDir, LoadMatchTable = False):

    UMscoresPathCV12 = os.path.join(SaveDir, 'UM Scores CV12.npz')
    UMScoresCV12 = dict(np.load(UMscoresPathCV12))
    UMscoresPathCV21 = os.path.join(SaveDir, 'UM Scores CV21.npz')
    UMScoresCV21 = dict(np.load(UMscoresPathCV21))

    # #Load ClusInfo
    ClusInfoPath = os.path.join(SaveDir, 'ClusInfo.pickle')
    with open(ClusInfoPath, 'rb') as fp:
        ClusInfo = pickle.load(fp)

    #Load param
    ParamPath = os.path.join(SaveDir, 'UMparam.pickle')
    with open(ParamPath, 'rb') as fp:
        param = pickle.load(fp)

    #Load output nUnit*nUnits probabilite array
    MatchProbPathCV12 = os.path.join(SaveDir, 'MatchProb CV12.npy')
    MatchProbCV12 = np.load(MatchProbPathCV12)
    MatchProbPathCV21 = os.path.join(SaveDir, 'MatchProb CV21.npy')
    MatchProbCV21 = np.load(MatchProbPathCV21)

    #Load Waveform info
    WavefromInfoPath = os.path.join(SaveDir, 'WaveformInfo.npz')
    WavefromInfo =dict(np.load(WavefromInfoPath))

    #save autimatuc matches
    MatchesPathCV12 = os.path.join(SaveDir, 'Matches CV12.npy')
    MatchesCV12 = np.load(MatchesPathCV12)
    MatchesPathCV21 = os.path.join(SaveDir, 'Matches CV21.npy')
    MatchesCV21 = np.load(MatchesPathCV21)

    if LoadMatchTable == True:
        MatchTablePath = os.path.join(SaveDir, 'MatchTable.csv') 
        MatchTable = pd.read_csv(MatchTablePath)

        return UMScoresCV12, UMScoresCV21, ClusInfo, param, MatchProbCV12, MatchProbCV21,MatchesCV12, MatchesCV21, WavefromInfo, MatchTable
    return UMScoresCV12, UMScoresCV21, ClusInfo, param, MatchProbCV12, MatchProbCV21,MatchesCV12, MatchesCV21, WavefromInfo