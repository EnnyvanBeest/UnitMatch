from tkinter import *
from tkinter import ttk
import matplotlib.pyplot as plt
from matplotlib.figure import Figure 
from matplotlib.backends.backend_tkagg import (FigureCanvasTkAgg,  
NavigationToolbar2Tk) 
import numpy as np


def run_GUI():
    global CVtkinter
    global root
    global EntryA
    global EntryB
    global SessionEntryA
    global SessionEntryB
    global MatchIdx
    global FrameTable
    global ScoreTable
    global AvgWaveformPlot
    global TrajectoryPlot
    global BayesLabel
    global OriginalIDLabel
    global RawWaveformPlot
    global HistPlot
    global IsMatch
    global NotMatch
    global OptionA
    global EntryFrame
    
    np.set_printoptions(suppress = True)
    IsMatch = []
    NotMatch = []

    root = Tk()
    # download theme from https://sourceforge.net/projects/tcl-awthemes/
    root.tk.call('lappend', 'auto_path', r'TkinterTheme\awthemes-10.4.0')
    root.tk.call('package', 'require', 'awdark')
    s = ttk.Style(root)
    s.theme_use('awdark')
    root.title('UMPy - Manual Curation')
    #root.geometry('800x800')
    root.iconbitmap(r'UM Logo.ico')
    background = ttk.Frame(root)
    background.place(x=0, y=0, relwidth=1.0, relheight=1.0)

    FrameTable = ttk.LabelFrame(root)
    ScoreTable = ttk.LabelFrame(root)
    AvgWaveformPlot = Canvas(root) 
    TrajectoryPlot = Canvas(root)
    BayesLabel = ttk.Label(root)
    OriginalIDLabel = ttk.Label(root)
    RawWaveformPlot = Canvas(root)
    HistPlot = Canvas(root)

    #Unit entry
    ######################################################################################
    #Will have Unit A - color green and unit B- color Blue
    EntryFrame = ttk.LabelFrame(root, text = 'Select Units')
    LabelA = ttk.Label(EntryFrame, text = 'Unit A')
    LabelA.configure(foreground='green')
    LabelB = ttk.Label(EntryFrame, text = 'Unit B')
    LabelB.configure(foreground='blue')

    #select the session
    SessionsList = np.arange(1,param['nSessions']+1).tolist()
    SessionEntryA = ttk.Combobox(EntryFrame, value = SessionsList, width = 2)
    SessionEntryB = ttk.Combobox(EntryFrame, value = SessionsList, width = 2)
    SessionEntryA.set(1) #Start wiht session1,2
    SessionEntryB.set(2)
    LabelSessionA = ttk.Label(EntryFrame, text = 'Session No.')
    LabelSessionB = ttk.Label(EntryFrame, text = 'Session No.')


    #select CV
    CVoptions = [('Avg', 0),  ('(1,2)', 1), ('(2,1)',2)]
    CVtkinter = IntVar()
    CVtkinter.set(0)
    LabelCV = ttk.Label(EntryFrame, text = 'Select the cv option')
    for i, option in enumerate(CVoptions):
        RadioCV = ttk.Radiobutton(EntryFrame, text = option[0], value =option[1], variable = CVtkinter, command = update_unit_cv).grid(row = i+1, column = 0)



    #selecting the unit
    SessionA = int(SessionEntryA.get())
    SessionB = int(SessionEntryB.get()) 
    CV = get_cv_option()
    CVoption = CVtkinter.get() - 1 #want 0 for cv (1,2) , want 1 for cv (2,1)
    if CV == 'Avg':
        #Use average CV values
        TmpIdxA = np.argwhere((SessionSwitch[SessionA -1] <= MatchesAvg[:,0]) * (MatchesAvg[:,0] < SessionSwitch[SessionA]) == True ).squeeze()
        TmpIdxB = np.argwhere((SessionSwitch[SessionB -1] <= MatchesAvg[:,1]) * (MatchesAvg[:,1] < SessionSwitch[SessionB]) == True).squeeze()
        InBoth = np.isin(TmpIdxA, TmpIdxB)
        OptionA = MatchesAvg[TmpIdxA[InBoth],:].tolist()

        EntryA = ttk.Combobox(EntryFrame, values = OptionA, width = 10)
        EntryA.set(OptionA[0][0])

        OptionB = np.flip(np.argsort(OutputAvg[int(EntryA.get()),SessionSwitch[SessionB-1]:SessionSwitch[SessionB]]) + SessionSwitch[SessionB-1])
        OptionB = OptionB.tolist()

    else :    
        TmpIdxA = np.argwhere((SessionSwitch[SessionA -1] <= MatchesGUI[CVoption][:,0]) * ( MatchesGUI[CVoption][:,0] < SessionSwitch[SessionA]) == True ).squeeze()
        TmpIdxB = np.argwhere((SessionSwitch[SessionB -1] <= MatchesGUI[CVoption][:,1]) * ( MatchesGUI[CVoption][:,1] < SessionSwitch[SessionB]) == True).squeeze()
        InBoth = np.isin(TmpIdxA, TmpIdxB)

        OptionA =  MatchesGUI[CVoption][TmpIdxA[InBoth],:].tolist()
        if CVoption == 0:
            OptionA = sorted(OptionA)

        EntryA = ttk.Combobox(EntryFrame, values = OptionA, width = 10)
        EntryA.set(OptionA[0][0])


        OptionB = np.flip(np.argsort(OutputGUI[CVoption][int(EntryA.get()),SessionSwitch[SessionB-1]:SessionSwitch[SessionB]]) + SessionSwitch[SessionB-1])
        OptionB = OptionB.tolist()

    ##################
    EntryB = ttk.Combobox(EntryFrame, values = OptionB, width = 10 )
    EntryB.set(OptionB[0])
    EntryA.bind('<<ComboboxSelected>>', Update_Units)
    EntryB.bind('<<ComboboxSelected>>', update)

    SessionEntryA.bind('<<ComboboxSelected>>', update_unit_entryA)
    SessionEntryB.bind('<<ComboboxSelected>>', update_unit_entryB)

    #adding a button which swaps unit A and B
    SwapButton = ttk.Button(EntryFrame, text = 'Swap Units', command = swap_units)

    #Clalculate the score histrograms
    #for each CV pair
    #make global varibel so the functions can acces the histograms.
    global HistNamesAvg
    global HistNames12
    global HistNames21
    global HistAvg
    global Hist12
    global Hist21
    global HistMatchesAvg
    global HistMatches12
    global HistMatches21

    HistNamesAvg, HistAvg, HistMatchesAvg =  get_score_histrograms(Scores2IncludeAvg, (OutputAvg > MatchThreshold))
    HistNames12, Hist12, HistMatches12 =  get_score_histrograms(Scores2IncludeGUI[0], (OutputGUI[0] > MatchThreshold))
    HistNames21, Hist21, HistMatches21 =  get_score_histrograms(Scores2IncludeGUI[1], (OutputGUI[1] > MatchThreshold))

    update(None)
    MatchIdx = 0
    #place the widgets on the EntryFrame
    LabelCV.grid(row = 0, column = 0)
    LabelA.grid(row = 0, column = 1)
    LabelB.grid(row = 0, column = 3)
    LabelSessionA.grid(row = 1, column = 1)
    SessionEntryA.grid(row = 1, column = 2, padx = 15)
    LabelSessionB.grid(row = 1, column = 3)
    SessionEntryB.grid(row = 1, column = 4, padx = 15)
    EntryA.grid(row = 2, column = 1, columnspan = 2, stick = 'WE', padx = 5)
    EntryB.grid(row = 2, column = 3, columnspan = 2, sticky = 'WE', padx = 5)
    SwapButton.grid(row = 3, column = 2, columnspan = 2, sticky = 'WE')
    ######################################################################################


    #MatchButtons
    ######################################################################################
    MatchButton = ttk.Button(root, text = 'Set as Match', command = set_match)
    NonMatchButton = ttk.Button(root, text='Set as Non Match', command = set_not_match)



    # Set up Key-Board shortcuts
    root.bind_all('u', update)
    root.bind_all('<Return>', update)
    root.bind_all('<Right>', next_pair)
    root.bind_all('<Left>', previous_pair)
    root.bind_all('q', set_match)
    root.bind_all('m', set_match)
    root.bind_all('e', set_not_match)
    root.bind_all('n', set_not_match)

    #Grid the units
    EntryFrame.grid(row = 0, column = 0, pady=5, padx = 5)
    MatchButton.grid(row = 4, column = 3, sticky = 'E', padx = 5)
    NonMatchButton.grid(row = 4, column = 4, sticky = 'W', padx = 5)

    root.mainloop()

    return IsMatch, NotMatch

def process_info_for_GUI(Output, MatchThresholdIn, Scores2Include, TotalScore, AmplitudeIn, SpatialDecayIn,
                         AvgCentroidIn, AvgWaveformIn, AvgWaveformPerTPIn, WaveIdxIn, MaxSiteIn, MaxSiteMeanIn, 
                         waveformIn, WithinSessionIn, ChannelPosIn, ClusInfoIn, paramIn):
    """ 
    This function has two jobs.
    1. Process data so it is in the correct form for the GUI e.g sperate matricies for CV options.
    2. Pass variable and make them global in the scope of GUI.py so all function in GUI.py can use them.
    """
    global MatchesAvg
    global MatchesGUI
    global OutputAvg
    global OutputGUI
    global Scores2IncludeAvg
    global Scores2IncludeGUI
    global Amplitude
    global AmplitudeAvg
    global SpatialDecay
    global SpatialDecayAvg
    global AvgCentroid
    global AvgCentroidAvg
    global AvgWaveform
    global AvgWaveformAvg
    global AvgWaveformPerTP
    global AvgWaveformPerTPAvg
    global WaveIdx
    global MaxSite
    global MaxSiteMean
    global waveform
    global ClusInfo
    global SessionSwitch
    global param
    global MatchThreshold
    global WithinSession
    global ChannelPos

    Amplitude = AmplitudeIn
    SpatialDecay = SpatialDecayIn
    AvgCentroid = AvgCentroidIn
    AvgWaveform = AvgWaveformIn
    AvgWaveformPerTP = AvgWaveformPerTPIn
    WaveIdx = WaveIdxIn
    MaxSite = MaxSiteIn
    MaxSiteMean = MaxSiteMeanIn
    waveform = waveformIn
    ClusInfo = ClusInfoIn
    SessionSwitch = ClusInfo['SessionSwitch']
    param = paramIn
    MatchThreshold = MatchThresholdIn
    WithinSession = WithinSessionIn
    ChannelPos = ChannelPosIn

    OutputThreshold = np.zeros_like(Output)
    OutputThreshold[Output > MatchThreshold] = 1

    matches = np.argwhere(OutputThreshold == 1) # need all matches including same sessoin for GUI

    # The correct output for different
    #for the code it is helpful for Mathces to contain both (session 1, session 2) and (session 2, session 1)
    # so when changin Unit A/B session it can find matches for all permutations
    Matches12part1 = np.argwhere(np.tril(Output) > MatchThreshold) 
    Matches12part2 = np.argwhere(np.tril(Output).T > MatchThreshold)
    Matches12 = np.unique(np.concatenate((Matches12part1,Matches12part2)), axis = 0)


    Matches21part1 = np.argwhere(np.triu(Output) > MatchThreshold) 
    Matches21part2 = np.argwhere(np.triu(Output).T > MatchThreshold)
    Matches21 = np.unique(np.concatenate((Matches21part1,Matches21part2)), axis = 0)

    MatchesGUI = [Matches12, Matches21]

    OutputGUI1part1 = np.tril(Output)
    OutputGUI1part2 = np.tril(Output).T
    np.fill_diagonal(OutputGUI1part2, 0)
    OutputGUI1 = OutputGUI1part1 + OutputGUI1part2

    OutputGUI2part1 = np.triu(Output)
    OutputGUI2part2 = np.triu(Output).T
    np.fill_diagonal(OutputGUI2part2, 0)
    OutputGUI2 = OutputGUI2part1 + OutputGUI2part2

    OutputGUI = [OutputGUI1, OutputGUI2]

    Scores2IncludeGUI = []
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

    Scores2IncludeGUI = [Scores2Include12, Scores2Include21]


    #gettin avg CV data
    # for the Scores where can do (X + X.T)/2 andtake upper triangular part
    TotalScoreAvg = np.triu( (TotalScore + TotalScore.T) / 2)

    Scores2IncludeAvg = {}
    for key, value in Scores2Include.items():
        Scores2IncludeAvg[key] = (value + value.T) / 2

    OutputAvg = (OutputGUI[0] + OutputGUI[1]) / 2

    MatchesAvgPart1 = np.argwhere( OutputAvg > MatchThreshold )
    MatchesAvgPart2 = np.argwhere( OutputAvg.T > MatchThreshold )
    MatchesAvg = np.unique(np.concatenate((MatchesAvgPart1,MatchesAvgPart2)), axis = 0)

    # or an simply average over both CV
    AmplitudeAvg = np.mean(Amplitude, axis = -1)
    SpatialDecayAvg = np.mean(SpatialDecay, axis = -1)
    AvgCentroidAvg = np.mean(AvgCentroid, axis = -1)
    AvgWaveformAvg = np.mean(AvgWaveform, axis = -1)
    AvgWaveformPerTPAvg = np.mean(AvgWaveformPerTP, axis = -1)



def update(event):
        
    UnitA = int(EntryA.get())
    UnitB = int(EntryB.get())

    CV = get_cv_option() # CV = 'Avg', if AVG is selcted else it equal [0,1] or [1,0]
    CVoption = CVtkinter.get() - 1

    table = get_table_data(UnitA, UnitB, CV)
    MakeTable(table)

    ScoresTable = get_unit_score_table(UnitA, UnitB, CVoption)
    make_unit_score_table(ScoresTable)

    plot_AvgWaveforms(UnitA, UnitB, CV)
    plot_trajectories(UnitA, UnitB, CV)
    plot_raw_waveforms(UnitA, UnitB, CV)
    add_probability_label(UnitA, UnitB, CVoption)
    add_original_ID(UnitA, UnitB)

    #plot histograms based of off the CV
    if CVoption == -1:
        plot_Histograms(HistNamesAvg, HistAvg, HistMatchesAvg, Scores2IncludeAvg, UnitA, UnitB)
    if CVoption == 0:
        plot_Histograms(HistNames12, Hist12, HistMatches12, Scores2IncludeGUI[CVoption], UnitA, UnitB)
    if CVoption == 1:
        plot_Histograms(HistNames21, Hist21, HistMatches21, Scores2IncludeGUI[CVoption], UnitA, UnitB)


#sort CV function as to be calleed as part of update
def get_cv_option():
    """
    Will read in the values of the radio button and assign an appropriate value to CV.
    In general is a list where itis [Unit A cv, UnitB cv], however it could be the string Avg 
    """
    global CVtkinter
    ChosenOption = CVtkinter.get()
    if ChosenOption == 0:
        CV = 'Avg'
    elif ChosenOption == 1:
        CV = [0, 1]
    elif ChosenOption == 2:
        CV = [1, 0] 
    return CV

#These arethe function used to selct units,including how the CV selctition radio buttons, session selction, unitselection andmoving left and rigthfor next units
def update_unit_cv():
    """
    When updating the CV we need to do /not do the following:
    - keep the same selected units,
    - update the options in boxes for unitA and unitB as matches and likely mathces units can change
    - make it so when scrolling the list #MATCHIDX AUTOMATICALLY UPDATES TO THE CORECT POINT INTHE NEW CV
    - update the screen to show the new CV
    """
    global EntryA
    global EntryB
    global OptionA
    global SessionEntryA
    global SessionEntryB
    global MatchIdx
    
    #selecting the unit
    SessionA = int(SessionEntryA.get())
    SessionB = int(SessionEntryB.get()) 

    #Keep track of the unit it was before as we dont want to change theunit viewed when changing the CV
    EntryAtmp = int(EntryA.get())
    EntryBtmp = int(EntryB.get())

    if EntryA.winfo_exists() == 1:
        EntryA.destroy()
    if EntryB.winfo_exists() == 1:
        EntryB.destroy()

    CV = get_cv_option()
    CVoption = CVtkinter.get() - 1 #want 0 for cv (1,2) , want 1 for cv (2,1)
    if CV == 'Avg':
        #Use average CV values
        TmpIdxA = np.argwhere((SessionSwitch[SessionA -1] <= MatchesAvg[:,0]) * (MatchesAvg[:,0] < SessionSwitch[SessionA]) == True ).squeeze()
        TmpIdxB = np.argwhere((SessionSwitch[SessionB -1] <= MatchesAvg[:,1]) * (MatchesAvg[:,1] < SessionSwitch[SessionB]) == True).squeeze()
        InBoth = np.isin(TmpIdxA, TmpIdxB)
        OptionA = MatchesAvg[TmpIdxA[InBoth],:].tolist()

        EntryA = ttk.Combobox(EntryFrame, values = OptionA, width = 10)
        EntryA.set(EntryAtmp)

        OptionB = np.flip(np.argsort(OutputAvg[int(EntryA.get()),SessionSwitch[SessionB-1]:SessionSwitch[SessionB]]) + SessionSwitch[SessionB-1])
        OptionB = OptionB.tolist()

    else :

        TmpIdxA = np.argwhere((SessionSwitch[SessionA -1] <= MatchesGUI[CVoption][:,0]) * ( MatchesGUI[CVoption][:,0] < SessionSwitch[SessionA]) == True).squeeze()
        TmpIdxB = np.argwhere((SessionSwitch[SessionB -1] <= MatchesGUI[CVoption][:,1]) * ( MatchesGUI[CVoption][:,1] < SessionSwitch[SessionB]) == True).squeeze()
        InBoth = np.isin(TmpIdxA, TmpIdxB)

        #So is orderby unit A
        OptionA =  MatchesGUI[CVoption][TmpIdxA[InBoth],:].tolist()
        if CVoption == 0:
            OptionA = sorted(OptionA)

        EntryA = ttk.Combobox(EntryFrame, values = OptionA, width = 10)
        EntryA.set(EntryAtmp)

        ##NOT SORTED FOR CV
        OptionB = np.flip(np.argsort(OutputGUI[CVoption][int(EntryA.get()),SessionSwitch[SessionB-1]:SessionSwitch[SessionB]]) + SessionSwitch[SessionB-1])
        OptionB = OptionB.tolist()

    ##################
    EntryB = ttk.Combobox(EntryFrame, values = OptionB, width = 10 )
    EntryB.set(EntryBtmp)
    
    EntryA.bind('<<ComboboxSelected>>', Update_Units)
    EntryB.bind('<<ComboboxSelected>>', update)
    EntryA.grid(row = 2, column = 1, columnspan = 2, stick = 'WE', padx = 5)
    EntryB.grid(row = 2, column = 3, columnspan = 2, sticky = 'WE', padx = 5)

    tmpList = [int(EntryAtmp), int(EntryBtmp)]
    if tmpList in OptionA:
        MatchIdx = OptionA.index(tmpList)

    update(None)


def update_unit_entryA(event):
    """
    This will be called when the user changes the sesion A, and the following will happen:
    - Keep the same B unit
    - Set A to the first match of the new session 
    - Set the Match IDX at 0, for the new list of match pairs
    - update the options in list A
    - update the options for list B
    - update what is in the screen
    """
    global EntryA
    global EntryB
    global OptionA
    global MatchIdx
    global SessionEntryB
    
    SessionA = int(SessionEntryA.get())
    SessionB = int(SessionEntryB.get()) 
    EntryBtmp = int(EntryB.get())

    if EntryA.winfo_exists() == 1:
        EntryA.destroy()
    if EntryB.winfo_exists() == 1:
        EntryB.destroy()

    CV = get_cv_option()
    CVoption = CVtkinter.get() - 1 #want 0 for cv (1,2) , want 1 for cv (2,1)
    if CV == 'Avg':
        #Use average CV values
        TmpIdxA = np.argwhere((SessionSwitch[SessionA -1] <= MatchesAvg[:,0]) * (MatchesAvg[:,0] < SessionSwitch[SessionA]) == True ).squeeze()
        TmpIdxB = np.argwhere((SessionSwitch[SessionB -1] <= MatchesAvg[:,1]) * (MatchesAvg[:,1] < SessionSwitch[SessionB]) == True).squeeze()
        InBoth = np.isin(TmpIdxA, TmpIdxB)
        OptionA = MatchesAvg[TmpIdxA[InBoth],:].tolist()

        EntryA = ttk.Combobox(EntryFrame, values = OptionA, width = 10)
        EntryA.set(OptionA[0][0])

        OptionB = np.flip(np.argsort(OutputAvg[int(EntryA.get()),SessionSwitch[SessionB-1]:SessionSwitch[SessionB]]) + SessionSwitch[SessionB-1])
        OptionB = OptionB.tolist()

    else :

        TmpIdxA = np.argwhere((SessionSwitch[SessionA -1] <= MatchesGUI[CVoption][:,0]) * ( MatchesGUI[CVoption][:,0] < SessionSwitch[SessionA]) == True ).squeeze()
        TmpIdxB = np.argwhere((SessionSwitch[SessionB -1] <= MatchesGUI[CVoption][:,1]) * ( MatchesGUI[CVoption][:,1] < SessionSwitch[SessionB]) == True).squeeze()
        InBoth = np.isin(TmpIdxA, TmpIdxB)
        OptionA =  MatchesGUI[CVoption][TmpIdxA[InBoth],:].tolist()
        if CVoption == 0:
            OptionA = sorted(OptionA)

        EntryA = ttk.Combobox(EntryFrame, values = OptionA, width = 10)
        EntryA.set(OptionA[0][0])

        OptionB = np.flip(np.argsort(OutputGUI[CVoption][int(EntryA.get()),SessionSwitch[SessionB-1]:SessionSwitch[SessionB]]) + SessionSwitch[SessionB-1])
        OptionB = OptionB.tolist()
    ##################
    EntryB = ttk.Combobox(EntryFrame, values = OptionB, width = 10 )
    EntryB.set(EntryBtmp)
    EntryA.bind('<<ComboboxSelected>>', Update_Units)
    EntryB.bind('<<ComboboxSelected>>', update)

    EntryA.grid(row = 2, column = 1, columnspan = 2, stick = 'WE', padx = 5)
    EntryB.grid(row = 2, column = 3, columnspan = 2, sticky = 'WE', padx = 5)
    MatchIdx = 0

    update(event)

def update_unit_entryB(event):
    """
    This will be called when the user changes the sesion B, and the following will happen:
    - Keep the same A unit
    - update the options in list A
    - update the options for list B
    - Set B to the highest match of unit A for the new session B 
    - Set the Match IDX at 0, for the new list of match pairs
    - update what is in the screen
    """
    global EntryA
    global EntryB
    global OptionA
    global MatchIdx
    global SessionEntryA
    
    SessionA = int(SessionEntryA.get())
    SessionB = int(SessionEntryB.get()) 
    EntryAtmp = int(EntryA.get())

    if EntryA.winfo_exists() == 1:
        EntryA.destroy()
    if EntryB.winfo_exists() == 1:
        EntryB.destroy()

    CV = get_cv_option()
    CVoption = CVtkinter.get() - 1 #want 0 for cv (1,2) , want 1 for cv (2,1)
    if CV == 'Avg':
        #Use average CV values
        TmpIdxA = np.argwhere((SessionSwitch[SessionA -1] <= MatchesAvg[:,0]) * (MatchesAvg[:,0] < SessionSwitch[SessionA]) == True ).squeeze()
        TmpIdxB = np.argwhere((SessionSwitch[SessionB -1] <= MatchesAvg[:,1]) * (MatchesAvg[:,1] < SessionSwitch[SessionB]) == True).squeeze()
        InBoth = np.isin(TmpIdxA, TmpIdxB)
        OptionA = MatchesAvg[TmpIdxA[InBoth],:].tolist()

        EntryA = ttk.Combobox(EntryFrame, values = OptionA, width = 10)
        EntryA.set(EntryAtmp)

        OptionB = np.flip(np.argsort(OutputAvg[int(EntryA.get()),SessionSwitch[SessionB-1]:SessionSwitch[SessionB]]) + SessionSwitch[SessionB-1])
        OptionB = OptionB.tolist()

    else :

        TmpIdxA = np.argwhere((SessionSwitch[SessionA -1] <= MatchesGUI[CVoption][:,0]) * ( MatchesGUI[CVoption][:,0] < SessionSwitch[SessionA]) == True ).squeeze()
        TmpIdxB = np.argwhere((SessionSwitch[SessionB -1] <= MatchesGUI[CVoption][:,1]) * ( MatchesGUI[CVoption][:,1] < SessionSwitch[SessionB]) == True).squeeze()
        InBoth = np.isin(TmpIdxA, TmpIdxB)
        OptionA =  MatchesGUI[CVoption][TmpIdxA[InBoth],:].tolist()
        if CVoption == 0:
            OptionA = sorted(OptionA)

        EntryA = ttk.Combobox(EntryFrame, values = OptionA, width = 10)
        EntryA.set(EntryAtmp)

        OptionB = np.flip(np.argsort(OutputGUI[CVoption][int(EntryAtmp),SessionSwitch[SessionB-1]:SessionSwitch[SessionB]]) + SessionSwitch[SessionB-1])
        OptionB = OptionB.tolist()
    
    ##################
    EntryB = ttk.Combobox(EntryFrame, values = OptionB, width = 10 )
    EntryB.set(OptionB[0])
    EntryA.bind('<<ComboboxSelected>>', Update_Units)
    EntryB.bind('<<ComboboxSelected>>', update)  

    EntryA.grid(row = 2, column = 1, columnspan = 2, stick = 'WE', padx = 5)
    EntryB.grid(row = 2, column = 3, columnspan = 2, sticky = 'WE', padx = 5)
    MatchIdx = 0    

    update(event)

def Update_Units(event):
    """ 
    This function is called when a option in the Unit A dropdown is slected (is a pair of units)
    - split the list of units into two pairs for unit labels for UnitA and UnitB
    - Update the unit B dropdown box to reflect matches for the new Unit A
    - Find the new Match Idx for the new pair which is selected
    """
    global MatchIdx
    global EntryA
    global EntryB
    global OptionA
    global SessionEntryB
    

    SessionB = int(SessionEntryB.get())
    selected = EntryA.get()

    if EntryB.winfo_exists() == 1:
        EntryB.destroy()
    # selected is s string which is xxx yyy , where xxx and yyy are the two units seperates by a space
    tmpA = selected.split()[0]
    tmpB = selected.split()[1]
    EntryA.set(tmpA)
    

    #need to make it so the list for B updates
    CV = get_cv_option()
    CVoption = CVtkinter.get() - 1 #want 0 for cv (1,2) , want 1 for cv (2,1)
    if CV == 'Avg':
        #Use average CV value
        OptionB = np.flip(np.argsort(OutputAvg[int(tmpA),SessionSwitch[SessionB-1]:SessionSwitch[SessionB]]) + SessionSwitch[SessionB-1])
        OptionB = OptionB.tolist()
    else :
        OptionB = np.flip(np.argsort(OutputGUI[CVoption][int(tmpA),SessionSwitch[SessionB-1]:SessionSwitch[SessionB]]) + SessionSwitch[SessionB-1])
        OptionB = OptionB.tolist()
    
    ##################
    EntryB = ttk.Combobox(EntryFrame, values = OptionB, width = 10 )
    EntryB.set(tmpB)
    EntryB.bind('<<ComboboxSelected>>', update)  

    EntryB.grid(row = 2, column = 3, columnspan = 2, sticky = 'WE', padx = 5)


    #need to get list index of the selected postion
    tmpList = [int(tmpA), int(tmpB)]
    MatchIdx = OptionA.index(tmpList)
    update(event)

def next_pair(event):
    """ 
    This function is called when moving to the next unit in the mathc list:
    - get the new unit idx for unit A and unit B
    - update the Units
    - update the dropdown box for unit b to reflect the new Unit A
    """
    global MatchIdx
    global EntryA
    global EntryB
    global OptionA
    global SessionEntryB

    SessionB = int(SessionEntryB.get())
    MatchIdx +=1
    tmpA, tmpB = OptionA[MatchIdx]
    EntryA.set(tmpA)

    if EntryB.winfo_exists() == 1:
        EntryB.destroy()

    #Update B, and change it's dropbox
    CV = get_cv_option()
    CVoption = CVtkinter.get() - 1 #want 0 for cv (1,2) , want 1 for cv (2,1)
    if CV == 'Avg':
        #Use average CV value
        OptionB = np.flip(np.argsort(OutputAvg[int(tmpA),SessionSwitch[SessionB-1]:SessionSwitch[SessionB]]) + SessionSwitch[SessionB-1])
        OptionB = OptionB.tolist()
    else :
        OptionB = np.flip(np.argsort(OutputGUI[CVoption][int(tmpA),SessionSwitch[SessionB-1]:SessionSwitch[SessionB]]) + SessionSwitch[SessionB-1])
        OptionB = OptionB.tolist()
    
    ##################
    EntryB = ttk.Combobox(EntryFrame, values = OptionB, width = 10 )
    EntryB.set(tmpB)
    EntryB.bind('<<ComboboxSelected>>', update)  

    EntryB.grid(row = 2, column = 3, columnspan = 2, sticky = 'WE', padx = 5)

    update(event)

def previous_pair(event):
    """ 
    This function is called when moving to the previous unit in the mathc list:
    - get the new unit idx for unit A and unit B
    - update the Units
    - update the dropdown box for unit b to reflect the new Unit A
    """
    global MatchIdx
    global EntryA
    global EntryB
    global OptionA
    global SessionEntryB
    
    SessionB = int(SessionEntryB.get())
    MatchIdx -=1
    tmpA, tmpB = OptionA[MatchIdx]
    EntryA.set(tmpA)


    if EntryB.winfo_exists() == 1:
        EntryB.destroy()

    #Update B, and change it's dropbox
    CV = get_cv_option()
    CVoption = CVtkinter.get() - 1 #want 0 for cv (1,2) , want 1 for cv (2,1)
    if CV == 'Avg':
        #Use average CV value
        OptionB = np.flip(np.argsort(OutputAvg[int(tmpA),SessionSwitch[SessionB-1]:SessionSwitch[SessionB]]) + SessionSwitch[SessionB-1])
        OptionB = OptionB.tolist()
    else :
        OptionB = np.flip(np.argsort(OutputGUI[CVoption][int(tmpA),SessionSwitch[SessionB-1]:SessionSwitch[SessionB]]) + SessionSwitch[SessionB-1])
        OptionB = OptionB.tolist()
    
    ##################
    EntryB = ttk.Combobox(EntryFrame, values = OptionB, width = 10 )
    EntryB.set(tmpB)
    EntryB.bind('<<ComboboxSelected>>', update)  

    EntryB.grid(row = 2, column = 3, columnspan = 2, sticky = 'WE', padx = 5)

    update(event)

def swap_units():
    global EntryFrame
    global EntryA
    global EntryB
    global SessionEntryB
    global SessionEntryB
    global MatchIdx
    global OptionA

    #get all initial info
    EntryATmp = int(EntryA.get())
    EntryBTmp = int(EntryB.get())
    SessionATmp = int(SessionEntryA.get())
    SesssionBTmp = int(SessionEntryB.get())

    #delte old box a and b
    if EntryA.winfo_exists() == 1:
        EntryA.destroy()
    if EntryB.winfo_exists() == 1:
        EntryB.destroy()

    #swap allpairs
    SessionEntryA.set(SesssionBTmp)
    SessionEntryB.set(SessionATmp)
    
    #Update A and B dropsown boxs
    CV = get_cv_option()
    CVoption = CVtkinter.get() - 1 #want 0 for cv (1,2) , want 1 for cv (2,1)
    if CV == 'Avg':
        #Use average CV values
        TmpIdxA = np.argwhere((SessionSwitch[SesssionBTmp -1] <= MatchesAvg[:,0]) * (MatchesAvg[:,0] < SessionSwitch[SesssionBTmp]) == True ).squeeze()
        TmpIdxB = np.argwhere((SessionSwitch[SessionATmp -1] <= MatchesAvg[:,1]) * (MatchesAvg[:,1] < SessionSwitch[SessionATmp]) == True).squeeze()
        InBoth = np.isin(TmpIdxA, TmpIdxB)
        OptionA = MatchesAvg[TmpIdxA[InBoth],:].tolist()

        EntryA = ttk.Combobox(EntryFrame, values = OptionA, width = 10)
        EntryA.set(EntryBTmp)

        OptionB = np.flip(np.argsort(OutputAvg[int(EntryA.get()),SessionSwitch[SessionATmp-1]:SessionSwitch[SessionATmp]]) + SessionSwitch[SessionATmp-1])
        OptionB = OptionB.tolist()

    else :

        TmpIdxA = np.argwhere((SessionSwitch[SesssionBTmp -1] <= MatchesGUI[CVoption][:,0]) * ( MatchesGUI[CVoption][:,0] < SessionSwitch[SesssionBTmp]) == True).squeeze()
        TmpIdxB = np.argwhere((SessionSwitch[SessionATmp -1] <= MatchesGUI[CVoption][:,1]) * ( MatchesGUI[CVoption][:,1] < SessionSwitch[SessionATmp]) == True).squeeze()
        InBoth = np.isin(TmpIdxA, TmpIdxB)

        #So is orderby unit A
        OptionA =  MatchesGUI[CVoption][TmpIdxA[InBoth],:].tolist()
        if CVoption == 0:
            OptionA = sorted(OptionA)

        EntryA = ttk.Combobox(EntryFrame, values = OptionA, width = 10)
        EntryA.set(EntryBTmp)

        ##NOT SORTED FOR CV
        OptionB = np.flip(np.argsort(OutputGUI[CVoption][int(EntryA.get()),SessionSwitch[SessionATmp-1]:SessionSwitch[SessionATmp]]) + SessionSwitch[SessionATmp-1])
        OptionB = OptionB.tolist()

    ##################
    EntryB = ttk.Combobox(EntryFrame, values = OptionB, width = 10 )
    EntryB.set(EntryATmp)
    
    EntryA.bind('<<ComboboxSelected>>', Update_Units)
    EntryB.bind('<<ComboboxSelected>>', update)
    EntryA.grid(row = 2, column = 1, columnspan = 2, stick = 'WE', padx = 5)
    EntryB.grid(row = 2, column = 3, columnspan = 2, sticky = 'WE', padx = 5)

    tmpList = [int(EntryBTmp), int(EntryATmp)]

    MatchIdx = OptionA.index(tmpList)
    update(None)

def get_score_histrograms(Scores2Include, OutputThreshold):
    """  
    Scores2Include is the dictionary of all scores.
    ProbThreshold is a nUnits*nUnits array where each index is 0 (Not Match) 1 (Match)
    """

    #are lsit of length 6 each list item is 2 np arrays(bins, values)
    HistNames = []
    Hist = []
    HistMatches = []
    for key, values in Scores2Include.items():
        HistNames.append(key)
        Hist.append(np.histogram(values, bins = 100, density = True))
        HistMatches.append(np.histogram(values[OutputThreshold.astype(bool)], bins = 100, density = True))

    return HistNames, Hist, HistMatches    

def add_original_ID(UnitA, UnitB):
    global OriginalIDLabel

    if OriginalIDLabel.winfo_exists():
        OriginalIDLabel.destroy()

    OriginalIDLabel = ttk.Label(root, text = f'The Original Unit IDs are:\nUnit A: {int(ClusInfo["OriginalID"][UnitA].squeeze())}   Unit B: {int(ClusInfo["OriginalID"][UnitB].squeeze())}', borderwidth = 2 , relief= 'groove')
    OriginalIDLabel.grid( row = 0, column = 2, ipadx = 5, ipady = 5)


def add_probability_label(UnitA, UnitB, CVoption):
    global BayesLabel

    if BayesLabel.winfo_exists():
        BayesLabel.destroy()

    if CVoption == -1:
        BayesLabel = ttk.Label(root, text = f'The UM probabilty for this match is:\n {np.round(OutputAvg[UnitA, UnitB],5)}', borderwidth = 2 , relief= 'groove' )
        BayesLabel.grid(row = 0, column = 1, ipadx = 5, ipady = 5)
    else:
        BayesLabel = ttk.Label(root, text = f'The UM probabilty for this match is:\n {np.round(OutputGUI[CVoption][UnitA, UnitB],5)}', borderwidth = 2 , relief= 'groove')
        BayesLabel.grid(row = 0, column = 1, ipadx = 5, ipady = 5)


def set_match(event = None):
    global IsMatch
    UnitA = int(EntryA.get())
    UnitB = int(EntryB.get())

    IsMatch.append( [UnitA, UnitB] )
    IsMatch.append( [UnitB, UnitA] )

    

def set_not_match(event = None):
    global NotMatch
    UnitA = int(EntryA.get())
    UnitB = int(EntryB.get())

    NotMatch.append( [UnitA, UnitB] )
    NotMatch.append( [UnitB, UnitA] )


def MakeTable(table):
    global FrameTable
    
    if FrameTable.winfo_exists() == 1:
        FrameTable.destroy()

    total_rows = len(table)
    total_columns = len(table[0])

    colors = ['black', 'green', 'blue']
    FrameTable = ttk.LabelFrame(root, text = 'UnitData')
    for i in range(total_rows):
        for j in range(total_columns):
                
            e = ttk.Entry(FrameTable, width=20)
            e.insert(END, table[i][j])
            e.configure(state='readonly')             
            e.grid(row=i, column=j)
    FrameTable.grid(row = 3, column = 0, padx = 10, pady = 10)

#get table data #ADD STABILTY - prob of unit with itself accros cv
def get_table_data(UnitA, UnitB, CV):
    template = [['Unit','A','B'],
        ['Avg Centroid','tmp','tmp'],
        ['Amplitude','tmp','tmp'],
        ['Spatial Decay','tmp','tmp'],
        ['Units Matches','tmp','tmp'],
        ['Stability','tmp','tmp']]

    total_rows = len(template)
    total_columns = len(template[0])
    UnitIdxTmp =[0 , UnitA, UnitB]
    table = template

    if CV == 'Avg':
        for i in range(total_rows):
            for j in range(total_columns):
                if j == 0 :
                    continue        
                if i == 0:
                    table[i][j] = str(UnitIdxTmp[j])
                if i == 1:
                    table[i][j] = str( np.round(AvgCentroidAvg[:,UnitIdxTmp[j]], 2))
                if  i == 2:
                    table[i][j] = str( np.round(AmplitudeAvg[UnitIdxTmp[j]], 2))
                if i == 3:
                    table[i][j] = str( np.round(SpatialDecayAvg[UnitIdxTmp[j]], 3))
                if i == 4:
                    table[i][j] = str(np.argwhere( (WithinSession*OutputAvg)[UnitIdxTmp[j],:] > MatchThreshold)).replace('[', '').replace(']', '')
                if i ==5:
                    table[i][j] = str( np.round(OutputGUI[0][UnitIdxTmp[j], UnitIdxTmp[j]], 3))
           
    else:
        for i in range(total_rows):
            for j in range(total_columns):
                if j == 0 :
                    continue        
                if i == 0:
                    table[i][j] = str(UnitIdxTmp[j])
                if i == 1:
                    table[i][j] = str( np.round(AvgCentroid[:,UnitIdxTmp[j],CV[j-1]], 2))
                if  i == 2:
                    table[i][j] = str( np.round(Amplitude[UnitIdxTmp[j],CV[j-1]], 2))
                if i == 3:
                    table[i][j] = str( np.round(SpatialDecay[UnitIdxTmp[j],CV[j-1]], 3))
                if i == 4:
                    table[i][j] = str(np.argwhere( (WithinSession*OutputGUI[CV[0]])[UnitIdxTmp[j],:] > MatchThreshold )).replace('[', '').replace(']', '')
                if i ==5:
                    table[i][j] = str( np.round(OutputGUI[0][UnitIdxTmp[j], UnitIdxTmp[j]], 3))
                    
    return table

def get_unit_score_table(UnitA, UnitB, CVoption):
  
    table = [['tmp'] * 2 for i in range((len(Scores2IncludeAvg) +1))]

    table[0] = ['Score', f'{UnitA} and {UnitB}']

    if CVoption == -1:
        for i in range(len(Scores2IncludeAvg)):
            for j in range(2):
                if j ==0:
                    table[i+1][j] = list(Scores2IncludeAvg.keys())[i]
                else:
                    table[i+1][j] = str( np.round(Scores2IncludeAvg[list(Scores2IncludeAvg.keys())[i]][UnitA, UnitB], 3))
   
    else:
        for i in range(len(Scores2IncludeAvg)):
            for j in range(2):
                if j ==0:
                    table[i+1][j] = list(Scores2IncludeGUI[CVoption].keys())[i]
                else:
                    table[i+1][j] = str( np.round(Scores2IncludeGUI[CVoption][list(Scores2IncludeGUI[CVoption].keys())[i]][UnitA, UnitB], 3))

    return table

def make_unit_score_table(table):
    global ScoreTable
    
    if ScoreTable.winfo_exists() == 1:
        ScoreTable.destroy()

    total_rows = len(table)
    total_columns = len(table[0])

    colors = ['black', 'Purple']
    ScoreTable = ttk.LabelFrame(root, text = 'UM Scores')
    for i in range(total_rows):
        for j in range(total_columns):
                
            e = ttk.Entry(ScoreTable, width = 30)
            e.insert(END, table[i][j])
            e.configure(state='readonly')             
            e.grid(row=i, column=j)

    ScoreTable.grid(row = 3, column = 1, columnspan=2, padx = 10, pady = 10)


def plot_AvgWaveforms(UnitA, UnitB, CV):
    global AvgWaveformPlot
    if AvgWaveformPlot.winfo_exists() == 1:
        AvgWaveformPlot.destroy()

    fig = Figure(figsize = (3, 3), 
                 dpi = 100) 
    fig.patch.set_facecolor('#33393b')
    
    plt1 = fig.add_subplot(111)
    #plt1.spines[["left", "bottom"]].set_position(("data", 0))
    plt1.spines[["bottom"]].set_position(("data", 0))
    plt1.spines[["top", "right"]].set_visible(False)
    #plt1.patch.set_facecolor('none')
    plt1.xaxis.set_label_coords(0.9,0)


    if CV =='Avg':
        plt1.plot(AvgWaveformAvg[:,UnitA], 'g', label=str(UnitA))
        plt1.plot(AvgWaveformAvg[:,UnitB], 'b', label=str(UnitB))
        plt1.set_xlabel('Time')
        plt1.set_ylabel('Amplitude')
        # plt1.set_xlim(left = 0)
        # plt1.set_xticks([])

    else:
        plt1.plot(AvgWaveform[:,UnitA,CV[0]], 'g', label=str(UnitA))
        plt1.plot(AvgWaveform[:,UnitB,CV[1]], 'b', label=str(UnitB))
        plt1.set_xlabel('Time')
        plt1.set_ylabel('Amplitude')
        # plt1.set_xlim(left = 0)

    AvgWaveformPlot=FigureCanvasTkAgg(fig, master = root)
    AvgWaveformPlot.draw()
    AvgWaveformPlot = AvgWaveformPlot.get_tk_widget()
    AvgWaveformPlot.grid(row = 1, column = 0)

def plot_trajectories(UnitA, UnitB, CV):
    global TrajectoryPlot
    if TrajectoryPlot.winfo_exists() == 1:
        TrajectoryPlot.destroy()

    fig = Figure(figsize = (3,3), 
                 dpi = 100) 
    fig.patch.set_facecolor('#33393b')

    
    plt2 = fig.add_subplot(111)
    #plt2.patch.set_facecolor('#33393b')
    plt2.set_aspect(0.5)
    plt2.spines[['right', 'top']].set_visible(False)
    

    if CV =='Avg':
        
        # AM not doing a time averaged WaveIDX (where you fins goodtimepoints), will just uses CV 0 for both
        plt2.plot(AvgWaveformPerTPAvg[1,UnitA,WaveIdx[UnitA,:,0].astype(bool)], AvgWaveformPerTPAvg[2,UnitA,WaveIdx[UnitA,:,0].astype(bool)], 'g')
        plt2.scatter(AvgCentroidAvg[1,UnitA], AvgCentroidAvg[2,UnitA], c = 'g')

        plt2.plot(AvgWaveformPerTPAvg[1,UnitB,WaveIdx[UnitB,:,0].astype(bool)], AvgWaveformPerTPAvg[2,UnitB,WaveIdx[UnitB,:,0].astype(bool)],'b')
        plt2.scatter(AvgCentroidAvg[1,UnitB], AvgCentroidAvg[2,UnitB], c = 'b')

        plt2.set_xlabel(r'Xpos ($\mu$m)')
        plt2.set_ylabel(r'Ypos ($\mu$m)')

    else:
        plt2.plot(AvgWaveformPerTP[1,UnitA,WaveIdx[UnitA,:,CV[0]].astype(bool), CV[0]], AvgWaveformPerTP[2,UnitA,WaveIdx[UnitA,:,CV[0]].astype(bool), CV[0]], 'g')
        plt2.scatter(AvgCentroid[1,UnitA,CV[0]], AvgCentroid[2,UnitA,CV[0]], c = 'g')

        plt2.plot(AvgWaveformPerTP[1,UnitB,WaveIdx[UnitB,:,CV[1]].astype(bool), CV[1]], AvgWaveformPerTP[2,UnitB,WaveIdx[UnitB,:,CV[1]].astype(bool), CV[1]],'b')
        plt2.scatter(AvgCentroid[1,UnitB,CV[1]], AvgCentroid[2,UnitB,CV[1]], c = 'b')

        plt2.set_xlabel(r'Xpos ($\mu$m)')
        plt2.set_ylabel(r'Ypos ($\mu$m)')

    TrajectoryPlot=FigureCanvasTkAgg(fig, master = root)
    TrajectoryPlot.draw()
    TrajectoryPlot = TrajectoryPlot.get_tk_widget()
#    TrajectoryPlot.configure(bg = '#33393b')
    TrajectoryPlot.grid(row = 1, column = 1, columnspan = 2)

def order_good_sites(GoodSites, ChannelPos, SessionNo):
    # make it so y goes from biggest to smallest
    ReorderedIdx = np.argsort(-ChannelPos[SessionNo][GoodSites,2].squeeze())
    ReOrderedGoodSites = GoodSites[ReorderedIdx]

    #re-arange x-axis so it goes (smaller x, bigger x)
    for i in range(9):
        a,b = ChannelPos[SessionNo][ReOrderedGoodSites[[2*i, 2*i+1]],1]

        if a > b:
            #swap order
            ReOrderedGoodSites[[2*i + 1, 2*i]] = ReOrderedGoodSites[[2*i, 2*i+1]]
    return ReOrderedGoodSites

def nearest_channels(MaxSite, MaxSiteMean, ChannelPos, ClusInfo, Unit, CV):

    SessionNo = ClusInfo['SessionID'][Unit]
    if CV == 'Avg':
        maxsite = MaxSiteMean[Unit].squeeze()
        __, x, y = ChannelPos[SessionNo][maxsite,:]

        GoodXsites = np.argwhere( np.logical_and((x-50 < ChannelPos[SessionNo][:,1]) == True, (ChannelPos[SessionNo][:,1] < x+50) == True))
        yValues = ChannelPos[SessionNo][GoodXsites,2]

        yDist2MaxSite = np.abs(yValues - ChannelPos[SessionNo][maxsite,2])
        GoodSites = np.argsort(yDist2MaxSite,axis = 0 )[:18]
        GoodSites = GoodXsites[GoodSites]
        ReOrderedGoodSites = order_good_sites(GoodSites, ChannelPos, SessionNo)

    else:
        maxsite = MaxSite[Unit,CV]
        __, x, y = ChannelPos[SessionNo][maxsite,:]

        GoodXsites = np.argwhere( np.logical_and((x-50 < ChannelPos[SessionNo][:,1]) == True, (ChannelPos[SessionNo][:,1] < x+50) == True))
        yValues = ChannelPos[SessionNo][GoodXsites,2]

        yDist2MaxSite = np.abs(yValues - ChannelPos[SessionNo][maxsite,2])
        GoodSites = np.argsort(yDist2MaxSite,axis = 0 )[:18]
        GoodSites = GoodXsites[GoodSites]
        ReOrderedGoodSites = order_good_sites(GoodSites, ChannelPos, SessionNo)
        
    return ReOrderedGoodSites

def plot_raw_waveforms(UnitA, UnitB, CV):

    SessionNoA = ClusInfo['SessionID'][UnitA]
    global RawWaveformPlot
    if RawWaveformPlot.winfo_exists() == 1:
        RawWaveformPlot.destroy()

    fig = Figure(figsize=(4,6), dpi = 100)
    fig.set_tight_layout(False)
    fig.patch.set_facecolor('#33393b')

    MainAx = fig.add_axes([0.2,0.2,0.8,0.8])
    #MainAx.set_facecolor('none')
    MainAxOffset = 0.2
    MainAxScale = 0.8

    if CV =='Avg':
        GoodChannels = nearest_channels(MaxSite, MaxSiteMean, ChannelPos, ClusInfo, UnitA, CV)

        minx, miny = ChannelPos[SessionNoA][GoodChannels[-2],[1,2]].squeeze()
        maxx, maxy = ChannelPos[SessionNoA][GoodChannels[1],[1,2]].squeeze()
        deltaX = (maxx - minx) / 2
        deltaY = (maxy - miny) / 18

        #may want to change so it find this for both units and selects the most extreme arguments
        #however i dont think tis will be necessary
        SubMiny = np.nanmin(waveform[UnitA,:,GoodChannels].mean(axis=-1))
        SubMaxy = np.nanmax(waveform[UnitA,:,GoodChannels].mean(axis=-1))
        # shift each waveform so 0 is at the channel site, 1/9 is width of a y waveform plot
        WaveformYoffset = (np.abs(SubMaxy) / (np.abs(SubMiny) + np.abs(SubMaxy)) ) * 1/9
    
    else:
        GoodChannels = nearest_channels(MaxSite, MaxSiteMean, ChannelPos, ClusInfo, UnitA, CV[0])

        minx, miny = ChannelPos[SessionNoA][GoodChannels[-2],[1,2]].squeeze()
        maxx, maxy = ChannelPos[SessionNoA][GoodChannels[1],[1,2]].squeeze()
        deltaX = (maxx - minx) / 2
        deltaY = (maxy - miny) / 18
        #may want to change so it find this for both units and selects the most extreme arguments
        #however i dont think this will be necessary
        SubMiny = np.nanmin(waveform[UnitA,:,GoodChannels, CV[0]])
        SubMaxy = np.nanmax(waveform[UnitA,:,GoodChannels, CV[0]])
        # shift each waveform so 0 is at the channel site, 1/9 is width of a y waveform plot
        WaveformYoffset = (np.abs(SubMaxy) / (np.abs(SubMiny) + np.abs(SubMaxy)) ) * 1/9


    #make the main scatter positiose site as scatter with opacity 
    MainAx.scatter(ChannelPos[SessionNoA][GoodChannels,1], ChannelPos[SessionNoA][GoodChannels,2], c = 'grey', alpha = 0.3)
    MainAx.set_xlim(minx - deltaX, maxx + deltaX)
    MainAx.set_ylim(miny - deltaY, maxy + deltaY)

    for i in range(9):
        for j in range(2):
            #may need to change this positioning if units sizes are irregular
            if j == 0:
                #The peak in the waveform is not half way, so maths says the x axis should be starting at
                #0.1 and 0.6 so the middle is at 0.25/0.76 however chosen these values so it loks better by eye

                #
                ax =  fig.add_axes([MainAxOffset + MainAxScale*0.25, MainAxOffset + MainAxScale*(i/9 - 1/18 + WaveformYoffset), MainAxScale*0.25, MainAxScale*1/9])
            if j == 1:
                ax = fig.add_axes([MainAxOffset + MainAxScale*0.75, MainAxOffset + MainAxScale*(i/9 - 1/18 + WaveformYoffset), MainAxScale*0.25, MainAxScale*1/9])

            if CV =='Avg':
                ax.plot(waveform[UnitA,:,GoodChannels[i*2 + j]].mean(axis=-1).squeeze(), color = 'g')             
                ax.plot(waveform[UnitB,:,GoodChannels[i*2 + j]].mean(axis = -1).squeeze(), color = 'b', lw=0.8)
            else:
                ax.plot(waveform[UnitA,:,GoodChannels[i*2 + j], CV[0]].squeeze(), color = 'g')             
                ax.plot(waveform[UnitB,:,GoodChannels[i*2 + j],CV[1]].squeeze(), color = 'b', lw=0.8)                
            ax.set_ylim(SubMiny,SubMaxy)
            ax.set_axis_off()


    MainAx.spines.right.set_visible(False)
    MainAx.spines.top.set_visible(False)
    MainAx.set_xticks([minx, maxx])
    MainAx.set_xlabel('Xpos ($\mu$m)', size = 14)
    MainAx.set_ylabel('Ypos ($\mu$m)', size = 14)


    RawWaveformPlot = FigureCanvasTkAgg(fig, master = root)
    RawWaveformPlot.draw()
    RawWaveformPlot = RawWaveformPlot.get_tk_widget()
    #RawWaveformPlot.configure(bg = '#33393b')

    RawWaveformPlot.grid(row = 0, column = 4, rowspan = 4, padx = 15, pady = 25, ipadx = 15)

def plot_Histograms(HistNames, Hist, HistMatched, Scores2Include, UnitA, UnitB):

    global HistPlot
    if HistPlot.winfo_exists() == 1:
        HistPlot.destroy()

    fig = Figure(figsize= (4,6), dpi = 100, layout = 'constrained')
    fig.patch.set_facecolor('#33393b')
    axs = fig.subplots(3,2, sharex='col')
    axs = axs.flat


    #loop over indexes..
    for i in range(len(Hist)):
        axs[i].step(Hist[i][1][:-1], Hist[i][0], color = 'orange')
        axs[i].step(HistMatched[i][1][:-1], HistMatched[i][0], color = 'magenta')
        axs[i].set_ylim(bottom=0)
        axs[i].set_title(HistNames[i], fontsize = 12)
        #axs[i].get_yaxis().set_visible(False)
        axs[i].axvline(Scores2Include[HistNames[i]][UnitA,UnitB], ls = '--', color = 'grey')


    HistPlot = FigureCanvasTkAgg(fig, master = root)
    HistPlot.draw()
    HistPlot = HistPlot.get_tk_widget()

    HistPlot.grid(row = 0, column = 3, rowspan = 4, padx = 5, pady = 20)

