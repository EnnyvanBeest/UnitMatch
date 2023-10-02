function UMparam = DefaultParametersUnitMatch(SaveDir,KSDir,channelpos,RawDataPaths)

%% Parameters for extracting raw waveforms
UMparam.sampleamount = 1000; % n raw waveforms to extract
UMparam.spikeWidth = 82; % width of spikes in samples (typically assuming 30KhZ sampling)
UMparam.RedoExtraction = 0; % Redoing raw average spike extraction --> this is time consuming

%% Parameters used in standard UnitMatch
UMparam.UseHistology = 0; % You can use actual Allen coordinates, if avaible and if using the example pipeline
UMparam.ProbabilityThreshold = 0.5; % Threshold for assigning as a match
UMparam.Scores2Include = {'CentroidDist','WavformSim','CentroidOverlord','spatialdecaySim','AmplitudeSim','LocTrajectorySim'};
UMparam.ApplyExistingBayesModel = 0; % Just apply an already existing Bayes Model (not recommended)
UMparam.GoodUnitsOnly = 1; % Recommended, to only use units that are good single units (e.g. use Bombcell)

%% Post processing
UMparam.AssignUniqueID = 1; % Use our method of Assigning Unique ID based on match probability (recommended)

%% Inspection
UMparam.MakePlotsOfPairs = 0; % Plots and saves matches for you to inspect
UMparam.GUI = 0; % Use our GUI to flip through and curate our matches
UMparam.RunPyKSChronicStitched = 0; % whether stitched KS was used (concatenated recordings) (it will compare KS labels to UM labels)

%% Parameters only used for compute functional scores
UMparam.ACGbinSize = 1E-03; % 
UMparam.ACGduration = 1; % in seconds
UMparam.binsize = 0.01; % in seconds

%% Path information
UMparam.SaveDir = SaveDir;
UMparam.KSDir = KSDir; % Cell array with path in every cell to KS output
UMparam.channelpos = channelpos;
UMparam.AllDecompPaths = RawDataPaths; % Assuming these are .bin files. You can use the general pipeline if you need to compress data
UMparam.RawDataPaths = RawDataPaths;  % Could be .cbin or .bin files
