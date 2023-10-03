function UMparam = DefaultParametersUnitMatch(UMparam)

%% Check if path information is given
if ~isfield(UMparam,'SaveDir')
    disp('Warning, no SaveDir given. Assigning current directory')
    UMparam.SaveDir = fullfile(cd,'UnitMatch');
end
if ~isfield(UMparam,'KSDir')
    error('Warning, no directories provided with sorted data and/or raw waveforms')
end

if ~isfield(UMparam,'RawDataPaths')
    warning('No RawDataPaths given... if raw waveforms are already extracted this may not be a problem')
    UMparam.RawDataPaths = [];
end
% if not given, assume decompPath (decompressed path) is same as raw
if ~isfield(UMparam,'AllDecompPaths')
    warning('No Decompressed data paths given, assuming RawDataPaths point at .bin files')
    UMparam.AllDecompPaths = UMparam.RawDataPaths;
end

if ~isfield(UMparam,'AllChannelPos')
    error('ChannelPositions not given. TIP: Use UMparam = ExtractKilosortData(KiloSortPaths, UMparam) to automatically extract these')
end

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

return
