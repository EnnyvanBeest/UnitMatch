function UMparam = DefaultParametersUnitMatch(UMparam)
% Checks if all required parameters are there, and will fill in default
% parameters if not wherever possible
%% Check if path information is given
if ~isfield(UMparam,'SaveDir')
    disp('Warning, no SaveDir given. Assigning current directory')
    UMparam.SaveDir = fullfile(cd,'UnitMatch');
end
if ~isfield(UMparam,'KSDir')
    error('Error, no directories provided with sorted data and/or raw waveforms')
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
if ~isfield(UMparam,'sampleamount')
    UMparam.sampleamount = 1000; % n raw waveforms to extract
end
if ~isfield(UMparam,'spikeWidth')
    UMparam.spikeWidth = 82; % width of spikes in samples (typically assuming 30KhZ sampling)
end
if ~isfield(UMparam,'RedoExtraction')
    UMparam.RedoExtraction = 0; % Redoing raw average spike extraction --> this is time consuming
end

%% Parameters used in standard UnitMatch
if ~isfield(UMparam,'RunPyKSChronicStitched')
    UMparam.RunPyKSChronicStitched = 0;
end

if ~isfield(UMparam,'UseHistology')
    UMparam.UseHistology = 0; % You can use actual Allen coordinates, if avaible and if using the example pipeline
end
if ~isfield(UMparam,'ProbabilityThreshold')
    UMparam.ProbabilityThreshold = 0.5; % Threshold for assigning as a match
end
if ~isfield(UMparam,'Scores2Include')
    UMparam.Scores2Include = {'CentroidDist','WavformSim','CentroidOverlord','spatialdecaySim','AmplitudeSim','LocTrajectorySim'};
end
if ~isfield(UMparam,'ApplyExistingBayesModel')
    UMparam.ApplyExistingBayesModel = 0; % Just apply an already existing Bayes Model (not recommended)
end
if ~isfield(UMparam,'GoodUnitsOnly')
    UMparam.GoodUnitsOnly = 1; % Recommended, to only use units that are good single units (e.g. use Bombcell)
end
if ~isfield(UMparam,'minGoodUnits')
    UMparam.minGoodUnits = 25; % Recommended
end

%% Post processing
if ~isfield(UMparam,'AssignUniqueID')
    UMparam.AssignUniqueID = 1; % Use our method of Assigning Unique ID based on match probability (recommended)
end
if ~isfield(UMparam,'UseDatadrivenProbThrs')
    UMparam.UseDatadrivenProbThrs = 0; % Use data driven probability threshold, if 0 use UMparam.ProbabilityThreshold
end
if ~isfield(UMparam,'min_angledist')
    UMparam.min_angledist = 0.1 % the minimum distance a centroid must move to take it angle, avoid getting angle of noise
end    
%% Inspection
if ~isfield(UMparam,'MakePlotsOfPairs')
    UMparam.MakePlotsOfPairs = 0; % Plots and saves matches for you to inspect
end
if ~isfield(UMparam,'GUI')
    UMparam.GUI = 0; % Use our GUI to flip through and curate our matches
end
if ~isfield(UMparam,'RunPyKSChronicStitched')
    UMparam.RunPyKSChronicStitched = 0; % whether stitched KS was used (concatenated recordings) (it will compare KS labels to UM labels)
end
%% Parameters only used for compute functional scores
if ~isfield(UMparam,'ACGbinSize')
    UMparam.ACGbinSize = 1E-03; %
end
if ~isfield(UMparam,'ACGduration')
    UMparam.ACGduration = 1; % in seconds
end
if ~isfield(UMparam,'binsize')
    UMparam.binsize = 0.01; % in seconds
end

return
