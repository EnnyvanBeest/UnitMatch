function [PrepareClusInfoparams, UMparam, UniqueIDConversion, MatchTable, WaveformInfo] = um_runUnitMatch(kilosort_dirs, ephys_dirs, SaveDir, tmpdatafolder, thisRecordingType, mouseName)
% Run unit match
% ------
% Inputs
% ------
% kilosort_dirs: cell array, with each cell containing the path to your raw .bin, .cbin or .dat data.
% ephys_dirs: cell array, with each cell containing the path to your kilosort output files.
% SaveDir: string, where to save unitmatch data
% tmpdatafolder: string, where to store any temporary files (eg decompressed ephys files, ...)
% thisRecordingType: string. either 'Chronic' or 'Acute'
% mouseName: string
% ------
% Outputs
% ------
% PrepareClusInfoparams: strcuture containing the paramaters used
% UMparam: strcuture containing the paramaters used
% UniqueIDConversion: structure containing a summary unit match's results:
%    - UniqueIDConversion.UniqueID returns, for each unit, the ID assigned
%    to them. Any duplicate values indicates units that were matched
%    - UniqueIDConversion.OriginalClusID is a 1 x number_of_units unit32 vector, gives each unit's original label
%    - UniqueIDConversion.recsesAll is a 1 x number_of_units double vector, returns from which recording (1:number_of_recordings)
%    each unit comes from
%    - UniqueIDConversion.GoodID is a 1 x number_of_units binary double
%    vector, returns for each unit whether it was classified as good or
%    not.
% MatchTable: structure detailing unit match's results, with scores for
% each metric used:
% ------

%% Parameters and settings
PrepareClusInfoparams.RunPyKSChronicStitched = 0; % Default 0. if 1, run PyKS chronic recordings stitched when same IMRO table was used
PrepareClusInfoparams.CopyToTmpFirst = 1; % If 1, copy data to local first, don't run from server (= advised!)
PrepareClusInfoparams.DecompressLocal = 1; % If 1, uncompress data first if it's currently compressed (= necessary for unitmatch and faster for QualityMetrics)
PrepareClusInfoparams.deNoise = 0;
PrepareClusInfoparams.nSavedChans = 385;
PrepareClusInfoparams.nSyncChans = 1;

% Storing preprocessed data?
PrepareClusInfoparams.ReLoadAlways = 1; % If 1, SP & Clusinfo are always loaded from KS output
PrepareClusInfoparams.saveSp = 1; % Save SP struct for easy loading of preprocessed data
PrepareClusInfoparams.binsz = 0.01; %Bin size for PSTHs in seconds
PrepareClusInfoparams.deNoise = 1; %whether to try and denoise data

% Quality Metrics
PrepareClusInfoparams.RunQualityMetrics = 1; % If 1, Run the quality metrics (Bombcell @JulieFabre)
PrepareClusInfoparams.RedoQM = 0; %if 1, redo quality metrics if it already exists
PrepareClusInfoparams.InspectQualityMetrics = 0; % If 1, Inspect the quality matrix/data set using the GUI (manual inspection)
PrepareClusInfoparams.loadPCs = 0; % Only necessary when computiong isoluation metrics/drift in QM. You save a lot of time keeping this at 0

% UnitMatch
PrepareClusInfoparams.UnitMatch = 1; % If 1, find identical units across sessions or oversplits in a fast and flexible way
PrepareClusInfoparams.RedoUnitMatch = 1; % if 1, Redo unitmatch
PrepareClusInfoparams.separateIMRO = 0; % Run for every IMRO separately (for memory reasons or when having multiple probes this might be a good idea)

% UnitMatch Parameters:
% All parameters to choose from: {'AmplitudeSim','spatialdecaySim','WavformMSE','WVCorr','CentroidDist','CentroidVar','CentroidDistRecentered','TrajAngleSim','TrajDistSim'};
% WavformSim is average of WVCorr and WavformMSE
% CentroidOverlord is average of CentroidDistRecentered and CentroidVar
PrepareClusInfoparams.Scores2Include = {'CentroidDist', 'WavformSim', 'CentroidOverlord', 'spatialdecaySim', 'AmplitudeSim', 'TrajAngleSim'}; %{'AmplitudeSim','spatialdecayfitSim','WavformSim','CentroidDist','CentroidVar','TrajAngleSim'}; %
PrepareClusInfoparams.ApplyExistingBayesModel = 0; %If 1, use probability distributions made available by us -
PrepareClusInfoparams.MakePlotsOfPairs = 0; % Plots pairs for inspection (UnitMatch)
PrepareClusInfoparams.AssignUniqueID = 1; % Assign UniqueID
PrepareClusInfoparams.GoodUnitsOnly = 1; % Include only good untis in the UnitMatch analysis - faster and more sensical
PrepareClusInfoparams.extractSync = 0;
PrepareClusInfoparams.plot = 0;
PrepareClusInfoparams.UseHistology =0;

PrepareClusInfoparams.SaveDir = SaveDir; % Save results here
PrepareClusInfoparams.tmpdatafolder = tmpdatafolder; % use this as a local directory (should be large enough to handle all sessions you want to combine)

PrepareClusInfoparams.loadMATsToSave = 1;
%% Loading data from kilosort/phy
subsesopt = kilosort_dirs;

if strcmp(thisRecordingType, 'Chronic')
    if ~PrepareClusInfoparams.RunPyKSChronicStitched %MatchUnitsAcrossDays
        disp('Unit matching in Matlab')
        subsesopt(cell2mat(cellfun(@(X) any(strfind(X, 'Chronic')), subsesopt, 'UniformOutput', 0))) = []; %Use separate days and match units via matlab script
    else
        disp('Using chronic pyks option')
        subsesopt = subsesopt(cell2mat(cellfun(@(X) any(strfind(X, 'Chronic')), subsesopt, 'UniformOutput', 0))); %Use chronic output from pyks
    end
end


%Create saving directories
clear params
thisIMRO = '';
thisdate = [];
if ~exist(fullfile(SaveDir))
    mkdir(fullfile(SaveDir))
end
PrepareClusInfoparams.SaveDir = fullfile(SaveDir);
if isempty(subsesopt)
    disp(['No data found for ', mouseName])

end
saveJF = 1;

%% Pre-process: quality metrics
PrepareClusInfoparams = PrepareClusInfo(subsesopt, PrepareClusInfoparams, ephys_dirs);

%% Run unit match
[UMparam, UniqueIDConversion, MatchTable, WaveformInfo] = RunUnitMatch(subsesopt, PrepareClusInfoparams, ephys_dirs);

%% Output plots
if PrepareClusInfoparams.plot
    % Evaluate (within unit ID cross-validation)
    EvaluatingUnitMatch(UMparam.SaveDir);

    % Function analysis
    ComputeFunctionalScores(UMparam.SaveDir, saveJF)
end
% Figures
if UMparam.MakePlotsOfPairs
    DrawBlind = 0; %1 for blind drawing (for manual judging of pairs)
    DrawPairsUnitMatch(UMparam.SaveDir, DrawBlind, saveJF);
end

% Quality metrics
try
    QualityMetricsROCs(UMparam.SaveDir);
catch ME
    disp(['Couldn''t do Quality metrics for ', mouseName])
end

end