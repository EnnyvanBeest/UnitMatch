%% Example pipeline going from raw Neuropixels data recorded with SpikeGLX to curated (with Bombcell) units identified with UnitMatch

%% User Input
%% Path information
DataDir = {'H:\OpenEphys_Example\RawData'}; %{'H:\MatchingUnits\RawDataMonthApart'};% ;%Raw data folders, typically servers were e.g. *.cbin files are stored
SaveDir = 'H:\OpenEphys_Example\SaveDir'; %'\\znas.cortexlab.net\Lab\Share\UNITMATCHTABLES_ENNY_CELIAN_JULIE\2ConsecutiveDays\Stitched';%'\\znas.cortexlab.net\Lab\Share\UNITMATCHTABLES_ENNY_CELIAN_JULIE\MonthApart\Stitched';%%'H:\MatchingUnits\Output\MonthApartStitched'% 'H:\MatchingUnits\Output\NotConcatenated';%'\\znas.cortexlab.net\Lab\Share\Celian\UnitMatch\MatchTables\NewSep27\MonthApart\Stitched'% %%;% %'H:\MatchingUnits\Output\ManyRecordings'%Folder where to store the results
tmpdatafolder = 'H:\OpenEphys_Example\Tmp'; % temporary folder for temporary decompression of data 
KilosortDir = 'H:\OpenEphys_Example\KilosortOutput'; % '\\znas.cortexlab.net\Lab\Share\Enny\UnitMatch\KSComparisonSubset';%'\\znas.cortexlab.net\Lab\Share\Enny\UnitMatch\KilosortOutputMonthApart';%'H:\MatchingUnits\KilosortOutputMonthApart';%'\\znas.cortexlab.net\Lab\Share\Celian\UnitMatch\KilosortOutputMonthApart';% Kilosort output folder
GithubDir = 'C:\Users\EnnyB\Documents\GitHub'; % Github directory
PythonEXE = 'C:\Users\EnnyB\anaconda3\envs\pyks2\pythonw.exe' % Python version to run python code in:

%% Information on experiments
MiceOpt = {'AL032','AV008','CB016','EB019','JF067'}; %'AL032', Add all mice you want to analyze
DataDir2Use = repmat(1,[1,length(MiceOpt)]); % In case you have multiple DataDir, index which directory is used for each mouse
RecordingType = repmat({'Chronic'},1,length(MiceOpt)); % And whether recordings were Chronic (default)
RecordingType(ismember(MiceOpt,{''}))={'Acute'}; %EB014', % Or maybe acute?

%% Parameters on how to prepare units/data for analysis
PipelineParams.RunPyKSChronicStitched = 1; % Default 0. if 1, run PyKS chronic recordings stitched when same IMRO table was used
PipelineParams.CopyToTmpFirst = 1; % If 1, copy data to local first, don't run from server (= advised!)
PipelineParams.DecompressLocal = 1; % If 1, uncompress data first if it's currently compressed (= necessary for unitmatch and faster for QualityMetrics)

% Storing preprocessed data?
PipelineParams.ReLoadAlways = 0; % If 1, SP & Clusinfo are always loaded from KS output
PipelineParams.saveSp = 1; % Save SP struct for easy loading of preprocessed data
PipelineParams.binsz = 0.01; %Bin size for PSTHs in seconds

% Quality Metrics
PipelineParams.RunQualityMetrics = 1; % If 1, Run the quality metrics (Bombcell @JulieFabre)
PipelineParams.RedoQM = 0; %if 1, redo quality metrics if it already exists
PipelineParams.InspectQualityMetrics = 0; % If 1, Inspect the quality matrix/data set using the GUI (manual inspection)
PipelineParams.loadPCs = 0; % Only necessary when computiong isoluation metrics/drift in QM. You save a lot of time keeping this at 0

% UnitMatch
PipelineParams.UnitMatch = 1; % If 1, find identical units across sessions or oversplits in a fast and flexible way
PipelineParams.RedoUnitMatch = 1; % if 1, Redo unitmatch
PipelineParams.separateIMRO = 0; % Run for every IMRO separately (for memory reasons or when having multiple probes this might be a good idea)
PipelineParams.UseHistology = 0; % Use real coordinates (3D space of tracked probes if available)

% UnitMatch Parameters:
% All parameters to choose from: {'AmplitudeSim','spatialdecaySim','WavformMSE','WVCorr','CentroidDist','CentroidVar','CentroidDistRecentered','TrajAngleSim','TrajDistSim','spatialdecayfitSim'};
% WavformSim is average of WVCorr and WavformMSE
% CentroidOverlord is average of CentroidDistRecentered and CentroidVar
% LocTrajectorySim is average of TrajAngleSim and TrajDistSim
PipelineParams.Scores2Include = {'CentroidDist','WavformSim','CentroidOverlord','spatialdecaySim','AmplitudeSim','LocTrajectorySim'}; %{'AmplitudeSim','spatialdecayfitSim','WavformSim','CentroidDist','CentroidVar','TrajAngleSim'}; % 
PipelineParams.ApplyExistingBayesModel = 0; %If 1, use probability distributions made available by us - 
PipelineParams.AssignUniqueID = 1; % Assign UniqueID 
PipelineParams.GoodUnitsOnly = 1; % Include only good untis in the UnitMatch analysis - faster and more sensical
PipelineParams.MakePlotsOfPairs = 0; % Plots pairs for inspection (UnitMatch)
PipelineParams.GUI = 1; % Flick through and do manual curation of matching - only works if MakePlotsofPairs = 1

%% Automatic from here
PipelineParams.SaveDir = SaveDir; % Save results here
PipelineParams.tmpdatafolder = tmpdatafolder; % use this as a local directory (should be large enough to handle all sessions you want to combine)

%% All dependencies you want to add (you may need to download these, all available via github)
addpath(genpath(cd))

% Required (for using UnitMatch):
addpath(genpath(fullfile(GithubDir,'spikes')))% Should work with normal spikes toolbox, but I use the forked version in https://github.com/EnnyvanBeest/spikes
addpath(genpath(fullfile(GithubDir,'npy-matlab'))) % https://github.com/kwikteam/npy-matlab
addpath(genpath(fullfile(GithubDir,'mtscomp'))) % https://github.com/int-brain-lab/mtscomp

% Advised (quality metrics for unit selection):
addpath(genpath(fullfile(GithubDir,'bombcell'))) % DOI: 10.5281/zenodo.8172822, https://github.com/Julie-Fabre/bombcell 

% Optional for histology:
addpath(genpath(fullfile(GithubDir,'AP_histology'))) % https://github.com/petersaj/AP_histology
addpath(genpath(fullfile(GithubDir,'allenCCF'))) % https://github.com/cortex-lab/allenCCF

% UNITMATCH - Move to top of paths 
addpath(genpath(fullfile(GithubDir,'UnitMatch'))) % Make sure to have this one fresh in the path (so run this last)

try
    % Only need to do this once:
    % Follow instructions on installing pykilosort in anaconda environment,
    % eg. https://github.com/int-brain-lab/pykilosort
    % Additially run (from within this environment):
    % pip install matlab
    % Python version to run python code in:    pyversion(PythonEXE) %Explanation on how to do this is provided in the README
catch ME
    disp(ME)
end

%% Actual pipeline
%% PyKS - run pykilosort from Matlab/Python integration
RunPyKS2_FromMatlab

%% Runs unitmatch across all data from a mouse to generate a table
RunUnitMatchAllDataPerMouse

%% Across Mice Graphs
% SummarizeAcrossMice

FromDate = datetime("2023-10-02 13:00:00");
UMFiles = cell(1,0); % Define your UMfolders here or use below:
groupvec = nan(1,0);
if ~exist('UMFiles') || isempty(UMFiles) % When using the example pipeline this may be useful:
    for midx = 1:length(MiceOpt)
        fprintf('Reference %s...\n', MiceOpt{midx})
        % Identify all UM tables
        tmpfile = dir(fullfile(SaveDir, MiceOpt{midx},'*','*','UnitMatch', 'UnitMatch.mat'));
        if isempty(tmpfile) 
            continue
        end
        for id = 1:length(tmpfile)
            if datetime(tmpfile(id).date) >FromDate 
%                AssignUniqueID(fullfile(tmpfile(id).folder,tmpfile(id).name));

                %             FolderParts = strsplit(tmpfile(id).folder,filesep);
                %             idx = find(ismember(FolderParts,MiceOpt{midx}));
                UMFiles = cat(2,UMFiles,fullfile(tmpfile(id).folder,tmpfile(id).name));
                groupvec = cat(2,groupvec,midx);
            end
        end
    end
end
summaryFunctionalPlots(UMFiles, 1, groupvec)