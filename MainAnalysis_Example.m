%% User Input
%% Path information
DataDir = {'H:\MatchingUnits\RawData'};% Check DataDir2Use
SaveDir = 'H:\MatchingUnits\Output\Concatenated1Day\'
tmpdatafolder = 'H:\MatchingUnits\Tmp'; % temporary folder 
KilosortDir = 'H:\MatchingUnits\KilosortOutput';% 'E:\Data\KiloSortOutput';%

%% Information on experiments
MiceOpt = {'AL032','AV008','CB016','EB019','JF067'};%{'AV009','AV015','EB014','EB019','AL032','AV008','JF067','CB016'}%,'AV008','JF067','CB016''EB019'}; %CB016 %AL032 'AV008' JF067 Add all mice you want to analyse
% nidq_sync_used = zeros(1,length(MiceOpt)); % Was an external nidq used for syncing (typically sync feeds directly into IMEC)
% nidq_sync_used(ismember(MiceOpt,{'EB001','CB007','CB008'}))=1; % Except for these mice...
DataDir2Use = repmat(1,[1,length(MiceOpt)]); % In case you have multiple DataDir, index which directory is used for each mouse
RecordingType = repmat({'Acute'},1,length(MiceOpt)); % And whether recordings were acute (default)
RecordingType(ismember(MiceOpt,{'AL032','EB019','CB016','AV008','JF067'}))={'Chronic'}; %EB014', % Or maybe Chronic?

%% Parameters on how to prepare units/data for analysis
PrepareClusInfoparams.RunPyKSChronicStitched = 0; % Default 0. if 1, run PyKS chronic recordings stitched when same IMRO table was used
PrepareClusInfoparams.CopyToTmpFirst = 1; % If 1, copy data to local first, don't run from server (= advised!)
PrepareClusInfoparams.DecompressLocal = 1; % If 1, uncompress data first if it's currently compressed (= necessary for unitmatch and faster for QualityMetrics))
PrepareClusInfoparams.RedoQM = 0; %if 1, redo quality metrics if it already exists
PrepareClusInfoparams.RunQualityMetrics = 1; % If 1, Run the quality metrics (Bombcell @JulieFabre)
PrepareClusInfoparams.loadPCs=1; % Do we need the PCs for data analysis (YES IF QM!)? If not you save a lot of time keeping this at 0
PrepareClusInfoparams.InspectQualityMetrics=0; % If 1, Inspect the quality matrix/data set using the GUI (manual inspection)
PrepareClusInfoparams.UnitMatch = 1; % If 1, find identical units across sessions or oversplits
PrepareClusInfoparams.RedoUnitMatch = 1; % if 1, Redo unitmatch
PrepareClusInfoparams.SaveDir = SaveDir; % Save results here
PrepareClusInfoparams.tmpdatafolder = tmpdatafolder; % use this as a local directory (should be large enough to handle all sessions you want to combine)
PrepareClusInfoparams.separateIMRO = 0; % Run for every IMRO separately (for memory reasons this might be a good idea)
PrepareClusInfoparams.ReLoadAlways = 1; % If 1, SP & Clusinfo are always loaded from KS output
PrepareClusInfoparams.binsz = 0.01; %Binsz for unitmatch PSTHs
PrepareClusInfoparams.saveSp = 1; % Save SP struct for easy loading of preprocessed data
% UnitMatch Parameters:
% PrepareClusInfoparams.Scores2Include =
% {'AmplitudeSim','spatialdecaySim','WavformMSE','WVCorr','CentroidDist','CentroidVar','CentroidDistRecentered','TrajAngleSim','TrajDistSim'};
% %Full set. WavformSim is average of WVCorr and WavformMSE
PrepareClusInfoparams.Scores2Include = {'CentroidDist','WavformSim','CentroidDistRecentered','spatialdecayfitSim','AmplitudeSim','TrajAngleSim'}; %{'AmplitudeSim','spatialdecayfitSim','WavformSim','CentroidDist','CentroidVar','TrajAngleSim'}; % 
PrepareClusInfoparams.ApplyExistingBayesModel = 0; %If 1, use probability distributions made available by us
PrepareClusInfoparams.MakePlotsOfPairs = 0; % Plots pairs for inspection (UnitMatch)
PrepareClusInfoparams.AssignUniqueID = 1; % Assign UniqueID 
PrepareClusInfoparams.GoodUnitsOnly = 1; % Include only good untis in the UnitMatch analysis
%% All dependencies you want to add (you may need to download these, all available via github)
addpath(genpath(cd))

% Optional (depends on precise analysis):
addpath(genpath('C:\Users\EnnyB\Documents\GitHub\AP_histology'))
addpath(genpath('C:\Users\EnnyB\Documents\Github\rastermap'))
addpath(genpath('C:\Users\EnnyB\Documents\GitHub\allenCCF')) 
addpath(genpath('C:\Users\EnnyB\Documents\GitHub\bombcell')) % https://github.com/Julie-Fabre/bombcell, branch 'enny'

% Required (for using UnitMatch):
addpath(genpath('C:\Users\EnnyB\Documents\GitHub\spikes'))% please use https://github.com/EnnyvanBeest/spikes
addpath(genpath('C:\Users\EnnyB\Documents\GitHub\npy-matlab'))
addpath(genpath('C:\Users\EnnyB\Documents\GitHub\mtscomp'))  % https://github.com/int-brain-lab/mtscomp
addpath(genpath('C:\Users\EnnyB\Documents\GitHub\UnitMatch')) % Make sure to have this one fresh in the path (so run this last)

try
    % Python version to run python code in:
    pyversion('C:\Users\EnnyB\anaconda3\envs\pyks2\pythonw.exe') %Explanation on how to do this is provided in the README
catch ME
    disp(ME)
end

%% Automatic from here
%% PyKS - run pykilosort from Matlab/Python integration
RunPyKS2_FromMatlab

%% Runs unitmatch across all data from a mouse to generate a table
RunUnitMatchAllDataPerMouse

%% Across Mice Graphs
SummarizeAcrossMice