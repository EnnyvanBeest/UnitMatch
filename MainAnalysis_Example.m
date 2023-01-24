%% User Input
%% Path information
DataDir = {'H:\MatchingUnits\RawData'};%{'\\znas\Subjects','\\128.40.198.18\Subjects','\\zaru.cortexlab.net\Subjects'}%%'\\znas\Subjects' %' Check DataDir2Use
SaveDir = 'H:\MatchingUnits\Output\'
tmpdatafolder = 'H:\MatchingUnits\Tmp\'; % temporary folder 
KilosortDir = 'H:\MatchingUnits\KilosortOutput\';% 'E:\Data\KiloSortOutput';%
AllenCCFPath = 'C:\Users\EnnyB\Documents\MATLAB\allenCCF'; % Path to allen common coordinate framework
storevideopath=fullfile(tmpdatafolder,'Videos');
HistoFolder = 'E:\Data\Histology\'; % Necessary when aligning to histology

%% Information on experiments
MiceOpt = {'AL032'}; % Add all mice you want to analyse
nidq_sync_used = zeros(1,length(MiceOpt)); % Was an external nidq used for syncing (typically sync feeds directly into IMEC)
nidq_sync_used(ismember(MiceOpt,{'EB001','CB007','CB008'}))=1; % Except for these mice...
DataDir2Use = repmat(1,[1,length(MiceOpt)]); % In case you have multiple DataDir, index which directory is used for each mouse
% DataDir2Use(ismember(MiceOpt,{'EB001','EB002','EB003','EB004','EB005','CB007','CB008','AL056'}))=2; 
ProbeType = repmat({'1_3b'},1,length(MiceOpt)); % WE need to know which probe type was used for each mice. Default:
ProbeType(ismember(MiceOpt,{'AL032'}))={'2_4S'}; % Change to 2_4Shank for these mice
RecordingType = repmat({'Acute'},1,length(MiceOpt)); % And whether recordings were acute (default)
RecordingType(ismember(MiceOpt,{'AL032','EB014'}))={'Chronic'}; %EB014', % Or maybe Chronic?

%% Parameters on how to prepare units/data for analysis
PrepareClusInfoparams.RunPyKSChronicStitched = 1; % if 1, run PyKS chronic recordings stitched when same IMRO table was used
PrepareClusInfoparams.CopyToTmpFirst = 1; % If 1, copy data to local first, don't run from server (= advised!)
PrepareClusInfoparams.DecompressLocal = 1; % If 1, uncompress data first if it's currently compressed (= necessary for unitmatch and faster for QualityMetrics))
PrepareClusInfoparams.RedoQM = 0; %if 1, redo quality metrics if it already exists
PrepareClusInfoparams.RunQualityMetrics = 1; % If 1, Run the quality metrics (Bombcell @JulieFabre)
PrepareClusInfoparams.loadPCs=1; % Do we need the PCs for data analysis (YES IF QM!)? If not you save a lot of time keeping this at 0
PrepareClusInfoparams.InspectQualityMatrix =0; % If 1, Inspect the quality matrix/data set using the GUI (manual inspection)
PrepareClusInfoparams.UnitMatch = 1; % If 1, find identical units across sessions or oversplits
PrepareClusInfoparams.RedoUnitMatch = 0; % if 1, Redo unitmatch
PrepareClusInfoparams.SaveDir = SaveDir; % Save results here
PrepareClusInfoparams.tmpdatafolder = tmpdatafolder; % use this as a local directory (should be large enough to handle all sessions you want to combine)

%% Parameters for further analysis
maxsessnr = 2; %max nr. sessions on a day (doesn't need to be accurate)
MinDist2Include = 85; %Including only trials in which mouse reached at least xcm for some statistics (past all stimuli?)
pretrialtime = 2; %take up to x seconds prior trial
posttrialtime = 2; % take up to x seconds post trial

%% All dependencies you want to add (you may need to download these, all available via github)
addpath(genpath(cd))

% Optional (depends on precise analysis):
addpath(genpath('C:\Users\EnnyB\Documents\GitHub\AP_histology'))
addpath(genpath('C:\Users\EnnyB\Documents\Github\rastermap'))
addpath(genpath('C:\Users\EnnyB\Documents\GitHub\allenCCF')) 

% Required (for using UnitMatch):
addpath(genpath('C:\Users\EnnyB\Documents\GitHub\spikes'))% please use https://github.com/EnnyvanBeest/spikes
addpath(genpath('C:\Users\EnnyB\Documents\GitHub\npy-matlab'))
addpath(genpath('C:\Users\EnnyB\Documents\GitHub\mtscomp'))  % https://github.com/int-brain-lab/mtscomp
addpath(genpath('C:\Users\EnnyB\Documents\GitHub\bombcell')) % https://github.com/Julie-Fabre/bombcell, branch 'enny'
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

%% Average Probe Activity - example script which shows how you can use prepareclusinfo & UnitMatch in data analysis
AverageProbeActivity