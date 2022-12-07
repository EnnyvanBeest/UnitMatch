%% User Input
DataDir = {'H:\MatchingUnits\RawData'};%{'\\znas\Subjects','\\128.40.198.18\Subjects','\\zaru.cortexlab.net\Subjects'}%%'\\znas\Subjects' %' Check DataDir2Use
SaveDir = 'H:\MatchingUnits\Output\'
tmpdatafolder = 'H:\MatchingUnits\Tmp\'; % temporary folder for acute recordings - preferably SSD
ChronicTmpFolder = 'H:\MatchingUnits\Tmp\' %Temporary folder for chronic, needs to be big!
LocalDir = 'H:\MatchingUnits\KilosortOutput\';% 'E:\Data\KiloSortOutput';%
AllenCCFPath = 'C:\Users\EnnyB\Documents\MATLAB\allenCCF'
storevideopath=fullfile(tmpdatafolder,'Videos');
HistoFolder = 'E:\Data\Histology\';
SaveFiguresTo = 'E:\Data\Figures'
MiceOpt = {'EB014'};%{'AL032'};
nidq_sync_used = zeros(1,length(MiceOpt));
nidq_sync_used(ismember(MiceOpt,{'EB001','CB007','CB008'}))=1;
DataDir2Use = repmat(1,[1,length(MiceOpt)]); %in case you have multiple locations with data
ProbeType = repmat({'1_3b'},1,length(MiceOpt));
ProbeType(ismember(MiceOpt,{'AL032'}))={'2_4S'};
RecordingType = repmat({'Acute'},1,length(MiceOpt));
RecordingType(ismember(MiceOpt,{'AL032'}))={'Chronic'}; %EB014',

RunPyKSChronic = 0
CopyToTmpFirst = 0
DecompressLocal = 1; %if 1, uncompress data first if it's currently compressed
RedoQM = 1; %if 1, redo quality matrix if it already exists
RunQualityMatrix = 1; % If 1, Run the quality matrix
InspectQualityMatrix =1; % Inspect the quality matrix/data set using the GUI
MatchUnitsAcrossDays = 1; % Matching chronic recording using QM instead of using pyks chronic output
maxsessnr = 2; %max nr. sessions on a day (doesn't need to be accurate)
MinDist2Include = 85; %Including only trials in which mouse reached at least xcm for some statistics (past all stimuli?)
pretrialtime = 2; %take up to x seconds prior trial
posttrialtime = 2; % take up to x seconds post trial
UseCluster = 0; %If 1 some of the steps will not be done in Matlab, but rather parameter scripts/ bash files will be created to submit to cluster computing

%% All dependencies
addpath(genpath(cd))
addpath(genpath('C:\Users\EnnyB\Documents\GitHub\spikes'))%please use https://github.com/EnnyvanBeest/spikes
addpath(genpath('C:\Users\EnnyB\Documents\GitHub\npy-matlab'))
addpath(genpath('C:\Users\EnnyB\Documents\GitHub\AP_histology'))
addpath(genpath('C:\Users\EnnyB\Documents\GitHub\UnitMatch'))
addpath(genpath('C:\Users\EnnyB\Documents\GitHub\mtscomp'))  % https://github.com/int-brain-lab/mtscomp
addpath(genpath('C:\Users\EnnyB\Documents\GitHub\bombcell'))
addpath(genpath('C:\Users\EnnyB\Documents\Github\rastermap'))
addpath(genpath('C:\Users\EnnyB\Documents\GitHub\allenCCF'))

try
    % Python version to run python code in:
    pyversion('C:\Users\EnnyB\anaconda3\envs\pyks2\pythonw.exe')
catch ME
    disp(ME)
end

%% PyKS
RunPyKS2_FromMatlab

%% Average Probe Activity
AverageProbeActivity