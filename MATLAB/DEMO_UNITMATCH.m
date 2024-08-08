%% DEMO UNIT MATCH 

% Example data can be found here: https://doi.org/10.6084/m9.figshare.24305758.v1

%% README

% UnitMatch matches neurons within and across days using only their spatiotemporal waveforms. 
% To run, it needs a folder called 'RawWaveforms' for each recording, where there should be a NPY file for every cluster
% containing the average waveform (recommended of at least 500 spikes) for every recording channel 
% for every half of a recording for that cluster (spikeWidth x nRecordingChannels x 2).

% If you use spikeGLX or OpenEphys, the UnitMatch pipeline will extract those waveforms for you (using ExtractAndSaveAverageWaveforms.m). 

%% -- Preprocessing ---

% First, let's preprocess the data. 

%% Add required and optional paths and subpaths

GithubDir = 'C:\Users\user_name\Documents\GitHub'; % Github directory

% Required (for using UnitMatch):
addpath(genpath(fullfile(GithubDir,'spikes'))) % https://github.com/cortex-lab/spikes
addpath(genpath(fullfile(GithubDir,'npy-matlab'))) % https://github.com/kwikteam/npy-matlab
addpath(genpath(fullfile(GithubDir,'mtscomp'))) % https://github.com/int-brain-lab/mtscomp

% Advised (quality metrics for unit selection):
addpath(genpath(fullfile(GithubDir,'bombcell'))) % DOI: 10.5281/zenodo.8172822, https://github.com/Julie-Fabre/bombcell 

% UNITMATCH - Move to top of paths 
addpath(genpath(fullfile(GithubDir,'UnitMatch'))) % Make sure to have this one fresh in the path (so run this last)

%% User input

% This is the path where the results will be saved ('\\path\to\save\UnitMatch'), e.g.:
UMparam.SaveDir = 'D:\MatchingUnits\Output\UnitMatch'; 

% This is a cell array with a path to each recording's Kilosort output directory, where there should be a subfolder called 'RawWaveforms'. 
% N.B. if you want to use the functional score evaluation of UnitMatch, 'KSDir' should also contain the Kilosort output (e.g. spike times etc.)/
% Takes the form of "{'\\path\to\firstrecording','\\path\to\secondrecording','\\path\to\nthrecording'};", e.g.:  
 UMparam.KSDir = {'D:\MatchingUnits\Data\tmp\Mouse1\AL032\2019-11-21\Probe0\1','D:\MatchingUnits\Data\tmp\Mouse1\AL032\2019-11-22\Probe0\1'};  

%% Get recording information

% In this part, you will provide the paths for the raw data, and the channel position, for each recording.

% -- IF USING KILOSORT -- 

% If you use Kilosort, you can use the following. It extracts KS data and do some noise removal, 
% and will generate a "PreparedData.mat" file that will be useful in the next section. 
% RawDataPaths and AllChannelPos will be found automatically.
% Optionally, it will decompresse cbin to bin data and can use BOMBCELL quality metric to define good single units. 
UMparam = ExtractKilosortData(UMparam.KSDir, UMparam);

% You can also force the raw data paths as a 3rd argument: UMparam = ExtractKilosortData(KiloSortPaths, UMparam, RawDataPaths);
% If RawDataPaths is not given, raw data path will be obtained from Kilosort's output file "params.py".

% -- IF NOT --
% You can generate those:
% % A cell array with info on where to find or store the compressed recording (.cbin files OR .bin files):
% UMparam.RawDataPaths = {'\\path\to\firstrecording.cbin','\\path\to\secondrecording.cbin','\\path\to\nthrecording.cbin'};  
% % A cell array with info on where to find the decompressed recording (.bin files) -- typically a temporary folder: 
% UMparam.AllDecompPaths = {'\\path\to\firstrecording.bin','\\path\to\secondrecording.bin','\\path\to\nthrecording.bin'};  
% % The coordinates of every recording channel on the probe (e.g. nRecordingChannels x 2), in the same order as the channels were recorded:
% UMparam.AllChannelPos = {[RecordingSites_Recording1],[RecordingSites_Recording2]}; 

%% Get cluster information

% -- IF USING KILOSORT -- 

% If you use kilosort and have used "ExtractKilosortData" above to generate the "PreparedData.mat" files, you can use:
clusinfo = getClusinfo(UMparam.KSDir); 

% -- IF NOT --
% In this part, you will have to generate "clusinfo", a structure that contains basic information about the clusters:
% * cluster_id (e.g. kilosort output clus_id)
% * Good_ID: ones for units that should be included in the analysis (typically well isolated single units)
% * RecSesID: Recording Session ID
% * Probe: Which probe (if just 1, ones of numel cluster_id)
% * Depth: depth on probe (optional)
% * Shank: Which shank (optional)

%% -- Running UnitMatch --

% Now, we can run UnitMatch.
% (Note that waveform extraction will take place in UnitMatch.m if not performed before)

%% Load default parameters:

UMparam = DefaultParametersUnitMatch(UMparam);

%% UnitMatch algorithm:

% Obtain the probability for each pair of units of being a "match"
[UniqueIDConversion, MatchTable, WaveformInfo, UMparam] = UnitMatch(clusinfo, UMparam);

% Apply tracking algoritms based on UnitMatch's output probabilities
if UMparam.AssignUniqueID
    [UniqueIDConversion, MatchTable] = AssignUniqueID(UMparam.SaveDir);
end

%% -- Evaluating the output -- (optional)

%% If you want, you can run additional scripts to evaluate the results:

% Within session cross-validation
EvaluatingUnitMatch(UMparam.SaveDir); 

% Evaluating the matching using functional scores (only works when having access to Kilosort output, e.g. spike times etc. )
ComputeFunctionalScores(UMparam.SaveDir)

%% You can also visualize the units on the probe:

PlotUnitsOnProbe(clusinfo,UMparam,UniqueIDConversion,WaveformInfo)

%% Curation:
if UMparam.MakePlotsOfPairs
    DrawPairsUnitMatch(UMparam.SaveDir);
    if UMparam.GUI
        FigureFlick(UMparam.SaveDir)
        pause
    end
end

%% If using bombcell

% Evaluating how much quality metrics predict "matchability" (only works in combination with bombcell)
QualityMetricsROCs(UMparam.SaveDir); 

%% The following scripts work with several animals and plot summaries of the matching evaluation.

UMFiles = {fullfile(UMparam.SaveDir,'UnitMatch.mat')}; % cell containing a list of all the UnitMatch outputs you want to combine -- let's use one for now.
groupvec = 1; % How do you want to group these outputs? E.g., group all the outputs from the same animal into one

% Plots a summary of the functional metrics results
summaryFunctionalPlots(UMFiles, 'Corr', groupvec); 

% Plots a summary of the matching probability
summaryMatchingPlots(UMFiles); 
