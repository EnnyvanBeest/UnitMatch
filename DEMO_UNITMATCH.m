%% DEMO UNIT MATCH 

%% READ ME
% If you do not use the suggested pipeline to extract raw waveforms (e.g. you don't use Neuropixels/SpikeGLX), make
% sure your 'KiloSortPaths' contains a subfolder called 'RawWaveforms'. There
% should be a NPY file for every cluster with the dimensions of
% UMparam.spikeWidth X nRecordingChannels X 2 (1 for each half of a
% recording). This should contain the average waveform (recommended of at
% least 500 spikes) for every recording channel for every half of a
% recording for that cluster.

%% User input: 
SaveDir = '\\path\to\save\UnitMatch'; % Recommended to use end this path with \Probe0\IMRO_1\ if more probes/IMRO tables were used or \AllProbes\AllIMRO\ otherwise
KiloSortPaths = {'\\path\to\firstrecording','\\path\to\secondrecording','\\path\to\nthrecording'};  % This is a cell array with a path, in the path there should be a subfolder called 'RawWaveforms'. 

% The following user input can also be automatically extracted using [[]]]
% If you want to use the functional score evaluation of UnitMatch this should also contain typical 'Kilosort output', (e.g. spike times etc.)
RawDataPaths = {'\\path\to\firstrecording','\\path\to\secondrecording','\\path\to\nthrecording'};  % This is a cell array with info on where to find the decompressed recording (.bin files) --> Necessary when you want UnitMatch to do waveform extraction
channelpos = {[RecordingSites_Recording1],[RecordingSites_Recording2]}; % These are coordinates of every recording channel on the probe (e.g. nRecordingChannels x 2)
clusinfo = struct; % Note, this can be kilosort input, 
% - clusinfo (this is a struct that contains per unit the following information):
% * cluster_id (e.g. kilosort output clus_id)
% * Good_ID: ones for units that should be included in the analysis
% * RecSesID: Recording Session ID
% * Coordinates: Typically 3D coordinates per unit
% * Depth: depth on probe
% * Shank: Which shank 
% * Probe: Which probe

% Params = ExtractKilosortData(KiloSortPaths, Params, RawDataPaths) 
RawDataPaths = Params.RawDataPaths; % This is a cell array with info on where to find the raw recording (.bin files) --> Necessary when you want UnitMatch to do waveform extraction
channelpos = Params.AllChannelPos; % These are coordinates of every recording channel on the probe (e.g. nRecordingChannels x 2)

%% Load default parameters
UMparam = DefaultParametersUnitMatch(SaveDir,KiloSortPaths,RawDataPaths,channelpos);

%% UnitMatch algorithm:
[UniqueIDConversion, MatchTable, WaveformInfo, UMparam] = UnitMatch(clusinfo, UMparam);
tmpfile = dir(fullfile(UMparam.SaveDir,'UnitMatch.mat'));
if UMparam.AssignUniqueID
    AssignUniqueID(fullfile(tmpfile.folder,tmpfile.name));
end

%% Automatic evaluation:
EvaluatingUnitMatch(UMparam.SaveDir); % Within session cross-validation
QualityMetricsROCs(UMparam.SaveDir); % Only works in combination with BOMBCELL
ComputeFunctionalScores(UMparam.SaveDir) % Only works when having access to Kilosort output (e.g. spike times etc.) 

%% Curation:
if UMparam.MakePlotsOfPairs
    DrawPairsUnitMatch(UMparam.SaveDir);
    if UMparam.GUI
        FigureFlick(UMparam.SaveDir)
        pause
    end
end