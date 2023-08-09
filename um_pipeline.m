%% ~~ Example unit match pipeline ~~
% 1. Adjust the paths in the following sections:
%         'Where to save data '
%         'Information on mice and recording types'
%         'get all raw ephys and kilosort directories'
%         and the parameters in um_runUnitMatch
%
% 2. Run!
%         This pipeline will:
%           (1) load your ephys data, 
%           (2) decompress your raw data if it is in .cbin format 
%           (3) run quality metrics (bombcell) on your data and save the output
%           (4) run unit match on your data and save the output
%           (5) bring up summary plots and a GUI to flip through classified cells.
%         The first time, this pipeline will be significantly slower (10-20' more)
%         than after because it extracts raw waveforms. Subsequent times these
%         pre-extracted waveforms are simply loaded in.
%
% 3. Check the output plots and output of `EvaluatingUnitMatch`. 
%         You may want to adjust/fine-tune the following:
%         - Which metrics are used in match. They are defined in PrepareClusInfoparams.Scores2Include. 
%         You want to remove any metrics with low AUC and/or add any metrics with a high AUC
%         - the prior used in the model (e.g. 'Acute' vs 'Chronic') 

%% Where to save data - CHANGE THESE PATHS
SaveDir = '/home/netshare/zinu/JF067/'; % Folder where to store the results
tmpdatafolder = '/media/julie/ExtraHD/data_temp'; % temporary folder for temporary decompression of data 
user = 'julie';

%% Information on mice and recording types - CHANGE THESE MOUSE NAMES AND RECORDING TYPES
MiceOpt = {'JF067'}; % Add all mice you want to analyze
RecordingType(ismember(MiceOpt,{'JF067'}))={'Chronic'}; % 'Acute' or 'Chronic'

for iMouse = 1:size(MiceOpt,2)

    %% Get all raw ephys and kilosort directories - CHANGE THESE PATHS
    ephys_dirs = {path_to_all_ephys_dirs}; % this should be a cell array, with each cell containing the path to your raw .bin, .cbin or .dat data. 
    kilosort_dirs = {path_to_all_kilosort_dirs}; % this should be a cell array, with each cell containing the path to your kilosort output files 
    thisRecordingType = RecordingType{iMouse};
    mouseName = MiceOpt{iMouse};

    %% Run UnitMatch
    [PrepareClusInfoparams, UMparam, UniqueIDConversion, MatchTable, WaveformInfo]  = um_runUnitMatch(kilosort_dirs, ephys_dirs, SaveDir, tmpdatafolder, thisRecordingType, mouseName);  

    %% Evaluate UnitMatch's output
    EvaluatingUnitMatch(SaveDir);
    ComputeFunctionalScores(SaveDir)

    %% GUI: look at matches and evaluate them
    DrawBlind = 0; %1 for blind drawing (for manual judging of pairs)
    DrawPairsUnitMatch(SaveDir, DrawBlind);
    
    % GUI guide:
    %   Right arrow: next pair
    %   Left arrow: previous pair
    %   Up arrow: label as match
    %   Down arrow: label as non-match
    %   o: label as I don't know (=uncurated)
    %   m: go to first uncurated pair
    %   p: go to pair number
    %   s: SAVE
    recompute = 1;

    FigureFlick(SaveDir,user,recompute)
    % your labels are saved in MatchTable.<user>
end