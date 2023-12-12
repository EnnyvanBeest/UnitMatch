function Params = DefaultParametersExtractKSData(Params,KiloSortPaths)
% Checks if all required Paramseters are there, and will fill in default
% Paramseters if not wherever possible
%% Check if path information is given
if ~isfield(Params,'loadPCs')
    Params.loadPCs = 0;
end
if ~isfield(Params,'RunPyKSChronicStitched')
    Params.RunPyKSChronicStitched = 0;
end
if ~isfield(Params,'DecompressionFlag')
    Params.DecompressionFlag = 0; %if 1, uncompress data first if it's currently compressed
end
if ~isfield(Params,'DecompressLocal')
    Params.DecompressLocal = 1; %if 1, uncompress data first if it's currently compressed
end
if ~isfield(Params,'CleanUpTemporary')
    Params.CleanUpTemporary = 0; % Clean up temporary data
end
if ~isfield(Params,'tmpdatafolder')
    Params.tmpdatafolder = KiloSortPaths(1); % Directory to temporarily decompress data --> must
% be large enough!
end
if ~isfield(Params,'saveSp')
    Params.saveSp = 1; % Save prepared data?
end
if ~isfield(Params,'ReLoadAlways')
    Params.ReLoadAlways = 0; % Always reprep data, even if we already have a saved set?
end
if ~isfield(Params,'deNoise')
    Params.deNoise = 1; % Calls RemoveNoiseAmplitudeBased for average channel based removal of noise
end
if ~isfield(Params, 'MinRecordingDuration')
    Params.MinRecordingDuration = 10; % Clean up temporary folder
end
%% Quality metrics?
if ~isfield(Params,'RedoQM')
    Params.RedoQM = 0; %if 1, redo quality metrics if it already exists
end
if ~isfield(Params,'RunQualityMetrics')
    Params.RunQualityMetrics = 1; % If 1, Run the quality metrics
end
if ~isfield(Params,'InspectQualityMetrics')
    Params.InspectQualityMetrics = 0; % Inspect the quality metrics/data set using the GUI
end

%% Channel information

if ~isfield(Params,'extractSync')
    Params.extractSync = 0; % Only necessary if not extracted and using function scores
end

if ~isfield(Params,'nSavedChans')
    Params.nSavedChans = 385; % For Neuropixels typically 385 channels of which 1 is syncc
end
if ~isfield(Params,'nSyncChans')
    Params.nSyncChans = 1; % Number of sync channels
end

%% Inspection
if ~isfield(Params,'MakePlotsOfPairs')
    Params.MakePlotsOfPairs = 0; % Plots and saves matches for you to inspect
end
if ~isfield(Params,'GUI')
    Params.GUI = 0; % Use our GUI to flip through and curate our matches
end
if ~isfield(Params,'RunPyKSChronicStitched')
    Params.RunPyKSChronicStitched = 0; % whether stitched KS was used (concatenated recordings) (it will compare KS labels to UM labels)
end
%% Paramseters only used for compute functional scores
if ~isfield(Params,'ACGbinSize')
    Params.ACGbinSize = 1E-03; %
end
if ~isfield(Params,'ACGduration')
    Params.ACGduration = 1; % in seconds
end
if ~isfield(Params,'binsize')
    Params.binsize = 0.01; % in seconds
end

return
