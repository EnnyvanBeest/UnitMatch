%% Example pipeline going from raw Neuropixels data recorded with SpikeGLX to curated (with Bombcell) units identified with UnitMatch

%% User Input
%% Path information
DataDir = {'H:\MatchingUnits\RawData'};%{'H:\MatchingUnits\RawDataMonthApart'};%  ;%Raw data folders, typically servers were e.g. *.cbin files are stored
SaveDir = 'H:\MatchingUnits\Output';%'H:\MatchingUnits\OutputMonthApart';%'\\znas.cortexlab.net\Lab\Share\UNITMATCHTABLES_ENNY_CELIAN_JULIE\Learning_Striatum_new';%'H:\MatchingUnits\OutputMonthApart'; %'\\znas.cortexlab.net\Lab\Share\UNITMATCHTABLES_ENNY_CELIAN_JULIE\2ConsecutiveDays\Stitched';%'\\znas.cortexlab.net\Lab\Share\UNITMATCHTABLES_ENNY_CELIAN_JULIE\MonthApart\Stitched';%%'H:\MatchingUnits\Output\MonthApartStitched'% 'H:\MatchingUnits\Output\NotConcatenated';%'\\znas.cortexlab.net\Lab\Share\Celian\UnitMatch\MatchTables\NewSep27\MonthApart\Stitched'% %%;% %'H:\MatchingUnits\Output\ManyRecordings'%Folder where to store the results
tmpdatafolder = 'D:\tmpdata'; % temporary folder for temporary decompression of data 
KilosortDir = 'H:\MatchingUnits\KilosortOutput';%'H:\MatchingUnits\KilosortOutputMonthApart';%'\\znas.cortexlab.net\Lab\Share\UNITMATCHTABLES_ENNY_CELIAN_JULIE\Learning_Striatum_new\KS4\';%'H:\MatchingUnits\KilosortOutput';% '\\znas.cortexlab.net\Lab\Share\Enny\UnitMatch\KSComparisonSubset';%'\\znas.cortexlab.net\Lab\Share\Enny\UnitMatch\KilosortOutputMonthApart';%'H:\MatchingUnits\KilosortOutputMonthApart';%'\\znas.cortexlab.net\Lab\Share\Celian\UnitMatch\KilosortOutputMonthApart';% Kilosort output folder
GithubDir = 'C:\Users\EnnyB\Documents\GitHub'; % Github directory
PythonEXE = 'C:\Users\EnnyB\anaconda3\envs\kilosort\pythonw.exe' % Python version to run python code in:

%% Information on experiments
MiceOpt = {'AL032'};%{'JF067'};%{'AL032','AV008','CB016','EB019','JF067'}; %'AL032', Add all mice you want to analyze
DataDir2Use = repmat(1,[1,length(MiceOpt)]); % In case you have multiple DataDir, index which directory is used for each mouse
RecordingType = repmat({'Chronic'},1,length(MiceOpt)); % And whether recordings were Chronic (default)
RecordingType(ismember(MiceOpt,{''}))={'Acute'}; %EB014', % Or maybe acute?

%% Parameters on how to prepare units/data for analysis
PipelineParams.RunPyKSChronicStitched = 0; % Default 0. if 1, run PyKS chronic recordings stitched when same IMRO table was used
PipelineParams.CopyToTmpFirst = 1; % If 1, copy data to local first, don't run from server (= advised!)
PipelineParams.DecompressLocal = 1; % If 1, uncompress data first if it's currently compressed (= necessary for unitmatch and faster for QualityMetrics)

% Storing preprocessed data?
PipelineParams.ExtractNewDataNow = 0; % If data is not (yet) extracted, don't bother for now if 0
PipelineParams.ReLoadAlways = 0;% If 1, SP & Clusinfo are always loaded from KS output, if 2 , only when date is from before FromDate
PipelineParams.saveSp = 1; % Save SP struct for easy loading of preprocessed data
PipelineParams.binsz = 0.01; %Bin size for PSTHs in seconds

% Quality Metrics
PipelineParams.RunQualityMetrics = 1; % If 1, Run the quality metrics (Bombcell @JulieFabre)
PipelineParams.RedoQM = 0; %if 1, redo quality metrics if it already exists (N.B. only if reloading data!)
PipelineParams.InspectQualityMetrics = 0; % If 1, Inspect the quality matrix/data set using the GUI (manual inspection)
PipelineParams.loadPCs = 0; % Only necessary when computiong isoluation metrics/drift in QM. You save a lot of time keeping this at 0

% UnitMatch
PipelineParams.UnitMatch = 1; % If 1, find identical units across sessions or oversplits in a fast and flexible way
PipelineParams.RedoUnitMatch = 1; % if 1, Redo unitmatch
PipelineParams.separateIMRO = 1; % Run for every IMRO separately (for memory reasons or when having multiple probes this might be a good idea)
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
PipelineParams.MakePlotsOfPairs =0; % Plots pairs for inspection (UnitMatch)
PipelineParams.GUI = 0; % Flick through and do manual curation of matching - only works if MakePlotsofPairs = 1

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
    % Follow instructions on installing kilosort in anaconda environment,
    % eg. https://github.com/MouseLand/Kilosort
    % Additionally run (from within this environment):
    % Python version to run python code in:
    pyversion(PythonEXE) %Explanation on how to do this is provided in the README

catch ME
    disp(ME)
end


%% Actual pipeline
%% PyKS - run pykilosort from Matlab/Python integration
if PipelineParams.ExtractNewDataNow
    RunKS4_FromMatlab
end

%% Runs unitmatch across all data from a mouse to generate a table
RunUnitMatchAllDataPerMouse

%% Across Mice Graphs
% SummarizeAcrossMice
FromDate = datetime("2024-02-26 09:00:00");
UMFiles = cell(1,0); % Define your UMfolders here or use below:
groupvec = nan(1,0);
if ~exist('UMFiles') || isempty(UMFiles) % When using the example pipeline this may be useful:
    countid = 0;
    for midx = 1:length(MiceOpt)
        fprintf('Reference %s...\n', MiceOpt{midx})
        % Identify all UM tables
        tmpfile = dir(fullfile(SaveDir, MiceOpt{midx},'*','*','UnitMatch', 'UnitMatch.mat'));
        if isempty(tmpfile) 
            continue
        end
        countid = countid+1;

        for id = 1:length(tmpfile)
            if datetime(tmpfile(id).date) > FromDate % && any(cell2mat(cellfun(@(X) any(strfind(fullfile(tmpfile(id).folder,tmpfile(id).name),X)),UMFiles2Take,'Uni',0)))
                %             FolderParts = strsplit(tmpfile(id).folder,filesep);
                %             idx = find(ismember(FolderParts,MiceOpt{midx}));
                UMFiles = cat(2,UMFiles,fullfile(tmpfile(id).folder,tmpfile(id).name));
                groupvec = cat(2,groupvec,countid);
%             elseif datetime(tmpfile(id).date)<FromDate
%                 rmdir(tmpfile(id).folder,'s')

            end
        end
    end
end
res = summaryFunctionalPlots(UMFiles, 'Corr', groupvec, 0, 0);
summaryMatchingPlots(UMFiles,{'UID1Liberal','UID1','UID1Conservative'},groupvec,1)

% resKS = summaryFunctionalPlots(UMFiles, 'Corr', groupvec, 1, 0);

%% Compare with KS stitched performance (Fig S3e/f)
figure('name','KS versus UM')
fnames = fieldnames(res.FPSum);
cols = distinguishable_colors(length(MiceOpt));

for fid = 1:numel(fnames)
    subplot(ceil(sqrt(numel(fnames))),round(sqrt(numel(fnames))),fid)
    hold on
    for midx = 1:length(MiceOpt)
        tmp1 = squeeze(res.FPSum.(fnames{fid}).AUC{midx}(1,:,:));
        tmp2 = squeeze(resKS.FPSum.(fnames{fid}).AUC{midx}(1,:,:));
        days = res.deltaDays{midx};
        if length(MiceOpt)==1
        scatter(tmp1(:),tmp2(:),35,days(:),'filled')
            colormap(copper)

        else
        scatter(tmp1(:),tmp2(:),35,cols(midx,:),'filled')
        colormap(cols)

        end
    end
    makepretty
    offsetAxes
    xlabel('UnitMatch')
    ylabel('Kilosort')
    xlim([0.5 1])
    ylim([0.5 1])
    line([0.5 1],[0.5 1],'color',[0.2 0.2 0.2])
    if length(MiceOpt)>1 & fid == 1
        hc = colorbar;

        hc.Ticks = linspace(0,1,length(MiceOpt));
        hc.TickLabels = MiceOpt;
    end

    title(fnames{fid})
end

%%
[unitPresence, unitProbaMatch, days, EPosAndNeg] = summaryMatchingPlots(UMFiles,{'UID1','ID1'},groupvec,1)

trackWithFunctionalMetrics(UMFiles)



%% 
