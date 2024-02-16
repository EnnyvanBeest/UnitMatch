DataDir = {'H:\MatchingUnits\RawData'};% ;%Raw data folders, typically servers were e.g. *.cbin files are stored
SaveDir = 'H:\MatchingUnits\Output'; %'\\znas.cortexlab.net\Lab\Share\UNITMATCHTABLES_ENNY_CELIAN_JULIE\2ConsecutiveDays\Stitched';%'\\znas.cortexlab.net\Lab\Share\UNITMATCHTABLES_ENNY_CELIAN_JULIE\MonthApart\Stitched';%%'H:\MatchingUnits\Output\MonthApartStitched'% 'H:\MatchingUnits\Output\NotConcatenated';%'\\znas.cortexlab.net\Lab\Share\Celian\UnitMatch\MatchTables\NewSep27\MonthApart\Stitched'% %%;% %'H:\MatchingUnits\Output\ManyRecordings'%Folder where to store the results
tmpdatafolder = 'H:\OpenEphys_Example\Tmp'; % temporary folder for temporary decompression of data
KilosortDir = 'H:\MatchingUnits\KilosortOutput'; % '\\znas.cortexlab.net\Lab\Share\Enny\UnitMatch\KSComparisonSubset';%'\\znas.cortexlab.net\Lab\Share\Enny\UnitMatch\KilosortOutputMonthApart';%'H:\MatchingUnits\KilosortOutputMonthApart';%'\\znas.cortexlab.net\Lab\Share\Celian\UnitMatch\KilosortOutputMonthApart';% Kilosort output folder
GithubDir = 'C:\Users\EnnyB\Documents\GitHub'; % Github directory
PythonEXE = 'C:\Users\EnnyB\anaconda3\envs\pyks2_debug\pythonw.exe' % Python version to run python code in:
%% Information on experiments
MiceOpt = {'AL032','AV008','CB016','EB019','JF067'}; %'AL032', Add all mice you want to analyze
DataDir2Use = repmat(1,[1,length(MiceOpt)]); % In case you have multiple DataDir, index which directory is used for each mouse
RecordingType = repmat({'Chronic'},1,length(MiceOpt)); % And whether recordings were Chronic (default)
RecordingType(ismember(MiceOpt,{''}))={'Acute'}; %EB014', % Or maybe acute?

%% Parameters on how to prepare units/data for analysis
PipelineParams.RunPyKSChronicStitched = 0; % Default 0. if 1, run PyKS chronic recordings stitched when same IMRO table was used
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
    % Follow instructions on installing pykilosort in anaconda environment,
    % eg. https://github.com/int-brain-lab/pykilosort
    % Additionally run (from within this environment):
    % pip install matlab
    % pip install pyqt5-tools
    % If it doesn't work, please try the forked version on my github repo:
    % https://github.com/EnnyvanBeest/pykilosort/tree/UnitMatchPipeline
    % Python version to run python code in:
    pyversion(PythonEXE) %Explanation on how to do this is provided in the README

catch ME
    disp(ME)
end

%%
clear DateOpt
% %dd = arrayfun(@(X) fullfile(DataDir{DataDir2Use(X)},MiceOpt{X},'*-*'),1:length(MiceOpt),'UniformOutput',0);
DateOpt = arrayfun(@(X) dir(fullfile(DataDir{DataDir2Use(X)},MiceOpt{X},'*-*')),1:length(MiceOpt),'UniformOutput',0); % DataDir2Use = server
DateOpt = cellfun(@(X) X([X.isdir]),DateOpt,'UniformOutput',0);
DateOpt = cellfun(@(X) {X.name},DateOpt,'UniformOutput',0);
FromDate = datetime("2023-10-03 09:00:00");

LogError = {}; % Keep track of which runs didn't work
if ~exist('PipelineParamsOri','var')
    PipelineParamsOri = PipelineParams;
end
DriftOptions = -45:3:45;
AllFPFN = nan(length(MiceOpt),length(DriftOptions),2); %FP/FN

for midx = 1:length(MiceOpt)
    close all % to not overcrowd the graphics card
    PipelineParams = PipelineParamsOri; % Reset
    %% Loading data from kilosort/phy easily
    if ~isempty(KilosortDir)
        myKsDir = fullfile(KilosortDir,MiceOpt{midx});
        subksdirs = dir(fullfile(myKsDir,'**','Probe*')); %This changed because now I suddenly had 2 probes per recording
        if length(subksdirs)<1
            clear subksdirs
            subksdirs.folder = myKsDir; %Should be a struct array
            subksdirs.name = 'Probe0';
        end
        ProbeOpt = (unique({subksdirs(:).name}));

        myKsDir = fullfile(KilosortDir,MiceOpt{midx});
        % Check for multiple subfolders?
        AllKiloSortPaths = dir(fullfile(myKsDir,'**','channel_positions.npy'));
        AllKiloSortPaths=arrayfun(@(X) AllKiloSortPaths(X).folder,1:length(AllKiloSortPaths),'Uni',0);
        % Remove anything that contains the name 'noise'
        AllKiloSortPaths(cellfun(@(X) contains(X,'NOISE'),AllKiloSortPaths)) = [];
    else
        myKsDir = fullfile(DataDir{DataDir2Use(midx)},MiceOpt{midx});
        AllKiloSortPaths = [];
        for did = 1:length(DateOpt{midx})
            disp(['Finding all pyKS directories in ' myKsDir ', ' DateOpt{midx}{did}])
            tmpfiles = dir(fullfile(myKsDir,DateOpt{midx}{did},'**','pyKS'));
            tmpfiles(cellfun(@(X) ismember(X,{'.','..'}),{tmpfiles(:).name})) = [];
            % Conver to string
            tmpfiles = arrayfun(@(X) fullfile(tmpfiles(X).folder,tmpfiles(X).name),1:length(tmpfiles),'uni',0);
            AllKiloSortPaths = [AllKiloSortPaths, tmpfiles];
        end
    end

    if isempty(AllKiloSortPaths)
        display(['No data found for ' MiceOpt{midx} ', continue...'])
        continue
    end

    if strcmp(RecordingType{midx},'Chronic')
        if ~PipelineParams.RunPyKSChronicStitched %MatchUnitsAcrossDays
            disp('Unit matching in Matlab')
            AllKiloSortPaths(cell2mat(cellfun(@(X) any(strfind(X,'Chronic')),AllKiloSortPaths,'UniformOutput',0))) = []; %Use separate days and match units via matlab script
        else
            disp('Using chronic pyks option')
            AllKiloSortPaths = AllKiloSortPaths(cell2mat(cellfun(@(X) any(strfind(X,'Chronic')),AllKiloSortPaths,'UniformOutput',0))); %Use chronic output from pyks
        end
    end

    %% Create saving directory
    clear params
    if ~exist(fullfile(SaveDir,MiceOpt{midx}))
        mkdir(fullfile(SaveDir,MiceOpt{midx}))
    end
    if isempty(AllKiloSortPaths)
        disp(['No data found for ' MiceOpt{midx}])
        continue
    end

    %% Prepare cluster information
    PipelineParams = ExtractKilosortData(AllKiloSortPaths,PipelineParams);

    PipelineParams.RecType = RecordingType{midx};%

    % Remove empty ones
    EmptyFolders = find(cellfun(@isempty,PipelineParams.AllChannelPos));
    AllKiloSortPaths(EmptyFolders) = [];
    PipelineParams.AllChannelPos(EmptyFolders) = [];
    PipelineParams.AllProbeSN(EmptyFolders) = [];
    PipelineParams.RawDataPaths(EmptyFolders) = [];
    PipelineParams.KSDir = AllKiloSortPaths;

    %% Might want to run UM for separate IMRO tables & Probes (although UM can handle running all at the same time and takes position into account)
    if ~PipelineParams.separateIMRO
        RunSet = ones(1,length(AllKiloSortPaths)); %Run everything at the same time
        nRuns = 1;
    else
        % Extract different IMRO tables
        channelpositionMatrix = cat(3,PipelineParams.AllChannelPos{:});
        [UCHanOpt,~,idIMRO] = unique(reshape(channelpositionMatrix,size(channelpositionMatrix,1)*size(channelpositionMatrix,2),[])','rows','stable');
        UCHanOpt = reshape(UCHanOpt',size(channelpositionMatrix,1),size(channelpositionMatrix,2),[]);

        % Extract unique probes used
        [ProbeOpt,~,idProbe]  = unique([PipelineParams.AllProbeSN{:}]);
        PosComb = combvec(1:length(ProbeOpt),1:size(UCHanOpt,3)); % Possible combinations Probe X IMRO
        % Assign a number to each KS path related to PosComb
        if strcmp(RecordingType{midx},'Chronic')

            RunSet = nan(1,length(AllKiloSortPaths));
            for ksid = 1:length(AllKiloSortPaths)
                RunSet(ksid) = find(PosComb(2,:)==idIMRO(ksid) & PosComb(1,:)==idProbe(ksid));
            end
            nRuns = length(PosComb);
        elseif strcmp(RecordingType{midx},'Acute') % Run every recording separately
            RunSet = 1:length(AllKiloSortPaths);
            PosComb = cat(1,idProbe',RunSet);
            nRuns = length(PosComb);
        else
            disp('Unknown recording type')
        end

    end

    ORIParams = PipelineParams; % RESET

    %% Run UnitMatch
    for runid = 1:nRuns
        try
            PipelineParams = ORIParams; % RESET
            idx = find(RunSet==runid);
            if isempty(idx)
                continue
            end
            if ~PipelineParams.separateIMRO
                PipelineParams.SaveDir = fullfile(SaveDir,MiceOpt{midx},'AllProbes','AllIMRO');
                PipelineParams.KSDir = AllKiloSortPaths;
            else
                PipelineParams.SaveDir = fullfile(SaveDir,MiceOpt{midx},['Probe' num2str(PosComb(1,runid)-1)],['IMRO_' num2str(PosComb(2,runid))]);
                PipelineParams.AllChannelPos = PipelineParams.AllChannelPos(idx);
                PipelineParams.AllProbeSN = PipelineParams.AllProbeSN(idx);
                PipelineParams.RawDataPaths = PipelineParams.RawDataPaths(idx);
                PipelineParams.KSDir = AllKiloSortPaths(idx);
            end

            %% Prepare personal save/directory and decompressed data paths
            PipelineParams.SaveDir = fullfile(PipelineParams.SaveDir ,'UnitMatch');
            if isstruct(PipelineParams.RawDataPaths{1})
                if length(PipelineParams.RawDataPaths)==1
                    PipelineParams.AllDecompPaths = arrayfun(@(X) fullfile(PipelineParams.tmpdatafolder, strrep(X.name, 'cbin', 'bin')), PipelineParams.RawDataPaths{1}, 'Uni', 0);
                else
                    PipelineParams.AllDecompPaths = cellfun(@(X) fullfile(PipelineParams.tmpdatafolder, strrep(X.name, 'cbin', 'bin')), PipelineParams.RawDataPaths, 'Uni', 0);
                end
            else
                allDecompPaths_dirs = arrayfun(@(X) dir(Params.RawDataPaths{X}), 1:length(Params.RawDataPaths), 'Uni', 0);
                PipelineParams.AllDecompPaths = arrayfun(@(X) fullfile(Params.tmpdatafolder, strrep(allDecompPaths_dirs{X}.name, 'cbin', 'bin')), 1:length(allDecompPaths_dirs), 'Uni', 0);
            end

            %% Load other (Default) parameters
            UMparam = DefaultParametersUnitMatch(PipelineParams);
            UnitMatchExist = dir(fullfile(UMparam.SaveDir,'UnitMatch.mat'));

            if isempty(UnitMatchExist) || PipelineParams.RedoUnitMatch || UnitMatchExist.date<FromDate

                %% Get clusinfo
                clusinfo = getClusinfo(UMparam.KSDir);
                if ~any(clusinfo.Good_ID) || sum(clusinfo.Good_ID)<UMparam.minGoodUnits
                    disp('No good units, continue')
                end


                %% Actual UnitMatch & Unique UnitID assignment
                % DriftOptions = -45:3:45;
                % AllFPFN = nan(length(MiceOpt),length(DriftOptions),2); %FP/FN

                for shufid = 1:length(DriftOptions)
                    UMparam.ImposeDrift = DriftOptions(shufid); % impose random drift % Can only realistically do multiples of 1.5micron
                    [UniqueIDConversion, MatchTable, WaveformInfo, UMparam] = UnitMatch(clusinfo, UMparam);
                    load(fullfile(UMparam.SaveDir,'UnitMatchModel.mat')) % Load best model
                    AllFPFN(midx,shufid,1) = nanmean(BestMdl.FalsePositiveEstimate);
                    AllFPFN(midx,shufid,2) = nanmean(BestMdl.FalseNegativeEstimate);
                    close all
                end
            end
        end
    end
end

figure; subplot(1,2,1); plot(DriftOptions,AllFPFN(:,:,1)); xlabel('Drift'); ylabel('False Positives'); makepretty; offsetAxes
subplot(1,2,2); h = plot(DriftOptions,AllFPFN(:,:,2)); xlabel('Drift'); ylabel('False negatives'); makepretty; offsetAxes
legend([h(:)],MiceOpt)
linkaxes

save(fullfile(SaveDir,'DriftvsFPFN.mat'),'DriftOptions','AllFPFN','MiceOpt')
