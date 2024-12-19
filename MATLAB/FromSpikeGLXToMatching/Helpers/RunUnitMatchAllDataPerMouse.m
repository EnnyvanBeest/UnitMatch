%% Automated
% Load all data
% Find available datasets (always using dates as folders)
clear DateOpt
% %dd = arrayfun(@(X) fullfile(DataDir{DataDir2Use(X)},MiceOpt{X},'*-*'),1:length(MiceOpt),'UniformOutput',0);
DateOpt = arrayfun(@(X) dir(fullfile(DataDir{DataDir2Use(X)},MiceOpt{X},'*-*')),1:length(MiceOpt),'UniformOutput',0); % DataDir2Use = server
DateOpt = cellfun(@(X) X([X.isdir]),DateOpt,'UniformOutput',0);
DateOpt = cellfun(@(X) {X.name},DateOpt,'UniformOutput',0);

RemoveOldCopies = 1;

LogError = {}; % Keep track of which runs didn't work
if ~exist('PipelineParamsOri','var')
    PipelineParamsOri = PipelineParams;
end
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
            disp(['Finding all KS4 directories in ' myKsDir ', ' DateOpt{midx}{did}])
            % First try KS4
            tmpfiles = dir(fullfile(myKsDir,DateOpt{midx}{did},'**','KS4'));
            tmpfiles(cellfun(@(X) ismember(X,{'.','..'}),{tmpfiles(:).name})) = [];
            % Conver to string
            tmpfiles = arrayfun(@(X) fullfile(tmpfiles(X).folder,tmpfiles(X).name),1:length(tmpfiles),'uni',0);


            if isempty(tmpfiles)
                disp(['Trying pyKS directories in ' myKsDir ', ' DateOpt{midx}{did}])
                tmpfiles = dir(fullfile(myKsDir,DateOpt{midx}{did},'**','pyKS'));
                tmpfiles(cellfun(@(X) ismember(X,{'.','..'}),{tmpfiles(:).name})) = [];
                % Conver to string
                tmpfiles = arrayfun(@(X) fullfile(tmpfiles(X).folder,tmpfiles(X).name),1:length(tmpfiles),'uni',0);
            end
            % if isempty(tmpfiles)
            %     disp(['Trying other kilosort directories in ' myKsDir ', ' DateOpt{midx}{did}])
            % 
            %     tmpfiles = dir(fullfile(myKsDir,DateOpt{midx}{did},'**','kilosort*','**','spike_clusters.npy'));
            %     tmpfiles = arrayfun(@(X) fullfile(tmpfiles(X).folder),1:length(tmpfiles),'uni',0);
            % 
            % end
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
            AllKiloSortPaths(cell2mat(cellfun(@(X) any(strfind(X,'Chrc2Sessions')),AllKiloSortPaths,'UniformOutput',0))) = []; %Use separate days and match units via matlab script

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

    %% Remove old copies?
    UnitMatchExist = dir(fullfile(SaveDir,MiceOpt{midx},'**','UnitMatch.mat'));
    if ~isempty(UnitMatchExist)
        MouseAllDone = 1;
    else 
        MouseAllDone = 0;
    end
    for id = 1:numel(UnitMatchExist)
        if RemoveOldCopies && UnitMatchExist(id).date<FromDate
            delete(fullfile(UnitMatchExist(id).folder,'**'))
            rmdir(UnitMatchExist(id).folder)

            MouseAllDone = 0;
        end
    end

    if MouseAllDone
        continue
    end
    %% Prepare cluster information
    if exist('FromDate','var')
        PipelineParams.FromDate = FromDate;
    end
    PipelineParams = ExtractKilosortData(AllKiloSortPaths,PipelineParams);
     
    PipelineParams.RecType = RecordingType{midx};%

    % Remove empty ones
    EmptyFolders = find(cellfun(@isempty,PipelineParams.AllChannelPos));
    PipelineParams.AllChannelPos(EmptyFolders) = [];
    PipelineParams.AllProbeSN(EmptyFolders) = [];
    PipelineParams.RawDataPaths(EmptyFolders) = [];
    PipelineParams.KSDir(EmptyFolders) = [];
    AllKiloSortPaths = PipelineParams.KSDir;

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
                PipelineParams.SaveDir = fullfile(SaveDir,MiceOpt{midx},['Probe' num2str(ProbeOpt(PosComb(1,runid)))],['IMRO_' num2str(PosComb(2,runid))]);
                PipelineParams.AllChannelPos = PipelineParams.AllChannelPos(idx);
                PipelineParams.AllProbeSN = PipelineParams.AllProbeSN(idx);
                PipelineParams.RawDataPaths = PipelineParams.RawDataPaths(idx);
                PipelineParams.KSDir = AllKiloSortPaths(idx);

                % Convert2Yuanetal(PipelineParams.KSDir,fullfile('H:\MatchingUnits\Yuan\AcrossManyDays\',MiceOpt{midx}))
                % % to compare to Yuan (makes the data according to their
                % format, then you need to download their pipeline to run
                % it)
            end

            %% Prepare personal save/directory and decompressed data paths
            PipelineParams.SaveDir = fullfile(PipelineParams.SaveDir ,'UnitMatch');

            if isstruct(PipelineParams.RawDataPaths{1})
               
                if length(PipelineParams.RawDataPaths)==1
                    PipelineParams.AllDecompPaths = arrayfun(@(X) fullfile(PipelineParams.tmpdatafolder, strrep(X.name, 'cbin', 'bin')), PipelineParams.RawDataPaths{1}, 'Uni', 0);
                    PipelineParams.RawDataPaths = arrayfun(@(X) fullfile(X.folder,X.name), PipelineParams.RawDataPaths{1}, 'Uni', 0);
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
                    continue
                end
                UMtime = tic;
                %% Actual UnitMatch & Unique UnitID assignment
                [UniqueIDConversion, MatchTable, WaveformInfo, UMparam] = UnitMatch(clusinfo, UMparam);
                if UMparam.AssignUniqueID
                    [UniqueIDConversion, MatchTable] = AssignUniqueID(UMparam.SaveDir);
                end
                UMtime = toc(UMtime)

                % Save phy similarity info for individual sessions
                SaveWithinSessionProba_forPhy(UMparam.SaveDir)
              
                % Evaluate (within unit ID cross-validation)
                EvaluatingUnitMatch(UMparam.SaveDir);
    
                %% Function analysis
                ComputeFunctionalScores(UMparam.SaveDir,1)
                  % Visualization
                PlotUnitsOnProbe(clusinfo,UMparam,UniqueIDConversion,WaveformInfo)

                %% Figures
                if UMparam.MakePlotsOfPairs
                    DrawBlind = 0; %1 for blind drawing (for manual judging of pairs)
                    DrawPairsUnitMatch(UMparam.SaveDir,DrawBlind);
                    if UMparam.GUI
                        FigureFlick(UMparam.SaveDir)
                        pause
                    end
                end

                %% QM
                try
                    QualityMetricsROCs(UMparam.SaveDir);
                catch ME
                    disp(['Couldn''t do Quality metrics for ' MiceOpt{midx}])
                end
            
            % else
            %     % Get clusinfo
            %     clusinfo = getClusinfo(UMparam.KSDir);
            %     if ~any(clusinfo.Good_ID) || sum(clusinfo.Good_ID)<UMparam.minGoodUnits
            %         disp('No good units, continue')
            %         continue
            %     end
            %     load(fullfile(UMparam.SaveDir,'UnitMatch.mat'),'UMparam','UniqueIDConversion','WaveformInfo')
            %     %% Visualization
            %     PlotUnitsOnProbe(clusinfo,UMparam,UniqueIDConversion,WaveformInfo)
            % 
            %     % Evaluate (within unit ID cross-validation)
            %     EvaluatingUnitMatch(UMparam.SaveDir);
            % 
                
            end
           
            %%
            disp(['Preprocessed data for ' MiceOpt{midx} ' run  ' num2str(runid) '/' num2str(nRuns)])
            close all
        catch ME
            disp([MiceOpt{midx} ' run  ' num2str(runid) '/' num2str(nRuns) ' crashed... continue with others'])

            LogError = {LogError{:} [MiceOpt{midx} '_run' num2str(runid)]};
        end


    end
end
%
