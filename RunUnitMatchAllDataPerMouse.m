%% Automated
% Load all data
% Find available datasets (always using dates as folders)
clear DateOpt
for idx = 1:length(DataDir)
    DateOpt{idx} = cellfun(@(X) dir(fullfile(DataDir{idx},X,'*-*')),MiceOpt(DataDir2Use==idx),'UniformOutput',0);
end
DateOpt = cat(2,DateOpt{:});
DateOpt = cellfun(@(X) X([X.isdir]),DateOpt,'UniformOutput',0);
DateOpt = cellfun(@(X) {X.name},DateOpt,'UniformOutput',0);
PlotLFP = 0;
plotUnitActivity = 0;
GoodUnits = cell(1,length(MiceOpt));
for midx = 1:length(MiceOpt)
    %% Loading data from kilosort/phy easily
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
    subsesopt = dir(fullfile(myKsDir,'**','channel_positions.npy'));
    subsesopt=arrayfun(@(X) subsesopt(X).folder,1:length(subsesopt),'Uni',0);
    if strcmp(RecordingType{midx},'Chronic')
        if ~PrepareClusInfoparams.RunPyKSChronicStitched %MatchUnitsAcrossDays
            disp('Unit matching in Matlab')
            subsesopt(cell2mat(cellfun(@(X) any(strfind(X,'Chronic')),subsesopt,'UniformOutput',0))) = []; %Use separate days and match units via matlab script
        else
            disp('Using chronic pyks option')
            subsesopt = subsesopt(cell2mat(cellfun(@(X) any(strfind(X,'Chronic')),subsesopt,'UniformOutput',0))); %Use chronic output from pyks
        end
    end

    % Copy file and then run pykilosort on it
    channelposition = cell(1,length(subsesopt));
    for id = 1:length(subsesopt)
        tmpfile = dir(fullfile(subsesopt{id},'channel_positions.npy'));
        if isempty(tmpfile)
            continue
        end
        channelposition{id} = readNPY(fullfile(tmpfile.folder,tmpfile.name));
    end
    subsesopt(cellfun(@isempty,channelposition))=[];
    channelposition(cellfun(@isempty,channelposition))=[];
    subsesoptAll = subsesopt;
    ImroCount = 1;
    IMROTableOpt = {};
    ChannelPosOpt = {};
    subsesoptGroups = {};
    for id = 1:length(channelposition)
        id1 = find(cell2mat(cellfun(@(X) all(channelposition{id}(:)==X(:)),ChannelPosOpt,'Uni',0)));
        if ~isempty(id1)
            subsesoptGroups{id1} = [subsesoptGroups{id1} id];
        else
            IMROTableOpt = {IMROTableOpt{:} ['IMRO' num2str(ImroCount)]};
            subsesoptGroups{ImroCount} = id;
            ChannelPosOpt{ImroCount} = channelposition{id};
            ImroCount = ImroCount+1;
        end
    end

    %% Create saving directoryed
    clear params
    thisIMRO = '';
    thisdate = [];
    if ~exist(fullfile(SaveDir,MiceOpt{midx}))
        mkdir(fullfile(SaveDir,MiceOpt{midx}))
    end
    PrepareClusInfoparams.SaveDir = fullfile(SaveDir,MiceOpt{midx});

    %% Prepare cluster information
    PrepareClusInfoparams = PrepareClusInfo(subsesoptAll,PrepareClusInfoparams);
    PrepareClusInfoparams.RecType = RecordingType{midx};%
    %% Run UnitMatch
    UnitMatchExist = dir(fullfile(PrepareClusInfoparams.SaveDir,'**','UnitMatch.mat'));
    if isempty(UnitMatchExist) || PrepareClusInfoparams.RedoUnitMatch
       UMparam = RunUnitMatch(subsesoptAll,PrepareClusInfoparams);
    else
        UMparam = PrepareClusInfoparams;
        UMparam.SaveDir = fullfile(PrepareClusInfoparams.SaveDir,'UnitMatch');
    end
    %% Function analysis
    ComputeFunctionalScores(UMparam.SaveDir)

    %% Figures
    if UMparam.MakePlotsOfPairs
        DrawBlind = 0; %1 for blind drawing (for manual judging of pairs)
        DrawPairsUnitMatch(UMparam.SaveDir,DrawBlind);
    end
    
    %% Evaluate (within unit ID cross-validation)
    EvaluatingUnitMatch(UMparam.SaveDir);

    %% QM
    QualityMetricsROCs(UMparam.SaveDir);


    %% 
    disp(['Preprocessed data for ' MiceOpt{midx}])
    

end

