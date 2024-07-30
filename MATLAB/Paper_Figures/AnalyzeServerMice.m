% 
SaveDir = '\\znas.cortexlab.net\Lab\Share\UNITMATCHTABLES_ENNY_CELIAN_JULIE\FullAnimal_KSChanMap'; %H:\Ongoing\'%'H:\SfN_2022'; %%'E:\Data\ResultsOngoing' %
% SaveDir = '\\znas.cortexlab.net\Lab\Share\UNITMATCHTABLES_ENNY_CELIAN_JULIE\2ConsecutiveDays_KSChanMap\Stitched\'
% SaveDir = 'H:\UnitMatch\'
FromDate = datetime("2024-03-08 09:00:00");
AssignUnitDate = datetime("2024-03-14 10:00:00");

UMFiles = cell(1,0); % Define your UMfolders here or use below:
groupvec = nan(1,0);
SpGLXV = cell(1,0);
if ~exist('UMFiles') || isempty(UMFiles) % When using the example pipeline this may be useful:
    MiceOpt = dir(SaveDir);
    MiceOpt = arrayfun(@(X) X.name,MiceOpt,'Uni',0);
    MiceOpt(ismember(MiceOpt,{'.','..'})) = [];

    for midx = 1:numel(MiceOpt)
        fprintf('Reference %s...\n', MiceOpt{midx})
        % Identify all UM tables
        % tmpfile = dir(fullfile(SaveDir, MiceOpt{midx}, 'UnitMatch.mat'));

        tmpfile = dir(fullfile(SaveDir, MiceOpt{midx},'*','*','UnitMatch', 'UnitMatch.mat'));

        if isempty(tmpfile)
            continue
        end
        for id = 1:length(tmpfile)
            % if datetime(tmpfile(id).date) > FromDate % && any(cell2mat(cellfun(@(X) any(strfind(fullfile(tmpfile(id).folder,tmpfile(id).name),X)),UMFiles2Take,'Uni',0)))
                
                if datetime(tmpfile(id).date) < AssignUnitDate
                AssignUniqueID(fullfile(tmpfile(id).folder)) % REDO
                end
                % Check that these data are not too noisy
              
                % load(fullfile(tmpfile(id).folder,tmpfile(id).name),'UMparam')
                % for rid = 1:numel(UMparam.RawDataPaths)
                %    meta = ReadMeta2(UMparam.RawDataPaths{rid}.folder);
                %    SpGLXV = {SpGLXV{:} meta.appVersion};
                % end


                %             FolderParts = strsplit(tmpfile(id).folder,filesep);
                %             idx = find(ismember(FolderParts,MiceOpt{midx}));
                UMFiles = cat(2,UMFiles,fullfile(tmpfile(id).folder,tmpfile(id).name));
                groupvec = cat(2,groupvec,midx);
            % else
                % keyboard
            % end


        end
    end
    close all
end

Info  = DataSetInfo(UMFiles)
Info.RecSes
nanmean(cat(1,Info.nGoodUnits{:})./cat(1,Info.nTotalUnits{:}).*100)
nanstd(cat(1,Info.nGoodUnits{:})./cat(1,Info.nTotalUnits{:}).*100)


AUCExtract(UMFiles)


summaryMatchingPlots(UMFiles,{'UID1Liberal','UID1','UID1Conservative'},groupvec,1)
summaryFunctionalPlots(UMFiles, 'Corr', groupvec)

%
summaryFunctionalPlots_Part2(UMFiles, groupvec, 0)


%% For figure S3
res = summaryFunctionalPlots(UMFiles, 'Corr', groupvec, 0, 0);
resKS = summaryFunctionalPlots(UMFiles, 'Corr', groupvec, 1, 0);

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


%% Redo
if 0
    for midx = 1:length(UMFiles)

        load(UMFiles{midx})
        % UMparam.SaveDir = strrep(UMparam.SaveDir,'\\znas.cortexlab.net\Lab\Share\UNITMATCHTABLES_ENNY_CELIAN_JULIE\MonthApartNew\NonStitched\','\\znas.cortexlab.net\Lab\Share\UNITMATCHTABLES_ENNY_CELIAN_JULIE\MonthApartKSMap\NonStitched\')
        clusinfo = getClusinfo(UMparam.KSDir);

        for ksid = 1:numel(UMparam.KSDir)
            myClusFile = dir(fullfile(UMparam.KSDir{ksid}, 'channel_map.npy'));
            channelmaptmp = readNPY(fullfile(myClusFile(1).folder, myClusFile(1).name));

            myClusFile = dir(fullfile(UMparam.KSDir{ksid}, 'channel_positions.npy'));
            channelpostmp = readNPY(fullfile(myClusFile(1).folder, myClusFile(1).name));

            if size(channelpostmp,1) ~=384
                channelpos = nan(UMparam.nSavedChans-UMparam.nSyncChans,2);
                channelpos(channelmaptmp+1,:) = channelpostmp; % 0 to 1-index

                % KS 2 output (JFxxx) has missing channels, interpollate these
                channelpos(:,2) = fillmissing(channelpos(:,2),'linear'); % Not ideal but c'est la vie
                channelpos(:,1) = fillmissing(channelpos(:,1),'linear'); % Not ideal but c'est la vie


                UMparam.AllChannelPos{ksid} = channelpos;
            else
                UMparam.AllChannelPos{ksid} = channelpostmp;


            end
        end

        % Actual UnitMatch & Unique UnitID assignment
        [UniqueIDConversion, MatchTable, WaveformInfo, UMparam] = UnitMatch(clusinfo, UMparam);
        if UMparam.AssignUniqueID
            [UniqueIDConversion, MatchTable] = AssignUniqueID(UMparam.SaveDir);
        end

        % Evaluate (within unit ID cross-validation)
        EvaluatingUnitMatch(UMparam.SaveDir);

        % Function analysis
        ComputeFunctionalScores(UMparam.SaveDir,1)
        % Visualization
        PlotUnitsOnProbe(clusinfo,UMparam,UniqueIDConversion,WaveformInfo)

    end



end