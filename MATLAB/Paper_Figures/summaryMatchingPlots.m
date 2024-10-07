function [unitPresence, unitProbaMatch, days, EPosAndNeg, DataSetInfo, pSinceAppearance, popCorr_Uni, popAUC_Uni] = summaryMatchingPlots(UMFiles,UIDtoUse,groupVector,computeFuncScores)

%% Settings
PlotIndividualMice = 0;
% Drifts2Exclude = 50; % more than 50 micron, expected to lose neurons. This has nothing to do with software but with quality of recordings
MinNumsUnit = -inf; % should have at least this amount of neurons, otherwise not a good session
try
    savedir = UMFiles{1}(1:find(UMFiles{1}~=UMFiles{end},1,'first')-1);
catch
    warning('can''t create save dir')
    savedir = [];
end
if isempty(savedir)
    savedir = fileparts(UMFiles{1});
end
RecSes = 0;
%% Initialize
if nargin<2
    UIDtoUse = {'UID1'};
end
if ~iscell(UIDtoUse)
    UIDtoUse = {UIDtoUse};
end
if any(ismember(UIDtoUse,'ID1'))
    warning('Assuming stitched recording?')
end
if nargin<3 || ~exist('groupVector','var') || isempty(groupVector)
    groupVector = 1:length(UMFiles);
end

if nargin<4
    computeFuncScores = 0;
end
groups = unique(groupVector);
groupColor = distinguishable_colors(max(groups)+1);

days = cell(1, length(UMFiles));
deltaDays = cell(1, length(UMFiles));
deltaDaysUni = cell(1, length(UMFiles));
unitPresence = cell(length(UIDtoUse), length(UMFiles));
unitProbaMatch = cell(length(UIDtoUse), length(UMFiles));
popProbaMatch = cell(length(UIDtoUse), length(UMFiles));

% WIthin session check
FalsePositives = nan(3,length(UMFiles)); % Liberal/intermedi/conserv
FalseNegatives = nan(3,length(UMFiles));
EPosAndNeg = nan(2,length(UMFiles));
nUnitsPerRec = cell(1,length(UMFiles));

% deltaDaysBinsOri = [0.1 1 2 3 4 5 10 20 50 100 inf];
% deltaDaysBinsOri = [0.1 1 2 3 4 5 10 15 20 30 40 50 75 100 125 150 inf];
deltaDaysBinsOri = 2.^(-1:8);
deltaDaysBins = [-deltaDaysBinsOri(end:-1:1)-0.1,deltaDaysBinsOri(1:end)];
deltaBinsVec =[-deltaDaysBinsOri(end:-1:1),deltaDaysBinsOri];
yTickLabels = arrayfun(@(X) num2str(round(X*10)/10),deltaBinsVec(1:end-1),'Uni',0);

% yTickLabels = cell(1,numel(deltaDaysBins)-1);
% for bb = 1:numel(deltaDaysBins)-1
%     if deltaDaysBins(bb)<0
%         yTickLabels{bb} = sprintf('%.0f',deltaDaysBins(bb+1));
%     elseif deltaDaysBins(bb)==0
%         yTickLabels{bb} = sprintf('%.0f',deltaDaysBins(bb));
%     else
%         yTickLabels{bb} = sprintf('%.0f',deltaDaysBins(bb));
%     end
%     % if deltaDaysBins(bb+1)<0
%     %     yTickLabels{bb} = sprintf('< %.0f',deltaDaysBins(bb+1));
%     % else
%     %     yTickLabels{bb} = sprintf('> %.0f',deltaDaysBins(bb));
%     % end
% end
% yTickLabels{ismember(yTickLabels,'> -0')} = '0';

UPrestoKeep = cell(1,length(UMFiles));
pSinceAppearance = nan(length(deltaDaysBins)-1,length(UMFiles), length(UIDtoUse));
popCorr_Uni.ISI = nan(length(deltaDaysBins)-1,length(UMFiles), length(UIDtoUse));
popCorr_Uni.natImResp = nan(length(deltaDaysBins)-1,length(UMFiles), length(UIDtoUse));
popCorr_Uni.refPop = nan(length(deltaDaysBins)-1,length(UMFiles), length(UIDtoUse));
popAUC_Uni.ISI = nan(length(deltaDaysBins)-1,length(UMFiles), length(UIDtoUse));
popAUC_Uni.natImResp = nan(length(deltaDaysBins)-1,length(UMFiles), length(UIDtoUse));
popAUC_Uni.refPop = nan(length(deltaDaysBins)-1,length(UMFiles), length(UIDtoUse));
for midx = 1:length(UMFiles)

    %% Load data

    fprintf('Reference %s...\n', UMFiles{midx})

    tmpfile = dir(UMFiles{midx});
    if isempty(tmpfile)
        continue
    end

    %% Get some baseline information:
    fprintf('Loading the data...\n')
    tic
    load(fullfile(tmpfile.folder, tmpfile.name), 'MatchTable', 'UMparam', 'UniqueIDConversion');
    toc

    %% Sort out day situation
    if ~isstruct(UMparam.RawDataPaths{1})
        UMparam.RawDataPaths = cellfun(@(x) dir(x), UMparam.RawDataPaths, 'uni', 0);
    end
    days{midx} = cellfun(@(y) datenum(y), cellfun(@(x) regexp(x.folder,'\\\d*-\d*-\d*\\','match'), UMparam.RawDataPaths, 'uni', 0), 'uni', 0);
    if all(cell2mat(cellfun(@(x) isempty(x), days{midx}, 'uni', 0)))
        days{midx} = cellfun(@(y) datenum(y), cellfun(@(x) regexp(x.folder,'\\\d*-\d*-\d*$','match'), UMparam.RawDataPaths, 'uni', 0), 'uni', 0);
    end
    if all(cell2mat(cellfun(@(x) isempty(x), days{midx}, 'uni', 0)))
        days{midx} = cellfun(@(y) datenum(y{1}(2:end-1)), cellfun(@(x) regexp(x.folder,'\/\d*-\d*-\d*\/','match'), UMparam.RawDataPaths, 'uni', 0), 'uni', 0);
    end
    if all(cell2mat(cellfun(@(x) isempty(x), days{midx}, 'uni', 0)))
        days{midx} = cellfun(@(y) datenum(y{1}(2:end-1)), cellfun(@(x) regexp(x.folder,'\/\d*-\d*-\d*$','match'), UMparam.RawDataPaths, 'uni', 0), 'uni', 0);
    end
    days{midx} = cell2mat(days{midx}) - days{midx}{1};
    deltaDays{midx} = days{midx} - days{midx}';
    deltaDaysUni{midx} = unique(deltaDays{midx});
    deltaDaysUniBins = [deltaDaysUni{midx}-0.5; deltaDaysUni{midx}(end)+0.5];


    %% Check if enough units?
    nclusPerSess = arrayfun(@(X) numel(unique(MatchTable.ID1(MatchTable.RecSes1 == X))),unique(MatchTable.RecSes1));
    if all(nclusPerSess<UMparam.minGoodUnits)
        continue
    end

    % Remove data for days that are having too little neurons
    if any(nclusPerSess<UMparam.minGoodUnits)
        if median(nclusPerSess)<UMparam.minGoodUnits
            continue
        end
    end

    RecSes = RecSes + numel(nclusPerSess);
    nUnitsPerRec{midx} = nclusPerSess;

    %% Extra Positive
    EPosAndNeg(1, midx) = sum(MatchTable.MatchProb(MatchTable.ID1 ~= MatchTable.ID2 & MatchTable.RecSes1 == MatchTable.RecSes2) > UMparam.ProbabilityThreshold)./sum(MatchTable.ID1 ~= MatchTable.ID2 & MatchTable.RecSes1 == MatchTable.RecSes2);

    %Extra Negative
    EPosAndNeg(2, midx) = sum(MatchTable.MatchProb(MatchTable.ID1 == MatchTable.ID2 & MatchTable.RecSes1 == MatchTable.RecSes2) < UMparam.ProbabilityThreshold)./sum(MatchTable.ID1 == MatchTable.ID2 & MatchTable.RecSes1 == MatchTable.RecSes2);

    if UMparam.RunPyKSChronicStitched
        EPosAndNeg(3, midx) = sum(MatchTable.MatchProb(MatchTable.ID1 == MatchTable.ID2 & MatchTable.RecSes1 < MatchTable.RecSes2) < UMparam.ProbabilityThreshold)./sum(MatchTable.ID1 == MatchTable.ID2 & MatchTable.RecSes1 < MatchTable.RecSes2);
    end

    %% Within session noise levels?
    % FalsePositives = nan(3,length(UMFiles));
    % FalseNegatives = nan(3,length(UMFiles));
    FalsePositives(1,midx) = sum(MatchTable.UID1Liberal == MatchTable.UID2Liberal & MatchTable.ID1 ~= MatchTable.ID2 & MatchTable.RecSes1 == MatchTable.RecSes2)./ sum(MatchTable.RecSes1 == MatchTable.RecSes2);
    FalsePositives(2,midx) = sum(MatchTable.UID1 == MatchTable.UID2 & MatchTable.ID1 ~= MatchTable.ID2 & MatchTable.RecSes1 == MatchTable.RecSes2)./ sum(MatchTable.RecSes1 == MatchTable.RecSes2);
    FalsePositives(3,midx) = sum(MatchTable.UID1Conservative == MatchTable.UID2Conservative & MatchTable.ID1 ~= MatchTable.ID2 & MatchTable.RecSes1 == MatchTable.RecSes2)./ sum(MatchTable.RecSes1 == MatchTable.RecSes2);

    FalseNegatives(1,midx) = sum(MatchTable.UID1Liberal ~= MatchTable.UID2Liberal & MatchTable.ID1 == MatchTable.ID2 & MatchTable.RecSes1 == MatchTable.RecSes2)./ sum(MatchTable.RecSes1 == MatchTable.RecSes2);
    FalseNegatives(2,midx) = sum(MatchTable.UID1 ~= MatchTable.UID2 & MatchTable.ID1 == MatchTable.ID2 & MatchTable.RecSes1 == MatchTable.RecSes2)./ sum(MatchTable.RecSes1 == MatchTable.RecSes2);
    FalseNegatives(3,midx) = sum(MatchTable.UID1Conservative ~= MatchTable.UID2Conservative & MatchTable.ID1 == MatchTable.ID2 & MatchTable.RecSes1 == MatchTable.RecSes2)./ sum(MatchTable.RecSes1 == MatchTable.RecSes2);
    for uidtype = 1:numel(UIDtoUse)


        %% For each cluster, find presence and proba of being matched in subsequent recordings
        [UIDuni,indx,~] = unique([MatchTable.(UIDtoUse{uidtype})]);
        RecSes = MatchTable.RecSes1(indx);
        RecSesOpt = unique(RecSes);
        nMatches = sum(MatchTable.(UIDtoUse{uidtype}) == MatchTable.(strrep(UIDtoUse{uidtype},'1','2')) & MatchTable.RecSes2>MatchTable.RecSes1);
        if nMatches < 20*numel(RecSesOpt)
            durationflag = 0;
            for recid = 1:length(UMparam.RawDataPaths)
                meta = ReadMeta2(UMparam.RawDataPaths{recid}.folder);
                Dur = str2num(meta.fileTimeSecs)./60;
                if Dur < UMparam.MinRecordingDuration
                    durationflag = 1;
                end
            end
            if durationflag == 1
                continue
            end
        end


        unitPresence{uidtype,midx} = zeros(numel(days{midx}), numel(UIDuni));
        unitProbaMatch{uidtype,midx} = zeros(numel(deltaDaysUni{midx}), numel(UIDuni));
        for uidx = 1:numel(UIDuni)
            sessUnitPresent = unique(MatchTable(find(MatchTable.(UIDtoUse{uidtype}) == UIDuni(uidx)),:).RecSes1);

            % Get unit presence
            unitPresence{uidtype,midx}(sessUnitPresent,uidx) = 1;

            % Get unit proba of being matched in prev and next days
            tmp = nan(numel(days{midx}),numel(days{midx}));
            tmp(sessUnitPresent,:) = unitPresence{uidtype,midx}(sessUnitPresent,uidx)*unitPresence{uidtype,midx}(:,uidx)';
            tmp(1 + (1+size(tmp,1))*[0:size(tmp,2)-1]) = nan; % remove diagonal
            unitProbaMatch{uidtype,midx}(:,uidx) = histcounts(deltaDays{midx}(tmp == 1),deltaDaysUniBins)./ ...
                histcounts(deltaDays{midx}(ismember(tmp, [0 1])),deltaDaysUniBins);
        end

        % Probability of finding a unit back as a function of minimum number of
        % units available
        UPres = (unitPresence{uidtype,midx}*unitPresence{uidtype,midx}');
        % Compute number of minimum units between each pair of recordings
        nUnits = diag(UPres);
        [nUnits,~] = meshgrid(nUnits);
        % Repeat with ExclFirst
        UPres = UPres./nUnits';
        UPres(logical(eye(size(UPres)))) = nan; % Exclude within
        UPrestoKeep{uidtype} = UPres;
        % Compute probability of a unit returning since it's appearance
        for binid = 1:length(deltaDaysBins)-1
            Idx = deltaDays{midx} > deltaDaysBins(binid) & deltaDays{midx} <= deltaDaysBins(binid+1);
            pSinceAppearance(binid,midx,uidtype) = nanmean(UPres(Idx));
        end


        % Takes way too long! :(
        if computeFuncScores
            tic
            fprintf('Computing probabilities. ')
            %     popProbaMatch{uidtype,midx} = nan(size(deltaDays{midx}));
            popCorr.ISI{midx} = nan(size(deltaDays{midx}));
            popCorr.natImResp{midx} = nan(size(deltaDays{midx}));
            popCorr.refPop{midx} = nan(size(deltaDays{midx}));
            
            popAUC.ISI{midx} = nan(size(deltaDays{midx}));
            popAUC.natImResp{midx} = nan(size(deltaDays{midx}));
            popAUC.refPop{midx} = nan(size(deltaDays{midx}));
            idxMatched = MatchTable.(UIDtoUse{uidtype}) == MatchTable.(strrep(UIDtoUse{uidtype},'1','2'));
            idxNonMatched = MatchTable.(UIDtoUse{uidtype}) ~= MatchTable.(strrep(UIDtoUse{uidtype},'1','2'));
            UID = MatchTable.(UIDtoUse{uidtype});
            for dd1 = 1:numel(days{midx})
                fprintf('Day %d.\n', dd1)
                for dd2 = 1:numel(days{midx})
                    sessIdx = MatchTable.RecSes1 == dd1 & MatchTable.RecSes2 == dd2;
                    unitIdx = sessIdx & ismember(MatchTable.(UIDtoUse{uidtype}), unique(UID(sessIdx))) & idxMatched;
                    unitNonIdx =  sessIdx & ismember(MatchTable.(UIDtoUse{uidtype}), unique(UID(sessIdx))) & idxNonMatched;
                    popCorr.ISI{midx}(dd1,dd2) = nanmedian(MatchTable(unitIdx, :).ISICorr);
                    popCorr.natImResp{midx}(dd1,dd2) = nanmedian(MatchTable(unitIdx, :).natImRespCorr);
                    popCorr.refPop{midx}(dd1,dd2) = nanmedian(MatchTable(unitIdx, :).refPopCorr);
                    % compute AUC
                    if sum(unitIdx)==0
                        continue
                    end
                    labels = [ones(1,sum(unitIdx)) zeros(1,sum(unitNonIdx))];

                    scores = [MatchTable(unitIdx, :).ISICorr; MatchTable(unitNonIdx, :).ISICorr]';
                    if any(~isnan(scores)) && length(unique(labels(~isnan(scores))))>1
                        [~, ~, ~, AUC] = perfcurve(labels(~isnan(scores)), scores(~isnan(scores)), 1);
                        popAUC.ISI{midx}(dd1,dd2) = AUC;
                    end
                    scores = [MatchTable(unitIdx, :).natImRespCorr; MatchTable(unitNonIdx, :).natImRespCorr]';
                    if any(~isnan(scores)) && length(unique(labels(~isnan(scores))))>1
                        [~, ~, ~, AUC] = perfcurve(labels(~isnan(scores)), scores(~isnan(scores)), 1);
                        popAUC.natImResp{midx}(dd1,dd2) = AUC;
                    end

                    scores = [MatchTable(unitIdx, :).refPopCorr; MatchTable(unitNonIdx, :).refPopCorr]';
                    if any(~isnan(scores)) && length(unique(labels(~isnan(scores))))>1
                        [~, ~, ~, AUC] = perfcurve(labels(~isnan(scores)), scores(~isnan(scores)), 1);
                        popAUC.refPop{midx}(dd1,dd2) = AUC;
                    end
                end
            end
            %     popProbaMatch{uidtype,midx}(eye(size(popProbaMatch{uidtype,midx}))==1) = nan; % remove diagonal

            %     popProbaMatch_Uni = nan(1,numel(deltaDaysUni{midx}));
            for binid = 1:length(deltaDaysBins)-1
                Idx = deltaDays{midx} > deltaDaysBins(binid) & deltaDays{midx} <= deltaDaysBins(binid+1);
                popCorr_Uni.ISI(binid,midx,uidtype) = nanmean(popCorr.ISI{midx}(Idx));
                popCorr_Uni.natImResp(binid,midx,uidtype) = nanmean(popCorr.natImResp{midx}(Idx));
                popCorr_Uni.refPop(binid,midx,uidtype) = nanmean(popCorr.refPop{midx}(Idx));
                %         popProbaMatch_Uni(dd) = nanmean(popProbaMatch{uidtype,midx}(Idx));
                popAUC_Uni.ISI(binid,midx,uidtype) = nanmean(popAUC.ISI{midx}(Idx));
                popAUC_Uni.natImResp(binid,midx,uidtype) = nanmean(popAUC.natImResp{midx}(Idx));
                popAUC_Uni.refPop(binid,midx,uidtype) = nanmean(popAUC.refPop{midx}(Idx));

            end
            toc
        end


    end
    %% Plots
    if PlotIndividualMice
        % Units lifetime
        figure('name',UMFiles{midx});
        for uidtype = 1:numel(UIDtoUse)
            subplot(1,3,uidtype)
            imagesc(unitPresence{uidtype,midx}')
            c = colormap('gray'); c = flipud(c); colormap(c)
            caxis([0 1])
            ylabel('Unit')
            xlabel('Days')
            xticks(1:numel(days{midx}));
            xticklabels(num2str(days{midx}'))
            title([UIDtoUse{uidtype}])
        end
        %
        % Probe of matching a unit
        %         figure;
        %         plot(deltaDaysUni{midx},popProbaMatch_Uni,'k')
        %         ylabel('P(match)')
        %         xlabel('Delta days')

        figure('name',UMFiles{midx});
        for uidtype = 1:numel(UIDtoUse)
            subplot(1,3,uidtype)

            [~,sortIdx] = sort(nanmean(unitProbaMatch{uidtype,midx},1),'descend');
            imagesc(deltaDaysUni{midx},1:numel(sortIdx),unitProbaMatch{uidtype,midx}(:,sortIdx)')
            colormap(c)
            hcb = colorbar; hcb.Title.String = 'P(track)';
            caxis([0 1])
            ylabel('Unit')
            xlabel('Delta days')
            title([UIDtoUse{uidtype}])
        end
        for uidtype = 1:numel(UIDtoUse)

            figure('name',[UIDtoUse{uidtype} ' ' UMFiles{midx} ]);

            subplot(2,2,1)
            h=imagesc(UPrestoKeep{uidtype},'AlphaData',~isnan(UPrestoKeep{uidtype}));
            % set(gca,'ydir','normal')
            c = colormap(flipud(gray));
            hcb = colorbar; hcb.Title.String = 'P(track)';
            caxis([0 0.6])
            ylabel('Recording')
            xlabel('Recording')
            makepretty
       

            subplot(2,2,2)
            plot(pSinceAppearance(:,midx,uidtype),'k')
            xlabel('delta Days')
            set(gca,'XTick',1:numel(yTickLabels),'XTickLabel',yTickLabels)
            ylabel('P(track)')
            makepretty

            subplot(2,2,3); hold all
            plot(popCorr_Uni.ISI(:,midx,uidtype),'k')
            plot(popCorr_Uni.natImResp(:,midx,uidtype),'r')
            plot(popCorr_Uni.refPop(:,midx,uidtype),'b')
            xlabel('delta Days')
            ylabel('Functional score')
            set(gca,'XTick',1:numel(deltaDaysBins)-1,'XTickLabel',yTickLabels)
            ylabel('Corr')
            ylim([0 1])
            makepretty

            subplot(2,2,4); hold all
            plot(popAUC_Uni.ISI(:,midx,uidtype),'k')
            plot(popAUC_Uni.natImResp(:,midx,uidtype),'r')
            plot(popAUC_Uni.refPop(:,midx,uidtype),'b')
            xlabel('delta Days')
            ylabel('Functional score')
            set(gca,'XTick',1:numel(deltaDaysBins)-1,'XTickLabel',yTickLabels)
            ylabel('AUC')
            ylim([0 1])
            makepretty
        end

        figure('name',[UMFiles{midx}])
        yyaxis left
        plot(nUnits(1,:),'k-');
        ylabel('total number Units')
        yyaxis right
        plot(days{midx},'r-');
        ylabel('Days')
        makepretty
        offsetAxes

    end
end

%%
UIDCols = [0 0.7 0.2; 0 0 0; 0.7 0.2 0];
IndivFig = figure('name','IndividualDatasets');
AvgFig = figure('name','AverageDatasets');
FunctionalFig = figure('name','Functional');
clear h
for uidtype = 1:numel(UIDtoUse)
    figure(IndivFig)
    subplot(3,3,uidtype)
    hold on
    for midx = 1:length(UMFiles)
        plot(pSinceAppearance(:,midx,uidtype),'color',groupColor(groupVector(midx),:))
    end
    xlabel('delta Days')
    set(gca,'XTick',1:numel(deltaDaysBins)-1,'XTickLabel',yTickLabels)
    ylabel('P(track)')
    title(UIDtoUse{uidtype})
    makepretty
    offsetAxes

    fnames = fieldnames(popCorr_Uni);
    for ff = 1:numel(fnames)
        subplot(3,3,3+ff)
        hold on
        for midx = 1:length(UMFiles)
            plot(popCorr_Uni.(fnames{ff})(:,midx,uidtype),'color',groupColor(groupVector(midx),:))
        end
        xlabel('delta Days')
        set(gca,'XTick',1:numel(deltaDaysBins)-1,'XTickLabel',yTickLabels)
        ylabel('corr')
        title(fnames{ff})
        makepretty
        offsetAxes

        subplot(3,3,6+ff)
        hold on
        for midx = 1:length(UMFiles)
            plot(popAUC_Uni.(fnames{ff})(:,midx,uidtype),'color',groupColor(groupVector(midx),:))
        end
        xlabel('delta Days')
        set(gca,'XTick',1:numel(deltaDaysBins)-1,'XTickLabel',yTickLabels)
        ylabel('AUC')
        title(fnames{ff})
        makepretty
        offsetAxes

    end

    figure(AvgFig)
    nonnanNr = sum(~isnan(pSinceAppearance(:,:,uidtype)),2);
    subplot(2,1,1)
    hold on
    h(uidtype) = errorbar(1:size(pSinceAppearance,1),nanmean(pSinceAppearance(:,:,uidtype),2),nanstd(pSinceAppearance(:,:,uidtype),[],2)./sqrt(nonnanNr-1),'linestyle','-','color',UIDCols(uidtype,:));

    % shadedErrorBar(1:size(pSinceAppearance,1),nanmean(pSinceAppearance,2),nanstd(pSinceAppearance,[],2)./sqrt(nonnanNr-1),'transparent',1)
    set(gca,'XTick',1:numel(deltaDaysBins)-1,'XTickLabel',yTickLabels)
    xlabel('delta Days (>)')
    ylabel('P(track)')
    makepretty
    offsetAxes

    subplot(2,1,2)
    hold on
    plot(1:size(pSinceAppearance(:,:,uidtype),1),nonnanNr,'.','color',UIDCols(uidtype,:))
    ylabel('nr Datasets')
    makepretty
    offsetAxes

   

    figure(FunctionalFig)
    for ff = 1:numel(fnames)
        subplot(2,numel(fnames),ff)
        hold on
        errorbar(1:size(popCorr_Uni.(fnames{ff})(:,:,uidtype),1),nanmean(popCorr_Uni.(fnames{ff})(:,:,uidtype),2),nanstd(popCorr_Uni.(fnames{ff})(:,:,uidtype),[],2)./sqrt(nonnanNr-1),'linestyle','-','color',UIDCols(uidtype,:));
        set(gca,'XTick',1:numel(deltaDaysBins)-1,'XTickLabel',yTickLabels)
        xlabel('delta Days (>)')
        ylabel('corr')
        title(fnames{ff})
        ylim([-0.5 1])
        makepretty
        offsetAxes

        subplot(2,numel(fnames),numel(fnames)+ff)
        hold on
        h(uidtype) = errorbar(1:size(popAUC_Uni.(fnames{ff})(:,:,uidtype),1),nanmean(popAUC_Uni.(fnames{ff})(:,:,uidtype),2),nanstd(popAUC_Uni.(fnames{ff})(:,:,uidtype),[],2)./sqrt(nonnanNr-1),'linestyle','-','color',UIDCols(uidtype,:));
        set(gca,'XTick',1:numel(deltaDaysBins)-1,'XTickLabel',yTickLabels)
        xlabel('delta Days (>)')
        ylabel('AUC')
        ylim([0 1])

        makepretty
        offsetAxes
        
    end   
end
legend([h(:)],UIDtoUse,'Location','best')

if ~exist(savedir)
    mkdir(savedir)
end
saveas(IndivFig,fullfile(savedir,['IndividData.fig']))
saveas(IndivFig,fullfile(savedir,['IndividData.bmp']))
 
saveas(AvgFig,fullfile(savedir,['TrackingProbabilities.fig']))
saveas(AvgFig,fullfile(savedir,['TrackingProbabilities.bmp']))

saveas(FunctionalFig,fullfile(savedir,['FunctionalScores.fig']))
saveas(FunctionalFig,fullfile(savedir,['FunctionalScores.bmp']))

%%
figure('name','FalsePositives')
subplot(1,2,1)
boxplot(FalsePositives'.*100,'boxstyle','filled','extrememode','clip')
set(gca,'XTick',1:3,'XTickLabel',{'Liberal','Intermediate','Conservative'})
makepretty
offsetAxes
ylim([0 10])
ylabel('FP (%)')

%% False Pos and negative
subplot(1,2,2)
scatter(EPosAndNeg(1,:).*100,(EPosAndNeg(2,:)).*100,35,distinguishable_colors(length(UMFiles)),'filled')
xlabel('False Positives (%)')
ylabel('False negative (%)')
title('Using Match Prob')
xlim([0 5])
ylim([0 10])
disp([num2str(nanmedian(EPosAndNeg(1,:).*100)) '+/- ' num2str(mad(EPosAndNeg(1,:).*100)) ' false positives'])
disp([num2str(nanmedian(EPosAndNeg(2,:).*100)) '+/- ' num2str(mad(EPosAndNeg(2,:).*100)) ' false negatives'])
if UMparam.RunPyKSChronicStitched
    disp([num2str(nanmedian(EPosAndNeg(3,:).*100)) '+/- ' num2str(mad(EPosAndNeg(3,:).*100)) ' KS negatives'])
end

%% Info on dataset
DataSetInfo.RecSes = RecSes;
DataSetInfo.nUnits = nUnitsPerRec;
