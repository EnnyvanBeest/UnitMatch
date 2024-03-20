function [unitPresence, unitProbaMatch, days, FalsePositives] = summaryMatchingPlots(UMFiles,groupVector,computeFuncScores)

%% Settings
PlotIndividualMice = 0;
% Drifts2Exclude = 50; % more than 50 micron, expected to lose neurons. This has nothing to do with software but with quality of recordings
MinNumsUnit = -inf; % should have at least this amount of neurons, otherwise not a good session
savedir = UMFiles{1}(1:find(UMFiles{1}~=UMFiles{end},1,'first')-1);

%% Initialize
if nargin<2 || ~exist('groupVector','var') || isempty(groupVector)
    groupVector = 1:length(UMFiles);
end

if nargin<3
    computeFuncScores = 0;
end
groups = unique(groupVector);
groupColor = distinguishable_colors(max(groups)+1);

days = cell(1, length(UMFiles));
deltaDays = cell(1, length(UMFiles));
deltaDaysUni = cell(1, length(UMFiles));
unitPresence = cell(1, length(UMFiles));
unitProbaMatch = cell(1, length(UMFiles));
popProbaMatch = cell(1, length(UMFiles));

% WIthin session check
FalsePositives = nan(3,length(UMFiles)); % Liberal/intermedi/conserv
FalseNegatives = nan(3,length(UMFiles));

% deltaDaysBinsOri = [0.1 1 2 3 4 5 10 20 50 100 inf];
% deltaDaysBinsOri = [0.1 1 2 3 4 5 10 15 20 30 40 50 75 100 125 150 inf];
deltaDaysBinsOri = 2.^(-1:8);
deltaDaysBins = [-deltaDaysBinsOri(end:-1:1)-0.1,deltaDaysBinsOri(1:end)];

deltaBinsVec =[-deltaDaysBinsOri(end:-1:1),deltaDaysBinsOri];
pSinceAppearance = nan(length(deltaDaysBins)-1,length(UMFiles));
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

popCorr_Uni.ISI = nan(length(deltaDaysBins)-1,length(UMFiles));
popCorr_Uni.natImResp = nan(length(deltaDaysBins)-1,length(UMFiles));
popCorr_Uni.refPop = nan(length(deltaDaysBins)-1,length(UMFiles));

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

    %% Within session noise levels?
    % FalsePositives = nan(3,length(UMFiles));
    % FalseNegatives = nan(3,length(UMFiles));
    FalsePositives(1,midx) = sum(MatchTable.UID1Liberal == MatchTable.UID2Liberal & MatchTable.ID1 ~= MatchTable.ID2 & MatchTable.RecSes1 == MatchTable.RecSes2)./ sum(MatchTable.RecSes1 == MatchTable.RecSes2);
    FalsePositives(2,midx) = sum(MatchTable.UID1 == MatchTable.UID2 & MatchTable.ID1 ~= MatchTable.ID2 & MatchTable.RecSes1 == MatchTable.RecSes2)./ sum(MatchTable.RecSes1 == MatchTable.RecSes2);
    FalsePositives(3,midx) = sum(MatchTable.UID1Conservative == MatchTable.UID2Conservative & MatchTable.ID1 ~= MatchTable.ID2 & MatchTable.RecSes1 == MatchTable.RecSes2)./ sum(MatchTable.RecSes1 == MatchTable.RecSes2);

    FalseNegatives(1,midx) = sum(MatchTable.UID1Liberal ~= MatchTable.UID2Liberal & MatchTable.ID1 == MatchTable.ID2 & MatchTable.RecSes1 == MatchTable.RecSes2)./ sum(MatchTable.RecSes1 == MatchTable.RecSes2);
    FalseNegatives(2,midx) = sum(MatchTable.UID1 ~= MatchTable.UID2 & MatchTable.ID1 == MatchTable.ID2 & MatchTable.RecSes1 == MatchTable.RecSes2)./ sum(MatchTable.RecSes1 == MatchTable.RecSes2);
    FalseNegatives(3,midx) = sum(MatchTable.UID1Conservative ~= MatchTable.UID2Conservative & MatchTable.ID1 == MatchTable.ID2 & MatchTable.RecSes1 == MatchTable.RecSes2)./ sum(MatchTable.RecSes1 == MatchTable.RecSes2);


    % 
    % if any(nansum(UMparam.drift,3)>Drifts2Exclude)
    %     keyboard
    % end

    % Remove the splits?
    %   UMparam = load(fullfile(tmpfile(id).folder,tmpfile(id).name),'UMparam','UniqueIDConversion');
    %             UniqueIDConversion = UMparam.UniqueIDConversion;
    %             UMparam = UMparam.UMparam;
    % 
    %             % Check for enough good units per session
    %             [nums,id1,id2] = unique(UniqueIDConversion.recsesAll(logical(UniqueIDConversion.GoodID)));
    %             numsperSes = diff([id1; length(id2)]);
    %             if any(numsperSes<MinNumsUnit)
    %                  continue
    %             end
    % 
    %             if any(nansum(UMparam.drift,3)>Drifts2Exclude)
    %                 continue
    %             end
    nclusPerSess = arrayfun(@(X) numel(unique(MatchTable.ID1(MatchTable.RecSes1 == X))),unique(MatchTable.RecSes1));
    if all(nclusPerSess<UMparam.minGoodUnits)
        continue
    end
    
  

    %% For each cluster, find presence and proba of being matched in subsequent recordings
    UDtoUse = 'UID1';

    [UIDuni,indx,~] = unique([MatchTable.(UDtoUse)]);
    RecSes = MatchTable.RecSes1(indx);
    RecSesOpt = unique(RecSes);
    nMatches = sum(MatchTable.UID1 == MatchTable.UID2 & MatchTable.RecSes2>MatchTable.RecSes1);
    if nMatches < 20*numel(RecSesOpt)
        durationflag = 0;
        for recid = 1:length(UMparam.AllRawPaths)
            meta = ReadMeta2(UMparam.AllRawPaths{recid}.folder)
            Dur = str2num(meta.fileTimeSecs)./60;
            if Dur<UMparam.MinRecordingDuration
                durationflag = 1;
            end
        end
        if durationflag == 1
            continue
        end
    end

    days{midx} = cellfun(@(y) datenum(y), cellfun(@(x) regexp(x.folder,'\\\d*-\d*-\d*\\','match'), UMparam.RawDataPaths, 'uni', 0), 'uni', 0);
    days{midx} = cell2mat(days{midx}) - days{midx}{1};
    deltaDays{midx} = days{midx} - days{midx}';

    % Remove data for days that are having too little neurons
    if any(nclusPerSess<UMparam.minGoodUnits)
        if median(nclusPerSess)<UMparam.minGoodUnits
            continue    
        end
    end

    unitPresence{midx} = zeros(numel(days{midx}), numel(UIDuni));
    deltaDaysUni{midx} = unique(deltaDays{midx});
    deltaDaysUniBins = [deltaDaysUni{midx}-0.5; deltaDaysUni{midx}(end)+0.5];
    unitProbaMatch{midx} = zeros(numel(deltaDaysUni{midx}), numel(UIDuni));
    for uidx = 1:numel(UIDuni)
        sessUnitPresent = unique(MatchTable(find(MatchTable.(UDtoUse) == UIDuni(uidx)),:).RecSes1);

        % Get unit presence
        unitPresence{midx}(sessUnitPresent,uidx) = 1;

        % Get unit proba of being matched in prev and next days
        tmp = nan(numel(days{midx}),numel(days{midx}));
        tmp(sessUnitPresent,:) = unitPresence{midx}(sessUnitPresent,uidx)*unitPresence{midx}(:,uidx)';
        tmp(1 + (1+size(tmp,1))*[0:size(tmp,2)-1]) = nan; % remove diagonal
        unitProbaMatch{midx}(:,uidx) = histcounts(deltaDays{midx}(tmp == 1),deltaDaysUniBins)./ ...
            histcounts(deltaDays{midx}(ismember(tmp, [0 1])),deltaDaysUniBins);
    end

    % Probability of finding a unit back as a function of minimum number of
    % units available
    UPres = (unitPresence{midx}*unitPresence{midx}');
    % Compute number of minimum units between each pair of recordings
    nUnits = diag(UPres);      
    [nUnits,~] = meshgrid(nUnits);
    % Repeat with ExclFirst   
    UPres = UPres./nUnits';
    UPres(logical(eye(size(UPres)))) = nan; % Exclude within
    % Compute probability of a unit returning since it's appearance
    for binid = 1:length(deltaDaysBins)-1
        Idx = deltaDays{midx} > deltaDaysBins(binid) & deltaDays{midx} <= deltaDaysBins(binid+1);
        pSinceAppearance(binid,midx) = nanmean(UPres(Idx));
    end


    % Takes way too long! :(
    if computeFuncScores
        tic
        fprintf('Computing probabilities. ')
        %     popProbaMatch{midx} = nan(size(deltaDays{midx}));
        popCorr.ISI{midx} = nan(size(deltaDays{midx}));
        popCorr.natImResp{midx} = nan(size(deltaDays{midx}));
        popCorr.refPop{midx} = nan(size(deltaDays{midx}));
        idxMatched = MatchTable.UID1 == MatchTable.UID2;
        UID = MatchTable.UID1;
        for dd1 = 1:numel(days{midx})
            fprintf('Day %d.\n', dd1)
            for dd2 = 1:numel(days{midx})
                sessIdx = MatchTable.RecSes1 == dd1 & MatchTable.RecSes2 == dd2;
                unitIdx = sessIdx & ismember(MatchTable.(UDtoUse), unique(UID(sessIdx)))& idxMatched;
                popCorr.ISI{midx}(dd1,dd2) = nanmedian(MatchTable(unitIdx, :).ISICorr);
                popCorr.natImResp{midx}(dd1,dd2) = nanmedian(MatchTable(unitIdx, :).natImRespCorr);
                popCorr.refPop{midx}(dd1,dd2) = nanmedian(MatchTable(unitIdx, :).refPopCorr);
                %             nMatched = numel(unique(UID(idxMatched & sessIdx))); % number of matched units from day 1
                %             nTot = numel(unique(UID(sessIdx))); % total number of units on day 1
                %             popProbaMatch{midx}(dd1,dd2) = nMatched/nTot;
            end
        end
        %     popProbaMatch{midx}(eye(size(popProbaMatch{midx}))==1) = nan; % remove diagonal

        %     popProbaMatch_Uni = nan(1,numel(deltaDaysUni{midx}));
        for binid = 1:length(deltaDaysBins)-1
            Idx = deltaDays{midx} > deltaDaysBins(binid) & deltaDays{midx} <= deltaDaysBins(binid+1);
            popCorr_Uni.ISI(binid,midx) = nanmean(popCorr.ISI{midx}(Idx));
            popCorr_Uni.natImResp(binid,midx) = nanmean(popCorr.natImResp{midx}(Idx));
            popCorr_Uni.refPop(binid,midx) = nanmean(popCorr.refPop{midx}(Idx));
            %         popProbaMatch_Uni(dd) = nanmean(popProbaMatch{midx}(Idx));
        end
        toc
    end

    %% Plots
    if PlotIndividualMice
        % Units lifetime
        figure;
        imagesc(unitPresence{midx}')
        c = colormap('gray'); c = flipud(c); colormap(c)
        caxis([0 1])
        ylabel('Unit')
        xlabel('Days')
        xticks(1:numel(days{midx}));
        xticklabels(num2str(days{midx}'))
        % 
        % Probe of matching a unit
%         figure;
%         plot(deltaDaysUni{midx},popProbaMatch_Uni,'k')
%         ylabel('P(match)')
%         xlabel('Delta days')

        figure;
        [~,sortIdx] = sort(nanmean(unitProbaMatch{midx},1),'descend');
        imagesc(deltaDaysUni{midx},1:numel(sortIdx),unitProbaMatch{midx}(:,sortIdx)')
        colormap(c)
        hcb = colorbar; hcb.Title.String = 'P(track)';
        caxis([0 1])
        ylabel('Unit')
        xlabel('Delta days')

        figure;
        subplot(2,2,1)
        h=imagesc(UPres,'AlphaData',~isnan(UPres));
        % set(gca,'ydir','normal')
        c = colormap(flipud(gray)); 
        hcb = colorbar; hcb.Title.String = 'P(track)';
        caxis([0 0.6])
        ylabel('Recording')
        xlabel('Recording')
        makepretty

        subplot(2,2,3)
        yyaxis left
        plot(nUnits(1,:),'k-');
        ylabel('total number Units')
        yyaxis right
        plot(days{midx},'r-');
        ylabel('Days')
        makepretty
        offsetAxes
  
        subplot(2,2,2)
        plot(pSinceAppearance(:,midx),'k')
        xlabel('delta Days')
        ylabel('P(track)')
        set(gca,'XTick',1:numel(yTickLabels),'XTickLabel',yTickLabels)
        ylabel('P(match)')
        makepretty

        subplot(2,2,4); hold all
        plot(popCorr_Uni.ISI(:,midx),'k')
        plot(popCorr_Uni.natImResp(:,midx),'r')
        plot(popCorr_Uni.refPop(:,midx),'b')
        xlabel('delta Days')
        ylabel('Functional score')
        set(gca,'XTick',1:numel(deltaDaysBins)-1,'XTickLabel',yTickLabels)
        ylabel('P(track)')
        makepretty

    end
end

%%
figure;
subplot(3,2,[1,3,5])
hold on
for midx = 1:length(UMFiles)
    plot(pSinceAppearance(:,midx),'color',groupColor(groupVector(midx),:))
end
xlabel('delta Days')
set(gca,'XTick',1:numel(deltaDaysBins)-1,'XTickLabel',yTickLabels)
ylabel('P(track)')
makepretty
offsetAxes

nonnanNr = sum(~isnan(pSinceAppearance),2);
subplot(3,2,[4,6])
hold on
errorbar(1:size(pSinceAppearance,1),nanmean(pSinceAppearance,2),nanstd(pSinceAppearance,[],2)./sqrt(nonnanNr-1),'linestyle','-')

% shadedErrorBar(1:size(pSinceAppearance,1),nanmean(pSinceAppearance,2),nanstd(pSinceAppearance,[],2)./sqrt(nonnanNr-1),'transparent',1)
set(gca,'XTick',1:numel(deltaDaysBins)-1,'XTickLabel',yTickLabels)
xlabel('delta Days (>)')
ylabel('P(track)')
makepretty
offsetAxes


subplot(3,2,[2])
plot(1:size(pSinceAppearance,1),nonnanNr,'.')
makepretty
offsetAxes

saveas(gcf,fullfile(savedir,[UDtoUse '_TrackingProbabilities.fig']))
saveas(gcf,fullfile(savedir,[UDtoUse '_TrackingProbabilities.bmp']))
figure
fnames = fieldnames(popCorr_Uni);
for ff = 1:numel(fnames)
    subplot(1,numel(fnames)+1,ff)
    hold on
    for midx = 1:length(UMFiles)
        plot(popCorr_Uni.(fnames{ff})(:,midx),'color',groupColor(groupVector(midx),:))
    end
    xlabel('delta Days')
    set(gca,'XTick',1:numel(deltaDaysBins)-1,'XTickLabel',yTickLabels)
    ylabel('Functional score')
    title(fnames{ff})
    makepretty
    offsetAxes
end

subplot(1,numel(fnames)+1,numel(fnames)+1)
hold on
h(1) = errorbar(1:size(popCorr_Uni.ISI,1),nanmean(popCorr_Uni.ISI,2),nanstd(popCorr_Uni.ISI,[],2)./sqrt(nonnanNr-1),'linestyle','-','markerfacecolor',[1 0 0]);
h(2) = errorbar(1:size(popCorr_Uni.natImResp,1),nanmean(popCorr_Uni.natImResp,2),nanstd(popCorr_Uni.natImResp,[],2)./sqrt(nonnanNr-1),'linestyle','-','markerfacecolor',[0 1 0]);
h(3) = errorbar(1:size(popCorr_Uni.refPop,1),nanmean(popCorr_Uni.refPop,2),nanstd(popCorr_Uni.refPop,[],2)./sqrt(nonnanNr-1),'linestyle','-','markerfacecolor',[0 0 1]);
set(gca,'XTick',1:numel(deltaDaysBins)-1,'XTickLabel',yTickLabels)
xlabel('delta Days (>)')
ylabel('Functional score')
legend([h(:)],{'ISI','NatIm','RefPop'})
makepretty
offsetAxes
saveas(gcf,fullfile(savedir,[UDtoUse '_FunctionalScores.fig']))
saveas(gcf,fullfile(savedir,[UDtoUse '_FunctionalScores.bmp']))

%% 
figure('name','FalsePositives')
boxplot(FalsePositives'.*100,'boxstyle','filled','extrememode','clip')
set(gca,'XTick',1:3,'XTickLabel',{'Liberal','Intermediate','Conservative'})
makepretty
offsetAxes
ylim([0 10])
ylabel('FP (%)')
% h=barwitherr(nanstd(FalsePositives,[],2),nanmean(FalsePositives,2));
% %% Summary plots
% probaBinned = nan(numel(deltaDaysBins)-1, length(UMFiles));
% for midx = 1:length(UMFiles)
%     for bb = 1:numel(deltaDaysBins)-1
%         idx = deltaDays{midx} > deltaDaysBins(bb) & deltaDays{midx} <= deltaDaysBins(bb+1);
%         if any(idx(:))
%             probaBinned(bb,midx) = nanmean(popProbaMatch{midx}(idx));
%         end
%     end
% end
% figure;
% hold all
% x = (1:numel(deltaDaysBins)-1)';
% y = nanmean(probaBinned,2);
% err = 2*nanstd(probaBinned,[],2)./sqrt(sum(~isnan(probaBinned),2));
% plot(x,y,'k')
% plot(x,nanmean(pSinceAppearance,2),'k')
% patch([x; flip(x)], [y-err; flip(y+err)], 'k', 'FaceAlpha',0.25, 'EdgeColor','none')
% xticks(1:numel(deltaDaysBins)-1)
% yTickLabels = cell(1,numel(deltaDaysBins)-1);
% for bb = 1:numel(deltaDaysBins)-1
%     yTickLabels{bb} = sprintf('%.0f< %cdays < %.0f',deltaDaysBins(bb), 916, deltaDaysBins(bb+1));
% end
% xticklabels(yTickLabels)
% ylabel('P(match)')