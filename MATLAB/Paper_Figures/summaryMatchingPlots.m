function [unitPresence, unitProbaMatch, days] = summaryMatchingPlots(UMFiles,groupVector)

%% Settings
PlotIndividualMice = 0;

%% Initialize
if nargin<2 || ~exist('groupVector','var') || isempty(groupVector)
    groupVector = 1:length(UMFiles);
end
groups = unique(groupVector);
groupColor = gray(length(groups)+1);

days = cell(1, length(UMFiles));
deltaDays = cell(1, length(UMFiles));
deltaDaysUni = cell(1, length(UMFiles));
unitPresence = cell(1, length(UMFiles));
unitProbaMatch = cell(1, length(UMFiles));
popProbaMatch = cell(1, length(UMFiles));

deltaDaysBinsOri = [0.1 2 5 10 20 50 100 inf];
deltaDaysBins = [-deltaDaysBinsOri(end:-1:1)-0.1, deltaDaysBinsOri(1:end)];
pSinceAppearance = nan(length(deltaDaysBins)-1,length(UMFiles));
yTickLabels = cell(1,numel(deltaDaysBins)-1);
for bb = 1:numel(deltaDaysBins)-1
    yTickLabels{bb} = sprintf('%.0f< %cdays < %.0f',deltaDaysBins(bb), 916, deltaDaysBins(bb+1));
end

for midx = 1:length(UMFiles)
    %% Load data

    fprintf('Reference %s...\n', UMFiles{midx})

    tmpfile = dir(UMFiles{midx});
    if isempty(tmpfile)
        continue
    end

    fprintf('Loading the data...\n')
    tic
    load(fullfile(tmpfile.folder, tmpfile.name), 'MatchTable', 'UMparam', 'UniqueIDConversion');
    toc
  

    %% Clean up data

    % Remove the splits?

    %% For each cluster, find presence and proba of being matched in subsequent recordings

    UIDuni = unique([MatchTable.UID1]);
    days{midx} = cellfun(@(y) datenum(y), cellfun(@(x) regexp(x.folder,'\\\d*-\d*-\d*\\','match'), UMparam.RawDataPaths, 'uni', 0), 'uni', 0);
    days{midx} = cell2mat(days{midx}) - days{midx}{1};
    deltaDays{midx} = days{midx} - days{midx}';

    unitPresence{midx} = zeros(numel(days{midx}), numel(UIDuni));
    deltaDaysUni{midx} = unique(deltaDays{midx});
    deltaDaysUniBins = [deltaDaysUni{midx}-0.5; deltaDaysUni{midx}(end)+0.5];
    unitProbaMatch{midx} = zeros(numel(deltaDaysUni{midx}), numel(UIDuni));
    for uidx = 1:numel(UIDuni)
        sessUnitPresent = unique(MatchTable(find(MatchTable.UID1 == UIDuni(uidx)),:).RecSes1);

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
    [nUnits,y] = meshgrid(nUnits);
     UPres = UPres./nUnits;
    UPres(logical(eye(size(UPres)))) = nan;
    % Compute probability of a unit returning since it's appearance
    for binid = 1:length(deltaDaysBins)-1
        Idx = deltaDays{midx}>deltaDaysBins(binid)&deltaDays{midx}<=deltaDaysBins(binid+1);
        pSinceAppearance(binid,midx) = nanmean(UPres(Idx));
    end


    % unitProbaMatch will contain the units only once, which means that
    % can't be used to compute the "average" proba.
    % Will have to do it day by day until I find something smarter...
    % Takes way too long! :(
    tic
    fprintf('Computing probabilities.')
    popProbaMatch{midx} = nan(size(deltaDays{midx}));
    idxMatched = MatchTable.UID1 == MatchTable.UID2;
    UID = MatchTable.UID1;
    for dd1 = 1:numel(days{midx})
        for dd2 = 1:numel(days{midx})
            sessIdx = MatchTable.RecSes1 == dd1 & MatchTable.RecSes2 == dd2;
            nMatched = numel(unique(UID(idxMatched & sessIdx))); % number of matched units from day 1
            nTot = numel(unique(UID(sessIdx))); % total number of units on day 1
            popProbaMatch{midx}(dd1,dd2) = nMatched/nTot;
        end
    end
    popProbaMatch{midx}(eye(size(popProbaMatch{midx}))==1) = nan; % remove diagonal

    popProbaMatch_Uni = nan(1,numel(deltaDaysUni{midx}));
    for dd = 1:numel(deltaDaysUni{midx})
        popProbaMatch_Uni(dd) = nanmean(popProbaMatch{midx}(deltaDays{midx} == deltaDaysUni{midx}(dd)));
    end
    toc

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
        figure;
        plot(deltaDaysUni{midx},popProbaMatch_Uni,'k')
        ylabel('P(match)')
        xlabel('Delta days')

        figure;
        [~,sortIdx] = sort(nanmean(unitProbaMatch{midx},1),'descend');
        imagesc(deltaDaysUni{midx},1:numel(sortIdx),unitProbaMatch{midx}(:,sortIdx)')
        colormap(c)
        hcb = colorbar; hcb.Title.String = 'P(match)';
        caxis([0 1])
        ylabel('Unit')
        xlabel('Delta days')

        figure;
        subplot(1,2,1)
        imagesc(UPres)
        set(gca,'XTick',1:length(days{midx}),'XTickLabel',days{midx},'yTick',1:length(days{midx}),'yTickLabel',days{midx})
        c = colormap('gray'); c = flipud(c); colormap(c)
        hcb = colorbar; hcb.Title.String = 'P(match)';
        caxis([0 1])
        ylabel('Delta days')
        xlabel('Delta days')
        makepretty

        subplot(1,2,2)
        plot(pSinceAppearance(:,midx),'k')
        xlabel('delta Days')
        ylabel('P(match)')
        set(gca,'XTick',1:numel(deltaDaysBins)-1,'XTickLabel',yTickLabels)
        ylabel('P(match)')
        makepretty

    end
end

figure;
subplot(1,2,1)
hold on
for midx = 1:length(UMFiles)
    plot(pSinceAppearance(:,midx),'color',groupColor(groupVector(midx),:))
end
xlabel('delta Days')
set(gca,'XTick',1:numel(deltaDaysBins)-1,'XTickLabel',yTickLabels)
ylabel('P(match)')
makepretty
offsetAxes

subplot(1,2,2)
hold on
shadedErrorBar(nanmean(pSinceAppearance,2),nanstd(pSinceAppearance,[],2),'transparent',1)
set(gca,'XTick',1:numel(deltaDaysBins)-1,'XTickLabel',yTickLabels)
xlabel('delta Days')
ylabel('P(match)')
makepretty
offsetAxes


%% Summary plots
probaBinned = nan(numel(deltaDaysBins)-1, length(UMFiles));
for midx = 1:length(UMFiles)
    for bb = 1:numel(deltaDaysBins)-1
        idx = deltaDays{midx} > deltaDaysBins(bb) & deltaDays{midx} <= deltaDaysBins(bb+1);
        if any(idx(:))
            probaBinned(bb,midx) = nanmean(popProbaMatch{midx}(idx));
        end
    end
end
figure;
hold all
x = (1:numel(deltaDaysBins)-1)';
y = nanmean(probaBinned,2);
err = 2*nanstd(probaBinned,[],2)./sqrt(sum(~isnan(probaBinned),2));
plot(x,y,'k')
patch([x; flip(x)], [y-err; flip(y+err)], 'k', 'FaceAlpha',0.25, 'EdgeColor','none')
xticks(1:numel(deltaDaysBins)-1)
yTickLabels = cell(1,numel(deltaDaysBins)-1);
for bb = 1:numel(deltaDaysBins)-1
    yTickLabels{bb} = sprintf('%.0f< %cdays < %.0f',deltaDaysBins(bb), 916, deltaDaysBins(bb+1));
end
xticklabels(yTickLabels)
ylabel('P(match)')