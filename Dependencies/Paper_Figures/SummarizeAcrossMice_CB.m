function summaryPlots(UMFiles, TakeRank, groupVector)

%% Define parameters

% Initialize
TakeRank = 0; %if 0 , take cross-correlation scores (these may be less informative than rank)
if TakeRank
    fprintf("Taking the rank!\n")
    FPNames = {'FRRankScore','ACGRankScore','RankScore','NImgRankScore'};
    stepsz = [1 1 1 1];
    minVal = [1 1 1 1];
    maxVal = [21 21 21 21];
    flipROC = [0 0 0 0];
else
    fprintf("Taking the correlation values!\n")
    FPNames = {'FRDiff','ACGCorr','FingerprintCor','NatImCorr'};
    stepsz = [0.1 0.1 0.1 0.1];
    minVal = [0 -1 -1 -1];
    maxVal = [15 1 1 1];
    flipROC = [0 1 1 1];
end

histBins = cell(1,numel(FPNames));
histBinsCenter = cell(1,numel(FPNames));
for fpIdx = 1:numel(FPNames)
    histBins{fpIdx} = minVal(fpIdx):stepsz(fpIdx):maxVal(fpIdx);
    histBinsCenter{fpIdx} = histBins{fpIdx}(1:end-1) + diff(histBins{fpIdx})/2;
end
UseKSLabels = PrepareClusInfoparams.RunPyKSChronicStitched;
ROCBins = 0:0.01:1;
minMatches = 20;

groups = unique(groupVector);
groupColor = gray(length(groups)+1);

%% Loop over mice to get all Distributions / ROCs / AUCs

clear FPSum
deltaDays = cell(1, length(UMFiles));
numMatchedUnits = cell(1, length(UMFiles));
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

    sessIDs = unique(MatchTable.RecSes1);

    % Initialize
    for fpIdx = 1:numel(FPNames)
        FPNameCurr = FPNames{fpIdx};
        FPSum.(FPNameCurr).Distr{midx} = nan(numel(histBins{fpIdx})-1, 3, numel(sessIDs)-1, numel(sessIDs));
        FPSum.(FPNameCurr).AUC{midx} = nan(3, numel(sessIDs)-1, numel(sessIDs));
        FPSum.(FPNameCurr).ROC{midx} = nan(numel(ROCBins), 3, numel(sessIDs)-1, numel(sessIDs));
    end

    %%% HACK -- FIGURE OUT A HOMOGENEIZED VAR TYPE
    if ~iscell(UMparam.AllRawPaths)
        for ii = 1:numel(UMparam.AllRawPaths)
            tmp{ii} = UMparam.AllRawPaths(ii);
        end
        UMparam.AllRawPaths = tmp;
    end

    %% Loop through pairs of sessions

    days = cellfun(@(y) datenum(y), cellfun(@(x) regexp(x.folder,'\\\d*-\d*-\d*\\','match'), UMparam.AllRawPaths, 'uni', 0), 'uni', 0);
    days = cell2mat(days) - days{1};
    deltaDays{midx} = nan(numel(sessIDs)-1,numel(sessIDs));
    numMatchedUnits{midx} = nan(numel(sessIDs)-1,numel(sessIDs));
    for sess1Idx = 1:numel(sessIDs)-1

        sess1 = sessIDs(sess1Idx);
        day1 = regexp(UMparam.AllRawPaths{sess1}.folder,'\d*-\d*-\d*','match'); day1 = datenum(day1{1});

        for sess2Idx = sess1Idx+1:numel(sessIDs)

            sess2 = sessIDs(sess2Idx);
            day2 = regexp(UMparam.AllRawPaths{sess2}.folder,'\d*-\d*-\d*','match'); day2 = datenum(day2{1});
            deltaDays{midx}(sess1Idx,sess2Idx) = abs(day2 - day1);

            %% Cut table to specific days

            MatchTable_pair = MatchTable(ismember(MatchTable.RecSes1, [sess1 sess2]) & ismember(MatchTable.RecSes2, [sess1 sess2]), :);

            %% Number of matches

            if ~UseKSLabels
                numMatchedUnits{midx}(sess1Idx,sess2Idx) = sum((MatchTable_pair.UID1 == MatchTable_pair.UID2) & ...
                    (MatchTable_pair.RecSes1 ~= MatchTable_pair.RecSes2))/2;
            else
                numMatchedUnits{midx}(sess1Idx,sess2Idx) = sum((MatchTable_pair.ID1 == MatchTable_pair.ID2) & ...
                    (MatchTable_pair.RecSes1 ~= MatchTable_pair.RecSes2))/2;
            end

            %% Looping through fingerprints
            for fpIdx = 1:numel(FPNames)
                FPNameCurr = FPNames{fpIdx};

                if ~any(ismember(MatchTable_pair.Properties.VariableNames, FPNameCurr))
                    fprintf('%s has no %s in table...', UMFiles{midx}, FPNameCurr)
                    continue
                end

                %% select pairs within and across days

                % Special criteria
                if contains(FPNameCurr,'NatI') || contains(FPNameCurr,'NImg') %%% clean
                    % Select pairs of which the first neuron has good
                    % test-retest reliability
                    Unit1ID = double(MatchTable_pair.ID1) + 10e9*MatchTable_pair.RecSes1; % identifier for each pair
                    uUnit1ID = unique(Unit1ID); % list of units
                    reliability = MatchTable_pair((MatchTable_pair.ID1 == MatchTable_pair.ID2) & (MatchTable_pair.RecSes1 == MatchTable_pair.RecSes2),:).NatImCorr; % test-retest reliability of each unit

                    validPairs = ismember(Unit1ID, uUnit1ID(reliability > 0.2));
                else
                    validPairs = ones(size(MatchTable_pair,1),1);
                end

                % Extract groups: "within", "match", "non-match"
                if ~UseKSLabels
                    WithinIdx = find((MatchTable_pair.UID1 == MatchTable_pair.UID2) & (MatchTable_pair.RecSes1 == MatchTable_pair.RecSes2) & validPairs); %Within session, same unit (cross-validation)
                    MatchIdx = find((MatchTable_pair.UID1 == MatchTable_pair.UID2) & (MatchTable_pair.RecSes1 ~= MatchTable_pair.RecSes2) & validPairs); %Across session, same unit (cross-validation)
                    NonMatchIdx = find((MatchTable_pair.UID1 ~= MatchTable_pair.UID2) & validPairs); % Not the same unit
                else
                    WithinIdx = find((MatchTable_pair.ID1 == MatchTable_pair.ID2) & (MatchTable_pair.RecSes1 == MatchTable_pair.RecSes2) & validPairs); %Within session, same unit (cross-validation)
                    MatchIdx = find((MatchTable_pair.ID1 == MatchTable_pair.ID2) & (MatchTable_pair.RecSes1 ~= MatchTable_pair.RecSes2) & validPairs); %Across session, same unit (cross-validation)
                    NonMatchIdx = find((MatchTable_pair.ID1 ~= MatchTable_pair.ID2) & validPairs); % Not the same unit
                end

                %% Check that passed the criteria

                % Condition to go ahead
                goAhead = numel(MatchIdx)/2 >= minMatches && ... % minimum number of matches
                    ~all(isnan(MatchTable_pair.(FPNameCurr)(MatchIdx))); % not all nans

                if goAhead
                    %% Compute fingerprint correlations and AUCs

                    FPCorr = MatchTable_pair.(FPNameCurr);

                    % Compute distributions
                    hw = histcounts(FPCorr(WithinIdx), histBins{fpIdx}) ./ length(WithinIdx);
                    hm = histcounts(FPCorr(MatchIdx), histBins{fpIdx}) ./ length(MatchIdx);
                    hn = histcounts(FPCorr(NonMatchIdx), histBins{fpIdx}) ./ length(NonMatchIdx);
                    % Save
                    FPSum.(FPNameCurr).Distr{midx}(:,:,sess1Idx,sess2Idx) = [hw', hm', hn'];

                    % Compute ROCs and AUCs
                    [ROC1, AUC1] = getAUC(FPCorr,MatchIdx,NonMatchIdx,ROCBins,flipROC(fpIdx));
                    [ROC2, AUC2] = getAUC(FPCorr,WithinIdx,MatchIdx,ROCBins,flipROC(fpIdx));
                    [ROC3, AUC3] = getAUC(FPCorr,WithinIdx,NonMatchIdx,ROCBins,flipROC(fpIdx));
                    % Save
                    FPSum.(FPNameCurr).ROC{midx}(:,:,sess1Idx,sess2Idx) = [ROC1,ROC2,ROC3];
                    FPSum.(FPNameCurr).AUC{midx}(:,sess1Idx,sess2Idx) = [AUC1, AUC2, AUC3];
                end
            end
        end
    end
end

%% Figure -- example mouse

% Best example mouse
% midx = find(cellfun(@(X) nansum(X(:)),numMatchedUnits) == max(cellfun(@(X)  nansum(X(:)),numMatchedUnits)),1,'first');

% Number of matched units (matrix)
for midx = 1:numel(UMFiles)
    figure('Position', [80 700 1700 170],'Name', fileparts(fileparts(UMFiles{midx})));
    s(1) = subplot(1,numel(FPNames)+2,1);
    imagesc(numMatchedUnits{midx})
    xticks(1:size(deltaDays{midx},2))
    xticklabels(days)
    yticks(1:size(deltaDays{midx},1))
    yticklabels(days(1:end-1))
    axis equal tight
    colorbar
    c = colormap("gray"); colormap(s(1), flipud(c));
    subplot(1,numel(FPNames)+2,2)
    scatter(deltaDays{midx}(:), mat2vec(numMatchedUnits{midx}),20,groupColor(groupVector(midx),:),'filled')
    ylabel('Number of matches')
    xlabel('\Deltadays')

    for fpIdx = 1:numel(FPNames)
        FPNameCurr = FPNames{fpIdx};
        s(fpIdx+2) = subplot(1,numel(FPNames)+2,fpIdx+2);
        imagesc(squeeze(FPSum.(FPNameCurr).AUC{midx}(1,:,:)))
        xticks(1:size(deltaDays{midx},2))
        xticklabels(days)
        yticks(1:size(deltaDays{midx},1))
        yticklabels(days(1:end-1))
        axis equal tight
        colormap( s(fpIdx+2), "RedBlue")
        clim([0 1])
        title(sprintf('Fingerprint %s', FPNameCurr))
    end
end

%% Figure -- all mice

% Build matrices across mice -- average in groups according to groupVector
distMatrix = cell(1,numel(FPNames));
ROCMatrix = cell(1,numel(FPNames));
AUCMatrix = cell(1,numel(FPNames));
for fpIdx = 1:numel(FPNames)
    FPNameCurr = FPNames{fpIdx};
    distMatrix{fpIdx} = nan(numel(histBinsCenter{fpIdx}), 3, length(groups));
    ROCMatrix{fpIdx} = nan(numel(ROCBins), 3, length(groups));
    AUCMatrix{fpIdx} = nan(3, length(groups));
    for gg = 1:length(groups)
        % Distributions
        tmp = cellfun(@(x) x(:,:,:), FPSum.(FPNameCurr).Distr(groupVector == gg), 'uni', 0);
        distMatrix{fpIdx}(:,:,gg) = nanmean(cat(3,tmp{:}),3);

        % ROC
        tmp = cellfun(@(x) x(:,:,:), FPSum.(FPNameCurr).ROC(groupVector == gg), 'uni', 0);
        ROCMatrix{fpIdx}(:,:,gg) = nanmean(cat(3,tmp{:}),3);

        % AUC
        tmp = cellfun(@(x) x(:,:), FPSum.(FPNameCurr).AUC(groupVector == gg), 'uni', 0);
        AUCMatrix{fpIdx}(:,gg) = nanmean(cat(2,tmp{:}),2);
    end
end

% Plot
distrCols = [0 0.7 0; 1 0 0; 0 0 0.7]; % Within / Match / Non-match
ROCCols = [1 0 0; 0.5 0.5 0.5; 0 0 0]; % Match vs Non-match / Within vs Match / Within vs Non-match
figure('Position', [400 270 800 700]);
for fpIdx = 1:numel(FPNames)
    FPNameCurr = FPNames{fpIdx};

    % Plot distribution
    subplot(4,numel(FPNames),0*numel(FPNames)+fpIdx); hold all
    for hid = 1:3
        h(hid) = shadedErrorBar(histBinsCenter{fpIdx}, nanmean(distMatrix{fpIdx}(:,hid,:),3), ...
            nanstd(distMatrix{fpIdx}(:,hid,:),[],3)./sqrt(sum(~isnan(distMatrix{fpIdx}(:,hid,:)),3)));
        h(hid).mainLine.Color = distrCols(hid,:);
        h(hid).patch.FaceColor = distrCols(hid,:);
        h(hid).edge(1).Color = 'none';
        h(hid).edge(2).Color = 'none';
    end
    title(sprintf('%s', FPNameCurr))
    if TakeRank; xlabel('Rank'); else; xlabel('Correlation'); end
    if TakeRank; xticks([1 10 maxVal]); xticklabels({'1','10',sprintf('>%d',maxVal-1)}); end
    ylabel('Proportion')
    offsetAxes
    makepretty

    % Plot ROC
    subplot(4,numel(FPNames),1*numel(FPNames)+fpIdx); hold all
    for hid = 1:3
        h(hid) = shadedErrorBar(ROCBins, nanmean(ROCMatrix{fpIdx}(:,hid,:),3), ...
            nanstd(ROCMatrix{fpIdx}(:,hid,:),[],3)./sqrt(sum(~isnan(ROCMatrix{fpIdx}(:,hid,:)),3)));
        h(hid).mainLine.Color = ROCCols(hid,:);
        h(hid).patch.FaceColor = ROCCols(hid,:);
        h(hid).edge(1).Color = 'none';
        h(hid).edge(2).Color = 'none';
    end
    plot([0 1], [0 1], 'k--')
    xlim([0 1])
    ylim([0 1])
    axis equal tight
    xlabel('False positives')
    ylabel('Hits')
    offsetAxes
    makepretty

    % Plot AUC
    subplot(4,numel(FPNames),2*numel(FPNames)+fpIdx); hold all
    for hid = 1:3
        scatter(hid + 0.1*randn(1,numel(groups)), AUCMatrix{fpIdx}(hid,:), 30, ROCCols(hid,:), 'filled')
    end
    ylim([0 1])
    hline(0.5, 'k--')
    xticks(1:3)
    xticklabels({'M v NM', 'W v M', 'W v NM'})
    xtickangle(45)
    ylabel('AUC')
    offsetAxes
    makepretty

    % Plot stability of AUC with delta days
    subplot(4,numel(FPNames),3*numel(FPNames)+fpIdx); hold all
    for midx = 1:length(UMFiles)
        if ~isempty(FPSum.(FPNameCurr).AUC{midx})
            scatter(deltaDays{midx}(:), mat2vec(FPSum.(FPNameCurr).AUC{midx}(1,:,:)),20,groupColor(groupVector(midx),:),'filled')
        end
    end
    xlabel('\Delta days')
    ylabel('AUC')
    ylim([0 1])
    hline(0.5)
    offsetAxes
    makepretty
end

%% Additional figures

% Plot number of matches as a function of delta days
figure;
hold all
for midx = 1:length(UMFiles)
    scatter(deltaDays{midx}(:), mat2vec(numMatchedUnits{midx}),20,groupColor(groupVector(midx),:),'filled')
end
ylabel('Number of matches')
xlabel('\Deltadays')

% AUC as a function of delta days
figure;
for aucIdx = 1:3
    for fpIdx = 1:numel(FPNames)
        FPNameCurr = FPNames{fpIdx};
        subplot(3,numel(FPNames),(aucIdx-1)*numel(FPNames)+fpIdx);
        hold all
        for midx = 1:length(UMFiles)
            if ~isempty(FPSum.(FPNameCurr).AUC{midx})
                scatter(deltaDays{midx}(:), mat2vec(FPSum.(FPNameCurr).AUC{midx}(aucIdx,:,:)),20,groupColor(groupVector(midx),:),'filled')
            end
        end
        xlabel('\Delta days')
        ylabel('AUC')
        title(sprintf('Fingerprint %s', FPNameCurr))
        ylim([0 1])
        hline(0.5)
    end
end

% Dependence of AUC on number of matched units
figure;
for aucIdx = 1:3
    for fpIdx = 1:numel(FPNames)
        FPNameCurr = FPNames{fpIdx};
        subplot(3,numel(FPNames),(aucIdx-1)*numel(FPNames)+fpIdx);
        hold all
        for midx = 1:length(UMFiles)
            if ~isempty(FPSum.(FPNameCurr).AUC{midx})
                scatter(numMatchedUnits{midx}(:), mat2vec(FPSum.(FPNameCurr).AUC{midx}(aucIdx,:,:)),20,groupColor(groupVector(midx),:),'filled')
            end
        end
        xlabel('Number of matched neurons')
        ylabel('AUC')
        title(sprintf('Fingerprint %s', FPNameCurr))
        ylim([0 1])
        hline(0.5)
    end
end